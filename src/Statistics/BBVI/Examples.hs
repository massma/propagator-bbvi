{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
module Statistics.BBVI.Examples
  ( genMixture
  , mixtureFit
  -- , genNormal
  -- , normalFit
  -- , genDirichlet
  -- , dirichletFit
  -- , genMixedMem
  -- , mixedMemFit
  )
where

import           Control.Monad.ST
import           Data.Maybe                     ( fromMaybe )
import           Data.Propagator
import           Data.Foldable                  ( concat )
import qualified Data.Vector                   as V
import           Numeric.AD                     ( grad
                                                , grad'
                                                , auto
                                                , diff
                                                )
import           System.Random.MWC              ( create
                                                , uniform
                                                , initialize
                                                , GenST
                                                )
import qualified System.Random.MWC.Distributions
                                               as MWCD
import qualified Statistics.Sample             as Samp
import           Statistics.BBVI

--- helper functions

logSum v1 = log . minimumProb . V.sum . V.map exp . V.zipWith (+) v1

minimumProb :: (Floating a, Ord a) => a -> a
minimumProb = max 1e-308

--- useful global params
globalMaxStep :: Int
globalMaxStep = 1000
globalDelta :: Double
globalDelta = 1e-16 -- 0.00001 --

nDimension :: Int
nDimension = 1
nTime :: Int
nTime = 1000
nState :: Int
nState = 3 -- 10

-----------------
--  mixture
------------------
mixtureLike :: (Floating a, Ord a) => a -> V.Vector a -> a -> V.Vector a -> a
mixtureLike std theta ob betas = logSum
  (V.map log theta)
  (V.map (\mu -> diffableNormalLogProb mu std ob) betas)

updateReparam
  :: (Dist b SampleVector, Differentiable c Double) -- Dist a SampleVector,
  => Int
  -> (  forall e
      . (Floating e, Ord e)
     => V.Vector e
     -> e
     -> V.Vector e
     -> e
     )
  -> GenST s
  -> V.Vector (V.Vector Double)
  -> GradientParams b
  -> GradientParams c
  -> PropNode b
  -> V.Vector (V.Vector (PropNode c))
  -> ST
       s
       (PropNode b, V.Vector (V.Vector (PropNode c)))
updateReparam nSamp likeF gen obss thetaG betaG thetaN betasN = do
  -- obss        <- V.replicateM nObs (resample xss gen)
  thetaSamp   <- V.replicateM nSamp (resample theta gen)
  betaSamples <- V.replicateM nSamp (epsilon ((betass V.! 0) V.! 0) gen)
  let
    (likes, gradLikes) = V.unzip $ V.zipWith
      (\eps th ->
        V.foldl1'
            (\(x1, y1) (x2, y2) -> (x1 + x2, V.zipWith (V.zipWith (+)) y1 y2))
          $ V.map
              (\obs ->
                let (likes', gradss) = V.unzip $ V.zipWith
                      (\ob betas -> grad'
                        (likeF (V.map auto th) (auto ob))
                        (V.map (\d -> transform d eps) betas)
                      )
                      obs
                      betass
                in  (V.sum likes', gradss)
              )
              (obss :: V.Vector (V.Vector Double))
      )
      betaSamples
      thetaSamp
  return
    $ ( gradientScore thetaG thetaN (nFactor, likes, thetaSamp)
      , V.imap
        (\dim v -> V.imap
          (\c d -> gradientReparam
            betaG
            d
            (nFactor, V.map ((V.! c) . (V.! dim)) gradLikes, betaSamples)
          )
          v
        )
        betasN
      )
 where
  nFactor = fromIntegral $ V.length obss
  theta   = dist thetaN
  betass  = V.map (V.map dist) betasN

mixtureLikeScore nSamp likeF gen obss thetaG betaG thetaN betasN = do
  -- obss        <- V.replicateM nObs (resample xs gen)
  thetaSamp   <- V.replicateM nSamp (resample theta gen)
  betaSamples <- V.replicateM nSamp
                              (V.mapM (V.mapM (\b -> resample b gen)) betass)
  let likes = V.zipWith
        (\theta' betass' ->
          V.foldl1' (V.zipWith (+))
            . V.map (\obs -> V.zipWith (likeF theta') obs betass')
            $ obss
        )
        thetaSamp
        betaSamples
  return
    $ ( gradientScore thetaG thetaN (nFactor, V.map V.sum likes, thetaSamp)
      , V.imap
        (\dim v -> V.imap
          (\c d -> gradientScore
            betaG
            d
            ( nFactor
            , V.map (V.! dim) likes
            , V.map ((V.! c) . (V.! dim)) betaSamples
            )
          )
          v
        )
        betasN
      )
 where
  -- xs     = dist xsN
  nFactor = fromIntegral $ V.length obss
  theta   = dist thetaN
  betass  = V.map (V.map dist) betasN

genMixture :: ST s (V.Vector (V.Vector Double))
genMixture = do
  gen <- create
  let theta' = MWCD.categorical (V.replicate nState 1.0) gen
  let std    = 1.0
  let
    mixtures = V.generate
      nState
      (\k ->
        (MWCD.normal
          ((fromIntegral k - (fromIntegral (nState - 1) * 0.5)) * 5)
          std
          gen
        )
      )
  xs <- V.replicateM nTime (V.replicateM nDimension . (mixtures V.!) =<< theta')
  return xs

mixtureFit xs = runST $ do
  genG <- create
  gen1 <- initialize =<< V.replicateM 256 (uniform genG)
  let priorTheta = dirichlet (V.replicate nState 1.0)
  let priorBeta = normalDistr 0 (5 * (fromIntegral nState - 1) :: Double)
  let nSamp      = 10
  let localStep  = 20
  let thetaGrad =
        GParams (fromIntegral $ V.length xs) priorTheta (rhoKuc defaultKucP) --
  qTheta <- cellWith $ mergeGeneric globalMaxStep globalDelta
  write qTheta $ defaultPropNode (dirichlet (V.replicate nState 1.0))
  let betaGrad =
        GParams (fromIntegral $ V.length xs) priorBeta (rhoKuc defaultKucP) --
  qBetas <- cellWith $ mergeGenericss globalMaxStep globalDelta
  write qBetas =<< V.replicateM
    nDimension
    (V.generateM
      nState
      (\_i -> do
        mu <- resample priorBeta gen1
        return $ defaultPropNode (normalDistr mu (1.0 :: Double))
      )
    )
  stepTogether
    (updateReparam nSamp (mixtureLike 1.0) gen1 xs thetaGrad betaGrad)
    qTheta
    qBetas
  -- stepSeparate localStep
  --              globalDelta
  --              (updateScore nSamp (mixtureLike 1.0) gen1 xs thetaGrad betaGrad)
  --              qTheta
  --              qBetas
  thetaF <- unsafeContent qTheta
  betaF  <- unsafeContent qBetas
  let betaDists = V.map (V.map dist) betaF
  return
    ( alphas (dist thetaF)
    , time thetaF
    , Samp.mean . V.map (Samp.mean . V.map (fromIntegral . time)) $ betaF
    , fmap (\c -> Samp.mean . V.map (mean . (V.! c)) $ betaDists)
           [0 .. nState - 1]
    -- , thetaF
    , fmap (\c -> Samp.mean . V.map (stdDev . (V.! c)) $ betaDists)
           [0 .. nState - 1]
    -- , betaF
    )
