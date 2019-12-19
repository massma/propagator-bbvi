{-# LANGUAGE FlexibleContexts #-}
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
mixtureLikeReparam
  :: (Dist b SampleVector, Differentiable c Double) -- Dist a SampleVector,
  => Int
  -> Double
  -> GenST s
  -> V.Vector (V.Vector Double)
  -> GradientParams b
  -> GradientParams c
  -> PropNode b
  -> V.Vector (V.Vector (PropNode c))
  -> ST
       s
       (PropNode b, V.Vector (V.Vector (PropNode c)))
mixtureLikeReparam nSamp std gen obss thetaG betaG thetaN betasN = do
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
                let
                  (likes', gradss) = V.unzip $ V.zipWith
                    (\ob betas -> grad'
                      (\bs -> logSum
                        (V.map (auto . log) th)
                        (V.map
                          (\mu -> diffableNormalLogProb mu (auto std) (auto ob))
                          bs
                        )
                      )
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
  -- xss    = dist xsN
  nFactor = fromIntegral $ V.length obss
  theta   = dist thetaN
  betass  = V.map (V.map dist) betasN

mixtureLikeScore nSamp std gen obss thetaG betaG thetaN betasN = do
  -- obss        <- V.replicateM nObs (resample xs gen)
  thetaSamp   <- V.replicateM nSamp (resample theta gen)
  betaSamples <- V.replicateM nSamp
                              (V.mapM (V.mapM (\b -> resample b gen)) betass)
  let likes = V.zipWith
        (\bss th -> V.foldl1' (V.zipWith (+)) $ V.map
          (\obs -> V.zipWith
            (\ob bs -> logSum
              (V.map log th)
              (V.map (\mu -> logProb (normalDistr mu std) ob) bs)
            )
            obs
            bss
          )
          (obss :: V.Vector (V.Vector Double))
        )
        betaSamples
        thetaSamp
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
  -- qTheta <- known $ defaultPropNode (dirichlet (V.replicate nState 1.0))
  qTheta <- cellWith $ mergeGeneric globalDelta globalMaxStep
  write qTheta $ defaultPropNode (dirichlet (V.replicate nState 1.0))
  let betaGrad =
        GParams (fromIntegral $ V.length xs) priorBeta (rhoKuc defaultKucP) --
  -- qBetas <- known =<< V.replicateM
  --   nDimension
  --   (V.generateM
  --     nState
  --     (\_i -> do
  --       mu <- resample priorBeta gen1
  --       return $ defaultPropNode (normalDistr mu (1.0 :: Double))
  --     )
  --   )
  qBetas <- cellWith $ mergeGenericss globalDelta globalMaxStep
  write qBetas =<< V.replicateM
    nDimension
    (V.generateM
      nState
      (\_i -> do
        mu <- resample priorBeta gen1
        return $ defaultPropNode (normalDistr mu (1.0 :: Double))
      )
    )
  stepTogether (mixtureLikeScore nSamp 1.0 gen1 xs thetaGrad betaGrad)
               qTheta
               qBetas
  -- (\tP bPs -> watch tP $ \theta' -> with bPs $ \betas' -> do
  --     (upT, upB) <- mixtureLikeScore nSamp
  --                                    1.0
  --                                    gen1
  --                                    xs
  --                                    thetaGrad
  --                                    betaGrad
  --                                    theta'
  --                                    betas'
  --     write bPs upB
  --     write tP  upT
  --   )
  --   qTheta
  --   qBetas
  -- (\tP bPs0 -> watch tP $ \theta' -> do
  --     bPs <-
  --       known =<< (V.map (V.map (initLocal localStep)) <$> unsafeContent bPs0)
  --     watch bPs $ \betas' -> do
  --       (_upTh, upB) <- mixtureLikeReparam nSamp 1.0 gen1 xs theta' betas'
  --       write bPs upB
  --     bPsNew <- unsafeContent bPs
  --     write bPs0 bPsNew
  --   )
  --   qTheta
  --   qBetas
  -- (\tP0 bPs -> watch bPs $ \betas' -> do
  --     tP <- known =<< (initLocal localStep <$> unsafeContent tP0)
  --     watch tP
  --       $ (\theta' -> do
  --           (upTh, _upB) <- mixtureLikeScore nSamp 1.0 gen1 xs theta' betas'
  --           write tP upTh
  --         )
  --     tPNew <- unsafeContent tP
  --     write tP0 tPNew
  --   )
  --   qTheta
  --   qBetas
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
