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
unsafeContent = (fromMaybe (error "impos") <$>) . content
initLocal step p = p { maxStep = time p + step }

logSum v1 = log . minimumProb . V.sum . V.map exp . V.zipWith (+) v1

minimumProb :: (Floating a, Ord a) => a -> a
minimumProb = max 1e-308

--- useful global params
globalMaxStep :: Int
globalMaxStep = 100000

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
mixtureLike
  :: (Dist a SampleVector, Dist b SampleVector, Differentiable c Double)
  => Int
  -> Int
  -> Double
  -> GenST s
  -> PropNode a
  -> PropNode b
  -> V.Vector (V.Vector (PropNode c))
  -> ST s (PropNode b, V.Vector (V.Vector (PropNode c)))
mixtureLike nSamp nObs std gen xsN thetaN betasN = do
  obss        <- V.replicateM nObs (resample xss gen)
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
    $ ( gradientScore thetaN (1, likes, thetaSamp)
      , V.imap
        (\dim v -> V.imap
          (\c d -> gradientReparam
            d
            (1, V.map ((V.! c) . (V.! dim)) gradLikes, betaSamples)
          )
          v
        )
        betasN
      )
 where
  xss    = dist xsN
  theta  = dist thetaN
  betass = V.map (V.map dist) betasN

mixtureLikeScore nSamp nObs std gen xsN thetaN betasN = do
  obss        <- V.replicateM nObs (resample xs gen)
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
    $ ( gradientScore thetaN (1, V.map V.sum likes, thetaSamp)
      , V.imap
        (\dim v -> V.imap
          (\c d -> gradientScore
            d
            (1, V.map (V.! dim) likes, V.map ((V.! c) . (V.! dim)) betaSamples)
          )
          v
        )
        betasN
      )
 where
    -- obs = xsN
  xs     = dist xsN
  theta  = dist thetaN
  betass = V.map (V.map dist) betasN

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
  let nObs       = 100
  let localStep  = 20
  qBetas <- known =<< V.replicateM
    nDimension
    (V.generateM
      nState
      (\i -> do
        mu <- resample priorBeta gen1
        return
          (defaultNormalDist { dist    = normalDistr mu (1.0 :: Double)
                             , maxStep = globalMaxStep
                             , delta   = globalDelta
                             , prior   = priorBeta
                             , weight  = 1
                             }
          )
      )
    )
  qThetas <-
    -- resample (dirichlet (V.replicate nClusters 1)) gen1 >>= \startTh ->
    ( known
    $ ((defaultDirichlet priorTheta) { maxStep = globalMaxStep
                                     , delta   = globalDelta -- 0.0000001
                                     , dist    = dirichlet
                                       (V.replicate nState 1.0) -- startTh
                                     , weight  = 1
                                     }
      )
    )
  xProp <- known $ defaultObs xs
  -- (\tP bPs xP ->
  --    watch tP $ \theta' ->
  --      with xP $ \xs' ->
  --        with bPs $ \betas' -> do
  --          (upTh, upB) <-
  --            mixtureLike nSamp nObs 1.0 gen1 xs' theta' betas'
  --          write bPs upB
  --          write tP upTh)
  --   qThetas
  --   qBetas
  --   xProp
  (\tP bPs0 xP -> watch tP $ \theta' ->
      with xP
        $ (\xs' -> do
            bPs <-
              known
                =<< (V.map (V.map (initLocal localStep)) <$> unsafeContent bPs0)
            watch bPs $ \betas' -> do
              (_upTh, upB) <- mixtureLikeScore nSamp
                                               nObs
                                               1.0
                                               gen1
                                               xs'
                                               theta'
                                               betas'
              write bPs upB
            bPsNew <- unsafeContent bPs
            write bPs0 bPsNew
          )
    )
    qThetas
    qBetas
    xProp
  (\tP0 bPs xP -> watch bPs $ \betas' ->
      with xP
        $ (\xs' -> do
            tP <- known =<< (initLocal localStep <$> unsafeContent tP0)
            watch tP
              $ (\theta' -> do
                  (upTh, _upB) <- mixtureLikeScore nSamp
                                                   nObs
                                                   1.0
                                                   gen1
                                                   xs'
                                                   theta'
                                                   betas'
                  write tP upTh
                )
            tPNew <- unsafeContent tP
            write tP0 tPNew
          )
    )
    qThetas
    qBetas
    xProp
  thetaF <- unsafeContent qThetas
  betaF  <- unsafeContent qBetas
  let betaDists = V.map (V.map dist) betaF
  return
    ( alphas (dist thetaF)
    , time thetaF
    , V.map (V.map time) betaF
    , V.map (V.map (\d -> (mean d, stdDev d))) betaDists
    )
