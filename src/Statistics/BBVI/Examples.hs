{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
module Statistics.BBVI.Examples
  ( genMixture
  , mixtureFit
  , genNormal
  , normalFit
  , genDirichlet
  , dirichletFit
  -- , genMixedMem
  -- , mixedMemFit
  )
where

import           Control.Monad.ST
import           Data.Propagator
import qualified Data.Vector                   as V
import           Numeric.AD                     ( grad'
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
import           Statistics.BBVI

--- helper functions
logSum :: (Floating c, Ord c) => V.Vector c -> V.Vector c -> c
logSum v1 = log . minimumProb . V.sum . V.map exp . V.zipWith (+) v1

minimumProb :: (Floating a, Ord a) => a -> a
minimumProb = max 1e-308

--- useful global params
-- for fit
globalMaxStep :: Int
globalMaxStep = 10000
globalDelta :: Double
globalDelta = 1e-16 -- 0.00001 --

-- for generating data
nDimension :: Int
nDimension = 1
nTime :: Int
nTime = 1000
nState :: Int
nState = 3 -- 10

-----------------
--  mixture
------------------
-- | log joint of a mixture model, with categorical marginalized out
mixtureJoint :: (Floating a, Ord a) => a -> V.Vector a -> a -> V.Vector a -> a
mixtureJoint std theta ob betas = logSum
  (V.map log theta)
  (V.map (\mu -> diffableNormalLogProb mu std ob) betas)

-- | transform a log joint to a gradient propagator, for use with
-- models with a mixture-like plate structure. Updates the mixture
-- components with the reparameterization gradient (Kucukelbir et al
-- 2017).
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
  -> DistInvariant b
  -> DistInvariant c
  -> DistCell b
  -> V.Vector (V.Vector (DistCell c))
  -> ST
       s
       (DistCell b, V.Vector (V.Vector (DistCell c)))
updateReparam nSamp jointF gen obss thetaG betaG thetaN betasN = do
  -- obss        <- V.replicateM nObs (resample xss gen)
  thetaSamp   <- V.replicateM nSamp (resample theta gen)
  betaSamples <- V.replicateM nSamp (epsilon ((betass V.! 0) V.! 0) gen)
  let
    (joints, gradJoints) = V.unzip $ V.zipWith
      (\e th ->
        V.foldl1'
            (\(x1, y1) (x2, y2) -> (x1 + x2, V.zipWith (V.zipWith (+)) y1 y2))
          $ V.map
              (\obs ->
                let (joints', gradss) = V.unzip $ V.zipWith
                      (\ob betas -> grad' (jointF (V.map auto th) (auto ob))
                                          (V.map (\d -> transform d e) betas)
                      )
                      obs
                      betass
                in  (V.sum joints', gradss)
              )
              (obss :: V.Vector (V.Vector Double))
      )
      betaSamples
      thetaSamp
  return
    $ ( gradientScore thetaG thetaN (nFactor, joints, thetaSamp)
      , V.imap
        (\dim v -> V.imap
          (\c d -> gradientReparam
            betaG
            d
            (nFactor, V.map ((V.! c) . (V.! dim)) gradJoints, betaSamples)
          )
          v
        )
        betasN
      )
 where
  nFactor = fromIntegral $ V.length obss
  theta   = dist thetaN
  betass  = V.map (V.map dist) betasN

-- | transform a log joint to a gradient propagator, for use with
-- models with a mixture-like plate structure. Updates the mixture
-- components with the score gradient (Ranganath et al 2014)
updateScore
  :: (Dist a1 c1, Dist a2 c2)
  => Int
  -> (c1 -> a3 -> V.Vector c2 -> Double)
  -> GenST s
  -> V.Vector (V.Vector a3)
  -> DistInvariant a1
  -> DistInvariant a2
  -> DistCell a1
  -> V.Vector (V.Vector (DistCell a2))
  -> ST
       s
       ( DistCell a1
       , V.Vector (V.Vector (DistCell a2))
       )
updateScore nSamp jointF gen obss thetaG betaG thetaN betasN = do
  -- obss        <- V.replicateM nObs (resample xs gen)
  thetaSamp   <- V.replicateM nSamp (resample theta gen)
  betaSamples <- V.replicateM nSamp
                              (V.mapM (V.mapM (\b -> resample b gen)) betass)
  let joints = V.zipWith
        (\theta' betass' ->
          V.foldl1' (V.zipWith (+))
            . V.map (\obs -> V.zipWith (jointF theta') obs betass')
            $ obss
        )
        thetaSamp
        betaSamples
  return
    $ ( gradientScore thetaG thetaN (nFactor, V.map V.sum joints, thetaSamp)
      , V.imap
        (\dim v -> V.imap
          (\c d -> gradientScore
            betaG
            d
            ( nFactor
            , V.map (V.! dim) joints
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

-- | generate data for testing mixture model
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

-- | fit a mixture model, given data
mixtureFit
  :: V.Vector (V.Vector Double) -> (Dirichlet, V.Vector (V.Vector NormalDist))
mixtureFit xs = runST $ do
  genG <- create
  gen1 <- initialize =<< V.replicateM 256 (uniform genG)
  let priorTheta = dirichlet (V.replicate nState 1.0)
  let priorBeta = normalDistr 0 (5 * (fromIntegral nState - 1) :: Double)
  let nSamp      = 10
  let localStep  = (20 :: Int)

  -- initialize distribution cells (both invariant and variants)
  let thetaGrad = DistInvariant (fromIntegral $ V.length xs)
                                priorTheta
                                (rhoKuc defaultKucP) --
  qTheta <- cellWith $ mergeGeneric globalMaxStep globalDelta
  write qTheta $ defaultDistCell (dirichlet (V.replicate nState 1.0))
  let betaGrad = DistInvariant (fromIntegral $ V.length xs)
                               priorBeta
                               (rhoKuc defaultKucP) --
  qBetas <- cellWith $ mergeGenericss globalMaxStep globalDelta
  write qBetas =<< V.replicateM
    nDimension
    (V.generateM
      nState
      (\_i -> do
        mu <- resample priorBeta gen1
        return $ defaultDistCell (normalDistr mu (1.0 :: Double))
      )
    )

  -- attach gradient propagators to cells
  stepTogether
    (updateReparam nSamp (mixtureJoint 1.0) gen1 xs thetaGrad betaGrad)
    qTheta
    qBetas
  -- stepSeparate localStep
  --              globalDelta
  --              (updateScore nSamp (mixtureJoint 1.0) gen1 xs thetaGrad betaGrad)
  --              qTheta
  --              qBetas

  -- pull content (run network to quiesence)
  thetaF <- unsafeContent qTheta
  betaF  <- unsafeContent qBetas
  let betaDists = V.map (V.map dist) betaF
  return (dist thetaF, betaDists)


--------------------------------------
-- fit simple normal distribution mean
--------------------------------------
-- | generate normally distributed data for testing
genNormal :: ST s (V.Vector Double)
genNormal = do
  gen <- create
  xs  <- V.replicateM 1000 (resample (normalDistr (5.0 :: Double) 3.0) gen)
  return xs

-- | log joint of normal distribution
jointNormal :: Floating a => a -> a -> a -> a
jointNormal std ob mu = diffableNormalLogProb mu std ob

-- | fit a gaussian using reparam gradient
normalFit :: Dist (Obs a) SampleDouble => V.Vector a -> (NormalDist, Time)
normalFit xs = runST $ do
  genG <- create
  gen1 <- initialize =<< V.replicateM 256 (uniform genG)
  let nSamp = 100
  xProp <- known $ defaultObs xs
  q     <- cellWith $ mergeGeneric globalMaxStep globalDelta
  write q (defaultDistCell $ normalDistr 0.0 2.0)
  let qInvar = DistInvariant (fromIntegral $ V.length xs)
                             (normalDistr 0.0 2.0)
                             (rhoKuc defaultKucP)

  (\qP xP -> watch qP $ \q' -> with xP $ \xs' ->
      singleUpdateReparam nSamp (jointNormal 1.0) gen1 qInvar xs' q' >>= write q
    )
    q
    xProp
  q' <- unsafeContent q
  return (dist q', time q')

-- | transform a log joint to a gradient propagator for a single
-- distribution cell, using reparameterization graidient. TODO: more
-- work on generalizing these types of functions!
singleUpdateReparam
  :: (Dist b SampleDouble, Differentiable c Double)
  => Int
  -> (forall  e . (Floating e, Ord e) => e -> e -> e)
  -> GenST s
  -> DistInvariant c
  -> DistCell b
  -> DistCell c
  -> ST s (DistCell c)
singleUpdateReparam nSamp joint gen qG xsN qN = do
  obs     <- V.replicateM nSamp (resample xs gen)
  samples <- V.replicateM nSamp (epsilon q gen)
  let gradJoint = V.map
        (\mu -> V.sum $ V.map (\ob -> diff (joint (auto ob)) mu) obs)
        (V.map (transform q) samples)
  return $ gradientReparam qG qN ((fromIntegral nSamp), gradJoint, samples)
 where
  xs = dist xsN
  q  = dist qN


---------------------------
-- simple dirichlet example
---------------------------
-- | generate data for testing the fit to a dirichlet
genDirichlet :: ST s (V.Vector Int)
genDirichlet = do
  gen <- create
  xs  <- V.replicateM
    1000
    (   resample (dirichlet (V.fromList [10.0, 20.0])) gen
    >>= \cat -> MWCD.categorical (cat :: V.Vector Double) gen
    )
  return (xs :: V.Vector Int)

-- | git a dirichlet to samples form a categorical
dirichletFit :: V.Vector Int -> (Dirichlet, Time)
dirichletFit xs = runST $ do
  genG <- create
  gen1 <- initialize =<< V.replicateM 256 (uniform genG)
  let priorTheta = dirichlet (V.fromList [1.0, 1.0])
  let nSamp      = 100
  qTheta <- cellWith $ mergeGeneric globalMaxStep globalDelta
  let qInvar = DistInvariant 1 priorTheta (rhoKuc defaultKucP)
  write qTheta (defaultDistCell priorTheta)
  (\tP -> watch tP $ \theta' -> do
      upTh <- dirichletGradProp nSamp gen1 xs qInvar theta'
      write tP upTh
    )
    qTheta
  thetaF <- unsafeContent qTheta
  return (dist thetaF, time thetaF)

dirichletGradProp
  :: Dist a (V.Vector Double)
  => Int
  -> GenST s
  -> V.Vector Int
  -> DistInvariant a
  -> DistCell a
  -> ST s (DistCell a)
dirichletGradProp nSamp gen xs qInvar diriQ =
  V.replicateM nSamp (resample diri gen) >>= \samples -> return $ gradientScore
    qInvar
    diriQ
    (1, V.map (\z -> V.sum $ V.map (\i -> log (z V.! i)) xs) samples, samples)
  where diri = dist diriQ
