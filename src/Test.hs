module Test () where

import Prelude
import Data.Propagator
import Control.Monad.ST
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import System.Random.MWC (create, uniform, initialize)
import qualified System.Random.MWC.Distributions as MWCD
import qualified Statistics.Sample as Samp
import Simpler

------------------
-- Normal
------------------

genNormal :: ST s (V.Vector Double)
genNormal = do
  gen <- create
  xs <- V.replicateM 1000 (resample (normalDistr (5.0 :: Double) 3.0) gen)
  return xs

normalFit xs =
  runST $ do
    genG <- create
    gen1 <- initialize =<< V.replicateM 256 (uniform genG)
    let nSamp = 100
    xProp <- known $ defaultObs xs
    q <-
      known $
      (defaultNormalDist
         { prior = normalDistr 0.0 2.0
         , dist = normalDistr 0.0 2.0
         , weight = fromIntegral nSamp
         })
    (\qP xP ->
       watch qP $ \q' ->
         with xP $ \xs' ->
           normalLikeAD2 nSamp 1.0 gen1 (xs', q') >>= write q)
      q
      xProp
    -- qPropAD2 (normalLikeAD2 nSamp nSamp 1.0 gen1) (xProp, q) -- (V.length xs)
    q' <- fromMaybe (error "impos") <$> content q
    xs' <- fromMaybe (error "impos") <$> content xProp
    return (dist q', time q', Samp.mean xs, time xs')

normalLikeAD nSamp std gen xs q =
  V.replicateM nSamp (epsilon q gen) >>= \samples ->
    return
      ( fromIntegral $ V.length xs
      , V.map
          (\eps ->
             V.sum $
             V.map (mean . paramGradOfLogQ (normalDistr (transform q eps) std)) xs)
          samples
      , samples)

normalLikeAD2 nSamp std gen (xsN, qN) = do
  obs <- V.replicateM nSamp (resample xs gen)
  samples <- V.replicateM nSamp (epsilon q gen)
  let gradLikes =
        V.map
          (\eps ->
             V.sum $
             V.map
               (mean . paramGradOfLogQ (normalDistr (transform q eps) std))
               (obs :: V.Vector Double))
          samples
  return $ gradientReparam qN ((fromIntegral nSamp), gradLikes, samples)
  where
    xs = dist xsN
    q = dist qN

------------------
--  Dirichlet
------------------
genDirichlet = do
  gen <- create
  xs <-
    V.replicateM
      1000
      (resample (dirichlet (V.fromList [10.0, 20.0])) gen >>= \cat ->
         MWCD.categorical (cat :: V.Vector Double) gen)
  return (xs :: V.Vector Int)

dirichletFit xs =
  runST $ do
    genG <- create
    gen1 <- initialize =<< V.replicateM 256 (uniform genG)
    let priorTheta = Diri (V.fromList [1.0, 1.0])
    let nSamp = 100
    let nClusters = 2
    qTheta <-
      -- resample (dirichlet (V.replicate nClusters 1)) gen1 >>= \startTh ->
        (known $
         ((defaultDirichlet priorTheta)
            { maxStep = globalMaxStep
            , delta = globalDelta -- 0.0000001
            , dist = dirichlet (V.fromList [(1.0 :: Double), 1.0]) -- startTh
            , weight = 1
            }))
    (\tP ->
       watch tP $ \theta' -> do
        upTh <- dirichletLike nSamp gen1 xs theta'
        write tP upTh)
      qTheta
    thetaF <- unsafeContent qTheta
    return (dist thetaF, time thetaF)

dirichletLike nSamp gen xs diriQ =
  V.replicateM nSamp (resample diri gen) >>= \samples ->
    return $
      gradientScore
        diriQ
        ( 1
        , V.map
          (\z ->
             V.sum $ V.map (\i -> log (z V.! i)) xs)
          samples
        , samples)
  where
    diri = dist diriQ
