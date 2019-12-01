{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Simpler
  ( someFunc
  ) where
import Prelude
import Data.Propagator
import Control.Monad.ST
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import System.Random.MWC (create, uniform, initialize, GenST)
import Statistics.Distribution
import Statistics.Distribution.Normal
import qualified System.Random.MWC.Distributions as MWCD
import qualified Statistics.Sample as Samp

class VariationalLogic a where
  difference :: a -> a -> Double
  logProb :: a -> Double -> Double
  paramVector :: a -> V.Vector Double
  fromParamVector :: V.Vector Double -> a
  gradNuLogQ :: a -> Double -> V.Vector Double -- nu are parameters to variational family
  nParams :: a -> Int
  transform :: a -> Double -> Double
  resample :: a -> GenST s -> ST s Double

class VariationalLogic a => Differentiable a where
  gradNuTransform :: a -> Double -> V.Vector Double
  gradZQ :: a -> Double -> Double


type Time = Int

maxStep :: Time
maxStep = 10000 -- 100000
-- | time can only increase in increments of 1. syntax is to write t
-- to 1, b/c that way an empty t starts at 1

type Memory = V.Vector Double

type Sample = V.Vector Double

type Weight = Double

data QDist a = QDist
  { time :: !Time
  , memory :: !Memory
  , weight :: !Weight
  , dist :: !a
  , prior :: !NormalDistribution
  , samples :: !Sample
  , rhoF :: !(V.Vector Double -> QDist a -> (Memory, V.Vector Double))
  }

instance VariationalLogic a => Propagated (QDist a) where
  merge q1 q2
    | difference (dist q1) (dist q2) < 0.00001 = Change False q1 -- 0.00001
    | time q1 >= maxStep = Change False q1
    | otherwise = Change True q2

instance VariationalLogic NormalDistribution where
  difference x1 x2 = sqrt (f mean + f stdDev)
    where
      f g = (g x1 - g x2) ^ (2 :: Int)
  logProb = logDensity
  paramVector d = V.fromList [mean d, stdDev d]
  fromParamVector v = normalDistr (v V.! 0) (v V.! 1)
  nParams _d = 2
  gradNuLogQ d x = V.fromList [(x - mu) / std, 1 / std ^ (3 :: Int) * (x - mu) ^ (2 :: Int) - 1 / std]
    where
      mu = mean d
      std = stdDev d
  transform d epsilon = mean d + stdDev d * epsilon
  resample _d gen = MWCD.standard gen

instance Differentiable NormalDistribution where
  gradZQ d z = -(z - mean d)/(stdDev d ** 2)
  gradNuTransform _d epsilon = V.fromList [1.0 , epsilon]


gradient :: VariationalLogic a => QDist a -> V.Vector Double -> V.Vector Double
gradient QDist{..} like = V.map (/ (fromIntegral $ V.length summed)) summed
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l -> V.map (* (l + weight * (logProb prior s - logProb dist s))) (gradNuLogQ dist s))
        samples
        like

gradientAD :: Differentiable a => QDist a -> V.Vector Double -> V.Vector Double
gradientAD QDist{..} like = V.map (/ (fromIntegral $ V.length summed)) summed
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l -> V.map (* (l + weight * (gradZQ prior (transform dist s) - gradZQ dist (transform dist s)))) (gradNuTransform dist s))
        samples
        like

rho alpha eta tau epsilon grad QDist{..} =
  ( memNew
  , V.map
      (\s ->
         eta * (fromIntegral time) ** (negate 0.5 + epsilon) *
         (1.0 / (tau + sqrt s)))
      memNew)
  where
    memNew = V.zipWith (\g s -> alpha * g ^ (2 :: Int) + (1.0 - alpha) * s) grad memory

qGeneric gradF likeFunc gen nSamp xs q@(QDist{..}) = do
  s <- V.replicateM nSamp (resample newQ gen)
  return $ q {time = time + 1, memory = memory', dist = newQ, samples = s} -- QDist (time + 1) memory' weight newQ prior s rhoF
  where
    like = V.map (\mu -> (sum $ V.map (likeFunc (transform dist mu)) xs)) samples
    grad = gradF q like
    (memory', rho') = rhoF grad q
    newQ = fromParamVector $ V.zipWith3 (\x g r -> x + r*g) (paramVector dist) grad rho'
-- >>> create >>= \gen -> updateQ 10 (Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty V.empty) 1

qPropGeneric ::
  VariationalLogic a =>
     (QDist a -> V.Vector Double -> V.Vector Double)
  -> (Double -> Double -> Double)
  -> GenST s
  -> Int
  -> V.Vector Double
  -> Cell s (QDist a)
  -> ST s ()
qPropGeneric gradF likeFunc gen nSamp xs q =
  watch q $ \q' -> do
    write q =<< qGeneric gradF likeFunc gen nSamp xs q'

qProp ::
  VariationalLogic a =>
  (Double -> Double -> Double)
  -> GenST s
  -> Int
  -> V.Vector Double
  -> Cell s (QDist a)
  -> ST s ()
qProp = qPropGeneric gradient

qPropAD ::
  Differentiable a =>
  (Double -> Double -> Double)
  -> GenST s
  -> Int
  -> V.Vector Double
  -> Cell s (QDist a)
  -> ST s ()
qPropAD = qPropGeneric gradientAD

qNormalProp :: VariationalLogic a => Double -> GenST s -> Int -> V.Vector Double -> Cell s (QDist a) -> ST s ()
qNormalProp std = qProp (normalLike std)

qNormalPropAD :: Differentiable a => Double -> GenST s -> Int -> V.Vector Double -> Cell s (QDist a) -> ST s ()
qNormalPropAD std = qPropAD (normalLikeAD std)

normalLike std z x = logDensity (normalDistr z std) x

-- | specialize all fmaps to vector???
normalLikeAD std z = (V.! 0) . gradNuLogQ (normalDistr z std)

propagator xs = runST $ do
  genG <- create
  seed <- V.replicateM 256 (uniform genG)
  gen1 <- initialize seed
  let prior = normalDistr 0.0 2.0
  let nSamp = 100
  let qDist = normalDistr 0.0 2.0
  stdNormal <- V.replicateM nSamp (MWCD.standard gen1)
  let alpha = 0.1 -- from kuckelbier et al
  let eta = 0.1 -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
  let tau = 1.0
  let epsilon = 1e-16
  q <- known $ QDist 1 (V.replicate 2 0) 1 qDist prior stdNormal (rho alpha eta tau epsilon)
  -- qNormalProp 1.0 gen1 nSamp xs q
  qNormalPropAD 1.0 gen1 nSamp xs q
  q' <- fromMaybe (error "impos") <$> content q
  return (dist q', time q', Samp.mean xs)

someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- V.replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
--
