{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Simpler
  ( someFunc
  ) where
import Prelude
import Data.Propagator
import Control.Monad.ST
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import System.Random.MWC (create, uniform, initialize)
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

class VariationalLogic a => Differentiable a where
  gradNuTransform :: a -> Double -> V.Vector Double
  gradZQ :: a -> Double -> Double
  transform :: a -> Double -> Double

type Time = Int

maxStep :: Time
maxStep = 100000
-- | time can only increase in increments of 1. syntax is to write t
-- to 1, b/c that way an empty t starts at 1

type Memory = V.Vector Double

type Sample = V.Vector Double

type Weight = Double

data QDist = QDist { time :: Time
                   , memory :: Memory
                   , weight :: Weight
                   , dist :: NormalDistribution
                   , prior :: NormalDistribution
                   , samples :: Sample
                   , standardNormal :: Sample}

instance Propagated QDist where
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

instance Differentiable NormalDistribution where
  gradZQ d z = -(z - mean d)/(stdDev d ** 2)
  gradNuTransform d epsilon = V.fromList [1.0 , epsilon]
  transform d epsilon = mean d + stdDev d * epsilon

gradient :: QDist -> V.Vector Double -> V.Vector Double
gradient QDist{..} like = V.map (/ (fromIntegral $ V.length summed)) summed
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l -> V.map (* (l + weight * (logDensity prior s - logProb dist s))) (gradNuLogQ dist s))
        samples
        like

gradientAD QDist{..} like = V.map (/ (fromIntegral $ V.length summed)) summed
  where
    summed =
      -- TODO: add prior!!
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith3
        (\tS sn l -> V.map (* (l + weight * (gradZQ prior tS - gradZQ dist tS))) (gradNuTransform dist sn))
        samples
        standardNormal
        like

rho alpha eta tau epsilon t grad s0 =
  ( s1
  , V.map
      (\s ->
         eta * (fromIntegral t) ** (negate 0.5 + epsilon) *
         (1.0 / (tau + sqrt s)))
      s1)
  where
    s1 = V.zipWith (\g s -> alpha * g ^ (2 :: Int) + (1.0 - alpha) * s) grad s0

qNormalGeneric likeFunc gradF gen nSamp std xs q@(QDist{..}) = do
  s <- V.replicateM nSamp (MWCD.standard gen)
  return $ QDist (time + 1) memory' weight newQ prior (fmap (transform newQ) s) s
  where
    like = likeFunc std xs samples
    alpha = 0.1 -- from kuckelbier et al
    eta = 0.1 -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
    tau = 1.0
    epsilon = 1e-16
    grad = gradF q like
    (memory', rho') = rho alpha eta tau epsilon time grad memory
    newQ = fromParamVector $ V.zipWith3 (\x g r -> x + r*g) (paramVector dist) grad rho'
-- >>> create >>= \gen -> updateQ 10 (Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty V.empty) 1

-- will content trick work, or should we do watch???? seems liek
-- content could repeat a bunch of identical computations
qNormalPropGeneric likeFunc gradF gen nSamp std xs q =
  watch q $ \q' -> do
    -- q' <- fromMaybe (error "impos") <$> content q
    write q =<< qNormalGeneric likeFunc gradF gen nSamp std xs q'

updateQNormalProp = qNormalPropGeneric normalLike gradient

updateQNormalPropAD = qNormalPropGeneric normalLikeAD gradientAD

normalLike std xs samples =
  fmap (\mu -> (sum $ fmap (logDensity (normalDistr mu std)) xs)) samples

-- | specialize all fmaps to vector???
normalLikeAD std xs samples =
  fmap
    (\mu -> (sum $ V.map ((V.! 0) . gradNuLogQ (normalDistr mu std)) xs))
    samples

-- propagator :: () -- Maybe VariationalProp
propagator xs = runST $ do
  genG <- create
  seed <- V.replicateM 256 (uniform genG)
  gen1 <- initialize seed
  let prior = normalDistr 0.0 2.0
  let nSamp = 100
  let qDist = normalDistr 0.0 2.0
  stdNormal <- V.replicateM nSamp (MWCD.standard gen1)
  let initSamp = V.map (transform qDist) stdNormal
  q <- known $ QDist 1 (V.replicate 2 0) 1 qDist prior initSamp stdNormal
  updateQNormalProp gen1 nSamp 1.0 xs q
  updateQNormalPropAD gen1 nSamp 1.0 xs q
  q' <- fromMaybe (error "impos") <$> content q
  return (dist q', time q', Samp.mean xs)

someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- V.replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
--
