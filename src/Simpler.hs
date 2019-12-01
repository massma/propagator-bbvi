{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Simpler
  ( someFunc
  ) where
import Prelude
import Data.Propagator
import GHC.Generics (Generic)
import Control.Monad (replicateM, when)
import Control.Monad.ST
import Control.Monad.Primitive
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import Data.Word (Word32)
import Numeric.AD (grad', grad, diff, Mode, auto)
import System.Random.MWC (create, GenST, uniform, initialize)
import qualified System.Random.MWC.Distributions as MWCD
import Statistics.Distribution
import Statistics.Distribution.Normal
import qualified Statistics.Sample as Samp

class VariationalLogic a where
  difference :: a -> a -> Double
  logProb :: a -> Double -> Double
  paramVector :: a -> V.Vector Double
  fromParamVector :: V.Vector Double -> a
  gradNuLogQ :: a -> Double -> V.Vector Double -- nu are parameters to variational family
  nParams :: a -> Int

type Time = Int

maxStep :: Time
maxStep = 100000

-- | time can only increase in increments of 1. syntax is to write t
-- to 1, b/c that way an empty t starts at 1

type Memory = V.Vector Double

norm = sqrt . V.sum . V.map (^2)

type Sample = V.Vector Double

type Weight = Double

data QDist = QDist { time :: Time
                   , memory :: Memory
                   , weight :: Weight
                   , dist :: NormalDistribution
                   , prior :: NormalDistribution
                   , samples :: Sample}

instance Propagated QDist where
  merge q1 q2
    | difference (dist q1) (dist q2) < 0.00001 = Change False q1 -- 0.00001
    | otherwise = Change True q2

instance VariationalLogic NormalDistribution where
  difference x1 x2 = sqrt (f mean + f stdDev)
    where
      f g = (g x1 - g x2) ^ 2
  logProb = logDensity
  paramVector d = V.fromList [mean d, stdDev d]
  fromParamVector v = normalDistr (v V.! 0) (v V.! 1)
  nParams _d = 2
  gradNuLogQ d x = V.fromList [(x - mu) / std, 1 / std ^ 3 * (x - mu) ^ 2 - 1 / std]
    where
      mu = mean d
      std = stdDev d

gradient QDist{..} like = V.map (/ (fromIntegral $ V.length summed)) summed
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l -> V.map (* (l + weight * (logDensity prior s - logProb dist s))) (gradNuLogQ dist s))
        samples
        like

rho alpha eta tau epsilon t grad s0 =
  ( s1
  , V.map
      (\s ->
         eta * (fromIntegral t) ** (negate 0.5 + epsilon) *
         (1.0 / (tau + sqrt s)))
      s1)
  where
    s1 = V.zipWith (\g s -> alpha * g ^ 2 + (1.0 - alpha) * s) grad s0

updateQNormal gen nSamp std xs q@(QDist{..}) =
  QDist (time + 1) memory' weight newQ prior <$> V.replicateM nSamp (genContinuous newQ gen)
  where
    like = normalLike std xs samples
    alpha = 0.1 -- from kuckelbier et al
    eta = 0.1 -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
    tau = 1.0
    epsilon = 1e-16
    dft = (V.replicate (nParams dist) 0.0)
    grad = gradient q like
    (memory', rho') = rho alpha eta tau epsilon time grad memory
    newQ = fromParamVector $ V.zipWith3 (\x g r -> x + r*g) (paramVector dist) grad rho'
-- >>> create >>= \gen -> updateQ 10 (Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty V.empty) 1

-- will content trick work, or should we do watch???? seems liek
-- content could repeat a bunch of identical computations
updateQNormalProp gen nSamp std xs q =
  watch q $ \q' -> do
    -- q' <- fromMaybe (error "impos") <$> content q
    write q =<<  updateQNormal gen nSamp std xs q'

normalLike std xs samples =
  fmap (\mu -> (sum $ fmap (logDensity (normalDistr mu std)) xs)) samples

-- propagator :: () -- Maybe VariationalProp
propagator xs = runST $ do
  genG <- create
  seed <- V.replicateM 256 (uniform genG)
  gen1 <- initialize seed
  gen2 <- initialize seed
  let prior = normalDistr 0.0 2.0
  let nSamp = 100
  let qDist = (normalDistr 0.0 2.0)
  initSamp <- V.replicateM nSamp (genContinuous qDist gen1)
  q <- known $ QDist 1 (V.replicate 2 0) 1 qDist prior initSamp
  updateQNormalProp gen1 nSamp 1.0 xs q
  q' <- fromMaybe (error "impos") <$> content q
  return (dist q', time q')

someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- V.replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
--
