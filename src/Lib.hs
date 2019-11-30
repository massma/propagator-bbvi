{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Lib
  ( someFunc
  ) where
import Prelude
import Data.Propagator
import GHC.Generics (Generic)
import Control.Monad (replicateM, when)
import Control.Monad.ST
import Control.Monad.Primitive
import qualified Data.Vector as V
import Data.Word (Word32)
import Numeric.AD (grad', grad, diff, Mode, auto)
import System.Random.MWC (create, GenST, uniform, initialize)
import Statistics.Distribution
import Statistics.Distribution.Normal

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

newtype Time = T Int deriving (Show, Eq, Ord, Read, Num)

time (T t) = t

maxStep :: Time
maxStep = T 100000

-- | time can only increase in increments of 1. syntax is to write t
-- to 1, b/c that way an empty t starts at 1
instance Propagated Time where
  merge t1 t2
    | t1 >= maxStep = Change False t1
    | otherwise = Change True (T (time t1 + 1))

type S = (V.Vector Double)

type Count = Int
type LogLikelihood = (Time, Count, Count, (V.Vector Double))

instance Propagated LogLikelihood where
  merge (t, cnt, maxCt, v1) (_, _, _, v2)
    | (cnt == 0) && (newCnt == maxCt)  = Change True ((t+1), 0, maxCt, v2)
    | (cnt == 0) = Change True ((t+1), newCnt, maxCt, v2)
    | cnt == maxCt = Change True (t+(T 1), 0, maxCt, (V.zipWith (+) v1 v2))
    | t >= maxStep = Change False (t, cnt, maxCt, v1)
    | otherwise = Change True (t, newCnt, maxCt, (V.zipWith (+) v1 v2))
    where
      newCnt = cnt + 1

-- | Q: should we store prior with QProp? - then would have to pass xs to updateQ
type QDistribution = NormalDistribution

data QProp s = Q { gradMemory :: !S
                 , distribution :: !QDistribution
                 , generator :: !(GenST s)
                 , stdSamples :: !(V.Vector Double)
                 , samples :: !(V.Vector Double)
                 }

instance Propagated (QProp s) where
  merge q1 q2
    | f q1 q2 < 0.00001 = Change False q1
    | otherwise = Change True q2
    where
      f x1 x2 = difference (distribution x1) (distribution x2)

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

instance Differentiable NormalDistribution where
  gradZQ d z = -(z - mean d)/(2 * stdDev d ** 2)
  gradNuTransform d epsilon = V.fromList [1.0 , epsilon]

gradient dist like samples = V.map (/ (fromIntegral $ V.length summed)) summed
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l -> V.map (* (l - logProb dist s)) (gradNuLogQ dist s))
        samples
        like

gradientAD dist like transformedSamples samples = V.map (/ (fromIntegral $ V.length summed)) summed
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith3
        (\tS s l -> V.map (* (l - gradZQ dist tS)) (gradNuTransform dist s))
        transformedSamples
        samples
        like
-- >>> gradient (normalDistr 0.0 2.0) V.empty V.empty
-- [0.0,0.0]

rho alpha eta tau epsilon t grad s0 =
  ( s1
  , V.map
      (\s ->
         eta * (fromIntegral t) ** (negate 0.5 + epsilon) *
         (1.0 / (tau + sqrt s)))
      s1)
  where
    s1 = V.zipWith (\g s -> alpha * g ^ 2 + (1.0 - alpha) * s) grad s0
-- >>> rho 0.1 0.01 1.0 1e-16 1 (V.fromList [0.0, 0.0]) (V.fromList [0.0, 0.0])
-- ([0.0,0.0],[1.0e-2,1.0e-2])

updateQ gradF nSamp Q{..} (T t) l =
  Q s' newDist generator <$> V.replicateM nSamp (genContinuous newDist generator)
  where
    alpha = 0.1 -- from kuckelbier et al
    eta = 0.1 -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
    tau = 1.0
    epsilon = 1e-16
    grad = gradF distribution l samples
    (s', rho') = rho alpha eta tau epsilon t grad gradMemory
    newDist = fromParamVector $ V.zipWith3 (\x g r -> x + r*g) (paramVector distribution) grad rho'
-- >>> create >>= \gen -> updateQ 10 (Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty V.empty) 1

updatePropQ nSamp t q lAD l =
  watch l $ \(tl, _cnt, _maxCnt, l') ->
    with t $ \tGlobal -> when (tl == tGlobal) $
      with q $ \q' ->
          updateQ gradient nSamp q' tGlobal l' >>= \newq -> write t (T 1) >> write q newq

updatePropQAD nSamp t q lAD l =
  watch l $ \(tl, _cnt, _maxCnt, l') ->
    with t $ \tGlobal -> when (tl == tGlobal) $
      with q $ \q' ->
          updateQ gradientAD nSamp q' tGlobal l' >>= \newq -> write t (T 1) >> write q newq

-- | don't forget prior with variational distribution - thin kmaybe we
-- should icnorporate prior into QProp and use it when we update q
-- TODO: consider setting explicit minimum on stddev

normalPropAD prior std xs q l = do
  watch q $ \qprop ->
    write
      l
      ( (T 1)
      , (1 :: Int)
      , (1 :: Int)
      , (fmap
           (\mu ->
              gradZQ prior mu + -- TODO: move prior updateQ
              (sum $ fmap (gradZQ (normalDistr mu std)) xs))
           (samples qprop)))

normalProp prior std xs q l = do
  watch q $ \qprop ->
    write
      l
      ( (T 1)
      , (1 :: Int)
      , (1 :: Int)
      , (fmap
           (\mu ->
              logDensity prior mu + -- TODO: move prior updateQ
              (sum $ fmap (logDensity (normalDistr mu std)) xs))
           (samples qprop)))

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
  q1 <- known $ Q (V.fromList [0.0, 0.0]) qDist gen1 initSamp
  t1 <- known $ T 1
  l1 <- cell
  normalProp prior 1.0 xs q1 l1
  updatePropQ nSamp t1 q1 l1
  q1' <- content q1
  t1' <- content t1
  return (distribution <$> q1', time <$> t1')

someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
--
