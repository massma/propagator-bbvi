module Lib
    ( someFunc
    ) where

import Prelude
import Data.Propagator
import GHC.Generics (Generic)
import Control.Monad (replicateM)
import Control.Monad.ST
import Control.Monad.Primitive
import qualified Data.Vector as V
import System.Random.MWC (create, GenST)
import Statistics.Distribution
import Statistics.Distribution.Normal

newtype Time = T Int deriving (Show, Eq, Ord, Read)

time (T t) = t

type S = (V.Vector Double)
type LogLikelihood = Double


-- | Q: should we store prior with QProp? - then would have to pass xs to updateQ
type QDistribution = NormalDistribution
data QProp s = Q { gradMemory :: !S
                 , distribution :: !QDistribution
                 , generator :: !(GenST s)
                 , loglike :: !(V.Vector Double)
                 , samples :: !(V.Vector Double)
                 }

class VariationalLogic a where
  difference :: a -> a -> Double

  logProb :: a -> Double -> Double

  paramVector :: a -> V.Vector Double

  fromParamVector :: V.Vector Double -> a

  gradLogProb :: a -> Double -> V.Vector Double

  nParams :: a -> Int

instance VariationalLogic NormalDistribution where
  difference x1 x2 = sqrt (f mean + f stdDev)
    where
      f g = (g x1 - g x2) ^ 2

  logProb = logDensity

  paramVector d = V.fromList [mean d, stdDev d]

  fromParamVector v = normalDistr (v V.! 0) (v V.! 1)

  nParams _d = 2

  gradLogProb d x = V.fromList [(x - mu) / std, 1 / std ^ 3 * (x - mu) ^ 2 - 1 / std]
    where
      mu = mean d
      std = stdDev d

-- see https://stats.stackexchange.com/questions/404191/what-is-the-log-of-the-pdf-for-a-normal-distribution
instance Propagated (QProp s) where
  merge q1 q2
--    | f q1 q2 < 0.01 = Change False q1
    | null (loglike q2) = Change True q2
    | null (loglike q1) = Change False (q1 {loglike = loglike q2})
    | otherwise = Change False (q1 {loglike = (V.zipWith (+) (loglike q1) (loglike q2))})
    where
      f x1 x2 = difference (distribution x1) (distribution x2)

maxStep :: Time
maxStep = T 10

instance Propagated Time where
  merge t1 t2
    | t1 > maxStep = Change False t1
    | t2 > t1 = Change True t2
    | otherwise = Change False t1

-- | don't forget prior with variational distribution
normalProp prior std xs t q = do
  watch q $ \qprop ->
    write
      q
      (qprop
         { loglike =
             fmap
               (\mu ->
                  logDensity prior mu +
                  (sum $ fmap (logDensity (normalDistr mu std)) xs))
               (samples qprop)
         }) >>
    with t (\(T t') -> write t (T (t' + 1)))

gradient dist like samples = V.map (/ (fromIntegral $ V.length summed)) summed
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l -> V.map (* (l - logProb dist s)) (gradLogProb dist s))
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


updateQ nSamp Q{..} t' =
  Q s' newDist generator V.empty <$> V.replicateM nSamp (genContinuous newDist generator)
  where
    alpha = 0.1 -- from kuckelbier et al
    eta = 0.01 -- coservative, from kuckelbiet et al
    tau = 1.0
    epsilon = 1e-16
    grad = gradient distribution loglike samples
    (s', rho') = rho alpha eta tau epsilon t' grad gradMemory
    newDist = fromParamVector $ V.zipWith3 (\x g r -> x + r*g) (paramVector distribution) grad rho'
-- >>> create >>= \gen -> updateQ 10 (Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty V.empty) 1

updatePropQ nSamp t q =
  watch t $ \(T t') ->
    -- with q $ \qprop -> write q qprop
    with q $ \qprop -> updateQ nSamp qprop t' >>= \newq -> write q newq

-- how to efficiently watch all of our t's????

-- propagator :: () -- Maybe VariationalProp
propagator xs = runST $ do
  gen <- create
  let prior = normalDistr 0.0 2.0
  q <- known $ Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty V.empty
  t <- known $ T 1
  updatePropQ 10 t q
  normalProp prior 1.0 xs t q
  q' <- content q
  t' <- content t
  return (distribution <$> q', time <$> t', samples <$> q', loglike <$> q')

newtype Test = Test Int deriving (Eq, Ord, Read, Show)

instance Propagated Test where
  merge t1 t2
    | t2 >= (Test 30) = Change False t2
    | t1 >= (Test 5) = Change False t2
    | otherwise = Change True t2

testProp = runST $ do
  x <- known $ Test 0
  lift1 (\(Test x) -> Test 100) x x
  content x
-- >>> testProp
-- Just (Test 0)


someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
-- <interactive>:1492:2-9: error: Variable not in scope: someFunc
