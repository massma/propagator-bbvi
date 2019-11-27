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

instance VariationalLogic NormalDistribution where
  difference x1 x2 = sqrt (f mean + f stdDev)
    where
      f g = (g x1 - g x2) ^ 2

  logProb = logDensity

  paramVector d = V.fromList [mean d, stdDev d]

  fromParamVector v = normalDistr (v V.! 0) (v V.! 1)

  gradLogProb d x = V.fromList [(x - mu) / std, 1 / std ^ 3 * (x - mu) ^ 2 - 1 / std]
    where
      mu = mean d
      std = stdDev d

-- see https://stats.stackexchange.com/questions/404191/what-is-the-log-of-the-pdf-for-a-normal-distribution
instance Propagated (QProp s) where
  merge q1 q2
    | f q1 q2 < 0.01 = Change False q1
    | null (loglike q2) = Change True q2
    | null (loglike q1) = Change False (q1 {loglike = loglike q2})
    | otherwise = Change False (q1 {loglike = (V.zipWith (+) (loglike q1) (loglike q2))})
    where
      f x1 x2 = difference (distribution x1) (distribution x2)

maxStep :: Time
maxStep = T 20000

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
      V.foldl1' (V.zipWith (+)) $
      V.zipWith
        (\s l -> V.map (* (l - logProb dist s)) (gradLogProb dist s))
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


updateQ Q{..} t' =
  Q s' newDist generator V.empty <$> V.replicateM (V.length samples) (genContinuous newDist generator)
  where
    alpha = 0.1 -- from kuckelbier et al
    eta = 0.01 -- coservative, from kuckelbiet et al
    tau = 1.0
    epsilon = 1e-16
    grad = gradient distribution loglike samples
    (s', rho') = rho alpha eta tau epsilon t' grad gradMemory
    newDist = fromParamVector $ V.zipWith3 (\x g r -> x + r*g) (paramVector distribution) grad rho'

updatePropQ t q = 0.0
  watch t $ \(T t') ->
    with q $ \qprop -> updateQ qprop t' >>= \newq -> write q newq

-- how to efficiently watch all of our t's????

-- propagator :: () -- Maybe VariationalProp
propagator xs = runST $ do
  gen <- create
  let prior = normalDistr 0.0 2.0
  q <- known $ Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty V.empty
  (distribution <$>) <$> content q

someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
-- Just (normalDistr 0.0 1.0)
