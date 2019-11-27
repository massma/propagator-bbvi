module Lib
    ( someFunc
    ) where

import Prelude
import Data.Propagator
import GHC.Generics (Generic)
import Control.Monad (replicateM)
import Control.Monad.ST
import Control.Monad.Primitive
import System.Random.MWC (create, GenST)
import Statistics.Distribution
import Statistics.Distribution.Normal

newtype Time = T Int deriving (Show, Eq, Ord, Read)

type S = Double
type LogLikelihood = Double

type QDistribution = NormalDistribution
data QProp s = Q { gradMemory :: !S
                 , distribution :: !QDistribution
                 , generator :: !(GenST s)
                 , loglike :: !Double
                 }

class VariationalLogic a where
  difference :: a -> a -> Double

  logProb :: a -> Double -> Double


instance VariationalLogic NormalDistribution where
  difference x1 x2 = sqrt (f mean + f stdDev)
    where
      f g = (g x1 - g x2) ^ 2

  logProb = logDensity

instance Propagated (QProp s) where
  merge q1 q2
    | f q1 q2 < 0.01 = Change False q1
    | loglike q2 == 0 = Change True q2
    | otherwise = Change False (q1 {loglike = (loglike q1 + loglike q2)})
    where
      f x1 x2 = difference (distribution x1) (distribution x2)

maxStep :: Time
maxStep = T 20000

instance Propagated Time where
  merge t1 t2
    | t1 > maxStep = Change False t1
    | t2 > t1 = Change True t2
    | otherwise = Change False t1

normalProp std prior xs t q = do
  watch q $ \qprop ->
    write q (qprop {loglike = sum (fmap (logProb (distribution qprop)) xs)}) >>
    with t (\(T t') -> write t (T (t' + 1)))

updateQ t q = 0.0  -- this is todo, make sure to zero log likelihood and calcualted needed things.

-- propagator :: () -- Maybe VariationalProp
propagator xs = runST $ do
  gen <- create
  let prior = normalDistr 0.0 2.0
  q <- known $ Q 0.0 (normalDistr 0.0 2.0) gen 0.0
  (distribution <$>) <$> content q

someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
-- Just (normalDistr 0.0 1.0)
