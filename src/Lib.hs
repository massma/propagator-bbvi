module Lib
    ( someFunc
    ) where

import Prelude
import Data.Propagator
import GHC.Generics (Generic)
import Control.Monad.ST

type Step = Int
type S = Double

newtype VariationalProp = V (Step, S, Distribution) deriving (Show, Read, Eq, Generic)

newtype Distribution = D NormalDistribution deriving (Show, Read, Eq, Generic)

data NormalDistribution = ND
  { mean :: !Double
  , std :: !Double
  , variance :: !Double
  } deriving (Show, Read, Eq, Generic)

normalDistr :: Double -> Double -> NormalDistribution
normalDistr m s = ND m s (s * s)

instance Propagated VariationalProp where
  merge (V (i1, s1, d1)) (V (i2, s2, d2))
    | i1 >= i2 = Change False (V (i2, s1, d1))
    | otherwise = Change True (V (i2, s2, d2))

-- propagator :: () -- Maybe VariationalProp
propagator = runST $ do
  x <- known $ V (0, 0.0, D (normalDistr 0.0 1.0))
  lift1 (\(V (i, s, d)) -> if i == 5 then V (i, s, d) else V (i + 1, s, d)) x x
  content x

someFunc :: IO ()
someFunc = putStrLn (show propagator)
-- >>> someFunc
-- Just (V (5,0.0,D (ND {mean = 0.0, std = 1.0, variance = 1.0})))
