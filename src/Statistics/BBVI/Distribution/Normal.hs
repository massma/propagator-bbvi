{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
module Statistics.BBVI.Distribution.Normal
  ( NormalDist(..)
  , normalDistr
  , diffableNormalLogProb
  )
where

import qualified Data.Vector                   as V
import           Numeric.MathFunctions.Constants
                                                ( m_sqrt_2_pi )
import           Statistics.BBVI.Class
import           Statistics.BBVI.StepSize
import           Statistics.BBVI.Propagator     ( PropNode(..)
                                                , SampleDouble
                                                )
import           System.Random.MWC.Distributions
                                                ( normal
                                                , standard
                                                )

data NormalDist = ND {mean :: Double, stdDev :: Double} deriving (Show, Eq, Ord, Read)

normalDistr mu std = ND mu (max 1e-10 std)

instance DistUtil NormalDist where
  fromParamVector xs = normalDistr (xs V.! 0) (xs V.! 1)
  toParamVector d = V.fromList [mean d, stdDev d]
  nParams _d = 2

instance Dist NormalDist SampleDouble where
  resample d gen = normal (mean d) (stdDev d) gen
  logProb d x = (-xm * xm / (2 * sd * sd)) - ndPdfDenom
   where
    xm         = x - mean d
    sd         = stdDev d
    ndPdfDenom = log $ m_sqrt_2_pi * sd
  paramGradOfLogQ d x = V.fromList
    [ xm / std ^ (2 :: Int)
    , 1 / std ^ (3 :: Int) * (xm ^ (2 :: Int) - std ^ (2 :: Int))
    ]
   where
    std = stdDev d
    xm  = x - mean d
-- >>> paramGradOfLogQ (normalDistr 0.0 3.0) (2.0 :: Double)
-- [0.2222222222222222,-0.18518518518518517]

diffableNormalLogProb mu sd x = -xm * xm / (2 * sd * sd) - ndPdfDenom
 where
  xm         = x - mu
  ndPdfDenom = log $ sqrt (2 * pi) * sd

instance Differentiable NormalDist SampleDouble where
  transform d samp = mean d + stdDev d * samp
  epsilon _d gen = standard gen
  sampleGradOfLogQ d z = -(z - mean d) / (stdDev d ** 2)
  gradTransform _d samp = V.fromList [1.0, samp] -- ND $ V.fromList [1.0 , eps * stdDev d] -- --
-- >>> (transform (normalDistr 0.0 2.0) 1.0 :: Double)
-- 2.0
-- >>> gradTransform (normalDistr 0.0 1.0) (2.0 :: Double)
-- [1.0,2.0]
