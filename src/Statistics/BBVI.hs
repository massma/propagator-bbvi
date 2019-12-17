module Statistics.BBVI
  ( DistUtil(..)
  , Dist(..)
  , Differentiable(..)
  , PropNode(..)
  , NormalDist(..)
  , normalDistr
  , defaultNormalDist
  , diffableNormalLogProb
  , defaultDirichlet
  , dirichlet
  , alphas
  , defaultObs
  , Obs(..)
  , gradientScore
  , gradientReparam
  , SampleDouble
  , SampleVector
  )
where

import           Statistics.BBVI.Class
import           Statistics.BBVI.Propagator
import           Statistics.BBVI.StepSize
import           Statistics.BBVI.Gradient
import           Statistics.BBVI.Observed
import           Statistics.BBVI.Distribution.Normal
import           Statistics.BBVI.Distribution.Dirichlet
