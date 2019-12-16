module Statistics.BBVI
  ( DistUtil(..)
  , Dist(..)
  , Differentiable(..)
  , PropNode(..)
  , NormalDist(..)
  , normalDistr
  , defaultNormalDist
  , defaultDirichlet
  , dirichlet
  , gradientScore
  , gradientReparam
  , globalMaxStep
  , globalDelta
  )
where

import           Statistics.BBVI.Class
import           Statistics.BBVI.Propagator
import           Statistics.BBVI.StepSize
import           Statistics.BBVI.Gradient
import           Statistics.BBVI.Defaults
import           Statistics.BBVI.Observed
import           Statistics.BBVI.Distribution.Normal
import           Statistics.BBVI.Distribution.Dirichlet
