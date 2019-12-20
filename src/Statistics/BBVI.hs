module Statistics.BBVI
  ( DistUtil(..)
  , Dist(..)
  , Differentiable(..)
  , PropNode(..)
  , PropNodes
  , PropNodess
  , NormalDist(..)
  , normalDistr
  , defaultPropNode
  , diffableNormalLogProb
  , dirichlet
  , alphas
  , defaultObs
  , Obs(..)
  , gradientScore
  , gradientReparam
  , SampleDouble
  , SampleVector
  , stepTogether
  , stepSeparate
  , unsafeContent
  , GradientParams(..)
  , mergeGeneric
  , mergeGenerics
  , mergeGenericss
  , rhoKuc
  , defaultKucP
  , KucP(..)
  , dist
  , time
  -- , initLocal
  )
where

import           Statistics.BBVI.Class
import           Statistics.BBVI.Propagator
import           Statistics.BBVI.StepSize
import           Statistics.BBVI.Gradient
import           Statistics.BBVI.Observed
import           Statistics.BBVI.Distribution.Normal
import           Statistics.BBVI.Distribution.Dirichlet
import           Statistics.BBVI.Scheduler
