module Statistics.BBVI
  ( -- * typeclasses for defining new distributions
    DistUtil(..)
  , Dist(..)
  , Differentiable(..)

    -- * types, acessors and merges for distribution cells
  , DistCell(..)
  , DistCells
  , DistCellss
  , DistInvariant(..)
  , dist
  , time
  , defaultDistCell
  , mergeGeneric
  , mergeGenerics
  , mergeGenericss

  -- * step size functions
  , rhoKuc
  , defaultKucP
  , KucP(..)


  -- * normal distribution types and functions
  , NormalDist(..)
  , normalDistr
  , diffableNormalLogProb

  -- * dirichlet distribution types and functions
  , dirichlet
  , alphas
  , defaultObs

  -- * resample observations with a propagator
  , Obs(..)

  -- * helper functions for transforming the log joint into gradient propagators
  , gradientScore
  , gradientReparam

  -- * sampling types
  , SampleDouble
  , SampleVector

  -- * gradient propagator attachment functions
  , stepTogether
  , stepSeparate
  , unsafeContent
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
