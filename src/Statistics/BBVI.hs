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

-- >>> (\(Change True x) -> x) $ mergeGeneric 1e-6 10 (Node 1 (V.replicate 2 0.0) (normalDistr 0.0 1.0)) (U (V.fromList [0.1, 0.9]) (V.fromList [0.2, 0.6]))
-- <interactive>:2314:3-23: warning: [-Wincomplete-uni-patterns]
--     Pattern match(es) are non-exhaustive
--     In a lambda abstraction:
--         Patterns not matched:
--             (Change False _)
--             (Contradiction _ _)
-- Node 2 [0.1,0.9] (ND {mean = 0.2, stdDev = 1.6})
-- >>> (\(Change True x) -> x) $ merge (Node 1 (V.replicate 2 0.0) (normalDistr 0.0 1.0)) (U (V.fromList [0.1, 0.9]) (V.fromList [0.2, 0.6]))
-- <interactive>:2315:3-23: warning: [-Wincomplete-uni-patterns]
--     Pattern match(es) are non-exhaustive
--     In a lambda abstraction:
--         Patterns not matched:
--             (Change False _)
--             (Contradiction _ _)
-- Node 2 [0.1,0.9] (ND {mean = 0.2, stdDev = 1.6})
