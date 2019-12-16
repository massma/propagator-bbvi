module Statistics.BBVI.Defaults
  ( globalMaxStep
  , globalDelta
  , globalEta
  )
where

import           Statistics.BBVI.Propagator

globalMaxStep :: Time
globalMaxStep = 100

globalDelta :: Double
globalDelta = 1e-16 -- 0.00001 --

globalEta :: Double
globalEta = 0.1 -- 0.1 --10.0
