{-# LANGUAGE RecordWildCards #-}
module Statistics.BBVI.StepSize
  ( KucP(..)
  , rhoKuc
  , defaultKucP
  )
where

import qualified Data.Vector                   as V
import           Statistics.BBVI.Propagator     ( DistCell(..)
                                                , Gradient
                                                , Memory
                                                )

-- | step size, vector as same length as n parameters in distribution
type Rho = V.Vector Double

-- | calculates a step size following Kucukelbir et al 2017
rhoKuc
  :: KucP -- ^ parameters
  -> Gradient -- ^ gradient at time t
  -> DistCell a -- ^ distribution cell (for memory)
  -> (Memory, Rho) -- ^ change in memory, stepsize
rhoKuc KucP {..} gra (Node time memory _dist) =
  ( deltaM
  , V.zipWith
    (\ds s ->
      eta
        *  (fromIntegral time)
        ** (negate 0.5 + eps)
        *  (1.0 / (tau + sqrt (s + ds)))
    )
    deltaM
    memory
  )
 where
  deltaM = V.zipWith (\g s -> alpha * g ^ (2 :: Int) - alpha * s) gra memory
rhoKuc _kp _gra (U{}) = error "called step size cell on an update cell!"

-- | parameters for 'rhoKuc'
data KucP = KucP
  { alpha :: Double
  , eta :: Double
  , tau :: Double
  , eps :: Double
  } deriving (Show, Eq, Ord, Read)

-- | default parameters for 'rhoKuc.' eta is what you probably want to
-- tune: kucukelbir et al try 0.01 0.1 1 10 100
defaultKucP :: KucP
defaultKucP = KucP { alpha = 0.1, eta = 0.1, tau = 1.0, eps = 1e-16 }
