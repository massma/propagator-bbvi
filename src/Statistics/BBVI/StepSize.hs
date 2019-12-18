{-# LANGUAGE RecordWildCards #-}
module Statistics.BBVI.StepSize
  ( KucP(..)
  , rhoKuc
  , defaultKucP
  )
where

import qualified Data.Vector                   as V
import           Statistics.BBVI.Propagator     ( PropNode(..)
                                                , Gradient
                                                , Memory
                                                )

type Rho = V.Vector Double

rhoKuc :: KucP -> Gradient -> PropNode a -> (Memory, Rho)
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

data KucP = KucP
  { alpha :: Double
  , eta :: Double
  , tau :: Double
  , eps :: Double
  } deriving (Show, Eq, Ord, Read)

-- | eta is what you probably want to tune: kucukelbir trys 0.01 0.1 1 10 100
defaultKucP :: KucP
defaultKucP = KucP { alpha = 0.1, eta = 0.1, tau = 1.0, eps = 1e-16 }
