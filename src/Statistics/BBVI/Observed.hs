{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Statistics.BBVI.Observed
  ( Obs(..)
  , defaultObs
  )
where

import qualified Data.Vector                   as V
import           Statistics.BBVI.Class
import           Statistics.BBVI.Propagator     ( PropNode(..) )
import           System.Random.MWC              ( uniformR )

newtype Obs a = O (V.Vector a) deriving (Show, Eq, Ord, Read)

defaultObs :: V.Vector Double -> PropNode (Obs Double)
defaultObs xs =
  (Node 1
        100000 -- global max step
        1e-5 -- delta
        (V.empty)
        0
        (O xs)
        (O V.empty)
        (\_x1 _x2 -> (V.empty, V.empty))
  )

instance Dist (Obs Double) Double where
  resample (O d) gen = return . (d V.!) =<< uniformR (0, (V.length d - 1)) gen
  logProb _d _x = 0
  paramGradOfLogQ _d _x = V.empty

instance DistUtil (Obs Double) where
  nParams _x = 0
  toParamVector _ = V.empty
  fromParamVector _ = O V.empty
