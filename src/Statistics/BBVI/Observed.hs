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
import           Statistics.BBVI.Propagator     ( PropNode(..)
                                                , SampleVector
                                                , SampleDouble
                                                )
import           System.Random.MWC              ( uniformR )

newtype Obs a = O (V.Vector a) deriving (Show, Eq, Ord, Read)

defaultObs :: V.Vector a -> PropNode (Obs a)
defaultObs d = (Node 1 V.empty (O d))

instance DistUtil (Obs Double) where
  nParams _x = 0
  toParamVector _ = V.empty
  fromParamVector _ = O V.empty

instance Dist (Obs Double) SampleDouble where
  resample (O d) gen = return . (d V.!) =<< uniformR (0, (V.length d - 1)) gen
  logProb _d _x = 0
  paramGradOfLogQ _d _x = V.empty

instance DistUtil (Obs SampleVector) where
  nParams _x = 0
  toParamVector _ = V.empty
  fromParamVector _ = O V.empty

instance Dist (Obs SampleVector) SampleVector where
  resample (O d) gen = return . (d V.!) =<< uniformR (0, (V.length d - 1)) gen
  logProb _d _x = 0
  paramGradOfLogQ _d _x = V.empty
