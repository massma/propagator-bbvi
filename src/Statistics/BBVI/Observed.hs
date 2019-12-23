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
import           Statistics.BBVI.Propagator     ( DistCell(..)
                                                , SampleVector
                                                , SampleDouble
                                                )
import           System.Random.MWC              ( uniformR )

-- | dummy distribution representing a vector of observations, for use
-- with building "observation" distribution cells.  these can be used
-- to easily/selectively subsample (and resample) observations using
-- existing typeclass methods for use with stochatic gradient updates.
newtype Obs a = O (V.Vector a) deriving (Show, Eq, Ord, Read)

-- | helper function to build a distribution cell of observations
defaultObs :: V.Vector a -> DistCell (Obs a)
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
