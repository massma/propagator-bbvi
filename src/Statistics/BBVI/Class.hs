{-# LANGUAGE MultiParamTypeClasses #-}
module Statistics.BBVI.Class
  ( DistUtil(..)
  , Dist(..)
  , Differentiable(..)
  )
where

import qualified Data.Vector                   as V
import           System.Random.MWC              ( GenST )
import           Control.Monad.ST               ( ST )

type ParamVector = V.Vector Double

-- | utility functions easing vector arithmetic on gradient updates
class DistUtil a where
  -- | build a distribution from a vector of parameters
  fromParamVector :: ParamVector -> a
  -- | convert a distrbution to a vector of parameters
  toParamVector :: a -> ParamVector
  -- | provide the number of parameters of a distribution
  nParams :: a -> Int

-- | minimum necessary functions to calculate score gradients
-- (Ranganath et al 2014)
class DistUtil a => Dist a c where
  -- | generate a sample from a distribution
  resample :: a -> GenST s -> ST s c
  -- | calculate the log probability (or density) of a sample under a distribution
  logProb :: a -> c -> Double
  -- | calculate the gradient of the log probability of a sample under
  -- a distribution, with respect to the parameters. THe location of
  -- parameters' partial derivatives in the vector should match the
  -- location of the parameters in 'fromParamVector' and 'toParamVector'.
  paramGradOfLogQ ::  a -> c -> ParamVector

-- | minimum necessary functions to calculate reparameterization
-- gradients (Kucukelbir et al 2017)
class (Dist a c) => Differentiable a c where
  -- | generate an un-transformed sample
  epsilon ::  a -> GenST s -> ST s c
  -- | transform a sample generated from 'epsilon' to a sample under the distribution
  transform :: a -> c -> c
  -- | the gradient of 'transform' with respect to the paramters of the distribution
  gradTransform :: a -> c -> ParamVector
  -- | the gradient of the log probability of a sample of \(z\) under
  -- the variational distribution, with respect to \(z\)
  sampleGradOfLogQ ::  a -> c -> c
