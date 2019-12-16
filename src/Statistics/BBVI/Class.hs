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

class DistUtil a where
  fromParamVector :: ParamVector -> a
  toParamVector :: a -> ParamVector
  nParams :: a -> Int

class DistUtil a => Dist a c where
  resample :: a -> GenST s -> ST s c
  logProb :: a -> c -> Double
  paramGradOfLogQ ::  a -> c -> ParamVector -- gradient of parameters evauated at some sample of x

class (Dist a c) => Differentiable a c where
  gradTransform :: a -> c -> ParamVector
  sampleGradOfLogQ ::  a -> c -> c -- gradient of a sample evaluate with params of q
  transform :: a -> c -> c
  epsilon ::  a -> GenST s -> ST s c
