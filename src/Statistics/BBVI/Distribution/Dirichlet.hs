{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
module Statistics.BBVI.Distribution.Dirichlet
  ( alphas
  , defaultDirichlet
  , dirichlet
  , Dirichlet(..)
  )
where

import qualified Data.Vector                   as V
import           Numeric.SpecFunctions          ( logGamma
                                                , digamma
                                                )
import           Statistics.BBVI.Class
import           Statistics.BBVI.StepSize
import           Statistics.BBVI.Propagator     ( PropNode(..)
                                                , SampleVector
                                                )
import qualified System.Random.MWC.Distributions
                                               as MWCD


newtype Dirichlet = Diri (V.Vector Double) deriving (Show, Eq, Ord, Read)

defaultDirichlet prior =
  (Node 1
        100000
        1e-6
        (fmap (const 0.0) (toParamVector prior))
        1
        prior
        prior
        (rhoKuc defaultKucP)
  )

dirichlet xs = Diri $ V.map (max 1e-100) xs

alphas :: Dirichlet -> V.Vector Double
alphas (Diri xs) = xs

logB :: Dirichlet -> Double
logB d = V.sum (V.map logGamma as) - logGamma (V.sum as) where as = alphas d

instance DistUtil Dirichlet where
  nParams         = V.length . alphas
  fromParamVector = dirichlet
  toParamVector   = alphas

instance Dist Dirichlet SampleVector where
  resample d gen = MWCD.dirichlet (alphas d) gen
  logProb d cat = V.sum (V.zipWith (\alpha x -> (alpha - 1) * log x) as cat)
    - logB d
    where as = alphas d
  paramGradOfLogQ d cat = V.zipWith (\a x -> summed - digamma a + log x) as cat
   where
    as     = alphas d
    summed = digamma (V.sum as)
