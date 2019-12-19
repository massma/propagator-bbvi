{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}
module Statistics.BBVI.Gradient
  ( gradientScore
  , gradientReparam
  , GradientParams(..)
  )
where

import qualified Data.Vector                   as V
import           Statistics.BBVI.Propagator     ( PropNode(..)
                                                , Gradient
                                                , Memory
                                                , dist
                                                )
import           Statistics.BBVI.Class

type Samples = V.Vector

data GradientParams a = GParams { weight :: !Double
                        , prior :: !a
                        , rhoF :: !(Gradient -> PropNode a -> (Memory, V.Vector Double))}

gradient
  :: (DistUtil a)
  => (PropNode a -> Double -> c -> Double -> V.Vector Double)
  -> GradientParams a
  -> PropNode a
  -> (Double, V.Vector Double, Samples c)
  -> PropNode a3
gradient f (GParams {..}) no (nFactors, like, samples) = U
  memory'
  (V.zipWith (*) rho' gr)
 where
  summed =
    V.foldl' (V.zipWith (+)) (V.replicate (nParams (dist no)) 0.0)
      $ V.zipWith (f no nFactors) samples like
  gr              = fmap (/ (fromIntegral $ V.length samples)) summed
  (memory', rho') = rhoF gr no

gradientScore
  :: Dist a c
  => GradientParams a
  -> PropNode a
  -> (Double, V.Vector Double, Samples c)
  -> PropNode a
gradientScore gp@(GParams {..}) = gradient f gp
 where
  f (Node _time _memory d) nFactors s l = fmap
    (* (l + nFactors / weight * (logProb prior s - logProb d s)))
    (paramGradOfLogQ d s)

gradientReparam
  :: (Differentiable a Double)
  => GradientParams a
  -> PropNode a
  -> (Double, V.Vector Double, Samples Double)
  -- ^ nFactor, likelihoods, samples (corresponding ot each likelihood)
  -> PropNode a
gradientReparam gp@(GParams {..}) = gradient f gp
 where
  f (Node _time _memory d) nFactors s l = fmap
    (* ( l
       + nFactors
       / weight
       * ( sampleGradOfLogQ prior (transform d s)
         - sampleGradOfLogQ d     (transform d s)
         )
       )
    )
    (gradTransform d s)
