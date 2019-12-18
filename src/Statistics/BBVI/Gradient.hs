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
                                                )
import           Statistics.BBVI.Class

type Samples = V.Vector

data GradientParams a = GParams { weight :: !Double
                        , prior :: !a
                        , rhoF :: !(Gradient -> PropNode a -> (Memory, V.Vector Double))}

-- weight (Node _ _ _ _ w _ _ _) = w
-- prior (Node _ _ _ _ _ _ p _) = p
-- rhoF (Node _ _ _ _ _ _ _ r) = r


gradient
  :: (DistUtil a)
  => (PropNode a -> Double -> c -> Double -> V.Vector Double)
  -> GradientParams a
  -> PropNode a
  -> (Double, V.Vector Double, Samples c)
  -> PropNode a3
gradient f (GParams {..}) no@(Node _time _memory dist) (nFactors, like, samples)
  = U memory' (V.zipWith (*) rho' gr)
 where
  summed =
    V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0)
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
  f (Node _time _memory dist) nFactors s l = fmap
    (* (l + nFactors / weight * (logProb prior s - logProb dist s)))
    (paramGradOfLogQ dist s)

gradientReparam
  :: (Differentiable a Double)
  => GradientParams a
  -> PropNode a
  -> (Double, V.Vector Double, Samples Double)
  -- ^ nFactor, likelihoods, samples (corresponding ot each likelihood)
  -> PropNode a
gradientReparam gp@(GParams {..}) = gradient f gp
 where
  f (Node _time _memory dist) nFactors s l = fmap
    (* ( l
       + nFactors
       / weight
       * ( sampleGradOfLogQ prior (transform dist s)
         - sampleGradOfLogQ dist  (transform dist s)
         )
       )
    )
    (gradTransform dist s)
