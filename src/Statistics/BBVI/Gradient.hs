{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}
module Statistics.BBVI.Gradient
  ( gradientScore
  , gradientReparam
  , DistInvariant(..)
  )
where

import qualified Data.Vector                   as V
import           Statistics.BBVI.Propagator     ( DistCell(..)
                                                , Gradient
                                                , Memory
                                                , dist
                                                )
import           Statistics.BBVI.Class

type Samples = V.Vector

-- | invariant data for a variational distribution cell
data DistInvariant a = DistInvariant { weight :: !Double
                                     -- ^ total occurences of z in log likelihood
                                     , prior :: !a
                                     -- ^ prior on latent variable
                                     , rhoF :: !(Gradient -> DistCell a -> (Memory, V.Vector Double))
                                     -- ^ calculates change in memory and setpsize
                                     }

-- | internal helper function for building gradient propagators
gradient
  :: (DistUtil a)
  => (DistCell a -> Double -> c -> Double -> V.Vector Double)
  -> DistInvariant a
  -> DistCell a
  -> (Double, V.Vector Double, Samples c)
  -> DistCell a3
gradient f (DistInvariant {..}) no (nFactors, like, samples) = U
  memory'
  (V.zipWith (*) rho' gr)
 where
  summed =
    V.foldl' (V.zipWith (+)) (V.replicate (nParams (dist no)) 0.0)
      $ V.zipWith (f no nFactors) samples like
  gr              = fmap (/ (fromIntegral $ V.length samples)) summed
  (memory', rho') = rhoF gr no


-- | helper function for transforming log joint to score gradient
-- propagator
gradientScore
  :: Dist a c
  => DistInvariant a
  -> DistCell a -- ^ current variational distribution
  -> (Double, V.Vector Double, Samples c)   -- ^ n factors in log
  -- joint, log joint, samples (corresponding ot each log joint)
  -> DistCell a -- ^ update to variational distribution
gradientScore gp@(DistInvariant {..}) = gradient f gp
 where
  f (Node _time _memory d) nFactors s l = fmap
    (* (l + nFactors / weight * (logProb prior s - logProb d s)))
    (paramGradOfLogQ d s)
  f (U _ _) _ _ _ = error "called gradient update on update cell"


-- | helper function for transforming log joint to reparameterization
-- gradient propagator
gradientReparam
  :: (Differentiable a Double)
  => DistInvariant a
  -> DistCell a -- ^ current variational distribution
  -> (Double, V.Vector Double, Samples Double)
  -- ^ n factors in log joint, gradient of joint wrt z, samples
  -- (corresponding ot each gradient)
  -> DistCell a -- ^ update to variational distribution
gradientReparam gp@(DistInvariant {..}) = gradient f gp
 where
  f (Node _time _memory d) nFactors s g = fmap
    (* ( g
       + nFactors
       / weight
       * ( sampleGradOfLogQ prior (transform d s)
         - sampleGradOfLogQ d     (transform d s)
         )
       )
    )
    (gradTransform d s)
  f (U _ _) _ _ _ = err

err :: a
err = error "called gradient update on update cell"
