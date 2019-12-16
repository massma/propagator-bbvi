{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}
module Statistics.BBVI.Gradient
  ( gradientScore
  , gradientReparam
  )
where

import qualified Data.Vector                   as V
import           Statistics.BBVI.Propagator     ( PropNode(..) )
import           Statistics.BBVI.Class

type Samples = V.Vector

gradientScore
  :: Dist a c
  => PropNode a
  -> (Double, V.Vector Double, Samples c)
  -> PropNode a
gradientScore = gradient f
 where
  f Node {..} nFactors s l = fmap
    (* (l + nFactors / weight * (logProb prior s - logProb dist s)))
    (paramGradOfLogQ dist s)

gradient f q@(Node {..}) (nFactors, like, samples) = U
  memory'
  (V.zipWith (*) rho' gr)
 where
  summed =
    V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0)
      $ V.zipWith (f q nFactors) samples like
  gr              = fmap (/ (fromIntegral $ V.length samples)) summed
  (memory', rho') = rhoF gr q

gradientReparam
  :: (Differentiable a Double)
  => PropNode a
  -> (Double, V.Vector Double, Samples Double)
  -> PropNode a
gradientReparam = gradient f
 where
  f Node {..} nFactors s l = fmap
    (* ( l
       + nFactors
       / weight
       * ( sampleGradOfLogQ prior (transform dist s)
         - sampleGradOfLogQ dist  (transform dist s)
         )
       )
    )
    (gradTransform dist s)
