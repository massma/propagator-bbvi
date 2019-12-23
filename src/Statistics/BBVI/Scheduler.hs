module Statistics.BBVI.Scheduler
  ( stepTogether
  , stepSeparate
  , unsafeContent
  )
where

import           Control.Monad.ST               ( ST )
import           Data.Propagator
import           Data.Maybe                     ( fromMaybe )
import qualified Data.Vector                   as V
import           Statistics.BBVI.Class
import           Statistics.BBVI.Propagator     ( Time
                                                , time
                                                , DistCell(..)
                                                , DistCellss
                                                , mergeGeneric
                                                , mergeGenericss
                                                )
-- | partial function to pull content from cell
unsafeContent :: Cell s a -> ST s a
unsafeContent =
  (fromMaybe (error "called unsafe content but no content in cell") <$>)
    . content

-- | update all propagators in one step until first argument reaches
-- quiesence; synonymous with traditional synchronous algorithms
stepTogether :: (a -> b -> ST s (a, b)) -> Cell s a -> Cell s b -> ST s ()
stepTogether f x y = watch x $ \x' -> with y $ \y' -> do
  (upX, upY) <- f x' y'
  -- careful: order below can create a bug. safer would be to use
  -- unsafe content instead of with y
  write y upY
  write x upX

-- | partially optimize each cell with the other held fixed, and
-- repeat until quiesence. Each local optimization will proceed some
-- number of steps, or until the local optimization converges, before
-- triggering the global cell to fire.
stepSeparate
  :: (DistUtil a1, DistUtil a2)
  => Time -- ^ number of local steps
  -> Double
  -> (  DistCell a2
     -> DistCellss a1
     -> ST s (DistCell a2, DistCellss a1)
     ) -- ^ gradient propagator
  -> Cell s (DistCell a2)
  -> Cell s (V.Vector (V.Vector (DistCell a1))) --
  -> ST s ()
stepSeparate nLocal localDelta f x0 ys0 = do
  watch x0 $ \x -> do
    yTemp <- unsafeContent ys0
    ys    <- cellWith $ mergeGenericss
      ((+ nLocal) . time . (V.! 0) . (V.! 0) $ yTemp)
      localDelta
    write ys yTemp
    watch ys $ \ys' -> do
      (_upX, upYs) <- f x ys'
      -- upYs <- f2 x ys'
      write ys upYs
    ysNew <- unsafeContent ys
    write ys0 ysNew
  watch ys0 $ \ys -> do
    xTemp <- unsafeContent x0
    x     <- cellWith $ mergeGeneric (nLocal + time xTemp) localDelta
    write x xTemp
    watch x $ \x' -> do
      (upX, _upYs) <- f x' ys
      -- upX <- f1 x' ys
      write x upX
    xNew <- unsafeContent x
    write x0 xNew
