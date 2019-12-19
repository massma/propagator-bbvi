module Statistics.BBVI.Scheduler
  ( stepTogether
  , unsafeContent
  )
where

import           Data.Propagator
import           Data.Maybe
import           Control.Monad.ST               ( ST )
import           Statistics.BBVI.Propagator     ( Time
                                                , PropNode(..)
                                                )

unsafeContent :: Cell s a -> ST s a
unsafeContent =
  (fromMaybe (error "called unsafe content but no content in cell") <$>)
    . content

-- | todo: make below a cellWith thing
-- initLocal :: Time -> Node a -> Node a
-- initLocal step p = p { maxStep = time p + step }

stepTogether :: (a -> b -> ST s (a, b)) -> Cell s a -> Cell s b -> ST s ()
{-# INLINE stepTogether #-}
stepTogether f x y = watch x $ \x' -> with y $ \y' -> do
  (upX, upY) <- f x' y'
  write y upY
  write x upX
