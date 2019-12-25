{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
module Statistics.BBVI.Propagator
  ( DistCell(..)
  , DistCells
  , DistCellss
  , SampleVector
  , SampleDouble
  , Gradient
  , Memory
  , Time
  , mergeGeneric
  , mergeGenerics
  , mergeGenericss
  , dist
  , time
  , defaultDistCell
  )
where

import           Statistics.BBVI.Class
import qualified Data.Vector                   as V
import           Data.Propagator

-- | a sample with a vector type e.g. sample from a dirichlet
type SampleVector = V.Vector Double

-- | a sample with a real type e.g. number sampled from a gaussian
type SampleDouble = Double

type Gradient = V.Vector Double

type Memory = V.Vector Double

-- | number of times a distribution cell has been updated
type Time = Int

-- | distribution cell
data DistCell a
  =
    U !Memory !Gradient -- ^ update to a distribution cell
  | Node !Time !Memory !a -- ^ distribution cell with a distribution
                          -- of type a
  deriving (Show, Eq, Ord, Read)

-- !(Gradient -> DistCell a -> (Memory, Gradient))

-- | helper function for initializing distribution cells
defaultDistCell
  :: DistUtil a
  => a -- ^ initial point for distribution
  -> DistCell a -- ^ distribution cell
defaultDistCell d = Node 1 (V.replicate (nParams d) 0) d

-- | time accessor for a distribution cell
time :: DistCell a -> Time
time (Node t _ _) = t
time (U{}       ) = error "called time on update propnode!"
-- memory (Node _ m _) = m

-- | distribution accessor for a distribution cell
dist :: DistCell a -> a
dist (Node _ _ d) = d
dist (U{}       ) = error "called dist on update propnode!"

-- | single cell representing a vector of distribution cells
type DistCells a = V.Vector (DistCell a)

-- | single cell representing an array (vector of vector) of distribution cells
type DistCellss a = V.Vector (V.Vector (DistCell a))

-- | generic merge for DistCells: useful for customizing
-- "quiesence" thresholds with cellWith
mergeGeneric
  :: DistUtil a
  => Time -- ^ maximum number of updates to a cell
  -> Double -- ^ threshold for information to be considered "new"; if
            -- the l2-norm of the change in the distribution's
            -- parameters is less than this threshold, the cell
            -- remains unchanged
  -> DistCell a -- ^ current cell
  -> DistCell a -- ^ proposed update to cell
  -> Change (DistCell a)
mergeGeneric maxStep delta !x1 !x2 = m x1 x2
 where
  m no@(Node t memory d) (U memUp gradUp)
    | norm gradUp < delta = Change False no
    | t >= maxStep        = Change False no
    | otherwise           = Change True updateNode
   where
    updateNode = Node (t + 1)
                      (V.zipWith (+) memory memUp)
                      (fromParamVector newQ)
      where newQ = V.zipWith (+) (toParamVector d) gradUp
  m (U _ _) _ = Contradiction mempty "Trying to update a gradient"
  -- CAREFUL: below is dangerous if I start doing the ideas i thought
  -- about: changing maxstep and elta node for local optmizations
  m no1@(Node t1 _m1 d1) no2@(Node t2 _m2 d2)
    | t1 >= maxStep
    = Change False no1
    | (t2 > t1)
      && (norm (V.zipWith (-) (toParamVector d1) (toParamVector d2)) >= delta)
    = Change True no2
    | otherwise
    = Change False no1

-- | generalization of 'mergeGeneric' to cells of vectors of distributions
mergeGenerics
  :: DistUtil a
  => Time -- ^ maximum number of updates to the cells
  -> Double -- ^ threshold for information to be considered "new"; if
            -- the l2-norm of change in the distribution's parameters
            -- is less than this threshold for all distributions, the
            -- cell remains unchanged
  -> DistCells a -- ^ current cells
  -> DistCells a -- ^ proposed update to cells
  -> Change (DistCells a)
mergeGenerics m d x1 x2 = V.sequence . V.zipWith (mergeGeneric m d) x1 $ x2

-- | generalization of 'mergeGeneric' to cells of arrays of distributions
mergeGenericss
  :: DistUtil a
  => Time -- ^ maximum number of updates to the cells
  -> Double -- ^ threshold for information to be considered "new"; if
            -- the l2-norm of change in the distribution's parameters
            -- is less than this threshold for all distributions, the
            -- cell remains unchanged
  -> DistCellss a -- ^ current cell
  -> DistCellss a -- ^ proposed update to cell
  -> Change (DistCellss a)
mergeGenericss m d v1 v2 = V.sequence . V.zipWith (mergeGenerics m d) v1 $ v2

instance DistUtil a => Propagated (DistCell a) where
  merge = mergeGeneric 1000000 1e-16

instance DistUtil a => Propagated (DistCells a) where
  merge ns updates = V.sequence $ V.zipWith merge ns updates

instance DistUtil a => Propagated (DistCellss a) where
  merge ns updates = V.sequence $ V.zipWith merge ns updates

norm :: V.Vector Double -> Double
norm = sqrt . V.sum . V.map (^ (2 :: Int))
