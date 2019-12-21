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


type SampleVector = V.Vector Double

type SampleDouble = Double

type Gradient = V.Vector Double

type Memory = V.Vector Double

type Time = Int

data DistCell a
  = U !Memory !Gradient
  | Node !Time !Memory !a deriving (Show, Eq, Ord, Read)

-- !(Gradient -> DistCell a -> (Memory, Gradient))

defaultDistCell :: DistUtil a => a -> DistCell a
defaultDistCell d = Node 1 (V.replicate (nParams d) 0) d

time :: DistCell a -> Time
time (Node t _ _) = t
-- time (U{}       ) = error "called time on update propnode!"
-- memory (Node _ m _) = m

dist :: DistCell a -> a
dist (Node _ _ d) = d
-- dist (U{}       ) = error "called dist on update propnode!"

type DistCells a = V.Vector (DistCell a)

type DistCellss a = V.Vector (V.Vector (DistCell a))

mergeGeneric
  :: DistUtil a
  => Int
  -> Double
  -> DistCell a
  -> DistCell a
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
  m no1@(Node t1 _m1 d1) no2@(Node t2 m2 d2)
    | t1 >= maxStep
    = Change False no1
    | (t2 > t1)
      && (norm (V.zipWith (-) (toParamVector d1) (toParamVector d2)) >= delta)
    = Change True no2
    | otherwise
    = Change False no1

mergeGenerics m d x1 x2 = V.sequence . V.zipWith (mergeGeneric m d) x1 $ x2
mergeGenericss m d v1 v2 = V.sequence . V.zipWith (mergeGenerics m d) v1 $ v2

instance DistUtil a => Propagated (DistCell a) where
  merge = mergeGeneric 1000000 1e-16

instance DistUtil a => Propagated (DistCells a) where
  merge ns updates = V.sequence $ V.zipWith merge ns updates

instance DistUtil a => Propagated (DistCellss a) where
  merge ns updates = V.sequence $ V.zipWith merge ns updates

norm :: V.Vector Double -> Double
norm = sqrt . V.sum . V.map (^ (2 :: Int))
