{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
module Statistics.BBVI.Propagator
  ( PropNode(..)
  , PropNodes(..)
  , PropNodess(..)
  , SampleVector
  , SampleDouble
  , Time
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

data PropNode a
  = U { memoryUpdate :: !Memory
      , gradientUpdate :: !Gradient }
  | Node { time :: !Time
         , maxStep :: !Time
         , delta :: !Double
         , memory :: !(V.Vector Double)
         , weight :: !Double
         , dist :: !a
         , prior :: !a
         , rhoF :: !(Gradient -> PropNode a -> (Memory, Gradient)) }

type PropNodes a = V.Vector (PropNode a)

type PropNodess a = V.Vector (V.Vector (PropNode a))

instance DistUtil a => Propagated (PropNode a) where
  merge node@(Node {..}) (U {..})
    | norm gradientUpdate < delta = Change False node
    | time >= maxStep             = Change False node
    | otherwise                   = Change True updateNode
   where
    updateNode =
      (node { time   = time + 1
            , memory = V.zipWith (+) memory memoryUpdate
            , dist   = fromParamVector newQ
            }
      )
      where newQ = V.zipWith (+) (toParamVector dist) gradientUpdate
  merge (U _ _) _ = Contradiction mempty "Trying to update a gradient"
-- CAREFUL: below is dangerous if I start doing the ideas i thought
-- about: changing maxstep and elta node for local optmizations
  merge node1@(Node{}) node2@(Node{})
    | time node1 >= maxStep node1
    = Change False node1
    | (time node2 > time node1)
      && (  norm
             (V.zipWith (-)
                        (toParamVector $ dist node1)
                        (toParamVector $ dist node2)
             )
         >= (delta node1)
         )
    = Change True (node2 { maxStep = maxStep node1, delta = delta node1 })
    | otherwise
    = Change False node1

instance DistUtil a => Propagated (PropNodes a) where
  merge nodes updates = V.sequence $ V.zipWith merge nodes updates

instance DistUtil a => Propagated (PropNodess a) where
  merge nodes updates = V.sequence $ V.zipWith merge nodes updates

norm = sqrt . V.sum . V.map (^ (2 :: Int))
