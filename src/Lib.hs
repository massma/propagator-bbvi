{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Lib
    ( someFunc
    ) where
-- idea for tomorrow: likelihood own cell,
-- have updateq watch t, when t advances then
-- grab likelihood and reset it, and update q

-- better ideas: 1. counter for n likelihood, t watches counter and
-- updates and then resets counter by subtracting value. have counter
-- just add whatever int it is given for merge 2. also test
-- paraleleization by just runnin ginference on two idential gaussians
-- (can be same sample and stuff)
import Prelude
import Data.Propagator
import GHC.Generics (Generic)
import Control.Monad (replicateM, when)
import Control.Monad.ST
import Control.Monad.Primitive
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import Data.Word (Word32)
import System.Random.MWC (create, GenST, uniform, initialize)
import Statistics.Distribution
import Statistics.Distribution.Normal

class VariationalLogic a where
  difference :: a -> a -> Double
  logProb :: a -> Double -> Double
  paramVector :: a -> V.Vector Double
  fromParamVector :: V.Vector Double -> a
  gradLogProb :: a -> Double -> V.Vector Double
  nParams :: a -> Int

newtype Time = T Int deriving (Show, Eq, Ord, Read, Num)

maxStep :: Time
maxStep = T 100000

-- | time can only increase in increments of 1. syntax is to write t
-- to 1, b/c that way an empty t starts at 1
instance Propagated Time where
  merge t1 t2
    | t1 >= maxStep = Change False t1
    | otherwise = Change True (t1 + (T 1))

newtype Counter = C (Int, Int) deriving (Show, Eq, Ord, Read)

counter (C c) = c

instance Propagated Counter where
  merge (C (cnt, maxCt))  c2
    | newCnt > maxCt = Change True (C (1, maxCt))
    | otherwise = Change True (C (newCnt, maxCt))
    where
      newCnt = cnt + 1

type S = (V.Vector Double)
type LogLikelihood = (V.Vector Double)

instance Propagated LogLikelihood where
  merge l1 l2
    | null l2 = Change True l2
    | null l1 = Change True l2
    | otherwise = Change True (V.zipWith (+) l1 l2)


-- | Q: should we store prior with QProp? - then would have to pass xs to updateQ
type QDistribution = NormalDistribution

data QProp s = Q { gradMemory :: !S
                 , distribution :: !QDistribution
                 , generator :: !(GenST s)
                 , time :: !Time
                 , samples :: !(V.Vector Double)
                 }

instance Propagated (QProp s) where
  merge q1 q2
    -- | f q1 q2 < 0.00001 = Change False q1
    | time q2 > maxStep = Change False q1
    | time q2 > time q1 = Change True q2
    | otherwise = Contradiction mempty ""
    where
      f x1 x2 = difference (distribution x1) (distribution x2)

instance VariationalLogic NormalDistribution where
  difference x1 x2 = sqrt (f mean + f stdDev)
    where
      f g = (g x1 - g x2) ^ 2
  logProb = logDensity
  paramVector d = V.fromList [mean d, stdDev d]
  fromParamVector v = normalDistr (v V.! 0) (v V.! 1)
  nParams _d = 2
  gradLogProb d x = V.fromList [(x - mu) / std, 1 / std ^ 3 * (x - mu) ^ 2 - 1 / std]
    where
      mu = mean d
      std = stdDev d

gradient dist like samples = V.map (/ (fromIntegral $ V.length summed)) summed
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l -> V.map (* (l - logProb dist s)) (gradLogProb dist s))
        samples
        like
-- >>> gradient (normalDistr 0.0 2.0) V.empty V.empty
-- [0.0,0.0]

rho alpha eta tau epsilon (T t) grad s0 =
  ( s1
  , V.map
      (\s ->
         eta * (fromIntegral t) ** (negate 0.5 + epsilon) *
         (1.0 / (tau + sqrt s)))
      s1)
  where
    s1 = V.zipWith (\g s -> alpha * g ^ 2 + (1.0 - alpha) * s) grad s0
-- >>> rho 0.1 0.01 1.0 1e-16 1 (V.fromList [0.0, 0.0]) (V.fromList [0.0, 0.0])
-- ([0.0,0.0],[1.0e-2,1.0e-2])

watchCnt c t =
  watch c $ \(C (c, m)) -> when (c == m) (write t (T 1))

updateQ nSamp Q{..} l =
  Q s' newDist generator (time + (T 1)) <$> V.replicateM nSamp (genContinuous newDist generator)
  where
    alpha = 0.1 -- from kuckelbier et al
    eta = 0.1 -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
    tau = 1.0
    epsilon = 1e-16
    grad = gradient distribution l samples
    (s', rho') = rho alpha eta tau epsilon time grad gradMemory
    newDist = fromParamVector $ V.zipWith3 (\x g r -> x + r*g) (paramVector distribution) grad rho'
-- >>> create >>= \gen -> updateQ 10 (Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty V.empty) 1

-- | still not sure this is proper
updatePropQ nSamp q l =
  with q $ \q' -> do
    l' <- fromMaybe (error "no content in likelohood") <$> content l
    newQ <- updateQ nSamp q' l'
    write l V.empty
    write q newQ


-- | don't forget prior with variational distribution - thin maybe we
-- should icnorporate prior into QProp and use it when we update q
-- TODO: consider setting explicit minimum on stddev
normalProp nSamp prior std xs q l = do
  watch q $ \qprop ->
    mapM_
      (\x -> write l (V.map (\mu -> logDensity (normalDistr mu std) x) (samples qprop)))
      xs >>
    updatePropQ nSamp q l
    -- write
    --   l
    --   (fmap
    --      (\mu ->
    --         logDensity prior mu + -- TODO: move prior updateQ
    --         (sum $ fmap (logDensity (normalDistr mu std)) xs))
    --      (samples qprop))

-- how to efficiently watch all of our t's????

-- propagator :: () -- Maybe VariationalProp
propagator xs = runST $ do
  genG <- create
  seed <- V.replicateM 256 (uniform genG)
  gen1 <- initialize seed
  let prior = normalDistr 0.0 2.0
  let nSamp = 100
  let qDist = (normalDistr 0.0 2.0)
  initSamp <- V.replicateM nSamp (genContinuous qDist gen1)
  q1 <- known $ Q (V.fromList [0.0, 0.0]) qDist gen1 (T 1) initSamp
  l1 <- known $ V.empty
  t1 <- known $ (T 1)
  normalProp nSamp prior 1.0 xs q1 l1
  q1' <- content q1
  return (distribution <$> q1', time <$> q1')

someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- (Just (normalDistr 4.669105916353021 0.1600423436219402),Just 2259)
-- >>> someFunc
-- (Just (normalDistr 4.669105916353021 0.1600423436219402),Just 2259)
