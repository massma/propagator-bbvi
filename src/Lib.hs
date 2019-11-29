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
import qualified Data.Vector as V
import System.Random.MWC (create, GenST)
import Statistics.Distribution
import Statistics.Distribution.Normal

class VariationalLogic a where
  difference :: a -> a -> Double
  logProb :: a -> Double -> Double
  paramVector :: a -> V.Vector Double
  fromParamVector :: V.Vector Double -> a
  gradLogProb :: a -> Double -> V.Vector Double
  nParams :: a -> Int

newtype Time = T Int deriving (Show, Eq, Ord, Read)

time (T t) = t

maxStep :: Time
maxStep = T 1000000

-- | time can only increase in increments of 1. syntax is to write t
-- to 1, b/c that way an empty t starts at 1
instance Propagated Time where
  merge t1 t2
    | t1 >= maxStep = Change False t1
    | otherwise = Change True (T (time t1 + 1))

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
                 , samples :: !(V.Vector Double)
                 }

instance Propagated (QProp s) where
  merge q1 q2
    | f q1 q2 < 0.0001 = Change False q1
    | otherwise = Change True q2
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

rho alpha eta tau epsilon t grad s0 =
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

updateQ nSamp Q{..} t l =
  Q s' newDist generator <$> V.replicateM nSamp (genContinuous newDist generator)
  where
    alpha = 0.1 -- from kuckelbier et al
    eta = 0.1 -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
    tau = 1.0
    epsilon = 1e-16
    grad = gradient distribution l samples
    (s', rho') = rho alpha eta tau epsilon t grad gradMemory
    newDist = fromParamVector $ V.zipWith3 (\x g r -> x + r*g) (paramVector distribution) grad rho'
-- >>> create >>= \gen -> updateQ 10 (Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty V.empty) 1

updatePropQ nSamp t q l =
  watch t $ \(T t') ->
    with q $ \q' ->
      with l $ \l' ->
        updateQ nSamp q' t' l' >>= \newq -> write l V.empty >> write q newq

-- | don't forget prior with variational distribution - thin kmaybe we
-- should icnorporate prior into QProp and use it when we update q
-- TODO: consider setting explicit minimum on stddev
normalProp prior std xs c q l = do
  watch q $ \qprop ->
    write
      l
      (fmap
         (\mu ->
            logDensity prior mu + -- TODO: move prior updateQ
            (sum $ fmap (logDensity (normalDistr mu std)) xs))
         (samples qprop)) >>
    write c (C (1, 1))


-- how to efficiently watch all of our t's????

-- propagator :: () -- Maybe VariationalProp
propagator xs = runST $ do
  gen <- create
  let prior = normalDistr 0.0 2.0
  let nSamp = 100
  let qDist = (normalDistr 0.0 2.0)
  initSamp <- V.replicateM nSamp (genContinuous qDist gen)
  q <- known $ Q (V.fromList [0.0, 0.0]) qDist gen initSamp
  t <- cell
  c <- cell
  l <- known $ V.empty
  updatePropQ nSamp t q l
  normalProp prior 1.0 xs c q l
  watchCnt  c t
  q' <- content q
  t' <- content t
  return (distribution <$> q', time <$> t')

someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
--
