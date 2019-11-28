module Lib
    ( someFunc
    ) where
-- idea for tomorrow: likelihood own cell,
-- have updateq watch t, when t advances then
-- grab likelihood and reset it, and update q

import Prelude
import Data.Propagator
import GHC.Generics (Generic)
import Control.Monad (replicateM)
import Control.Monad.ST
import Control.Monad.Primitive
import qualified Data.Vector as V
import System.Random.MWC (create, GenST)
import Statistics.Distribution
import Statistics.Distribution.Normal

newtype Time = T Int deriving (Show, Eq, Ord, Read)

time (T t) = t

type S = (V.Vector Double)
type LogLikelihood = (V.Vector Double)

-- | Q: should we store prior with QProp? - then would have to pass xs to updateQ
type QDistribution = NormalDistribution
data QProp s = Q { gradMemory :: !S
                 , distribution :: !QDistribution
                 , generator :: !(GenST s)
                 , samples :: !(V.Vector Double)
                 }

class VariationalLogic a where
  difference :: a -> a -> Double

  logProb :: a -> Double -> Double

  paramVector :: a -> V.Vector Double

  fromParamVector :: V.Vector Double -> a

  gradLogProb :: a -> Double -> V.Vector Double

  nParams :: a -> Int

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

-- see https://stats.stackexchange.com/questions/404191/what-is-the-log-of-the-pdf-for-a-normal-distribution
instance Propagated (QProp s) where
  merge q1 q2
    -- | f q1 q2 < 0.01 = Change False q1
    | otherwise = Change True q2
    where
      f x1 x2 = difference (distribution x1) (distribution x2)

instance Propagated LogLikelihood where
  merge l1 l2
    | null l2 = Change True l2
    | null l1 = Change True l2
    | otherwise = Change True (V.zipWith (+) l1 l2)

maxStep :: Time
maxStep = T 5

instance Propagated Time where
  merge t1 t2
    | t1 > maxStep = Change False t1
    | t2 > t1 = Change True t2
    | otherwise = Change False t1

-- | don't forget prior with variational distribution - thin kmaybe we
-- should icnorporate prior into QProp and use it when we update q
normalProp prior std xs t q l = do
  watch q $ \qprop ->
    write
      l
      (fmap
         (\mu ->
            logDensity prior mu + -- TODO: move prior updateQ
            (sum $ fmap (logDensity (normalDistr mu std)) xs))
         (samples qprop)) >>
    with t (\(T t') -> write t (T (t' + 1)))

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


updateQ nSamp Q{..} t l =
  Q s' newDist generator <$> V.replicateM nSamp (genContinuous newDist generator)
  where
    alpha = 0.1 -- from kuckelbier et al
    eta = 0.01 -- coservative, from kuckelbiet et al
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

-- how to efficiently watch all of our t's????

-- propagator :: () -- Maybe VariationalProp
propagator xs = runST $ do
  gen <- create
  let prior = normalDistr 0.0 2.0
  q <- known $ Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty
  t <- known $ T 1
  l <- known $ V.empty
  updatePropQ 10 t q l
  normalProp prior 1.0 xs t q l
  q' <- content q
  t' <- content t
  return (distribution <$> q', time <$> t', samples <$> q')

newtype Test = Test Int deriving (Eq, Ord, Read, Show)

instance Propagated Test where
  merge t1 t2
    | t2 >= (Test 30) = Change False t2
    | t1 >= (Test 5) = Change False t2
    | otherwise = Change True t2

testProp = runST $ do
  x <- known $ Test 0
  lift1 (\(Test x) -> Test 100) x x
  content x
-- >>> testProp
-- Just (Test 0)


someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
-- (Just (normalDistr 4.698429664574176e-2 1.9510165133417525),Just 6,Just [-4.175372991091955,-0.6927839199544554,-0.7017981683893668,-1.560056375926866,-0.4274962807947022,-2.006133086056551,2.655967899924271,-0.7725469560239188,0.22084589204472316,-1.5117093518251352])
