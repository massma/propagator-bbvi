{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Lib
  ( someFunc
  ) where
import Prelude
import Data.Propagator
import GHC.Generics (Generic)
import Control.Monad (replicateM, when)
import Control.Monad.ST
import Control.Monad.Primitive
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import Data.Word (Word32)
import Numeric.AD (grad', grad, diff, Mode, auto)
import System.Random.MWC (create, GenST, uniform, initialize)
import qualified System.Random.MWC.Distributions as MWCD
import Statistics.Distribution
import Statistics.Distribution.Normal
import qualified Statistics.Sample as Samp

class VariationalLogic a where
  difference :: a -> a -> Double
  logProb :: a -> Double -> Double
  paramVector :: a -> V.Vector Double
  fromParamVector :: V.Vector Double -> a
  gradNuLogQ :: a -> Double -> V.Vector Double -- nu are parameters to variational family
  nParams :: a -> Int

class VariationalLogic a => Differentiable a where
  gradNuTransform :: a -> Double -> V.Vector Double
  gradZQ :: a -> Double -> Double
  transform :: a -> Double -> Double

newtype Time = T Int deriving (Show, Eq, Ord, Read, Num)

time (T t) = t

maxStep :: Time
maxStep = T 100000

-- | time can only increase in increments of 1. syntax is to write t
-- to 1, b/c that way an empty t starts at 1
instance Propagated Time where
  merge t1 t2
    | t1 >= maxStep = Change False t1
    | otherwise = Change True (T (time t1 + 1))

newtype Memory = M (V.Vector Double) deriving (Show, Eq, Ord, Read)

norm = sqrt . V.sum . V.map (^2)

instance Propagated Memory where
  merge (M s1) (M s2) = Change True (M s2) -- always increasing, max at max double for compu

newtype StandardNormal = SN (Time, V.Vector Double) deriving (Show, Eq, Ord, Read)

fromStandard (SN x) = x

instance Propagated (StandardNormal) where
  merge (SN (t1, v1)) (SN (t2, v2))
    | t1 >= maxStep = Change False (SN (t1, v1))
    | t2 > t1 = Change True (SN (t2, v2))
    | otherwise = Change False (SN (t1, v1))

type Sample = V.Vector Double

instance Propagated Sample where
  merge s1 s2
    | (s1 V.! 0) /= (s2 V.! 0) = Change True s2
    | otherwise = Change False s1

type GradLike = V.Vector Double
type Like = V.Vector Double
type Count = Int
type Weight = Double
data LogLikelihood = LL { timeL :: !Time
                        , cnt :: Count
                        , maxCnt :: Count
                        , gradLike :: !(Maybe (Weight, (V.Vector Double)))
                        , like :: !(Maybe (Weight, (V.Vector Double)))
                        }

instance Propagated LogLikelihood where
  merge l1 l2
    | (cnt l1 == 0) && (newCnt == maxCnt l1)  = Change True (LL (timeL l1 + T 1) 0 (maxCnt l1) (gradLike l2) (like l2))
    | (cnt l1 == 0) = Change True (LL (timeL l1) newCnt (maxCnt l1) (gradLike l2) (like l2))
    | cnt l1 == maxCnt l1 = Change True (LL (timeL l1 + T 1) 0 (maxCnt l1) (m gradLike) (m like))
    | timeL l1 >= maxStep = Change False l1
    | otherwise = Change True (LL (timeL l1) newCnt (maxCnt l1) (m gradLike) (m like))
    where
      newCnt = cnt l1 + 1
      m f =
        case (f l1, f l2) of
          (Nothing, Nothing) -> Nothing
          (Just (w1, x1), Just (w2, x2)) -> Just $ (w1 + w2, V.zipWith (+) x1 x2)
          (Just x1, Nothing) -> Just x1
          (Nothing, Just x2) -> Just x2

-- | Q: should we store prior with QProp? - then would have to pass xs to updateQ
instance Propagated NormalDistribution where
  merge q1 q2
    | difference q1 q2 < 0.00001 = Change False q1 -- 0.00001
    | otherwise = Change True q2

instance VariationalLogic NormalDistribution where
  difference x1 x2 = sqrt (f mean + f stdDev)
    where
      f g = (g x1 - g x2) ^ 2
  logProb = logDensity
  paramVector d = V.fromList [mean d, stdDev d]
  fromParamVector v = normalDistr (v V.! 0) (v V.! 1)
  nParams _d = 2
  gradNuLogQ d x = V.fromList [(x - mu) / std, 1 / std ^ 3 * (x - mu) ^ 2 - 1 / std]
    where
      mu = mean d
      std = stdDev d

instance Differentiable NormalDistribution where
  gradZQ d z = -(z - mean d)/(stdDev d ** 2)
  gradNuTransform d epsilon = V.fromList [1.0 , epsilon]
  transform d epsilon = mean d + stdDev d * epsilon

gradient prior total dist Nothing samples = Nothing
gradient prior total dist (Just (cnt, like)) samples = Just $ V.map (/ (fromIntegral $ V.length summed)) summed
  where
    factor = cnt / total
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l -> V.map (* (l + factor * (logDensity prior s - logProb dist s))) (gradNuLogQ dist s))
        samples
        like

gradientAD prior total dist Nothing transformedSamples samples = Nothing
gradientAD prior total dist (Just (cnt, like)) transformedSamples samples = Just $ V.map (/ (fromIntegral $ V.length summed)) summed
  where
    factor = cnt / total
    summed =
      -- TODO: add prior!!
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith3
        (\tS s l -> V.map (* (l + factor * (gradZQ prior tS - gradZQ dist tS))) (gradNuTransform dist s))
        transformedSamples
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

resample q stdNorm derived =
  watch q $ \q' ->
    with stdNorm $ \(SN (_t', s')) ->
      write derived $ V.map (transform q') s'

updateQ prior (T t) (M gradMemory) q s stdNorm l =
  ((M s'), newQ)
  where
    alpha = 0.1 -- from kuckelbier et al
    eta = 0.1 -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
    tau = 1.0
    epsilon = 1e-16
    dft = (V.replicate (nParams q) 0.0)
    gradAD = fromMaybe dft  $ gradientAD prior (fromIntegral (maxCnt l)) q (gradLike l) s stdNorm
    gradS = fromMaybe dft $ gradient prior (fromIntegral (maxCnt l)) q (like l) s
    grad = V.zipWith (+) gradAD gradS
    (s', rho') = rho alpha eta tau epsilon t grad gradMemory
    newQ = fromParamVector $ V.zipWith3 (\x g r -> x + r*g) (paramVector q) grad rho'
-- >>> create >>= \gen -> updateQ 10 (Q (V.fromList [0.0, 0.0]) (normalDistr 0.0 2.0) gen V.empty V.empty) 1

updateQProp gen nSamp prior l stdNorm s q memory =
  watch l $ \l' ->
    with stdNorm $ \(SN (tGlobal, stdNorm')) ->
      when (timeL l' == tGlobal) $
      with q $ \q' ->
        with s $ \s' ->
          with memory $ \mem' ->
            let (memNew, qNew) = updateQ prior tGlobal mem' q' s' stdNorm' l'
             in V.replicateM nSamp (MWCD.standard gen) >>= \newSamp ->
                  write stdNorm (SN (timeL l' + (T 1), newSamp)) >>
                  write memory memNew >>
                  write q qNew

-- | don't forget prior with variational distribution - think maybe we
-- should icnorporate prior into QProp and use it when we update q
-- TODO: consider setting explicit minimum on stddev

normalPropAD std xs s l = do
  watch s $ \_s' ->
    (fromMaybe (error "impos") <$> content  s) >>= \s' ->
      write
        l
        (LL
           (T 1)
           (1 :: Int)
           (1 :: Int)
           (Just
              ( 1.0
              , (fmap
                   (\mu ->
                      (sum $ fmap ((V.! 0) . gradNuLogQ (normalDistr mu std)) xs))
                   s')))
           Nothing)

normalProp std xs s l = do
  watch s $ \_s' ->
    (fromMaybe (error "impos") <$> content s) >>= \s' ->
      write
        l
        (LL
           (T 1)
           (1 :: Int)
           (1 :: Int)
           Nothing
           (Just
              ( 1.0
              , (fmap
                   (\mu -> (sum $ fmap (logDensity (normalDistr mu std)) xs))
                   s'))))

-- propagator :: () -- Maybe VariationalProp
propagator xs = runST $ do
  genG <- create
  seed <- V.replicateM 256 (uniform genG)
  gen1 <- initialize seed
  gen2 <- initialize seed
  let prior = normalDistr 0.0 2.0
  let nSamp = 100
  let qDist = (normalDistr 0.0 2.0)
  initSamp <- V.replicateM nSamp (genContinuous qDist gen1)
  initStandard <- V.replicateM nSamp (MWCD.standard gen1)
  q <- known $ qDist
  stdNorm <- known $ SN (T 1, initSamp)
  l <- cell
  samp <- cell
  mem <- known $ M (V.replicate(nParams  qDist) 0.0)
  resample q stdNorm samp
  normalProp 1.0 xs samp l
  -- normalPropAD 1.0 xs samp l
  updateQProp gen1 nSamp prior l stdNorm samp q mem
  q' <- content q
  t' <- content stdNorm
  mem' <- content mem
  return (q', fst . fromStandard <$> t', mem', Samp.mean xs)

someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- V.replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
--
