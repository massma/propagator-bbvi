{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Simpler
  ( someFunc
  ) where
import Prelude
import Data.Propagator
import Control.Monad.ST
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import Numeric.SpecFunctions (logGamma, digamma)
import System.Random.MWC (create, uniform, initialize, GenST, uniformR)
import Statistics.Distribution
import Statistics.Distribution.Normal
import qualified System.Random.MWC.Distributions as MWCD
import qualified Statistics.Sample as Samp

-- from HMM:
-- splash belief propagation (BP)
type SampleVector = V.Vector Double

type SampleDouble = Double

-- data Sample = V SampleVector | D SampleDouble

-- fromRD (D x) = x

class Dist a b where
  -- | transform <$> resample = true sample from dist
  transform :: a -> b -> b
  -- inv transform?
  resample :: a -> GenST s -> ST s b

class VecParam a where
  paramVector :: a -> V.Vector Double
  fromParamVector :: V.Vector Double -> a
  nParams :: a -> Int

class (VecParam a, Dist a b) => VariLogic a b where
  logProb :: a -> b -> Double
  paramGradOfLogQ :: a -> b -> V.Vector Double -- gradient of parameters evauated at some sample of x

class VariLogic a b => Differentiable a b where
  gradTransform :: a -> b -> V.Vector Double
  sampleGradOfLogQ :: a -> b -> Double -- gradient of a sample evaluate with params of q

type Time = Int

maxStep :: Time
maxStep = 100000 -- 4664 -- 100000

type Memory = V.Vector Double

type Samples = V.Vector

type Likelihood = V.Vector Double

type Weight = Double

type Gradient = V.Vector Double

type GradUpdate = (Memory, Gradient)

data PropNode a b = N (Node a) | U GradUpdate

data Node a = Node
  { time :: !Time
  , memory :: !Memory
  , weight :: !Weight
  , dist :: !a
  , prior :: !a
  , rhoF :: !(V.Vector Double -> Node a -> (Memory, V.Vector Double))
  }

norm :: V.Vector Double -> Double
norm = sqrt . V.sum . V.map (^ (2 :: Int))

newtype Obs = O (V.Vector Double) deriving (Show, Eq, Ord, Read)

fromO (O v) = v

instance Dist Obs Double where
  transform _d eps = eps
  resample (O d) gen =
    return . (d V.!) =<< uniformR (0, (V.length d - 1)) gen

instance VecParam Obs where
  paramVector (O d) = d
  fromParamVector v = (O v)
  nParams _d = 0

instance VariLogic Obs Double where
  logProb _d _x = 0.0
  paramGradOfLogQ _d _x = V.empty

instance Differentiable Obs Double where
  gradTransform _d _x = V.empty
  sampleGradOfLogQ _d _x = 0.0

newtype Dirichlet = Diri (V.Vector Double) deriving (Show, Eq, Ord, Read)

logB :: Dirichlet -> Double
logB (Diri alphas) = V.sum (V.map logGamma alphas) - logGamma (V.sum alphas)

instance Dist Dirichlet SampleVector where
  transform _d x = x
  resample (Diri diri) gen = MWCD.dirichlet diri gen

instance VecParam Dirichlet where
  paramVector (Diri diri) = diri
  fromParamVector = Diri
  nParams (Diri diri) = V.length diri

instance VariLogic Dirichlet SampleVector where
  logProb (Diri diri) cat =
    V.sum (V.zipWith (\alpha x -> (alpha - 1) * log x) diri cat) -
    logB (Diri diri)
  paramGradOfLogQ (Diri diri) cat =
    V.zipWith (\a x -> summed - digamma a + a * log x) diri cat
    where
      summed = digamma (V.sum cat)

instance (VecParam a, VariLogic a b) => Propagated (PropNode a b) where
  merge (N node) (U (deltaM, g))
    | norm g < 0.00001 = Change False (N node) -- 0.00001
    | time node >= maxStep = Change False (N node)
    | otherwise = Change True updateNode
    where
      updateNode =
        N
          (node
             { time = (time node) + 1
             , memory = V.zipWith (+) (memory node) deltaM
             , dist = newQ
             })
        where
          newQ = fromParamVector $ V.zipWith (+) (paramVector (dist node)) g
  merge (U _) _ = Contradiction mempty "Trying to update a gradient"
  merge (N _) (N _) = Contradiction mempty "Trying overwrite a node"

instance VecParam NormalDistribution where
  paramVector d = V.fromList [mean d, stdDev d]
  fromParamVector v = normalDistr (v V.! 0) ((v V.! 1))
  nParams _d = 2

instance Dist NormalDistribution SampleDouble where
  transform d epsilon = mean d + stdDev d * epsilon
  resample _d gen = MWCD.standard gen

instance VariLogic NormalDistribution SampleDouble where
  logProb d x= logDensity d x
  paramGradOfLogQ d x = V.fromList [(x - mu) / std, 1 / std ^ (3 :: Int) * (x - mu) ^ (2 :: Int) - 1 / std]
    where
      mu = mean d
      std = stdDev d

instance Differentiable NormalDistribution SampleDouble where
  sampleGradOfLogQ d z = -(z - mean d)/(stdDev d ** 2)
  gradTransform _d epsilon = V.fromList [1.0 , epsilon]

gradient ::
     (VecParam a, VariLogic a b)
  => Node a
  -> (Weight, V.Vector Double, Samples b)
  -> PropNode a b
gradient q@(Node{..}) (nFactors, like, samples) =
  U (memory', V.zipWith (*) rho' grad)
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l -> V.map (* (l + nFactors / weight * (logProb prior s - logProb dist s))) (paramGradOfLogQ dist s))
        samples
        like
    grad = V.map (/ (fromIntegral $ V.length summed)) summed
    (memory', rho') = rhoF grad q

-- | TODO: speed up by calc length in one pass
gradientAD ::
     Differentiable a b
  => Node a
  -> (Weight, V.Vector Double, Samples b)
  -> PropNode a b
gradientAD q@(Node {..}) (nFactors, like, samples) =
  U (deltaMem, V.zipWith (*) rho' grad)
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l ->
           V.map
             (* (l +
                 nFactors / weight *
                 (sampleGradOfLogQ prior (transform dist s) -
                  sampleGradOfLogQ dist (transform dist s))))
             (gradTransform dist s))
        samples
        like
    grad = V.map (/ (fromIntegral $ V.length summed)) summed
    (deltaMem, rho') = rhoF grad q

rho alpha eta tau epsilon grad Node{..} =
  ( deltaM
  , V.zipWith
      (\ds s ->
         eta * (fromIntegral time) ** (negate 0.5 + epsilon) *
         (1.0 / (tau + sqrt (s + ds))))
      deltaM memory)
  where
    deltaM = V.zipWith (\g s -> alpha * g ^ (2 :: Int) - alpha * s) grad memory

-- | TODO: try making obs a Node and see if it still performant
-- implement weight on each propgator, and then divide by each
-- variationa distributions' factor (I was doing this wrong).
qProp ::
     (VecParam a, VariLogic a b)
  => (a -> ST s (Weight, V.Vector Double, Samples b))
  -> Cell s (PropNode a b)
  -> ST s ()
qProp likeFunc q =
  watch q $ \(N q') -> write q =<< (gradient q' <$> likeFunc (dist q'))

qPropAD ::
     (Differentiable a b)
  => (a -> ST s (Weight, V.Vector Double, Samples b))
  -> Cell s (PropNode a b)
  -> ST s ()
qPropAD likeFunc q =
  watch q $ \(N q') -> write q =<< (gradientAD q' <$> likeFunc (dist q'))

qProp2 ::
     (VariLogic a c, VariLogic b d)
  => ((a, b) -> ST s (Weight, V.Vector Double, (Samples c, Samples d)))
  -> (Cell s (PropNode a c), Cell s (PropNode b d))
  -> ST s ()
qProp2 likeFunc (q1, q2) =
  watch q1 $ \(N q1') ->
    watch q2 $ \(N q2') -> do
      (w, l, (s1, s2)) <- likeFunc (dist q1', dist q2')
      write q1 (gradient q1' (w, l, s1))
      write q2 (gradient q2' (w, l, s2))

qPropAD2 ::
     (Differentiable a c, Differentiable b d)
  => ((a, b) -> ST s ((Weight, V.Vector Double, Samples c), (Weight, V.Vector Double, Samples d)))
  -> (Cell s (PropNode a c), Cell s (PropNode b d))
  -> ST s ()
qPropAD2 likeFunc (q1, q2) =
  watch q1 $ \(N q1') ->
    watch q2 $ \(N q2') -> do
      (l1, l2) <- likeFunc (dist q1', dist q2')
      write q1 (gradientAD q1' l1)
      write q2 (gradientAD q2' l2)

normalLike nSamp std gen xs q =
  V.replicateM nSamp (resample q gen) >>= \samples ->
    return
      ( fromIntegral $ V.length xs
      , V.map
          (\eps ->
             V.sum $ V.map (logProb (normalDistr (transform q eps) std)) xs)
          samples
      , samples)

-- normalLike2 nSamp nObs std gen (xs, q) = do
--   samples <- V.replicateM nSamp (resample q gen)
--   obs <- V.replicateM nSamp (resample xs gen)
--   return
--     ( fromIntegral nObs
--     , V.map
--         (\eps ->
--            V.sum $ V.map (logProb (normalDistr (transform q eps) std) . fromRD) obs)
--         samples
--     , (V.empty, samples))

-- | specialize all fmaps to vector???
-- CAREFUL: with rao-blackweixation and my wieghting approach, all terms
-- in log likelihood functions must contain all nodes q
normalLikeAD nSamp std gen xs q =
  V.replicateM nSamp (resample q gen) >>= \samples ->
    return
      ( fromIntegral $ V.length xs
      , V.map
          (\eps ->
             V.sum $
             V.map ((V.! 0) . paramGradOfLogQ (normalDistr (transform q eps) std)) xs)
          samples
      , samples)

normalLikeAD2 nSamp nObs std gen (xs, q) = do
  obs <- V.replicateM nObs (resample xs gen)
  samples <- V.replicateM nSamp (resample q gen)
  return
    ( ( (fromIntegral nObs)
      , V.empty
      , V.empty)
    , ( (fromIntegral nObs)
      , V.map
          (\eps ->
             V.sum $
             V.map
               ((V.! 0) . paramGradOfLogQ (normalDistr (transform q eps) std))
               (obs :: V.Vector Double))
          samples
      , samples))

normalFit xs =
  runST $ do
    genG <- create
    gen1 <- initialize =<< V.replicateM 256 (uniform genG)
    let prior = normalDistr 0.0 2.0
    let nSamp = 100
    let qDist = normalDistr 0.0 2.0
    let alpha = 0.1 -- from kuckelbier et al
    let eta = 0.1 -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
    let tau = 1.0
    let epsilon = 1e-16 -- (fromIntegral $ V.length xs)
    let xDist = (O xs)
    q <-
      known $
      (N (Node
            1
            (V.replicate 2 0)
            (fromIntegral nSamp)
            qDist
            prior
            (rho alpha eta tau epsilon)) :: PropNode NormalDistribution Double)
    xProp <-
      known $
      (N (Node 1 (V.empty) 0 xDist (O V.empty) (\_x1 _x2 -> (V.empty, V.empty))) :: PropNode Obs Double)
  -- qNormalProp 1.0 gen1 nSamp xs q
  -- qPropAD (normalLikeAD nSamp 1.0 gen1 xs) q
    qPropAD2 (normalLikeAD2 nSamp nSamp 1.0 gen1) (xProp, q) -- (V.length xs)
    (N q') <- fromMaybe (error "impos") <$> content q
    (N xs') <- fromMaybe (error "impos") <$> content xProp
    return (dist q', time q', Samp.mean xs, time xs')

genNormal = do
  gen <- create
  xs <- V.replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  return xs

someFunc :: IO ()
someFunc = do
  xs <- genNormal
  putStrLn (show $ normalFit xs)
-- >>> someFunc
--
