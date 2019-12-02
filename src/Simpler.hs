{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Simpler
  ( someFunc
  ) where
import Prelude
import Data.Propagator
import Control.Monad.ST
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import System.Random.MWC (create, uniform, initialize, GenST, uniformR)
import Statistics.Distribution
import Statistics.Distribution.Normal
import qualified System.Random.MWC.Distributions as MWCD
import qualified Statistics.Sample as Samp

class VariationalLogic a where
  logProb :: a -> Double -> Double
  paramVector :: a -> V.Vector Double
  fromParamVector :: V.Vector Double -> a
  gradNuLogQ :: a -> Double -> V.Vector Double -- nu are parameters to variational family
  nParams :: a -> Int
  transform :: a -> Double -> Double
  -- inv transform?
  resample :: a -> GenST s -> ST s Double

class VariationalLogic a => Differentiable a where
  gradNuTransform :: a -> Double -> V.Vector Double
  gradZQ :: a -> Double -> Double

type Time = Int

maxStep :: Time
maxStep = 100000 -- 4664 -- 100000
-- | time can only increase in increments of 1. syntax is to write t
-- to 1, b/c that way an empty t starts at 1

type Memory = V.Vector Double

type Sample = V.Vector Double

type Likelihood = V.Vector Double

type Weight = Double

type Gradient = V.Vector Double

type GradUpdate = (Memory, Gradient)

data PropNode a = N (Node a) | U GradUpdate

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

instance VariationalLogic Obs where
  logProb _d _x = 0.0
  paramVector (O d) = d
  fromParamVector v = (O v)
  gradNuLogQ _d _x = V.empty
  nParams _d = 0
  transform _d eps = eps
  resample (O d) gen =
    return . (d V.!) =<< uniformR (0, (V.length d - 1)) gen

instance Differentiable Obs where
  gradNuTransform _d _x = V.empty
  gradZQ _d _x = 0.0

instance VariationalLogic a => Propagated (PropNode a) where
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

instance VariationalLogic NormalDistribution where
  logProb = logDensity
  paramVector d = V.fromList [mean d, stdDev d]
  fromParamVector v = normalDistr (v V.! 0) ((v V.! 1))
  nParams _d = 2
  gradNuLogQ d x = V.fromList [(x - mu) / std, 1 / std ^ (3 :: Int) * (x - mu) ^ (2 :: Int) - 1 / std]
    where
      mu = mean d
      std = stdDev d
  transform d epsilon = mean d + stdDev d * epsilon
  resample _d gen = MWCD.standard gen

instance Differentiable NormalDistribution where
  gradZQ d z = -(z - mean d)/(stdDev d ** 2)
  gradNuTransform _d epsilon = V.fromList [1.0 , epsilon]

gradient :: VariationalLogic a => Node a -> (Weight, V.Vector Double, Sample) -> PropNode a
gradient q@(Node{..}) (nFactors, like, samples) =
  U (memory', V.zipWith (*) rho' grad)
  where
    summed =
      V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0) $
      V.zipWith
        (\s l -> V.map (* (l + nFactors / weight * (logProb prior s - logProb dist s))) (gradNuLogQ dist s))
        samples
        like
    grad = V.map (/ (fromIntegral $ V.length summed)) summed
    (memory', rho') = rhoF grad q

-- | TODO: speed up by calc length in one pass
gradientAD :: Differentiable a => Node a -> (Weight, V.Vector Double, Sample) -> PropNode a
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
                 (gradZQ prior (transform dist s) -
                  gradZQ dist (transform dist s))))
             (gradNuTransform dist s))
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
     VariationalLogic a
  => (a -> ST s (Weight, V.Vector Double, Sample))
  -> Cell s (PropNode a)
  -> ST s ()
qProp likeFunc q =
  watch q $ \(N q') -> write q =<< (gradient q' <$> likeFunc (dist q'))

qPropAD ::
     Differentiable a
  => (a -> ST s (Weight, V.Vector Double, Sample))
  -> Cell s (PropNode a)
  -> ST s ()
qPropAD likeFunc q =
  watch q $ \(N q') -> write q =<< (gradientAD q' <$> likeFunc (dist q'))

qProp2 ::
     (VariationalLogic a, VariationalLogic b)
  => ((a, b) -> ST s (Weight, V.Vector Double, (Sample, Sample)))
  -> (Cell s (PropNode a), Cell s (PropNode b))
  -> ST s ()
qProp2 likeFunc (q1, q2) =
  watch q1 $ \(N q1') ->
    watch q2 $ \(N q2') -> do
      (w, l, (s1, s2)) <- likeFunc (dist q1', dist q2')
      write q1 (gradient q1' (w, l, s1))
      write q2 (gradient q2' (w, l, s2))

qPropAD2 ::
     (Differentiable a, Differentiable b)
  => ((a, b) -> ST s ((Weight, V.Vector Double, Sample), (Weight, V.Vector Double, Sample)))
  -> (Cell s (PropNode a), Cell s (PropNode b))
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
             V.sum $ V.map (logDensity (normalDistr (transform q eps) std)) xs)
          samples
      , samples)

normalLike2 nSamp nObs std gen (xs, q) = do
  samples <- V.replicateM nSamp (resample q gen)
  obs <- V.replicateM nSamp (resample xs gen)
  return
    ( fromIntegral nObs
    , V.map
        (\eps ->
           V.sum $ V.map (logDensity (normalDistr (transform q eps) std)) obs)
        samples
    , (V.empty, samples))

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
             V.map ((V.! 0) . gradNuLogQ (normalDistr (transform q eps) std)) xs)
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
               ((V.! 0) . gradNuLogQ (normalDistr (transform q eps) std))
               obs)
          samples
      , samples))

propagator xs = runST $ do
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
  q <- known $ N (Node 1 (V.replicate 2 0) (fromIntegral nSamp) qDist prior (rho alpha eta tau epsilon))
  xProp <- known $ N (Node 1 (V.empty) 0 xDist (O V.empty)  (\_x1 _x2 -> (V.empty, V.empty)))
  -- qNormalProp 1.0 gen1 nSamp xs q
  -- qPropAD (normalLikeAD nSamp 1.0 gen1 xs) q
  qPropAD2 (normalLikeAD2 nSamp nSamp 1.0 gen1) (xProp, q) -- (V.length xs)
  (N q') <- fromMaybe (error "impos") <$> content q
  (N xs') <- fromMaybe (error "impos") <$> content xProp
  return (dist q', time q', Samp.mean xs, time xs', V.length . fromO . dist $ xs')

someFunc :: IO ()
someFunc = do
  gen <- create
  xs <- V.replicateM 1000 (genContinuous (normalDistr 5.0 3.0) gen)
  putStrLn (show $ propagator xs)
-- >>> someFunc
--
