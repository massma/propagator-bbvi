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
import Numeric.AD (grad, auto)
import Numeric.MathFunctions.Constants (m_sqrt_2_pi)
import Numeric.SpecFunctions (logGamma, digamma)
import System.Random.MWC (create, uniform, initialize, GenST, uniformR)
-- import Statistics.Distribution
-- import Statistics.Distribution.Normal
import qualified System.Random.MWC.Distributions as MWCD
import qualified Statistics.Sample as Samp

-- from HMM:
-- splash belief propagation (BP)
type SampleVector = V.Vector Double

type SampleDouble = Double

-- data Sample = V SampleVector | D SampleDouble

-- fromRD (D x) = x

class Functor a => Dist a b where
  -- | transform <$> resample = true sample from dist
  transform :: a Double -> b -> b
  -- inv transform?
  resample :: a Double -> GenST s -> ST s b

class (Functor a) => VariUtil a where
  zipDist :: (Double -> Double -> Double) -> a Double -> a Double -> a Double
  norm :: a Double -> Double

class (VariUtil a, Dist a b) => VariLogic a b where
  logProb :: a Double -> b -> Double
  paramGradOfLogQ :: a Double -> b -> a Double -- gradient of parameters evauated at some sample of x

class VariLogic a b => Differentiable a b where
  gradTransform :: a Double -> b -> a Double
  sampleGradOfLogQ :: a Double -> b -> Double -- gradient of a sample evaluate with params of q

instance VariUtil V.Vector where
  zipDist = V.zipWith
  norm = sqrt . V.sum . V.map (^ (2 :: Int))

type Time = Int

maxStep :: Time
maxStep = 100000 -- 4664 -- 100000

type Memory = V.Vector Double

type Samples = V.Vector

type Likelihood = V.Vector Double

type Weight = Double

type Gradient = V.Vector Double

data PropNode a b = N (Node a) | U (a Double, a Double) -- memory, gradient

data Node a = Node
  { time :: !Time
  , memory :: !(a Double)
  , weight :: !Weight
  , dist :: !(a Double)
  , prior :: !(a Double)
  , rhoF :: !(a Double -> Node a -> (a Double, a Double)) -- memory, gradient
  }

newtype Obs a = O (V.Vector a) deriving (Show, Eq, Ord, Read)

fromO (O v) = v

instance Dist Obs Double where
  transform _d eps = eps
  resample (O d) gen =
    return . (d V.!) =<< uniformR (0, (V.length d - 1)) gen

instance Functor Obs where
  fmap _f _o = O V.empty -- WARNING: does not satisfy functor laws

instance VariUtil Obs where
  zipDist _f _o1 _o2 = O V.empty
  norm _ = 0.0

instance VariLogic Obs Double where
  logProb _d _x = 0.0
  paramGradOfLogQ _d _x = O V.empty

instance Differentiable Obs Double where
  gradTransform _d _x = O V.empty
  sampleGradOfLogQ _d _x = 0.0

newtype Dirichlet a = Diri (V.Vector a) deriving (Show, Eq, Ord, Read)

logB :: Dirichlet Double -> Double
logB (Diri alphas) = V.sum (V.map logGamma alphas) - logGamma (V.sum alphas)

instance Dist Dirichlet SampleVector where
  transform _d x = x
  resample (Diri diri) gen = MWCD.dirichlet diri gen

instance Functor Dirichlet where
  fmap f (Diri diri) = Diri $ V.map f diri

instance VariUtil Dirichlet where
  zipDist f (Diri x1) (Diri x2) = Diri $ V.zipWith f x1 x2
  norm (Diri x1) = norm x1

instance VariLogic Dirichlet SampleVector where
  logProb (Diri diri) cat =
    V.sum (V.zipWith (\alpha x -> (alpha - 1) * log x) diri cat) -
    logB (Diri diri)
  paramGradOfLogQ (Diri diri) cat =
    Diri $ V.zipWith (\a x -> summed - digamma a + a * log x) diri cat
    where
      summed = digamma (V.sum cat)

instance VariLogic a b => Propagated (PropNode (a) b) where
  merge (N node) (U (deltaM, g))
    | norm g < 0.00001 = Change False (N node) -- 0.00001
    | time node >= maxStep = Change False (N node)
    | otherwise = Change True updateNode
    where
      updateNode =
        N
          (node
             { time = (time node) + 1
             , memory = zipDist (+) (memory node) deltaM
             , dist = newQ
             })
        where
          newQ = zipDist (+) (dist node) g
  merge (U _) _ = Contradiction mempty "Trying to update a gradient"
  merge (N _) (N _) = Contradiction mempty "Trying overwrite a node"

newtype NormalDist a = ND (V.Vector a) deriving (Show, Eq, Ord, Read)

mean (ND xs) = xs V.! 0
stdDev (ND xs) = xs V.! 1
normalDistr mu std = ND (V.fromList [mu, std])

instance Functor NormalDist where
  fmap f (ND v) = ND $ V.map f v

instance Dist NormalDist SampleDouble where
  transform d epsilon = mean d + stdDev d * epsilon
  resample _d gen = MWCD.standard gen

instance VariUtil NormalDist where
  zipDist f (ND x1) (ND x2) = ND $ V.zipWith f x1 x2
  norm (ND x) = norm x

instance VariLogic NormalDist SampleDouble where
  logProb d x = (-xm * xm / (2 * sd * sd)) - ndPdfDenom
    where
      xm = x - mean d
      sd = stdDev d
      ndPdfDenom = log $ m_sqrt_2_pi * sd

  paramGradOfLogQ d x = ND $ V.fromList [(x - mu) / std, 1 / std ^ (3 :: Int) * (x - mu) ^ (2 :: Int) - 1 / std]
    where
      mu = mean d
      std = stdDev d

instance Differentiable NormalDist SampleDouble where
  sampleGradOfLogQ d z = -(z - mean d)/(stdDev d ** 2)
  gradTransform _d epsilon = ND $ V.fromList [1.0 , epsilon] -- grad (\d' -> transform d' (auto epsilon)) d -- ND $ V.fromList [1.0 , epsilon] --

gradientScore ::
     VariLogic a b
  => Node (a)
  -> (Weight, V.Vector Double, Samples b)
  -> PropNode (a) b
gradientScore = gradient f
  where
    f Node {..} nFactors s l =
      fmap
        (* (l + nFactors / weight * (logProb prior s - logProb dist s)))
        (paramGradOfLogQ dist s)

gradient f q@(Node{..}) (nFactors, like, samples) =
  U (memory', zipDist (*) rho' grad)
  where
    summed =
      V.foldl' (zipDist (+)) (fmap (const 0.0) dist) $
      V.zipWith
        (f q nFactors)
        samples
        like
    grad = fmap (/ (fromIntegral $ V.length samples)) summed
    (memory', rho') = rhoF grad q

-- | TODO: speed up by calc length in one pass
gradientReparam ::
     Differentiable a b
  => Node (a)
  -> (Weight, V.Vector Double, Samples b)
  -> PropNode (a) b
gradientReparam = gradient f
  where
    f Node {..} nFactors s l =
      fmap
        (* (l +
            nFactors / weight *
            (sampleGradOfLogQ prior (transform dist s) -
             sampleGradOfLogQ dist (transform dist s))))
        (gradTransform dist s)


rho alpha eta tau epsilon grad Node{..} =
  ( deltaM
  , zipDist
      (\ds s ->
         eta * (fromIntegral time) ** (negate 0.5 + epsilon) *
         (1.0 / (tau + sqrt (s + ds))))
      deltaM memory)
  where
    deltaM = zipDist (\g s -> alpha * g ^ (2 :: Int) - alpha * s) grad memory

-- | TODO: try making obs a Node and see if it still performant
-- implement weight on each propgator, and then divide by each
-- variationa distributions' factor (I was doing this wrong).
qProp ::
     VariLogic a b
  => ((a Double) -> ST s (Weight, V.Vector Double, Samples b))
  -> Cell s (PropNode (a) b)
  -> ST s ()
qProp likeFunc q =
  watch q $ \(N q') -> write q =<< (gradientScore q' <$> likeFunc (dist q'))

qPropAD ::
     (Differentiable a b)
  => ((a Double) -> ST s (Weight, V.Vector Double, Samples b))
  -> Cell s (PropNode (a) b)
  -> ST s ()
qPropAD likeFunc q =
  watch q $ \(N q') -> write q =<< (gradientReparam q' <$> likeFunc (dist q'))

qProp2 ::
     (VariLogic a c, VariLogic b d)
  => ((a Double, b Double) -> ST s (Weight, V.Vector Double, (Samples c, Samples d)))
  -> (Cell s (PropNode (a) c), Cell s (PropNode (b) d))
  -> ST s ()
qProp2 likeFunc (q1, q2) =
  watch q1 $ \(N q1') ->
    watch q2 $ \(N q2') -> do
      (w, l, (s1, s2)) <- likeFunc (dist q1', dist q2')
      write q1 (gradientScore q1' (w, l, s1))
      write q2 (gradientScore q2' (w, l, s2))

qPropAD2 ::
     (Differentiable a c, Differentiable b d)
  => (((a Double), (b Double)) -> ST s ((Weight, V.Vector Double, Samples c), (Weight, V.Vector Double, Samples d)))
  -> (Cell s (PropNode (a) c), Cell s (PropNode (b) d))
  -> ST s ()
qPropAD2 likeFunc (q1, q2) =
  watch q1 $ \(N q1') ->
    watch q2 $ \(N q2') -> do
      (l1, l2) <- likeFunc (dist q1', dist q2')
      write q1 (gradientReparam q1' l1)
      write q2 (gradientReparam q2' l2)

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
             V.map (mean . paramGradOfLogQ (normalDistr (transform q eps) std)) xs)
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
               (mean . paramGradOfLogQ (normalDistr (transform q eps) std))
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
            (normalDistr 0 0)
            (fromIntegral nSamp)
            qDist
            prior
            (rho alpha eta tau epsilon)) :: PropNode NormalDist Double)
    xProp <-
      known $
      (N (Node 1 (O V.empty) 0 xDist (O V.empty) (\_x1 _x2 -> (O V.empty, O V.empty))) :: PropNode Obs Double)
  -- qNormalProp 1.0 gen1 nSamp xs q
  -- qPropAD (normalLikeAD nSamp 1.0 gen1 xs) q
    qPropAD2 (normalLikeAD2 nSamp nSamp 1.0 gen1) (xProp, q) -- (V.length xs)
    (N q') <- fromMaybe (error "impos") <$> content q
    (N xs') <- fromMaybe (error "impos") <$> content xProp
    return (dist q', time q', Samp.mean xs, time xs')

genNormal = do
  gen <- create
  xs <- V.replicateM 1000 (transform (normalDistr 5.0 3.0) <$> resample (normalDistr 5.0 3.0) gen)
  return xs

genMixture = do
  gen <- create
  let theta' = MWCD.categorical (V.fromList [5.0, 10.0]) gen
  let std = 1.0
  let mixtures =
        V.fromList
          [MWCD.normal 2.0 std gen :: IO Double, MWCD.normal (-2.0) std gen]
  xs <- V.replicateM 1000 ((mixtures V.!) =<< theta')
  return xs

-- mixtureFit xs =
--   runST $ do
--     genG <- create
--     gen1 <- initialize =<< V.replicateM 256 (uniform genG)
--     let priorTheta = Diri (V.fromList [1.0, 1.0])
--     let priorBeta = V.fromList [normalDistr 0.0 1.0, normalDistr 0.0 1.0]
--     let nSamp = 100
--     let qDist = normalDistr 0.0 2.0
--     let alpha = 0.1 -- from kuckelbier et al
--     let eta = 0.1 -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
--     let tau = 1.0
--     let epsilon = 1e-16 -- (fromIntegral $ V.length xs)
--     let xDist = (O xs)
--     q <-
--       known $
--       (N (Node
--             1
--             (V.replicate 2 0)
--             (fromIntegral nSamp)
--             qDist
--             prior
--             (rho alpha eta tau epsilon)) :: PropNode NormalDist Double)
--     xProp <-
--       known $
--       (N (Node 1 (V.empty) 0 xDist (O V.empty) (\_x1 _x2 -> (V.empty, V.empty))) :: PropNode Obs Double)
--   -- qNormalProp 1.0 gen1 nSamp xs q
--   -- qPropAD (normalLikeAD nSamp 1.0 gen1 xs) q
--     qPropAD2 (normalLikeAD2 nSamp nSamp 1.0 gen1) (xProp, q) -- (V.length xs)
--     (N q') <- fromMaybe (error "impos") <$> content q
--     (N xs') <- fromMaybe (error "impos") <$> content xProp
--     return (dist q', time q', Samp.mean xs, time xs')

someFunc :: IO ()
someFunc = do
  let xs = runST $ genNormal
  putStrLn (show $ normalFit xs)
-- >>> someFunc
--
