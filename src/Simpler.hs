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
import Numeric.AD (grad, grad', auto, diff)
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

class Sampleable a c where
  resample :: a Double -> GenST s -> ST s c

class (Functor a) => DistUtil a where
  zipDist :: (Double -> Double -> Double) -> a Double -> a Double -> a Double
  norm :: a Double -> Double

class (Functor a, Sampleable a c, DistUtil a) => Dist a c where
  logProb :: a Double -> c -> Double
  paramGradOfLogQ ::  a Double -> c -> a Double -- gradient of parameters evauated at some sample of x

class (Functor a) => Differentiable a where
  diffableLogProb :: (Floating b) => a b -> b -> b
  gradTransform :: (Floating b) => a b -> b -> a b
  sampleGradOfLogQ :: (Floating b) => a b -> b -> b -- gradient of a sample evaluate with params of q
  transform :: (Floating b) => a b -> b -> b
  epsilon ::  (Floating b) => a b -> GenST s -> ST s b

normVec = sqrt . V.sum . V.map (^ (2 :: Int))

type Time = Int

maxStep :: Time
maxStep = 100000 -- 4664 -- 100000

type Samples = V.Vector

data PropNode a b = N (Node a b) | U { memoryUpdate :: !(a b)
                                     , gradientUpdate :: !(a b)
                                     } -- memory, gradient

data Node a b = Node
  { time :: !Time
  , memory :: !(a b)
  , weight :: !b
  , dist :: !(a b)
  , prior :: !(a b)
  , rhoF :: !(a b -> Node a b -> (a b, a b)) -- memory, gradient
  }

newtype Obs a = O (V.Vector a) deriving (Show, Eq, Ord, Read)

fromO (O v) = v

instance Sampleable Obs Double where
  resample (O d) gen =
    return . (d V.!) =<< uniformR (0, (V.length d - 1)) gen

instance Functor Obs where
  fmap f (O d) = O $ V.map f d

instance Dist Obs Double where
  logProb _d _x = 0
  paramGradOfLogQ _d _x = O V.empty

instance DistUtil Obs where
  zipDist _f _o1 _o2 = O V.empty
  norm _ = 0.0

instance Differentiable Obs where

newtype Dirichlet a = Diri (V.Vector a) deriving (Show, Eq, Ord, Read)

dirichlet = Diri

logB :: Dirichlet Double -> Double
logB (Diri alphas) = V.sum (V.map logGamma alphas) - logGamma (V.sum alphas)

instance Sampleable Dirichlet (V.Vector Double) where
  resample (Diri diri) gen = MWCD.dirichlet diri gen

instance Functor Dirichlet where
  fmap f (Diri diri) = Diri $ V.map f diri

instance Foldable Dirichlet where
  foldr f x0 (Diri diri) = V.foldr f x0 diri

instance Traversable Dirichlet where
  traverse f (Diri diri) = Diri <$> traverse f diri

instance Dist Dirichlet (V.Vector Double) where
  logProb (Diri diri) cat =
    V.sum (V.zipWith (\alpha x -> (alpha - 1) * log x) diri cat) -
    logB (Diri diri)
  paramGradOfLogQ (Diri diri) cat =
    Diri $ V.zipWith (\a x -> summed - digamma a + log x) diri cat
    where
      summed = digamma (V.sum cat)

instance DistUtil Dirichlet where
  zipDist f (Diri x1) (Diri x2) = Diri $ V.zipWith f x1 x2
  norm (Diri x1) = normVec x1

instance DistUtil a => Propagated (PropNode a Double) where
  merge (N node) (U { .. })
    | norm gradientUpdate < 0.00000001 = Change False (N node) -- 0.00001
    | time node >= maxStep = Change False (N node)
    | otherwise = Change True updateNode
    where
      updateNode =
        N
          (node
             { time = (time node) + 1
             , memory = zipDist (+) (memory node) memoryUpdate
             , dist = newQ
             })
        where
          newQ = zipDist (+) (dist node) gradientUpdate
  merge (U _ _) _ = Contradiction mempty "Trying to update a gradient"
  merge (N _) (N _) = Contradiction mempty "Trying overwrite a node"

newtype NormalDist a = ND (V.Vector a) deriving (Show, Eq, Ord, Read)

fromPropNode (N node) = node
fromPropNode (U _ _) = error "called from prop node on update"

mean (ND xs) = xs V.! 0
stdDev (ND xs) = xs V.! 1
normalDistr mu std = ND (V.fromList [mu, std])

instance Functor NormalDist where
  fmap f (ND v) = ND $ V.map f v

instance Foldable NormalDist where
  foldr f x0 (ND diri) = V.foldr f x0 diri

instance Traversable NormalDist where
  traverse f (ND diri) = ND <$> traverse f diri

instance Sampleable NormalDist Double where
  resample d gen = MWCD.normal (mean d) (stdDev d) gen -- transform d <$> epsilon d gen

instance Differentiable NormalDist where
  diffableLogProb d x = (-xm * xm / (2 * sd * sd)) - ndPdfDenom
    where
      xm = x - mean d
      sd = stdDev d
      ndPdfDenom = log $ sqrt (2*pi) * sd
  transform d eps = mean d + stdDev d * eps
  epsilon _d gen = realToFrac <$> MWCD.standard gen
  sampleGradOfLogQ d z = -(z - mean d)/(stdDev d ** 2)
  gradTransform _d eps = ND $ V.fromList [1.0 , eps] -- grad (\d' -> transform d' (auto eps)) d

instance DistUtil NormalDist where
  zipDist f (ND x1) (ND x2) = ND $ V.zipWith f x1 x2
  norm (ND x) = normVec x

instance Dist NormalDist Double where
  logProb d x = (-xm * xm / (2 * sd * sd)) - ndPdfDenom
    where
      xm = x - mean d
      sd = stdDev d
      ndPdfDenom = log $ m_sqrt_2_pi * sd

  paramGradOfLogQ d x = ND $ V.fromList [(x - mu) / std, 1 / std ^ (3 :: Int) * (x - mu) ^ (2 :: Int) - 1 / std]
    where
      mu = mean d
      std = stdDev d

gradientScore ::
     Dist a c
  => Node a Double
  -> (Double, V.Vector Double, Samples c)
  -> PropNode a Double
gradientScore = gradient f
  where
    f Node {..} nFactors s l =
      fmap
        (* (l + nFactors / weight * (logProb prior s - logProb dist s)))
        (paramGradOfLogQ dist s)

gradient f q@(Node{..}) (nFactors, like, samples) =
  U memory' (zipDist (*) rho' grad)
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
     (DistUtil a, Differentiable a)
  => Node a Double
  -> (Double, V.Vector Double, Samples Double)
  -> PropNode a Double
gradientReparam = gradient f
  where
    f Node {..} nFactors s l =
      fmap
        (* (l +
            nFactors / weight *
            (sampleGradOfLogQ prior (transform dist s) -
             sampleGradOfLogQ dist (transform dist s))))
        (gradTransform dist s)

rho alpha eta tau eps gra Node{..} =
  ( deltaM
  , zipDist
      (\ds s ->
         eta * (fromIntegral time) ** (negate 0.5 + eps) *
         (1.0 / (tau + sqrt (s + ds))))
      deltaM memory)
  where
    deltaM = zipDist (\g s -> alpha * g ^ (2 :: Int) - alpha * s) gra memory

normalLike nSamp std gen xs q =
  V.replicateM nSamp (resample q gen) >>= \samples ->
    return
      ( fromIntegral $ V.length xs
      , V.map
          (\z ->
             V.sum $ V.map (logProb (normalDistr z std)) xs)
          samples
      , samples)

-- normalLike2 nSamp std gen (xs, q) = do
--   samples <- V.replicateM nSamp (resample q gen)
--   obs <- V.replicateM nSamp (resample xs gen)
--   let like = V.map
--                (\eps ->
--                   V.sum $ V.map (logProb (normalDistr (transform q eps) std)) obs)
--                samples
--   return
--     ( U (O V.empty, O V.empty)
--     , gradientScore q (fromIntegral nObs, like, samples))

-- | specialize all fmaps to vector???
-- CAREFUL: with rao-blackweixation and my wieghting approach, all terms
-- in log likelihood functions must contain all nodes q
normalLikeAD nSamp std gen xs q =
  V.replicateM nSamp (epsilon q gen) >>= \samples ->
    return
      ( fromIntegral $ V.length xs
      , V.map
          (\eps ->
             V.sum $
             V.map (mean . paramGradOfLogQ (normalDistr (transform q eps) std)) xs)
          samples
      , samples)

normalLikeAD2 nSamp std gen (xsN, qN) = do
  obs <- V.replicateM nSamp (resample xs gen)
  samples <- V.replicateM nSamp (epsilon q gen)
  let gradLikes =
        V.map
          (\eps ->
             V.sum $
             V.map
               (mean . paramGradOfLogQ (normalDistr (transform q eps) std))
               (obs :: V.Vector Double))
          samples
  return $ gradientReparam qN ((fromIntegral nSamp), gradLikes, samples)
  where
    xs = dist xsN
    q = dist qN

normalFit xs =
  runST $ do
    genG <- create
    gen1 <- initialize =<< V.replicateM 256 (uniform genG)
    let prior = normalDistr 0.0 2.0
    let nSamp = 100
    let qDist = normalDistr 0.0 2.0
    let alpha = 0.1 -- from kuckelbier et al
    let eta = 0.01 -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
    let tau = 1.0
    let eps = 1e-16 -- (fromIntegral $ V.length xs)
    let xDist = (O xs)
    q <-
      known $
      N
        (Node
           1
           (normalDistr 0 0)
           (fromIntegral nSamp)
           qDist
           prior
           (rho alpha eta tau eps))
    xProp <-
      known $
      N
        (Node
           1
           (O V.empty)
           0
           xDist
           (O V.empty)
           (\_x1 _x2 -> (O V.empty, O V.empty)))
    (\qP xP ->
       watch qP $ \(N q') ->
         with xP $ \(N xs') -> normalLikeAD2 nSamp 1.0 gen1 (xs', q') >>= write q)
      q
      xProp
    -- qPropAD2 (normalLikeAD2 nSamp nSamp 1.0 gen1) (xProp, q) -- (V.length xs)
    (N q') <- fromMaybe (error "impos") <$> content q
    (N xs') <- fromMaybe (error "impos") <$> content xProp
    return (dist q', time q', Samp.mean xs, time xs')

mixedFit xs =
  runST $ do
    genG <- create
    gen1 <- initialize =<< V.replicateM 256 (uniform genG)
    let priorTheta = Diri (V.fromList [0.1, 0.1])
    let priorBeta = normalDistr 0.0 4.0
    let nSamp = 10
    let alpha = 0.1 -- from kuckelbier et al
    let eta = 0.01  -- 1 -- 10 -- 100 -- 0.01 -- this needs tuning
    let tau = 1.0
    let eps = 1e-16 -- (fromIntegral $ V.length xs)
    let xDist = (O xs)
    let nClusters = 2
    qBetas <- V.generateM nClusters
      (\_i -> do
          mu <- resample priorBeta gen1
          known $
           N
            (Node
               1
               (normalDistr 0 0)
               (fromIntegral nSamp)
               (normalDistr mu 2.0)
               priorBeta
               (rho alpha eta tau eps)))
    qThetas <- resample (dirichlet (V.replicate nClusters 1)) gen1 >>= \startTh ->
      (known $
           N
            (Node
               1
               (dirichlet $ V.replicate nClusters 0.0)
               (fromIntegral (nSamp))
               (dirichlet startTh) -- (dirichlet $ V.replicate nClusters 0.1)
               priorTheta
               (rho alpha eta tau eps)))
    xProp <-
      known $
      N
        (Node
           1
           (O V.empty)
           0
           xDist
           (O V.empty)
           (\_x1 _x2 -> (O V.empty, O V.empty)))
    (\tP bPs xP ->
       watch tP $ \(N theta') ->
        with xP $ \(N xs') -> do
           betasP <- V.mapM ((fromPropNode . fromMaybe (error "impos") <$>) . content) bPs
           (upTh, upB) <- mixedLike nSamp 1.0 gen1 xs' theta' betasP -- xs'
           V.mapM_ (uncurry write) $ V.zip bPs upB
           write tP upTh)
      qThetas
      qBetas
      xProp
    -- qPropAD2 (normalLikeAD2 nSamp nSamp 1.0 gen1) (xProp, q) -- (V.length xs)
    (N thetaF) <- fromMaybe (error "impos") <$> content qThetas
    betaF <- V.mapM ((fromPropNode . fromMaybe (error "impos") <$>) . content) qBetas
    return (dist thetaF, time thetaF, V.map time betaF, V.map dist betaF)

mixedLike nSamp std gen xsN thetaN betasN = do
  obs <- V.replicateM nSamp (resample xs gen)
  thetaSamp <- V.replicateM nSamp (resample theta gen)
  betaSamples <- V.replicateM nSamp (epsilon (betas V.! 0) gen)
  let (likes, gradLikes) =
        V.unzip $
        V.zipWith
          (\eps th ->
             V.foldl1' (\(x1, y1) (x2, y2) -> (x1 + x2, V.zipWith (+) y1 y2)) $
             V.map
               (\x ->
                  grad'
                    (\bs ->
                       logSum
                         (V.map (auto . log) th)
                         (V.map
                            (\mu ->
                               diffableLogProb
                                 (normalDistr mu (auto std))
                                 (auto x))
                            bs))
                    (V.map (\d -> transform d eps) betas))
               (obs :: V.Vector Double))
          betaSamples
          thetaSamp
  return $
    ( gradientScore thetaN (fromIntegral nSamp, likes, thetaSamp)
    , V.imap
        (\i d ->
           gradientReparam
             d
             (fromIntegral nSamp, V.map (V.! i) gradLikes, betaSamples))
        betasN)
  where
    logSum v1 = V.sum . V.map exp . V.zipWith (+) v1
    -- obs = xsN
    xs = dist xsN
    theta = dist thetaN
    betas = V.map dist betasN

genNormal :: ST s (V.Vector Double)
genNormal = do
  gen <- create
  xs <- V.replicateM 1000 (resample (normalDistr (5.0 :: Double) 3.0) gen)
  return xs

genMixture :: ST s (V.Vector Double)
genMixture = do
  gen <- create
  let theta' = MWCD.categorical (V.fromList [5.0, 10.0]) gen
  let std = 1.0
  let mixtures =
        V.fromList
          [MWCD.normal 2.0 std gen, MWCD.normal (-2.0) std gen]
  xs <- V.replicateM 1000 ((mixtures V.!) =<< theta')
  return xs

someFunc :: IO ()
someFunc = do
  -- let xs = runST $ genNormal
  -- putStrLn (show $ normalFit xs)
  let xs = runST $ genMixture
  putStrLn $ unlines $ V.toList $ V.map show xs
  putStrLn (show $ mixedFit xs)
-- >>> someFunc
--
