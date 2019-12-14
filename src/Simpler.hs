{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Simpler
  ( someFunc
  , Differentiable(..)
  , Dist(..)
  , Dirichlet(..)
  , defaultDirichlet
  , dirichlet
  , NormalDist(..)
  , defaultNormalDist
  , normalDistr
  , Obs(..)
  , defaultObs
  , PropNode(..)
  , PropNodes
  , gradientScore
  , gradientReparam
  , unsafeContent

  -- * probably ones we want to remove
  , globalMaxStep
  , globalDelta
  )
where
import           Prelude
import           Data.Propagator
import           Control.Monad.ST
import           Data.Maybe                     ( fromMaybe )
import qualified Data.Vector                   as V
import           Numeric.AD                     ( grad
                                                , grad'
                                                , auto
                                                , diff
                                                )
import           Numeric.MathFunctions.Constants
                                                ( m_sqrt_2_pi )
import           Numeric.SpecFunctions          ( logGamma
                                                , digamma
                                                )
import           System.Random.MWC              ( create
                                                , uniform
                                                , initialize
                                                , GenST
                                                , uniformR
                                                )
-- import Statistics.Distribution
-- import Statistics.Distribution.Normal
import qualified Statistics.Sample             as StSa
import qualified System.Random.MWC.Distributions
                                               as MWCD

type SampleVector = V.Vector Double

type SampleDouble = Double

type ParamVector = V.Vector Double

class DistUtil a where
  fromParamVector :: ParamVector -> a
  toParamVector :: a -> ParamVector
  nParams :: a -> Int

class DistUtil a => Dist a c where
  resample :: a -> GenST s -> ST s c
  logProb :: a -> c -> Double
  paramGradOfLogQ ::  a -> c -> ParamVector -- gradient of parameters evauated at some sample of x

class (Dist a c) => Differentiable a c where
  gradTransform :: a -> c -> ParamVector
  sampleGradOfLogQ ::  a -> c -> c -- gradient of a sample evaluate with params of q
  transform :: a -> c -> c
  epsilon ::  a -> GenST s -> ST s c

norm = sqrt . V.sum . V.map (^ (2 :: Int))

type Time = Int

globalMaxStep :: Time
globalMaxStep = 1000

globalDelta :: Double
globalDelta = 1e-16 -- 0.00001 --

globalEta :: Double
globalEta = 1.0 -- 0.1 --10.0

type Gradient = V.Vector Double
type Memory = V.Vector Double

type Samples = V.Vector

data PropNode a
  = U { memoryUpdate :: !(V.Vector Double)
      , gradientUpdate :: !(V.Vector Double) }
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
-- | CAREFUL: below is dangerous if I start doing the ideas i thought
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

newtype Obs a = O (V.Vector a) deriving (Show, Eq, Ord, Read)

defaultObs :: V.Vector Double -> PropNode (Obs Double)
defaultObs xs =
  (Node 1
        globalMaxStep
        globalDelta
        (V.empty)
        0
        (O xs)
        (O V.empty)
        (\_x1 _x2 -> (V.empty, V.empty))
  )

instance Dist (Obs Double) Double where
  resample (O d) gen = return . (d V.!) =<< uniformR (0, (V.length d - 1)) gen
  logProb _d _x = 0
  paramGradOfLogQ _d _x = V.empty

instance DistUtil (Obs Double) where
  nParams _x = 0
  toParamVector _ = V.empty
  fromParamVector _ = O V.empty

newtype Dirichlet = Diri (V.Vector Double) deriving (Show, Eq, Ord, Read)

defaultDirichlet prior =
  (Node 1
        globalMaxStep
        globalDelta
        (fmap (const 0.0) (toParamVector prior))
        1
        prior
        prior
        (rhoKuc defaultKucP)
  )

dirichlet xs = Diri $ xs

alphas :: Dirichlet -> V.Vector Double
alphas (Diri xs) = xs

logB :: Dirichlet -> Double
logB d = V.sum (V.map logGamma as) - logGamma (V.sum as) where as = alphas d

instance DistUtil Dirichlet where
  nParams = V.length . alphas
  fromParamVector xs = dirichlet $ V.map exp xs
  toParamVector (Diri xs) = V.map log xs

instance Dist Dirichlet SampleVector where
  resample d gen = MWCD.dirichlet (alphas d) gen
  logProb d cat = V.sum (V.zipWith (\alpha x -> (alpha - 1) * log x) as cat)
    - logB d
    where as = alphas d
  paramGradOfLogQ d cat = V.zipWith
    (\a x -> a * (summed - digamma a + log x))
    as
    cat
   where
    as     = alphas d
    summed = digamma (V.sum as)

defaultNormalDist =
  (Node 1
        globalMaxStep
        globalDelta
        (V.replicate (nParams dflt) 0)
        1
        dflt
        dflt
        (rhoKuc defaultKucP)
  )
  where dflt = (normalDistr 0 1)

data NormalDist = ND {mean :: Double, stdDev :: Double} deriving (Show, Eq, Ord, Read)

normalDistr mu std = ND mu std

instance DistUtil NormalDist where
  fromParamVector xs = normalDistr (xs V.! 0) (exp (xs V.! 1))
  toParamVector d = V.fromList [mean d, log (stdDev d)]
  nParams _d = 2

instance Dist NormalDist SampleDouble where
  resample d gen = MWCD.normal (mean d) (stdDev d) gen
  logProb d x = (-xm * xm / (2 * sd * sd)) - ndPdfDenom
   where
    xm         = x - mean d
    sd         = stdDev d
    ndPdfDenom = log $ m_sqrt_2_pi * sd
  paramGradOfLogQ d x = V.fromList
    [(x - mu) / std, 1 / std ^ (2 :: Int) * (x - mu) ^ (2 :: Int) - 1]
   where
    mu  = mean d
    std = stdDev d
-- >>> paramGradOfLogQ (normalDistr 0.0 1.0) (2.0 :: Double)
-- [2.0,3.0]
-- >>> grad (\[mu, omega] -> diffableNormalLogProb mu (exp omega) (auto 2.0)) [0.0, log 1.0] :: [Double]
-- <interactive>:129:8-70: warning: [-Wincomplete-uni-patterns]
--     Pattern match(es) are non-exhaustive
--     In a lambda abstraction:
--         Patterns not matched:
--             []
--             [_]
--             (_:_:_:_)
-- [2.0,3.0]

diffableNormalLogProb mu sd x = (-xm * xm / (2 * sd * sd)) - ndPdfDenom
 where
  xm         = x - mu
  ndPdfDenom = log $ sqrt (2 * pi) * sd

instance Differentiable NormalDist SampleDouble where
  transform d eps = mean d + stdDev d * eps
  epsilon _d gen = MWCD.standard gen
  sampleGradOfLogQ d z = -(z - mean d) / (stdDev d ** 2)
  gradTransform d eps = V.fromList [1.0, eps * stdDev d] -- ND $ V.fromList [1.0 , eps * stdDev d] -- --
-- >>> (transform (normalDistr 0.0 2.0) 1.0 :: Double)
-- 2.0
-- >>> gradTransform (normalDistr 0.0 1.0) (2.0 :: Double)
-- [1.0,2.0]

gradientScore
  :: Dist a c
  => PropNode a
  -> (Double, V.Vector Double, Samples c)
  -> PropNode a
gradientScore = gradient f
 where
  f Node {..} nFactors s l = fmap
    (* (l + nFactors / weight * (logProb prior s - logProb dist s)))
    (paramGradOfLogQ dist s)

gradient f q@(Node {..}) (nFactors, like, samples) = U
  memory'
  (V.zipWith (*) rho' gr)
 where
  summed =
    V.foldl' (V.zipWith (+)) (V.replicate (nParams dist) 0.0)
      $ V.zipWith (f q nFactors) samples like
  gr              = fmap (/ (fromIntegral $ V.length samples)) summed
  (memory', rho') = rhoF gr q

-- | TODO: speed up by calc length in one pass
gradientReparam
  :: (Differentiable a Double)
  => PropNode a
  -> (Double, V.Vector Double, Samples Double)
  -> PropNode a
gradientReparam = gradient f
 where
  f Node {..} nFactors s l = fmap
    (* ( l
       + nFactors
       / weight
       * ( sampleGradOfLogQ prior (transform dist s)
         - sampleGradOfLogQ dist  (transform dist s)
         )
       )
    )
    (gradTransform dist s)

rhoKuc KucP {..} gra Node {..} =
  ( deltaM
  , V.zipWith
    (\ds s ->
      eta
        *  (fromIntegral time)
        ** (negate 0.5 + eps)
        *  (1.0 / (tau + sqrt (s + ds)))
    )
    deltaM
    memory
  )
 where
  deltaM = V.zipWith (\g s -> alpha * g ^ (2 :: Int) - alpha * s) gra memory

data KucP = KucP
  { alpha :: Double
  , eta :: Double
  , tau :: Double
  , eps :: Double
  } deriving (Show, Eq, Ord, Read)

-- | eta is what you probably want to tune: kucukelbir trys 0.01 0.1 1 10 100
defaultKucP = KucP { alpha = 0.1, eta = globalEta, tau = 1.0, eps = 1e-16 }

initLocal step p = p { maxStep = time p + step }
initLocalDefault p = initLocal (maxStep p) p

unsafeContent = (fromMaybe (error "impos") <$>) . content


-- mixed membership
mixedLike nSamp nObs std gen obss thetasN betassN = do
  thetaSamples <- V.replicateM nSamp (V.mapM (\th -> resample th gen) thetas)
  betaSamples  <- V.replicateM nSamp (epsilon ((betass V.! 0) V.! 0) gen)
  let
    (thetaLikes, gradLikesTranspose) = V.unzip $ V.zipWith
      (\eps ths -> V.unzip $ V.zipWith
        (\obs th ->
          let aThs = V.map log th
          in
            -- here grad is vector of loc x cluster, we want cluster x loc
            V.foldl1' (\(l1, g1) (l2, g2) -> (l1 + l2, V.zipWith (+) g1 g2))
              $ V.imap
                  (\i ob -> grad'
                    (\bs -> logSum
                      (V.map auto aThs)
                      (V.map
                        (\mu -> diffableNormalLogProb mu (auto std) (auto ob))
                        bs
                      )
                    )
                    (V.map (\v -> transform (v V.! i) eps) betass)
                  )
                  obs
        )
        (obss :: V.Vector (V.Vector Double)) -- Maybe Double
        ths
      )
      betaSamples
      thetaSamples
  return
    $ ( V.imap
        (\i d -> gradientScore
          d
          (1, V.map (V.! i) thetaLikes, (V.map (V.! i) thetaSamples))
        )
        thetasN
      , V.imap
        (\c clusters -> V.imap
          (\loc d -> gradientReparam
            d
            (1, V.map ((V.! c) . (V.! loc)) gradLikesTranspose, betaSamples)
          )
          clusters
        )
        betassN
      )
 where
  thetas = V.map dist thetasN
  betass = V.map (V.map dist) betassN

logSum v1 = log . V.sum . V.map exp . V.zipWith (+) v1

-- mixedLikeScore nSamp nObs std gen obss thetasN betassN = do
--   thetaSamples <- V.replicateM nSamp (V.mapM (\th -> resample th gen) thetas)
--   betaSamples  <- V.replicateM nSamp (epsilon ((betass V.! 0) V.! 0) gen)
--   let
--     (thetaLikes, gradLikesTranspose) = V.unzip $ V.zipWith
--       (\eps ths -> V.unzip $ V.zipWith
--         (\obs th ->
--           let aThs = V.map log th
--           in
--             -- here grad is vector of loc x cluster, we want cluster x loc
--             V.foldl1' (\(l1, g1) (l2, g2) -> (l1 + l2, V.zipWith (+) g1 g2))
--               $ V.imap
--                   (\i ob -> grad'
--                     (\bs -> logSum
--                       (V.map auto aThs)
--                       (V.map
--                         (\mu -> diffableNormalLogProb mu (auto std) (auto ob))
--                         bs
--                       )
--                     )
--                     (V.map (\v -> transform (v V.! i) eps) betass)
--                   )
--                   obs
--         )
--         (obss :: V.Vector (V.Vector Double)) -- Maybe Double
--         ths
--       )
--       betaSamples
--       thetaSamples
--   return
--     $ ( V.imap
--         (\i d -> gradientScore
--           d
--           (1, V.map (V.! i) thetaLikes, (V.map (V.! i) thetaSamples))
--         )
--         thetasN
--       , V.imap
--         (\c clusters -> V.imap
--           (\loc d -> gradientReparam
--             d
--             (1, V.map ((V.! c) . (V.! loc)) gradLikesTranspose, betaSamples)
--           )
--           clusters
--         )
--         betassN
--       )
--  where
--   thetas = V.map dist thetasN
--   betass = V.map (V.map dist) betassN

nLocs :: Int
nLocs = 1 -- 9 * 36
nDays :: Int
nDays = 1000 -- 118 * 31
nStates :: Int
nStates = 2 -- 10

genMixed :: ST s (V.Vector (V.Vector (Double))) -- Maybe Double
genMixed = do
  gen <- create
  -- initialize dirichlet that favors sparse categoricals
  let diriParam  = 0.01
  let thetaParam = MWCD.dirichlet (V.replicate nStates diriParam) gen
  thetasTrue <- V.replicateM nDays thetaParam
  let thetas' = V.map (\theta' -> MWCD.categorical theta' gen) thetasTrue
  let betaStd = 0.001
  let betas' = V.generate
        nStates
        (\k -> V.replicate
          nLocs
          (MWCD.normal
            ((fromIntegral k - (fromIntegral (nStates - 1) * 0.5)) * 5)
            betaStd
            gen
          )
        )
  xs <- V.generateM
    nDays
    (\day -> V.generateM
      nLocs
      (\loc -> (thetas' V.! day) >>= \z -> ((betas' V.! z) V.! loc)) -- Just <$>
    )
  return xs

mixedFit xs = runST $ do
  genG <- create
  gen1 <- initialize =<< V.replicateM 256 (uniform genG)
  let priorTheta = dirichlet (V.replicate nStates 0.1)
  -- let priorBeta  = normalDistr 0.0 4.0
  let nSamp      = 10
  let nObs       = (100 :: Int)
  -- let localStep  = 20
  qBetas <- known =<< V.generateM
    nStates
    (\i ->
      let mu' = if i == 0 then 2 else 3
      in  V.replicateM
            nLocs
            (do
                                                          -- mu <- resample priorBeta gen1
              return
                (defaultNormalDist { dist    = normalDistr mu' (1.0 :: Double)
                                   , maxStep = globalMaxStep
                                   , delta   = globalDelta -- 0.0000001
                                   , prior   = normalDistr mu' (1.0 :: Double) -- priorBeta
                                   , weight  = 1
                                   }
                )
            )
    )
  qThetas <- known $ V.replicate
    nDays
    ((defaultDirichlet priorTheta) { maxStep = globalMaxStep
                                   , delta   = globalDelta -- 0.0000001
                                   , weight  = 1
                                   }
    )
  (\tP bPs -> watch tP $ \theta' -> with bPs $ \betas' -> do
      (upTh, upB) <- mixedLike nSamp nObs 1.0 gen1 xs theta' betas'
      write bPs upB
      write tP  upTh
    )
    qThetas
    qBetas
  -- (\tP bPs0 xP -> watch tP $ \theta' ->
  --     with xP
  --       $ (\xs' -> do
  --           bPs <- known =<< (V.map (initLocal localStep) <$> unsafeContent bPs0)
  --           watch bPs $ \betas' -> do
  --             (_upTh, upB) <- mixedLikeScore nSamp nObs 1.0 gen1 xs' theta' betas'
  --             write bPs upB
  --           bPsNew <- unsafeContent bPs
  --           write bPs0 bPsNew
  --         )
  --   )
  --   qThetas
  --   qBetas
  --   xProp
  -- (\tP0 bPs xP -> watch bPs $ \betas' ->
  --     with xP
  --       $ (\xs' -> do
  --           tP <- known =<< (initLocal localStep <$> unsafeContent tP0)
  --           watch tP
  --             $ (\theta' -> do
  --                 (upTh, _upB) <- mixedLikeScore nSamp
  --                                                nObs
  --                                                1.0
  --                                                gen1
  --                                                xs'
  --                                                theta'
  --                                                betas'
  --                 write tP upTh
  --               )
  --           tPNew <- unsafeContent tP
  --           write tP0 tPNew
  --         )
  --   )
    -- qThetas
    -- qBetas
    -- xProp
  thetaF <- unsafeContent qThetas
  betaF  <- unsafeContent qBetas
  let betaDists = V.map (V.map dist) betaF
  return
    ( StSa.mean . V.map (fromIntegral . time) $ thetaF
    , V.map (StSa.mean . V.map (fromIntegral . time)) betaF
    , V.map (StSa.mean . V.map mean) betaDists
    , V.map (StSa.mean . V.map stdDev) betaDists
    )

someFunc :: IO ()
someFunc = do
  -- let xs = runST $ genNormal
  -- putStrLn (show $ normalFit xs)
  -- let xs = runST $ genMixture
  -- putStrLn $ unlines $ V.toList $ V.map show xs
  -- putStrLn (show $ mixtureFit xs)
  -- let xs = runST $ genDirichlet
  -- putStrLn (show $ dirichletFit xs)
  let xs = runST $ genMixed
  putStrLn (show $ mixedFit xs)
-- >>> someFunc
--
