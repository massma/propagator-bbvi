module Test () where

import Prelude
import Data.Propagator
import Control.Monad.ST
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import System.Random.MWC (create, uniform, initialize)
import qualified System.Random.MWC.Distributions as MWCD
import qualified Statistics.Sample as Samp
import Simpler

------------------
-- Normal
------------------

genNormal :: ST s (V.Vector Double)
genNormal = do
  gen <- create
  xs <- V.replicateM 1000 (resample (normalDistr (5.0 :: Double) 3.0) gen)
  return xs

normalFit xs =
  runST $ do
    genG <- create
    gen1 <- initialize =<< V.replicateM 256 (uniform genG)
    let nSamp = 100
    xProp <- known $ defaultObs xs
    q <-
      known $
      (defaultNormalDist
         { prior = normalDistr 0.0 2.0
         , dist = normalDistr 0.0 2.0
         , weight = fromIntegral nSamp
         })
    (\qP xP ->
       watch qP $ \q' ->
         with xP $ \xs' ->
           normalLikeAD2 nSamp 1.0 gen1 (xs', q') >>= write q)
      q
      xProp
    -- qPropAD2 (normalLikeAD2 nSamp nSamp 1.0 gen1) (xProp, q) -- (V.length xs)
    q' <- fromMaybe (error "impos") <$> content q
    xs' <- fromMaybe (error "impos") <$> content xProp
    return (dist q', time q', Samp.mean xs, time xs')

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

------------------
--  Dirichlet
------------------
genDirichlet = do
  gen <- create
  xs <-
    V.replicateM
      1000
      (resample (dirichlet (V.fromList [10.0, 20.0])) gen >>= \cat ->
         MWCD.categorical (cat :: V.Vector Double) gen)
  return (xs :: V.Vector Int)

dirichletFit xs =
  runST $ do
    genG <- create
    gen1 <- initialize =<< V.replicateM 256 (uniform genG)
    let priorTheta = Diri (V.fromList [1.0, 1.0])
    let nSamp = 100
    let nClusters = 2
    qTheta <-
      -- resample (dirichlet (V.replicate nClusters 1)) gen1 >>= \startTh ->
        (known $
         ((defaultDirichlet priorTheta)
            { maxStep = globalMaxStep
            , delta = globalDelta -- 0.0000001
            , dist = dirichlet (V.fromList [(1.0 :: Double), 1.0]) -- startTh
            , weight = 1
            }))
    (\tP ->
       watch tP $ \theta' -> do
        upTh <- dirichletLike nSamp gen1 xs theta'
        write tP upTh)
      qTheta
    thetaF <- unsafeContent qTheta
    return (dist thetaF, time thetaF)

dirichletLike nSamp gen xs diriQ =
  V.replicateM nSamp (resample diri gen) >>= \samples ->
    return $
      gradientScore
        diriQ
        ( 1
        , V.map
          (\z ->
             V.sum $ V.map (\i -> log (z V.! i)) xs)
          samples
        , samples)
  where
    diri = dist diriQ


------------------
--  mixture
------------------
mixtureLike nSamp nObs std gen xsN thetaN betasN = do
  obs         <- V.replicateM nObs (resample xs gen)
  thetaSamp   <- V.replicateM nSamp (resample theta gen)
  betaSamples <- V.replicateM nSamp (epsilon (betas V.! 0) gen)
  let
    (likes, gradLikes) = V.unzip $ V.zipWith
      (\eps th ->
        V.foldl1' (\(x1, y1) (x2, y2) -> (x1 + x2, V.zipWith (+) y1 y2)) $ V.map
          (\x -> grad'
            (\bs -> logSum
              (V.map (auto . log) th)
              (V.map (\mu -> diffableNormalLogProb mu (auto std) (auto x)) bs)
            )
            (V.map (\d -> transform d eps) betas)
          )
          (obs :: V.Vector Double)
      )
      betaSamples
      thetaSamp
  return
    $ ( gradientScore thetaN (1, likes, thetaSamp)
      , V.imap
        (\i d -> gradientReparam d (1, V.map (V.! i) gradLikes, betaSamples))
        betasN
      )
 where
    -- obs = xsN
  xs    = dist xsN
  theta = dist thetaN
  betas = V.map dist betasN

logSum v1 = log . V.sum . V.map exp . V.zipWith (+) v1

mixtureLikeScore nSamp nObs std gen xsN thetaN betasN = do
  obs         <- V.replicateM nObs (resample xs gen)
  thetaSamp   <- V.replicateM nSamp (resample theta gen)
  betaSamples <- V.replicateM nSamp (V.mapM (\b -> resample b gen) betas)
  let
    likes = V.zipWith
      (\bs th -> V.sum $ V.map
        (\x -> logSum (V.map log th)
                      (V.map (\mu -> logProb (normalDistr mu std) x) bs)
        )
        (obs :: V.Vector Double)
      )
      betaSamples
      thetaSamp
  return
    $ ( gradientScore thetaN (1, likes, thetaSamp)
      , V.imap (\i d -> gradientScore d (1, likes, V.map (V.! i) betaSamples))
               betasN
      )
 where
    -- obs = xsN
  xs    = dist xsN
  theta = dist thetaN
  betas = V.map dist betasN

genMixture :: ST s (V.Vector Double)
genMixture = do
  gen <- create
  let theta' = MWCD.categorical (V.fromList [5.0, 10.0]) gen
  let std    = 1.0
  let mixtures =
        V.fromList [MWCD.normal 5.0 std gen, MWCD.normal (-5.0) std gen]
  xs <- V.replicateM 1000 ((mixtures V.!) =<< theta')
  return xs

mixtureFit xs = runST $ do
  genG <- create
  gen1 <- initialize =<< V.replicateM 256 (uniform genG)
  let priorTheta = dirichlet (V.fromList [0.1, 0.1])
  let priorBeta  = normalDistr 0.0 4.0
  let nSamp      = 10
  let nObs       = 100
  let localStep  = 20
  let nClusters  = 2
  qBetas <- known =<< V.generateM
    nClusters
    (\i -> do
         -- mu <- resample priorBeta gen1
         -- really I am doing empiracle bayes here ( would do mu - 1std, mu + 1 std), this approach makes things easier for testing
      let mu = if i == 0 then 2 else (negate 2)
      return
        (defaultNormalDist { dist    = normalDistr mu (1.0 :: Double)
                           , maxStep = globalMaxStep
                           , delta   = globalDelta -- 0.0000001
                           , prior   = normalDistr mu (1.0 :: Double) -- priorBeta
                           , weight  = 1
                           }
        )
    )
  qThetas <-
    -- resample (dirichlet (V.replicate nClusters 1)) gen1 >>= \startTh ->
    ( known
    $ ((defaultDirichlet priorTheta) { maxStep = globalMaxStep
                                     , delta   = globalDelta -- 0.0000001
                                     , dist    = dirichlet (V.replicate 2 1.0) -- startTh
                                     , weight  = 1
                                     }
      )
    )
  xProp <- known $ defaultObs xs
  -- (\tP bPs xP ->
  --    watch tP $ \theta' ->
  --      with xP $ \xs' ->
  --        with bPs $ \betas' -> do
  --          (upTh, upB) <-
  --            mixtureLike nSamp nObs 1.0 gen1 xs' theta' betas'
  --          write bPs upB
  --          write tP upTh)
  --   qThetas
  --   qBetas
  --   xProp
  (\tP bPs0 xP -> watch tP $ \theta' ->
      with xP
        $ (\xs' -> do
            bPs <- known =<< (V.map (initLocal localStep) <$> unsafeContent bPs0)
            watch bPs $ \betas' -> do
              (_upTh, upB) <- mixtureLikeScore nSamp nObs 1.0 gen1 xs' theta' betas'
              write bPs upB
            bPsNew <- unsafeContent bPs
            write bPs0 bPsNew
          )
    )
    qThetas
    qBetas
    xProp
  (\tP0 bPs xP -> watch bPs $ \betas' ->
      with xP
        $ (\xs' -> do
            tP <- known =<< (initLocal localStep <$> unsafeContent tP0)
            watch tP
              $ (\theta' -> do
                  (upTh, _upB) <- mixtureLikeScore nSamp
                                                 nObs
                                                 1.0
                                                 gen1
                                                 xs'
                                                 theta'
                                                 betas'
                  write tP upTh
                )
            tPNew <- unsafeContent tP
            write tP0 tPNew
          )
    )
    qThetas
    qBetas
    xProp
  thetaF <- unsafeContent qThetas
  betaF  <- unsafeContent qBetas
  let betaDists = V.map dist betaF
  return
    ( alphas (dist thetaF)
    , time thetaF
    , V.map time betaF
    , V.map (\d -> (mean d, stdDev d)) betaDists
    )
