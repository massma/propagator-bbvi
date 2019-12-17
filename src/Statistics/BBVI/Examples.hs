{-# LANGUAGE FlexibleContexts #-}
module Statistics.BBVI.Examples
  ( genNormal
  , normalFit
  , genDirichlet
  , dirichletFit
  , genMixture
  , mixtureFit
  , genMixedMem
  , mixedMemFit
  )
where

import           Control.Monad.ST
import           Data.Maybe                     ( fromMaybe )
import           Data.Propagator
import qualified Data.Vector                   as V
import           Numeric.AD                     ( grad
                                                , grad'
                                                , auto
                                                , diff
                                                )
import           System.Random.MWC              ( create
                                                , uniform
                                                , initialize
                                                , GenST
                                                )
import qualified System.Random.MWC.Distributions
                                               as MWCD
import qualified Statistics.Sample             as Samp
import           Statistics.BBVI

--- helper functions
unsafeContent = (fromMaybe (error "impos") <$>) . content
initLocal step p = p { maxStep = time p + step }

------------------
-- Normal
------------------
genNormal :: ST s (V.Vector Double)
genNormal = do
  gen <- create
  xs  <- V.replicateM 1000 (resample (normalDistr (5.0 :: Double) 3.0) gen)
  return xs

normalFit xs = runST $ do
  genG <- create
  gen1 <- initialize =<< V.replicateM 256 (uniform genG)
  let nSamp = 100
  xProp <- known $ defaultObs xs
  q     <-
    known
      $ (defaultNormalDist { prior  = normalDistr 0.0 2.0
                           , dist   = normalDistr 0.0 2.0
                           , weight = fromIntegral nSamp
                           }
        )
  (\qP xP -> watch qP
      $ \q' -> with xP $ \xs' -> normalLike nSamp 1.0 gen1 (xs', q') >>= write q
    )
    q
    xProp
  q'  <- fromMaybe (error "impos") <$> content q
  xs' <- fromMaybe (error "impos") <$> content xProp
  return (dist q', time q', Samp.mean xs, time xs')

normalLike nSamp std gen (xsN, qN) = do
  obs     <- V.replicateM nSamp (resample xs gen)
  samples <- V.replicateM nSamp (epsilon q gen)
  let gradLikes = V.map
        (\eps -> V.sum $ V.map
          -- paramGradOfLogQ should be faster than using AD, but hacky
          ((V.! 0) . paramGradOfLogQ (normalDistr (transform q eps) std))
          (obs :: V.Vector Double)
        )
        samples
  return $ gradientReparam qN ((fromIntegral nSamp), gradLikes, samples)
 where
  xs = dist xsN
  q  = dist qN


------------------
--  Dirichlet
------------------
genDirichlet = do
  gen <- create
  xs  <- V.replicateM
    1000
    (   resample (dirichlet (V.fromList [10.0, 20.0])) gen
    >>= \cat -> MWCD.categorical (cat :: V.Vector Double) gen
    )
  return (xs :: V.Vector Int)

dirichletFit xs = runST $ do
  genG <- create
  gen1 <- initialize =<< V.replicateM 256 (uniform genG)
  let priorTheta = dirichlet (V.fromList [1.0, 1.0])
  let nSamp      = 100
  let nClusters  = 2
  qTheta <-
    -- resample (dirichlet (V.replicate nClusters 1)) gen1 >>= \startTh ->
    ( known
    $ ((defaultDirichlet priorTheta)
        { maxStep = globalMaxStep
        , delta   = globalDelta -- 0.0000001
        , dist    = dirichlet (V.fromList [(1.0 :: Double), 1.0]) -- startTh
        , weight  = 1
        }
      )
    )
  (\tP -> watch tP $ \theta' -> do
      upTh <- dirichletLike nSamp gen1 xs theta'
      write tP upTh
    )
    qTheta
  thetaF <- unsafeContent qTheta
  return (dist thetaF, time thetaF)

dirichletLike nSamp gen xs diriQ =
  V.replicateM nSamp (resample diri gen) >>= \samples -> return $ gradientScore
    diriQ
    (1, V.map (\z -> V.sum $ V.map (\i -> log (z V.! i)) xs) samples, samples)
  where diri = dist diriQ

-----------------
--  mixture
------------------
mixtureLike
  :: (Dist b SampleVector, Differentiable c Double)
  => Int
  -> Int
  -> Double
  -> GenST s
  -> PropNode (Obs Double)
  -> PropNode b
  -> V.Vector (PropNode c)
  -> ST s (PropNode b, V.Vector (PropNode c))
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
              (_upTh, upB) <- mixtureLikeScore nSamp
                                               nObs
                                               1.0
                                               gen1
                                               xs'
                                               theta'
                                               betas'
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


-- mixed membership
mixedMemLike
  :: (Dist b SampleVector, Differentiable c Double)
  => Int
  -> Int
  -> Double
  -> GenST s
  -> V.Vector (V.Vector Double)
  -> V.Vector (PropNode b)
  -> V.Vector (V.Vector (PropNode c))
  -> ST
       s
       ( V.Vector (PropNode b)
       , V.Vector (V.Vector (PropNode c))
       )
mixedMemLike nSamp nObs std gen obss thetasN betassN = do
  thetaSamples <- V.replicateM nSamp (V.mapM (\th -> resample th gen) thetas)
  betaSamples  <- V.replicateM nSamp (epsilon ((betass V.! 0) V.! 0) gen)
  let
    (thetaLikes, gradLikesTranspose) = V.unzip $ V.zipWith
      (\eps ths ->
        let
          (fulLikes, unSummedGradTs) = V.unzip $ V.zipWith
            (\obs th ->
              let
                aThs           = V.map log th
                (likes, gradT) = V.unzip $ V.imap
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
              in
                (V.sum likes, gradT)
            )
            (obss :: V.Vector (V.Vector Double)) -- Maybe Double
            ths
        in  (fulLikes, V.foldl1' (V.zipWith (V.zipWith (+))) unSummedGradTs)
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

logSum v1 = log . minimumProb . V.sum . V.map exp . V.zipWith (+) v1

-- | should dynamically run below to figure out the range
-- >>> minimumDouble
-- -323.0
minimumDouble =
  ( (+ 1)
    . head
    . dropWhile (not . isInfinite . log . ((10 :: Double) **))
    $ intList
  , (+ 1)
    . head
    . dropWhile (not . isInfinite . ((1 :: Double) /) . (10 **))
    $ intList
  )
  where intList = [-100, -101 ..]

minimumProb :: (Floating a, Ord a) => a -> a
minimumProb = max 1e-308

mixedMemLikeScore nSamp nObs std gen obss thetasN betassN = do
  thetaSamples <- V.replicateM nSamp (V.mapM (\th -> resample th gen) thetas)
  betaSamples  <- V.replicateM nSamp
                               (V.mapM (V.mapM (\b -> resample b gen)) betass)
  let
    (thetaLikes, betaLikes) = V.unzip $ V.zipWith
      (\bss ths -> V.unzip $ V.zipWith
        (\obs th ->
          let aThs   = V.map log th
              logVec = V.imap
                (\loc ob -> logSum
                  aThs
                  (V.map (\mu -> logProb (normalDistr mu std) ob)
                         (V.map (V.! loc) bss)
                  )
                )
                obs
          in  (V.sum logVec, logVec)
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
          (\loc d -> gradientScore
            d
            -- trusting compiler to recognize below is repeated for each c
            ( 1
            , V.map (V.sum . V.map (V.! loc)) betaLikes
            , V.map ((V.! loc) . (V.! c)) betaSamples
            )
          )
          clusters
        )
        betassN
      )
 where
  thetas = V.map dist thetasN
  betass = V.map (V.map dist) betassN

nLocs :: Int
nLocs = 1 -- 9 * 36
nDays :: Int
nDays = 1000 -- 118 * 31
nStates :: Int
nStates = 2 -- 10

genMixedMem :: ST s (V.Vector (V.Vector (Double))) -- Maybe Double
genMixedMem = do
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

mixedMemFit xs = runST $ do
  genG <- create
  gen1 <- initialize =<< V.replicateM 256 (uniform genG)
  let priorTheta = dirichlet (V.replicate nStates 1.0)
  -- let priorBeta  = normalDistr 0.0 4.0
  let nSamp      = 100
  let nObs       = (100 :: Int)
  -- let localStep  = 20
  let std        = 1.0
  qBetas <- known =<< V.generateM
    nStates
    (\i ->
      let mu' = if i == 0 then -1 else 1
      in  V.replicateM
            nLocs
            (do
              return
                (defaultNormalDist { dist    = normalDistr mu' 1.0 -- (std :: Double)
                                   , maxStep = globalMaxStep
                                   , delta   = globalDelta -- 0.0000001
                                   , prior   = normalDistr mu' (4.0 :: Double) -- priorBeta
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
      (upTh, upB) <- mixedMemLike nSamp nObs std gen1 xs theta' betas'
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
  --             (_upTh, upB) <- mixedMemLikeScore nSamp nObs std gen1 xs' theta' betas'
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
  --                 (upTh, _upB) <- mixedMemLikeScore nSamp
  --                                                nObs
  --                                                strd
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
    ( Samp.mean . V.map (fromIntegral . time) $ thetaF
    , V.map (Samp.mean . V.map (fromIntegral . time)) betaF
    , V.map (Samp.mean . V.map mean) betaDists
    , V.map (Samp.mean . V.map stdDev) betaDists
    )
