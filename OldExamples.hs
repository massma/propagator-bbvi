----------------------
-- mixed membership
----------------------
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

genMixedMem :: ST s (V.Vector (V.Vector (Double))) -- Maybe Double
genMixedMem = do
  gen <- create
  -- initialize dirichlet that favors sparse categoricals
  let diriParam  = 0.01
  let thetaParam = MWCD.dirichlet (V.replicate nState diriParam) gen
  thetasTrue <- V.replicateM nTime thetaParam
  let thetas' = V.map (\theta' -> MWCD.categorical theta' gen) thetasTrue
  let betaStd = 0.001
  let betas' = V.generate
        nState
        (\k -> V.replicate
          nDimension
          (MWCD.normal
            ((fromIntegral k - (fromIntegral (nState - 1) * 0.5)) * 5)
            betaStd
            gen
          )
        )
  xs <- V.generateM
    nTime
    (\day -> V.generateM
      nDimension
      (\loc -> (thetas' V.! day) >>= \z -> ((betas' V.! z) V.! loc)) -- Just <$>
    )
  return xs

mixedMemFit xs = runST $ do
  genG <- create
  gen1 <- initialize =<< V.replicateM 256 (uniform genG)
  let priorTheta = dirichlet (V.replicate nState 1.0)
  let priorBeta = normalDistr 0 (5 * (fromIntegral nState - 1) :: Double)
  let nSamp      = 100
  let nObs       = (100 :: Int)
  -- let localStep  = 20
  let std        = 1.0
  qBetas <- known =<< V.generateM
    nState
    (\_i -> V.replicateM
      nDimension
      (do
        mu <- resample priorBeta gen1
        return
          (defaultNormalDist { dist    = normalDistr mu 1.0 -- (std :: Double)
                             , maxStep = globalMaxStep
                             , delta   = globalDelta -- 0.0000001
                             , prior   = priorBeta
                             , weight  = 1
                             }
          )
      )
    )
  qThetas <- known $ V.replicate
    nTime
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
