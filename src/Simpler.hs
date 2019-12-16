{-# LANGUAGE FlexibleContexts #-}
module Simpler
  ( someFunc
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
import qualified Statistics.Sample             as StSa
import qualified System.Random.MWC.Distributions
                                               as MWCD
import           Text.Printf                    ( printf )
import           Statistics.BBVI

type Time = Int

type Samples = V.Vector

initLocal step p = p { maxStep = time p + step }
initLocalDefault p = initLocal (maxStep p) p

unsafeContent = (fromMaybe (error "impos") <$>) . content

-- mixed membership
-- mixedLike nSamp nObs std gen obss thetasN betassN = do
--   thetaSamples <- V.replicateM nSamp (V.mapM (\th -> resample th gen) thetas)
--   betaSamples  <- V.replicateM nSamp (epsilon ((betass V.! 0) V.! 0) gen)
--   let
--     (thetaLikes, gradLikesTranspose) = V.unzip $ V.zipWith
--       (\eps ths ->
--         let
--           (fulLikes, unSummedGradTs) = V.unzip $ V.zipWith
--             (\obs th ->
--               let
--                 aThs           = V.map log th
--                 (likes, gradT) = V.unzip $ V.imap
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
--               in
--                 (V.sum likes, gradT)
--             )
--             (obss :: V.Vector (V.Vector Double)) -- Maybe Double
--             ths
--         in  (fulLikes, V.foldl1' (V.zipWith (V.zipWith (+))) unSummedGradTs)
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

mixedLikeScore nSamp nObs std gen obss thetasN betassN = do
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
            -- something is wrong belwo, before this was indexed by cluster but that came from grad
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
      (upTh, upB) <- mixedLikeScore nSamp nObs std gen1 xs theta' betas'
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
  --             (_upTh, upB) <- mixedLikeScore nSamp nObs std gen1 xs' theta' betas'
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
  putStrLn . unlines . fmap show . V.toList $ xs
  putStrLn (show $ mixedFit xs)
-- >>> someFunc
--
