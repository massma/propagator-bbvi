module Main where

import           Control.Monad.ST               ( runST )
import qualified Data.Vector                   as V
import           Statistics.BBVI.Examples

main :: IO ()
main = -- let xs = runST $ genMixedMem in putStrLn (show $ mixedMemFit xs)
  -- let xs = runST $ genMixture
  -- in  (putStrLn . unlines . fmap show . take 100 . V.toList) xs
  --       >> putStrLn (showTuple2 $ mixtureFit xs)
  let xs = runST $ genNormal
  in  (putStrLn . unlines . fmap show . take 100 . V.toList) xs
        >> putStrLn (showTuple2 $ normalFit xs)

showTuple2 (x1, x2) = show x1 <> "\n" <> show x2 <> "\n"
showTuple5 (x1, x2, x3, x4, x5) =
  show x1
    <> "\n"
    <> show x2
    <> "\n"
    <> show x3
    <> "\n"
    <> show x4
    <> "\n"
    <> show x5
