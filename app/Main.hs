module Main where

import           Control.Monad.ST               ( runST )
import qualified Data.Vector                   as V
import           Statistics.BBVI.Examples

main :: IO ()
main = -- let xs = runST $ genMixedMem in putStrLn (show $ mixedMemFit xs)
  let xs = runST $ genMixture
  in  (putStrLn . unlines . fmap show . take 100 . V.toList) xs
        >> putStrLn (show $ mixtureFit xs)
