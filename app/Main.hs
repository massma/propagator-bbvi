module Main where

import           Statistics.BBVI.Examples
import           Control.Monad.ST               ( runST )

main :: IO ()
main = -- let xs = runST $ genMixedMem in putStrLn (show $ mixedMemFit xs)
  let xs = runST $ genMixture in putStrLn (show $ mixtureFit xs)
