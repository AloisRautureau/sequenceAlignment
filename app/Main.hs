module Main where

import Lib
import System.Environment

main :: IO ()
main = do
    (seq1:seq2:_) <- getArgs
    let (aligned1, aligned2) = align seq1 seq2
    putStrLn "Result:"
    putStrLn aligned1
    putStrLn aligned2
