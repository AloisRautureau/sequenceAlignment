module Main where

import Lib
import System.Environment

main :: IO ()
main = do
    (seq1:seq2:_) <- getArgs
    let (galigned1, galigned2) = globalAlign seq1 seq2
    let (laligned1, laligned, laligned2) = localAlign seq1 seq2
    putStrLn "Global alignment:"
    putStrLn galigned1
    putStrLn galigned2
    putStrLn "Local alignment:"
    putStrLn laligned1
    putStrLn laligned
    putStrLn laligned2
