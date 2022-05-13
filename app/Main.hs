module Main where

import Alignment
import Rendering
import System.Environment
import Graphics.Gloss
import Prelude hiding (Left)

main :: IO ()
main = do
    (seq1:seq2:_) <- getArgs
    let (galigned1, galigned2) = globalAlign seq1 seq2
    let (laligned1, laligned2) = localAlign seq1 seq2
    putStrLn "Global alignment:"
    putStrLn galigned1
    putStrLn galigned2
    putStrLn "Local alignment:"
    putStrLn laligned1
    putStrLn laligned2

    -- display (InWindow "Matrices" (200, 200) (10, 10)) white (showGlobalMatrix seq1 seq2)
