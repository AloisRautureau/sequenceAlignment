module Lib (alignRecursive, align) where
import Data.List (sortBy)
import Data.Array

type Sequence = String

-- naive recursive version

hammingDistance :: Eq a => [a] -> [a] -> Int
hammingDistance x y = sum [1 | (x', y') <- x `zip` y, x' /= y']

alignRecursive :: Sequence -> Sequence -> (Sequence, Sequence)
alignRecursive [] seq2 = (replicate (length seq2) '-', seq2)
alignRecursive seq1 [] = (seq1, replicate (length seq1) '-')
alignRecursive seq1@(x:xs) seq2@(y:ys) = head $ sortBy alignementQuality possibilities
    where possibilities = [fuse (x, y) (alignRecursive xs ys), fuse ('-', y) (alignRecursive seq1 ys), fuse (x, '-') (alignRecursive xs seq2)]
          alignementQuality (a1, a2) (b1, b2) = compare (hammingDistance a1 a2) (hammingDistance b1 b2)
          fuse (x, y) (x', y') = (x:x', y:y')

-- Needleman-Wunsch algorithm

-- This could have been a similarity matrix, but it felt like forcibly added
-- complexity. This is a basic implementation after all.
type Cost = Int
matchCost = 3
mismatchCost = -1
indelCost = -2

data Alignement = Match | Mismatch | Up | Left deriving (Eq, Show, Ord)
type Path = [Alignement]

-- Given two sequences, returns a cost and an optimal path following
-- the Needleman-Wunsch algorithm
needlemanWunschPath :: Sequence -> Sequence -> Path
needlemanWunschPath seq1 seq2 = snd $ matrix ! (length seq2) ! (length seq1)
    where matrix = listArray (0, length seq2) [row i | i <- [0..length seq2]]
          -- In order to compute the value of a cell, we only need the row
          -- above it, and the left neighboring cell.
          -- Therefore, we only need to initialize first row, and first cell
          -- of each row after that.
          row 0 = listArray (0, length seq1) [(j * indelCost, replicate j Lib.Left) | j <- [0..length seq1]]
          row i = listArray (0, length seq1) $ (i * indelCost, replicate i Up) : [cell i j | j <- [1..length seq1]]
          cell i j = head $ sortBy scoreSort [
                update (matrix ! (i-1) ! (j-1)) (score i j), 
                update (matrix ! (i-1) ! j) (indelCost, Lib.Up),
                update (matrix ! i ! (j-1)) (indelCost, Lib.Left)
            ]
          score i j | (seq1 !! (j-1)) == (seq2 !! (i-1)) = (matchCost, Match) | otherwise = (mismatchCost, Mismatch)
          -- We sort possible move given possible score, and prioritize matches
          -- when the possible scores are equal.
          scoreSort (s,m) (s',m')
            | s == s' = compare m m'
            | otherwise = compare s' s
          update (score, path) (cost, move) = (score + cost, path ++ [move])

-- Align from a sequence of moves in the Needleman-Wunsch matrix
align :: Sequence -> Sequence -> (Sequence, Sequence)
align seq1 seq2 = withAlignement seq1 seq2 (needlemanWunschPath seq1 seq2)
    where withAlignement (x:xs) (y:ys) (alignment:as)
            | alignment == Up = appendTuple ('-', y) (withAlignement (x:xs) ys as)
            | alignment == Lib.Left = appendTuple (x, '-') (withAlignement xs (y:ys) as)
            | otherwise = appendTuple (x, y) (withAlignement xs ys as)
          withAlignement [] s _ = (replicate (length s) '-', s)
          withAlignement s [] _ = (s, replicate (length s) '-')
          appendTuple (a, b) (la, lb) = (a:la,b:lb)
