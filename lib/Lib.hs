module Lib (alignRecursive, globalAlign, localAlign) where
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
matchCost = 1
mismatchCost = -1
indelCost = -2

data Alignement = Match | Mismatch | Up | Left deriving (Eq, Show, Ord)
type Path = [Alignement]

-- Given two sequences, returns an optimal path following
-- the Needleman-Wunsch algorithm
needlemanWunschAlignment :: Sequence -> Sequence -> Path
needlemanWunschAlignment seq1 seq2 = snd $ matrix ! (length seq2) ! (length seq1)
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

-- Given two sequences, returns an optimal local alignment using the
-- Smith-Waterman algorithm
smithWatermanAlignment :: Sequence -> Sequence -> (Path, (Int, Int))
smithWatermanAlignment seq1 seq2 = ((snd $ matrix ! maxCellI ! maxCellJ), maxCellIx)
    where matrix = listArray (0, length seq2) [row i | i <- [0..length seq2]]
          row 0 = listArray (0, length seq1) (replicate (length seq1 + 1) (0, []))
          row i = listArray (0, length seq1) $ (0, []) : [cell i j | j <- [1..length seq1]]
          cell i j = head $ sortBy scoreSort [
                update (matrix ! (i-1) ! (j-1)) (score i j),
                update (matrix ! (i-1) ! j) (indelCost, Lib.Up),
                update (matrix ! i ! (j-1)) (indelCost, Lib.Left),
                (0, [])
            ]
          score i j | (seq1 !! (j-1)) == (seq2 !! (i-1)) = (matchCost, Match) | otherwise = (mismatchCost, Mismatch)
          scoreSort (s, m) (s', m')
            | s == s' = compare m m'
            | otherwise = compare s' s
          update (score, path) (cost, move) = (score + cost, path ++ [move])
          maxCellIx@(maxCellI, maxCellJ) = fst $ foldl1 f [((i, j), cell) | i <- [0..length seq2], let (j, cell) = maxCellRow i]
            where f x@(_, a) y@(_, b) | scoreSort a b == LT = x | otherwise = y
                  maxCellRow i = foldl1 f (assocs (matrix ! i))

-- Functions returning the aligned sequences
traceback :: Sequence -> Sequence -> Path -> Char -> (Sequence, Sequence)
traceback seq1@(x:xs) seq2@(y:ys) (a:as) c
    | a == Up = appendTuple (c, y) (traceback seq1 ys as c)
    | a == Lib.Left = appendTuple (x, c) (traceback xs seq2 as c)
    | otherwise = appendTuple (x, y) (traceback xs ys as c)
    where appendTuple (a, b) (la, lb) = (a:la, b:lb)
traceback [] s _ c = (replicate (length s) c, s)
traceback s [] _ c = (s, replicate (length s) c)

globalAlign :: Sequence -> Sequence -> (Sequence, Sequence)
globalAlign seq1 seq2 = traceback seq1 seq2 (needlemanWunschAlignment seq1 seq2) '-'

localAlign :: Sequence -> Sequence -> (Sequence, Sequence, Sequence)
localAlign seq1 seq2 = (seq1'', reverse seq2'', alignedSeq2)
    where seq1' = drop (length seq1 - ixSeq1) seq1
          seq2' = reverse $ take ixSeq2 seq2
          (path, (ixSeq2, ixSeq1)) = smithWatermanAlignment seq1 seq2
          (seq1'', seq2'') = traceback seq1' seq2' path ' '
          alignedSeq2 = replicate (length seq1 - length seq2') ' ' ++ seq2
