module Lib (alignRecursive, align) where
import Data.List (sortBy)
import Data.Array

type Sequence = String

hammingDistance :: Eq a => [a] -> [a] -> Int
hammingDistance x y = sum [1 | (x', y') <- x `zip` y, x' /= y']

alignRecursive :: Sequence -> Sequence -> (Sequence, Sequence)
alignRecursive [] seq2 = (replicate (length seq2) '-', seq2)
alignRecursive seq1 [] = (seq1, replicate (length seq1) '-')
alignRecursive seq1@(x:xs) seq2@(y:ys) = head $ sortBy alignementQuality possibilities
    where possibilities = [fuse (x, y) (alignRecursive xs ys), fuse ('-', y) (alignRecursive seq1 ys), fuse (x, '-') (alignRecursive xs seq2)]
          alignementQuality (a1, a2) (b1, b2) = compare (hammingDistance a1 a2) (hammingDistance b1 b2)
          fuse (x, y) (x', y') = (x:x', y:y')

matchCost = 3
mismatchCost = -1
indelCost = -2

data Alignement = Match | DelX | DelY | Mismatch | Init deriving (Eq, Show)

needlemanWunschPath :: Sequence -> Sequence -> [Alignement]
needlemanWunschPath seq1 seq2 = reverse $ pathFrom (length seq1) (length seq2)
    where cell 0 0 = (0, Init)
          cell i 0 = (i * indelCost, DelX)
          cell 0 j = (j * indelCost, DelY)
          cell i j = head $ sortBy (\(a, _) (b, _) -> compare b a) possibilities 
            where possibilities = [
                    (fst (cell (i-1) j) + indelCost, DelX),
                    (fst (cell i (j-1)) + indelCost, DelY),
                    (fst (cell (i-1) (j-1)) + matchScore i j, matchStatus i j)
                    ]
                  max3 a b = max a . max b
                  charAt n = last . take (n+1)
                  matchStatus i j | charAt (i-1) seq1 == charAt (j-1) seq2 = Match | otherwise = Mismatch
                  matchScore i j | matchStatus i j == Match = matchCost | otherwise = mismatchCost
          pathFrom i j = 
            let direction = snd (cell i j) 
            in direction : case direction of
                DelX -> pathFrom (i-1) j
                DelY -> pathFrom i (j-1)
                Match -> pathFrom (i-1) (j-1)
                Mismatch -> pathFrom (i-1) (j-1)
                Init -> []

align :: Sequence -> Sequence -> (Sequence, Sequence)
align seq1 seq2 = _align seq1 seq2 (needlemanWunschPath seq1 seq2)
    where _align (x:xs) (y:ys) (a:as) = case a of
                                            Match -> fuse (x, y) (_align xs ys as)
                                            Mismatch -> fuse (x, y) (_align xs ys as)
                                            DelY -> fuse ('-', y) (_align (x:xs) ys as)
                                            DelX -> fuse (x, '-') (_align xs (y:ys) as)
                                            Init -> _align (x:xs) (y:ys) as
          _align [] seq2 _ = alignRecursive [] seq2
          _align seq1 [] _ = alignRecursive seq1 []
          fuse (x, y) (x', y') = (x:x', y:y')
