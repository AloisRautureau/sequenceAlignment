module Rendering (showGlobalMatrix) where

import Graphics.Gloss
import Alignment
import Data.Array
import Prelude hiding (Left)

arrow :: Picture
arrow = pictures [tail, head]
    where tail = rectangleSolid 3 10
          head = polygon [(-5, 5), (5, 5), (0, 10)]

cell :: (Cost, Alignment.Path) -> Picture
cell (c, p) = pictures [costText, border, pathArrow]
    where costText = textPos (show c) $ scale 0.125 0.125 $ text (show c)
          textPos t = translate ((-5) + (-2.5) * (fromIntegral $ length t)) (-7)
          border = rectangleWire 40 40
          pathArrow
            | p /= [] = translate (-10) (10) $ arrowRotation arrow
            | otherwise = blank
          arrowRotation = case last p of
                Match -> rotate (-45)
                Mismatch -> rotate (-45)
                Up -> rotate 0
                Left -> rotate (-90)

row :: [(Cost, Alignment.Path)] -> Picture
row [] = blank
row (c:cs) = pictures [cell c, translate 40 0 (row cs)]

seqHorizontal :: Sequence -> Picture
seqHorizontal [] = blank
seqHorizontal (c:cs) = pictures [textPic, translate 40 0 (seqHorizontal cs)] 
      where textPic = textPos (show c) $ scale 0.125 0.125 $ text (show c)
            textPos t = translate ((-5) + (-2.5) * (fromIntegral $ length t)) (-7)

seqVertical :: Sequence -> Picture
seqVertical [] = blank
seqVertical (c:cs) = pictures [textPic, translate 0 (-40) (seqVertical cs)] 
      where textPic = textPos (show c) $ scale 0.125 0.125 $ text (show c)
            textPos t = translate ((-5) + (-2.5) * (fromIntegral $ length t)) (-7)

showGlobalMatrix :: Sequence -> Sequence -> Picture
showGlobalMatrix seq1 seq2 = pictures (descript : [translate 0 ((-40)*fromIntegral i) (row (elems r)) | (i, r) <- rows])
      where matrix = needlemanWunschMatrix seq1 seq2
            rows = assocs matrix
            descript = pictures [translate 40 40 (seqHorizontal seq1), translate (-40) (-40) (seqVertical seq2)]