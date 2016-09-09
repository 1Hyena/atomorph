#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

# TEST: Transformer
TITLE="transformer"
../../atomorph -f test_${TITLE} -i ../data --blobs-as-texture transformer_1.png transformer_2.png -F 32 --verbose -b 128 -M 5 -O 30 -D 2 -p 3 -c 4 -z 2 --finite
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x4 -geometry 64x64+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten -quality 50 montage_${TITLE}.jpg
rm montage_${TITLE}.png

