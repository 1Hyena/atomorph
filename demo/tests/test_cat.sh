#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

# TEST: Walking Cat
TITLE="cat"
../../atomorph -f test_${TITLE} -i ../data --blobs-as-texture cat_1.png cat_2.png cat_3.png cat_4.png cat_5.png cat_6.png -F 48 -t 5 --verbose -m 4 -M 5 -O 10 -D 2 -B 32 -c 2 -p 5 -z 1 --blend-blobs
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x6 -geometry 64x64+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten -quality 50 montage_${TITLE}.jpg
rm montage_${TITLE}.png

