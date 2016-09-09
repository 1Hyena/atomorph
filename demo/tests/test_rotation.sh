#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

# TEST 1: Rotation
# Intuitively, this test should generate a
# perfectly rotating star. However, instead
# it generates interestingly anomalous yet
# smooth rotation.
TITLE="rotation"
../../atomorph -f test_${TITLE} -i ../data --blobs-as-texture shape_1.png shape_2.png shape_3.png shape_4.png shape_5.png shape_6.png -F 48 --verbose -t 255 -M 1 -O 10 -D 2
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x6 -geometry 64x64+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten montage_${TITLE}.jpg
rm montage_${TITLE}.png

