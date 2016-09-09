#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

# TEST: Perlin Noise Dependent Colour Interpolation
TITLE="perlin"
../../atomorph -f test_${TITLE} -i ../data --blobs-as-texture RGB_1.png RGB_2.png RGB_3.png RGB_4.png -F 32 --verbose -c 0 -p 1 -z 0 -O 1 -D 1
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x4 -geometry +1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten montage_${TITLE}.jpg
rm montage_${TITLE}.png

