#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

# TEST: Battle Lord
TITLE="battlelord"
../../atomorph -f test_${TITLE} -i ../data --blobs-as-texture battlelord_1.png battlelord_2.png battlelord_3.png battlelord_4.png battlelord_5.png battlelord_6.png -F 48 --verbose -m -1 -O 60 -D 2
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x6 -geometry 64x64+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten -quality 50 montage_${TITLE}.jpg
rm montage_${TITLE}.png

