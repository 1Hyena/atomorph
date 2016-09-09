#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

# TEST: Many Armors
TITLE="armors"
../../atomorph -f test_${TITLE} -i ../data/armors --blobs-as-texture -F 120 --verbose -m -1 -O 60 -D 2 a-0.png a-1.png a-2.png a-3.png a-4.png a-5.png a-6.png a-7.png a-8.png a-9.png a-10.png a-11.png a-12.png a-13.png a-14.png a-15.png a-16.png a-17.png a-18.png a-19.png
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 10x12 -geometry 64x64+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten -quality 50 montage_${TITLE}.jpg
rm montage_${TITLE}.png

