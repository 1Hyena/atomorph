#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

# TEST 1: Volatile Blobs 
# AtoMorph does not require all blobs to
# have correspondences on all key frames.
TITLE="volatile"
../../atomorph -f test_${TITLE} -i ../data --blobs-as-texture volatile_1.png volatile_2.png volatile_3.png volatile_4.png -F 32 --verbose -c 1 -p 0 -z 0 -O 1 -D 1
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x4 -geometry +1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten montage_${TITLE}.jpg
rm montage_${TITLE}.png


# TEST 2: Empty Key Frames
# AtoMorph intuitively handles the situations
# where some key frames are completely empty.
TITLE="empty"
../../atomorph -f test_${TITLE} -i ../data --blobs-as-texture RGB_1.png empty.png squares_4.png volatile_1.png volatile_4.png -F 40 --verbose -c 1 -p 1 -z 1 -O 1 -D 1
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x5 -geometry +1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten montage_${TITLE}.jpg
rm montage_${TITLE}.png


# TEST 3: Simple Fluid
TITLE="fluid_simple"
../../atomorph -f test_${TITLE} -i ../data --blobs-as-texture volatile_1.png volatile_4.png -F 48 --verbose -c 1 -p 0 -z 0 -M 1 -O 1 -D 1 -L 10
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x6 -geometry +1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten montage_${TITLE}.jpg
rm montage_${TITLE}.png

# TEST 4: Chaotic Fluid
TITLE="fluid_chaotic"
../../atomorph -f test_${TITLE} -i ../data --blobs-as-texture volatile_1.png RGB_1.png squares_1.png -F 48 --verbose -c 0 -p 1 -z 0 -M 1 -O 1 -D 1 -L 10 --blob-feather 1
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x6 -geometry +1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten montage_${TITLE}.jpg
rm montage_${TITLE}.png

