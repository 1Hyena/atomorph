#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

# TEST 1: fruit
TITLE="fruit"
../../atomorph -f test_${TITLE} -i ../data/food --blobs-as-texture fruit_1.png fruit_2.png fruit_3.png -F 24 --verbose -m -1 -O 60 -D 2
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x3 -geometry 64x64+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten -quality 90 montage_${TITLE}.jpg
rm montage_${TITLE}.png

# TEST 2: berry
TITLE="berry"
../../atomorph -f test_${TITLE} -i ../data/food --blobs-as-texture berry_1.png berry_2.png berry_3.png berry_4.png berry_5.png berry_6.png -F 48 --verbose -m -1 -O 60 -D 2
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x6 -geometry 64x64+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten -quality 90 montage_${TITLE}.jpg
rm montage_${TITLE}.png

# TEST 3: lettuce
TITLE="lettuce"
../../atomorph -f test_${TITLE} -i ../data/food --blobs-as-texture lettuce_1.png lettuce_2.png lettuce_3.png lettuce_4.png lettuce_5.png lettuce_6.png -F 48 --verbose -m -1 -O 60 -D 2
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x6 -geometry 64x64+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten -quality 90 montage_${TITLE}.jpg
rm montage_${TITLE}.png

# TEST 4: spice
TITLE="spice"
../../atomorph -f test_${TITLE} -i ../data/food --blobs-as-texture spice_1.png spice_2.png spice_3.png spice_4.png spice_5.png spice_6.png -F 48 --verbose -m -1 -O 60 -D 2
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x6 -geometry 64x64+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten -quality 90 montage_${TITLE}.jpg
rm montage_${TITLE}.png
