#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

TITLE="combine"

# TEST: Combine Layers
# This test shows that individual morphs
# can be combined together with ImageMagick.
../../atomorph -f test_body -i ../data/battlelord --blobs-as-texture body_1.png body_2.png body_3.png body_4.png body_5.png body_6.png -F 48 --verbose -m -1 -O 60 -D 2
convert +append test_body_*.png strip_body.png
rm test_body_*.png

../../atomorph -f test_head -i ../data/battlelord --blobs-as-texture head_1.png head_2.png head_3.png head_4.png head_5.png head_6.png -F 48 --verbose -m -1 -O 30 -D 2
convert +append test_head_*.png strip_head.png
rm test_head_*.png

../../atomorph -f test_lleg -i ../data/battlelord --blobs-as-texture lleg_1.png lleg_2.png lleg_3.png lleg_4.png lleg_5.png lleg_6.png -F 48 --verbose -m -1 -O 5 -D 2
convert +append test_lleg_*.png strip_lleg.png
rm test_lleg_*.png

../../atomorph -f test_lfoot -i ../data/battlelord --blobs-as-texture lfoot_1.png lfoot_2.png lfoot_3.png lfoot_4.png lfoot_5.png lfoot_6.png -F 48 --verbose -m -1 -O 5 -D 2
convert +append test_lfoot_*.png strip_lfoot.png
rm test_lfoot_*.png

../../atomorph -f test_rleg -i ../data/battlelord --blobs-as-texture rleg_1.png rleg_2.png rleg_3.png rleg_4.png rleg_5.png rleg_6.png -F 48 --verbose -m -1 -O 20 -D 2
convert +append test_rleg_*.png strip_rleg.png
rm test_rleg_*.png

../../atomorph -f test_rfoot -i ../data/battlelord --blobs-as-texture rfoot_1.png rfoot_2.png rfoot_3.png rfoot_4.png rfoot_5.png rfoot_6.png -F 48 --verbose -m -1 -O 5 -D 2
convert +append test_rfoot_*.png strip_rfoot.png
rm test_rfoot_*.png

../../atomorph -f test_rarm -i ../data/battlelord --blobs-as-texture rarm_1.png rarm_2.png rarm_3.png rarm_4.png rarm_5.png rarm_6.png -F 48 --verbose -m -1 -O 30 -D 2
convert +append test_rarm_*.png strip_rarm.png
rm test_rarm_*.png

../../atomorph -f test_rhand -i ../data/battlelord --blobs-as-texture rhand_1.png rhand_2.png rhand_3.png rhand_4.png rhand_5.png rhand_6.png -F 48 --verbose -m -1 -O 30 -D 2
convert +append test_rhand_*.png strip_rhand.png
rm test_rhand_*.png

convert strip_lfoot.png strip_rfoot.png strip_lleg.png strip_rleg.png strip_body.png strip_rarm.png strip_rhand.png strip_head.png -background none -flatten strip_${TITLE}.png
rm strip_body.png
rm strip_head.png
rm strip_lleg.png
rm strip_lfoot.png
rm strip_rleg.png
rm strip_rfoot.png
rm strip_rarm.png
rm strip_rhand.png
convert strip_${TITLE}.png -crop 128x128 +repage test_${TITLE}_%04d.png
rm strip_${TITLE}.png
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
montage test_${TITLE}_*.png -tile 8x6 -geometry 64x64+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten montage_${TITLE}.jpg
rm montage_${TITLE}.png

