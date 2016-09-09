#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

# TEST: Doom Monsters
TITLE="monsters"
../../atomorph -f test_${TITLE} -i ../data/monsters --blobs-as-texture -F 112 --verbose -m -1 -O 60 -D 4 hero.png d.png cg.png cd.png baron.png imp.png mc.png pain.png spider.png vile.png sg.png revenant.png zm.png hk.png
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif
rm test_${TITLE}_*.png
../../atomorph -f test_${TITLE} -i ../data/monsters --blobs-as-texture -F 24 --verbose -m -1 -O 60 -D 4 hero.png d.png cg.png cd.png baron.png imp.png mc.png pain.png spider.png vile.png revenant.png hk.png
montage test_${TITLE}_*.png -tile 4x6 -geometry 128x128+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten -quality 50 montage_${TITLE}.jpg
rm montage_${TITLE}.png

