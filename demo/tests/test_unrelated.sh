#!/bin/bash
# ImageMagick is needed for these tests to run properly.

TITLE="unrelated"

cd output

# TEST: A Morph of Unrelated Photos using Fluid Dynamics
../../atomorph -f test_${TITLE} -i ../data flowers.png merilin.png --blobs-as-texture --verbose --blend-blobs -c 2 -p 3 -z 1 -M 10 -O 50 --finite --blob-feather 4 --keep-background -F 48 -B 256
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose background test_${TITLE}_*.png test_${TITLE}.gif

convert +append test_${TITLE}_*.png strip_${TITLE}.png
rm test_${TITLE}_*.png
convert strip_${TITLE}.png -crop 256x256 +repage test_${TITLE}_%04d.png

for i in 1 2
do    
    rm test_${TITLE}_*1.png
    rm test_${TITLE}_*3.png
    rm test_${TITLE}_*5.png
    rm test_${TITLE}_*7.png
    rm test_${TITLE}_*9.png
    convert +append test_${TITLE}_*.png strip_${TITLE}.png
    rm test_${TITLE}_*.png    
    convert strip_${TITLE}.png -crop 256x256 +repage test_${TITLE}_%04d.png
done
rm strip_${TITLE}.png

montage test_${TITLE}_*.png -tile 4x3 -geometry 128x128+1+1 -frame 3 -background white -bordercolor white montage_${TITLE}.png
rm test_${TITLE}_*.png
convert montage_${TITLE}.png -background white -flatten -quality 50 montage_${TITLE}.jpg
rm montage_${TITLE}.png


