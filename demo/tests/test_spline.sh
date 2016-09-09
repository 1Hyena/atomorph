#!/bin/bash
# ImageMagick is needed for these tests to run properly.

cd output

# TEST: Spline Interpolation
../../atomorph -f test_motion_none -i ../data --blobs-as-texture dots_1.png dots_2.png dots_3.png dots_4.png -F 80 --verbose -c 0 -p 0 -z 1 -O 1 -D 1 --fading-linear --motion-none
../../atomorph -f test_motion_linear -i ../data --blobs-as-texture dots_1.png dots_2.png dots_3.png dots_4.png -F 80 --verbose -c 0 -p 0 -z 1 -O 1 -D 1 --fading-linear --motion-linear
../../atomorph -f test_motion_spline -i ../data --blobs-as-texture dots_1.png dots_2.png dots_3.png dots_4.png -F 80 --verbose -c 0 -p 0 -z 1 -O 1 -D 1 --fading-linear --motion-spline
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose none test_motion_none_*.png test_motion_none.gif
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose none test_motion_linear_*.png test_motion_linear.gif
convert -channel A -threshold 99% -delay 10 -loop 0 -dispose none test_motion_spline_*.png test_motion_spline.gif
convert test_motion_none_*.png -layers flatten montage_motion_none.png
convert test_motion_linear_*.png -layers flatten montage_motion_linear.png
convert test_motion_spline_*.png -layers flatten montage_motion_spline.png
rm test_motion_none_*.png
rm test_motion_linear_*.png
rm test_motion_spline_*.png
convert montage_motion_none.png -background white -flatten montage_motion_none.jpg
convert montage_motion_linear.png -background white -flatten montage_motion_linear.jpg
convert montage_motion_spline.png -background white -flatten montage_motion_spline.jpg
rm montage_motion_none.png
rm montage_motion_linear.png
rm montage_motion_spline.png


