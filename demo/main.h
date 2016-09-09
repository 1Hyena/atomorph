/*
 * AtoMorph - Command line program for morphing images with AtoMorph Library.
 * See Copyright Notice at the end of this file.
 */

#include "../atomorph.h"
#include "options.h"

bool init        (int argc, char **argv);
bool fill_morph  (am::morph *morph, size_t frame, std::vector<unsigned char> *image, unsigned width);
void write_image (am::morph *morph, size_t frame, unsigned width, unsigned height);
bool load_files  (am::morph *morph);
void save_files  (am::morph *morph);
void main_loop   (am::morph *morph);

/*
The MIT License (MIT)

Copyright (c) 2013-2014 Erich Erstu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
