/*
 * AtoMorph Demo - Simple Demo showing what AtoMorph is capable of doing.
 * See Copyright Notice at the end of this file.
 */

#include "../../atomorph.h"

#include <allegro5/allegro.h>
#include <allegro5/allegro_font.h>
#include <allegro5/allegro_image.h>
#include <allegro5/allegro_color.h>
#include <allegro5/allegro_primitives.h>

extern const float FPS;
extern const int SCREEN_W;
extern const int SCREEN_H;

extern ALLEGRO_DISPLAY     *display     ;
extern ALLEGRO_EVENT_QUEUE *event_queue ;
extern ALLEGRO_TIMER       *timer       ;
extern ALLEGRO_FONT        *font        ;

bool init(int argc, char **argv);
int calculate_fps();
void draw(AM_SCENE *scene, ALLEGRO_BITMAP *to, double t, double weight);
void render_morph(AM_IMAGE *img, ALLEGRO_BITMAP *to);
void blend_morphs(AM_BLENDER *b, ALLEGRO_BITMAP *to);
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
