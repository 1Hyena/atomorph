/*
 * AtoMorph - Simple Library for Morphing 2D Particle Clouds
 * See Copyright Notice at the end of this file.
 */

#ifndef _ATOMORPH_H_
#define _ATOMORPH_H_

#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <future>
#include <random>
#include <map>
#include <set>
#include <limits>
#include <assert.h>

#ifdef ATOMORPH_OPENCV
#include <opencv/cv.h>
#endif

#include "spline.h"
#include "perlin.h"
#include "fluidmodel.h"
#include "color.h"

#ifdef ATOMORPH_DEPRECATED
#define AM_NONE                  0
#define AM_LINEAR                1
#define AM_COSINE                2
#define AM_PERLIN                3
#define AM_SPLINE                4

#define AM_WARN_POINTER_SIZE     1 // Set when sizeof(void *) is less than 64 bits.
#define AM_WARN_ATOM_SIZE        2 // Set when sizeof(AM_ATOM) is more than 64 bits.

typedef struct AM_COLOR {
    unsigned char r;
    unsigned char g;
    unsigned char b;
    unsigned char a;
} AM_COLOR;

typedef struct AM_ATOM {
    uint16_t x;
    uint16_t y;
    unsigned char  r;
    unsigned char  g;
    unsigned char  b;
    unsigned char  a;
} AM_ATOM;

class AM_SCENE {
    public:
    AM_SCENE();
    ~AM_SCENE();
    
    void clear();    
    bool init(size_t atoms, size_t frames);
    bool push_atom(size_t frame, AM_ATOM atom);    
    void renew_splines();
    double get_certainty(double t) const;
    void get_xy  (size_t atom, double t, double *x, double *y, unsigned method) const;   
    AM_COLOR get_rgba(size_t atom, double t, double lag, double slope, unsigned interpolation) const;    
    double get_current_path_length(size_t atom, double t) const;    
    AM_ATOM get_atom(size_t atom, size_t frame) const;
    size_t atom_count() const  {return atoms;}
    size_t frame_count() const {return frames;}
    size_t get_sorted_atom_at(size_t position);
    size_t get_current_frame(double t) const;
    bool copy_map_from(const AM_SCENE *scene);
    bool copy_candidates_from(const AM_SCENE *scene);
    void shuffle();
    void sort_atoms();
    double get_cost() const;
    double get_path_length(size_t atom) const;
    double get_path_color(size_t atom) const;    
    bool elect_atoms();        
    bool swap_atoms(size_t frame, size_t atom1, size_t atom2);
    const std::vector<AM_ATOM> *get_candidates(size_t frame) const;
    
    private:
    size_t atoms;
    size_t frames;
    
    std::vector<size_t> sorted_atoms;
    std::vector<AM_ATOM> *candidates;
    glnemo::CRSpline     *splines;
    
    AM_ATOM **map;
};

class AM_THREAD {
    public:
    AM_THREAD();
    ~AM_THREAD();
    
    bool clear();
    bool init(const AM_SCENE *scene);
    
    void set_seed(unsigned seed);
    void set_step_size(int step_size);
    void set_gradient_importance(double weight);
    void set_magic_exponent(double exponent);    
    
    bool is_running() const {return running;}
    bool is_paused()  const {return paused;}

    void stop()   { signal_stop  = true; while(running) std::this_thread::sleep_for(std::chrono::milliseconds(0)); step_thread.join();}
    void start()  { running      = true; step_thread = std::thread(&AM_THREAD::run, this);}
    void pause()  { signal_pause = true; while(!paused) std::this_thread::sleep_for(std::chrono::milliseconds(0)); signal_pause = false;}
    void resume() { paused       = false;}    
    
    double get_cost() const {return cost;}
    bool   fetch_scene(AM_SCENE *target) const;
    
    private:
    void run();
    void step();
    double chain_length(AM_ATOM a1, AM_ATOM a2, AM_ATOM a3);
    double chain_gradient(AM_ATOM a1, AM_ATOM a2, AM_ATOM a3);
    AM_SCENE scene;
    double *subcost;
    double  cost;
    std::default_random_engine e1;
    int step_size;
    double magic_exponent;
    double gradient_importance;
        
    std::atomic<bool> signal_stop;
    std::atomic<bool> running;
    std::atomic<bool> signal_pause;
    std::atomic<bool> paused;
    std::thread step_thread;    
};

class AM_IMAGE {
    public:
    AM_IMAGE();
    ~AM_IMAGE();
    
    bool set_scene(const AM_SCENE *scene);
    bool set_resolution(size_t width, size_t height);
    bool set_time(double time);
    bool set_seed(unsigned seed);
    bool set_color_interpolation(unsigned method);
    bool set_path_interpolation(unsigned method);    
        
    bool is_running() const {return running;}
    bool is_paused()  const {return paused;}
    size_t pixel_count() const  {return atoms.size();}    
    
    bool get_xy  (size_t pixel, int *x, int *y) const;
    bool get_rgba(size_t pixel, unsigned char *r, unsigned char *g, unsigned char *b, unsigned char *a) const;    
    size_t get_pixel_count() const {return atoms.size();}
    
    bool fetch_pixels(std::vector<AM_ATOM> *to) const;

    void stop()   { signal_stop  = true; while(running) std::this_thread::sleep_for(std::chrono::milliseconds(0)); step_thread.join();}
    void start()  { running      = true; step_thread = std::thread(&AM_IMAGE::run, this);}
    void pause()  { signal_pause = true; while(!paused) std::this_thread::sleep_for(std::chrono::milliseconds(0)); signal_pause = false;}
    void resume() { paused       = false;}            
        
    private:
    void run();
    void render();    
    std::vector<AM_ATOM> atoms;
    size_t w;
    size_t h;
    double t;
    AM_SCENE scene;    
    bool done;
    unsigned seed;
    unsigned color_interpolation;
    unsigned path_interpolation;
    PerlinNoise lag_map;
    PerlinNoise slope_map;    
    
    std::atomic<bool> signal_stop;
    std::atomic<bool> running;
    std::atomic<bool> signal_pause;
    std::atomic<bool> paused;
    std::thread step_thread;        
};

class AM_BLENDER {
    public:
    AM_BLENDER();
    ~AM_BLENDER();
    
    bool set_resolution(size_t width, size_t height);
    bool set_median_combining(bool value);
    bool clear();
    bool add_image(const AM_IMAGE *img);
        
    bool get_xy  (size_t pixel, int *x, int *y) const;
    bool get_rgba(size_t pixel, unsigned char *r, unsigned char *g, unsigned char *b, unsigned char *a) const;            
    size_t pixel_count() const  {return atoms.size();}        
        
    bool is_running() const {return running;}
    bool is_paused()  const {return paused;}        
        
    void stop()   { signal_stop  = true; while(running) std::this_thread::sleep_for(std::chrono::milliseconds(0)); step_thread.join();}
    void start()  { running      = true; step_thread = std::thread(&AM_BLENDER::run, this);}
    void pause()  { signal_pause = true; while(!paused) std::this_thread::sleep_for(std::chrono::milliseconds(0)); signal_pause = false;}
    void resume() { paused       = false;}            
        
    private:
    void run();
    void render();        
    std::vector<AM_ATOM> atoms;
    std::vector<size_t> layers;
    size_t w;
    size_t h;
    bool done;
    bool median_combining;
    
    std::atomic<bool> signal_stop;
    std::atomic<bool> running;
    std::atomic<bool> signal_pause;
    std::atomic<bool> paused;
    std::thread step_thread;         
};

AM_ATOM      am_create_atom(double x, double y, unsigned char r, unsigned char g, unsigned char b, unsigned char a);
double       am_atom_distance(AM_ATOM a1, AM_ATOM a2);
double       am_atom_gradient(AM_ATOM a1, AM_ATOM a2);
const char * am_get_version();
size_t       am_get_warning();
#endif

namespace am {

const unsigned RGB    = 0; // Don't change the color space.
const unsigned HSP    = 1; // Convert from RGB to HSP.
const unsigned NONE   = 2; // Color/motion interpolation is off. 
const unsigned LINEAR = 3; // Uses linear interpolation.
const unsigned SPLINE = 4; // Spline interpolation (motion only).
const unsigned COSINE = 5; // Cosine interpolation (colors only).
const unsigned PERLIN = 6; // Perlin noise dependent interpolation (colors only).
const unsigned STATE_BLOB_DETECTION   = 0;
const unsigned STATE_BLOB_UNIFICATION = 1;
const unsigned STATE_BLOB_MATCHING    = 2;
const unsigned STATE_ATOM_MORPHING    = 3;
const unsigned STATE_DONE             = 4;
const unsigned char HAS_PIXEL         = 1;
const unsigned char HAS_FLUID         = 2;
typedef struct pixel {
    uint16_t x;
    uint16_t y;

    color c;
} pixel;

typedef struct blob { 
    size_t index;             // Position in the nest vector.
    std::set<size_t> surface; // Pixel positions for surface.      
    std::set<size_t> border;  // Pixel positions around surface.
    size_t group  =0;
    bool   unified=false; 
    double x,y,r,g,b,a;
} blob;

typedef struct key_frame {
    std::vector<blob*>      blobs;  // Vector of blobs.
    std::map<size_t, pixel> pixels; // Map of pixels by position.
    std::map<size_t, blob*> owners; // Map of blobs by position.
    double                  x,y,r,g,b,a;
    size_t index          =0; // Position in the nest container.
    size_t first_expansion=0; // Blobs before this index cannot expand.
    size_t first_dust     =0; // First under-sized blob to be unified.
} frame;

typedef union key_point {
    struct {
    uint16_t x;
    uint16_t y;
    uint8_t  x_fract;
    uint8_t  y_fract;
    uint8_t  flags;
    } s;
    uint64_t word;
} point;

typedef struct blob_chain {
    size_t   width =0;
    size_t   height=0;
    point  **points=nullptr;
    size_t **places=nullptr;
    double   energy=0.0;
    size_t   max_surface=0;
    glnemo::CRSpline *splines=nullptr;
#ifdef ATOMORPH_OPENCV
    void            **kdtrees=nullptr;
    void            **feature=nullptr;
#endif
} chain;

const size_t WARN_POINTER_SIZE = 1; // Set when sizeof(void *) is less than 64 bits.
const size_t WARN_PIXEL_SIZE   = 2; // Set when sizeof(pixel) is more than 64 bits.
const size_t WARN_POINT_SIZE   = 3; // Set when sizeof(point) is more than 64 bits.

const unsigned TEXTURE  = 0;
const unsigned AVERAGE  = 1;
const unsigned DISTINCT = 2;

void clear_chain(chain *c);
bool renew_chain(chain *c, size_t width, size_t height);

pixel create_pixel(uint16_t x, uint16_t y, unsigned char r, unsigned char g, unsigned char b, unsigned char a);
pixel create_pixel(uint16_t x, uint16_t y, color c);
inline size_t xy2pos (uint16_t x, uint16_t y) {return (y*(UINT16_MAX+1)+x);}

inline void point2xy(point pt, float *x, float *y) {
    *x = pt.s.x + pt.s.x_fract/float(UINT8_MAX+1);
    *y = pt.s.y + pt.s.y_fract/float(UINT8_MAX+1);
}

inline uint32_t pixel_distance (pixel p1, pixel p2 ) {
    int32_t xd = p1.x-p2.x;
    int32_t yd = p1.y-p2.y;
        
    return xd*xd+yd*yd;
}

inline uint32_t approx_point_distance (point p1, point p2) {
    int32_t xd = p1.s.x-p2.s.x;
    int32_t yd = p1.s.y-p2.s.y;
        
    return xd*xd + yd*yd;
}

inline uint64_t point_distance (point p1, point p2) {
    uint64_t xd = std::abs(((UINT8_MAX+1)*p1.s.x+p1.s.x_fract) - ((UINT8_MAX+1)*p2.s.x+p2.s.x_fract));
    uint64_t yd = std::abs(((UINT8_MAX+1)*p1.s.y+p1.s.y_fract) - ((UINT8_MAX+1)*p2.s.y+p2.s.y_fract));
    
    return (xd*xd + yd*yd);
}

inline double distance(double x1, double y1, double x2, double y2) {
    double xd = x1-x2;
    double yd = y1-y2;
    
    return xd*xd+yd*yd;
}

inline bool point_has_pixel(point p) {
    return (p.s.flags & HAS_PIXEL);
}

inline bool point_has_fluid(point p) {
    return (p.s.flags & HAS_FLUID);
}

const char * get_version();
size_t get_warning();
bool uses_opencv();

}

#include "thread.h"
#include "morph.h"

#endif
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

