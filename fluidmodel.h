/*
 * See Copyright Notice at the end of this file.
 */
 
#ifndef FLUIDMODEL_H
#define FLUIDMODEL_H

#include <algorithm>

struct Line {
    double x1, y1;
    double x2, y2;
};

struct Material {

    Material(double m_,
             double rd_,
             double k_,
             double v_,
             double d_,
             double g_);

    double m;  //mass
    double rd;
    double k;
    double v;
    double d;
    double g;
};


struct Particle {

    Particle();
    Particle(Material *mat_,
             double x_,
             double y_,
             double u_,
             double v_);

    void clear();
    void copy_from(Particle *p);

    Material *mat;

    double x;
    double y;
    double u;
    double v;

    double dudx;
    double dudy;
    double dvdx;
    double dvdy;
    unsigned cx;
    unsigned cy;

    double px[3];
    double py[3];
    double gx[3];
    double gy[3];

    // CUSTOM VARIABLES:
    double gravity_x;
    double gravity_y;
    double freedom_r; // freedom radius multiplier
    bool   active;
    bool   mature;
    double R,G,B,A; // ideal color   (does not blend)
    double r,g,b,a; // current color (blends over time)
    size_t destination_pos;
    size_t source_pos;
    size_t frame_key;
    double strength; // how easily it regains its color
    bool   hack;
    bool   source_owner; // true when it was first one to occupy the source location.
    // END OF CUSTOM VARIABLES
};


struct Node {

    Node();

    void clear();

    double m;
    double d;
    double gx;
    double gy;
    double u;
    double v;
    double ax;
    double ay;
    bool   active;
    size_t pos;
    double r,g,b,a;
    double weight;
};


class FluidModel {

    public:

        FluidModel(unsigned gsizeX_,
                   unsigned gsizeY_,
                   unsigned particle_count);
        ~FluidModel();

        void setPressed(bool b);
        void setMovePos(double x, double y);

        void step(size_t steps_left, double freedom_radius, double t);
        const Line * getLines();
        Material * getMaterial();
        Particle * getParticles();

        unsigned get_particle_count(void) {return particle_count;}
    private:

        Particle *   particles;
        unsigned     gsizeX;
        unsigned     gsizeY;
        unsigned     particle_count;
        Node **      active;
        unsigned     activeCount;
        Material     water;

        bool         pressed;
        bool         pressedprev;

        double       mx;
        double       my;
        double       mxprev;
        double       myprev;

        Node **      grid;

        Line *       lines;
};

#endif

/* This version:
 * Copyright Xueqiao Xu ( http://code.google.com/p/mycodeplayground )
 * MIT License ( http://www.opensource.org/licenses/mit-license.php )
 * Download from: http://code.google.com/p/mycodeplayground

 * Python version:
 * Copyright Xueqiao Xu ( http://code.google.com/p/mycodeplayground )
 * MIT License ( http://www.opensource.org/licenses/mit-license.php )
 * Download from: http://code.google.com/p/mycodeplayground

 * Javascript version:
 * Copyright Stephen Sinclair (radarsat1) ( http://www.music.mcgill.ca/~sinclair )
 * MIT License ( http://www.opensource.org/licenses/mit-license.php )
 * Download from: http://www.music.mcgill.ca/~sinclair/blog

 * Flash version:
 * Copyright iunpin ( http://wonderfl.net/user/iunpin )
 * MIT License ( http://www.opensource.org/licenses/mit-license.php )
 * Download from: http://wonderfl.net/c/6eu4

 * Original Java version:
 * http://grantkot.com/MPM/Liquid.html
 */

