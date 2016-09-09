/*
 * See Copyright Notice in fluidmodel.h
 * See Copyright Notice in atomorph.h
 *
 */
 
#include <cmath>
#include <algorithm>

#include "fluidmodel.h"

Material::Material(double m_,
             double rd_,
             double k_,
             double v_,
             double d_,
             double g_)
    :m(m_),
     rd(rd_),
     k(k_),
     v(v_),
     d(d_),
     g(g_) {

}


Particle::Particle() {

}

Particle::Particle(Material *mat_,
             double x_,
             double y_,
             double u_,
             double v_)
    :mat(mat_),
     x(x_),
     y(y_),
     u(u_),
     v(v_),
     dudx(0.0),
     dudy(0.0),
     dvdx(0.0),
     dvdy(0.0),
     cx(0.0),
     cy(0.0) {

    std::fill(px, px + 3, 0.0);
    std::fill(gx, gx + 3, 0.0);
    std::fill(py, py + 3, 0.0);
    std::fill(gy, gy + 3, 0.0);

}


Node::Node()
    :m(0.0),
     d(0.0),
     gx(0.0),
     gy(0.0),
     u(0.0),
     v(0.0),
     ax(0.0),
     ay(0.0),
     active(false) {

}

void Node::clear() {
    m = d = gx = gy = u = v = ax = ay = 0.0;
    r = 0.0;
    g = 0.0;
    b = 0.0;
    a = 0.0;
    weight = 0.0;
    active = false;
}

void Particle::clear() {
    u = 0.0;
    v = 0.0;
    dudx = 0.0;
    dudy = 0.0;
    dvdx = 0.0;
    dvdy = 0.0;
    cx   = 0.0;
    cy   = 0.0;
    strength = 1.0;
    std::fill(px, px + 3, 0.0);
    std::fill(gx, gx + 3, 0.0);
    std::fill(py, py + 3, 0.0);
    std::fill(gy, gy + 3, 0.0);
    source_owner = false;
}

void Particle::copy_from(Particle *p) {
    x = p->x;
    y = p->y;
    gravity_x = p->gravity_x;
    gravity_y = p->gravity_y;
    freedom_r = p->freedom_r;
    R = p->R; G = p->G; B = p->B; A = p->A;
    r = p->r; g = p->g; b = p->b; a = p->a;    
}

FluidModel::FluidModel(unsigned gsizeX_,
                       unsigned gsizeY_,
                       unsigned particle_count_)
    :particles(NULL),
     gsizeX(gsizeX_),
     gsizeY(gsizeY_),
     particle_count(particle_count_),
     active(NULL),
     activeCount(0),
     water(1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
     pressed(false),
     pressedprev(false),
     mx(0.0),
     my(0.0),
     mxprev(0.0),
     myprev(0.0),
     grid(NULL),
     lines(NULL) {

    grid = new Node *[gsizeX];
    for (unsigned i = 0; i < gsizeX; ++i) {
        grid[i] = new Node[gsizeY];
    }
    for (size_t j=0; j<gsizeY; ++j) {
        for (size_t i=0; i<gsizeX; ++i) {
            grid[i][j].pos = (j * gsizeX)+i;
        }
    }

    active = new Node*[gsizeX * gsizeY];

    particles = new Particle[particle_count];

    for (size_t i=0; i< particle_count; i++) {
        particles[i] = Particle(&water, gsizeX/2,
                                        gsizeY/2,
                                        0.0,
                                        0.0);
    }

    lines = new Line[particle_count];
}


FluidModel::~FluidModel() {
    for (unsigned i = 0; i < gsizeX; ++i) {
        delete [] grid[i];
    }
    delete [] grid;

    delete [] active;

    delete [] particles;

    delete [] lines;
}


void FluidModel::step(size_t steps_left, double freedom_radius, double morph_time) {
    bool drag = false;
    double mdx = 0.0, mdy = 0.0;

    if (pressed and pressedprev) {
        drag = true;
        mdx = (mx - mxprev);
        mdy = (my - myprev);
    }

    pressedprev = pressed;
    mxprev = mx;
    myprev = my;

    for (unsigned i = 0; i < activeCount; ++i) {
        active[i]->clear();
    }
    activeCount = 0;

    double phi;
    double x, y;
    double fx = 0.0, fy = 0.0;

    Particle *p;
    Node *n;

    for (unsigned i = 0; i < particle_count; ++i) {
        p = particles + i;
        if (p->active == false) continue; // CUSTOMIZATION
        
        p->cx = static_cast<int>(p->x - 0.5);
        p->cy = static_cast<int>(p->y - 0.5);

        x = p->cx - p->x;

        p->px[0] = (0.5 * x * x + 1.5 * x + 1.125);
        p->gx[0] = (x + 1.5);
        x += 1.0;
        p->px[1] = (-x * x + 0.75);
        p->gx[1] = (-2.0 * x);
        x += 1.0;
        p->px[2] = (0.5 * x * x - 1.5 * x + 1.125);
        p->gx[2] = (x - 1.5);

        y = p->cy - p->y;

        p->py[0] = (0.5 * y * y + 1.5 * y + 1.125);
        p->gy[0] = (y + 1.5);
        y += 1.0;
        p->py[1] = (-y * y + 0.75);
        p->gy[1] = (-2.0 * y);
        y += 1.0;
        p->py[2] = (0.5 * y * y - 1.5 * y + 1.125);
        p->gy[2] = (y - 1.5);
        
        for (unsigned i = 0; i < 3; ++i) {
            for (unsigned j = 0; j < 3; ++j) {
                n = grid[p->cx + i] + p->cy + j;
                if (not n->active) {
                    active[activeCount++] = n;
                    n->active = true;
                }

                phi = p->px[i] * p->py[j];
                n->m += phi * p->mat->m;
                n->d += phi;
                n->gx += p->gx[i] * p->py[j];
                n->gy += p->px[i] * p->gy[j];

                // CUSTOM:
                if (p->mature && p->strength > 0.0) {
                    double new_w = p->strength;
                    double old_w = n->weight;
                    double sum_w = new_w + old_w;
                    n->r = (old_w * n->r + new_w * p->r)/sum_w;
                    n->g = (old_w * n->g + new_w * p->g)/sum_w;
                    n->b = (old_w * n->b + new_w * p->b)/sum_w;
                    n->a = (old_w * n->a + new_w * p->a)/sum_w;
                    n->weight = sum_w;
                }
                // END OF CUSTOM
            }
        }

    }

    unsigned cx, cy;
    unsigned cxi, cyi;

    double pdx, pdy;
    double C20, C02, C30, C03;
    double C21, C31, C12, C13, C11;

    double csum1, csum2;

    double u, u2, u3;
    double v, v2, v3;

    double vx, vy;

    double density;
    double pressure;
    double weight;

    Node *n01, *n02, *n11, *n12;

    for (unsigned i = 0; i < particle_count; ++i) {
        p = particles + i;
        if (p->active == false) continue; // CUSTOMIZATION
        
        cx = static_cast<int>(p->x);
        cy = static_cast<int>(p->y);
        cxi = cx + 1;
        cyi = cy + 1;

        n01 = grid[cx] + cy;
        n02 = grid[cx] + cyi;
        n11 = grid[cxi] + cy;
        n12 = grid[cxi] + cyi;

        pdx = n11->d - n01->d;
        pdy = n02->d - n01->d;
        C20 = 3.0 * pdx - n11->gx - 2.0 * n01->gx;
        C02 = 3.0 * pdy - n02->gy - 2.0 * n01->gy;
        C30 = -2.0 * pdx + n11->gx + n01->gx;
        C03 = -2.0 * pdy + n02->gy + n01->gy;

        csum1 = n01->d + n01->gy + C02 + C03;
        csum2 = n01->d + n01->gx + C20 + C30;

        C21 = 3.0 * n12->d - 2.0 * n02->gx - n12->gx -
                     3.0 * csum1 - C20;
        C31 = -2.0 * n12->d + n02->gx + n12->gx +
                     2.0 * csum1 - C30;
        C12 = 3.0 * n12->d - 2.0 * n11->gy - n12->gy -
                     3.0 * csum2 - C02;
        C13 = -2.0 * n12->d + n11->gy + n12->gy +
                     2.0 * csum2 - C03;
        C11 = n02->gx - C13 - C12 - n01->gx;

        u = p->x - cx;
        u2 = u * u;
        u3 = u * u2;
        v = p->y - cy;
        v2 = v * v;
        v3 = v * v2;

        density = n01->d + n01->gx * u +
                  n01->gy * v + C20 * u2 +
                  C02 * v2 + C30 * u3 +
                  C03 * v3 + C21 * u2 * v +
                  C31 * u3 * v + C12 * u * v2 +
                  C13 * u * v3 + C11 * u * v;

        pressure = density - 1.0;
        if (pressure > 2.0) {
            pressure = 2.0;
        }

        fx = fy = 0.0;

        if (p->x < 4.0) {
            fx += p->mat->m * (4.0 - p->x);
        } else if (p->x > gsizeX - 5) {
            fx += p->mat->m * (gsizeX - 5 - p->x);
        }

        if (p->y < 4.0) {
            fy += p->mat->m * (4.0 - p->y);
        } else if (p->y > gsizeY - 5) {
            fy += p->mat->m * (gsizeY - 5 - p->y);
        }

        if (drag) {
            vx = fabs(p->x - mx);
            vy = fabs(p->y - my);
            if (vx < 10.0 and vy < 10.0) {
                weight = p->mat->m * (1.0 - vx * 0.1) *
                         (1.0 - vy * 0.1);
                fx += weight * (mdx - p->u);
                fy += weight * (mdy - p->v);
            }
        }

        for (unsigned i = 0; i < 3; ++i) {
            for (unsigned j = 0; j < 3; ++j) {
                n = grid[p->cx + i] + p->cy + j;
                phi = p->px[i] * p->py[j];
                n->ax += -((p->gx[i] * p->py[j]) * pressure) + fx * phi;
                n->ay += -((p->px[i] * p->gy[j]) * pressure) + fy * phi;
            }
        }
    }

    for (unsigned i = 0; i < activeCount; ++i) {
        n = active[i];
        if (n->m > 0.0) {
            n->ax /= n->m;
            n->ay /= n->m;

            //ORIGINALLY GRAVITY WAS INDICATED BY THESE:
            //n->ay += 0.03;
            //n->ax += 0.03;
        }
    }

    double mu, mv;

    for (unsigned pi = 0; pi < particle_count; ++pi) {
        p = particles + pi;
        if (p->active == false) continue; // CUSTOMIZATION

        for (unsigned i = 0; i < 3; ++i) {
            for (unsigned j = 0; j < 3; ++j) {
                n = grid[p->cx + i] + p->cy + j;
                phi = p->px[i] * p->py[j];

                /* CUSTOM STUFF */
                double custom_ax=0.0;
                double custom_ay=0.0;
                {
                    double x1 = p->x;
                    double y1 = p->y;
                    double x2 = p->gravity_x;
                    double y2 = p->gravity_y;
                    double a = 0.03;                   

                    // Gravitate towards (x2, y2).
                    double A = std::abs(y1 - y2);
                    double B = std::abs(x1 - x2);
                    double C = std::sqrt(A*A + B*B);
                    
                    if (a >= C) a = C;
                    
                    if (B <= 0.0) {
                        if (y2 <= y1) custom_ay -= a;
                        else          custom_ay += a;
                    }
                    else if (C > 0.0) {
                        double dx = (a * B) / C;
                        double dy = (A *dx) / B;
                        
                        if (x1 <= x2) custom_ax += dx;
                        else          custom_ax -= dx;

                        if (y1 <= y2) custom_ay += dy;
                        else          custom_ay -= dy;
                    }
                }
                /* END OF CUSTOM STUFF */
                
                p->u += phi * (n->ax+custom_ax);
                p->v += phi * (n->ay+custom_ay);

            }
        }
        mu = p->mat->m * p->u;
        mv = p->mat->m * p->v;
        if (p->mature == false) {
            mu *= 0.0;
            mv *= 0.0;
        }
        for (unsigned i = 0; i < 3; ++i) {
            for (unsigned j = 0; j < 3; ++j) {
                n = grid[p->cx + i] + p->cy + j;
                phi = p->px[i] * p->py[j];
                n->u += phi * mu;
                n->v += phi * mv;
            }
        }
    }

    for (unsigned i = 0; i < activeCount; ++i) {
        n = active[i];
        if (n->m > 0.0) {
            n->u /= n->m;
            n->v /= n->m;
        }
    }

    double gu, gv;

    for (unsigned i = 0; i < particle_count; ++i) {
        p = particles + i;
        
        // CUSTOMIZATION:
        if (p->active == false) continue; 
        double nR = 0.0, nG = 0.0, nB = 0.0, nA = 0.0;
        double weight = 0.0;
        // END OF CUSTOMIZATION
        
        gu = 0.0;
        gv = 0.0;
        for (unsigned ii = 0; ii < 3; ++ii) {
            for (unsigned j = 0; j < 3; ++j) {
                n = grid[p->cx + ii] + p->cy + j;
                phi = p->px[ii] * p->py[j];
                gu += phi * n->u;
                gv += phi * n->v;
                
                if (n->weight > 0.0) {
                    weight += n->weight;
                    nR     += n->r * n->weight;
                    nG     += n->g * n->weight;
                    nB     += n->b * n->weight;
                    nA     += n->a * n->weight;
                }
            }
        }
        
        // CUSTOMIZATION:
        if (weight > 0.0) {
            double wr = 1.0;
            double wg = 1.0;
            double wb = 1.0;
            double wa = 1.0;
            nR /= weight; wr = std::abs(nR - p->R); 
            nG /= weight; wg = std::abs(nG - p->G); 
            nB /= weight; wb = std::abs(nB - p->B); 
            nA /= weight; wa = std::abs(nA - p->A); 
            
            if (p->mature == false) {
                p->r = nR;
                p->g = nG;
                p->b = nB;
                p->a = nA;
            }
            else {               
                p->r = (1.0 - wr) * p->r + wr * nR;
                p->g = (1.0 - wg) * p->g + wg * nG; 
                p->b = (1.0 - wb) * p->b + wb * nB; 
                p->a = (1.0 - wa) * p->a + wa * nA;
            }            
        }
        // END OF CUSTOMIZATION
        
        p->x += gu;
        p->y += gv;

        /* CUSTOM STUFF */
        {
            double x1 = p->x;
            double y1 = p->y;
            double x2 = p->gravity_x;
            double y2 = p->gravity_y;

            // Gravitate towards (x2, y2) if the distance
            // is greater than freedom radius.
            double A = std::abs(y1 - y2);
            double B = std::abs(x1 - x2);
            double C = std::sqrt(A*A + B*B);            

            double r = freedom_radius * p->freedom_r;
            if (C > r) {
                double mx = 0.0;
                double my = 0.0;
                
                double a = C - r;
                if (a >= C) a = C;
                
                if (B <= 0.0) {
                    if (y2 <= y1) my -= a;
                    else          my += a;
                }
                else if (C > 0.0) {
                    double dx = (a * B) / C;
                    double dy = (A *dx) / B;
                    
                    if (x1 <= x2) mx += dx;
                    else          mx -= dx;

                    if (y1 <= y2) my += dy;
                    else          my -= dy;
                }
            
                double w = 1.0 / (steps_left + 1);                                
                
                mx *= w;
                my *= w;
                
                p->x += mx;
                p->y += my;                
            }
        }
        /* END OF CUSTOM STUFF */

        p->u += gu - p->u;
        p->v += gv - p->v;
        if (p->x < 1.0) {
            p->x = 1.0 + static_cast<double>(rand()) / RAND_MAX * 0.01;
            p->u = 0.0;
        } else if (p->x > gsizeX - 2) {
            p->x = gsizeX - 2 - static_cast<double>(rand()) / RAND_MAX * 0.01;
            p->u = 0.0;
        }
        if (p->y < 1.0) {
            p->y = 1.0 + static_cast<double>(rand()) / RAND_MAX * 0.01;
            p->v = 0.0;
        } else if (p->y > gsizeY - 2) {
            p->y = gsizeY - 2 - static_cast<double>(rand()) / RAND_MAX * 0.01;
            p->v = 0.0;
        }
    }

    // ORIGINALLY THIS PIECE OF CODE WAS HERE
    // BUT FOR FLUID MORPHING IT IS NOT NEEDED:
    /*for (unsigned i = 0; i < particle_count; ++i) {
        p = particles + i;
        if (p->inactive) continue; // CUSTOMIZATION
                
        lines[i].x1 = p->x;
        lines[i].y1 = p->y;
        lines[i].x2 = p->x - p->u;
        lines[i].y2 = p->y - p->v;
    }*/
}

void FluidModel::setPressed(bool b) {
    pressed = b;
}

void FluidModel::setMovePos(double x, double y) {
    mx = x;
    my = y;
}

const Line * FluidModel::getLines() {
    return lines;
}

Particle * FluidModel::getParticles() {
    return particles;
}

Material * FluidModel::getMaterial() {
    return &water;
}

