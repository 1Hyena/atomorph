/*
 * See Copyright Notice in atomorph.h
 */
 
namespace am {

typedef struct color {
    uint8_t r;
    uint8_t g;
    uint8_t b;
    uint8_t a;
} color;

color create_color(unsigned char r, unsigned char g, unsigned char b, unsigned char a);
color create_color(double r, double g, double b, double a);

inline double color_distance(color c1, color c2) {
    int16_t rd = c1.r-c2.r;
    int16_t gd = c1.g-c2.g;
    int16_t bd = c1.b-c2.b;
    int16_t ad = c1.a-c2.a;
        
    return sqrt(rd*rd+gd*gd+bd*bd+ad*ad)/510.0;
}

color  rgb_to_hsp(color c);
color  hsp_to_rgb(color c);

void RGBtoHSP(double  R, double  G, double  B, double *H, double *S, double *P);
void HSPtoRGB(double  H, double  S, double  P, double *R, double *G, double *B);

}

