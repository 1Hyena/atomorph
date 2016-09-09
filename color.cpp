/*
 * See Copyright Notice in atomorph.h
 */
 
#include <math.h>
#include <stdlib.h>
#include <algorithm>

#include "color.h"

namespace am {

color create_color(unsigned char r, unsigned char g, unsigned char b, unsigned char a) {   
    color c;
    c.r = r;
    c.g = g;
    c.b = b;
    c.a = a;
    return c;
}

color create_color(double r, double g, double b, double a) {
    color c;
    c.r = std::round(r*255.0);
    c.g = std::round(g*255.0);
    c.b = std::round(b*255.0);
    c.a = std::round(a*255.0);
    return c;
}

color rgb_to_hsp(color c) {
    color hsp = c;    
    
    double r,g,b,h,s,p;
	r = c.r/255.0; 
	g = c.g/255.0; 
	b = c.b/255.0; 
	RGBtoHSP(r,g,b,&h,&s,&p);
	hsp.r = std::round(h*255.0);
	hsp.g = std::round(s*255.0);
	hsp.b = std::round(p*255.0);
    
    return hsp;
}

color hsp_to_rgb(color c) {
    color rgb = c;

    double r,g,b,h,s,p;
	h = c.r/255.0; 
	s = c.g/255.0; 
	p = c.b/255.0; 
	HSPtoRGB(h,s,p,&r,&g,&b);
	rgb.r = std::min(std::round(r*255.0), 255.0);
	rgb.g = std::min(std::round(g*255.0), 255.0);
	rgb.b = std::min(std::round(b*255.0), 255.0);    
    
    return rgb;
}

const double Pr = 0.299;
const double Pg = 0.587;
const double Pb = 0.114;

//  public domain function by Darel Rex Finley, 2006
//
//  This function expects the passed-in values to be on a scale
//  of 0 to 1, and uses that same scale for the return values.
//
//  See description/examples at alienryderflex.com/hsp.html
void RGBtoHSP(double  R, double  G, double  B, double *H, double *S, double *P) {
    //  Calculate the Perceived brightness.
    *P=sqrt(R*R*Pr+G*G*Pg+B*B*Pb);
    
    //  Calculate the Hue and Saturation.  (This part works
    //  the same way as in the HSV/B and HSL systems???.)
    if (R==G && R==B) {
        *H=0.; 
        *S=0.; 
        return; 
    }
    
    if (R>=G && R>=B) {// R is largest
        if (B>=G) {
            *H=6./6.-1./6.*(B-G)/(R-G); 
            *S=1.-G/R; 
        }
        else {
            *H=0./6.+1./6.*(G-B)/(R-B); 
            *S=1.-B/R; 
        }
    }
    else if (G>=R && G>=B) {// G is largest
        if (R>=B) {
            *H=2./6.-1./6.*(R-B)/(G-B); *S=1.-B/G; 
        }
        else {
            *H=2./6.+1./6.*(B-R)/(G-R); *S=1.-R/G; 
        }
    }
    else {// B is largest
        if (G>=R) {
            *H=4./6.-1./6.*(G-R)/(B-R); *S=1.-R/B; 
        }
        else {
            *H=4./6.+1./6.*(R-G)/(B-G); *S=1.-G/B;
        }
    }
}



//  public domain function by Darel Rex Finley, 2006
//
//  This function expects the passed-in values to be on a scale
//  of 0 to 1, and uses that same scale for the return values.
//
//  Note that some combinations of HSP, even if in the scale
//  0-1, may return RGB values that exceed a value of 1.  For
//  example, if you pass in the HSP color 0,1,1, the result
//  will be the RGB color 2.037,0,0.
//
//  See description/examples at alienryderflex.com/hsp.html
void HSPtoRGB(double  H, double  S, double  P, double *R, double *G, double *B) {
    double  part, minOverMax=1.-S ;
    
    if (minOverMax>0.) {
        if ( H<1./6.) {   //  R>G>B
            H= 6.*( H-0./6.); part=1.+H*(1./minOverMax-1.);
            *B=P/sqrt(Pr/minOverMax/minOverMax+Pg*part*part+Pb);
            *R=(*B)/minOverMax; *G=(*B)+H*((*R)-(*B));
        }
        else if ( H<2./6.) {   //  G>R>B
            H= 6.*(-H+2./6.); part=1.+H*(1./minOverMax-1.);
            *B=P/sqrt(Pg/minOverMax/minOverMax+Pr*part*part+Pb);
            *G=(*B)/minOverMax; *R=(*B)+H*((*G)-(*B));
        }
        else if ( H<3./6.) {   //  G>B>R
            H= 6.*( H-2./6.); part=1.+H*(1./minOverMax-1.);
            *R=P/sqrt(Pg/minOverMax/minOverMax+Pb*part*part+Pr);
            *G=(*R)/minOverMax; *B=(*R)+H*((*G)-(*R));
        }
        else if ( H<4./6.) {   //  B>G>R
            H= 6.*(-H+4./6.); part=1.+H*(1./minOverMax-1.);
            *R=P/sqrt(Pb/minOverMax/minOverMax+Pg*part*part+Pr);
            *B=(*R)/minOverMax; *G=(*R)+H*((*B)-(*R));
        }
        else if ( H<5./6.) {   //  B>R>G
            H= 6.*( H-4./6.); part=1.+H*(1./minOverMax-1.);
            *G=P/sqrt(Pb/minOverMax/minOverMax+Pr*part*part+Pg);
            *B=(*G)/minOverMax; *R=(*G)+H*((*B)-(*G)); }
        else {   //  R>B>G
            H= 6.*(-H+6./6.); part=1.+H*(1./minOverMax-1.);
            *G=P/sqrt(Pr/minOverMax/minOverMax+Pb*part*part+Pg);
            *R=(*G)/minOverMax; *B=(*G)+H*((*R)-(*G)); 
        }
    }
    else {
        if ( H<1./6.) {   //  R>G>B
            H= 6.*( H-0./6.); *R=sqrt(P*P/(Pr+Pg*H*H)); *G=(*R)*H; *B=0.; 
        }
        else if ( H<2./6.) {   //  G>R>B
            H= 6.*(-H+2./6.); *G=sqrt(P*P/(Pg+Pr*H*H)); *R=(*G)*H; *B=0.; 
        }
        else if ( H<3./6.) {   //  G>B>R
            H= 6.*( H-2./6.); *G=sqrt(P*P/(Pg+Pb*H*H)); *B=(*G)*H; *R=0.; 
        }
        else if ( H<4./6.) {   //  B>G>R
            H= 6.*(-H+4./6.); *B=sqrt(P*P/(Pb+Pg*H*H)); *G=(*B)*H; *R=0.; 
        }
        else if ( H<5./6.) {   //  B>R>G
            H= 6.*( H-4./6.); *B=sqrt(P*P/(Pb+Pr*H*H)); *R=(*B)*H; *G=0.; 
        }
        else {   //  R>B>G
            H= 6.*(-H+6./6.); *R=sqrt(P*P/(Pr+Pb*H*H)); *B=(*R)*H; *G=0.;
        }
    }
}      

}      

