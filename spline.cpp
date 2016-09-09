// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2013
// e-mail:   Jean-Charles.Lambert@oamp.fr
// address:  Dynamique des galaxies
//           Laboratoire d'Astrophysique de Marseille
//           Pôle de l'Etoile, site de Château-Gombert
//           38, rue Frédéric Joliot-Curie
//           13388 Marseille cedex 13 France
//           CNRS U.M.R 7326
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".
// ============================================================================
#include "spline.h"

namespace glnemo {

CRSpline::CRSpline() : vp(), delta_t(0.0) {}

CRSpline::CRSpline(const CRSpline& s) {
    for (int i = 0; i < (int)s.vp.size(); i++)
        vp.push_back(s.vp[i]);
    
    delta_t = s.delta_t;
}

CRSpline::~CRSpline() {}

// Solve the Catmull-Rom parametric equation for a given time(t) and vector quadruple (p1,p2,p3,p4)
Vec3D CRSpline::Eq(double t, const Vec3D& p1, const Vec3D& p2, const Vec3D& p3, const Vec3D& p4) {
    double t2 = t  * t;
    double t3 = t2 * t;
    double b1 = 0.5 * (    -t3 + 2.0*t2 -   t);
    double b2 = 0.5 * ( 3.0*t3 - 5.0*t2 + 2.0);
    double b3 = 0.5 * (-3.0*t3 + 4.0*t2 +   t);
    double b4 = 0.5 * (     t3 -     t2      );
    
    return (p1*b1 + p2*b2 + p3*b3 + p4*b4);
}

void CRSpline::AddSplinePoint(const Vec3D& v) {
    vp.push_back(v);
    delta_t = 1.0 / vp.size();
}

Vec3D CRSpline::GetInterpolatedSplinePoint(double t) {
    // Find out in which interval we are on the spline
    int p = (int)(t / delta_t);
    // Compute local control point indices
    int p0 = p - 1; p0 = (p0 < 0 ? vp.size()-1 : (p0 >= (int)vp.size() ? p0 - (int)vp.size() : p0));
    int p1 = p;     p1 = (p1 < 0 ? vp.size()-1 : (p1 >= (int)vp.size() ? p1 - (int)vp.size() : p1));
    int p2 = p + 1; p2 = (p2 < 0 ? vp.size()-1 : (p2 >= (int)vp.size() ? p2 - (int)vp.size() : p2));
    int p3 = p + 2; p3 = (p3 < 0 ? vp.size()-1 : (p3 >= (int)vp.size() ? p3 - (int)vp.size() : p3));
    // Relative (local) time
    double lt = (t - delta_t*(double)p) / delta_t;
    // Interpolate
    return CRSpline::Eq(lt, vp[p0], vp[p1], vp[p2], vp[p3]);
}

int CRSpline::GetNumPoints() {
    return vp.size();
}

Vec3D& CRSpline::GetNthPoint(int n) {
    return vp[n];
}

}

