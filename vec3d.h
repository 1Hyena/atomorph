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
#ifndef GLNEMOVEC3D_H
#define GLNEMOVEC3D_H
#include <math.h>
#include <iostream>
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/
namespace glnemo {

class Vec3D{
    public:
    ~Vec3D() {};

	double x, y, z;
    Vec3D( double InX, double InY, double InZ ) : x( InX ), y( InY ), z( InZ ) {}
    Vec3D( const Vec3D& V) : x( V.x ), y( V.y ), z( V.z )                      {}
	Vec3D( ) : x(0), y(0), z(0)                                                {}
    
    inline void set( const double InX, const double InY, const double InZ ) {
        x = InX; y = InY; z = InZ;
    } 

    inline bool   operator== (const Vec3D& V2) const { return (x == V2.x && y == V2.y && z == V2.z); }
    inline Vec3D  operator+  (const Vec3D& V2) const { return Vec3D( x + V2.x,  y + V2.y,  z + V2.z);}
    inline Vec3D  operator-  (const Vec3D& V2) const { return Vec3D( x - V2.x,  y - V2.y,  z - V2.z);}
    inline Vec3D  operator-  (               ) const { return Vec3D(-x, -y, -z);                     }
    inline Vec3D  operator/  (const Vec3D& V2) const { return Vec3D (x / V2.x,  y / V2.y,  z / V2.z);}
    inline Vec3D  operator*  (const Vec3D& V2) const { return Vec3D (x * V2.x,  y * V2.y,  z * V2.z);}
    inline Vec3D  operator*  (double S       ) const { return Vec3D (x * S,  y * S,  z * S);         }
    inline Vec3D  operator/  (double S       ) const { double f=1.0/S; return Vec3D(x*f,y*f,z*f);    }
    inline double operator[] (int i          )       { return (i == 0 ? x : (i == 1 ? y : z));       }
    inline Vec3D& operator=  (const Vec3D& V2)       { x=V2.x; y=V2.y; z=V2.z; return *this;         }
    inline void   operator+= (const Vec3D& V2)       { x += V2.x; y += V2.y; z += V2.z;              }
    inline void   operator-= (const Vec3D& V2)       { x -= V2.x; y -= V2.y; z -= V2.z;              }

    inline double Dot( const Vec3D &V1 ) const {
        return V1.x*x + V1.y*y + V1.z*z;
    }

    inline Vec3D CrossProduct( const Vec3D &V2 ) const {
        return Vec3D( y * V2.z  -  z * V2.y,
                      z * V2.x  -  x * V2.z,
                      x * V2.y  -  y * V2.x );
    }

    Vec3D RotByMatrix( const double m[16] ) const {
        return Vec3D( x*m[0] + y*m[4] + z*m[8],
                      x*m[1] + y*m[5] + z*m[9],
                      x*m[2] + y*m[6] + z*m[10] );
    }

	// These require math.h for the sqrtf function
    inline double Magnitude( ) const {
        return sqrt( x*x + y*y + z*z );
    }

    inline double Distance( const Vec3D &V1 ) const {
        return ( *this - V1 ).Magnitude();
    }

    inline void Normalize() {
        double fMag = ( x*x + y*y + z*z );
        if (fMag == 0) return;
        
        double fMult = 1.0/sqrtf(fMag);
        x *= fMult;
        y *= fMult;
        z *= fMult;
        return;
    }
};

}

#endif

