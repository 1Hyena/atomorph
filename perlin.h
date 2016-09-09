/*
 * See Copyright Notice at the end of this file.
 */
 
class PerlinNoise {
    public:
    PerlinNoise( unsigned seed = 1 );
    
    double noise( double x           ) const { return noise(x,0.0,0.0); }
    double noise( double x, double y ) const { return noise(x,y,0.0);   }
        
    double noise( double x, double y, double z ) const;
    double octaveNoise( double x, int octaves ) const;
    double octaveNoise( double x, double y, int octaves ) const;
    double octaveNoise( double x, double y, double z, int octaves ) const;

    private:   
    double fade( double t                     ) const { return t*t*t*(t*(t*6-15)+10); }
    double lerp( double t, double a, double b ) const { return a + t * (b - a);       }    
    
    double grad( int hash, double x, double y, double z ) const;
    int p[512];
};

/*
The MIT License (MIT)

Copyright (c) 2013 Reputeless

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

