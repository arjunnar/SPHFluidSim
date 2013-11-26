#ifndef SPHUTILITIES_H
#define SPHUTILITIES_H

#include "vecmath.h"

class KernelUtilities
{
public:
    static float spikyKernel(Vector3f r, float h);
    static Vector3f gradSpikyKernel(Vector3f r, float h);

    static float polySixKernel(Vector3f r, float h);
    //static Vector3f gradPolySixKernel(Vector3f r, float h);

    //static float viscosityKernel(Vector3f r, float h);
    static float laplacianViscosityKernel(Vector3f r, float h);
};

#endif // SPHUTILITIES_H
