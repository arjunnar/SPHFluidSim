#ifndef SPHUTILITIES_H
#define SPHUTILITIES_H

#include "vecmath.h"

    static const float h = 0.04;

    static const float SPIKY_KERNEL_GRAD_CONSTANT = -45.0 / (M_PI * pow(h, 6.0));

    static const float LAPLACIAN_VISCOSITY_KERNEL_CONSTANT = 45.0 / (M_PI * pow(h, 6.0));

    static const float POLY_SIX_KERNEL_CONSTANT = 315.0/(64.0 * M_PI * pow(h, 9.0) );
    static const float GRAD_POLY_SIX_KERNEL_CONSTANT =  -945.0 / (32.0 * M_PI * pow(h, 9.0));
    static const float LAPLACIAN_POLY_SIX_KERNEL_CONSTANT = 945.0 / (8 * M_PI * pow(h, 9.0));


class KernelUtilities
{
public:
    static Vector3f gradSpikyKernel(Vector3f r);

    static float polySixKernel(Vector3f r);
    static Vector3f gradPolySixKernel(Vector3f r);
    static float laplacianPolySixKernel(Vector3f r);

    //static float viscosityKernel(Vector3f r, float h);
    static float laplacianViscosityKernel(Vector3f r);





};

#endif // SPHUTILITIES_H
