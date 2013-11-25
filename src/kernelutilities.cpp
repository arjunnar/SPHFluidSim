#include "kernelutilities.h"
#define _USE_MATH_DEFINES

float KernelUtilities::spikyKernel(Vector3f r, float h)
{
    float rmag = r.abs();
    float constant = 15.0 / (M_PI * pow(h, 6.0));
    float result;

    if (0 <= rmag && rmag <= h)
    {
        result = constant * pow(h - rmag, 3.0);
    }

    else
    {
        result = 0.0;
    }

    return result;
}

Vector3f KernelUtilities::gradSpikyKernel(Vector3f r, float h)
{
    float rmag = r.abs();
    float constant = -45.0 / (M_PI * pow(h, 6.0));
    Vector3f grad;

    if (0 <= rmag && rmag <= h)
    {
        grad = constant * pow(h - rmag, 2.0) / rmag * r;
    }

    else
    {
        grad = Vector3f::ZERO;
    }

    return grad;
}