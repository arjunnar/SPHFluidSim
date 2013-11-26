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
    // Gradient of Eq(21)
    float rmag = r.abs();
    float constant = -45.0 / (M_PI * pow(h, 6.0));
    Vector3f grad;

    if (0 <= rmag && rmag <= h)
    {
        grad = constant * pow(h - rmag, 2.0) * r.normalized();
    }

    else
    {
        grad = Vector3f::ZERO;
    }

    return grad;
}

float KernelUtilities::polySixKernel(Vector3f r, float h)
{
    // Eq(20)
	float result;
	float rmag = r.abs();

    if (0 <= rmag && rmag <= h)
	{
		float constant = 315.0/(64.0 * M_PI * pow(h, 9.0) );
        float innerPart = pow( (pow(h, 2.0) - r.absSquared()), 3.0 );
		result = constant * innerPart;
	}

	else
	{
		result = 0.0;
	}

	return result;
}

float KernelUtilities::laplacianViscosityKernel(Vector3f r, float h)
{
    // Laplacian of Eq(22); computed for us in the paper
    float rmag = r.abs();
    float result;

    if (0 <= rmag && rmag <= h)
    {
        float constant = 45.0 / (M_PI * pow(h, 6.0));
        return constant * (h - r.abs());
    }

    else
    {
        result = 0.0f;
    }
}
