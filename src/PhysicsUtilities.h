/*
 * ForceUtilities.h
 *
 *  Created on: Oct 19, 2013
 *      Author: arjunnar
 */

#ifndef FORCEUTILITIES_H_
#define FORCEUTILITIES_H_
#include <vector>
#include <vecmath.h>

using namespace std;

class PhysicsUtilities {
public:
    static Vector3f getSpringForce(Vector3f xi, Vector3f xj,float springConstant, float springRestLength);
    static Vector3f getGravityForce(float mass, float gravityConstant);
    static Vector3f getDragForce(float dragConstant, Vector3f velocityVector);
    static Vector3f getPositionOfParticle(vector<Vector3f> &state, int particleNum);
    static Vector3f getVelocityOfParticle(vector<Vector3f> &state, int particleNum);

    static void setPositionOfParticle(vector<Vector3f> &state, int particleNum, Vector3f position);
    static void setVelocityOfParticle(vector<Vector3f> &state, int particleNum, Vector3f velocity);
    static vector<Vector3f> getParticlePositions(vector<Vector3f> &state);

    static float getPressureAtLocation(float desntiyAtLoc, float restDensity, float gasConstant);
    static Vector3f getPressureForce(float mi,
                                     float pi,
                                     float pj,
                                     float densityj,
                                     Vector3f gradKernel);
    static Vector3f getViscosityForce(float mi,
                                      float viscosityConstant,
                                      float densityj,
                                      float viscosityKernelLaplacian,
                                      Vector3f vj,
                                      Vector3f vi);

    // For surface tension:
    static Vector3f calculateSurfaceNormal(float mj, float densityj, Vector3f gradPolySixKernel);
    static float calculateLaplacianColorField(float mj, float densityj, float laplacianPolySixKernel);

  };

#endif /* FORCEUTILITIES_H_ */
