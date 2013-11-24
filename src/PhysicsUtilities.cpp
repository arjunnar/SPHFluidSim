/*
 * ForceUtilities.cpp
 *
 *  Created on: Oct 19, 2013
 *      Author: arjunnar
 */

#include "PhysicsUtilities.h"
#include <vector>

using namespace std;

Vector3f PhysicsUtilities::getSpringForce(Vector3f xi, Vector3f xj, float springConstant, float springRestLength)
{
    Vector3f d = xi - xj;

    float partialCalc = -1.0 * springConstant * (d.abs() - springRestLength);
    Vector3f force = ( partialCalc/d.abs() ) * d;
    return force;
}

Vector3f PhysicsUtilities::getGravityForce(float mass, float gravityConstant)
{
    return Vector3f(0.0, -1.0 * mass * gravityConstant, 0.0);
}

Vector3f PhysicsUtilities::getDragForce(float dragConstant, Vector3f velocityVector)
{
    return -1.0 * dragConstant * velocityVector;
}

Vector3f PhysicsUtilities::getPressureForce(float mi,
                                            float pi,
                                            float pj,
                                            float densityj,
                                            Vector3f gradKernel)
{
    float avgPressure = (pi + pj) / 2.0;
    return (-1.0 * mi * avgPressure / densityj) * gradKernel;
}

Vector3f PhysicsUtilities::getPositionOfParticle(vector<Vector3f> &state, int particleNum)
{
    return state[2 * particleNum];
}

Vector3f PhysicsUtilities::getVelocityOfParticle(vector<Vector3f> &state, int particleNum)
{
    return state[2 * particleNum + 1];
}

void PhysicsUtilities::setPositionOfParticle(vector<Vector3f> &state, int particleNum, Vector3f position)
{
    state[2 * particleNum] = position;
}

void PhysicsUtilities::setVelocityOfParticle(vector<Vector3f> &state, int particleNum, Vector3f velocity)
{
    state[2 * particleNum + 1] = velocity;
}

float PhysicsUtilities::getPressureAtLocation(float densityAtLoc, float restDensity, float gasConstant)
{
    return gasConstant * (densityAtLoc - restDensity);
}




