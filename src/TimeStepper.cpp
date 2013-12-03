#include "TimeStepper.h"

void ForwardEuler::takeStep(ParticleSystem* particleSystem, float stepSize)
{
    vector<Vector3f> origState = particleSystem->getState();
    vector<Vector3f> derivativeAtState = particleSystem->evalF(origState);

    vector<Vector3f> newStateVectorArray;

    for (int i = 0; i < origState.size(); ++i)
    {
        Vector3f newStateVector3f = origState.at(i) + stepSize * derivativeAtState.at(i);
        newStateVectorArray.push_back(newStateVector3f);
    }

    particleSystem->setState(newStateVectorArray);
}

void Trapezoidal::takeStep(ParticleSystem* particleSystem, float stepSize)
{
    vector<Vector3f> origState = particleSystem->getState();
    vector<Vector3f> f0 = particleSystem->evalF(origState);

    vector<Vector3f> afterEulerStepState;

    for (int i = 0; i < origState.size(); ++i)
    {
        Vector3f newStateVector3f = origState.at(i) + stepSize * f0.at(i);
        afterEulerStepState.push_back(newStateVector3f);
    }

    vector<Vector3f> f1 = particleSystem->evalF(afterEulerStepState);

    vector<Vector3f> newStateVectorArray;
    for (int i = 0; i < origState.size(); ++i)
    {
        Vector3f newStateVector3f = origState.at(i) + stepSize/2.0 * (f0.at(i) + f1.at(i));
        newStateVectorArray.push_back(newStateVector3f);
    }

    particleSystem->setState(newStateVectorArray);
}

void LeapFrog::takeStep(ParticleSystem *particleSystem, float stepSize)
{
    vector<Vector3f> currentState = particleSystem->getState();

    // Compute new positions using current velocity
    vector<Vector3f> newPositionsCurrentVelocities = vector<Vector3f>();
    for (int i = 0; i < currentState.size(); i += 2)
    {
        Vector3f currentPos = currentState[i];
        Vector3f currentVel = currentState[i + 1];
        Vector3f newPos = currentPos + stepSize * currentVel;
        newPositionsCurrentVelocities.push_back(newPos);
        newPositionsCurrentVelocities.push_back(currentVel);
    }

    // Compute forces are updated positions
    vector<Vector3f> deriv = particleSystem->evalF(newPositionsCurrentVelocities);

    // Update velocities using computed forces
    vector<Vector3f> newPositionsNewVelocities = vector<Vector3f>();
    for (int i = 0; i < newPositionsCurrentVelocities.size(); i += 2)
    {
        Vector3f newPos = newPositionsCurrentVelocities[i];
        Vector3f currentVel = newPositionsCurrentVelocities[i + 1];
        Vector3f accelFromUpdatedPositions = deriv[i + 1];
        Vector3f newVel = currentVel + stepSize * accelFromUpdatedPositions;
        newPositionsNewVelocities.push_back(newPos);
        newPositionsNewVelocities.push_back(newVel);
    }

    particleSystem->setState(newPositionsNewVelocities);
}

