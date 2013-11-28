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
