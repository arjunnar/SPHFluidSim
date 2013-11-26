#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#include <vector>
#include <vecmath.h>
#include <iostream>

using namespace std;

class ParticleSystem
{
public:

    int m_numParticles;

	ParticleSystem(int numParticles = 0);
	
	// For a given state, evaluate derivative f(X,t)
	virtual vector<Vector3f> evalF(vector<Vector3f> state) = 0;
	
	// Getter method for the system's state
    vector<Vector3f> getState(){ return m_vVecState; }
	
	// Setter method for the system's state
    virtual void setState(const vector<Vector3f>  & newState) { m_vVecState = newState; }
	
	virtual void reinitializeSystem() = 0;

	virtual void draw() = 0;
	
protected:

	vector<Vector3f> m_vVecState;
	
};

#endif
