#ifndef PENDULUMSYSTEM_H
#define PENDULUMSYSTEM_H

#include "extra.h"
#include <vector>

#include "particleSystem.h"

class PendulumSystem: public ParticleSystem
{

private:
    static constexpr float PARTICLE_MASS = 0.10;
    static constexpr float GRAVITY_CONSTANT = 1.0;

    static constexpr float SPRING_CONSTANT = 8.0;
    static constexpr float SPRING_REST_LENGTH = 0.30;
    static constexpr float DRAG_CONSTANT = 0.005;

    vector<int> getParticlesConnectedTo(int i);
    void setInitState();

public:
	PendulumSystem(int numParticles);
	vector<Vector3f> evalF(vector<Vector3f> state);
	void draw();
	void reinitializeSystem();
};

#endif
