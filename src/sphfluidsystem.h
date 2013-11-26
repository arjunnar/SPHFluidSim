#ifndef SPHFLUIDSYSTEM_H
#define SPHFLUIDSYSTEM_H

#include "particleSystem.h"
#include "ParticleGrid.h"
#include "extra.h"
#include <map>
#include <utility>
#include <vector>
#include "debugutilities.h"

using namespace std;

class SPHFluidSystem : public ParticleSystem
{
   public:
        SPHFluidSystem();
        SPHFluidSystem(int numParticles);
        virtual ~SPHFluidSystem();
        vector<Vector3f> evalF(vector<Vector3f> state);
        void draw();
        void reinitializeSystem();
        void checkCollisions();
   private:
        ParticleGrid particleGrid;
        map<int, float> densityCache;
        float calcDensity(int particleIndex, vector<int> &neighborIndexes, vector<Vector3f> &state);
        void initConstants();
        float PARTICLE_MASS;
        float GRAVITY_CONSTANT;
        float REST_DENSITY;
        float GAS_CONSTANT;
        float GRID_DIMENSION;
        float H_CONSTANT;
        float VISCOSITY_CONSTANT;
};

#endif // SPHFLUIDSYSTEM_H
