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
        static constexpr float PARTICLE_MASS = 1.0;
        static constexpr float GRAVITY_CONSTANT = 1.0;
        static constexpr float REST_DENSITY = 1.0;
        static constexpr float GAS_CONSTANT = 2.0;
        static constexpr float GRID_DIMENSION = 0.2;
        static constexpr float H_CONSTANT = GRID_DIMENSION * 2 * sqrt(3);
};

#endif // SPHFLUIDSYSTEM_H
