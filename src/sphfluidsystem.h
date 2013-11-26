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

        // Functions inherited from ParticleSystem
        vector<Vector3f> evalF(vector<Vector3f> state);
        void draw();
        void reinitializeSystem();

   private:
        // Instance variables
        ParticleGrid particleGrid;
        map<int, float> densityCache;

        // Physics constants
        float PARTICLE_MASS;
        float GRAVITY_CONSTANT;
        float REST_DENSITY;
        float GAS_CONSTANT;
        float GRID_DIMENSION;
        float H_CONSTANT;
        float VISCOSITY_CONSTANT;

        // Helper functions
        float calcDensity(int particleIndex, vector<int> &neighborIndexes, vector<Vector3f> &state);
        void initConstants();

        // Different system initializations
        void testOneInitializeSystem();
        void buildTwoParticleSystemNeighbors();
        void buildTwoParticleSystemNotNeighbors();

};

#endif // SPHFLUIDSYSTEM_H
