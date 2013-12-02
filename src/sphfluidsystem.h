#ifndef SPHFLUIDSYSTEM_H
#define SPHFLUIDSYSTEM_H

#include "particleSystem.h"
#include "ParticleGrid.h"
#include "extra.h"
#include <map>
#include <utility>
#include <vector>
#include "debugutilities.h"
#include "assert.h"

using namespace std;

class SPHFluidSystem : public ParticleSystem
{
   public:
        SPHFluidSystem();
        SPHFluidSystem(float boxSizeX, float boxSizeY, float boxSizeZ);
        virtual ~SPHFluidSystem();

        // Functions inherited from ParticleSystem
        vector<Vector3f> evalF(vector<Vector3f> state);
        void draw();
        void reinitializeSystem();
        void advanceState();

   private:
        // Instance variables
        ParticleGrid particleGrid;
        //map<int, float> densityCache;
        vector<float> vecParticleDensities;
        vector<float> vecParticlePressures;

        // For leapfrog integration
        float timeStep;
        vector<Vector3f> vHalfStepBefore;

        // Physics constants
        float PARTICLE_MASS;
        float GRAVITY_CONSTANT;
        float REST_DENSITY;
        float GAS_CONSTANT;
        float VISCOSITY_CONSTANT;
        float TENSION_CONSTANT;
        float TENSION_THRESHOLD;
        float SELF_DENSITY_CONSTANT;
        float SELF_LAPLACIAN_COLOR_FIELD;

        // Other constants
        float MIN_DENSITY;
        float MAX_DENSITY;

        // Helper functions
        void calculateDensitiesAndPressures(vector<Vector3f> &state);
        void initConstants();
        bool isNan(float val);
        bool isNan(Vector3f vec);
        Vector3f getParticleColor(float densityAtParticleLoc);

        // Different system initializations
        void build2DTestSystem();
        void build2DTestSystem2();
        void build3DTestSystem();

        void buildTwoParticleSystemNeighbors();
        void buildTwoParticleSystemNotNeighbors();
        void buildTestSystem2();
};

#endif // SPHFLUIDSYSTEM_H
