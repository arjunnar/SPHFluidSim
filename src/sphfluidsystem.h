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
#include "Utilities.cpp"

using namespace std;

class SPHFluidSystem : public ParticleSystem
{
   public:
        SPHFluidSystem();
        SPHFluidSystem(float boxSizeX, float boxSizeY, float boxSizeZ, FluidSystemType fluidSystemType);
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

        FluidSystemType typeOfSystem;

        // Physics constants
        float PARTICLE_MASS;
        float GRAVITY_CONSTANT;
        float REST_DENSITY;
        float REST_DENSITY_FIRST;
        float REST_DENSITY_SECOND;
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

        // 2D system initializations
        void build2DTestSystemSimple();
        void build2DTestSystemLarge();
        void build2DTestSystemTwoDensities();

        // 3D system initializations
        void build3DTestSystemSimple();
        void build3DTestSystemLarge();

        // Not really used anymore
        void buildTwoParticleSystemNeighbors();
        void buildTwoParticleSystemNotNeighbors();
        void buildTestSystem2();
};

#endif // SPHFLUIDSYSTEM_H
