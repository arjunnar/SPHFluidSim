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
#include "MarchingCubes.h"
#include "mpVector.h"

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
        void emitParticle();

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
        float calcDensity(Vector3f pos);
        void initConstants();
        bool isNan(float val);
        bool isNan(Vector3f vec);
        Vector3f getParticleColor(float densityAtParticleLoc);
        void drawFace(TRIANGLE face, bool reverse);

        // 2D system initializations
        void build2DTestSystemSimple();
        void build2DTestSystemLarge();
        void build2DTestSystemTwoDensities();

        // 3D system initializations
        void build3DTestSystemSimple();
        void build3DTestSystemLarge();

        // For emitting particles
        float angleToEmit;
        float angleToEmit2;
        float emitAngleIncrement;
        float emitAngleIncrement2;
        float emitVelocityConstant;
        Vector3f locOfCannon;
        Vector3f locOfCannon2;
        Vector3f locOfCannon3;


        // Not really used anymore
        void buildTwoParticleSystemNeighbors();
        void buildTwoParticleSystemNotNeighbors();
        void buildTestSystem2();
};

#endif // SPHFLUIDSYSTEM_H
