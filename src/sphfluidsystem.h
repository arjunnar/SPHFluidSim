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
        void toggleRotation()
        {
            rotEnabled = !rotEnabled;
            rotSpeed = 0.2;
        }

        void moveToFasterRotSpeed()
        {
            if (rotEnabled)
            {
                rotSpeed = 0.35;
            }
        }

        void toggleMarchingCubes()
        {
            renderWithMarchingCubes = !renderWithMarchingCubes;
        }

        void setMovingHorizontally(bool right);

        // For water droplet simulation
        void addSphereOfParticles(Vector3f origin, float radius);

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

        // Other variables
        float MIN_DENSITY;
        float MAX_DENSITY;
        float rotSpeed;
        bool rotEnabled;
        bool renderWithMarchingCubes;
        Vector3f directionOfGravity;

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
        void build3DWaveSystem();

        // For emitting particles
        float angleToEmit;
        float angleToEmit2;
        float emitAngleIncrement;
        float emitAngleIncrement2;
        float emitVelocityConstant;
        Vector3f locOfCannon;
        Vector3f locOfCannon2;
        Vector3f locOfCannon3;

        // For adding a droplet of water
        void addCircleOfParticles(Vector3f origin, float stepSize, float radius);
        void addDiskOfParticles(Vector3f origin, float radius);

        // Not really used anymore
        void buildTwoParticleSystemNeighbors();
        void buildTwoParticleSystemNotNeighbors();
        void buildTestSystem2();
};

#endif // SPHFLUIDSYSTEM_H
