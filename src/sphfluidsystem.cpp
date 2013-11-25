#include "sphfluidsystem.h"
#include "PhysicsUtilities.h"
#include "kernelutilities.h"

using namespace std;

SPHFluidSystem::SPHFluidSystem()
{
    SPHFluidSystem(10000); // Use 10000 particles by default
}

SPHFluidSystem::SPHFluidSystem(int numParticles) : ParticleSystem(numParticles)
{
    PARTICLE_MASS = 1.0;
    GRAVITY_CONSTANT = 1.0;
    REST_DENSITY = 1.0;
    GAS_CONSTANT = 2.0;
    GRID_DIMENSION = 0.2;
    H_CONSTANT = GRID_DIMENSION * 2 * sqrt(3);




    Vector3f pointOne(1.5, 1.5, 1.5);
    Vector3f pointTwo(1.51, 1.5, 1.5);
    Vector3f pointThree(1.5, 1.5, 1.51);


    m_vVecState.push_back(pointOne);
    m_vVecState.push_back(Vector3f::ZERO);
    m_vVecState.push_back(pointTwo);
    m_vVecState.push_back(Vector3f::ZERO);
    m_vVecState.push_back(pointThree);
    m_vVecState.push_back(Vector3f::ZERO);


    for (int i = 0; i < 10; ++i)
    {
    	for (int k = 0; k < 10; ++k) {
            Vector3f point(1.0 + i * .01, 1.0, 1.0 +  k * .01);
            m_vVecState.push_back(point);
            m_vVecState.push_back(Vector3f::ZERO);
    	}

    }



    /*
    for (int i = 0; i < 5; ++i)
    {
        Vector3f point(1.0 + i * .01, 1.005, 1.0);
        m_vVecState.push_back(point);
        m_vVecState.push_back(Vector3f::ZERO);
    }
    */

    Vector3f origin = Vector3f::ZERO;
    particleGrid = ParticleGrid(origin, 2.0);
    vector<Vector3f> particlePositions = PhysicsUtilities::getParticlePositions(m_vVecState);
    particleGrid.initializeGrid(particlePositions);
}

SPHFluidSystem::~SPHFluidSystem()
{

}

vector<Vector3f> SPHFluidSystem::evalF(vector<Vector3f> state)
{
    //vector<Vector3f> particlePositions = PhysicsUtilities::getParticlePositions(state);
    //particleGrid.initializeGrid(particlePositions);

    vector<Vector3f> derivative;
    for (int particleIndex = 0; particleIndex < m_numParticles; ++particleIndex)
    {
        Vector3f position = PhysicsUtilities::getPositionOfParticle(state, particleIndex);

        vector<int> neighborIndexes = particleGrid.getNeighborParticleIndexes(particleIndex, position);

        float densityAtParticleLoc = calcDensity(particleIndex, neighborIndexes, state); // densityi

        //cout << "Density at particle loc: " << densityAtParticleLoc << "\n";

        // CALCULATE PRESSURE FORCE
        Vector3f totalPressureForce = Vector3f::ZERO;
        float pressureAtParticleLoc = PhysicsUtilities::getPressureAtLocation(densityAtParticleLoc,
                                                                              REST_DENSITY,
                                                                              GAS_CONSTANT); // pi
        //cout << "Neighbor indexes: ";
        //DebugUtilities::printVectorInt(neighborIndexes);
        //cout << endl;

        for (int neighborI : neighborIndexes)
        {
            Vector3f neighborPos = PhysicsUtilities::getPositionOfParticle(state, neighborI);

            // Calculate density at neighbor location
            vector<int> neighborsOfNeighbor = particleGrid.getNeighborParticleIndexes(neighborI, neighborPos);

            float densityAtNeighborLoc = calcDensity(neighborI,
                                                     neighborsOfNeighbor,
                                                     state); // densityj

            // Calculate pressure at neighbor location
            float pressureAtNeighborLoc = PhysicsUtilities::getPressureAtLocation(densityAtNeighborLoc,
                                                                                  REST_DENSITY,
                                                                                  GAS_CONSTANT); // pj

            Vector3f spikyKernelGrad = KernelUtilities::gradSpikyKernel(position - neighborPos, H_CONSTANT);
            Vector3f pressureContribution = PhysicsUtilities::getPressureForce( PARTICLE_MASS,
                                                                                pressureAtParticleLoc,
                                                                                pressureAtNeighborLoc,
                                                                                densityAtNeighborLoc,
                                                                                spikyKernelGrad );

            totalPressureForce += pressureContribution;
        }

        // CALCULATE VISCOSITY FORCE
        Vector3f totalViscosityForce = Vector3f::ZERO;

        // CALCULATE SURFACE TENSION FORCE
        Vector3f surfaceTensionForce = Vector3f::ZERO;

        // Add the results to the derivative
        Vector3f totalForce = totalPressureForce + totalViscosityForce + surfaceTensionForce;
        //cout << "Total Force on particle " << particleIndex << " : " << " < " << totalForce.x() << " , " << totalForce.y()
          //   << " , " << totalForce.z() << " > " << endl;

        Vector3f velocity = PhysicsUtilities::getVelocityOfParticle(state, particleIndex);
        Vector3f acceleration = totalForce / densityAtParticleLoc;

        // HACKY
        float forceEpsilon = 0.005;
        if (totalForce.abs() > forceEpsilon)
        {
            acceleration = totalForce / densityAtParticleLoc;

        }

        else
        {
            acceleration = Vector3f::ZERO;
        }



        derivative.push_back(velocity);
        derivative.push_back(acceleration);
    }

    densityCache.clear();

    return derivative;
}

void SPHFluidSystem::checkCollisions()
{
    Vector3f leftFaceNormal(-1.0, 0.0, 0.0);
    Vector3f rightFaceNormal = -1.0 * leftFaceNormal;

    Vector3f bottomFaceNormal(0.0, -1.0, 0.0);
    Vector3f topFaceNormal = -1.0 * bottomFaceNormal;

    Vector3f backFaceNormal(0.0, 0.0, -1.0);
    Vector3f frontFaceNormal = -1.0 * frontFaceNormal;

    for (int i = 0; i < m_numParticles; ++i)
    {
        Vector3f pos = PhysicsUtilities::getPositionOfParticle(m_vVecState, i);
        Vector3f vel = PhysicsUtilities::getVelocityOfParticle(m_vVecState, i);

        //if ( !(0 <= pos.x() && pos.x <= particleGrid.getGridSize()) )
        //{

        //}
    }
}

void SPHFluidSystem::draw()
{
    for (int i = 0; i < m_numParticles; i++)
    {
        // Draw the particles
        Vector3f posParticle = PhysicsUtilities::getPositionOfParticle(m_vVecState, i);
        glPushMatrix();
        glTranslatef(posParticle[0], posParticle[1], posParticle[2] );
        glutSolidSphere(0.0250f,10.0f,10.0f);
        glPopMatrix();
    }
}

void SPHFluidSystem::reinitializeSystem()
{

}

// Helper functions
float SPHFluidSystem::calcDensity(int particleIndex, vector<int> &neighborIndexes, vector<Vector3f> &state)
{
    float density;
    //cout << "In method calcDensity" << endl;
    if (densityCache.find(particleIndex) != densityCache.end())
    {
        density = densityCache[particleIndex];
    }

    else
    {
        float densitySum = 0.0;
        Vector3f particleLoc = PhysicsUtilities::getPositionOfParticle(state, particleIndex);
       // cout << "Not already in cache" << endl;
        for (int neighborI : neighborIndexes)
        {
           Vector3f neighborLoc = PhysicsUtilities::getPositionOfParticle(state, neighborI);
           densitySum += PARTICLE_MASS * KernelUtilities::spikyKernel(particleLoc - neighborLoc, H_CONSTANT);
         //  cout << "Spiky kernel val: "  << KernelUtilities::spikyKernel(particleLoc - neighborLoc, 100) << endl;
        }

        density = densitySum;
        densityCache.insert(std::pair<int, float>(particleIndex, density));
    }

    return density;
}

