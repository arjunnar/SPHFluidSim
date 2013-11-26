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
    initConstants();
    buildTwoParticleSystemNeighbors();

	Vector3f origin = Vector3f::ZERO;
    particleGrid = ParticleGrid(origin, 0.5, 0.9, 0.5);
}

void SPHFluidSystem::initConstants()
{
    PARTICLE_MASS = 0.00065;
    GRAVITY_CONSTANT = 0.5;
	REST_DENSITY = 0.0;
    GAS_CONSTANT = 0.0060;
	GRID_DIMENSION = 0.015;
	H_CONSTANT = 0.02289;
	VISCOSITY_CONSTANT = 0.0072;
}

SPHFluidSystem::~SPHFluidSystem()
{

}

vector<Vector3f> SPHFluidSystem::evalF(vector<Vector3f> state)
{
    vector<Vector3f> particlePositionsInState = PhysicsUtilities::getParticlePositions(state);
    particleGrid.initializeGrid(particlePositionsInState);

    vector<Vector3f> derivative;
    for (int particleIndex = 0; particleIndex < m_numParticles; ++particleIndex)
    {
        Vector3f positionOfParticle = PhysicsUtilities::getPositionOfParticle(state, particleIndex);
        Vector3f velocityOfParticle = PhysicsUtilities::getVelocityOfParticle(state, particleIndex);

        // Get the neighbors of the particle in the passed in state
        vector<int> neighborIndexes = particleGrid.getNeighborParticleIndexes(particleIndex, positionOfParticle);

        float densityAtParticleLoc = calcDensity(particleIndex, neighborIndexes, state); // densityi

        // CALCULATE ACCELERATION FROM PRESSURE AND VISCOSITY
        Vector3f totalPressureForce = Vector3f::ZERO;
        Vector3f totalViscosityForce = Vector3f::ZERO;

        float pressureAtParticleLoc = PhysicsUtilities::getPressureAtLocation(densityAtParticleLoc,
                                                                              REST_DENSITY,
                                                                              GAS_CONSTANT); // pi
        for (int neighborI : neighborIndexes)
        {
            Vector3f neighborPos = PhysicsUtilities::getPositionOfParticle(state, neighborI);

            vector<int> neighborsOfNeighbor = particleGrid.getNeighborParticleIndexes(neighborI, neighborPos);
            float densityAtNeighborLoc = calcDensity(neighborI,neighborsOfNeighbor, state); // densityj

            // CALCULATE PRESSURE CONTRIBUTION
            Vector3f pressureContribution;
            if (neighborsOfNeighbor.size() == 0)
            {
            	pressureContribution = Vector3f::ZERO;
            }

            else
            {
            	float pressureAtNeighborLoc = PhysicsUtilities::getPressureAtLocation(densityAtNeighborLoc,
                                                                                      REST_DENSITY,
                                                                                      GAS_CONSTANT); // pj

                Vector3f spikyKernelGrad = KernelUtilities::gradSpikyKernel(positionOfParticle - neighborPos, H_CONSTANT);

                pressureContribution = PhysicsUtilities::getPressureForce( PARTICLE_MASS,
                                                                           pressureAtParticleLoc,
                                                                           pressureAtNeighborLoc,
                                                                           densityAtNeighborLoc,
                                                                           spikyKernelGrad );
            }

            totalPressureForce += pressureContribution;

            // CALCULATE VISCOSITY CONTRIBUTION
            Vector3f viscosityContribution;
            if (neighborsOfNeighbor.size() == 0)
            {
            	viscosityContribution = Vector3f::ZERO;
            }

            else
            {
                float laplacianKernel = KernelUtilities::laplacianViscosityKernel(positionOfParticle - neighborPos, H_CONSTANT);
            	Vector3f velocityNeighbor = PhysicsUtilities::getVelocityOfParticle(state, neighborI);
            	viscosityContribution = PhysicsUtilities::getViscosityForce( PARTICLE_MASS,
            																 VISCOSITY_CONSTANT,
            																 densityAtNeighborLoc,
            																 laplacianKernel,
            																 velocityNeighbor,
                                                                             velocityOfParticle);
            }

            totalViscosityForce += viscosityContribution;
        }

        Vector3f accelPressure;
        Vector3f accelViscosity;

        float densityEpsilon = 0.00005;

        // NOT SURE IF THIS IS NEEDED, SEEMS HACKY
        if (densityAtParticleLoc < densityEpsilon)
        {
            accelPressure = Vector3f::ZERO;
            accelViscosity = Vector3f::ZERO;
        }

        else
        {
			accelPressure = totalPressureForce / densityAtParticleLoc;
			accelViscosity = totalViscosityForce / densityAtParticleLoc;
        }

        // CALCULATE ACCELERATION FROM SURFACE TENSION
        Vector3f accelSurfaceTension = Vector3f::ZERO;

        // CALCULATE ACCELERATION FROM GRAVITY
        Vector3f accelGravity = PhysicsUtilities::getGravityForce(PARTICLE_MASS, GRAVITY_CONSTANT) / PARTICLE_MASS;

        // Add the results to the derivative
        Vector3f totalAccel = accelPressure + accelViscosity + accelSurfaceTension + accelGravity;

        derivative.push_back(velocityOfParticle);
        derivative.push_back(totalAccel);
    }

    densityCache.clear();

    return derivative;
}

void SPHFluidSystem::draw()
{
    for (int i = 0; i < m_numParticles; i++)
    {
        // Draw the particles
        Vector3f posParticle = PhysicsUtilities::getPositionOfParticle(m_vVecState, i);
        glPushMatrix();
        glTranslatef(posParticle[0], posParticle[1], posParticle[2] );
        glutSolidSphere(0.0150f,10.0f,10.0f);
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

    if (densityCache.find(particleIndex) != densityCache.end())
    {
        // Already computed the density for the particle in the given state. Just get it from the cache.
        density = densityCache[particleIndex];
    }

    else
    {
        // Compute the density using Eq(3)
        density = 0.0;
        Vector3f particleLoc = PhysicsUtilities::getPositionOfParticle(state, particleIndex);

        for (int neighborI : neighborIndexes)
        {
           Vector3f neighborLoc = PhysicsUtilities::getPositionOfParticle(state, neighborI);
           density += PARTICLE_MASS * KernelUtilities::polySixKernel(particleLoc - neighborLoc, H_CONSTANT);
        }

        densityCache.insert(std::pair<int, float>(particleIndex, density));
    }

    return density;
}

// Different system initializtions
void SPHFluidSystem::buildTwoParticleSystemNotNeighbors()
{
    Vector3f pos1(0.3, 0.5, 0.1);
    Vector3f pos2(0.1, 0.5, 0.1);

    m_vVecState.push_back(pos1);
    m_vVecState.push_back(Vector3f::ZERO);
    m_vVecState.push_back(pos2);
    m_vVecState.push_back(Vector3f::ZERO);
}

void SPHFluidSystem::buildTwoParticleSystemNeighbors()
{
    Vector3f pos1(0.3, 0.5, 0.1);
    Vector3f pos2(0.305, 0.5, 0.1);

    m_vVecState.push_back(pos1);
    m_vVecState.push_back(Vector3f::ZERO);
    m_vVecState.push_back(pos2);
    m_vVecState.push_back(Vector3f::ZERO);
}

void  SPHFluidSystem::testOneInitializeSystem()
{
    for (int k = 0; k < 1; k++)
	{
		for (int i = 0; i < 20; i++ )
		{
			for (int j = 0; j < 40; j++) {
	            Vector3f point(0.14 + .015 * i, 0.24 + j * .015, 0.2 +  k * .015);
	            m_vVecState.push_back(point);
	            m_vVecState.push_back(Vector3f::ZERO);
			}
		}
	}
}
