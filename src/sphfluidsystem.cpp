#include "sphfluidsystem.h"
#include "PhysicsUtilities.h"
#include "kernelutilities.h"

SPHFluidSystem::SPHFluidSystem()
{
    SPHFluidSystem(10000); // Use 10000 particles by default
}

SPHFluidSystem::SPHFluidSystem(int numParticles) : ParticleSystem(numParticles)
{

}

SPHFluidSystem::~SPHFluidSystem()
{

}

vector<Vector3f> SPHFluidSystem::evalF(vector<Vector3f> state)
{
    vector<Vector3f> derivative;
    for (int particleIndex = 0; particleIndex < m_numParticles; ++particleIndex)
    {
        Vector3f position = PhysicsUtilities::getPositionOfParticle(state, particleIndex);

        vector<int> neighborIndexes = particleGrid.getNeighborParticleIndexes(particleIndex, position);
        float densityAtParticleLoc = calcDensity(particleIndex, neighborIndexes, state); // densityi

        // CALCULATE PRESSURE FORCE
        Vector3f totalPressureForce = Vector3f::ZERO;
        float pressureAtParticleLoc = PhysicsUtilities::getPressureAtLocation(densityAtParticleLoc,
                                                                              REST_DENSITY,
                                                                              GAS_CONSTANT); // pi
        for (int neighborI : neighborIndexes)
        {
            Vector3f neighborPos = PhysicsUtilities::getPositionOfParticle(state, neighborI);

            // Calculate density at neighbor location
            vector<int> neighborsOfNeighbor = particleGrid.getNeighborParticleIndexes(neighborI, neighborPos);

            float densityAtNeighborLoc = calcDensity(neighborI,
                                                     neighborsOfNeighbor,
                                                     state); // densityj

            // Calculate pressure at neighbor location
            float pressureAtNeighborLoc = PhysicsUtilities::getPressureAtLocation(densityAtNeighborLoc, REST_DENSITY, GAS_CONSTANT); // pj

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
        Vector3f velocity = PhysicsUtilities::getVelocityOfParticle(state, particleIndex);
        Vector3f acceleration = (1.0 / densityAtParticleLoc) * totalForce;

        derivative.push_back(velocity);
        derivative.push_back(acceleration);
    }

    return derivative;
}


void SPHFluidSystem::draw()
{
    glDisable(GL_LIGHTING);
    glColor3f(245/255.0, 15.0/255.0, 222.0/255.0);

    for (int i = 0; i < m_numParticles; i++)
    {
        // Draw the particles
        Vector3f posParticle = PhysicsUtilities::getPositionOfParticle(this->m_vVecState, i);
        glPushMatrix();
        glTranslatef(posParticle[0], posParticle[1], posParticle[2] );
        glutSolidSphere(0.075f,10.0f,10.0f);
        glPopMatrix();

    }

    glEnable(GL_LIGHTING);
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
        density = densityCache[particleIndex];
    }

    else
    {
        float densitySum = 0.0;
        Vector3f particleLoc = PhysicsUtilities::getPositionOfParticle(state, particleIndex);

        for (int neighborI : neighborIndexes)
        {
           Vector3f neighborLoc = PhysicsUtilities::getPositionOfParticle(state, neighborI);
           densitySum += PARTICLE_MASS * KernelUtilities::spikyKernel(particleLoc - neighborLoc, H_CONSTANT);
        }

        density = densitySum;
        densityCache.insert(std::pair<int, float>(particleIndex, density));
    }

    return density;
}

