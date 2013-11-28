#include "sphfluidsystem.h"
#include "PhysicsUtilities.h"
#include "kernelutilities.h"

using namespace std;

SPHFluidSystem::SPHFluidSystem()
{
    SPHFluidSystem(5000);
}

SPHFluidSystem::SPHFluidSystem(int numParticles) : ParticleSystem(numParticles)
{
    initConstants();
    testOneInitializeSystem();

	Vector3f origin = Vector3f::ZERO;
    particleGrid = ParticleGrid(origin, 0.5 , 0.9, 0.5);
}

void SPHFluidSystem::initConstants()
{
    PARTICLE_MASS = 0.02;
    GRAVITY_CONSTANT = 6.5;
    REST_DENSITY = 1000.0;
    GAS_CONSTANT = 1.0;
    VISCOSITY_CONSTANT = 6.0; // Was previously 6.0;
    TENSION_CONSTANT = 0.10;
}

SPHFluidSystem::~SPHFluidSystem()
{

}

vector<Vector3f> SPHFluidSystem::evalF(vector<Vector3f> state)
{
    vector<Vector3f> particlePositionsInState = PhysicsUtilities::getParticlePositions(state);
    particleGrid.initializeGrid(particlePositionsInState);
    calculateDensitiesAndPressures();

    vector<Vector3f> derivative;
    for (int particleIndex = 0; particleIndex < m_numParticles; ++particleIndex)
    {
        Vector3f positionOfParticle = PhysicsUtilities::getPositionOfParticle(state, particleIndex);
        Vector3f velocityOfParticle = PhysicsUtilities::getVelocityOfParticle(state, particleIndex);

        // Get the neighbors of the particle
        vector<int> neighborIndexes = particleGrid.getNeighborParticleIndexes(particleIndex, positionOfParticle);

        float densityAtParticleLoc = vecParticleDensities[particleIndex];
        float pressureAtParticleLoc = vecParticlePressures[particleIndex];

        // CALCULATE ACCELERATION FROM PRESSURE AND VISCOSITY
        // ALSO CALCULATE PIECES OF THE SURFACE TENSION FORCE
        Vector3f totalPressureForce = Vector3f::ZERO;
        Vector3f totalViscosityForce = Vector3f::ZERO;

        for (int neighborI : neighborIndexes)
        {
            Vector3f neighborPos = PhysicsUtilities::getPositionOfParticle(state, neighborI);

            float densityAtNeighborLoc = vecParticleDensities[neighborI];

            Vector3f spikyKernelGradForDebugging;

            // CALCULATE PRESSURE CONTRIBUTION FROM THE NEIGHBOR
            Vector3f pressureContribution;
            if (false)
            {
            	pressureContribution = Vector3f::ZERO;
            }

            else
            {
                float pressureAtNeighborLoc = vecParticlePressures[neighborI];
                Vector3f rForPressureKernel = positionOfParticle - neighborPos;

                // We don't want particles to collapse to the same location. This would give us infinite pressure.
                // So, if the particles are too close, just separate them a little bit.
                float rPressureKernelEpsilon = 0.00005;
                if (rForPressureKernel.abs() < rPressureKernelEpsilon)
                {
                     // If 2D simulation, leave the z component 0.0. If 3D simulation, make the z component 0.003 as well.
                    rForPressureKernel = Vector3f(0.002, 0.002, 0.0);
                }

                Vector3f spikyKernelGrad = KernelUtilities::gradSpikyKernel(rForPressureKernel, H_CONSTANT);
                spikyKernelGradForDebugging = spikyKernelGrad;

                if (isNan(spikyKernelGrad))
                {
                    cout << "NAN spiky kernel gradient" << endl;
                    cout << "Spiky kernel gradient: "; DebugUtilities::printVector3f(spikyKernelGrad); cout << endl;
                    cout << "ri - rj magnitude: " << (positionOfParticle - neighborPos).abs() << endl;
                }

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

            // CALCULATE SURFACE TENSION CONTRIBUTION
            Vector3f surfaceTensContribution;
            if (neighborsOfNeighbor.size() == 0)
            {

                // If the neighbor has no neighbors then the density at the neighbor location will be zero.
                // We don't want to divide by zero, so we just set the surface tension contribution to be zero.
                surfaceTensContribution = Vector3f::ZERO;
            }

            else
            {
                Vector3f gradPolySixKernel = KernelUtilities::gradPolySixKernel(positionOfParticle - neighborPos, H_CONSTANT);

                Vector3f surfaceNormal = PhysicsUtilities::calculateSurfaceNormal(PARTICLE_MASS, densityAtNeighborLoc, gradPolySixKernel);

                float surfaceNormalEpsilon = 0.05;

                // Only evaluate the force if the magnitude of n exceeds a certain threshold
                if (surfaceNormal.abs() < surfaceNormalEpsilon)
                {
                    surfaceTensContribution = Vector3f::ZERO;
                }

                else
                {
                    float laplacianPolySixKernel = KernelUtilities::laplacianPolySixKernel(positionOfParticle - neighborPos, H_CONSTANT);
                    float laplacianColorField = PhysicsUtilities::calculateLaplacianColorField(PARTICLE_MASS,
                                                                                               densityAtNeighborLoc,
                                                                                               laplacianPolySixKernel);

                    surfaceTensContribution = -1.0 * TENSION_CONSTANT * laplacianColorField * surfaceNormal.normalized();
                }
            }

            totalSurfaceTensionForce += surfaceTensContribution;

            if(isNan(viscosityContribution))
            {
                cout << "Encountered NAN viscosity contribution" << endl;
                cout << "Viscosity contribution: "; DebugUtilities::printVector3f(viscosityContribution);

                cout << "Particle loc: "; DebugUtilities::printVector3f(positionOfParticle);
                cout << "Particle velocity: "; DebugUtilities::printVector3f(velocityOfParticle);

                cout << "NeighborI: " << neighborI << endl;
                cout << "Neighbor loc: "; DebugUtilities::printVector3f(neighborPos);
                Vector3f velocityNeighbor = PhysicsUtilities::getVelocityOfParticle(state, neighborI);
                cout << "Neighbor velocity: "; DebugUtilities::printVector3f(velocityNeighbor);
                cout << "Density at neighbor loc: " << densityAtNeighborLoc << endl;

                float laplacianKernel = KernelUtilities::laplacianViscosityKernel(positionOfParticle - neighborPos, H_CONSTANT);
                cout << "Laplacian Kernel: " << laplacianKernel << endl;

                assert(false); // Kill execution of the program
            }

            if(isNan(pressureContribution))
            {
                cout << "Encountered NAN pressure contribution" << endl;
                cout << "Pressure contribution: "; DebugUtilities::printVector3f(pressureContribution);

                cout << "Particle loc: "; DebugUtilities::printVector3f(positionOfParticle);
                cout << "Pressure at particle loc: " << pressureAtParticleLoc << endl;

                cout << "NeighborI: " << neighborI << endl;
                cout << "Neighbor loc: "; DebugUtilities::printVector3f(neighborPos);
                cout << "Density at neighbor loc: " << densityAtNeighborLoc << endl;
                float pressureAtNeighborLoc = PhysicsUtilities::getPressureAtLocation(densityAtNeighborLoc, REST_DENSITY, GAS_CONSTANT);
                cout << "Pressure at neighbor loc: " << pressureAtNeighborLoc << endl;

                cout << "Spiky kernel grad: "; DebugUtilities::printVector3f(spikyKernelGradForDebugging);

                assert(false); // Kill execution of the program
            }
        }

        Vector3f accelPressure;
        Vector3f accelViscosity;
        Vector3f accelSurfaceTension;

        float densityEpsilon = 0.0005;

        // NOT SURE IF THIS IS NEEDED, SEEMS HACKY
        if (densityAtParticleLoc < densityEpsilon)
        {
            accelPressure = Vector3f::ZERO;
            accelViscosity = Vector3f::ZERO;
            accelSurfaceTension = Vector3f::ZERO;
        }

        else
        {
			accelPressure = totalPressureForce / densityAtParticleLoc;
			accelViscosity = totalViscosityForce / densityAtParticleLoc;
            accelSurfaceTension = totalSurfaceTensionForce / densityAtParticleLoc;
        }

        //cout << "Accel surface tension: "; DebugUtilities::printVector3f(accelSurfaceTension);

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
    int numNanPositions = 0;
    for (int i = 0; i < m_numParticles; i++)
    {
        // Draw the particles
        Vector3f posParticle = PhysicsUtilities::getPositionOfParticle(m_vVecState, i);
        if (isNan(posParticle))
        {
            cout << "Encountered NAN position: " << "( " << posParticle.x() << " , " << posParticle.y() << " ," << posParticle.z() << " )" << endl;
            //++numNanPositions;
        }
        glPushMatrix();
        glTranslatef(posParticle[0], posParticle[1], posParticle[2] );
        glutSolidSphere(0.007550f,10.0f,10.0f);
        glPopMatrix();
    }

    //cout << "Num nans: " << numNanPositions << endl;
}

void SPHFluidSystem::reinitializeSystem()
{

}

// Helper functions
float SPHFluidSystem::calculateDensitiesAndPressures(vector<Vector3f> &state)
{
    vecParticleDensities = vector<float>();
    vecParticlePressures = vector<float>();

    for (int i = 0; i < m_numParticles; ++i)
    {
        float density = 0;
        Vector3f pos = PhysicsUtilities::getPositionOfParticle(state, i);
        vector<int> neighborIndexes = particleGrid.getNeighborParticleIndexes(i, pos);

        for (int neigborI : neighborIndexes)
        {
            Vector3f neighborPos = PhysicsUtilities::getPositionOfParticle(state, neighborI);
            density += mj * KernelUtilities::polySixKernel(pos - neighborPos);
        }

        density += SELF_DENSITY_CONSTANT;
        vecParticleDensities.push_back(density);

        float pressure = PhysicsUtilities::getPressureAtLocation(density, REST_DENSITY, GAS_CONSTANT);
        vecParticleDensities.push_back(pressure);
    }
}
/*
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

        density += PARTICLE_MASS * KernelUtilities::polySixKernel(particleLoc - particleLoc, H_CONSTANT);

        densityCache.insert(std::pair<int, float>(particleIndex, density));
    }

    return density;
}
*/

bool SPHFluidSystem::isNan(float val)
{
    return val != val;
}

bool SPHFluidSystem::isNan(Vector3f vec)
{
    return isNan(vec.x()) || isNan(vec.y()) || isNan(vec.z());
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

    m_numParticles = 2;
}

void SPHFluidSystem::buildTwoParticleSystemNeighbors()
{
    Vector3f pos1(0.3, 0.5, 0.1);
    Vector3f pos2(0.305, 0.5, 0.1);

    m_vVecState.push_back(pos1);
    m_vVecState.push_back(Vector3f::ZERO);
    m_vVecState.push_back(pos2);
    m_vVecState.push_back(Vector3f::ZERO);

    m_numParticles = 2;
}

void SPHFluidSystem::testOneInitializeSystem()
{
    for (int k = 0; k < 15; k++)
    {
        for (int i = 0; i < 15; i++ )
        {
            for (int j = 0; j < 20; j++) {
                Vector3f point(0.14 + .01 * i + .005 * (i %2), 0.24 + j * .01 +  .005 * (j %2), 0.1 +  k * .01 + .005 * (k %2));
                m_vVecState.push_back(point);
                m_vVecState.push_back(Vector3f::ZERO);
            }
        }
    }

    m_numParticles = 4500;
}

void SPHFluidSystem::build2DTestSystem()
{
    float k = 0.4;
    for (int i = 0; i < 10; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            Vector3f point(0.04 * (i + 1), 0.04 * (j + 1), k);
            m_vVecState.push_back(point);
            m_vVecState.push_back(Vector3f::ZERO);

        }
    }

    m_numParticles = 100;

}
