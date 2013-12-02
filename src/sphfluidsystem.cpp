#include "sphfluidsystem.h"
#include "PhysicsUtilities.h"
#include "kernelutilities.h"

using namespace std;

SPHFluidSystem::SPHFluidSystem()
{

}

SPHFluidSystem::SPHFluidSystem(float boxSizeX, float boxSizeY, float boxSizeZ)
{
    initConstants();

    build3DTestSystem();
    //build2DTestSystem();

    Vector3f origin = Vector3f::ZERO;
    particleGrid = ParticleGrid(origin, boxSizeX , boxSizeY,  boxSizeZ);
    vecParticleDensities = vector<float>();
    vecParticlePressures = vector<float>();

    vHalfStepBefore = vector<Vector3f>();
    for (int i = 0; i < m_numParticles; ++i)
    {
        vHalfStepBefore.push_back(Vector3f::ZERO);
    }
}

void SPHFluidSystem::initConstants()
{
    // Physics constants
    PARTICLE_MASS = 0.02;
    GRAVITY_CONSTANT = 6.5;
    REST_DENSITY = 920;
    GAS_CONSTANT = 1.0;
    VISCOSITY_CONSTANT = 6.0;
    TENSION_CONSTANT = 0.0004;
    TENSION_THRESHOLD = 6.0;
    SELF_DENSITY_CONSTANT = PARTICLE_MASS * KernelUtilities::polySixKernel(Vector3f::ZERO);
    SELF_LAPLACIAN_COLOR_FIELD = PARTICLE_MASS * KernelUtilities::laplacianPolySixKernel(Vector3f::ZERO);

    // Other constants
    MIN_DENSITY = 700;
    MAX_DENSITY = 5000;

    timeStep = 0.003;
}

SPHFluidSystem::~SPHFluidSystem()
{

}

void SPHFluidSystem::advanceState()
{
    vector<Vector3f> deriv = evalF(m_vVecState);
    vector<Vector3f> newState = vector<Vector3f>();

    vector<Vector3f> viPlusHalf = vector<Vector3f>();
    for (int i = 1; i < deriv.size(); i += 2)
    {
        Vector3f nextVal = vHalfStepBefore[i] + deriv[i] * timeStep;
        viPlusHalf.push_back(nextVal);
    }


    for (int i = 0; i < m_vVecState.size(); ++i)
    {
        if (i % 2 == 0)
        {
            Vector3f prevPosition = m_vVecState[i];
            Vector3f newPosition = prevPosition + vHalfStepBefore[i / 2] * timeStep;
            newState.push_back(newPosition);
        }

        else
        {
            newState.push_back(viPlusHalf[i / 2]);
        }
    }

    vHalfStepBefore = viPlusHalf;
    this->setState(newState);
}

vector<Vector3f> SPHFluidSystem::evalF(vector<Vector3f> state)
{
    vector<Vector3f> particlePositionsInState = PhysicsUtilities::getParticlePositions(state);
    particleGrid.initializeGrid(particlePositionsInState);
    calculateDensitiesAndPressures(state);

    vector<Vector3f> derivative;
    for (int particleIndex = 0; particleIndex < m_numParticles; ++particleIndex)
    {
        Vector3f positionOfParticle = PhysicsUtilities::getPositionOfParticle(state, particleIndex);
        Vector3f velocityOfParticle = PhysicsUtilities::getVelocityOfParticle(state, particleIndex);
        float densityAtParticleLoc = vecParticleDensities[particleIndex];
        float pressureAtParticleLoc = vecParticlePressures[particleIndex];


        // CALCULATE PRESSURE AND VISCOSITY FORCES
        // ALSO CALCULATE PIECES OF THE SURFACE TENSION FORCE
        Vector3f totalPressureForce = Vector3f::ZERO;
        Vector3f totalViscosityForce = Vector3f::ZERO;
        Vector3f totalGradColorField = Vector3f::ZERO;
        float totalLaplacianColorField = 0.0f;

        vector<int> neighborIndexes = particleGrid.getNeighborParticleIndexes(particleIndex, positionOfParticle);
        for (int neighborI : neighborIndexes)
        {
            Vector3f positionOfNeighbor = PhysicsUtilities::getPositionOfParticle(state, neighborI);

            float densityAtNeighborLoc = vecParticleDensities[neighborI];
            float pressureAtNeighborLoc = vecParticlePressures[neighborI];

            Vector3f spikyKernelGradForDebugging;

            Vector3f pressureContribution;
            Vector3f viscosityContribution;
            Vector3f gradColorFieldContribution;
            float laplacianColorFieldContribution;

            Vector3f rForKernel = positionOfParticle - positionOfNeighbor;

            float rEpsilon = 0.0005;
            if (rForKernel.abs() < rEpsilon || rForKernel.abs() > h)
            {

                pressureContribution = Vector3f::ZERO;
                viscosityContribution = Vector3f::ZERO;
                gradColorFieldContribution = Vector3f::ZERO;
                laplacianColorFieldContribution = 0.0f;

            }

            else
            {
                // CALCULATE PRESSURE CONTRIBUTION FROM THE NEIGHBOR
                Vector3f spikyKernelGrad = KernelUtilities::gradSpikyKernel(rForKernel);
                spikyKernelGradForDebugging = spikyKernelGrad;

                if (isNan(spikyKernelGrad))
                {
                    cout << "NAN spiky kernel gradient" << endl;
                    cout << "Spiky kernel gradient: "; DebugUtilities::printVector3f(spikyKernelGrad);
                    cout << "ri - rj magnitude: " << rForKernel.abs() << endl;
                }

                pressureContribution = PhysicsUtilities::getPressureForce( PARTICLE_MASS,
                                                                           pressureAtParticleLoc,
                                                                           pressureAtNeighborLoc,
                                                                           densityAtNeighborLoc,
                                                                           spikyKernelGrad );

                //cout << "Pressure contribution: "; DebugUtilities::printVector3f(pressureContribution);

                //cout << "Particle loc: "; DebugUtilities::printVector3f(positionOfParticle);
                //cout << "Density at particle loc: " << densityAtParticleLoc;
                //cout << "Pressure at particle loc: " << pressureAtParticleLoc << endl;

                //cout << "NeighborI: " << neighborI << endl;
                //cout << "Neighbor loc: "; DebugUtilities::printVector3f(positionOfNeighbor);
                //cout << "Density at neighbor loc: " << densityAtNeighborLoc << endl;
                float pressureAtNeighborLoc = PhysicsUtilities::getPressureAtLocation(densityAtNeighborLoc, REST_DENSITY, GAS_CONSTANT);
                //cout << "Pressure at neighbor loc: " << pressureAtNeighborLoc << endl;

                //cout << "Spiky kernel grad: "; DebugUtilities::printVector3f(spikyKernelGradForDebugging);

                // CALCULATE VISCOSITY CONTRIBUTION FROM THE NEIGHBOR
                float laplacianKernel = KernelUtilities::laplacianViscosityKernel(rForKernel);
                Vector3f velocityOfNeighbor = PhysicsUtilities::getVelocityOfParticle(state, neighborI);
                viscosityContribution = PhysicsUtilities::getViscosityForce(PARTICLE_MASS,
                                                                            VISCOSITY_CONSTANT,
                                                                            densityAtNeighborLoc,
                                                                            laplacianKernel,
                                                                            velocityOfNeighbor,
                                                                            velocityOfParticle);

                // CALCULATE GRAD COLOR FIELD CONTRIBUTION FROM NEIGHBOR
                Vector3f gradPolySixKernel = KernelUtilities::gradPolySixKernel(rForKernel);
                gradColorFieldContribution = (PARTICLE_MASS / densityAtNeighborLoc) * gradPolySixKernel; // Eq(15);

                // CALCULATE LAPLACIAN COLOR FIELD CONTRIBUTION FROM NEIGHBOR
                float laplacianPolySixKernel = KernelUtilities::laplacianPolySixKernel(rForKernel);
                laplacianColorFieldContribution = (PARTICLE_MASS / densityAtNeighborLoc) * laplacianPolySixKernel; // Eq(15)

                totalViscosityForce += viscosityContribution;
                totalPressureForce += pressureContribution;
                totalGradColorField += gradColorFieldContribution;

                //cout << "Grad color field contribution: "; DebugUtilities::printVector3f(gradColorFieldContribution);

                // For debugging
                if(isNan(viscosityContribution))
                {
                    cout << "Encountered NAN viscosity contribution" << endl;
                    cout << "Viscosity contribution: "; DebugUtilities::printVector3f(viscosityContribution);

                    cout << "Particle loc: "; DebugUtilities::printVector3f(positionOfParticle);
                    cout << "Particle velocity: "; DebugUtilities::printVector3f(velocityOfParticle);

                    cout << "NeighborI: " << neighborI << endl;
                    cout << "Neighbor loc: "; DebugUtilities::printVector3f(positionOfNeighbor);
                    Vector3f velocityNeighbor = PhysicsUtilities::getVelocityOfParticle(state, neighborI);
                    cout << "Neighbor velocity: "; DebugUtilities::printVector3f(velocityNeighbor);
                    cout << "Density at neighbor loc: " << densityAtNeighborLoc << endl;

                    float laplacianKernel = KernelUtilities::laplacianViscosityKernel(positionOfParticle - positionOfNeighbor);
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
                    cout << "Neighbor loc: "; DebugUtilities::printVector3f(positionOfNeighbor);
                    cout << "Density at neighbor loc: " << densityAtNeighborLoc << endl;
                    float pressureAtNeighborLoc = PhysicsUtilities::getPressureAtLocation(densityAtNeighborLoc, REST_DENSITY, GAS_CONSTANT);
                    cout << "Pressure at neighbor loc: " << pressureAtNeighborLoc << endl;

                    cout << "Spiky kernel grad: "; DebugUtilities::printVector3f(spikyKernelGradForDebugging);

                    assert(false); // Kill execution of the program
                }

                if (isNan(gradColorFieldContribution))
                {
                    cout << "Encountered NAN grad color field contribution" << endl;
                    cout << "Grad color field contribution: "; DebugUtilities::printVector3f(gradColorFieldContribution);
                    assert(false);
                }

                if (isNan(laplacianColorFieldContribution))
                {
                    cout << "Encountered NAN laplacian color field contribution" << endl;
                    cout << "Laplacian color field contribution: " << laplacianColorFieldContribution << endl;
                    assert(false);
                }
            }

            totalViscosityForce += viscosityContribution;
            totalPressureForce += pressureContribution;
            totalGradColorField += gradColorFieldContribution;
            totalLaplacianColorField += laplacianColorFieldContribution;

        }

        totalLaplacianColorField += (SELF_LAPLACIAN_COLOR_FIELD / densityAtParticleLoc);

        // CALCULATE TOTAL SURFACE TENSION FORCE FROM COMPONENTS
        Vector3f surfaceNormal = totalGradColorField;
        float surfaceNormalMag = surfaceNormal.abs();

        Vector3f totalSurfaceTensionForce;

        if (surfaceNormalMag > TENSION_THRESHOLD)
        {
            float constant = (-1.0 * TENSION_CONSTANT * totalLaplacianColorField) / surfaceNormalMag;
            totalSurfaceTensionForce = constant * surfaceNormal;
        }

        else
        {
            //cout << "Surface normal mag: " << surfaceNormalMag << endl;
            totalSurfaceTensionForce = Vector3f::ZERO;
        }

        //cout << "Particle: " << particleIndex << endl;
        //cout << "Total surface tension: "; DebugUtilities::printVector3f(totalSurfaceTensionForce);
        //cout << "Total viscosity: "; DebugUtilities::printVector3f(totalViscosityForce);
        //cout << "Total pressure: "; DebugUtilities::printVector3f(totalPressureForce);

        //if (totalSurfaceTensionForce.abs() > 0.0f)
        //{
            //cout << "Total surf tension force: "; DebugUtilities::printVector3f(totalSurfaceTensionForce);
            //cout << "Total grad color field: "; DebugUtilities::printVector3f(totalGradColorField);
            //cout << "Total laplacian color field: " << totalLaplacianColorField << endl;
        //}

        Vector3f totalForce = totalPressureForce + totalSurfaceTensionForce + totalViscosityForce;

        //cout << "Total force: "; DebugUtilities::printVector3f(totalForce);

        // COMPUTE TOTAL ACCELERATION OF PARTICLE


        Vector3f accelPressure = totalPressureForce / densityAtParticleLoc;
        Vector3f accelViscosity = totalViscosityForce / densityAtParticleLoc;
        Vector3f accelSurfaceTension = totalSurfaceTensionForce / densityAtParticleLoc;
        Vector3f accelGravity = PhysicsUtilities::getGravityForce(PARTICLE_MASS, GRAVITY_CONSTANT) / PARTICLE_MASS;

        Vector3f accelTotal = accelPressure + accelViscosity + accelGravity + accelSurfaceTension;

        //cout << "Accel total: "; DebugUtilities::printVector3f(accelTotal);

        derivative.push_back(velocityOfParticle);
        derivative.push_back(accelTotal);
    }


    vector<Vector3f> deriv;
    /*
    for (int i = 0; i < m_numParticles; ++i)
    {
        deriv.push_back(Vector3f::ZERO);
        deriv.push_back(Vector3f::ZERO);
    }

    return deriv;
    */
    return derivative;


}

void SPHFluidSystem::draw()
{
    for (int i = 0; i < m_numParticles; i++)
    {
        // Draw the particles
        Vector3f posParticle = PhysicsUtilities::getPositionOfParticle(m_vVecState, i);

        if (isNan(posParticle))
        {
            cout << "Encountered NAN position for particle: " << i; DebugUtilities::printVector3f(posParticle);
        }

        glPushMatrix();
        Vector3f colorParticle = getParticleColor(vecParticleDensities[i]);
        GLfloat particleColor[] = { colorParticle.x(), colorParticle.y(), colorParticle.z(), 1.0f};
        glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, particleColor );
        glTranslatef(posParticle[0], posParticle[1], posParticle[2] );
        glutSolidSphere(0.02, 10, 10);
        glPopMatrix();
    }

    //cout << "Num nans: " << numNanPositions << endl;
}

void SPHFluidSystem::reinitializeSystem()
{

}

// Helper functions
void SPHFluidSystem::calculateDensitiesAndPressures(vector<Vector3f> &state)
{
    vecParticleDensities = vector<float>();
    vecParticlePressures = vector<float>();

    // Get an idea of the min and max density values so we can color the particles based on density.
    float minDensity = INFINITY;
    float maxDensity = -1 * INFINITY;

    for (int i = 0; i < m_numParticles; ++i)
    {
        float density = 0;
        Vector3f pos = PhysicsUtilities::getPositionOfParticle(state, i);
        vector<int> neighborIndexes = particleGrid.getNeighborParticleIndexes(i, pos);

        for (int neighborI : neighborIndexes)
        {
            Vector3f neighborPos = PhysicsUtilities::getPositionOfParticle(state, neighborI);
            density += PARTICLE_MASS * KernelUtilities::polySixKernel(pos - neighborPos);
        }

        density += SELF_DENSITY_CONSTANT;
        vecParticleDensities.push_back(density);

        minDensity = min(density, minDensity);
        maxDensity = max(density, maxDensity);

        float pressure = PhysicsUtilities::getPressureAtLocation(density, REST_DENSITY, GAS_CONSTANT);
        vecParticlePressures.push_back(pressure);
    }

    //cout << "Min density: " << minDensity << endl;
    //cout << "Max density: " << maxDensity << endl;
}

Vector3f SPHFluidSystem::getParticleColor(float densityAtParticleLoc)
{
    densityAtParticleLoc = max(MIN_DENSITY, densityAtParticleLoc);
    densityAtParticleLoc = min(MAX_DENSITY, densityAtParticleLoc);

    float percentIntoInterval = (densityAtParticleLoc - MIN_DENSITY) / (MAX_DENSITY - MIN_DENSITY);
    return Vector3f(percentIntoInterval, 0.0f,  1 - percentIntoInterval);


}

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
    Vector3f pos1(0.3, 0.35, 0.1);
    Vector3f pos2(0.1, 0.35, 0.1);

    m_vVecState.push_back(pos1);
    m_vVecState.push_back(Vector3f::ZERO);
    m_vVecState.push_back(pos2);
    m_vVecState.push_back(Vector3f::ZERO);

    m_numParticles = 2;
}

void SPHFluidSystem::buildTwoParticleSystemNeighbors()
{
    Vector3f pos1(0.3, 0.35, 0.1);
    Vector3f pos2(0.33, 0.35, 0.1);

    m_vVecState.push_back(pos1);
    m_vVecState.push_back(Vector3f::ZERO);
    m_vVecState.push_back(pos2);
    m_vVecState.push_back(Vector3f::ZERO);

    m_numParticles = 2;
}

/*
void SPHFluidSystem::build3DTestSystem()
{
    for (int k = 0; k < 10; k++)
    {
        for (int i = 0; i < 10; i++ )
        {
            for (int j = 0; j < 10; j++) {
                //Vector3f point(0.05 + .02 * i + .005 * (), 0.24 + j * .01 +  .005 * (j %2), 0.1 +  k * .01 + .005 * (k %2));
                Vector3f point(0.05 + 0.03 * i, 0.05 + 0.03 * j, 0.1 + k * 0.03 );
                m_vVecState.push_back(point);
                m_vVecState.push_back(Vector3f::ZERO);
            }
        }
    }

    m_numParticles = 1000;
}
*/

void SPHFluidSystem::build2DTestSystem()
{
    float k = 0.2;
    for (int i = 0; i < 20; ++i)
    {
        for (int j = 0; j < 20; ++j)
        {
            Vector3f point(0.02 * (i + 1), 0.02 * (j + 1), k);
            m_vVecState.push_back(point);
            m_vVecState.push_back(Vector3f::ZERO);

        }
    }

    m_numParticles = 400;

}

void SPHFluidSystem::build2DTestSystem2()
{
    for (int k = 0; k < 15; +k)
    {
        for (int i = 0; i < 15; ++i)
        {
            for (int j = 0; j < 15; ++j)
            {
                Vector3f point(0.03 * (i + 1), 0.03 * (j + 1), 0.03 + k * 0.03);
                m_vVecState.push_back(point);
                m_vVecState.push_back(Vector3f::ZERO);

            }
        }
    }
    /*
    TENSION_CONSTANT = 0.1;
    TENSION_THRESHOLD = 6.0;
    GAS_CONSTANT = 2.0;
    VISCOSITY_CONSTANT = 2.6;
    REST_DENSITY = 300;
    */
    m_numParticles = 3375;
}

void SPHFluidSystem::build3DTestSystem()
{
    for (int k = 0; k < 10; k++)
    {
        for (int i = 0; i < 10; i++ )
        {
            for (int j = 0; j < 10; j++) {
                //Vector3f point(0.05 + .02 * i + .005 * (), 0.24 + j * .01 +  .005 * (j %2), 0.1 +  k * .01 + .005 * (k %2));
                Vector3f point(0.05 + 0.03 * i, 0.05 + 0.03 * j, 0.05 + k * 0.03 );
                m_vVecState.push_back(point);
                m_vVecState.push_back(Vector3f::ZERO);
            }
        }
    }

    m_numParticles = m_vVecState.size() / 2;
}





