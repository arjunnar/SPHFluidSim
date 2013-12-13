#include "sphfluidsystem.h"
#include "PhysicsUtilities.h"
#include "kernelutilities.h"

using namespace std;

SPHFluidSystem::SPHFluidSystem()
{

}

SPHFluidSystem::SPHFluidSystem(float boxSizeX, float boxSizeY, float boxSizeZ, FluidSystemType fluidSystemType)
{
    initConstants();

    typeOfSystem = fluidSystemType;
    m_numParticles = 0;

    switch (typeOfSystem)
    {
         case TwoDensitySystem2D:
            build2DTestSystemTwoDensities();
            break;

         case SystemSimple2D:
            build2DTestSystemSimple();
            break;

         case SystemLarge2D:
            build2DTestSystemLarge();
            break;

         case SystemSimple3D:
            build3DTestSystemSimple();
            break;

         case SystemLarge3D:
            build3DTestSystemLarge();
            break;

         case System2DEmitter:
            m_numParticles = 0;
            break;

         case SystemDroplet3D:
            m_numParticles = 0;
            addSphereOfParticles(Vector3f(0.35, 0.7, 0.35), 0.2);
            break;

         case System3DWaves:
            build3DWaveSystem();
            break;
    }

    Vector3f origin = Vector3f::ZERO;
    particleGrid = ParticleGrid(origin, boxSizeX , boxSizeY,  boxSizeZ);
    vecParticleDensities = vector<float>();
    vecParticlePressures = vector<float>();

    locOfCannon = Vector3f(boxSizeX/5.0f, boxSizeY/2.5f, 0.2);
    locOfCannon2 = Vector3f(boxSizeX * 4.0/5.0, boxSizeY/1.9f, 0.2);
    locOfCannon3 = Vector3f(boxSizeX * 2.5 /5.0f, boxSizeY/3.5f, 0.2);
    directionOfGravity = Vector3f(0.0f, -1.0f, 0.0f);

    renderWithMarchingCubes = false;
    rotEnabled = false;
}

void SPHFluidSystem::initConstants()
{
    // Physics constants
    PARTICLE_MASS = 0.02;
    GRAVITY_CONSTANT = 6.5;
    REST_DENSITY = 920;
    REST_DENSITY_FIRST = 920;
    REST_DENSITY_SECOND = 300;
    GAS_CONSTANT = 1.0;
    VISCOSITY_CONSTANT = 6.0;
    TENSION_CONSTANT = 0.0004;
    TENSION_THRESHOLD = 6.0;
    SELF_DENSITY_CONSTANT = PARTICLE_MASS * KernelUtilities::polySixKernel(Vector3f::ZERO);
    SELF_LAPLACIAN_COLOR_FIELD = PARTICLE_MASS * KernelUtilities::laplacianPolySixKernel(Vector3f::ZERO);

    // Other constants
    MIN_DENSITY = 700;
    MAX_DENSITY = 5000;


    emitAngleIncrement = M_PI/12;
    emitAngleIncrement2 = M_PI/24;
    angleToEmit = emitAngleIncrement;
    angleToEmit2 = emitAngleIncrement2;
    emitVelocityConstant = 2.7f;
}

SPHFluidSystem::~SPHFluidSystem()
{

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

        Vector3f rotForce(-1.0 * (positionOfParticle.z() - (particleGrid.getSideLengthZ() / 2.0)),
                          0.0f,
                          positionOfParticle.x() - particleGrid.getSideLengthX() / 2.0);

        Vector3f totalForce = totalPressureForce + totalSurfaceTensionForce + totalViscosityForce;

        //cout << "Total force: "; DebugUtilities::printVector3f(totalForce);

        // COMPUTE TOTAL ACCELERATION OF PARTICLE


        Vector3f accelPressure = totalPressureForce / densityAtParticleLoc;
        Vector3f accelViscosity = totalViscosityForce / densityAtParticleLoc;
        Vector3f accelSurfaceTension = totalSurfaceTensionForce / densityAtParticleLoc;
        Vector3f accelGravity = GRAVITY_CONSTANT * directionOfGravity;

        Vector3f accelTotal = accelPressure + accelViscosity + accelGravity + accelSurfaceTension;

        if (rotEnabled)
        {
            accelTotal += rotSpeed * rotForce / PARTICLE_MASS;
        }

        //cout << "Accel total: "; DebugUtilities::printVector3f(accelTotal);

        derivative.push_back(velocityOfParticle);
        derivative.push_back(accelTotal);
    }

    /*
    vector<Vector3f> deriv;

    for (int i = 0; i < m_numParticles; ++i)
    {
        deriv.push_back(Vector3f::ZERO);
        deriv.push_back(Vector3f::ZERO);
    }

    return deriv;
    */

    return derivative;
}/*
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

void SPHFluidSystem::drawFace(TRIANGLE face, bool reverse)
{
    //cout << "Trying to draw face " << endl;
    mpVector pos1 = face.p[0];
    mpVector pos2 = face.p[1];
    mpVector pos3 = face.p[2];
    mpVector normal1 = face.norm[0];
    mpVector normal2 = face.norm[1];
    mpVector normal3 = face.norm[2];

    if (reverse)
    {
        pos1 = face.p[2];
        pos3 = face.p[0];
        normal1 = face.norm[2];
        normal3 = face.norm[0];
    }

    glBegin(GL_TRIANGLES);

    // First normal and first vertex
    glNormal3d(normal1.x, normal1.y, normal1.z);
    glVertex3d(pos1.x, pos1.y, pos1.z);

    // Second normal and second vertex
    glNormal3d(normal2.x, normal2.y, normal2.z);
    glVertex3d(pos2.x, pos2.y, pos2.z);

    // Third normal and third vertex
    glNormal3d(normal3.x, normal3.y, normal3.z);
    glVertex3d(pos3.x, pos3.y, pos3.z);

    glEnd();
}

void SPHFluidSystem::setMovingHorizontally(bool right)
{
    /*
    Vector3f velocity = right ? Vector3f(3.0, 0.0, 0.0) : Vector3f(-3.0, 0.0, 0.0);
    for (int particleIndex = 0; particleIndex < m_numParticles; ++particleIndex)
    {
        PhysicsUtilities::setVelocityOfParticle(m_vVecState, particleIndex, velocity);
    }
    */

    directionOfGravity = Vector3f(-1.0, -1.0, 0);
}


float SPHFluidSystem::calcDensity(Vector3f pos)
{
    float density = 0;

    vector<int> neighborIndexes = particleGrid.getNeighborParticleIndexes(-1, pos);
   // cout << "neighbor indexes size: " << neighborIndexes.size() << endl;
    for (int neighborI : neighborIndexes)
    {
        Vector3f neighborPos = PhysicsUtilities::getPositionOfParticle(m_vVecState, neighborI);
        density += PARTICLE_MASS * KernelUtilities::polySixKernel(pos - neighborPos);
    }

    density += SELF_DENSITY_CONSTANT;
    return density;
}




void SPHFluidSystem::draw()
{
    if (renderWithMarchingCubes)
    {
        int numTrianglesCreated;

        int nX = particleGrid.getSideLengthX() / 0.01;
        int nY = particleGrid.getSideLengthY() / 0.01;
        int nZ = particleGrid.getSideLengthZ() / 0.01;
        float stepSize = 0.01;

        mp4Vector *mcPoints = new mp4Vector[(nX+1)*(nY+1)*(nZ+1)];

        for(int i=0; i < nX+1; i++)
            for(int j=0; j < nY+1; j++)
                for(int k=0; k < nZ+1; k++)
                {
                    mp4Vector vert(i*stepSize, j*stepSize, k*stepSize, 0);
                    vert.val = calcDensity(Vector3f(vert.x, vert.y, vert.z));
                    //cout << vert.val << endl;
                    mcPoints[i*(nY+1)*(nZ+1) + j*(nZ+1) + k] = vert;
                }

        TRIANGLE* triangleMeshPointer = MarchingCubes(nX, nY, nZ, 1.0f, 1.0f, 1.0f, 1000, mcPoints, numTrianglesCreated);
        //MarchingCubesCross(minX, maxX, minY, maxY, minZ, maxZ, numCellsInDimenX, numCellsInDimenY, numCellsInDimenZ, MIN_DENSITY, &calcDensity, numTrianglesCreated);

        //cout << "Num triangles created: " << numTrianglesCreated << endl;


        GLfloat particleColor[] = { 0.0f, 1.0f, 1.0f, 1.0f, };
        glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, particleColor );


        for (int triangleI = 0; triangleI < numTrianglesCreated; ++triangleI)
        {
            TRIANGLE face = *triangleMeshPointer;
            drawFace(face, false);
            drawFace(face, true);
            ++triangleMeshPointer;
            //cout << "Drawing triangle num: " << triangleI + 1 << endl;
        }
    }

    else
    {

    for (int i = 0; i < m_numParticles; i++)
        {
            // Draw the particles
            Vector3f posParticle = PhysicsUtilities::getPositionOfParticle(m_vVecState, i);

            if (isNan(posParticle))
            {
                cout << "Encountered NAN position for particle: " << i; DebugUtilities::printVector3f(posParticle);
            }


            if (typeOfSystem == FluidSystemType::TwoDensitySystem2D)
            {
                glPushMatrix();
                glDisable(GL_LIGHTING);
                glColor3f(0.0, 0.0, 1.0);
                if (i >= 200)
                {
                    //GLfloat particleColor[] = { 0.0, 1.0, 1.0, 1.0f};
                    //glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, particleColor );
                    glColor3f(1.0f, 1.0f, 0.0f);
                }

                else
                {

                    //GLfloat particleColor[] = { 0.0, 0.0, 1.0, 1.0f};
                    //glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, particleColor );
                    glColor3f(0.0, 0.0, 1.0);
                }


                glTranslatef(posParticle[0], posParticle[1], posParticle[2] );
                glutSolidSphere(0.02, 10, 10);
                glEnable(GL_LIGHTING);
                glPopMatrix();
            }

            else if (typeOfSystem == FluidSystemType::System2DEmitter)
            {
                glPushMatrix();
                glDisable(GL_LIGHTING);
                glTranslatef(posParticle[0], posParticle[1], posParticle[2] );

                if (i % 3 == 0)
                {
                    glColor3f(0.0, 1.0f, 1.0f);
                }

                else if (i % 3 == 1)
                {
                    glColor3f(0.0f, 1.0f, 0.0f);
                }

                else
                {
                    glColor3f(1.0f, 0.0f, 127.0/255.0f);
                }

                glutSolidSphere(0.02, 10, 10);
                glEnable(GL_LIGHTING);
                glPopMatrix();
            }

            else
            {
                glPushMatrix();
                Vector3f colorParticle = getParticleColor(vecParticleDensities[i]);
                GLfloat particleColor[] = { colorParticle.x(), colorParticle.y(), colorParticle.z(), 1.0f};
                glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, particleColor );

                glTranslatef(posParticle[0], posParticle[1], posParticle[2] );
                glutSolidSphere(0.02, 10, 10);
                glPopMatrix();
            }
        }
    }
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

        if ( typeOfSystem == FluidSystemType::TwoDensitySystem2D)
        {

            if (i >= 200 )
            {
                REST_DENSITY = REST_DENSITY_SECOND;
            }

            else
            {
                REST_DENSITY = REST_DENSITY_FIRST;
            }
        }

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

void SPHFluidSystem::build2DTestSystemSimple()
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

void SPHFluidSystem::build2DTestSystemLarge()
{
    float k = 0.2;
    for (int i = 0; i < 20; ++i)
    {
        for (int j = 0; j < 40; ++j)
        {
            Vector3f point(0.02 * (i + 1), 0.02 * (j + 1), k);
            m_vVecState.push_back(point);
            m_vVecState.push_back(Vector3f::ZERO);

        }
    }

    m_numParticles = 1000;
}

void SPHFluidSystem::build3DTestSystemSimple()
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

void SPHFluidSystem::build3DTestSystemLarge()
{
    for (int k = 0; k < 15; k++)
    {
        for (int i = 0; i < 15; i++ )
        {
            for (int j = 0; j < 40; j++) {
                //Vector3f point(0.05 + .02 * i + .005 * (), 0.24 + j * .01 +  .005 * (j %2), 0.1 +  k * .01 + .005 * (k %2));
                Vector3f point(0.075 + 0.03 * i, 0.10 + 0.03 * j, 0.075 + k * 0.03 );
                m_vVecState.push_back(point);
                m_vVecState.push_back(Vector3f::ZERO);
            }
        }
    }

    m_numParticles = m_vVecState.size() / 2;
}

void SPHFluidSystem::build2DTestSystemTwoDensities()
{
    float k = 0.2;
    for (int i = 0; i < 20; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            Vector3f point(0.02 * (i + 1), 0.9 + 0.02 * (j + 1), k);
            m_vVecState.push_back(point);
            m_vVecState.push_back(Vector3f::ZERO);

        }
    }

    for (int i = 0; i < 20; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            Vector3f point(0.02 * (i + 1), 0.3 + 0.02 * (j + 1), k);
            m_vVecState.push_back(point);
            m_vVecState.push_back(Vector3f::ZERO);

        }
    }

    m_numParticles = 400;
}

void SPHFluidSystem::build3DWaveSystem()
{
    for (int k = 0; k < 10; k++)
    {
        for (int i = 0; i < 39; i++ )
        {
            for (int j = 0; j < 9; j++) {
                Vector3f point(0.02 + 0.03 * i, 0.1 + 0.03 * j, 0.10 + k * 0.03 );
                m_vVecState.push_back(point);
                m_vVecState.push_back(Vector3f::ZERO);
            }
        }
    }

    m_numParticles = m_vVecState.size() / 2;
}

void SPHFluidSystem::addCircleOfParticles(Vector3f origin, float stepSize, float radius)
{
    float angle = 0.0f;

    while (angle <= 2 * M_PI)
    {
        float x = radius * cos(angle) + origin.x();
        float z = radius * sin(angle) + origin.z();
        m_vVecState.push_back(Vector3f(x, origin.y(), z));
        m_vVecState.push_back(Vector3f::ZERO);
        ++m_numParticles;
        angle += stepSize;
    }
}

void SPHFluidSystem::addDiskOfParticles(Vector3f origin, float radiusOfDisk)
{
    float stepSize = 0.05;

    m_vVecState.push_back(Vector3f(origin.x(), origin.y(), origin.z()));
    m_vVecState.push_back(Vector3f::ZERO);
    ++m_numParticles;

    float angularStep = M_PI/10;
    for (float radius = stepSize; radius <= radiusOfDisk; radius += stepSize)
    {
        addCircleOfParticles(origin, angularStep, radius);
    }
}

void SPHFluidSystem::addSphereOfParticles(Vector3f origin, float radiusOfSphere)
{
    float yStep = 0.05;
    for (float yDist = -1.0 * radiusOfSphere; yDist <= radiusOfSphere; yDist += yStep)
    {
        float yVal = origin.y() + yDist;

        Vector3f originOfDisk(origin.x(), yVal, origin.z());

        float radiusOfDisk = sqrt(pow(radiusOfSphere, 2) - pow(yDist, 2));

        addDiskOfParticles(originOfDisk, radiusOfDisk);
    }
}

void SPHFluidSystem::emitParticle()
{
    float xVel = emitVelocityConstant * cos(angleToEmit);
    float yVel = emitVelocityConstant * sin(angleToEmit);

    float xVel2 = emitVelocityConstant * cos(angleToEmit2);
    float yVel2 = emitVelocityConstant * sin(angleToEmit2);

    m_vVecState.push_back(locOfCannon);
    m_vVecState.push_back(Vector3f(xVel, yVel, 0.0f));

    m_vVecState.push_back(locOfCannon2);
    m_vVecState.push_back(Vector3f(-1.0 * xVel, yVel, 0.0f));

    m_vVecState.push_back(locOfCannon3);
    m_vVecState.push_back(Vector3f(xVel2, yVel2, 0.0f));

    m_numParticles += 3;

    angleToEmit += emitAngleIncrement;
    angleToEmit2 += emitAngleIncrement2;

    if (angleToEmit >= M_PI)
    {
        angleToEmit = M_PI - emitAngleIncrement;
        emitAngleIncrement *= -1;
    }

    if (angleToEmit <= 0.0f)
    {
        emitAngleIncrement *= -1;
        angleToEmit = emitAngleIncrement;
    }

    if (angleToEmit2 >= M_PI)
    {
        angleToEmit2 = M_PI - emitAngleIncrement2;
        emitAngleIncrement2 *= -1;
    }

    if (angleToEmit2 <= 0.0f)
    {
        emitAngleIncrement2 *= -1;
        angleToEmit2 = emitAngleIncrement2;
    }
}


