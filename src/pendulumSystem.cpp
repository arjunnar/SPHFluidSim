
#include "pendulumSystem.h"
#include "PhysicsUtilities.h"

PendulumSystem::PendulumSystem(int numParticles):ParticleSystem(numParticles)
{
	this->m_numParticles = numParticles;
	this->setInitState();
}

void PendulumSystem::setInitState()
{
    this->m_vVecState.clear();

    float y = 0.0;
    float zIncrement = 0.30;

    float z = 0.5;
    float x = 0.0;

    for (int i = 0; i < m_numParticles; ++i)
    {
        Vector3f initPos(0, y, z);
        Vector3f initVelocity(0.0, 0.0, 0.0);

        this->m_vVecState.push_back(initPos);
        this->m_vVecState.push_back(initVelocity);
        z += zIncrement;
    }
}

void PendulumSystem::reinitializeSystem()
{
    this->setInitState();
}

vector<Vector3f> PendulumSystem::evalF(vector<Vector3f> state)
{
    vector<Vector3f> derivative;

    for (int i = 0; i < m_numParticles; ++i)
    {
        Vector3f velocity;
        Vector3f totalForceOnParticle;

        if (i == 0)
        {
           velocity = Vector3f(0.0, 0.0, 0.0);
           totalForceOnParticle = Vector3f(0.0, 0.0, 0.0);
        }

        else
        {
            Vector3f gravityForce = PhysicsUtilities::getGravityForce(PendulumSystem::PARTICLE_MASS, PendulumSystem::GRAVITY_CONSTANT);
            Vector3f position = PhysicsUtilities::getPositionOfParticle(state, i);
            velocity = PhysicsUtilities::getVelocityOfParticle(state, i);
            Vector3f dragForce = PhysicsUtilities::getDragForce(PendulumSystem::DRAG_CONSTANT, velocity);

            // Calculate the total spring force on the particle
            vector<int> indexesConnectedTo = this->getParticlesConnectedTo(i);

            Vector3f totalSpringForce = Vector3f::ZERO;
            for (int connectedIndex : indexesConnectedTo)
            {
                Vector3f posConnectedTo = PhysicsUtilities::getPositionOfParticle(state, connectedIndex);
                Vector3f springForce = PhysicsUtilities::getSpringForce(position,
                                                                        posConnectedTo,
                                                                        PendulumSystem::SPRING_CONSTANT,
                                                                        PendulumSystem::SPRING_REST_LENGTH);

                totalSpringForce += springForce;
            }

            totalForceOnParticle += gravityForce + dragForce + totalSpringForce;
        }

        derivative.push_back(velocity);
        derivative.push_back(totalForceOnParticle/PendulumSystem::PARTICLE_MASS);
    }

    return derivative;
}

vector<int> PendulumSystem::getParticlesConnectedTo(int i)
{
    vector<int> indexesConnectedTo;

    if (i != 0)
    {
        indexesConnectedTo.push_back(i - 1);
    }

    if (i != this->m_numParticles - 1)
    {
        indexesConnectedTo.push_back(i + 1);
    }

    return indexesConnectedTo;
}

void PendulumSystem::draw()
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

		// Draw the springs connected to that particle
		vector<int> indexesConnectedTo = this->getParticlesConnectedTo(i);
		for (int index : indexesConnectedTo)
		{
		    Vector3f posOtherParticle = PhysicsUtilities::getPositionOfParticle(this->m_vVecState, index);
		    glBegin(GL_LINES);
		    glVertex3d(posParticle.x(), posParticle.y(), posParticle.z());
	        glVertex3d(posOtherParticle.x(), posOtherParticle.y(), posOtherParticle.z());
	        glEnd();
		}
	}

	glEnable(GL_LIGHTING);
}
