
#include "simpleSystem.h"

using namespace std;

SimpleSystem::SimpleSystem()
{
    m_numParticles = 1;
    this->setInitState();
}

vector<Vector3f> SimpleSystem::evalF(vector<Vector3f> state)
{
    Vector3f stateVector = state.at(0);
    Vector3f deriv(-1.0 * stateVector.y(), stateVector.x(), 0.0);

    vector<Vector3f> f;
    f.push_back(deriv);
    return f;
}

void SimpleSystem::reinitializeSystem()
{
    this->setInitState();
}

void SimpleSystem::setInitState()
{
    m_vVecState.clear();

    Vector3f initialPosition(0.5, 1.0, 0.5);
    m_vVecState.push_back(initialPosition);
}

void SimpleSystem::draw()
{
    glDisable(GL_LIGHTING);
    glColor3f(245.0/255.0, 15.0/255.0, 222.0/255.0);
    Vector3f pos = this->getState()[0];

    glPushMatrix();
    glTranslatef(pos[0], pos[1], pos[2] );
    glutSolidSphere(0.075f,10.0f,10.0f);
    glPopMatrix();
    glEnable(GL_LIGHTING);
}

