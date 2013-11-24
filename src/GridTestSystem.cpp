/*
 * GridTestSystem.cpp
 *
 *  Created on: Nov 23, 2013
 *      Author: lusmyk
 */

#include "GridTestSystem.h"

#include "PhysicsUtilities.h"

GridTestSystem::GridTestSystem(int numParticles) : ParticleSystem(numParticles)
{
	Vector3f origin(0,0,0);
	float gridSize = 4.0;
	int numCellsPerDimension = 10;
	particleGrid = ParticleGrid(origin, gridSize, numCellsPerDimension);
	//test 0.2, 0.2, 0.2
	Vector3f pointOne(.2, .2, .2);
	Vector3f pointTwo(1.0, .2, .2);

	vector<Vector3f> points;
	points.push_back(pointOne);
	points.push_back(pointTwo);
	particleGrid.initializeGrid(points);

	m_vVecState.push_back(pointOne);
	m_vVecState.push_back(Vector3f::ZERO);
	m_vVecState.push_back(pointTwo);
	m_vVecState.push_back(Vector3f::ZERO);

}

GridTestSystem::~GridTestSystem() {

}


vector<Vector3f> GridTestSystem::evalF(vector<Vector3f> state)
{
    vector<Vector3f> derivative;
    for (int i = 0; i < state.size(); i++)  {
    	derivative.push_back(Vector3f::ZERO);
    }

    return derivative;
}


void GridTestSystem::draw()
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

void GridTestSystem::reinitializeSystem() {

}
