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

	// World space to grid coords test 1
	//Vector3f pointOne(.2, .2, .2);
	//Vector3f pointTwo(1.0, .2, .2);

	// World space to grid coords test 2
	// Vector3f pointOne(2.5, 2.5, 2.5); // Expect (6, 6, 6)
	// Vector3f pointTwo(3.5, 1.2, 3.9); // Expect (8, 3, 9)

	// Neighbor test 1
	// Vector3f pointOne(1.2, 3.4, .2); // Expect (3, 8, 0)
	// Vector3f pointTwo(1.2, 3.5, .6); // Expect (3, 8, 1)
	// Vector3f pointThree(1.3, 3.3, .3);

	// Neighbor test 2
	/*
	Vector3f pointZero(2.5, 2.5, 2.5);
	Vector3f pointOne(2.5, 2.5, 2.1);
	Vector3f pointTwo(2.5, 2.1, 2.5);
	Vector3f pointThree(2.1, 2.1, 2.1);
	Vector3f pointFour(2.5, 2.5, 3.0);
	Vector3f pointFive(3.0, 3.0, 3.0);

	vector<Vector3f> points;
	points.push_back(pointZero);
	points.push_back(pointOne);
	points.push_back(pointTwo);
	points.push_back(pointThree);
	points.push_back(pointFour);
	points.push_back(pointFive);

	particleGrid.initializeGrid(points);

	m_vVecState.push_back(pointZero);
	m_vVecState.push_back(Vector3f::ZERO);
	m_vVecState.push_back(pointOne);
	m_vVecState.push_back(Vector3f::ZERO);

	Tuple::tuple<int, 3> gridCoords1 = particleGrid.getGridCoordinates(pointZero);
	Tuple::tuple<int, 3> gridCoords2 = particleGrid.getGridCoordinates(pointOne);
	Tuple::tuple<int, 3> gridCoords3 = particleGrid.getGridCoordinates(pointTwo);
	Tuple::tuple<int, 3> gridCoords4 = particleGrid.getGridCoordinates(pointThree);
	Tuple::tuple<int, 3> gridCoords5 = particleGrid.getGridCoordinates(pointFour);
	Tuple::tuple<int, 3> gridCoords6 = particleGrid.getGridCoordinates(pointFive);


	cout << "0 - i: " << gridCoords1[0] << " j: " << gridCoords1[1] << " k: " << gridCoords1[2] << "\n";
	cout << "1 - i: " << gridCoords2[0] << " j: " << gridCoords2[1] << " k: " << gridCoords2[2] << "\n";
	cout << "2 - i: " << gridCoords3[0] << " j: " << gridCoords3[1] << " k: " << gridCoords3[2] << "\n";
	cout << "3 - i: " << gridCoords4[0] << " j: " << gridCoords4[1] << " k: " << gridCoords4[2] << "\n";
	cout << "4 - i: " << gridCoords5[0] << " j: " << gridCoords5[1] << " k: " << gridCoords5[2] << "\n";
	cout << "5 - i: " << gridCoords6[0] << " j: " << gridCoords6[1] << " k: " << gridCoords6[2] << "\n";
	vector<int> neighborsOne = particleGrid.getNeighborParticleIndexes(0, pointZero);
	cout << "Neighbors list of particles 0" << "\n";
	for (int i : neighborsOne) {
		cout << i << " \n";
	}


	vector<int> neighborsTwo = particleGrid.getNeighborParticleIndexes(1, pointOne);

	cout << "Neighbors list of particles 1" << "\n";
	for (int i : neighborsTwo) {
		cout << i << " \n";
	}

	vector<int> neighborsThree = particleGrid.getNeighborParticleIndexes(2, pointTwo);

	cout << "Neighbors list of particles 2" << "\n";
	for (int i : neighborsThree) {
		cout << i << " \n";
	}

	vector<int> neighborsFour= particleGrid.getNeighborParticleIndexes(3, pointThree);

	cout << "Neighbors list of particles 3" << "\n";
	for (int i : neighborsFour) {
		cout << i << " \n";
	}

	vector<int> neighborsFive= particleGrid.getNeighborParticleIndexes(4, pointFour);

	cout << "Neighbors list of particles 4" << "\n";
	for (int i : neighborsFive) {
		cout << i << " \n";
	}

	vector<int> neighborsSix= particleGrid.getNeighborParticleIndexes(5, pointFive);

	cout << "Neighbors list of particles 5" << "\n";
	for (int i : neighborsSix) {
		cout << i << " \n";
	}
	*/
}

GridTestSystem::~GridTestSystem()
{

}

vector<Vector3f> GridTestSystem::evalF(vector<Vector3f> state)
{
    vector<Vector3f> derivative;
    for (int i = 0; i < state.size(); i++)
    {
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

void GridTestSystem::reinitializeSystem()
{

}
