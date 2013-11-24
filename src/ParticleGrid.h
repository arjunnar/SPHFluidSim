/*
 * ParticleGrid.h
 *
 *  Created on: Nov 23, 2013
 *      Author: arjunnar
 */

#ifndef PARTICLEGRID_H_
#define PARTICLEGRID_H_
#include <vector>
#include <vecmath.h>
#include "tuple.h"

typedef std::vector<std::vector<std::vector<std::vector<int> > > > Grid3D;
using namespace std;

class ParticleGrid
{


public:
	ParticleGrid();
	ParticleGrid(Vector3f origin, float size, int numCellsPerDimension);
	virtual ~ParticleGrid();
	std::vector<int> getNeighborParticleIndexes(int ParticleIndex, Vector3f &particleLoc);
	void initializeGrid(std::vector<Vector3f> &particleLocations);

	// Made this public for testing purposes only
	Tuple::tuple<int, 3> getGridCoordinates(Vector3f &particleLoc);

private:
	// Instance variables
	float sideLength;
	Vector3f origin;
	Vector3f topRightCorner;
	Grid3D grid;
	int numCellsPerDimension;
	float gridSideLength;

	// Helper functions
	void baseInitGrid();
	inline bool isCoordValid(int val);
};

#endif /* PARTICLEGRID_H_ */
