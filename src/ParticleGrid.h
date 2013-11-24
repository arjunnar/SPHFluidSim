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

class ParticleGrid
{


public:
	ParticleGrid(Vector3f origin, int size, int gridSideLength);
	virtual ~ParticleGrid();
	std::vector<int> getNeighborParticles(Vector3f &particleLoc);
	void initializeGrid(std::vector<Vector3f> &particleLocations);

private:
	int sideLength;
	Vector3f origin;
	Vector3f topRightCorner;
	Grid3D grid;
	int numCellsPerDimension;
	float gridSideLength;
	void initGrid();
	tuple<int, 3> getGridCoordinates(Vector3f &particleLoc);
};

#endif /* PARTICLEGRID_H_ */