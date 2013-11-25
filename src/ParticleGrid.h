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
#include  <map>
#include <list>

typedef std::vector<std::vector<std::vector<std::vector<int> > > > Grid3D;

using namespace std;

class ParticleGrid
{

public:
	ParticleGrid();
    ParticleGrid(Vector3f origin, float size);
	virtual ~ParticleGrid();
	std::vector<int> getNeighborParticleIndexes(int ParticleIndex, Vector3f &particleLoc);
	void initializeGrid(std::vector<Vector3f> &particleLocations);

	// Made this public for testing purposes only
	Tuple::tuple<int, 3> getGridCoordinates(Vector3f &particleLoc);

    float getGridSize()
    {
        return gridSideLength;
    };

private:
	// Instance variables
	float sideLength;
	Vector3f origin;
	Vector3f topRightCorner;
    static const int NUM_CELLS_PER_DIMEN = 100;
    std::vector<std::list<int>> grid;
	float gridSideLength;
    map<int, Tuple::tuple<int, 3>> indexesToGridCoords;

	// Helper functions
	inline bool isCoordValid(int val);
    inline int getGridIndex(int i, int j, int k);
    inline std::list<int> getGridListAt(int i, int j, int k);
};

#endif /* PARTICLEGRID_H_ */
