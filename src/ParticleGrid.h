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
#include "kernelutilities.h"

typedef std::vector<std::vector<std::vector<std::vector<int> > > > Grid3D;

using namespace std;

class ParticleGrid
{

public:
	ParticleGrid();
    ParticleGrid(Vector3f origin, float sizeX, float sizeY, float sizeZ);
	virtual ~ParticleGrid();

	std::vector<int> getNeighborParticleIndexes(int ParticleIndex, Vector3f &particleLoc);
	void initializeGrid(std::vector<Vector3f> &particleLocations);
	Tuple::tuple<int, 3> getGridCoordinates(Vector3f &particleLoc);

    // Getters
    float getSideLengthX()
    {
        return sideLengthX;
    }

    float getSideLengthY()
    {
        return sideLengthY;
    }

    float getSideLengthZ()
    {
        return sideLengthZ;
    }

    int getNumCellsX()
    {
        return numCellsPerDimenX;
    }

    int getNumCellsY()
    {
        return numCellsPerDimenY;
    }

    int getNumCellsZ()
    {
        return numCellsPerDimenZ;
    }

private:
	// Instance variables
	Vector3f origin;
	Vector3f topRightCorner;
    int numCellsPerDimenX;
    int numCellsPerDimenY;
    int numCellsPerDimenZ;
    std::vector<std::list<int>> grid;

    float sideLengthX;
	float sideLengthY;
	float sideLengthZ;

	float gridSideLengthX;
	float gridSideLengthY;
	float gridSideLengthZ;

    map<int, Tuple::tuple<int, 3>> mapIndexesToGridCoords;

	// Helper functions
    inline bool isCoordValidX(int val);
    inline bool isCoordValidY(int val);
    inline bool isCoordValidZ(int val);
    inline int getGridIndex(int i, int j, int k);
    inline std::list<int> getGridListAt(int i, int j, int k);
};

#endif /* PARTICLEGRID_H_ */
