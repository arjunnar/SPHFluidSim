/*
 * ParticleGrid.cpp
 *
 *  Created on: Nov 23, 2013
 *      Author: arjunnar
 */

#include "ParticleGrid.h"


ParticleGrid::ParticleGrid()
{
    ParticleGrid(Vector3f::ZERO, 0.6, 0.9, 0.6);
}


ParticleGrid::ParticleGrid(Vector3f origin, float sizeX, float sizeY, float sizeZ)
{
    grid = std::vector<std::list<int>>(100 * 100 * 100);
    for (int i = 0; i < grid.size(); ++i)
    {
        grid[i] = std::list<int>();
    }

    indexesToGridCoords = map<int, Tuple::tuple<int, 3>>();

	sideLengthX = sizeX;
	sideLengthY = sizeY;
	sideLengthZ = sizeZ;
	this->origin = origin;
	topRightCorner = origin + Vector3f(sideLengthX, sideLengthY, sideLengthZ);
    gridSideLengthX = sideLengthX / NUM_CELLS_PER_DIMEN;
    gridSideLengthY = sideLengthY / NUM_CELLS_PER_DIMEN;
    gridSideLengthZ = sideLengthZ / NUM_CELLS_PER_DIMEN;
}

ParticleGrid::~ParticleGrid()
{
}

void ParticleGrid::initializeGrid(std::vector<Vector3f> &particleLocations)
{
    for (int pIndex = 0; pIndex < particleLocations.size(); ++pIndex)
	{
        Vector3f loc = particleLocations[pIndex];

        Tuple::tuple<int, 3> gridCoords = getGridCoordinates(loc);
		int i = gridCoords[0];
		int j = gridCoords[1];
		int k = gridCoords[2];

		if ( !(isCoordValid(i) && isCoordValid(j) && isCoordValid(k)) )
		{
			continue;
		}

        bool oldCoordsInMap = indexesToGridCoords.find(pIndex) != indexesToGridCoords.end();

        int oldi;
        int oldj;
        int oldk;

        if (oldCoordsInMap)
        {
            Tuple::tuple<int, 3> oldGridCoords = indexesToGridCoords[pIndex];
            oldi = oldGridCoords[0];
            oldj = oldGridCoords[1];
            oldk = oldGridCoords[2];
        }

        if (oldCoordsInMap && i == oldi && j == oldj && k == oldk)
        {
            continue;
        }

        else
        {
            if (!oldCoordsInMap)
            {
                grid[getGridIndex(i, j, k)].push_back(pIndex);
            }

            else
            {
                grid[getGridIndex(oldi, oldj, oldk)].remove(pIndex);
                grid[getGridIndex(i, j, k)].push_back(pIndex);
                indexesToGridCoords.erase(pIndex);
            }

            indexesToGridCoords.insert(std::pair<int, Tuple::tuple<int, 3>>(pIndex, Tuple::tuple<int, 3>(i, j, k)));
        }
	}
}




Tuple::tuple<int, 3> ParticleGrid::getGridCoordinates(Vector3f &particleLoc)
{
    int i = (int) ( particleLoc.x() / gridSideLengthX );
    int j = (int) ( particleLoc.y() / gridSideLengthY );
    int k = (int) ( particleLoc.z() / gridSideLengthZ );
	return Tuple::tuple<int, 3>(i, j, k);
}

std::vector<int> ParticleGrid::getNeighborParticleIndexes(int particleIndex, Vector3f &particleLoc)
{
	std::vector<int> neighborParticleIndexes;
	Tuple::tuple<int, 3> gridCoordsOfParticle = getGridCoordinates(particleLoc);

	int iParticle = gridCoordsOfParticle[0];
	int jParticle = gridCoordsOfParticle[1];
	int kParticle = gridCoordsOfParticle[2];

	for (int iIncr = -1; iIncr <= 1; ++iIncr)
	{
        int iNeighbor = iParticle + iIncr;https://github.com/
		if (!isCoordValid(iNeighbor)) { continue; }

		for (int jIncr = -1; jIncr <= 1; ++jIncr)
		{
			int jNeighbor = jParticle + jIncr;
			if (!isCoordValid(jNeighbor)) { continue; }

			for (int kIncr = -1; kIncr <= 1; ++kIncr)
			{
				int kNeighbor = kParticle + kIncr;
				bool inSameCell = iNeighbor == iParticle && jNeighbor == jParticle && kNeighbor == kParticle;

				if (!isCoordValid(kNeighbor)) { continue; }

                std::list<int> neighborsInCell = getGridListAt(iNeighbor, jNeighbor, kNeighbor);

				for (int neighborIndex : neighborsInCell)
				{
					if (!inSameCell || neighborIndex != particleIndex)
                    {
						neighborParticleIndexes.push_back(neighborIndex);
					}
				}
			}
		}
	}

	return neighborParticleIndexes;
}

// Helper functions
/*
void ParticleGrid::baseInitGrid()
{
    for (int i = 0; i < numCellsPerDimension; ++i)
    {
        grid.push_back(std::vector<std::vector<std::vector<int> > >());
        for (int j = 0; j < numCellsPerDimension; ++j)
        {
            grid[i].push_back(std::vector<std::vector<int> >());
            for (int k = 0; k < numCellsPerDimension; ++k)
            {
                grid[i][j].push_back(std::vector<int>());
            }
        }
    }
}
*/

inline bool ParticleGrid::isCoordValid(int val)
{
    return 0 <= val && val < NUM_CELLS_PER_DIMEN;
}

inline int ParticleGrid::getGridIndex(int i, int j, int k)
{
    return NUM_CELLS_PER_DIMEN * NUM_CELLS_PER_DIMEN * i + NUM_CELLS_PER_DIMEN * j + k;
}
inline list<int> ParticleGrid::getGridListAt(int i, int j, int k)
{
    return grid[NUM_CELLS_PER_DIMEN * NUM_CELLS_PER_DIMEN * i + NUM_CELLS_PER_DIMEN * j + k];
}


