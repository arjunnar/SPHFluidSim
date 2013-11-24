/*
 * ParticleGrid.cpp
 *
 *  Created on: Nov 23, 2013
 *      Author: arjunnar
 */

#include "ParticleGrid.h"


ParticleGrid::ParticleGrid()
{
	ParticleGrid(Vector3f::ZERO, 4, 100);
}


ParticleGrid::ParticleGrid(Vector3f origin, float size, int numCellsPerDimension)
{
	sideLength = size;
	this->origin = origin;
	topRightCorner = origin + Vector3f(size, size, size);
	this->numCellsPerDimension = numCellsPerDimension;
	this->gridSideLength = (sideLength * 1.0) / numCellsPerDimension;
	baseInitGrid();
}

ParticleGrid::~ParticleGrid()
{
}

void ParticleGrid::initializeGrid(std::vector<Vector3f> &particleLocations)
{
	grid = Grid3D();
	baseInitGrid();
	for (int index = 0; index < particleLocations.size(); ++index)
	{
		Vector3f loc = particleLocations[index];
		Tuple::tuple<int, 3> gridCoords = getGridCoordinates(loc);
		int i = gridCoords[0];
		int j = gridCoords[1];
		int k = gridCoords[2];

		grid[i][j][k].push_back(index);
	}
}

Tuple::tuple<int, 3> ParticleGrid::getGridCoordinates(Vector3f &particleLoc)
{
	int i = (int) ( particleLoc.x() / gridSideLength );
	int j = (int) ( particleLoc.y() / gridSideLength );
	int k = (int) ( particleLoc.z() / gridSideLength );
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
		int iNeighbor = iParticle + iIncr;
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

				std::vector<int> neighborsInCell = grid[iNeighbor][jNeighbor][kNeighbor];

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

inline bool ParticleGrid::isCoordValid(int val)
{
	return 0 <= val && val < numCellsPerDimension;
}


