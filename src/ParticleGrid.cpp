/*
 * ParticleGrid.cpp
 *
 *  Created on: Nov 23, 2013
 *      Author: arjunnar
 */

#include "ParticleGrid.h"

ParticleGrid::ParticleGrid(Vector3f origin, int size, int numCellsPerDimension)
{
	sideLength = size;
	this->origin = origin;
	topRightCorner = origin + Vector3f(size, size, size);
	this->numCellsPerDimension = numCellsPerDimension;
	this->gridSideLength = (sideLength * 1.0) / numCellsPerDimension;
}

ParticleGrid::~ParticleGrid()
{
}

void ParticleGrid::initializeGrid(std::vector<Vector3f> &particleLocations)
{
	grid = Grid3D();
	initGrid();
	for (int index = 0; index < particleLocations.size(); ++index)
	{
		Vector3f loc = particleLocations[index];
		tuple<int, 3> gridCoords = getGridCoordinates(loc);
		int i = gridCoords[0];
		int j = gridCoords[1];
		int k = gridCoords[2];

		grid[i][j][k].push_back(index);
	}
}

std::vector<int> ParticleGrid::getNeighborParticles(Vector3f &particleLoc)
{
	tuple<int, 3> gridCoords = getGridCoordinates(particleLoc);
	int i = gridCoords[0];
	int j = gridCoords[1];
	int k = gridCoords[2];
	return grid[i][j][k];
}

// Helper functions
void ParticleGrid::initGrid()
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

tuple<int, 3> ParticleGrid::getGridCoordinates(Vector3f &particleLoc)
{
	int i = (int) ( particleLoc.x() / gridSideLength );
	int j = (int) ( particleLoc.y() / gridSideLength );
	int k = (int) ( particleLoc.z() / gridSideLength );
	return tuple<int, 3>(i, j, k);
}


