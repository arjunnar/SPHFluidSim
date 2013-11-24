/*
 * GridTestSystem.h
 *
 *  Created on: Nov 23, 2013
 *      Author: lusmyk
 */

#ifndef GRIDTESTSYSTEM_H_
#define GRIDTESTSYSTEM_H_


#include "extra.h"
#include <vector>
#include "ParticleGrid.h"
#include "particleSystem.h"

class GridTestSystem: public ParticleSystem
{

public:
	GridTestSystem(int numParticles);
    virtual ~GridTestSystem();
	vector<Vector3f> evalF(vector<Vector3f> state);
	void draw();
	void reinitializeSystem();
private:
    ParticleGrid particleGrid;



};





#endif /* GRIDTESTSYSTEM_H_ */
