#ifndef CLOTHSYSTEM_H
#define CLOTHSYSTEM_H

#include "extra.h"
#include <vector>
#include "PhysicsUtilities.h"
#include "particleSystem.h"

class ClothSystem: public ParticleSystem
{

public:
	ClothSystem(int height, int width);
	vector<Vector3f> evalF(vector<Vector3f> state);
	void reinitializeSystem();
	void applyVelocityToCorners(bool negative);
	void stopClothMovement();
	Vector3f getTopCornerPosition(bool topLeft);

	void setWind(bool windOn)
    {
        this->windOn = windOn;
    }

	bool getWind()
	{
	    return this->windOn;
	}

	void setDrawShaded(bool shadeCloth)
	{
	    this->shadeCloth = shadeCloth;
	}

	bool getDrawShaded()
	{
	    return this->shadeCloth;
	}

	void draw();


private:
	int height;
	int width;
	bool windOn;
	bool shadeCloth;
	vector<Vector3f> normals;
	vector<vector<int>> faces;

	static constexpr float GRAVITY_CONSTANT = 1.0f;
	static constexpr float PARTICLE_MASS = 0.1f;
	static constexpr float DRAG_CONSTANT = 0.035f;

	static constexpr float STRUC_SPRING_CONSTANT = 10.0f;
	static constexpr float STRUC_SPRING_REST_LENGTH = 0.3f;

	static constexpr float SHEAR_SPRING_CONSTANT = 2.5f;
	static constexpr float SHEAR_SPRING_REST_LENGTH = 0.5f;

	static constexpr float FLEXION_SPRING_CONSTANT = 2.5f;
	static constexpr float FLEXION_SPRING_REST_LENGTH = 0.6f;

	static const int NULL_INDEX = -1;

	void setInitState();

	int getIndex(int i, int j);
	vector<int> getIndexesForStructuralSprings(int i, int j);
	vector<int> getIndexesForShearSprings(int i, int j);
	vector<int> getIndexesForFlexionSprings(int i, int j);

	void setVelocityOfCorners(Vector3f velocity);

	bool isTopCorner(int i, int j);

	// For drawing the shaded cloth
	void extractVertexPositions();
	void calculateNormals();
	void calculateFaces();
	Vector3f getVectorBetweenParticles(int index1, int index2);
	void drawShadedCloth();
	void drawFace(vector<int> faceVector, bool negateNormals);

	// For drawing the wireframe
	void drawParticles();
	void drawAllSprings();
	void drawLinesForSprings(Vector3f posParticles, vector<int> connectedIndexes);
};


#endif
