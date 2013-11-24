#include "ClothSystem.h"
#include "PhysicsUtilities.h"
#include <assert.h>
#include <math.h>

ClothSystem::ClothSystem(int height, int width)
{
    this->height = height;
    this->width = width;
    this->m_numParticles = height * width;
    this->windOn = false;
    this->shadeCloth = false;
    this->calculateFaces();
    this->setInitState();
}

void ClothSystem::applyVelocityToCorners(bool negative)
{
    // Move back and forth in the z direction.
    Vector3f velocityOfCorners(0, 0, (negative ? -1.0 : 1.0) * 0.75);
    this->setVelocityOfCorners(velocityOfCorners);
}

void ClothSystem::stopClothMovement()
{
    this->setVelocityOfCorners(Vector3f::ZERO);
}

void ClothSystem::setVelocityOfCorners(Vector3f velocity)
{
    int cornerIndex1 = this->getIndex(0, 0);
    int cornerIndex2 = this->getIndex(0, this->width - 1);

    PhysicsUtilities::setVelocityOfParticle(this->m_vVecState, cornerIndex1, velocity);
    PhysicsUtilities::setVelocityOfParticle(this->m_vVecState, cornerIndex2, velocity);
}

Vector3f ClothSystem::getTopCornerPosition(bool topLeft)
{
    int cornerParticleNum = topLeft ? this->getIndex(0, 0) : this->getIndex(0, this->width - 1);
    return PhysicsUtilities::getPositionOfParticle(this->m_vVecState, cornerParticleNum);
}

vector<Vector3f> ClothSystem::evalF(vector<Vector3f> state)
{
    vector<Vector3f> derivative;

    Vector3f gravityForce = PhysicsUtilities::getGravityForce(ClothSystem::PARTICLE_MASS, ClothSystem::GRAVITY_CONSTANT);

    // Include a random force in the Z direction
    Vector3f windForce = this->windOn ? Vector3f(0.0, 0.0, -0.30 * (rand() % 2 + 0.5)) : Vector3f::ZERO;

    for (int i = 0; i < this->height; ++i)
    {
        for (int j = 0; j < this->width; ++j)
        {
            int indexOfParticle = this->getIndex(i, j);
            Vector3f velocityOfParticle = PhysicsUtilities::getVelocityOfParticle(state, indexOfParticle);
            Vector3f totalForce;

            if (this->isTopCorner(i, j))
            {
                totalForce = Vector3f::ZERO;
            }

            else
            {
                Vector3f posOfParticle = PhysicsUtilities::getPositionOfParticle(state, indexOfParticle);

                velocityOfParticle = PhysicsUtilities::getVelocityOfParticle(state, indexOfParticle);

                Vector3f dragForce = PhysicsUtilities::getDragForce(ClothSystem::DRAG_CONSTANT, velocityOfParticle);

                Vector3f totalStructuralForce = Vector3f::ZERO;
                Vector3f totalShearForce = Vector3f::ZERO;
                Vector3f totalFlexionForce = Vector3f::ZERO;


                // Calculate the force on the particle from structural springs
                vector<int> structuralSpringIndexes = this->getIndexesForStructuralSprings(i, j);
                for (int otherParticleIndex : structuralSpringIndexes)
                {
                    Vector3f posOfOtherParticle = PhysicsUtilities::getPositionOfParticle(state, otherParticleIndex);
                    Vector3f structuralForce = PhysicsUtilities::getSpringForce(posOfParticle,
                                                                                posOfOtherParticle,
                                                                                ClothSystem::STRUC_SPRING_CONSTANT,
                                                                                ClothSystem::STRUC_SPRING_REST_LENGTH);
                    totalStructuralForce += structuralForce;
                }

                // Calculate the force on the particle from shear springs
                vector<int> shearSpringIndexes = this->getIndexesForShearSprings(i, j);
                for (int otherParticleIndex : shearSpringIndexes)
                {
                    Vector3f posOtherParticle = PhysicsUtilities::getPositionOfParticle(state, otherParticleIndex);
                    Vector3f shearForce = PhysicsUtilities::getSpringForce(posOfParticle,
                                                                           posOtherParticle,
                                                                           ClothSystem::SHEAR_SPRING_CONSTANT,
                                                                           ClothSystem::SHEAR_SPRING_REST_LENGTH);
                    totalShearForce += shearForce;
                }

                // Calculate the force on the particle from flexion springs
                vector<int> flexionSpringIndexes = this->getIndexesForFlexionSprings(i, j);
                for (int otherParticleIndex : flexionSpringIndexes)
                {
                    Vector3f posOtherParticle = PhysicsUtilities::getPositionOfParticle(state, otherParticleIndex);
                    Vector3f flexionForce = PhysicsUtilities::getSpringForce(posOfParticle,
                                                                             posOtherParticle,
                                                                             ClothSystem::FLEXION_SPRING_CONSTANT,
                                                                             ClothSystem::FLEXION_SPRING_REST_LENGTH);
                    totalFlexionForce += flexionForce;
                }

                totalForce = gravityForce + windForce + dragForce + totalStructuralForce + totalShearForce + totalFlexionForce;
            }

            derivative.push_back(velocityOfParticle);
            derivative.push_back(totalForce/ClothSystem::PARTICLE_MASS);
        }
    }

	return derivative;
}

void ClothSystem::setInitState()
{
    this->m_vVecState.clear();
    this->windOn = false;

    float xIncrement = 0.3f;
    float zIncrement = 0.3f;
    float initZ = 0.3f;
    float initX = 0.3f;

    float x = initX;
    float y = 0.0f;
    float z = initZ;

    for (int i = 0; i < this->height; ++i)
    {
        for (int j = 0; j < this->width; ++j)
        {
            Vector3f initPosition(x, y, z);
            Vector3f initVelocity = Vector3f::ZERO;

            this->m_vVecState.push_back(initPosition);
            this->m_vVecState.push_back(initVelocity);

            x += xIncrement;
        }

        x = initX ;
        z += initZ;
    }
}

void ClothSystem::reinitializeSystem()
{
    this->setInitState();
}

bool ClothSystem::isTopCorner(int i, int j)
{
    return (i == 0) &&  (j == 0 || j == this->width - 1);

}
int ClothSystem::getIndex(int i, int j)
{
    int index;

    if (i >= 0 && i < this->height &&
        j >= 0 && j < this->width)
    {
        index = i * this->width + j;
    }

    else
    {
        index = ClothSystem::NULL_INDEX;
    }

    return index;
}

vector<int> ClothSystem::getIndexesForStructuralSprings(int i, int j)
{
    vector<int> indexes;

    int springIndex1 = this->getIndex(i + 1, j);
    int springIndex2 = this->getIndex(i, j + 1);
    int springIndex3 = this->getIndex(i - 1, j);
    int springIndex4 = this->getIndex(i, j - 1);

    if (springIndex1 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(springIndex1);
    }

    if (springIndex2 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(springIndex2);
    }

    if (springIndex3 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(springIndex3);
    }

    if (springIndex4 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(springIndex4);
    }

    return indexes;
}

vector<int> ClothSystem::getIndexesForShearSprings(int i, int j)
{
    vector<int> indexes;

    int shearIndex1 = this->getIndex(i + 1, j + 1);
    int shearIndex2 = this->getIndex(i + 1, j - 1);
    int shearIndex3 = this->getIndex(i - 1, j + 1);
    int shearIndex4 = this->getIndex(i - 1, j - 1);

    if (shearIndex1 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(shearIndex1);
    }

    if (shearIndex2 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(shearIndex2);
    }

    if (shearIndex3 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(shearIndex3);
    }

    if (shearIndex4 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(shearIndex4);
    }

    return indexes;
}

vector<int> ClothSystem::getIndexesForFlexionSprings(int i, int j)
{
    vector<int> indexes;

    int flexionIndex1 = this->getIndex(i + 2, j);
    int flexionIndex2 = this->getIndex(i, j + 2);
    int flexionIndex3 = this->getIndex(i - 2, j);
    int flexionIndex4 = this->getIndex(i, j - 2);

    if (flexionIndex1 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(flexionIndex1);
    }

    if (flexionIndex2 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(flexionIndex2);
    }

    if (flexionIndex3 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(flexionIndex3);
    }

    if (flexionIndex4 != ClothSystem::NULL_INDEX)
    {
        indexes.push_back(flexionIndex4);
    }

    return indexes;
}

void ClothSystem::draw()
{
    if (this->shadeCloth)
    {
        this->drawShadedCloth();
    }

    else
    {
        this->drawParticles();
        this->drawAllSprings();
    }
}

void ClothSystem::drawParticles()
{
    glColor3f(245/255.0, 15.0/255.0, 222.0/255.0);
    for (int i = 0; i < m_numParticles; i++)
    {
        Vector3f pos = PhysicsUtilities::getPositionOfParticle(this->m_vVecState, i);
        glDisable(GL_LIGHTING);
        glPushMatrix();
        glTranslatef(pos[0], pos[1], pos[2] );
        glutSolidSphere(0.075f,10.0f,10.0f);
        glPopMatrix();
        glEnable(GL_LIGHTING);
    }
}

void ClothSystem::drawAllSprings()
{
    for (int i = 0; i < this->height; ++i)
    {
        for (int j = 0; j < this->width; ++j)
        {
            int indexOfParticle = this->getIndex(i, j);
            Vector3f posParticle = PhysicsUtilities::getPositionOfParticle(m_vVecState, indexOfParticle);

            // Draw the structural springs
            glColor3f(15.0/255.0, 245.0/255.0, 214.0/255.0);
            vector<int> structuralSpringIndexes = this->getIndexesForStructuralSprings(i, j);
            this->drawLinesForSprings(posParticle, structuralSpringIndexes);

            // Draw the shear springs
            glColor3f(237.0/255.0, 245.0/255.0, 15.0/255.0);
            vector<int> shearSpringIndexes = this->getIndexesForShearSprings(i, j);
            this->drawLinesForSprings(posParticle, shearSpringIndexes);

            // Draw the flexion springs
            glColor3f(183.0/255.0, 18.0/255.0, 243.0/255.0);
            vector<int> flexionSpringIndexes = this->getIndexesForFlexionSprings(i, j);
            this->drawLinesForSprings(posParticle, flexionSpringIndexes);
        }
    }
}

void ClothSystem::drawLinesForSprings(Vector3f posParticle, vector<int> connectedIndexes)
{
    for (int otherParticleIndex : connectedIndexes)
    {
        glDisable(GL_LIGHTING);
        Vector3f posOtherParticle = PhysicsUtilities::getPositionOfParticle(m_vVecState, otherParticleIndex);
        glBegin(GL_LINES);
        glVertex3d(posParticle.x(), posParticle.y(), posParticle.z());
        glVertex3d(posOtherParticle.x(), posOtherParticle.y(), posOtherParticle.z());
        glEnd();
        glEnable(GL_LIGHTING);
    }
}

void ClothSystem::drawShadedCloth()
{
    this->calculateNormals();

    for (unsigned int i = 0; i < this->faces.size(); i++)
    {
        drawFace(this->faces[i], false);
        drawFace(this->faces[i], true);
    }

    this->normals.clear();
}

void ClothSystem::drawFace(vector<int> faceVector, bool reverse)
{
    int i1;
    int i2;
    int i3;

    if (reverse)
    {
        i1 = faceVector[2];
        i2 = faceVector[1];
        i3 = faceVector[0];
    }

    else
    {
        i1 = faceVector[0];
        i2 = faceVector[1];
        i3 = faceVector[2];
    }

    Vector3f pos1 = PhysicsUtilities::getPositionOfParticle(m_vVecState, i1);
    Vector3f pos2 = PhysicsUtilities::getPositionOfParticle(m_vVecState, i2);
    Vector3f pos3 = PhysicsUtilities::getPositionOfParticle(m_vVecState, i3);

    int normalMult = reverse ? -1.0 : 1.0;
    Vector3f normal1 = normalMult * normals[i1];
    Vector3f normal2 = normalMult * normals[i2];
    Vector3f normal3 = normalMult * normals[i3];

    glBegin(GL_TRIANGLES);

    // First normal and first vertex
    glNormal3d(normal1.x(), normal1.y(), normal1.z());
    glVertex3d(pos1.x(), pos1.y(), pos1.z());

    // Second normal and second vertex
    glNormal3d(normal2.x(), normal2.y(), normal2.z());
    glVertex3d(pos2.x(), pos2.y(), pos2.z());

    // Third normal and third vertex
    glNormal3d(normal3.x(), normal3.y(), normal3.z());
    glVertex3d(pos3.x(), pos3.y(), pos3.z());

    glEnd();

}

void ClothSystem::calculateNormals()
{
    for (int i = 0; i < this->height; ++i)
    {
        for (int j = 0; j < this->width; ++j)
        {
            int indexOfParticle = this->getIndex(i, j);

            int index1 = this->getIndex(i, j - 1);
            int index2 = this->getIndex(i - 1, j);
            int index3 = this->getIndex(i, j + 1);
            int index4 = this->getIndex(i + 1, j);


            Vector3f vec1;
            Vector3f vec2;
            Vector3f vec3;
            Vector3f vec4;

            int countVectors = 0;
            if (index1 != ClothSystem::NULL_INDEX)
            {
                vec1 = this->getVectorBetweenParticles(indexOfParticle, index1);
                ++countVectors;
            }

            if (index2 != ClothSystem::NULL_INDEX)
            {
                vec2 = this->getVectorBetweenParticles(indexOfParticle, index2);
                ++countVectors;
            }

            if (index3 != ClothSystem::NULL_INDEX)
            {
                vec3 = this->getVectorBetweenParticles(indexOfParticle, index3);
                ++countVectors;
            }

            if (index4 != ClothSystem::NULL_INDEX)
            {
                vec4 = this->getVectorBetweenParticles(indexOfParticle, index4);
                ++countVectors;
            }

            assert(countVectors != 0);

            Vector3f cross1 = Vector3f::ZERO;
            Vector3f cross2 = Vector3f::ZERO;
            Vector3f cross3 = Vector3f::ZERO;
            Vector3f cross4 = Vector3f::ZERO;

            if (index1 != ClothSystem::NULL_INDEX && index4 != ClothSystem::NULL_INDEX)
            {
                cross4 = Vector3f::cross(vec4, vec1);
            }

            if (index1 != ClothSystem::NULL_INDEX && index2 != ClothSystem::NULL_INDEX)
            {
                cross4 = Vector3f::cross(vec1, vec2);
            }

            if (index2 != ClothSystem::NULL_INDEX && index3 != ClothSystem::NULL_INDEX)
            {
                cross2 = Vector3f::cross(vec2, vec3);
            }

            if (index3 != ClothSystem::NULL_INDEX && index4 != ClothSystem::NULL_INDEX)
            {
                cross3 = Vector3f::cross(vec3, vec4);
            }

            Vector3f finalNormal = (cross1 + cross2 + cross3 + cross4)/(1.0 * countVectors);

            this->normals.push_back(finalNormal.normalized());
        }
    }
}

void ClothSystem::calculateFaces()
{
    for (int i = 0; i < this->height - 1; ++i)
    {
        for (int j = 0; j < this->width - 1; ++j)
        {
            vector<int> face1;
            face1.push_back(this->getIndex(i, j));
            face1.push_back(this->getIndex(i, j + 1));
            face1.push_back(this->getIndex(i + 1, j));

            vector<int>face2;
            face2.push_back(this->getIndex(i, j + 1));
            face2.push_back(this->getIndex(i + 1, j + 1));
            face2.push_back(this->getIndex(i + 1, j));

            faces.push_back(face1);
            faces.push_back(face2);
        }
    }
}

Vector3f ClothSystem::getVectorBetweenParticles(int index1, int index2) //pos(index2) - pos(index1)
{
    Vector3f pos1 = PhysicsUtilities::getPositionOfParticle(this->m_vVecState, index1);
    Vector3f pos2 = PhysicsUtilities::getPositionOfParticle(this->m_vVecState, index2);
    return pos2 - pos1;
}

