#ifndef DEBUGUTILITIES_H
#define DEBUGUTILITIES_H

#include <vector>
#include "vecmath.h"
#include <iostream>

using namespace std;

class DebugUtilities
{
public:
    static void printVectorInt(vector<int> &vec);
    static void printVector3f(Vector3f &vec);
};

#endif // DEBUGUTILITIES_H
