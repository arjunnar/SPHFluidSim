#include "debugutilities.h"

using namespace std;

void DebugUtilities::printVectorInt(vector<int> &vec)
{
    for (int i : vec)
    {
        cout << i << " ";
    }
}

void DebugUtilities::printVector3f(Vector3f &vec)
{
    cout << " ( " << vec.x() << " , " << vec.y() << " , " << vec.z() << " ) ";
    cout << endl;
}


