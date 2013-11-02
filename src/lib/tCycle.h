#ifndef TCYCLE_H
#define TCYCLE_H

#include <map>
#include <vcg/space/point3.h>

using namespace std;

class tCycle {

public:
	map< int, map< Point3<float>, float > > m;
	map< int, map< Point3<float>, float > > p;

};

#endif
