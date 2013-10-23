#include "VCG.h"

#ifndef DEF
#define DEF

#define CNAME mfes

struct Triangle{
	Point3<float> normal;
	Point3<float> t1;
	Point3<float> t2;
	Point3<float> t3;
};

#define CF_LAMBDA 0.5
#define TAUBIN_LAMBDA 0.5
#define TAUBIN_MU -0.53
// standard: mu = lambda = 0.5
// meshlab: mu = -0.53, lambad = 0.5

#endif
