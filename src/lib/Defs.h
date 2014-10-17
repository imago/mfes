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

/// 1 kJ/mol in e^2/(A*mol)
#define E2A 7.197585e-4;

#define CF_LAMBDA 0.5
#define TAUBIN_LAMBDA 0.5
#define TAUBIN_MU -0.53
// standard: mu = lambda = 0.5
// meshlab: mu = -0.53, lambad = 0.5

// kappa squared times I/eps_r [1/(eps_0)] but in Angstroem
#define DL 8.43249149e-22

#define CONVERT 5.74346052632
		//5.74342695

#define VOL_MESH_TRIES 5

#define NOT_IN_ST 111
#define HEAPSIZE 4000000

const double MEADUNITS = 1.3806505E-23 * 4.336663E17; /* e^2/(A * K * mol) */
const double T = 300; /* K */

#endif
