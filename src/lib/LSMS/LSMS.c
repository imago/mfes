#include "molTypes.h"
#include "gridTypes.h"
#include "grid.h"
#include "pdbLoader.h"
#include "sceneBuilder.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Timing variables
clock_t m_tTime;
int m_iNbFrames;
clock_t ttime;


short WIDTH = 512;
short HEIGHT = 512;
short N = 128;
float PR = 1.4;
Molecule mol = NULL;
Molecule catrace = NULL;
Grid grid = NULL;

void init(void)
{
	clock_t t;

	t = clock();
	buildScene(grid);
	printf("Scene building time %f seconds\n",(clock()-t)/(float)CLOCKS_PER_SEC);
	t = clock();
	printf("Total Molecular Surface Generation Time : %f seconds\n",(clock()-ttime)/(float)CLOCKS_PER_SEC);

}

float centerMolecule(Molecule mol, Molecule caTrace)
{
	// scale and translate the molecule so that the atom coordinates are between -250 and 250
	float minx,miny,minz,maxx,maxy,maxz;
	float s,sx,sy,sz;
	int i;

	minx = mol->xpoints[0];
	maxx = mol->xpoints[0];
	miny = mol->ypoints[0];
	maxy = mol->ypoints[0];
	minz = mol->zpoints[0];
	maxz = mol->zpoints[0];

	for (i=1;i<mol->npoints;i++){
		if (mol->xpoints[i]<minx)	minx=mol->xpoints[i];
		if (mol->xpoints[i]>maxx)	maxx=mol->xpoints[i];
		if (mol->ypoints[i]<miny)	miny=mol->ypoints[i];
		if (mol->ypoints[i]>maxy)	maxy=mol->ypoints[i];
		if (mol->zpoints[i]<minz)	minz=mol->zpoints[i];
		if (mol->zpoints[i]>maxz)	maxz=mol->zpoints[i];
	}
	sx = 400.0f/(maxx-minx);
	sy = 400.0f/(maxy-miny);
	sz = 400.0f/(maxz-minz);
	s = 0.0f;
	if (sx<sy && sx<sz)	s=sx;
	else if (sy<sx && sy<sz)s=sy;
	else s=sz;
	//printf("%f %f %f\n",sx,sy,sz);
	//printf("s = %f\n",s);
	printf("SPDBV quality = %f\n",(PR*s)/(512/N));

	for (i=0;i<mol->npoints;i++){
		mol->xpoints[i]=(s*(mol->xpoints[i]-((minx+maxx)/2)));
		mol->ypoints[i]=(s*(mol->ypoints[i]-((miny+maxy)/2)));
		mol->zpoints[i]=(s*(mol->zpoints[i]-((minz+maxz)/2)));
		mol->rpoints[i]=mol->rpoints[i]*s;
	}
	for (i=0;i<catrace->npoints;i++){
		catrace->xpoints[i]=(s*(caTrace->xpoints[i]-((minx+maxx)/2)));
		catrace->ypoints[i]=(s*(caTrace->ypoints[i]-((miny+maxy)/2)));
		catrace->zpoints[i]=(s*(caTrace->zpoints[i]-((minz+maxz)/2)));
		catrace->rpoints[i]=caTrace->rpoints[i]*s;
	}
	return s;
}

int main(int argc, char** argv)
{
	float s;
	clock_t t;
	char inner = 0;
	float k = 0;

	if (argc!=5)
	{
		printf("Usage: levelSet <pdb File> <gridSize> <probeSize> <inner caves(1) or outer surface(0)\n");
		return 1;
	}
	t = clock();
	mol = loadMolecule(argv[1]);
	catrace = loadCATrace(argv[1]);
	printf("Number of CA atoms = %d\n",catrace->npoints);
	printf("Molecule loading time %f seconds\n",(clock()-t)/(float)CLOCKS_PER_SEC);
	N = (short)atoi(argv[2]);
	PR = (float)atof(argv[3]);
	inner = argv[4][0]-'0';
	printf("Number of atoms: %d\n",mol->npoints);
	t = clock();
	ttime = clock();
	s = centerMolecule(mol, catrace);
	printf("Molecule centering time %f seconds\n",(clock()-t)/(float)CLOCKS_PER_SEC);
	PR = PR*s;

	t = clock();
	grid = createGrid(N);
	printf("Grid creation time %f seconds\n",(clock()-t)/(float)CLOCKS_PER_SEC);
	t = clock();

	signDistanceGridMol(grid,mol,PR);
	printf("Grid initialization time %f seconds\n",(clock()-t)/(float)CLOCKS_PER_SEC);
	
	t = clock();

	shrink(grid,PR);
/*	c = findProbesMol(grid,PR);
	printf("There are %d probes.\n",c);
	printf("ProbesMol time %f seconds.\n",(clock()-t)/(float)CLOCKS_PER_SEC);
*/
	printf("Shrinking time %f seconds.\n",(clock()-t)/(float)CLOCKS_PER_SEC);
	
	
	t = clock();
	k = s/(512/N);
	k = k * k * k;
	printf("Total cavity/molecule volume = %f\n",fastMarching(grid,inner)/k);
	printf("Fast marching time %f seconds\n",(clock()-t)/(float)CLOCKS_PER_SEC);
	
	init();

	return 0;
}
