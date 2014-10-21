#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

#include "VCG.h"

typedef boost::property_tree::ptree INI;


#include "LSMS/molTypes.h"
#include "LSMS/gridTypes.h"
#include "LSMS/grid.h"
#include "LSMS/pdbLoader.h"
#include "LSMS/sceneBuilder.h"


Molecule mol = NULL; // molecule
float minx,miny,minz,maxx,maxy,maxz;
double s;

class STLFacet
{
public:
  Point3f n;
  Point3f v[3];
//  short attr;
};

typedef typename mMesh::FaceIterator FaceIterator;
typedef typename mMesh::VertexIterator VertexIterator;

class LSMS {
public:
	int calcMC(mMesh &mSurface, vector<Atom> &atomList, INI& ini, string mode = "protein") {

		string debug = ini.get<string>("model.debug");
		double probeRadius = atof(ini.get_optional<string>("experiment.probe_radius").get_value_or("1.4").c_str());
		double ionExclR = atof(ini.get_optional<string>("experiment.ionr").get_value_or("2").c_str());
		bool calc_cavity = false;
		unsigned int gridSize = 512;
		float extendR = 0;

		if (mode == "protein") {
		    gridSize = ini.get<unsigned int>("model.grid_resolution");
		} else if (mode == "residue")  {
			gridSize = ini.get<unsigned int>("model.grid_residue_resolution");
		} else if (mode == "exclusion")  {
			gridSize = ini.get<unsigned int>("model.grid_resolution");
			extendR = ionExclR; // ion exclusion layer distance [Angstroem]
		} else {
			// cavity calculation is turned on
			calc_cavity = true;
		}

		cout << "grid size: " << gridSize << endl;
		cout << "probe radius: " << probeRadius << endl;
		cout << "atom list size: " << atomList.size() << endl;

		Molecule mol = (Molecule)malloc(sizeof(struct prot));
		mol->npoints = atomList.size();
		mol->xpoints = (float *)malloc(mol->npoints*sizeof(float));
		mol->ypoints = (float *)malloc(mol->npoints*sizeof(float));
		mol->zpoints = (float *)malloc(mol->npoints*sizeof(float));
		mol->rpoints = (float *)malloc(mol->npoints*sizeof(float));

		for (unsigned int i = 0; i<atomList.size() ; i++) {
			Atom currentAtom = atomList.at(i);
			Point3<float> currentCoord = currentAtom.getCoord();
			mol->xpoints[i] = currentCoord.X();
			mol->ypoints[i] = currentCoord.Y();
			mol->zpoints[i] = currentCoord.Z();
			mol->rpoints[i] = currentAtom.getRadius() + extendR; 
			if (debug == "yes")
				cout << "x, y, z, r: " << mol->xpoints[i] << ", " << mol->ypoints[i] << ", " << mol->zpoints[i] << ", " << mol->rpoints[i] << endl;
		}


		char inner = 0;
		float k = 0;

		if (calc_cavity)
			inner = 1; //calc_cavity-'0';
		else
			inner = 0;

		s = centerMolecule(mol,gridSize);

		probeRadius = probeRadius*s;

		//grid = createGrid(N);
		short s2;
		//		int i = N;
		int i = gridSize;
		int halfLength = 256;
		short m2,n2,k2;
		lGrid g = (lGrid)malloc(sizeof(struct g));
		//		g->N = N;
		g->N = gridSize;

		if (gridSize <= 512){
		  s2 = 512/i;
		  halfLength = 256;
		} else {
		  s2 = gridSize/i;
		  halfLength = gridSize/2;
		}
		g->stepSize = s2;

		g->matrix = (GridPoint***)malloc(i*sizeof(GridPoint**));
		for (m2=0;m2<i;m2++)
			g->matrix[m2] = (GridPoint**)malloc(i*sizeof(GridPoint*));
		for (m2=0;m2<i;m2++)
			for (n2=0;n2<i;n2++)
				g->matrix[m2][n2] = (GridPoint*)malloc(i*sizeof(GridPoint));

		for(m2=0;m2<i;m2++)
			for (n2=0;n2<i;n2++)
				for (k2=0;k2<i;k2++)
				{
					g->matrix[m2][n2][k2].point.x = -halfLength+m2*s2;
					g->matrix[m2][n2][k2].point.y = -halfLength+n2*s2;
					g->matrix[m2][n2][k2].point.z = -halfLength+k2*s2;
					g->matrix[m2][n2][k2].phi = 1; // everything is outside surface
					g->matrix[m2][n2][k2].from=-2;
					g->matrix[m2][n2][k2].dist= 0;
				}

		signDistanceGridMol(g,mol,probeRadius,gridSize);

		// Don't do shrinking, if probe radius is  0
		if (probeRadius == 0)
			printf("Not shrinking because probe radius is 0.\n");
		else
		  shrink(g,probeRadius,gridSize);

		//		k = s/(512/N);
		if (gridSize <= 512)
		  k = s/(512/gridSize);
		else 
		  k = s;

		k = k * k * k;

		vector<Triangle> result;

		float volume = fastMarching(g,inner)/k;
		printf("\nTotal cavity/molecule volume = %f A^3\n", volume);

		if (calc_cavity && volume == 0){
			cout << "no cavity found!" << endl;
			return 0;
		}
		init(result, g);

		cout << "Counting " << result.size() << " triangles." << endl;

		float minMaxX = minx + maxx;
		float minMaxY = miny + maxy;
		float minMaxZ = minz + maxz;

		for (unsigned int i = 0; i < result.size(); i++){
			Triangle currentTriangle = result.at(i);
			STLFacet f;

			f.n = currentTriangle.normal;
			currentTriangle.t1.X() = (currentTriangle.t1.X()*1/(s) + (minMaxX/2));
			currentTriangle.t1.Y() = (currentTriangle.t1.Y()*1/(s) + (minMaxY/2));
			currentTriangle.t1.Z() = (currentTriangle.t1.Z()*1/(s) + (minMaxZ/2));
			currentTriangle.t2.X() = (currentTriangle.t2.X()*1/(s) + (minMaxX/2));
			currentTriangle.t2.Y() = (currentTriangle.t2.Y()*1/(s) + (minMaxY/2));
			currentTriangle.t2.Z() = (currentTriangle.t2.Z()*1/(s) + (minMaxZ/2));
			currentTriangle.t3.X() = (currentTriangle.t3.X()*1/(s) + (minMaxX/2));
			currentTriangle.t3.Y() = (currentTriangle.t3.Y()*1/(s) + (minMaxY/2));
			currentTriangle.t3.Z() = (currentTriangle.t3.Z()*1/(s) + (minMaxZ/2));

			f.v[0] = currentTriangle.t1;
			f.v[1] = currentTriangle.t2;
			f.v[2] = currentTriangle.t3;

			FaceIterator fi=Allocator<mMesh>::AddFaces(mSurface,1);
			VertexIterator vi=Allocator<mMesh>::AddVertices(mSurface,3);
			for(int k=0;k<3;++k)
			{
				(*vi).P().Import(f.v[k]);
				(*fi).V(k)=&*vi;
				++vi;
			}
		}

		free(mol);
		//		for (int i = 0; i<N; i++) {
		for (int i = 0; i<gridSize; i++) {
		  //			for (int j = 0; j<N; j++) {
			for (int j = 0; j<gridSize; j++) {
				free(g->matrix[i][j]);
			}
		}
		//		for (int i = 0; i<N; i++) {
		for (int i = 0; i<gridSize; i++) {
			free(g->matrix[i]);
		}
		free(g);

		return 1;
	}

private:

	void init(vector<Triangle>& result, lGrid &grid){
		buildScene(result, grid);
	}

	double centerMolecule(Molecule mol, int gridSize)
	{
	  double sx,sy,sz;
	  double max=0;
	  int i;
	  double rmax = 0;

	  // set starting point
	  minx = mol->xpoints[0];
	  maxx = mol->xpoints[0];
	  miny = mol->ypoints[0];
	  maxy = mol->ypoints[0];
	  minz = mol->zpoints[0];
	  maxz = mol->zpoints[0];
	  
	    
	  for (i=0;i<mol->npoints;i++){
	    if (mol->xpoints[i]-mol->rpoints[i]<minx)	minx=mol->xpoints[i]-mol->rpoints[i];
	    if (mol->xpoints[i]+mol->rpoints[i]>maxx)	maxx=mol->xpoints[i]+mol->rpoints[i];
	    if (mol->ypoints[i]-mol->rpoints[i]<miny)	miny=mol->ypoints[i]-mol->rpoints[i];
	    if (mol->ypoints[i]+mol->rpoints[i]>maxy)	maxy=mol->ypoints[i]+mol->rpoints[i];
	    if (mol->zpoints[i]-mol->rpoints[i]<minz)	minz=mol->zpoints[i]-mol->rpoints[i];
	    if (mol->zpoints[i]+mol->rpoints[i]>maxz)	maxz=mol->zpoints[i]+mol->rpoints[i];
	    if (rmax < mol->rpoints[i])  rmax=ceil(mol->rpoints[i]);
	  }
	    
	  sx = gridSize/(maxx-minx+2*rmax);
	  sy = gridSize/(maxy-miny+2*rmax);
	  sz = gridSize/(maxz-minz+2*rmax);
	  
	  s = 0.0f;
	  if (sx<sy && sx<sz){ max=(maxx-minx); s=sx; }
	  else if (sy<sx && sy<sz){ max=(maxy-miny); s=sy;}
	  else { max=(maxz-minz); s=sz; }
	  
	  
	  for (i=0;i<mol->npoints;i++){
	    //		printf("x y z r before %f %f %f %f\n", mol->xpoints[i], mol->ypoints[i], mol->zpoints[i], mol->rpoints[i]);
	    mol->xpoints[i]=(s*(mol->xpoints[i]-((minx+maxx)/2)));
	    mol->ypoints[i]=(s*(mol->ypoints[i]-((miny+maxy)/2)));
	    mol->zpoints[i]=(s*(mol->zpoints[i]-((minz+maxz)/2)));
	    //		printf("x y z r after %f %f %f %f\n", mol->xpoints[i], mol->ypoints[i], mol->zpoints[i], mol->rpoints[i]);
	    mol->rpoints[i]=mol->rpoints[i]*s;
	  }
	  
	  cout << "min: " << minx << ", " << miny << ", " << minz << "; max: " << maxx << ", " << maxy << ", " << maxz << endl;
	  cout << "abs (x, y, z): (" << (maxx-minx) << ", " << (maxy-miny) << ", " << (maxz-minz) << ")" << endl;
	  cout << "scaling: " << s << endl;
	  cout << "grid size: " << gridSize << endl;
	  cout << "r_max: " << rmax << endl;

	  /*	  double h_eff = 512/gridSize* max/gridSize;
	  if (gridSize > 512)
	    h_eff = max/gridSize;
	  cout << "h_eff: " << h_eff << endl;
	  */

	  return s;
	}

};
