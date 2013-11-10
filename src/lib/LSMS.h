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

// Timing variables
clock_t m_tTime;
int m_iNbFrames;
clock_t ttime;


short WIDTH = 512;
short HEIGHT = 512;
short N = 32;
double PR = 1.4;
Molecule mol = NULL;
float SIZEMAX = 400.0f;
float minx,miny,minz,maxx,maxy,maxz;
double s;
bool usePredef = false;
float _extendedRadius = 0;

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

		bool calc_cavity = false;
		unsigned int gridSize;
		if (mode == "protein") {
		    gridSize = ini.get<unsigned int>("model.grid_resolution");
		} else if (mode == "residue")  {
			gridSize = ini.get<unsigned int>("model.grid_residue_resolution");
		} else {
			// cavity calculation is turned on
			gridSize = 512;
			calc_cavity = true;
		}

		double probeRadius = ini.get<double>("experiment.probe_radius");

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
			mol->rpoints[i] = currentAtom.getRadius(); //+extendedRadius;
			if (debug == "yes")
				cout << "x, y, z, r: " << mol->xpoints[i] << ", " << mol->ypoints[i] << ", " << mol->zpoints[i] << ", " << mol->rpoints[i] << endl;
		}


		char inner = 0;
		float k = 0;

		N = gridSize;
		PR = probeRadius;
		if (calc_cavity)
			inner = 1; //calc_cavity-'0';
		else
			inner = 0;

		s = centerMolecule(mol);

		PR = PR*s;

		//grid = createGrid(N);
		short s2;
		int i = N;
		short m2,n2,k2;
		lGrid g = (lGrid)malloc(sizeof(struct g));
		g->N = N;
		s2 = 512/i;
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
					g->matrix[m2][n2][k2].point.x = -256+m2*s2;
					g->matrix[m2][n2][k2].point.y = -256+n2*s2;
					g->matrix[m2][n2][k2].point.z = -256+k2*s2;
					g->matrix[m2][n2][k2].phi = 1; // everything is outside surface
					g->matrix[m2][n2][k2].from=-2;
					g->matrix[m2][n2][k2].dist= 0;
				}

		signDistanceGridMol(g,mol,PR);

		// Don't do shrinking, if PR == 0
		if (PR == 0)
			printf("Not shrinking because PR == 0.\n");
		else
			shrink(g,PR);

		k = s/(512/N);
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
		for (int i = 0; i<N; i++) {
			for (int j = 0; j<N; j++) {
				free(g->matrix[i][j]);
			}
		}
		for (int i = 0; i<N; i++) {
			free(g->matrix[i]);
		}
		free(g);

		return 1;
	}

private:

	void init(vector<Triangle>& result, lGrid &grid){
		buildScene(result, grid);
	}

	double centerMolecule(Molecule mol)
	{
		// scale and translate the molecule so that the atom coordinates are between -250 and 250
		double sx,sy,sz;
		double max=0;
		int i;

		if (!usePredef){
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
			}

			if (mol->npoints==1)
			{
	/*				sx = 40.0f;
					sy = 40.0f;
					sz = 40.0f;
					*/

				sx = SIZEMAX/(maxx-minx+2);
				sy = SIZEMAX/(maxy-miny+2);
				sz = SIZEMAX/(maxz-minz+2);

			}
			else
			{
				sx = SIZEMAX/(maxx-minx+2);
				sy = SIZEMAX/(maxy-miny+2);
				sz = SIZEMAX/(maxz-minz+2);

			}
			s = 0.0f;
			if (sx<sy && sx<sz){ max=(maxx-minx); s=sx; }
			else if (sy<sx && sy<sz){ max=(maxy-miny); s=sy;}
			else { max=(maxz-minz); s=sz; }
		}

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
		cout << "sizemax: " << SIZEMAX << endl;
		double h_eff = 512/N* max/SIZEMAX*0.5;
		cout << "h_eff: " << h_eff << endl;
		return s;
	}

/*	void writeSTL(string fileName){
		float maxx = 0, maxy = 0, maxz = 0, minx = 0, miny = 0, minz = 0;
		ofstream stlFile;
		stlFile.open (fileName.c_str());
		stlFile << "#Volume: " << getVolume() << endl;
		stlFile << "#SASA: " << getSASA() << endl;
		stlFile << "solid object" << endl;
		for (int i = 0; i < triangles.size(); i++) {
			Triangle currentTriangle = triangles.at(i);
			Point3D currentNormal;
			Point3D p1 = currentTriangle.p1;
			Point3D p2 = currentTriangle.p2;
			Point3D p3 = currentTriangle.p3;
			currentNormal = currentTriangle.calcNormal();
			currentNormal.normalize();
			stlFile << "	facet normal " << currentNormal.x << " " << currentNormal.y << " " << currentNormal.z << endl;
			stlFile << "		outer loop" << endl;
			stlFile << "			vertex	" << p1.x << "	" << p1.y << "	" << p1.z << endl;
			stlFile << "			vertex	" << p2.x << "	" << p2.y << "	" << p2.z << endl;
			stlFile << "			vertex	" << p3.x << "	" << p3.y << "	" << p3.z << endl;
			stlFile << "		endloop" << endl;
			stlFile << "	endfacet" << endl;
			if (p1.x < minx )
				minx = p1.x;
			if (p2.x < minx )
				minx = p2.x;
			if (p3.x < minx )
				minx = p3.x;
			if (p1.y < miny )
				miny = p1.y;
			if (p2.y < miny )
				miny = p2.y;
			if (p3.y < miny )
				miny = p3.y;
			if (p1.z < minz )
				minz = p1.z;
			if (p2.z < minz )
				minz = p2.z;
			if (p3.z < minz )
				minz = p3.z;
			if (p1.x > maxx )
				maxx = p1.x;
			if (p2.x > maxx )
				maxx = p2.x;
			if (p3.x > maxx )
				maxx = p3.x;
			if (p1.y > maxy )
				maxy = p1.y;
			if (p2.y > maxy )
				maxy = p2.y;
			if (p3.y > maxy )
				maxy = p3.y;
			if (p1.z > maxz )
				maxz = p1.z;
			if (p2.z > maxz )
				maxz = p2.z;
			if (p3.z > maxz )
				maxz = p3.z;
		}
		stlFile << "endsolid object" << endl;
		stlFile.close();
		cout << "STL-File: (minx, miny, minz) = (" << minx << ", " << miny << ", " << minz << ");" << endl;
		cout << "(maxx, maxy, maxz) = (" << maxx << ", " << maxy << ", " << maxz << ")" << endl;

	}
*/

/*	void setTriangles(vector<Triangle> &newtriangles){
		triangles = newtriangles;
	}

	float getArea(Triangle t){
		float x1 = t.p2.x - t.p1.x;
		float x2 = t.p2.y - t.p1.y;
		float x3 = t.p2.z - t.p1.z;

		float y1 = t.p3.x - t.p1.x;
		float y2 = t.p3.y - t.p1.y;
		float y3 = t.p3.z - t.p1.z;

		float area = 0.5*sqrt(pow(x2*y3-x3*y2, 2)+pow(x3*y1-x1*y3, 2)+pow(x1*y2-x2*y1,2));

		return area;
	}

	void readSurface(float scaling){
		string currentLine;
		float x1, x2, x3, y1, y2, y3, z1, z2, z3;


		vector<string> lines;
		vector<Triangle> result;

		ifstream ifile ("triangles.dat");
		if (!ifile)
		{
			cout << "Error reading triangles.dat." << endl;
		}
		else
		{
			while (getline(ifile, currentLine)){
				lines.push_back(currentLine);
			}
		}

		Point3D center;

		float area = 0;

		for (int i=0; i<lines.size(); i++) {
			istringstream iss(lines.at(i));
			iss >> x1 >> x2 >> x3 >> y1 >> y2 >> y3 >> z1 >> z2 >> z3;
			Triangle newTriangle;
			float minMaxX = min.x + max.x;
			float minMaxY = min.y + max.y;
			float minMaxZ = min.z + max.z;
			newTriangle.p1.x = (x1*1/(scalingFactor) + (minMaxX/2));
			newTriangle.p1.y = (x2*1/(scalingFactor) + (minMaxY/2));
			newTriangle.p1.z = (x3*1/(scalingFactor) + (minMaxZ/2));
			newTriangle.p2.x = (y1*1/(scalingFactor) + (minMaxX/2));
			newTriangle.p2.y = (y2*1/(scalingFactor) + (minMaxY/2));
			newTriangle.p2.z = (y3*1/(scalingFactor) + (minMaxZ/2));
			newTriangle.p3.x = (z1*1/(scalingFactor) + (minMaxX/2));
			newTriangle.p3.y = (z2*1/(scalingFactor) + (minMaxY/2));
			newTriangle.p3.z = (z3*1/(scalingFactor) + (minMaxZ/2));
			center.x += (newTriangle.p1.x+newTriangle.p2.x+newTriangle.p3.x)/3;
			center.y += (newTriangle.p1.y+newTriangle.p2.y+newTriangle.p3.y)/3;
			center.z += (newTriangle.p1.z+newTriangle.p2.z+newTriangle.p3.z)/3;
			//newTriangle.print();
			result.push_back(newTriangle);
			area += getArea(newTriangle);
		}

		center = center / result.size();

		if (scaling != 1){
			for (int i = 0; i < result.size(); i++){
				Triangle currentTriangle = result.at(i);
				Point3D x, y, z;
				x = currentTriangle.p1;
				y = currentTriangle.p2;
				z = currentTriangle.p3;
				x = x - center;
				y = y - center;
				z = z - center;
				x = x * scaling;
				y = y * scaling;
				z = z * scaling;
				x = x + center;
				y = y + center;
				z = z + center;
				currentTriangle.p1 = x;
				currentTriangle.p2 = y;
				currentTriangle.p3 = z;
				result.at(i) = currentTriangle;
			}
		}

		setTriangles(result);
		setSASA(area);
		cout << "SASA: " << area << " A^2" << endl;
	}

	void readSTL(string stlFileName){
		string currentLine;
		string temp;
		float x1, x2, x3, y1, y2, y3, z1, z2, z3;

		vector<string> lines;
		vector<Triangle> result;

		ifstream ifile (stlFileName.c_str());
		if (!ifile)
		{
			cout << "Error reading: " << stlFileName.c_str() << endl;
		}
		else
		{
			while (getline(ifile, currentLine)){
				string test = trim(currentLine);
				if(test.substr(0, 8) == "#Volume:"){
					float volume = atof(test.substr(8, 20).c_str());
					setVolume(volume);
					cout << "Surface volume: " << volume << " A^3" << endl;
				}
				if(test.substr(0, 6) == "#SASA:"){
					float sasa = atof(test.substr(6, 20).c_str());
					setSASA(sasa);
					cout << "SASA: " << sasa << " A^2" << endl;
				}
				if(test.substr(0, 6) == "vertex"){
					lines.push_back(currentLine);
				}
			}
		}

		for (int i=0; i<lines.size(); i=i+3) {
			istringstream is1(lines.at(i));
			is1 >> temp >> x1 >> x2 >> x3;
			istringstream is2(lines.at(i+1));
			is2 >> temp >> y1 >> y2 >> y3;
			istringstream is3(lines.at(i+2));
			is3 >> temp >> z1 >> z2 >> z3;
			Triangle newTriangle;
			newTriangle.p1.x = x1;
			newTriangle.p1.y = x2;
			newTriangle.p1.z = x3;
			newTriangle.p2.x = y1;
			newTriangle.p2.y = y2;
			newTriangle.p2.z = y3;
			newTriangle.p3.x = z1;
			newTriangle.p3.y = z2;
			newTriangle.p3.z = z3;
			result.push_back(newTriangle);
		}

		setTriangles(result);
	}


	void readSurface(string fileName){
		string currentLine;
		float x1, x2, x3, y1, y2, y3, z1, z2, z3;

		vector<string> lines;
		vector<Triangle> result;

		ifstream ifile (fileName.c_str());
		if (!ifile)
		{
			cout << "Error reading: " << fileName << endl;
		}
		else
		{
			while (getline(ifile, currentLine)){
				lines.push_back(currentLine);
			}
		}

		for (int i=0; i<lines.size(); i++) {
			istringstream iss(lines.at(i));
			iss >> x1 >> x2 >> x3 >> y1 >> y2 >> y3 >> z1 >> z2 >> z3;
			Triangle newTriangle;
			newTriangle.p1.x = x1;
			newTriangle.p1.y = x2;
			newTriangle.p1.z = x3;
			newTriangle.p2.x = y1;
			newTriangle.p2.y = y2;
			newTriangle.p2.z = y3;
			newTriangle.p3.x = z1;
			newTriangle.p3.y = z2;
			newTriangle.p3.z = z3;
			result.push_back(newTriangle);
		}

		setTriangles(result);
	}

	vector<Triangle> getTriangles(){
		return triangles;
	}


	void setVolume(float newVolume){
		volume = newVolume;
	}

	float getVolume(){
		return volume;
	}

	void setSASA(float newSASA){
		sasa = newSASA;
	}

	float getSASA(){
		return sasa;
	}

	void setFileName(string newFileName){
		fileName = newFileName;
	}

	string getFileName(){
		return fileName;
	}
	*/
};
