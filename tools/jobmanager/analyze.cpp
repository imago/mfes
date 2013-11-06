#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

class Parameter {
public:
	string jobname;
	string plocalh;
	float pmaxh;
	float pminh;
	unsigned int lsmsN;
	float pgrading;
	int poptsteps_3d;
	int poptsteps_2d;
	float pfineness;

	int taubin;
	int hc;
	int lap;

	float blocalh;
	float bmaxh;
	float bgrading;
	int boptsteps_3d;
	int boptsteps_2d;
	float bfineness;

	string solver;
	int maxsteps;
	int order;

	string pqr;
	string boundary;

  string getIdentifier(){
	  stringstream ss;
	  ss << jobname << "_plocalh-" << plocalh <<"_pmaxh-" << pmaxh << "_pminh-"<< pminh << "_N-" << lsmsN << "_pgrading-" << pgrading << "_poptsteps3d-" << poptsteps_3d;
	  ss << "_poptsteps2d-" << poptsteps_2d << "_pfineness-"<< pfineness << "_t-" << taubin << "_hc-" << hc << "_lap-";
	  ss << lap << "_blocalh-" << blocalh << "_bmaxh-"<< bmaxh << "_bgrading-" << bgrading << "_boptsteps3d-" << boptsteps_3d << "_boptsteps2d-" << boptsteps_2d;
	  ss << "_bfineness-" << bfineness << "_solver-" << solver << "_maxsteps-" << maxsteps << "_order-" << order;
	  return ss.str();
  }
};


class Result {
public:

	Result(){
		timeTotal = -1;
		timeInverse = -1;
		timeSolve = -1;
		timeProteinStl = -1;
		timeProteinNg = -1;
		timeResidueStl = -1;
		timeResidueNg = -1;
	}

	float deltaG;
	int gpTotal;
	int gpOnSurface;
	int gpOnBoundary;
	int gpProtein;
	int gpBetween;

	float areaProtein;
	float areaBoundary;

	float volumeTotal;
	float volumeProtein;
	float volumeBetween;

	float densityTotal;
	float densityProtein;
	float densityBoundary;
	float densityInProtein;
	float densityBetween;

	int ndof;

	float timeTotal;
	float timeInverse;
	float timeSolve;
	float timeResidueStl;
	float timeResidueNg;
	float timeProteinStl;
	float timeProteinNg;

	string getCSVLine(Parameter &p){
		stringstream result;
		result << p.jobname << "\t" << p.lsmsN << "\t" << p.plocalh << "\t" <<  p.pmaxh << "\t" << p.pminh << "\t" << p.pgrading
				<<"\t" <<  p.pfineness << "\t" << p.poptsteps_2d << "\t" << p.poptsteps_3d
				<<"\t" <<  p.taubin << "\t" << p.hc << "\t" << p.lap
				<<"\t" <<  p.blocalh << "\t" << p.bmaxh << "\t" << p.bgrading << "\t" << p.bfineness
				<<"\t" <<  p.boptsteps_2d << "\t" << p.boptsteps_3d
				<<"\t" <<  p.solver << "\t" << p.maxsteps << "\t" << p.order
				<<"\t" <<  deltaG << "\t" << gpTotal << "\t" << gpOnSurface << "\t" << gpOnBoundary << "\t" << gpProtein << "\t" << gpBetween
				<<"\t" <<  areaProtein << "\t" << areaBoundary
				<<"\t" <<  volumeTotal << "\t" << volumeProtein << "\t" << volumeBetween
				<<"\t" <<  densityTotal << "\t" << densityProtein << "\t" << densityBoundary << "\t" << densityInProtein << "\t" << densityBetween
				<<"\t" <<  ndof << "\t" << timeTotal << "\t" << timeInverse << "\t" << timeSolve << "\t" << timeProteinStl << "\t" << timeProteinNg
				<<"\t" <<  timeResidueStl << "\t" << timeResidueNg;

		// "\t" <<

		return result.str();
	}

	void print(){
		cout << "deltaG: " << deltaG << endl;
		cout << "gpTotal: "<< gpTotal << endl;
		cout << "gpOnSurface: "<< gpOnSurface << endl;
		cout << "gpOnBoundary: "<< gpOnBoundary << endl;
		cout << "gpProtein: "<< gpProtein << endl;
		cout << "gpBetween: " << gpBetween << endl;
		cout << endl;
		cout << "areaProtein: "<< areaProtein << endl;
		cout << "areaBoundary: "<< areaBoundary << endl;
		cout << endl;
		cout << "volumeTotal: "<< volumeTotal << endl;
		cout << "volumeProtein: "<< volumeTotal << endl;
		cout << "volumeBetween: "<< volumeBetween << endl;
		cout << endl;
		cout << "densityTotal: "<< densityTotal << endl;
		cout << "densityProtein: "<< densityProtein << endl;
		cout << "densityBoundary: "<< densityBoundary << endl;
		cout << "densityInProtein: "<< densityInProtein << endl;
		cout << "densityBetween: "<< densityBetween << endl;
		cout << endl;
		cout << "ndof: "<< ndof << endl;
		cout << endl;
		cout << "timeTotal: "<< timeTotal << endl;
		cout << "timeInverse: "<< timeInverse << endl;
		cout << "timeSolve: "<< timeSolve << endl;
		cout << "timeProteinStl: "<< timeProteinStl << endl;
		cout << "timeProteinNg: "<< timeProteinNg << endl;
		cout << "timeResidueStl: "<< timeResidueStl << endl;
		cout << "timeResidueNg: "<< timeResidueNg << endl;
	}

};


class data // 3 vertices of each triangle
{
public:
  float x1,y1,z1;
  float x2,y2,z2;
  float x3,y3,z3;
  float nx,ny,nz;
};
class Point
{
private:
    double x, y, z;

public:
    Point(double dX=0.0, double dY=0.0, double dZ=0.0)
    {
    x = dX;
    y = dY;
    z = dZ;
    }

    friend bool operator== (const Point &cP1, const Point &cP2);
    friend bool operator!= (const Point &cP1, const Point &cP2);
    friend bool operator< (const Point &cP1, const Point &cP2);
    friend bool operator> (const Point &cP1, const Point &cP2);

    double GetX() { return x; }
    double GetY() { return y; }
    double GetZ() { return z; }
    void print(){ std::cout << "(" << x << ", " << y << ", " << z << ")" << endl; }
};

inline bool operator<(const Point &cP1, const Point &cP2)
{

  if (cP1.x == cP2.x)
    if (cP1.y == cP2.y)
      return cP1.z < cP2.z;

  if (cP1.x == cP2.x)
    return cP1.y < cP2.y;

  return cP1.x < cP2.x;
}

inline bool operator>(const Point &cP1, const Point &cP2)
{

  if (cP1.x == cP2.x)
    if (cP1.y == cP2.y)
      return cP1.z > cP2.z;

  if (cP1.x == cP2.x)
    return cP1.y > cP2.y;

  return cP1.x > cP2.x;

}

inline bool operator==(const Point &cP1, const Point &cP2)
{
    return (cP1.x == cP2.x &&
            cP1.y == cP2.y &&
            cP1.z == cP2.z);
}

inline bool operator!=(const Point &cP1, const Point &cP2)
{
    return !(cP1 == cP2);
}


float getArea(data t){
  float x1 = t.x2 - t.x1;
  float x2 = t.y2 - t.y1;
  float x3 = t.z2 - t.z1;

  float y1 = t.x3 - t.x1;
  float y2 = t.y3 - t.y1;
  float y3 = t.z3 - t.z1;

  float area = 0.5*sqrt(pow(x2*y3-x3*y2, 2)+pow(x3*y1-x1*y3, 2)+pow(x1*y2-x2*y1,2));

  return area;
}

class element {
public:
	string surfnr;
	int bcnr;
	int domin;
	int domout;
	int np;
	int p1;
	int p2;
	int p3;
	void print(){
		cout << "surfnr: " << surfnr << ", bcnr: " << bcnr << ", domin: " << domin << ", domout: " << domout << ", np: " << np << ", p1: " << p1 << ", p2: " << p2 << ", p3: " << p3 << endl;
	};
};

class VOL {
public:
	void readNGVol(string fileName){
		cout << "reading ng vol ..." << endl;
		ifstream in(fileName.c_str());

		if (!in) {
			in.close();
			cout << "Cannot open ng vol:" << fileName << endl;
			exit(0);
		}
		string pqrline, temp;

		int length = 0;
		stringstream ss;

		while( !in.eof() ) {
			getline(in, pqrline);
			stringstream ss(pqrline.c_str());
			ss >> temp;
			if (temp == "surfaceelements"){
				getline(in, pqrline);
				length = atoi(pqrline.c_str());
				for (unsigned int i = 0; i < length; i++){
					getline(in, pqrline);
					stringstream ss(pqrline.c_str());
					element newElement;
					ss >> newElement.surfnr >> newElement.bcnr >> newElement.domin >> newElement.domout >> newElement.np >> newElement.p1 >> newElement.p2 >> newElement.p3;
//					if (newElement.surfnr == "#")
//						break;
					elementList.push_back(newElement);
//					newElement.print();
				}
			}
			if (temp == "points"){
				getline(in, pqrline);
				length = atoi(pqrline.c_str());
				for (unsigned int i = 0; i < length; i++){
					getline(in, pqrline);
					stringstream ss(pqrline.c_str());
					double x, y, z;
					ss >> x >> y >> z;
//					if (newPoint.surfnr == "#")
//						break;
					Point newPoint(x, y, z);
					pointList.push_back(newPoint);
//					newPoint.print();
				}
			}
		}
	}

	double getTotalArea(){
		double totalArea = 0;

		for (unsigned int i = 0; i < elementList.size(); i++){
			element currentElement = elementList.at(i);
			Point p1 = pointList.at(currentElement.p1-1);
			Point p2 = pointList.at(currentElement.p2-1);
			Point p3 = pointList.at(currentElement.p3-1);

			float x1 = p2.GetX() - p1.GetX(); //t.x2 - t.x1;
			float x2 = p2.GetY() - p1.GetY(); //t.y2 - t.y1;
			float x3 = p2.GetZ() - p1.GetZ(); //t.z2 - t.z1;

			float y1 = p3.GetX() - p1.GetX(); //.x3 - t.x1;
			float y2 = p3.GetY() - p1.GetY(); //t.y3 - t.y1;
			float y3 = p3.GetZ() - p1.GetZ(); //t.z3 - t.z1;

			totalArea += 0.5*sqrt(pow(x2*y3-x3*y2, 2)+pow(x3*y1-x1*y3, 2)+pow(x1*y2-x2*y1,2));

		}
		return totalArea;
	}

	int getNrPoints(){
		return pointList.size();
	}

	double getTotalVolume(){
		double totalVolume = 0;
		for (unsigned int i = 0; i < elementList.size(); i++){
			element currentElement = elementList.at(i);
			Point p1 = pointList.at(currentElement.p1-1);
			Point p2 = pointList.at(currentElement.p2-1);
			Point p3 = pointList.at(currentElement.p3-1);

			totalVolume += (p1.GetX()*p2.GetY()*p3.GetZ() -
					p1.GetX()*p3.GetY()*p2.GetZ() -
					p2.GetX()*p1.GetY()*p3.GetZ() +
					p2.GetX()*p3.GetY()*p1.GetZ() +
					p3.GetX()*p1.GetY()*p2.GetZ() -
					p3.GetX()*p2.GetY()*p1.GetZ()) / 6;

		}

		return totalVolume;
	}


	vector<element> elementList;
	vector<Point> pointList;
};


string getName(string path){
	std::string filename;

	size_t pos = path.find_last_of("//");
	if(pos != std::string::npos)
	    filename.assign(path.begin() + pos + 1, path.end());
	else
	    filename = path;

	return filename;
}

void writeMeshFiles(Parameter& p){
	// molecule
	ofstream c;
	c.open ("mesh_molecule.opt");
	c << "options.localh  " << p.plocalh << endl;
	c << "options.meshsize  " << p.pmaxh << endl;
	c << "options.minmeshsize  " << p.pminh << endl;
	c << "meshoptions.fineness  " << p.pfineness << endl;
	c << "options.grading  " << p.pgrading << endl;
	c << "options.optsteps2d  " << p.poptsteps_2d << endl;
	c << "options.optsteps3d  " << p.poptsteps_3d << endl;
	c.close();

	// boundary
	c.open ("mesh_boundary.opt");
	c << "options.localh  " << p.blocalh << endl;
	c << "options.meshsize  " << p.bmaxh << endl;
	c << "options.minmeshsize  0"  << endl;
	c << "meshoptions.fineness  " << p.bfineness << endl;
	c << "options.grading  " << p.bgrading << endl;
	c << "options.optsteps2d  " << p.boptsteps_2d << endl;
	c << "options.optsteps3d  " << p.boptsteps_3d << endl;
	c.close();

}


void writeConfigMFES(Parameter& p){
	ofstream c;
	c.open ("config.ini");
	c << "[general]" << endl;
	c << "jobname = " << p.jobname << endl;
	c << "molecule = " << getName(p.pqr) << endl;
	c << "mode = energy" << endl;
	c << endl;
	c << "[experiment]" << endl;
	c << "eps_in = 4" << endl;
	c << "eps_out = 80" << endl;
	c << "probe_radius = 1.4" << endl;
	c << endl;
	c << "[model]" << endl;
	c << "grid_resolution = " << p.lsmsN << endl;
	c << "# t: taubin, lap: laplace, hc: hc laplace, aw: laplace angle weighted" << endl;
	if (p.taubin > 0)
		c << "smoothing = t " << p.taubin << endl;
	else if (p.hc > 0)
		c << "smoothing = hc " << p.hc << endl;
	else if (p.lap > 0)
		c << "smoothing = l " << p.lap << endl;
	else
		c << "smoothing =  " << endl;
	c << "boundary = boundary.vol" << endl;
	c << "options_molecule = mesh_molecule.opt" << endl;
	c << "options_boundary = mesh_boundary.opt" << endl;
	c << "# refine file: x y z h [..]" << endl;
	c << "refine_file = " << endl;
	c << "#refine.ref" << endl;
	c << "surface_stl =" << endl;
	c << "# born_vcg.stl" << endl;
	c << "volume_vol  = protein.vol" << endl;
	c << "debug = analyze" << endl;
	c << endl;
	c << "[solver]" << endl;
	c << "solver = " << p.solver << endl;
	c << "solution_order = " << p.order << endl;
	c << "maxsteps = " << p.maxsteps << endl;

	c.close();
}

void parse(string fileName, Result &result){
	cout << "reading result files ..." << endl;
	ifstream in(fileName.c_str());

	if (!in) {
		in.close();
		cout << "Cannot open results:" << fileName << endl;
		exit(0);
	}
	string pqrline, temp;

	int length = 0;
	stringstream ss;

	while( !in.eof() ) {
		getline(in, pqrline);
		stringstream ss(pqrline.c_str());
		while(ss >> temp){
			if (temp == "ndof"){
				ss >> temp >> result.ndof;
				cout << result.ndof;
			}
			if (temp == "The"){
				ss >> temp >> temp >> temp >> result.deltaG;
			}
			if (temp == "execution"){
				ss >> temp >> result.timeTotal;
			}
		}

	}
	in.close();

	in.open("ng.prof");
	if (!in) {
		in.close();
		cout << "Cannot open file: ng.prof" << endl;
		exit(0);
	}

	string times, temp2;
	double time = 0;
	while( !in.eof() ) {
		getline(in, pqrline);
		stringstream ss(pqrline.c_str());
		ss >> temp >> temp >> temp >> times >> temp >> time >> temp >> temp;
		if (temp == "Mumps"){ // last mumps in ng.prof is total
			result.timeInverse = time/atoi(times.substr(0, 1).c_str());
		}
		if (temp == "Equation"){
			result.timeSolve = time/atoi(times.substr(0, 1).c_str());
		}

	}
	in.close();

	in.open("times");
	if (!in) {
		in.close();
		cout << "Cannot open file: times" << endl;
		exit(0);
	}

	time = 0;
	int count = 0;
	result.timeResidueStl = 0;
	result.timeResidueNg = 0;
	while( !in.eof() ) {
		getline(in, pqrline);
		stringstream ss(pqrline.c_str());
		ss >> temp >> time;
		if (temp == "protein_surface"){
			result.timeProteinStl = time;
		}
		if (temp == "protein_volume"){
			result.timeProteinNg = time;
		}

		if (temp == "model_surface"){
			result.timeResidueStl += time;
			count++;
		}

		if (temp == "model_volume"){
			result.timeResidueNg += time;
			count++;
		}

	}
	result.timeResidueStl /= count;
	result.timeResidueNg /= count;
	in.close();
}

int main(int argc, char* argv[]) {

	string jobName, pqrFile, jobFile, workDir;

	if ( argc != 4 ) {
	    cout<<"usage: "<< argv[0] <<" <jobname> <pqr file> <job file> <working directory>\n";
	    exit(0);
	}
	else {
		jobName  = argv[1];
		pqrFile  = argv[2];
		jobFile	 = argv[3];
	    workDir  = argv[4];
	}

	cout << "Job name: " << jobName << endl;
	cout << "PQR file: " << jobFile << endl;
	cout << "Parsing job file: " << jobFile << endl;
	cout << "Working dir     : " << workDir << endl;

	ifstream in(jobFile.c_str());

	if (!in) {
		in.close();
		cout << "Cannot open job file:" << jobFile << endl;
		exit(0);
	}
	string pqrline;

	vector<Parameter> pList;
	vector<Result> results;

	while( !in.eof() ) {
		getline(in, pqrline);
		Parameter p;
		p.pqr = pqrFile;
		p.jobname = jobName;
		stringstream ss(pqrline);
		ss >> p.plocalh;
		if (p.plocalh.substr(0,1) == " " || p.plocalh.substr(0.1) == "#")
			continue;

		ss >> p.pmaxh >> p.pminh >> p.lsmsN;
		ss >> p.pgrading >> p.poptsteps_3d >> p.poptsteps_2d;
		ss >> p.pfineness;

		ss >> p.taubin >> p.hc >> p.lap;

		ss >> p.blocalh >> p.bmaxh >> p.bgrading >> p.boptsteps_3d;
		ss >> p.boptsteps_2d >> p.bfineness >> p.solver;
		ss >> p.maxsteps >> p.order >> p.boundary ;
		pList.push_back(p);

		string uid = p.getIdentifier();
		string cmd;
		cmd = "mkdir -p "+workDir+"/"+uid;
		system(cmd.c_str());
		cmd = workDir+"/"+uid;
		chdir(cmd.c_str());
		cmd = "cp "+p.pqr+" "+p.boundary+" .";
		system(cmd.c_str());

		writeConfigMFES(p);
		writeMeshFiles(p);

		cmd = "mfes --ini config.ini | tee result.log";
		system(cmd.c_str());

		Result currentResult;

		VOL proteinSurface;
		proteinSurface.readNGVol("proteinSurface.vol");
//		cout << proteinSurface.elementList.size() << " elements read in." << endl;
//		cout << proteinSurface.pointList.size() << " points read in." << endl;


//		cout << endl << "proteinSurface.vol: " << endl;
		currentResult.gpOnSurface = proteinSurface.getNrPoints();
//		cout << "points on surface: " << currentResult.gpOnSurface << endl;
		currentResult.areaProtein = proteinSurface.getTotalArea();
//		cout << "current area: " << currentResult.areaProtein << endl;


//		cout << endl << "proteinVolume.vol: " << endl;
		VOL proteinVolume;
		proteinVolume.readNGVol("proteinVolume.vol");
		float total = proteinVolume.getNrPoints();
//		cout << "total points: " << total << endl;
		currentResult.gpProtein = total - currentResult.gpOnSurface;
//		cout << "points inner protein: " << currentResult.gpProtein << endl;

		currentResult.volumeProtein = proteinVolume.getTotalVolume();
//		cout << "current volume: " << currentResult.volumeProtein << endl;
//		cout << "current area: " << proteinVolume.getTotalArea() << endl;

//		cout << endl << "boundary.vol: " << endl;
		VOL boundarySurface;
		boundarySurface.readNGVol("boundary.vol");
		currentResult.gpOnBoundary = boundarySurface.getNrPoints();
//		cout << "points on surface: " << currentResult.gpOnBoundary << endl;
		currentResult.areaBoundary = boundarySurface.getTotalArea();
//		cout << "current area: " << currentResult.areaBoundary << endl;
		currentResult.volumeTotal = boundarySurface.getTotalVolume();
//		cout << "current total volume: " << currentResult.volumeTotal << endl;

//		cout << endl << "protein.vol: " << endl;
		VOL protein;
		protein.readNGVol("protein.vol");
		currentResult.gpTotal = protein.getNrPoints();
//		cout << "total points: " << currentResult.gpTotal << endl;

		currentResult.volumeBetween = currentResult.volumeTotal - currentResult.volumeProtein;
//		cout << "volume between protein <> boundary: " << currentResult.volumeBetween << endl;

		currentResult.densityTotal = currentResult.gpTotal/currentResult.volumeTotal;
		currentResult.densityProtein = currentResult.gpOnSurface/currentResult.areaProtein;
		currentResult.densityBoundary = currentResult.gpOnBoundary/currentResult.areaBoundary;
		currentResult.densityInProtein = currentResult.gpProtein/currentResult.volumeProtein;
		currentResult.densityBetween = (currentResult.gpTotal - currentResult.gpProtein) / (currentResult.volumeTotal - currentResult.volumeProtein);

		currentResult.gpBetween = currentResult.gpTotal - currentResult.gpProtein - currentResult.gpOnSurface - currentResult.gpOnBoundary;

//		cout << endl << "total density: " << currentResult.densityTotal << endl;
//		cout << "density on protein surface: " << currentResult.densityProtein << endl;
//		cout << "density on boundary surface: " << currentResult.densityBoundary << endl;
//		cout << "density in protein volume: " << currentResult.densityInProtein << endl;
//		cout << "density between protein <> boundary: " << currentResult.densityBetween << endl;

		parse("result.log", currentResult);
		results.push_back(currentResult);

		currentResult.print();

		cmd = "../../";
		chdir(cmd.c_str());
	}

	in.close();
	cout << pList.size() << " job(s) processed." << endl;
	cout << results.size() << " result(s) collected." << endl;

	for(unsigned int i = 0; i < results.size(); i++){
		Result currentResult = results.at(i);
		cout << currentResult.getCSVLine(pList.at(i)) << endl;
	}

	return 0;
}
