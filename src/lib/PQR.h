#include <iostream>
#include <string>
#include <vector>

#include "Atom.h"
#include "Residue.h"
#include "Model.h"
#include "ST.h"
#include "tCycle.h"

using namespace std;


class PQR {
public:
	PQR(string _fileName):
	    fileName(_fileName)
	{
//	      parse();
	}

	static bool STparsed;
	static bool patchNTE;
	static bool patchCTE;
	static bool explicitModels;

private:
  string fileName;
  string stFilesFolder;
  Model molecule;
  vector<Atom> atomList;
  vector<Residue> residueList;
  vector<Residue> titGroupList;
  vector<ST> stList;
  tCycle cycle;
  vector< vector< vector< vector<double> > > > wmatrix;
};

bool PQR::STparsed = false;
bool PQR::explicitModels = false;
bool PQR::patchCTE = false;
bool PQR::patchNTE = false;
