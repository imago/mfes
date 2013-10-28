#include <solve.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include "Atom.h"
#include "Charge.h"
#include "Residue.h"
#include "Model.h"
#include "ST.h"
#include "tCycle.h"

using namespace std;


namespace netgen
{
  int h_argc;
  char ** h_argv;
}

char **argv;
ngsolve::MyMPI mympi(0, argv);


typedef boost::property_tree::ptree INI;

class PQR {
public:
	PQR(string _fileName):
	    fileName(_fileName)
	{
	      if(parse())
	    	  cout << atomList.size() << " atom(s) found." << endl;
	      else {
	    	  cout << "An error occured while parsing file..." << endl;
	    	  exit(0);
	      }
	}

	void addInfo(INI &ini){
		string _stFilesFolder = ini.get<string>("pka.st_folder");
		stFilesFolder = _stFilesFolder;

		string _calcCTE = ini.get<string>("pka.calc_cte");
		if (_calcCTE == "yes")
			calcCTE = true;

		string _calcNTE = ini.get<string>("pka.calc_nte");
		if (_calcNTE == "yes")
			calcNTE = true;

		string _explicitModels = ini.get<string>("pka.explicit_models");
		if (_explicitModels == "yes")
			explicitModels = true;

	}

	void parseSTFolder(){
		boost::filesystem::path full_path( stFilesFolder );

		for ( boost::filesystem::directory_iterator it( full_path );
				it != boost::filesystem::directory_iterator(); ++it )
		{
			string titGroupName = it->path().filename().string();
			string ext = it->path().extension().string();
			boost::to_lower(ext);

			if ( boost::filesystem::is_regular_file( it->status() )
				&& ( ext == ".st" ))
			{
				cout << "Parsing... " << it->path().string() << endl;
				ifstream in(it->path().string().c_str());
				if (!in) {
					in.close();
					cout << "Cannot open st file:" << it->path().string() << endl;
					exit(0);
				}
				string currentLine, temp;

				char state; float deltaG; float preCalcShift;
				string atomName; double charge;
				vector<Charge> rules;
				unsigned int stateNr = 0;
				while( !in.eof() ) {
					getline(in, currentLine);
					if ( currentLine.find("ATOM", 0)==string::npos ) {
						if (stateNr > 0){
							 ST newST(titGroupName, state, stateNr, deltaG, rules);
							 stList.push_back(newST);
							 rules.clear();
						}
						istringstream ss(currentLine);
						ss >> deltaG >> temp>> state >> preCalcShift;
						stateNr++;
					} else {

						temp = currentLine.substr(12,4);
						istringstream ss(temp);
						ss >> atomName;

						temp = currentLine.substr(54,6);
						if( sscanf(temp.c_str(), "%lg", &(charge))!=1 ){
							cout << "could not load charge in: " << endl;
							cout << currentLine << endl;
							exit(0);
						}
						Charge newCharge(atomName, charge);
						rules.push_back(newCharge);
					}

				}

			}
		}

		cout << stList.size() << " patch rules added." << endl;

	}

	void parseSTFiles(){
		parseSTFolder(); // we have all ST file information
	    STparsed = true;
	}

	int parse(){
		vector<Atom> atomList;
		string temp;
		double x, y, z;
 
		ifstream in(fileName.c_str());

		if (!in) {
		  in.close();
		  cout << "Cannot open molecule file:" << fileName << endl;
		  exit(0);
		}
		string pqrline;

		while( !in.eof() ) {
			getline(in, pqrline);
			if ( pqrline.find("ATOM", 0)==string::npos ) continue;

			 if ( pqrline.size() < 72 ){
				 cout << "Molecule file is not well structured. Example: " << endl;
				 cout << "ATOM      1  N   VAL A   3      16.875  37.901  43.478-0.470 1.850      A" << endl;
				 exit(0);
			 }

			int atomNumber; string atomName; string residueName; char confID; char chainID;
			int residueNumber; double charge; double radius; string segName;

			// atom number
			temp = pqrline.substr(6,5);
			if( sscanf(temp.c_str(), "%i", &(atomNumber))!=1 )
				return 0;

			// atom name
			temp = pqrline.substr(12,4);
			istringstream ss(temp);
			ss >> atomName;

			// conformer id
			temp = pqrline.substr(16,1);
			if( sscanf(temp.c_str(), "%c", &(confID))!=1 )
				confID = ' ';

			// residue name
			residueName = pqrline.substr(17,3);

			// chain id
			temp = pqrline.substr(21,1);
			if( sscanf(temp.c_str(), "%c", &(chainID))!=1 )
			chainID = ' ';

			// residue number
			temp = pqrline.substr(22,4);
			if( sscanf(temp.c_str(), "%i", &(residueNumber))!=1 )
				return 0;

			// coordinates
			double x, y, z;
			temp = pqrline.substr(30,8);
			if( sscanf(temp.c_str(), "%lg", &x)!=1 )
				return 0;

			temp = pqrline.substr(38,8);
			if( sscanf(temp.c_str(), "%lg", &y)!=1 )
				return 0;
			temp = pqrline.substr(46,8);
			if( sscanf(temp.c_str(), "%lg", &z)!=1 )
				return 0;

			Point3<float> coord(x, y, z);

			// charge
			temp = pqrline.substr(54,6);
			if( sscanf(temp.c_str(), "%lg", &(charge))!=1 )
				return 0;
 
			// radius
			temp = pqrline.substr(60,6);
			if( sscanf(temp.c_str(), "%lg", &(radius))!=1 )
				return 0;

			// segment name
			try
			{
				temp = pqrline.substr(72,4);
				istringstream ss(temp);
				ss >> segName;
			}
			catch(out_of_range& x)
			{
				segName = "    ";
			}
			Atom newAtom(atomNumber, atomName, confID, residueName, chainID, residueNumber, coord, charge, radius, segName);
			this->atomList.push_back(newAtom);
		}

		return 1;

	}

	void calcModel(INI& ini){
		molecule.calcModel(atomList, ini);
	}

	void calcDeltaG(){
		cout << "Calculating deltaG ..." << endl;
		ngsolve::PDE pde;

	//	using namespace nglib;

		string pdeFile = "test.pde";

		try {
		    pde.LoadPDE (pdeFile.c_str());
		    pde.Solve();
		  } catch(ngstd::Exception & e) {
		    std::cout << "Caught exception: " << std::endl
			      << e.What() << std::endl;
		  }
	}

	void calcResidues(){
	  Residue newResidue;
	  vector<Atom> currentAtomList;
	  Atom lastAtom;
	  int oldResID = -1;
	  int currentResID = oldResID;

	  for (int i = 0; i < atomList.size(); i++){
		  Atom currentAtom = atomList.at(i);
		  currentResID = currentAtom.getResidueNumber();
		  if (oldResID == -1){ // First element
			  oldResID = currentResID;
			  lastAtom = currentAtom;
		  }
		  //if (currentAtom is last element){ // Last element
		  if (i == (atomList.size()-1)){
			  oldResID = -2;
			  currentAtomList.push_back(currentAtom);
		  }

		  if (oldResID != currentResID){
			  newResidue.setResidueName(lastAtom.getResidueName());
			  newResidue.setResidueNumber(lastAtom.getResidueNumber());
		      newResidue.setChainID(lastAtom.getChainID());
		      newResidue.setAtomList(currentAtomList);
		      residueList.push_back(newResidue);
		      if (newResidue.isTitGroup()){
		    	  titGroupList.push_back(newResidue);
		      }
		      currentAtomList.clear();
		  }
		  currentAtomList.push_back(currentAtom);
		  lastAtom = currentAtom;
		  oldResID = currentResID;
	  }

	   // Patching also N-terminus?
	   if (calcNTE){
	     Residue firstResidue = residueList.at(0);
	     firstResidue.setResidueName("NTE");
	     int residueNumber = firstResidue.getResidueNumber();
	     firstResidue.setResidueNumber(residueNumber--);
	  //   residueList.insert(residueList.begin(), firstResidue);
	     titGroupList.insert(titGroupList.begin(), firstResidue);
	   }

	   // Patching also ending C-terminus?
	   if (calcCTE){
	     Residue lastResidue = residueList.at(residueList.size()-1);
	     lastResidue.setResidueName("CTE");
	     int residueNumber = lastResidue.getResidueNumber();
	     lastResidue.setResidueNumber(residueNumber++);
	   //  residueList.insert(residueList.end(), lastResidue);
	     titGroupList.insert(titGroupList.end(), lastResidue);
	   }

	   cout << residueList.size() << " residue(s) with " << titGroupList.size() << " titratable group(s) found." << endl;

	}
	static bool STparsed;
	static bool calcNTE;
	static bool calcCTE;
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
bool PQR::calcCTE = false;
bool PQR::calcNTE = false;
