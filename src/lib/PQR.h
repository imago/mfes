#include <solve.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <set>

#include "Atom.h"
#include "Charge.h"
#include "Residue.h"
#include "Model.h"
#include "ST.h"
#include "tCycle.h"
#include "PDE.h"

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
			string titGroupName = it->path().stem().string();
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
		if (ini.get<string>("experiment.cavity") == "yes")
			molecule.calcModel(atomList, ini, "cavity.vol", true);

		molecule.calcModel(atomList, ini);

	}

	void calcDeltaG(){
		cout << "Calculating deltaG ..." << endl;
		ngsolve::PDE pde;

	//	using namespace nglib;

		string pdeFile = "energy.pde";

		try {
		    pde.LoadPDE (pdeFile.c_str());
		    pde.Solve();
		} catch(ngstd::Exception & e) {
			std::cout << "Caught exception: " << std::endl
				<< e.What() << std::endl;
		}

	}

	void calcResidues(INI &ini){
	  Residue newResidue;
	  vector<Atom> currentAtomList;
	  Atom lastAtom;
	  int oldResID = -1;
	  int currentResID = oldResID;

	  for (unsigned int i = 0; i < atomList.size(); i++){
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
//		      if (newResidue.isTitGroup()){
//		    	  titGroupList.push_back(newResidue);
//		      }
		      currentAtomList.clear();
		  }
		  currentAtomList.push_back(currentAtom);
		  lastAtom = currentAtom;
		  oldResID = currentResID;
	  }

	  // Read in sites file
	  string sitesFile = ini.get<string>("pka.sites_file");
	  set<string> sitesData;
	  ifstream in(sitesFile.c_str());
	  if (!in) {
		  in.close();
		  cout << "Cannot open sites file:" << sitesFile << endl;
		  exit(0);
	  }
	  string currentLine, _chainID, _resNumber, _resName;

	  while( !in.eof() ) {
		  getline(in, currentLine);
		  stringstream ss(currentLine);
		  ss >> _chainID >> _resNumber >> _resName ;
		  string name = _resName+"-"+_resNumber+"_"+_chainID;
		  sitesData.insert(name);

	  }
	  // Insert residues into titGroupList
	  for (unsigned int i = 0; i < residueList.size(); i++){
		  Residue currentResidue = residueList.at(i);
		  if (currentResidue.isTitGroup() && sitesData.find(currentResidue.getIdentifier()) != sitesData.end()){
			  vector<Atom> currentAtomList = currentResidue.getAtomList();
			  if (i >= 1){
				  Residue before = residueList.at(i-1);
				  vector<Atom> newAtomList = before.getAtomList();
				  unsigned int atomListSize = newAtomList.size();
				  Atom C = newAtomList.at(atomListSize-2);
				  Atom O = newAtomList.at(atomListSize-1);;
				  currentAtomList.insert(currentAtomList.begin(), O);
				  currentAtomList.insert(currentAtomList.begin(), C);
			  }

			  if (i+1 < residueList.size()){
				  Residue after = residueList.at(i+1);
				  vector<Atom> newAtomList = after.getAtomList();
				  Atom N = newAtomList.at(0);
				  Atom HN = newAtomList.at(1);
				  Atom CA = newAtomList.at(2);
				  currentAtomList.insert(currentAtomList.end(), N);
				  currentAtomList.insert(currentAtomList.end(), HN);
				  currentAtomList.insert(currentAtomList.end(), CA);
			  }
			  currentResidue.setAtomList(currentAtomList);
			  titGroupList.push_back(currentResidue);
		  }
	  }

	   // Patching also N-terminus?
	   if (calcNTE){
	     Residue firstResidue = residueList.at(0);
	     firstResidue.setResidueName("NTE");
	     int residueNumber = firstResidue.getResidueNumber();
	     firstResidue.setResidueNumber(residueNumber--);
	     vector<Atom> currentAtomList = firstResidue.getAtomList();

	     Residue after = residueList.at(1);
	     vector<Atom> newAtomList = after.getAtomList();
	     Atom N = newAtomList.at(0);
	     Atom HN = newAtomList.at(1);
	     Atom CA = newAtomList.at(2);
	     for (unsigned int i = 0; i < currentAtomList.size(); i++)
	    	 currentAtomList.at(i).setResidueName("NTE");
	     currentAtomList.insert(currentAtomList.end(), N);
	     currentAtomList.insert(currentAtomList.end(), HN);
	     currentAtomList.insert(currentAtomList.end(), CA);
	     firstResidue.setAtomList(currentAtomList);
	     titGroupList.insert(titGroupList.begin(), firstResidue);
	   }

	   // Patching also ending C-terminus?
	   if (calcCTE){
	     Residue lastResidue = residueList.at(residueList.size()-1);
	     lastResidue.setResidueName("CTE");
	     int residueNumber = lastResidue.getResidueNumber();
	     lastResidue.setResidueNumber(residueNumber++);
	     vector<Atom> currentAtomList = lastResidue.getAtomList();

	     Residue before = residueList.at(residueList.size()-2);
	     vector<Atom> newAtomList = before.getAtomList();
	     unsigned int atomListSize = newAtomList.size();
	     Atom C = newAtomList.at(atomListSize-2);
	     Atom O = newAtomList.at(atomListSize-1);
	     for (unsigned int i = 0; i < currentAtomList.size(); i++)
	    	 currentAtomList.at(i).setResidueName("CTE");
	     currentAtomList.insert(currentAtomList.begin(), O);
	     currentAtomList.insert(currentAtomList.begin(), C);
	     lastResidue.setAtomList(currentAtomList);
	     titGroupList.insert(titGroupList.end(), lastResidue);
	   }

	   cout << residueList.size() << " residue(s) with " << titGroupList.size() << " titratable group(s) found." << endl;

	}

	void calcExplicitModels(INI &ini){
		for (unsigned int i = 0; i < titGroupList.size(); i++){
			Residue currentTitGroup = titGroupList.at(i);
			string fname = currentTitGroup.getIdentifier()+".vol";
			vector<Atom> currentAtomList = currentTitGroup.getAtomList();
			molecule.calcModel(currentAtomList, ini, fname);
			// Write out PQR files
			for (unsigned int j = 1; j <= currentTitGroup.getNrStates(); j++){
				writePQR(currentTitGroup, j);
			}
		}

	}

	void writePQR(Residue &titGroup, int stateNr){
		for (unsigned int i = 0; i < stList.size(); i++){
			ST currentST = stList.at(i);
			if (currentST.getTitGroupName() == titGroup.getResidueName() &&
					currentST.getStateNr() == stateNr){
					vector<Atom> currentAtomList = titGroup.getAtomList();
					patch(currentAtomList, currentST);
					stringstream ss;
					ss << titGroup.getIdentifier() << ".state" << stateNr << ".pqr";
					ofstream pqr;
					pqr.open (ss.str());
					for (unsigned int j = 0; j < currentAtomList.size(); j++){
						Atom currentAtom = currentAtomList.at(j);
						pqr << currentAtom.pqrLine() << endl;
					}
					pqr.close();
				break;
			}
		}
	}

	void switchState(Residue &titGroup, int stateNr){
		for (unsigned int i = 0; i < stList.size(); i++){
			ST currentST = stList.at(i);
			if (currentST.getTitGroupName() == titGroup.getResidueName() &&
					currentST.getStateNr() == stateNr){
					vector<Atom> currentAtomList = titGroup.getAtomList();
					patch(currentAtomList, currentST);
					titGroup.setAtomList(currentAtomList);
				break;
			}
		}

	}

	void patch(vector<Atom> &atomList, ST &st){
		vector<Charge> rules = st.getRules();
		for (unsigned int i = 0; i < atomList.size(); i++){
			Atom currentAtom = atomList.at(i);
			bool found = false;
			for (unsigned int j = 0; j < rules.size(); j++){
				Charge currentRule = rules.at(j);
				if (currentRule.getAtomName() == currentAtom.getAtomName() &&
						st.getTitGroupName() == currentAtom.getResidueName()){
					currentAtom.setCharge(currentRule.getCharge());
					found = true;
					break;
				}
			}
			if (!found){
				currentAtom.setCharge(0);
			}
			atomList.at(i) = currentAtom;
		}
	}

	void writeOutW(string fileName){
		cout << "Writing out W matrix .. " << endl;
		 ofstream out( fileName.c_str() );
		 out.setf(ios_base::scientific);
		 out.precision(6);

		 cout << "Writing W matrix into '" << fileName << "' ... ";
		 cout.flush();

		 if (!out)
		 {
		  out.close();
		  cout << "error writing: " << fileName << endl;
		  exit(0);
//		  throw FileIOException(filename, "PQR::writeOutW");
		 }

		 // run over all matrix elements
		 for(unsigned int i=0; i<titGroupList.size(); i++)
		 {
		  for(unsigned int j=0; j<titGroupList.size(); j++)
		  {
		   for(unsigned int state1 = 1; state1<titGroupList.at(i).getNrStates(); state1++)
			{
			 for(unsigned int state2 = 1; state2<titGroupList.at(j).getNrStates(); state2++)
			 {
		     out.width(4);
		     out << i+1 << " " << state1+1 << " ";
		     out.width(4);
		     out << j+1 << " " << state2+1 << "    ";
			  out << wmatrix[i][j][state1-1][state2-1] << endl;
			 }
			}
		  }
		 }

		 out.close();

		 cout << "finished." << endl;

	}

	double calcWmv(int state, Residue &mstat, Residue &mref, string cycle)
	{

		double wmv = 0;
		double wmv_c0 = 0;
		double wmv_c1 = 0;

		tCycle pp_s0 = mref.getTCycle();
		tCycle pp_s1 = mstat.getTCycle();
		Atom charge_s0, charge_s1;

		switchState(mstat, state);

		vector<Atom> mrefAtomList = mref.getAtomList();
		vector<Atom> mstatAtomList = mstat.getAtomList();

		// mstat size equals to mref size
		unsigned int i = 0;
		while ((mrefAtomList.size() > i) & (mstatAtomList.size() > i)){
			Atom charge_s0 = mrefAtomList.at(i);
			Atom charge_s1 = mstatAtomList.at(i);
			i++;

			if (charge_s0.getCharge() == 0 && charge_s1.getCharge() == 0)
				continue;

			if (cycle == "cycle0" || cycle == "diff")
				wmv_c0 += (charge_s1.getCharge() * (
					// which cycle?
						pp_s1.m[state][charge_s1.getCoord()] - pp_s0.m[0][charge_s1.getCoord()]
					)
					- charge_s0.getCharge() * (
							pp_s1.m[state][charge_s0.getCoord()] - pp_s0.m[0][charge_s0.getCoord()]
					)
				);
			if (cycle == "cycle1" || cycle == "diff")
				wmv_c1 += (charge_s1.getCharge() * (
					// which cycle?
						pp_s1.p[state][charge_s1.getCoord()] - pp_s0.p[0][charge_s1.getCoord()]
					)
					- charge_s0.getCharge() * (
							pp_s1.p[state][charge_s0.getCoord()] - pp_s0.p[0][charge_s0.getCoord()]
					)
				);


			wmv += wmv_c0 + wmv_c1;

		}

		return wmv;
	}

	void calcW(string cycle){
		cout << "Calculation of W matrix started.. " << endl;
		double G;
		double convert = T*MEADUNITS;

		 // create wmatrix
		 for(unsigned int i = 0; i<titGroupList.size(); i++)
		 {
		  wmatrix.push_back( vector< vector< vector<double> > >() );
		  for(unsigned int j = 0; j<titGroupList.size(); j++)
		  {
		   wmatrix.at(i).push_back( vector< vector<double> >() );
			/* there are no matrix entries for reference states
			 that's why we take only the number of states-1 */
			for(unsigned int state1 = 0; state1<titGroupList.at(i).getNrStates()-1; state1++)
			{
			 wmatrix.at(i).at(j).push_back( vector<double>() );
			 for(unsigned int state2 = 0; state2<titGroupList[j].getNrStates()-1; state2++)
			 {
			  wmatrix.at(i).at(j).at(state1).push_back(0.0);
			 }
			}
		  }
		 }

		 // have to make diff of both cycles or just read one cycle
		 string currentCycle = "cycle0";
		 for(unsigned int i = 0; i<titGroupList.size(); i++)
		 {
			 Residue currentTitGroup_i = titGroupList.at(i);

		  // read potentials of residue 1
			 stringstream currentPotat;
			 currentPotat << currentCycle << "." << currentTitGroup_i.getIdentifier() << ".potat";
			 cout << currentPotat.str() << endl;
		  for(unsigned int j = 0; j<titGroupList.size(); j++)
		  {
			  Residue currentTitGroup_j = titGroupList.at(j);
			// diagonal elements to be zero
			if ( currentTitGroup_i.getResidueName() == currentTitGroup_j.getResidueName() )
			{
		    for(unsigned int state1 = 1; state1<currentTitGroup_i.getNrStates(); state1++)
		    {
		     for(unsigned int state2 = 1; state2<currentTitGroup_j.getNrStates(); state2++)
		     {
		      wmatrix[i][j][state1-1][state2-1] = 0.0;
			  }
			 }
			 continue;
			}

			// Switch titGroup i to reference state
			switchState(currentTitGroup_i, 0);

/*			// charges of residue 2
			residues[j]->select(pqr);
			pref  = pqr->fromSelection();
		   pstat = pqr->fromSelection();
			pqr->select_clear();

			// reference state
			pref->clear_all_charges();
			residues[j]->switchToReferenceState(pref);

		   // charged state
		   pstat->clear_all_charges();
*/
		   // run over all non-reference states
		   for(unsigned int state1 = 1; state1<currentTitGroup_i.getNrStates(); state1++)
			{
		    for(unsigned int state2 = 1; state2<currentTitGroup_j.getNrStates(); state2++)
		 	 {
		     // get charge of current state in residue 2
		    	// Switch titGroup j to current state2 (!= (reference which is 0))
		    	switchState(currentTitGroup_j, state2);
//			  residues[j]->switchToState(state2, pstat);


			  // Wmv
			  G = calcWmv(state1, currentTitGroup_j, currentTitGroup_i, cycle);

			  G = G * convert; // kt->e^2/A

			  // update matrix
			  wmatrix[i][j][state1-1][state2-1] = G;
			 }
			}

		  }
		 }

	    // symmetrize matrix
		 double dev=0.0, w1=0.0, w2=0.0, ave=0.0;

		 cout << "Symmetrizing W matrix ... ";

		 // setup output stream for symmetry check
		 ostringstream log;
		 log.setf(ios_base::scientific);
		 log.precision(6);

		 for(unsigned int i=0; i<titGroupList.size(); i++)
		 {
		  for(unsigned int j=i+1; (j>i) && (j<titGroupList.size()); j++)
		  {
		   for(unsigned int state1 = 1; state1<titGroupList.at(i).getNrStates(); state1++)
			{
			 for(unsigned int state2 = 1; state2<titGroupList.at(j).getNrStates(); state2++)
			 {
				 w1  = wmatrix[i][j][state1-1][state2-1];
				 w2  = wmatrix[j][i][state2-1][state1-1];
				 dev = (w1-w2);
				 ave = 0.5*(w1+w2);
				 wmatrix[i][j][state1-1][state2-1] = ave;
				 wmatrix[j][i][state2-1][state1-1] = ave;
				 // write symmetry check to log file
				 log << titGroupList.at(i).getResidueName() << "[" << state1+1 << "]/";
				 log << titGroupList.at(j).getResidueName() << "[" << state2+1 << "]: ";
				 log << "(" << w1 << "," << w2 << "), diff = " << dev;
				 log << ", ave = " << ave << endl;
			 }
			}
		  }
		 }

		 cout << "finished." << endl;
//		 cout << log.str() << endl;
	}


	void calcPkint(string cycleName){

		for (unsigned int i = 0; i < titGroupList.size(); i++){
			Residue currentTitGroup = titGroupList.at(i);
			string fileName = cycleName+"."+currentTitGroup.getIdentifier()+".potat";
			cout << "Reading " << fileName << " .... ";

			tCycle cycle;
			unsigned int nrStates = 0;
			unsigned int nrAtoms = 0;
			float potential = 0;
			unsigned int potNr = 0;

			ifstream in(fileName.c_str());
			if (!in) {
				in.close();
				cout << "Cannot open potat file:" << fileName << endl;
				exit(0);
			}
			string currentLine;
			while( !in.eof() ) {
				getline(in, currentLine);
				nrStates = atoi(currentLine.c_str());
				getline(in, currentLine);
				nrAtoms  = atoi(currentLine.c_str());
				for(unsigned int j = 0; j < nrStates; j++){
					for(unsigned int k = 0; k < nrAtoms; k++){
						Point3<float> newPoint;
						getline(in, currentLine);
						newPoint.X() = atof(currentLine.c_str());
						getline(in, currentLine);
						newPoint.Y() = atof(currentLine.c_str());
						getline(in, currentLine);
						newPoint.Z() = atof(currentLine.c_str());
						getline(in, currentLine);
						potential = atof(currentLine.c_str());
						cycle.p[j].insert( {newPoint, potential} );
						potNr++;
					}
					getline(in, currentLine);
				}
				nrStates = 0;
				nrAtoms = 0;
				nrStates = atoi(currentLine.c_str());
				getline(in, currentLine);
				nrAtoms  = atoi(currentLine.c_str());
				for(unsigned int i = 0; i < nrStates; i++){
					for(unsigned int j = 0; j < nrAtoms; j++){
						Point3<float> newPoint;
						getline(in, currentLine);
						newPoint.X() = atof(currentLine.c_str());
						getline(in, currentLine);
						newPoint.Y() = atof(currentLine.c_str());
						getline(in, currentLine);
						newPoint.Z() = atof(currentLine.c_str());
						getline(in, currentLine);
						potential = atof(currentLine.c_str());
						cycle.m[i].insert( {newPoint, potential} );
						potNr++;
					}
					getline(in, currentLine);

				}
				cout << potNr << " potentials read in." << endl;

			}
			currentTitGroup.setTCycle(cycle);
			titGroupList.at(i) = currentTitGroup;
			in.close();

		}

		float bornEner = 0;
		float backEner = 0;

		for (unsigned int i = 0; i < titGroupList.size(); i++){
			Residue currentOrigTitGroup = titGroupList.at(i);
			Residue currentTitGroup = titGroupList.at(i);
			unsigned int nrStates = currentTitGroup.getNrStates();
			for (unsigned int stateNr = 1; stateNr < nrStates; stateNr++){
				tCycle cycle = currentTitGroup.getTCycle();

				vector<Atom> refAtomList = getAtoms(currentTitGroup, 0);
				vector<Atom> compareAtomList = getAtoms(currentTitGroup, stateNr);

				float born0 = 0;
				float born1 = 0;

			    for (unsigned int j = 0; j < refAtomList.size(); j++){
			    	Atom currentRefAtom = refAtomList.at(j);
			    	Atom currentCompareAtom = compareAtomList.at(j);
//			    	cout << "currentRefAtom :";
//			    	currentRefAtom.print();
//			    	cout << "currentCompareAtom :";
//			    	currentCompareAtom.print();

			    	if (currentRefAtom.getCharge() != 0){
//						cout << "born0 old: " << born0 << endl;
						born0 += currentRefAtom.getCharge() * (cycle.p[0][currentRefAtom.getCoord()] - cycle.m[0][currentRefAtom.getCoord()]);
//						cout << currentRefAtom.getCharge() << " * (" << cycle.p[0][currentRefAtom.getCoord()] << " - " << cycle.m[0][currentRefAtom.getCoord()] << ") = " << currentRefAtom.getCharge() * (cycle.p[0][currentRefAtom.getCoord()] - cycle.m[0][currentRefAtom.getCoord()]) << endl;
//						cout << "born0 new is:" << born0 << endl;
			    	}

			    	if (currentCompareAtom.getCharge() != 0){
//			    		cout << "born1 old: " << born1 << endl;
						born1 += currentCompareAtom.getCharge() * (cycle.p[stateNr][currentCompareAtom.getCoord()] - cycle.m[stateNr][currentCompareAtom.getCoord()]);
//						cout << currentCompareAtom.getCharge() << " * (" << cycle.p[stateNr][currentCompareAtom.getCoord()] << " - " << cycle.m[stateNr][currentCompareAtom.getCoord()] << ") = " << currentCompareAtom.getCharge() * (cycle.p[stateNr][currentCompareAtom.getCoord()] - cycle.m[stateNr][currentCompareAtom.getCoord()]) << endl;
//						cout << "born1 new is:" << born1 << endl;
			    	}
			    }
		        bornEner = 0.5*(born1-born0);
//			    cout << "Returning 0.5 * (" << born1 << " - " << born0 << ") = " << bornEner << endl;

		        float back0 = 0;
		        float back1 = 0;

		        //Residue titGroup = titGroup !
		        vector<Atom> pAtomList = deleteAtoms(atomList, refAtomList); // from whole protein without patched atoms in currentTitGroup
		        vector<Atom> mAtomList = deleteAtoms(currentOrigTitGroup.getAtomList(), refAtomList); // not patched



		        for (unsigned int j = 0; j < pAtomList.size(); j++){
		        	Atom currentAtom = pAtomList.at(j);
//		        	cout << "currentAtom (p)";
//		        	currentAtom.print();

//		    		cout << "back1 old: " << born1 << endl;
		            back1 += currentAtom.getCharge() * (cycle.p[0][currentAtom.getCoord()] - cycle.p[stateNr][currentAtom.getCoord()]);
//					cout << currentAtom.getCharge() << " * (" << cycle.p[0][currentAtom.getCoord()] << " - " << cycle.p[stateNr][currentAtom.getCoord()] << ")" << endl;
//					cout << "back1 new is:" << born1 << endl;

		        }

		        for (unsigned int j = 0; j < mAtomList.size(); j++){
		        	Atom currentAtom = mAtomList.at(j);
//		        	cout << "currentAtom (m)";
//		        	currentAtom.print();

//		        	cout << "back0 old: " << born1 << endl;
		        	back0 += currentAtom.getCharge() * (cycle.m[0][currentAtom.getCoord()] - cycle.m[stateNr][currentAtom.getCoord()]);
//					cout << currentAtom.getCharge() << " * (" << cycle.m[0][currentAtom.getCoord()] << " - " <<  cycle.m[stateNr][currentAtom.getCoord()] << ")" << endl;
//					cout << "back0 new is:" << born1 << endl;
		        }

		        backEner = (back1-back0);
//		        cout << "Returning (" << back1 << " - " << back0 << ") = " << backEner << endl;

		        unsigned int cycleNr = 0;
		        if (cycleName == "cycle1")
		        	cycleNr = 1;

			    currentTitGroup.setBorn(bornEner, cycleNr, stateNr); // real titGroupList object
			    currentTitGroup.setBack(backEner, cycleNr, stateNr);
//			    cout << currentTitGroup.getResidueName() << " form R to current state " << stateNr << " has born energy of " <<
//			    		bornEner << " and back energy of " << backEner << "kJ/mol" << endl;
			    titGroupList.at(i) = currentTitGroup;

			}


		}



	}


	vector<Atom> deleteAtoms(vector<Atom> result, vector<Atom>& compareAtomList){
		// compareAtomList is << currentAtomList
		for (unsigned int i = 0; i < compareAtomList.size(); i++){
			Atom currentCompareAtom = compareAtomList.at(i);
			if (currentCompareAtom.getCharge() != 0){
				for (unsigned int j = 0; j < result.size(); j++){
					Atom currentAtom = result.at(j);
					if (currentAtom.getAtomNumber() == currentCompareAtom.getAtomNumber()){
						result.erase(result.begin()+j);
					}
				}
			}
		}
		return result;
	}

	vector<Atom> getAtoms(Residue currentTitGroup, int stateNr){
		vector<Atom> currentAtomList = currentTitGroup.getAtomList();

		for (unsigned int i = 0; i < stList.size(); i++){
			ST currentST = stList.at(i);
			if (currentST.getStateNr() == stateNr+1 && currentST.getTitGroupName() == currentTitGroup.getResidueName()){
				patch(currentAtomList, currentST);
			}
		}

		return currentAtomList;
	}

	void calcPotat(string cycleName){
		string pdeFile;
		try {
			if (cycleName == "cycle1"){
				ngsolve::PDE pde;

				pdeFile = "pka_"+cycleName+".pde";
				boost::filesystem::wpath file(pdeFile);

				if(boost::filesystem::exists(file)){
					pde.LoadPDE (pdeFile.c_str());
					pde.Solve();
				}

			} else {
				for (unsigned int i = 0; i < titGroupList.size(); i++){
					ngsolve::PDE pde;

					string prefix = titGroupList.at(i).getIdentifier();
					pdeFile = "pka_"+cycleName+"_"+prefix+".pde";
					boost::filesystem::wpath file(pdeFile);
					string potatFile = cycleName+"."+prefix+".potat";
					boost::filesystem::wpath potfile(potatFile);

					if(boost::filesystem::exists(file) && !boost::filesystem::exists(potfile)){
						pde.LoadPDE (pdeFile.c_str());
						pde.Solve();
					} else {
						cout << potatFile << ".. already calculated." << endl;
					}
				}
			}

		} catch(ngstd::Exception & e) {
			std::cout << "Caught exception: " << std::endl
					<< e.What() << std::endl;
		}



	}

	void writePDE(INI &ini, string mode = "explicit"){
		PDE currentPDE;
		if (mode == "explicit"){
			currentPDE.writeCycle0(titGroupList, ini);
			currentPDE.writeCycle1(fileName, titGroupList, ini);
		} else {
			currentPDE.writeEnergy(fileName, ini);
		}
	}

	void writeOutPkint(string cycleName, string fileName){
		 ofstream out( fileName.c_str() );
		 out.setf(ios_base::scientific);
		 out.precision(6);

		 if (!out)
		 {
			 out.close();
			 cout << "could not open: " << fileName << endl;
		 }


		  unsigned int cycleNr = 0;
		  if (cycleName == "cycle1")
			  cycleNr = 1;

		  for (unsigned int i = 0; i < titGroupList.size(); i++){
			  Residue currentTitGroup = titGroupList.at(i);
			  // write reference energy [kJ/mol] to pkint file
			  out  << 0.0 << " " << getState(currentTitGroup.getResidueName(), 0) << " ";

			  unsigned int nrStates = currentTitGroup.getNrStates();
			  for (unsigned int stateNr = 1; stateNr < nrStates; stateNr++){
				  float born, back, shift, result;
				  if (cycleName == "cycle0"){
					  born = currentTitGroup.getBorn(cycleNr, stateNr);
					  back = currentTitGroup.getBack(cycleNr, stateNr);
					  shift = getShift(currentTitGroup.getResidueName(), stateNr)*CONVERT;
					  result = shift - ( born - back );
				  } else if (cycleName == "cycle1") {
					  born = currentTitGroup.getBorn(cycleNr, stateNr);
					  back = currentTitGroup.getBack(cycleNr, stateNr);
					  shift = 0;
					  result = shift - ( born - back );
				  } else {
					  born = currentTitGroup.getBorn(1, stateNr) - currentTitGroup.getBorn(0, stateNr);
					  back = currentTitGroup.getBack(1, stateNr) - currentTitGroup.getBack(0, stateNr);
					  shift = getShift(currentTitGroup.getResidueName(), stateNr)*CONVERT;;
					  result = shift+(born - back);
				  }

				  cout << currentTitGroup.getIdentifier() << ": R -> " << getState(currentTitGroup.getResidueName(), stateNr) << endl;
				  cout << "Gborn: " << born << endl;
				  cout << "Gback: " << back << endl;
				  cout << "intrinsic pka = "<< shift << " + ( " << born << " - " << back << ") = " << result << " [kJ/mol]" << endl;
				  cout << "-> " << result/CONVERT << " []" << endl;

				  // write pkint [kJ/mol] and transition identifier to pkint file
				  out << result << " " << getState(currentTitGroup.getResidueName(), stateNr) << " ";

			  }

			  out << currentTitGroup.getResidueName() << endl;

		  }

	}

	char getState(string residueName, int stateNr){
		for (unsigned int i = 0; i < stList.size(); i++){
			ST currentST = stList.at(i);
			if (currentST.getTitGroupName() == residueName && currentST.getStateNr() == (stateNr+1)){
				return currentST.getState();
			}
		}
		return '?';
	}

	float getShift(string residueName, int stateNr){
		for (unsigned int i = 0; i < stList.size(); i++){
			ST currentST = stList.at(i);
			if (currentST.getTitGroupName() == residueName && currentST.getStateNr() == (stateNr+1)){
				return currentST.getShift();
			}
		}
		return 0;
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
