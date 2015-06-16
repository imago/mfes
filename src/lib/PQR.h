/**
 * @file PQR.h
 *
 * @brief      PQR model, energy and pKA computation management. Helper methods.
 * @details    PQR file (karlsberg style) is read in. Parsing atoms. Organizing into
 *	       vector<Atom> atomList, vector<Residue> residueList and vector<Residue> 
 *             titGroupList. Managing volume model generation. Managing potential 
 *             energy or pKA computations.
 *             
 * @authors    Ilkay Sakalli
 *
 */


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

/// namespace defined for NETGEN
namespace netgen
{
  int h_argc;
  char ** h_argv;
}

char **argv;
ngsolve::MyMPI mympi(0, argv);

/// Using INI files with the boost library (property_tree)
typedef boost::property_tree::ptree INI;

/**
 * @class PQR
 *
 * @brief PQR class model to handle proteins (atoms and titratable groups).
 *
 * Reading options in config file. Parsing ST files (description of titratable 
 * groups in proteins). Titrating N- and/or C-terminus. Patching groups. Computing 
 * potential energies. pKA values (W-matrix, pkint PKINT files) using one- or 
 * two-cycle approach. Managing model generation (protein and titratable groups). 
 *
 * @author Ilkay Sakalli
 *
 */ 
class PQR {
 public:
  /** @brief Constructor
   *
   *  Parsing the PQR file and adding atoms to atomList.
   *
   *  @param string fileName of the PQR file
   */  
  
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
  
  /** @brief Setting up some options 
   *
   *  Here, some options are defined: ST files folder, calc termini?, 
   *  use two-cycle approach?
   *
   *  @param INI Options INI file.
   *  @return Void.
   */  
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
  

  /** @brief Parsing ST files folder.
   *
   *  ST files given in CHARMM format are parsed to know charges of 
   *  different titration states of residues.
   *
   *  @return Void.
   */  
  void parseSTFolder(){
    boost::filesystem::path full_path( stFilesFolder );
    
    for ( boost::filesystem::directory_iterator it( full_path );
	  it != boost::filesystem::directory_iterator(); ++it )
      {
	string titGroupName = it->path().stem().string();
	string ext = it->path().extension().string();
	boost::to_lower(ext);
	
	/// For every file in the ST files folder defined in INI config file
	/// If the file has extension "st"
	if ( boost::filesystem::is_regular_file( it->status() )
	     && ( ext == ".st" ))
	  {
	    /// Parse current ST file
	    cout << "Parsing... " << it->path().string() << endl;
	    ifstream in(it->path().string().c_str());
	    if (!in) {
	      in.close();
	      cout << "Cannot open st file:" << it->path().string() << endl;
	      exit(0);
	    }
	    string currentLine, temp;
	    
	    char state; float deltaG; float shiftCycle0 = NOT_IN_ST;
	    string atomName; double charge;
	    vector<Charge> rules;
	    unsigned int stateNr = 0;
	    while( !in.eof() ) {
	      getline(in, currentLine);
	      if ( currentLine.find("ATOM", 0)==string::npos ) {
		if (stateNr > 0){
		  ST newST(titGroupName, state, stateNr, deltaG, rules, shiftCycle0);
		  stList.push_back(newST);
		  rules.clear();
		}
		istringstream ss(currentLine);
		ss >> deltaG >> temp>> state >> temp >> shiftCycle0;
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
  
  
  /** @brief Parsing ST files folder wrapper.
   *
   *  ST files given in CHARMM format are parsed to know charges of 
   *  different titration states of residues. Afterwards static variable
   *  STparsed is set to true.
   *
   *  @return Void.
   */  
  void parseSTFiles(){
    parseSTFolder(); // we have all ST file information
    STparsed = true;
  }
  

  /** @brief Parsing the PQR file.
   *
   *  Adding all atoms of PQR file to atomList.
   *
   *  @return int 1 if parsing was successful.
   */  
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
  
  /** @brief Wrapper method to compute molecular volume model of the protein.
   *
   *  Compute molecular volume model of protein depending on the 
   *  configuration of the INI file.
   *
   *  @param INI Configuration INI file.
   *  @return Void.
   */  
  void calcModel(INI& ini){
    /// Compute the cavity
    if (ini.get<string>("experiment.cavity") == "yes")
      molecule.calcModel(atomList, ini, "cavity.vol", true);
    
    /// Compute the ion exclusion layer if ionc != 0
    boost::optional<string> ionc = ini.get_optional<string>("experiment.ionc");
    if(ionc){
      if(ini.get<float>("experiment.ionc") != 0){
	molecule.calcIonLayer(atomList, ini, "exclusion.vol");
      }
    }

    /// Compute the molecular surface of the protein
    molecule.calcModel(atomList, ini);
  }
  
  /** @brief Compute energy difference using NETGEN and nglib libraries.
   *
   *  An electrostatic solvation energy difference is computed. NETGEN and
   *  MUMPS is used.
   *
   *  @return Void.
   */  
  void calcDeltaG(){
    cout << "Calculating deltaG ..." << endl;
    ngsolve::PDE pde;
   
    string pdeFile = "energy.pde"; 

    try {
      pde.LoadPDE (pdeFile.c_str());
      pde.Solve();
    } catch(ngstd::Exception & e) {
      std::cout << "Caught exception: " << std::endl
		<< e.What() << std::endl;
    }
    
  }
  
  /** @brief Using the atom list from PQR file parse all residues.
   *
   *  All residues which define a molecular model and for every titratable group
   *  are computed in the spirit of the pKA framework of Knapp et al. Also,
   *  C- and N-termini are added if defined in config INI file.
   *  
   *  @param INI Configuration INI file.
   *  @return Void.
   */  
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
    /*	   cout << "titration groups ordering" << endl;
	   for (unsigned int i = 0; i < titGroupList.size(); i++){
	   Residue currentTitGroup = titGroupList.at(i);
	   cout << currentTitGroup.getResidueName() << endl;
	   }*/
    

  }

  /** @brief Molecular volume surfaces for titratable groups are computed if using
   *         the two-cycle approach.
   *
   *  Molecular volume models are computed for titratable groups if the two-cycle 
   *  approach (explicit model computation) is chosen. PQR files are written out 
   *  for every titratable group in every state.
   *  
   *  @param INI Configuration INI file.
   *  @return Void.
   */  
  void calcExplicitModels(INI &ini){
    for (unsigned int i = 0; i < titGroupList.size(); i++){
      Residue currentTitGroup = titGroupList.at(i);
      // Write out PQR files
      for (unsigned int j = 1; j <= currentTitGroup.getNrStates(); j++){
	writePQR(currentTitGroup, j);
      }
      
      if (getCycle0Shift(currentTitGroup.getResidueName(), 1) == NOT_IN_ST || explicitModels){
	string fname = currentTitGroup.getIdentifier()+".vol";
	vector<Atom> currentAtomList = currentTitGroup.getAtomList();
	boost::optional<string> ionc = ini.get_optional<string>("experiment.ionc");
	if(ionc){
	  if(ini.get<float>("experiment.ionc") != 0){
	    molecule.calcIonLayer(atomList, ini, currentTitGroup.getIdentifier()+".vol_exclusion");
	  }
	}
		
	molecule.calcModel(currentAtomList, ini, fname);
      } 
      
    }
    
  }


  /** @brief Write out a PQR file using a given atom list.
   *
   *  A PQR file is written out using a given atom list to a 
   *  file as fileName.
   *  
   *  @param vector<Atom> Atom List to write out to PQR file.
   *  @param string File name of the resulting PQR file.
   *  @return Void.
   */  
  void writePQR(vector<Atom> &aList, string fileName){
    ofstream pqr;
    pqr.open (fileName);
    for (unsigned int i = 0; i < aList.size(); i++){
      Atom currentAtom = aList.at(i);
      pqr << currentAtom.pqrLine() << endl;
    }
    pqr.close();
    
  }
  
  /** @brief Overloading method to write out a given residue with a given
   *         state to a PQR file.
   *
   *  A PQR file is written out using a given residue with state.
   *  
   *  @param Residue Current resdiue which holds a list of atoms.
   *  @param int State of residue to patch to.
   *  @return Void.
   */  
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
  
  /** @brief Switches the state of a residue to a given state nr.
   *
   *  A residue which has an arbitrary state is switched into another state
   *  with this one using previous ST file information and atom list information
   *  from given residue.
   *  
   *  @param Residue Current resdiue which holds a list of atoms.
   *  @param int State of residue to patch to.
   *  @return Void.
   */  
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
  
  /** @brief Switches a given residue into neutral state.
   *
   *  A residue which has an arbitrary state is switched into neutral/
   *  reference state using previous ST file information and atom list 
   *  information from given residue.
   *  
   *  @param Residue Current resdiue which holds a list of atoms.
   *  @param int State of residue to patch to neutral state.
   *  @return Void.
   */  
  void neutralState(Residue &titGroup, int stateNr){
    for (unsigned int i = 0; i < stList.size(); i++){
      ST currentST = stList.at(i);
      if (currentST.getTitGroupName() == titGroup.getResidueName() &&
	  currentST.getStateNr() == stateNr){
	vector<Atom> currentAtomList = titGroup.getAtomList();
	patchNeutral(currentAtomList, currentST);
	titGroup.setAtomList(currentAtomList);
	break;
      }
    }
    
  }
  

  /** @brief Switches a given atom list (of a residue) to another state (using ST
   *         data structure which is already in a given state). 
   *
   *  An atom list (normally from a residue) is patched. Atoms in the list are 
   *  patched using ST which has charge information. ST is normally one state
   *  of a titratable group.
   *  
   *  @param vector<Atom> An atom list.
   *  @param ST ST structure with information about charges.
   *  @return Void.
   */  
  void patch(vector<Atom> &atomList, ST &st){
    vector<Charge> rules = st.getRules();
    for (unsigned int i = 0; i < atomList.size(); i++){
      Atom currentAtom = atomList.at(i);
      bool found = false;
      for (unsigned int j = 0; j < rules.size(); j++){
	Charge currentRule = rules.at(j);
	if (currentRule.getAtomName() == currentAtom.getAtomName() &&
	    st.getTitGroupName() == currentAtom.getResidueName()){
	  float currentCharge = currentRule.getCharge();
	  if ( currentCharge == 0){
	    currentCharge = 0.01;
	  }
	  
	  currentAtom.setCharge(currentCharge);
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
  
  /** @brief Switches a given atom list (of a residue) to neutral state (using ST
   *         data structure). 
   *
   *  An atom list (normally from a residue) is patched. Atoms in the list are 
   *  patched using ST which has charge information.
   *  
   *  @param vector<Atom> An atom list.
   *  @param ST ST structure with information about charges.
   *  @return Void.
   */  
  void patchNeutral(vector<Atom> &atomList, ST &st){
    vector<Charge> rules = st.getRules();
    for (unsigned int i = 0; i < atomList.size(); i++){
      Atom currentAtom = atomList.at(i);
      bool found = false;
      for (unsigned int j = 0; j < rules.size(); j++){
	Charge currentRule = rules.at(j);
	if (currentRule.getAtomName() == currentAtom.getAtomName() &&
	    st.getTitGroupName() == currentAtom.getResidueName()){
	  float currentCharge = currentRule.getCharge();
	  // There should be a better way than this.
	  if ( currentCharge == 0){
	    currentCharge = 0.01;
	  }
	  
	  currentAtom.setCharge(currentCharge);
	  found = true;
	  break;
	}
      }
      if (!found){
	//currentAtom.setCharge(0);
      }
      //currentAtom.print();
      atomList.at(i) = currentAtom;
    }
  }
  

  /** @brief Writing out W-matrix as flat G file.
   *
   *  A symmetric W-matrix which is computed before is written out.
   *  
   *  @param string Prefix of filename for W-matrix
   *  @return Void.
   */  
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
//  throw FileIOException(filename, "PQR::writeOutW");
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
  
  /** @brief Computing the W-matrix element given a residue in reference 
   *         and current state and cycle information.
   *
   *  A W-matrix element is computed using the pKA framework defined by Knapp et al.
   *  
   *  @param int Current state number.
   *  @param Residue Residue in reference state (protein environment).
   *  @param Residue Residue in reference state (water environment).
   *  @param Residue Residue in current state (water environment).
   *  @param string Name of current cycle.
   *  @return double W-matrix element.
   */  
  double calcWmv(int state, Residue &ref, Residue &mref, Residue &mstat, string cycle)
  {
    
    // double wmv = 0;
    double wmv_c0 = 0;
    double wmv_c1 = 0;
    
    Atom charge_s0, charge_s1;
    
    vector<Atom> mstatAtomList = mstat.getAtomList();
    vector<Atom> mrefAtomList = mref.getAtomList();
    
    tCycle pp_s1 = ref.getTCycle();

    // mstat size equals to mref size
    unsigned int i = 0;
    
    
    cycle = "cycle1";
    while (mrefAtomList.size() > i){	   
      Atom charge_s1 = mstatAtomList.at(i);
      Atom charge_s0 = mrefAtomList.at(i);
      i++;
      
      if (charge_s0.getCharge() == 0 && charge_s1.getCharge() == 0 )
	continue;
      
      if (cycle == "cycle0" || cycle == "diff"){
	wmv_c0 += (charge_s1.getCharge() * (
					    pp_s1.m[state][charge_s1.getCoord()] 
					    - pp_s1.m[1][charge_s1.getCoord()]
					    ) 
		   - charge_s0.getCharge() * (
					      pp_s1.m[state][charge_s0.getCoord()] 
					      - pp_s1.m[1][charge_s0.getCoord()] 
					      )
		   );
      }
      
      if (cycle == "cycle1" || cycle == "diff"){
	//  cout << "Lookint at charge_s1 (" << charge_s1.getCoord().X() << ", " << charge_s1.getCoord().Y() << ", " << charge_s1.getCoord().Z() << ")" << endl;
	// cout << "Lookint at charge_s0 (" << charge_s0.getCoord().X() << ", " << charge_s0.getCoord().Y() << ", " << charge_s0.getCoord().Z() << ")" << endl;
	float temp = (charge_s1.getCharge() - charge_s0.getCharge()) * (
									pp_s1.p[state][charge_s1.getCoord()] - pp_s1.p[0][charge_s1.getCoord()]
									);
	wmv_c1 += (charge_s1.getCharge() - charge_s0.getCharge()) * (
								     pp_s1.p[state][charge_s1.getCoord()] - pp_s1.p[0][charge_s1.getCoord()]
								     ); 
	
	// cout << "(" << charge_s1.getCharge() << " - " << charge_s0.getCharge() << ") * (" << pp_s1.p[state][charge_s1.getCoord()] << " - " << pp_s1.p[0][charge_s1.getCoord()] << ") = " << temp << endl;
	
      }
      
    }
    // cout << "wmv_c1 = " << wmv_c1 << ", wmv_c0 = " << wmv_c0 << endl;
    // cout << "returning wmv = " << wmv_c1 << endl;
    return (wmv_c1 - wmv_c0);
  }
  
  /** @brief Computing the whole W-matrix for every pairwise titratable group
   *         interaction.
   *
   *  A multidimensional W-matrix is generated. At first diagonal elements are 
   *  set to zero. Decide which cycle is computed and compute every W element 
   *  using calcWmv of one to another group for every state. Then the matrix is
   *  symmetrized.
   *  
   *  @param string Name of current cycle.
   *  @return Void.
   */  
  void calcW(string cycle){
    cout << "Calculation of W matrix started.. " << endl;
    double G;
    
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
    for(unsigned int i = 0; i<titGroupList.size(); i++)
      {
	Residue currentTitGroup_i = titGroupList.at(i);
	// cout << "i has name: " << currentTitGroup_i.getResidueName() << endl;
	// read potentials of residue 1
	for(unsigned int j = 0; j<titGroupList.size(); j++)
	  {
	    // cout << "at i,j: " << i << ", " << j << endl;
	    Residue currentTitGroup_j = titGroupList.at(j);
	    Residue currentTitGroup_j_Ref = titGroupList.at(j);
	    //  cout << "j has name: " << currentTitGroup_j.getResidueName() << endl;
	    // diagonal elements to be zero
	    if ( currentTitGroup_i.getResidueName() == currentTitGroup_j.getResidueName() && i == j)
	      {
		//  cout << "both have same name: " << currentTitGroup_i.getResidueName() << endl;
		
		for(unsigned int state1 = 1; state1<currentTitGroup_i.getNrStates(); state1++)
		  {
		    for(unsigned int state2 = 1; state2<currentTitGroup_j.getNrStates(); state2++)
		      {
			wmatrix[i][j][state1-1][state2-1] = 0.0;
		      }
		  }
		continue;
	      }
	    
	    
	    // run over all non-reference states
	    for(unsigned int state1 = 1; state1<currentTitGroup_i.getNrStates(); state1++)
	      {
		for(unsigned int state2 = 1; state2<currentTitGroup_j.getNrStates(); state2++)
		  {
		    
		    //  cout << "state1 is " << state1 << endl;
		    // cout << "state2 is " << state2 << endl;
		    
		    // get charge of current state in residue 2
		    // Switch titGroup j to current state2 (!= (reference which is 0))
		    // Switch titGroup i to reference state
		    switchState(currentTitGroup_i, 1);
		    switchState(currentTitGroup_j, state2+1);
		    switchState(currentTitGroup_j_Ref, 1);
		    
		    // Wmv
		    //						 G = calcWmv(state1, currentTitGroup_j, currentTitGroup_i, cycle);
		    
		    //			       cout << "G cycle0 is " << calcWmv(state1, pstat, pref, "cycle0");
		    //			       cout << "G cycle1 is " << calcWmv(state1, pstat, pref, "cycle1");
		    
		    G = calcWmv(state1, currentTitGroup_i, currentTitGroup_j_Ref, currentTitGroup_j, cycle); // in kJ/mol
		    
		    G *= E2A;
		    
		    //			       cout << currentTitGroup_i.getResidueName() << " " << i << " with " << currentTitGroup_j.getResidueName()  << " " << j ;
		    
		    //			       cout << "G: cycle1 - cycl0 difference is " << G << endl; 
		    
		    //			       cout << "writing at " << i << ", " << j << ", " << state1-1 << ", " << state2-1 << endl;
		    // update matrix
		    wmatrix[i][j][state1-1][state2-1] = G; 
		    
		  }
	      }
	  }
      }
    cout << "finished." << endl;
    
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
		    
		    
		    //				 cout << titGroupList.at(i).getResidueName() << " " << i << " with " << titGroupList.at(j).getResidueName() << " " << j;
		    //				 cout << "w1: " << w1 << ", w2: " << w2 << endl;
		    
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
  
  
  
  /** @brief Computing intrinsic pKA values.
   *
   *  Intrinsic pKA values (pkint) are computed parsing PKINT files computed before.
   *  Back and Born energy parts are assigned to titratable groups.
   *  
   *  @param string Name of current cycle.
   *  @param string A reference PQR definition is needed (written out before).
   *  @return Void.
   */  
  void calcPkint(string cycleName, string refPQRFileName){
    
    for (unsigned int i = 0; i < titGroupList.size(); i++){
      Residue currentTitGroup = titGroupList.at(i);
      if (cycleName == "cycle0" && (getCycle0Shift(titGroupList.at(i).getResidueName(), 1) != NOT_IN_ST && !explicitModels) ){
	tCycle cycle;
	currentTitGroup.setTCycle(cycle);
      } else {
	
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
    }
    
    vector<Atom> refAllAtomList;
    for (unsigned int i = 0; i < residueList.size(); i++){
      Residue currentResidue = residueList.at(i);
      
      if (i == 0 && calcNTE){
	Residue nteResidue = residueList.at(0);
	vector<Atom> list = nteResidue.getAtomList();
	for (unsigned int k = 0; k < list.size(); k++){ 
	  Atom currentAtom = list.at(k);
	  currentAtom.setResidueName("NTE");
	  list.at(k) = currentAtom;
	}
	nteResidue.setResidueName("NTE");
	nteResidue.setAtomList(list);
	neutralState(nteResidue, 1);
	list = nteResidue.getAtomList();
	for (unsigned int k = 0; k < list.size(); k++){ 
	  Atom currentAtom = list.at(k);
	  currentAtom.setResidueName(currentResidue.getResidueName());
	  list.at(k) = currentAtom;	      
	}
	currentResidue.setAtomList(list);
      } else if ( i == residueList.size()-1  && calcCTE ){
	Residue cteResidue = residueList.at(residueList.size()-1);
	vector<Atom> list = cteResidue.getAtomList();
	for (unsigned int k = 0; k < list.size(); k++){ 
	  Atom currentAtom = list.at(k);
	  currentAtom.setResidueName("CTE");
	  list.at(k) = currentAtom;
	}
	cteResidue.setResidueName("CTE");
	cteResidue.setAtomList(list);
	neutralState(cteResidue, 1);
	list = cteResidue.getAtomList();
	for (unsigned int k = 0; k < list.size(); k++){ 
	  Atom currentAtom = list.at(k);
	  currentAtom.setResidueName(currentResidue.getResidueName());
	  list.at(k) = currentAtom;	      
	}
	currentResidue.setAtomList(list);
	
      }
      
      bool found = false;
      for (unsigned int j = 0; j < titGroupList.size(); j++){
	Residue currentTGroup = titGroupList.at(j);
	if (currentTGroup.getResidueName() == currentResidue.getResidueName() && currentTGroup.getResidueNumber() == currentResidue.getResidueNumber()){
	  found = true;
	  cout << "reference: setting " << currentTGroup.getIdentifier() << " to neutral state" << endl;
	  
	  break;
	} 
      }
      if (found){
	neutralState(currentResidue, 1);
      }
      vector<Atom> residueAtoms = currentResidue.getAtomList();
      refAllAtomList.insert(refAllAtomList.end(), residueAtoms.begin(), residueAtoms.end());
    }
    
    writePQR(refAllAtomList, refPQRFileName);
    
    
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
	  if (currentRefAtom.getCharge() != 0){
	    born0 += currentRefAtom.getCharge() * (cycle.p[0][currentRefAtom.getCoord()] - cycle.m[0][currentRefAtom.getCoord()]);
	  }
	  
	  if (currentCompareAtom.getCharge() != 0){
	    born1 += currentCompareAtom.getCharge() * (cycle.p[stateNr][currentCompareAtom.getCoord()] - cycle.m[stateNr][currentCompareAtom.getCoord()]);
	  }
	}
	bornEner = 0.5*(born1-born0);
	
	float back0 = 0;
	float back1 = 0;
	
	//				cout << "Calculating BACK" << endl;
	//				cout << "================" << endl;
	
	//Residue titGroup = titGroup !
	//				vector<Atom> pAtomList = deleteAtoms(refAllAtomList, refAtomList); // from whole protein without patched atoms in currentTitGroup
	vector<Atom> pAtomList = deleteAtoms(refAllAtomList, refAtomList);
	stringstream ss;
	ss << "pAtomList_" << currentTitGroup.getIdentifier() << ".pqr"; 
	writePQR(pAtomList, ss.str());
	
	
	vector<Atom> mAtomList = deleteAtoms(currentOrigTitGroup.getAtomList(), refAtomList); // not patched
	stringstream ss2;
	ss2 << "mAtomList_" << currentTitGroup.getIdentifier() << ".pqr"; 
	writePQR(mAtomList, ss2.str());
	
	
	for (unsigned int j = 0; j < pAtomList.size(); j++){
	  Atom currentAtom = pAtomList.at(j);
	  if (currentAtom.getCharge() != 0){
	    back1 += currentAtom.getCharge() * (cycle.p[0][currentAtom.getCoord()] - cycle.p[stateNr][currentAtom.getCoord()]);
	    if (cycle.p[0][currentAtom.getCoord()] != 0){
	      Point3<float> c = currentAtom.getCoord();
	      //		    cout << "Looking at pcharge (" << c.X() << ", " << c.Y() << ", " << c.Z() << endl;
	      // cout << "back1 = " << currentAtom.getCharge() << " * ( " << cycle.p[0][currentAtom.getCoord()] << " - " << cycle.p[stateNr][currentAtom.getCoord()] << " ) = " << back1 << endl;
	    }
	  }
	}
	
	for (unsigned int j = 0; j < mAtomList.size(); j++){
	  Atom currentAtom = mAtomList.at(j);
	  if (currentAtom.getCharge() != 0){
	    back0 += currentAtom.getCharge() * (cycle.m[0][currentAtom.getCoord()] - cycle.m[stateNr][currentAtom.getCoord()]);
	    if (cycle.m[0][currentAtom.getCoord()] != 0){
	      Point3<float> c = currentAtom.getCoord();
	      
	      //  cout << "Looking at mcharge (" << c.X() << ", " << c.Y() << ", " << c.Z() << endl;
	      //		    cout << "back2 = " << currentAtom.getCharge() << " * ( " << cycle.m[0][currentAtom.getCoord()] << " - " << cycle.m[stateNr][currentAtom.getCoord()] << " ) = " << back0 << endl;
	    }
	  }
	}
	backEner = (back1-back0);
	// cout << "returning (" << back1 << "-" << back0 << ") =" << backEner << endl; 
	
	unsigned int cycleNr = 0;
	if (cycleName == "cycle1")
	  cycleNr = 1;
	
	currentTitGroup.setBorn(bornEner, cycleNr, stateNr); // real titGroupList object
	currentTitGroup.setBack(backEner, cycleNr, stateNr);
	titGroupList.at(i) = currentTitGroup;
      }
    }
  }
  
  /** @brief Deletes from atom List "result" all atoms which are also in atom 
   *         list "compareAtomList" via their atom id in PQR file.
   *         
   *  Returns the vector "result" which contains all atoms which were in list
   *  "result" before and are not same as in "compareAtomList". 
   *  If same atoms are found in "result" and "compareAtomList" this aroms are 
   *  erased from "result".
   *  
   *  @param vector<Atom> Atom list where same atoms with compareAtomList will 
   *                      be erased.
   *  @param vector<Atom> Atom list where atoms from this list are erased in 
   *                      result list.
   *  @return vector<Atom> Returns an atom list where same atoms with compareAtomList
   *                       are erased from result list.
   *
   */   
  vector<Atom> deleteAtoms(vector<Atom> result, vector<Atom>& compareAtomList){
    // compareAtomList is << currentAtomList
    for (unsigned int i = 0; i < compareAtomList.size(); i++){
      Atom currentCompareAtom = compareAtomList.at(i);
      //			currentCompareAtom.print();
      if (currentCompareAtom.getCharge() != 0){
	unsigned int pos = 0;
	for (unsigned int j = 0; j < result.size(); j++){
	  Atom currentAtom = result.at(j);
	  if (currentAtom.getAtomNumber() == currentCompareAtom.getAtomNumber()){
	    result.erase(result.begin()+j-pos);
	    pos++;
	  }
	}
      }
    }
    return result;
  }
  /** @brief Returns a list of atoms from a given residue which is patched to
   *         a given state nr.
   *         
   *  A residue is patched to a given state by stateNr. Its atoms are returned
   *  as an atom list. 
   *  
   *  @param Residue A titratable group which will be patched.
   *  @param int State number to patch given residue.
   *  @return vector<Atom> Returns an atom list where given residue currentTitGroup
   *                       is patched to state stateNr.
   *
   */     
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


  /** @brief Returns the averaged pKA value which is defined for the one-cycle 
   *         approach to skip a second thermodynamic cycle computation.
   *         
   *  Averaged pKA values are given by Sakalli et al. and are stored in ST files.
   *  
   *  @param string Name of residue receive an averaged pKA value.
   *  @param int State of given residue to return an averaged pKA value.
   *  @return float Averaged pKA value to skip second thermodynamic cycle computation.
   */    
  float getCycle0Shift(string residueName, int state){
    float shift = NOT_IN_ST;
    for (unsigned int i=0; stList.size(); i++){
      ST currentST = stList.at(i);
      if (currentST.getTitGroupName() == residueName && currentST.getStateNr() == state){
	shift = currentST.getShiftCycle0();
	return currentST.getShiftCycle0();
      }
    }
    return shift;
  }
  
  /** @brief Computes POTAT files which are needed to calculate intrinsic pKA values.
   *         
   *  A PDE file is written out before to control NGSolve which is a module of NETGEN.
   *  With this, POTAT files are computed which store all electrostatic potentials at 
   *  every atom position of every titratable goup in every state.
   *  
   *  @param string Name of current cycle.
   *  @return Void.
   */    
  void calcPotat(string cycleName){
    string pdeFile;
    try {
      if (cycleName == "cycle1"){
	ngsolve::PDE pde;
	bool skip = false;
	pdeFile = "pka_"+cycleName+".pde";
	boost::filesystem::wpath file(pdeFile);
	
	if(boost::filesystem::exists(file) && !skip){
	  pde.LoadPDE (pdeFile.c_str());
	  pde.Solve();
	} else {
	  cout << "cycle1 computation skipped manually." << endl;		 		  }
	
	
      } else {
	for (unsigned int i = 0; i < titGroupList.size(); i++){
	  ngsolve::PDE pde;
	  if (getCycle0Shift(titGroupList.at(i).getResidueName(), 1) == NOT_IN_ST || explicitModels){
	    string prefix = titGroupList.at(i).getIdentifier();
	    pdeFile = "pka_"+cycleName+"_"+prefix+".pde";
	    boost::filesystem::wpath file(pdeFile);
	    string potatFile = cycleName+"."+prefix+".potat";
	    boost::filesystem::wpath potfile(potatFile);
	    
	    if(boost::filesystem::exists(file)){
	      //if(boost::filesystem::exists(file)
	      // && !boost::filesystem::exists(potfile)){
	      
	      pde.LoadPDE (pdeFile.c_str());
	      pde.Solve();
	    } 
	    //						else {
	    //	cout << potatFile << ".. already calculated." << endl;
	    //}
	  }
	}
      }
      
    } catch(ngstd::Exception & e) {
      std::cout << "Caught exception: " << std::endl
		<< e.What() << std::endl;
    }    
  }

  
  /** @brief Writes out PDE files to compute POTAT files with NGSolve.
   *         
   *  A PDE file is written out to control NGSolve which is a module of NETGEN.
   *  
   *  @param INI Config INI file.
   *  @param string Approach used to compute pKA values (two-cycle approach: 
   *                implicit; one-cycle approach: explicit).
   *  @return Void.
   */    
  void writePDE(INI &ini, string mode = "explicit"){
    PDE currentPDE;
    if (mode == "explicit"){
      currentPDE.writeCycle0(titGroupList, ini);
      currentPDE.writeCycle1(fileName, titGroupList, ini);
    } else {
      currentPDE.writeEnergy(fileName, ini);
    }
  }
  
  /** @brief Writes out PKINT files to use with Monte-Carlo method, Karlsberg2.
   *         
   *  Depending on the current cycle a PKINT file is written out. Before, all
   *  pkint values where computed using POTAT files. A PKINT file contains all
   *  intrinsic pKA values for every titratable group from reference to its other
   *  states.
   *  
   *  @param string Current cycle name.
   *  @param string Prefix of file name to write out PKINT file.
   *  @param INI Config INI file.
   *  @return Void.
   */    
  void writeOutPkint(string cycleName, string fileName, INI &ini){
    string eps_in = ini.get<string>("experiment.eps_in");
    string eps_out = ini.get<string>("experiment.eps_out");
    bool eps_same = false;
    if (eps_in == eps_out)
      eps_same = true;
    
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
	float cycle0Shift = getCycle0Shift(currentTitGroup.getResidueName(), stateNr+1);
	if ( cycle0Shift != NOT_IN_ST && !explicitModels){
	  if (cycleName == "cycle0"){
	    cout << "using precalculated values." << endl;
	    cycle0Shift *= CONVERT;
	    born = 0; //currentTitGroup.getBorn(cycleNr, stateNr);
	    back = 0; //currentTitGroup.getBack(cycleNr, stateNr);
	    shift = getShift(currentTitGroup.getResidueName(), stateNr)*CONVERT;
	    if (eps_same)
	      result = 0;
	    else
	      result = cycle0Shift + ( born - back );
	    
	    cout << currentTitGroup.getIdentifier() << ": R -> " << getState(currentTitGroup.getResidueName(), stateNr) << endl;
	    cout << "Gborn: " << born << endl;
	    cout << "Gback: " << back << endl;
	    cout << "intrinsic pka = "<< shift << " + " << (cycle0Shift) << " + ( " << born << " - " << back << ") = " << result << " [kJ/mol]" << endl;
	    cout << "-> " << result/CONVERT << " []" << endl;
	    
	  } else if (cycleName == "cycle1") {
	    born = currentTitGroup.getBorn(cycleNr, stateNr);
	    back = currentTitGroup.getBack(cycleNr, stateNr);
	    shift = 0;
	    result = shift + ( born - back );
	    
	    cout << currentTitGroup.getIdentifier() << ": R -> " << getState(currentTitGroup.getResidueName(), stateNr) << endl;
	    cout << "Gborn: " << born << endl;
	    cout << "Gback: " << back << endl;
	    cout << "intrinsic pka = "<< shift << " + ( " << born << " - " << back << ") = " << result << " [kJ/mol]" << endl;
	    cout << "-> " << result/CONVERT << " []" << endl;
	    
	  } else {
	    cout << "using precalculated values." << endl;
	    cycle0Shift *= CONVERT;
	    born = currentTitGroup.getBorn(1, stateNr) - 0; //currentTitGroup.getBorn(0, stateNr);
	    back = currentTitGroup.getBack(1, stateNr) - 0; //currentTitGroup.getBack(0, stateNr);
	    shift = getShift(currentTitGroup.getResidueName(), stateNr)*CONVERT;
	    if (eps_same)
	      result = 0;
	    else 
	      result = shift+(born - back)-(cycle0Shift);
	    
	    cout << currentTitGroup.getIdentifier() << ": R -> " << getState(currentTitGroup.getResidueName(), stateNr) << endl;
	    cout << "Gborn: " << born << endl;
	    cout << "Gback: " << back << endl;
	    cout << "intrinsic pka = "<< shift << " + ( " << born << " - " << back << ") - " << (cycle0Shift) << " = " << result << " [kJ/mol]" << endl;
	    cout << "-> " << result/CONVERT << " []" << endl;
	  }
	  
	  // write pkint [kJ/mol] and transition identifier to pkint file
	  out << result << " " << getState(currentTitGroup.getResidueName(), stateNr) << " ";
	  
	} else {
	  
	  if (cycleName == "cycle0"){
	    born = currentTitGroup.getBorn(cycleNr, stateNr);
	    back = currentTitGroup.getBack(cycleNr, stateNr);
	    shift = 0;
	    result = shift + ( born - back );
	  } else if (cycleName == "cycle1") {
	    born = currentTitGroup.getBorn(cycleNr, stateNr);
	    back = currentTitGroup.getBack(cycleNr, stateNr);
	    shift = 0;
	    result = shift + ( born - back );
	  } else {
	    born = currentTitGroup.getBorn(1, stateNr) - currentTitGroup.getBorn(0, stateNr);
	    back = currentTitGroup.getBack(1, stateNr) - currentTitGroup.getBack(0, stateNr);
	    shift = getShift(currentTitGroup.getResidueName(), stateNr);
	    shift *= CONVERT;
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
	
      }
      
      out << currentTitGroup.getIdentifier() << endl;
      
    }
    
  }

  /** @brief Every state has a char identifier which is returned from ST file.
   *         
   *  Returns a char character which is stored in ST file for every titratable group
   *  and every state.
   *  
   *  @param string Name of the residue.
   *  @param int State of the residue to get it char identifier.
   *  @return char Char identifier of residue in given state.
   */    
  char getState(string residueName, int stateNr){
    for (unsigned int i = 0; i < stList.size(); i++){
      ST currentST = stList.at(i);
      if (currentST.getTitGroupName() == residueName && currentST.getStateNr() == (stateNr+1)){
	return currentST.getState();
      }
    }
    return '?';
  }
  
  /** @brief Returns the pKA shift from reference to given state.
   *         
   *  Returns a pKA shift value from experimental values (defined in ST file) 
   *  which is defined from reference residue to same residue to given state.
   *  
   *  @param string Name of the residue
   *  @param int State of the residue to get it char identifier.
   *  @return float Returns the pKA shift from reference to given state.
   */    
  float getShift(string residueName, int stateNr){
    for (unsigned int i = 0; i < stList.size(); i++){
      ST currentST = stList.at(i);
      if (currentST.getTitGroupName() == residueName && currentST.getStateNr() == (stateNr+1)){
	return currentST.getShift();
      }
    }
    return 0;
  }
  
  /** @brief Returns the averaged pKA value to skip second thermodynamic cycle.
   *         
   *  Returns an averaged pKA value which mimics the electrostatics part of the 
   *  energy shift in the thermodynamic cycle from reference to given state of 
   *  every titratable group/residue. Values are given from Sakalli et al.
   *  
   *  @param string Name of the residue
   *  @param int State of the residue to get it char identifier.
   *  @return float Returns the average pKA value for a given residue in 
   *                a given state
   *
   */   
  float getShiftCycle0(string residueName, int stateNr){
    for (unsigned int i = 0; i < stList.size(); i++){
      ST currentST = stList.at(i);
      if (currentST.getTitGroupName() == residueName && currentST.getStateNr() == (stateNr+1)){
	return currentST.getShiftCycle0();
      }
    }
    return 0;
  }
  
  /// Are ST files already parsed.
  static bool STparsed;
  /// Do we have to compute the N-terminus
  static bool calcNTE;
  /// Do we have to compute the C-terminus
  static bool calcCTE;
  /// Are we using the two-cycle approach (compute models explicitly)
  static bool explicitModels;
  
  
 private:
  /// File name of the PQR file which has to be parsed (protein).
  string fileName;
  /// Folder of ST files, which define charge pattern for different states
  /// of every titratable group in CHARMM format.
  string stFilesFolder;
  /// Placeholder for molecular volume mesh computation.
  Model molecule;
  /// Atom list of all atoms with all PQR information.
  vector<Atom> atomList;
  /// A list of all residues in the protein.
  vector<Residue> residueList;
  /// A list of all titratable groups in the protein (definition via framework of
  /// pKA computations from Knapp et al.).
  vector<Residue> titGroupList;
  /// A list of all states in every ST file for every titratable group.
  vector<ST> stList;
  /// An abstraction of the thermodynamic cycle
  tCycle cycle;
  /// A multidensional W-matrix vector to store pairwise interaction energies of all titratable groups with each other.
  vector< vector< vector< vector<double> > > > wmatrix;
};

bool PQR::STparsed = false;
bool PQR::explicitModels = false;
bool PQR::calcCTE = false;
bool PQR::calcNTE = false;
