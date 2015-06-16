/** @file Residue.h
 *  @brief Holds the general description of a residue in proteins.
 *
 *  A datastructure for residues is built, knowing if they are titratable groups,
 *  their atom lists and other information.
 *
 *  @author Ilkay Sakalli
 */

#ifndef RESIDUE_H
#define RESIDUE_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <string>
#include "Atom.h"
#include "tCycle.h"

using namespace std;

/**
 * @class Residue
 *
 * @brief This class is a datastructure to represent residues in proteins. 
 *
 * Residues in proteins may be titratable groups and have different states, 
 * have a residue name, chain ID, an atom list, are assigned to one of both
 * titratable group approaches to compute pKA values and hold born or back
 * energy terms.
 *
 * @author Ilkay Sakalli
 *
 */ 

class Residue {
public:
	Residue():
          resName(""),
          resNumber(-1),
 	  hasCycle(false)
	{
		born[0][0] = born[0][1] = born[0][2] = born[1][0] = born[1][1] = born[1][2] = 0.0f;
		back[0][0] = back[0][1] = back[0][2] = back[1][0] = back[1][1] = back[1][2] = 0.0f;

	}

	unsigned int getNrStates(){
	  map<string, int> nrStates  = { {"ARG", 2}, {"CTE", 2}, {"DPP", 3}, {"EPP", 3}, {"HSP", 3}, {"LYS", 2}, {"NTE", 2}, {"TYR", 2}, {"CYS", 2}};
	  return nrStates.at(resName);
	}

	bool isTitGroup(){
		set<string> validGroups = {"ARG","DPP", "HSP", "NTE", "CTE", "EPP", "LYS", "TYR", "CYS"};
		if(validGroups.find(resName.c_str())!=validGroups.end()){
	      return true;
		}
		return false;
	}

	void setResidueName(string _resName){
	  resName = _resName;
	}

	string getResidueName(){
		return resName;
	}

	void setResidueNumber(int _resNumber){
	    resNumber = _resNumber;
	}

	int getResidueNumber(){
		return resNumber;
	}

	void setAtomList(vector<Atom> &_atomList){
	  atomList = _atomList;
	}

	vector<Atom>& getAtomList(){
		return atomList;
	}

	void setChainID(char _chainID){
	    chainID = _chainID;
	}

	string getIdentifier(){
		stringstream result;
		result << resName << "-" << resNumber;
		if (chainID != ' ')
			result << "_" << chainID;
		return result.str();
	}

	void setTCycle(tCycle &_cycle){
		cycle = _cycle;
		hasCycle = true;
	}
	
	bool haveCycle(){
	  return hasCycle;
	}
	
	tCycle getTCycle(){
	    return cycle;
	}

	void setBorn(float bornEner, unsigned int cycleNr, unsigned int stateNr){
		born[cycleNr][stateNr] = bornEner;
	}

	float getBorn(unsigned int cycleNr, unsigned int stateNr){
		return born[cycleNr][stateNr];
	}

	void setBack(float backEner, unsigned int cycleNr, unsigned int stateNr){
		back[cycleNr][stateNr] = backEner;
	}

	float getBack(unsigned int cycleNr, unsigned int stateNr){
		return back[cycleNr][stateNr];
	}

	void print(){
	  for (unsigned int i = 0; i < atomList.size(); i++){
	    Atom currentAtom = atomList.at(i);
	    cout << currentAtom.pqrLine() << endl;
	  }
	}

private:
	string resName;
	int resNumber;
	char chainID;
	vector<Atom> atomList;
	tCycle cycle;
	bool hasCycle;
	float born[2][3];
	float back[2][3];
};

#endif
