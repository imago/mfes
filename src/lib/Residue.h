#ifndef RESIDUE_H
#define RESIDUE_H

#include <string>
#include <vector>
#include <set>
#include "Atom.h"
#include "tCycle.h"

using namespace std;

class Residue {
public:
	Residue():
	    resName(""),
	    resNumber(-1)
	{}

	bool isTitGroup(){
		string groups[] = {"ARG","DPP", "HSP", "NTE", "CTE", "EPP", "LYS", "TYR"};
		set<string> validGroups(groups, groups + 8);
		if(validGroups.find(resName.c_str())!=validGroups.end()){
	      return true;
		}
		return false;
	}

	void setResidueName(string _resName){
	  resName = _resName;
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

	void setChainID(char _chainID){
	    chainID = _chainID;
	}

private:
	string resName;
	int resNumber;
	char chainID;
	vector<Atom> atomList;
	tCycle cycle;
	float back;
	float born;
};

#endif
