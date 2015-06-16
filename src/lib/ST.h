/** @file ST.h
 *  @brief ST file parser and titratable group builder
 *
 *  Using this file, titratable groups are build as residues, have a state
 *  and their state may be changed depending on the framework to compute
 *  pKA values.
 *
 *  @author Ilkay Sakalli
 */

#include "Charge.h"
#include <vector>
#include <string>

using namespace std;
/**
 * @class ST
 *
 * @brief This class is a datastructure to represent ST files. These files hold
 *        all information about titratable groups in the protein.
 *
 * Titratable groups are generated and have a state. They can ge shifted to another 
 * state depending on the framework of pKA compuations. Every titratable group type
 * may have an average energy shift value for one-cycle approach computations.
 *
 * @author Ilkay Sakalli
 *
 */ 


class ST {
public:

	ST(string _titGroupName, char _state, int _stateNr, float _shift, vector<Charge>_rules, float _shiftCycle0):
    titGroupName(_titGroupName),
    state(_state),
    stateNr(_stateNr),
    shift(_shift),
    rules(_rules),
	shiftCycle0(_shiftCycle0){
	}

	string getTitGroupName(){
	    return titGroupName;
	}

	char getState(){
	   return state;
	}

	int getStateNr(){
	   return stateNr;
	}

	float getShift(){
		return shift;
	}

	float getShiftCycle0(){
		return shiftCycle0;
	}

	vector<Charge> getRules(){
	    return rules;
	}

	void print(){
		cout << "residue name: " << titGroupName << endl;
		cout << "state:        " << state << endl;
		cout << "state nr:     " << stateNr << endl;
		cout << "shift:        " << shift << endl;
		cout << "cycle0 shift: " << shiftCycle0 << endl;
	}

private:
  string titGroupName;
  char state;
  int stateNr;
  float shift;
  vector<Charge> rules;
  float shiftCycle0;
};
