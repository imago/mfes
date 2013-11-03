#include "Charge.h"
#include <vector>
#include <string>

using namespace std;

class ST {
public:

	ST(string _titGroupName, char _state, int _stateNr, float _shift, vector<Charge>_rules):
    titGroupName(_titGroupName),
    state(_state),
    stateNr(_stateNr),
    shift(_shift),
    rules(_rules){

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

	vector<Charge> getRules(){
	    return rules;
	}

private:
  string titGroupName;
  char state;
  int stateNr;
  float shift;
  vector<Charge> rules;
};
