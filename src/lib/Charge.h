#ifndef CHARGE_H
#define CHARGE_H

#include <string>
#include <vcg/space/point3.h>

using namespace std;

class Charge  {
public:
	Charge(string _atomName, float _charge):
	    atomName(_atomName),
	    charge(_charge) {
	}
	void print(){
		cout << "Atom name   : " << atomName << endl;
		cout << "Charge      : " << charge << endl;
	}

	string getAtomName(){
		return atomName;
	}

	float getCharge(){
		return charge;
	}

private:
	string atomName;
	float charge;

};

#endif
