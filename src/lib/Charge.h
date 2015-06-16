/** @file Charge.h
 *  @brief Defining a point charge in the implicit solvent model.
 *
 *  A file holds the definition of a point charge using the location 
 * of the atom point charge. A charge has a name and a charge in Coulomb.
 *
 *  @author Ilkay Sakalli
 */

#ifndef CHARGE_H
#define CHARGE_H

#include <string>
#include <vcg/space/point3.h>

using namespace std;

/**
 * @class Charge
 *
 * @brief An atom point charge is defined.
 *
 * This small class uses a three point class Point3f from VCGlib
 * and defines a face on the molecular surface.
 *
 * @author Ilkay Sakalli
 *
 */ 


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
