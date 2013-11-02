#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <iomanip>
#include <vcg/space/point3.h>

using namespace vcg;
using namespace std;

class Atom {
public:
	Atom(){};

	Atom(int _atomNumber, string _atomName, char _confID, string _residueName, char _chainID, int _residueNumber, Point3<float> _coord, float _charge, float _radius, string _segName):
		    atomNumber(_atomNumber),
		    atomName(_atomName),
		    confID(_confID),
		    residueName(_residueName),
		    chainID(_chainID),
		    residueNumber(_residueNumber),
		    coord(_coord),
		    charge(_charge),
		    radius(_radius),
		    segName(_segName)
		{

		}

	void print(){
		cout << "Atom number : " << atomNumber << endl;
		cout << "Atom name   : " << atomName << endl;
		cout << "Conform. nr : " << confID << endl;
		cout << "Residue name: " << residueName << endl;
		cout << "Chain nr    : " << chainID << endl;
		cout << "Residue nr  : " << residueNumber << endl;
		cout << "Coordinates : (" << coord.X() << ", " << coord.Y() << ", " << coord.Z() << ")" << endl;
		cout << "Charge      : " << charge << endl;
		cout << "Radius      : " << radius << endl;
		cout << "Segment name: " << segName << endl;
	}

	Point3<float> getCoord(){
		return coord;
	}

	float getRadius(){
		return radius;
	}

	int getResidueNumber(){
		return residueNumber;
	}

	string getResidueName(){
		return residueName;
	}

	void setResidueName(string _residueName){
		residueName = _residueName;
	}

	char getChainID(){
		return chainID;
	}

	string getAtomName(){
		return atomName;
	}

	void setCharge(float _charge){
		charge = _charge;
	}

	float getCharge(){
		return charge;
	}

	int getAtomNumber(){
		return atomNumber;
	}

	string pqrLine(){
		stringstream ss;
		ss << "ATOM" << setw(7) << atomNumber << setw(5) << atomName << setw(1)
				<< confID << setw(3) << residueName << setw(2) << chainID << setw(4)
				<< residueNumber << setw(12) << coord.X() << setw(8) << coord.Y()
				<< setw(8) << coord.Z() << setw(6) << charge << setw(6)
				<< radius << setw(7) << segName;
		return ss.str();
	}

private:
	int atomNumber;
	string atomName;
	char confID;
	string residueName;
	char chainID;
	int residueNumber;
	Point3<float> coord;
	float charge;
	float radius;
	string segName;
};

#endif
