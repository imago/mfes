#ifndef ATOM_H
#define ATOM_H

#include <string>
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

	char getChainID(){
		return chainID;
	}


private:
   string fieldName;
   char confID;
   int atomNumber;
   string residueName;
   char chainID;
   int residueNumber;
   Point3<float> coord;
   float radius;
   string segName;
   string atomName;
   float charge;

};

#endif
