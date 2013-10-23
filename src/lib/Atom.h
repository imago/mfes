#include <string>
#include "Point3.h"

using namespace std;

class Atom {
public:
	Atom(int _atomNumber, string _atomName, char _confID, string _residueName, char _chainID, int _residueNumber, Point3 _coord, float _charge, float _radius, string _segName):
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
		cout << "Coordinates : " << coord << endl;
		cout << "Charge      : " << charge << endl;
		cout << "Radius      : " << radius << endl;
		cout << "Segment name: " << segName << endl;


	}
protected:
   string atomName;
   float charge;

private:
   string fieldName;
   char confID;
   int atomNumber;
   string residueName;
   char chainID;
   int residueNumber;
   Point3 coord;
   float radius;
   string segName;

};
