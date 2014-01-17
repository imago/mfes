// precalculate atom position in netgen volume

#include <solve.hpp>
#include <vector>
using namespace ngsolve;


class NumProcPreCalculation : public NumProc
{
protected:
  string pqrfile;

  class Atom
   {
   public:
     string fieldName;
     int atomNo;
     string atomName;
     string residueName;
     string chainName;
     int residueNo;
     float x;
     float y;
     float z;
     float charge;
     float radius;
   };
   
   vector<Atom> molecule;
  
public:
    
  NumProcPreCalculation (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
      pqrfile = flags.GetStringFlag ("pqrfile","");      
      readMolecule();
  }

  virtual ~NumProcPreCalculation() 
  { ; }
  
  void readMolecule()
  {
	  string line;
	  ifstream in(pqrfile.c_str());

	  if (!in) {
		  in.close();
		  cout << "Cannot open pqr file:" << pqrfile << endl;
		  exit(0);
	  }

	  int i = 0;
	  while( !in.eof() ) {
		getline(in, line);
		if ( line.find("ATOM", 0)==string::npos ) continue;

		Atom t;
		i++;

		// atom name
		t.atomNo = atoi(line.substr(5,6).c_str());

		// residue name
		t.residueName = line.substr(17,3).c_str();

		// residue number
		t.residueNo = atoi(line.substr(22,4).c_str());

		// coordinates
		t.x = atof(line.substr(30,8).c_str());
		t.y = atof(line.substr(38,8).c_str());
		t.z = atof(line.substr(46,8).c_str());

		// charge
		t.charge = atof(line.substr(54,6).c_str());
			   
		// radius
		t.radius = atof(line.substr(60,6).c_str());

		molecule.push_back(t);
		//		cout << t.atomNo << ": " << t.fieldName << " " << t.x << " " << t.y  << " " << t.z  << " " << t.charge  << " " << t.radius <<  endl;
	}
  }
  


  virtual void Do(LocalHeap & lh)
  {
    cout << IM(1) << "Preprocessing calculations" << endl;
 
    int elnr;
    Array<int> dnums;
    Vec<3> point;


    if (MyMPI_GetNTasks() == 1 || MyMPI_GetId() != 0)
      for(unsigned int i=0; i<molecule.size(); i++)
	{
	  Atom currentAtom = molecule.at(i);
	  point(0) = currentAtom.x;
	  point(1) = currentAtom.y;
	  point(2) = currentAtom.z;

	  IntegrationPoint PtOnRef(0,0,0,1);

	  elnr = ma.FindElementOfPoint(point,PtOnRef,true);
	  double p[3] = { PtOnRef(0), PtOnRef(1), PtOnRef(2) };
	  ma.addElIndex(elnr);
	  ma.addIpPoint(p);
	}
        cout << ma.getElIndexSize() << " atoms preprocessed." << endl;
  }

  virtual string GetClassName () const
  {
    return "precalculation";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl;
  }

  static void PrintDoc (ostream & ost)
  {
    ost << 
      "Precalculat elnr for atom position coordinates.\n" 
	<< endl;
  }
};



static RegisterNumProc<NumProcPreCalculation> init_np_pts("precalculation");


