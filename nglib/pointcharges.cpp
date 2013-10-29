// introduce point sources in the linearform


#include <solve.hpp>
#include <vector>
using namespace ngsolve;


const double q0 = 1.60217646e-19;


class NumProcPointCharges : public NumProc
{
protected:
  LinearForm* lff;
  GridFunction * gfu;
  string pqrfile;
  bool interpolate;
  
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
    
  /*
    In the constructor, the solver class gets the flags from the pde - input file.
    the PDE class apde constains all bilinear-forms, etc...
  */
  NumProcPointCharges (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
      lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", "f"));
      pqrfile = flags.GetStringFlag ("pqrfile","");
      interpolate = flags.GetDefineFlag ("interpolate");
      
      readMolecule();
  }

  virtual ~NumProcPointCharges() 
  { ; }
  
  
  void reset()
  {
	  Vec<3> p;
	 
	  for (int j = 0; j < ma.GetNP(); j++)
	  {
		  lff->GetVector().FVDouble() (j)  =  0;
	  }

  }
  
  void readMolecule()
  {
	  string line;
	  ifstream in(pqrfile.c_str());
	  int i = 0;
	  while( !in.eof() ) {
		getline(in, line);
		if ( line.find("ATOM", 0)==string::npos ) continue;

		Atom t;
		i++;

		// atom name
		t.atomNo = atoi(line.substr(12,4).c_str());

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
		cout << i << ": " << t.fieldName << " " << t.x << " " << t.y  << " " << t.z  << " " << t.charge  << " " << t.radius <<  endl;
	}
  }
  


  virtual void Do(LocalHeap & lh)
  {
    cout << IM(1) << "Add point sources" << endl;
 
    // Decide if interpolate or add whole charge to a point

    // Set all values to zero
     if (MyMPI_GetNTasks() == 1 || MyMPI_GetId() != 0)
      reset();

    if (!interpolate){
    	Array<int> dnums;
		Vec<3> point;
		double charge;
		Vec<3> p;
		for(unsigned int i=0; i<molecule.size(); i++)
		{
			Atom currentAtom = molecule.at(i);
		  point(0) = currentAtom.x;
		  point(1) = currentAtom.y;
		  point(2) = currentAtom.z;
		  charge = currentAtom.charge;
		  if (charge != 0) {
			  double mindist = 1e99;
			  int minpi = -1;

			  for (int j = 0; j < ma.GetNP(); j++)
			  {
				  p = ma.GetPoint<3> (j);
				  //if (L2Norm (p-point) < mindist && L2Norm (p-point) != 0)
				  if (L2Norm (p-point) < mindist )
				  {
					  mindist = L2Norm (p-point);
					  minpi = j;
					 // cout << "Debug: p: " << p << ", point: " << point << ", i: " << i << ", mindist: " << mindist << ", minpi: " << minpi << ", L2Norm: " << L2Norm (p-point) << endl;
				  }
			  }
			  p = ma.GetPoint<3> (minpi);
			  cout << i << ": closest point to (" << point(0) << ", " << point(1) << ", " << point(2) << ") is point nr " << minpi << "(" << p(0) << ", " << p(1) << ", " << p(2) << "), dist = " << mindist <<", adding charge: " << charge << endl;

			  lff->GetVector().FVDouble() (minpi)  +=  charge * q0;
 
		  }
		}
    } else {
      cout << IM(1) << "Interpolating point charges" << endl;
		int elnr;
		Array<int> dnums;
		Vec<3> point;
		double charge;

		if (MyMPI_GetNTasks() == 1 || MyMPI_GetId() != 0)
		for(unsigned int i=0; i<molecule.size(); i++)
		{
			Atom currentAtom = molecule.at(i);
		 //  cout << "atom " << i << endl;
			point(0) = currentAtom.x;
			point(1) = currentAtom.y;
			point(2) = currentAtom.z;
			charge = currentAtom.charge;

			if (charge == 0){
			//cout << "skipping atom id "<< i << endl;
			  continue;
			}

			// cout << "adding charge (" << charge << ") for position (" << point(0) <<", " << point(1) <<", "<< point(2) << ")" << endl;
			IntegrationPoint PtOnRef(0,0,0,1);

			// cout << "point = " << point << endl;
			elnr = ma.FindElementOfPoint(point,PtOnRef,true);
			// cout << "id = " << MyMPI_GetId() << "elnr = " << elnr << endl;

			if (elnr == -1) continue;
			const FESpace & fes = lff->GetFESpace();
			fes.GetDofNrs (elnr, dnums);

			const ScalarFiniteElement<3> & fel = dynamic_cast<const ScalarFiniteElement<3>&> (fes.GetFE (elnr, lh));

			FlatVector<> shape(dnums.Size(), lh);

			fel.CalcShape(PtOnRef,shape);

			for(int i=0; i<dnums.Size(); i++)
			{
			  FlatVector<double> vector = lff->GetVector().FV<double>();
			  vector (dnums[i])  +=  charge * shape(i) * q0;
			}

		}
    }
  }

  virtual string GetClassName () const
  {
    return "pointcharges";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
	<< "Linear-form     = " << lff->GetName() << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << 
      "Insert point source\n" 
	<< endl;
  }
};



static RegisterNumProc<NumProcPointCharges> init_np_pts("pointcharges");


