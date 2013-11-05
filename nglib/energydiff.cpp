// implement an energy difference between two calculations
 
#include <solve.hpp>
using namespace ngsolve;
#include <vector>
#include <limits>

typedef std::numeric_limits< double > dbl;

class NumProcEnergyCalc : public NumProc
{
  protected:
	///
    BilinearForm * bfa;
    ///
    LinearForm * lff;
    ///
    GridFunction * gfu;
    ///
    GridFunction * gfu0;
    ///
    GridFunction * gfv;
    ///
    Vector<double> point;
    ///
    Array<int> domains;
    ///
    bool integrateonplanes;
    ///
    bool usepoint3and4;
    ///
    int variabledirection;
    ///
    int n[3];
    ///
    string filename, text;
    ///
    string variablename;
    ///
    bool applyd;
    ///
    bool hermitsch;
    ///
    int component;
    ///
    int outputprecision;
    
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
    ///
    NumProcEnergyCalc (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcEnergyCalc();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcEnergyCalc (pde, flags);
    }
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "energydiff";
    }
    
    void readMolecule();

  };


NumProcEnergyCalc :: NumProcEnergyCalc (PDE & apde, const Flags & flags)
: NumProc (apde), point(1)
{
	bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""), 1); 
    lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", ""), 1);
    gfu  = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""), 0);
    gfu0 = pde.GetGridFunction (flags.GetStringFlag ("gridfunction0", ""), 0);
    gfv = pde.GetGridFunction (flags.GetStringFlag ("gridfunction2", ""), 1); 
    pqrfile = flags.GetStringFlag ("pqrfile","pqr");
      
    readMolecule();

    variablename = flags.GetStringFlag ("resultvariable", "");
    
    point.SetSize(3);

    if (flags.NumListFlagDefined ("point"))
    {
    	const Array<double> & p = flags.GetNumListFlag ("point");
    	point.SetSize(p.Size());
    	for (int i = 0; i < p.Size(); i++)
    		point(i) = p[i];
          
    	cout << "point = " << point << endl;
    }

      
    integrateonplanes = flags.GetDefineFlag("integrateonplanes");

    variabledirection = static_cast<int>(flags.GetNumFlag("variabledirection",0))-1;
    
    n[0] = static_cast<int>(flags.GetNumFlag("n1",0));
    n[1] = static_cast<int>(flags.GetNumFlag("n2",0));
    n[2] = static_cast<int>(flags.GetNumFlag("n3",0));

    text = flags.GetStringFlag ("text","energydiff");

    if(flags.StringFlagDefined("filename"))
      filename = pde.GetDirectory() + dirslash + flags.GetStringFlag("filename","");
    else
      filename = "err.out";

    applyd = flags.GetDefineFlag ("applyd");
    hermitsch = flags.GetDefineFlag ("hermitsch");

    outputprecision = (pde.ConstantUsed("outputprecision")) ? int(pde.GetConstant("outputprecision")) : -1;
    if(flags.NumFlagDefined("outputprecision"))
      outputprecision = int(flags.GetNumFlag("outputprecision",-1));

    component = static_cast<int>(flags.GetNumFlag("cachecomp",1))-1;
}

NumProcEnergyCalc :: ~NumProcEnergyCalc()
{
    ;
}
  
void NumProcEnergyCalc::readMolecule()
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
	//	cout << i << ": " << t.fieldName << " " << t.x << " " << t.y  << " " << t.z  << " " << t.charge  << " " << t.radius <<  endl;
	}    
}



void NumProcEnergyCalc :: PrintDoc (ostream & ost)
{
	ost << "(documentation not updated!)";
}

void NumProcEnergyCalc :: Do(LocalHeap & lh)
{
	  
	double result = 0;
	double potential = 0;
	double conversion = 1.602176565e-19 * 6.02214129e23 * 0.5 / 1000.0;

	const BilinearFormIntegrator & bfi = (bfa) ? *bfa->GetIntegrator(0) : *gfu->GetFESpace().GetIntegrator();
	const BilinearFormIntegrator & bfi2 =  (bfa) ? *bfa->GetIntegrator(0) : *gfu->GetFESpace().GetIntegrator();
	
	cout.precision(dbl::digits10);

	for(unsigned int kk=0; kk<molecule.size(); kk++)	
	{
		Atom currentAtom = molecule.at(kk);
		point(0) = currentAtom.x;
		point(1) = currentAtom.y;
		point(2) = currentAtom.z;
		
		FlatVector<double> pflux(bfi.DimFlux(), lh);
		CalcPointFlux (ma, *gfu, point, domains,
				pflux, bfi, applyd, lh, component);
	      
		result = pflux(0);
		cout << "(" << pflux(0)*conversion;
	      
		CalcPointFlux (ma, *gfu0, point, domains,
				pflux, bfi2, applyd, lh, component);

		result -= pflux(0); 
		cout << " - " << pflux(0)*conversion << ")*" << currentAtom.charge;

	      
		potential += currentAtom.charge * result;

		cout << " = " << currentAtom.charge * result*conversion << endl;

	}

	double energy;
	energy = (1.602176565e-19 * 6.02214129e23 * 0.5 * potential)/1000;

	if (MyMPI_GetNTasks() == 1 || MyMPI_GetId() != 0){
		ofstream resultFile;
		resultFile.open("result.out");
		resultFile << "The energy difference is " << energy << " [kJ/mol].\n";
		resultFile.close();
		cout << "The energy difference is " << setprecision(12) << energy << " [kJ/mol].\n";
	}
        
	pde.GetVariable(variablename,true) = result;
}

void NumProcEnergyCalc :: PrintReport (ostream & ost)
{
	ost << "NumProcEnergyCalc:" << endl;
}

static RegisterNumProc<NumProcEnergyCalc> init_np_pts("energydiff");
