// writing out results to formatted potential files
 
#include <solve.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <list>

using namespace ngsolve;

 
///
class NumProcWritePotat : public NumProc
{
protected:
	///
    BilinearForm * bfa;
    ///
    LinearForm * lff;
    ///
    GridFunction * gfu;
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
    ///
    string pqrfile;
    ///
    string potatfile;
    ///
    int states;
    ///
    int noOfAtoms;
    ///
    string createNewFile;
    ///
    bool initsearchtree;
    ///
    bool predef;

    class Atom
     {
     public:
       string fieldName;
       int atomNo; // has to begin with 1!
       string atomName;
       string residueName;
       string chainName;
       int residueNo;
       float x;
       float y;
       float z;
       float charge;
       float radius;
       int operator<(const Atom &rhs) const
       {
    	   if( this->x < rhs.x){ return 1;
    	   } else if ( this->x == rhs.x){
    		   if (this->y < rhs.y){
    			   return 1;
    		   } else if (this->y == rhs.y){
    			   if (this->z < rhs.z){
    				   return 1;
    			   }
    		   }
    	   }
    	   return 0;
       }
    };
    
    list<Atom> molecule;
  
public:
    NumProcWritePotat (PDE & apde, const Flags & flags)
    	: NumProc (apde), point(1)
    {
    	bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""), 1);
    	lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", ""), 1);
    	gfu  = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""), 0);
    	gfv = pde.GetGridFunction (flags.GetStringFlag ("gridfunction2", ""), 1);
    	pqrfile = flags.GetStringFlag ("pqrfile","pqr");
    	potatfile = flags.GetStringFlag ("potatfile","potat");
	initsearchtree = flags.GetDefineFlag("initsearchtree");
    	states = flags.GetNumFlag ("statenr", 1);
	predef = flags.GetDefineFlag("predef");

    	createNewFile = flags.GetStringFlag ("file","");


    	noOfAtoms = readMolecule();

    	variablename = flags.GetStringFlag ("resultvariable", "");

    	point.SetSize(3);

    	if (flags.NumListFlagDefined ("point"))
    	{
    		const Array<double> & p = flags.GetNumListFlag ("point");
    		point.SetSize(p.Size());
    		for (int i = 0; i < p.Size(); i++)
    			point(i) = p[i];

//    		cout << "point = " << point << endl;
    	}


    	integrateonplanes = flags.GetDefineFlag("integrateonplanes");
    	variabledirection = static_cast<int>(flags.GetNumFlag("variabledirection",0))-1;

    	n[0] = static_cast<int>(flags.GetNumFlag("n1",0));
    	n[1] = static_cast<int>(flags.GetNumFlag("n2",0));
    	n[2] = static_cast<int>(flags.GetNumFlag("n3",0));

    	text = flags.GetStringFlag ("text","WritePotat");

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

    ~NumProcWritePotat()
    {
      ;
    }

    
  
    int readMolecule()
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
		//	        cout << t.atomNo << ": " << t.fieldName << " " << t.x << " " << t.y  << " " << t.z  << " " << t.charge  << " " << t.radius <<  endl;
    	}
    	return i;
    }

  void Do(LocalHeap & lh)
  {

    double conversion = 1.602176565e-19 * 6.02214129e23 / 1000.0;
    
    stringstream str(stringstream::out|stringstream::binary);
    
    if (MyMPI_GetNTasks() == 1 || MyMPI_GetId() != 0) {
      const BilinearFormIntegrator & bfi = (bfa) ? *bfa->GetIntegrator(0) : *gfu->GetFESpace().GetIntegrator();
     
      unsigned int elnr, idx;
	      
      list<Atom>::iterator it;
      int position = 1;
      molecule.sort();

      bool firstone(true);

      for (it=molecule.begin(); it!=molecule.end(); ++it)
	{
	  Atom line;
	  line = *it;
	  point(0) = line.x;
	  point(1) = line.y;
	  point(2) = line.z;

	  double* p = new double[3];
	  if(firstone)
	    {
	      // normally done in adding pointcharges
	      // just generated if no searchtree generated before
	      IntegrationPoint PtOnRef(0,0,0,1);
       	      elnr = ma.FindElementOfPoint(point,PtOnRef,true);
	      firstone = false;
	    }

	  FlatVector<double> pflux(bfi.DimFlux(), lh);

	  // predef = false;

	  if (!predef){
	    CalcPointFlux (ma, *gfu, point, domains, pflux, bfi, false, lh, component);
	  }
	  else {
	    idx = line.atomNo-1;
	    elnr     = ma.getElIndex(idx);
	    ma.setIpPoint(idx, p);

	    IntegrationPoint ip (p[0], p[1], p[2], 1);
	    
	    CalcElFlux (ma, *gfu, point, domains, pflux, bfi, false, lh, elnr, ip, component);
	    
	  }
 
	  str << point(0) << endl;
	  str << point(1) << endl;
	  str << point(2) << endl;
	  str << conversion * pflux(0) << endl;

	  position++;
	  
	}
	
    }
  
    
    // Writing out results
    cout << IM(1) << "looking for " << potatfile.c_str() << endl;

    ofstream potfile;
   
    struct stat st;
    if(stat(potatfile.c_str(),&st) != 0){
      cout << IM(1)  << "file not found. creating and writing states" << endl;
      potfile.open( potatfile.c_str());
    } else if (createNewFile == "create") {
      cout << IM(1)  << "file found but creating new one and writing states" << endl;
      potfile.open( potatfile.c_str());
    } else {
      potfile.open( potatfile.c_str(), ios::app );
    }
    if (states != -1)
      potfile << states << endl;
      
    cout << IM(1)  << "Writing noOfAtoms for next state: " << noOfAtoms << endl;
    potfile << noOfAtoms << endl;
    potfile << str.str();
    potfile.close();

  }

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcWritePotat (pde, flags);
    }

    virtual string GetClassName () const
    {
      return "writepotat";
    }

    virtual void PrintReport (ostream & ost)
    {
    	ost << GetClassName() << endl
    		  << "Linear-form     = " << lff->GetName() << endl;
    }

    static void PrintDoc (ostream & ost)
    {
      ost <<
        "Write out structured potential file\n"
  	<< endl;
    }


};


static RegisterNumProc<NumProcWritePotat> init_np_pts("writepotat");
