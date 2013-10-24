 // implement an energy difference between two calculations
 
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
    
    string pqrfile;

    string potatfile;

    int states;

    int nrOfValues;

    // TODO Should be removed 
    int noOfAtoms;
    
    
    class pqrLine
  {
  public:
    string fieldName;
    int atomSerialNumber;
    string atomName;
    string residueName;
    string chainIdentifier;
    int residueSequenceNumber;
    float x;
    float y;
    float z;
    float charge;
    float radius;
    string segmentIdentifier;
    float ener;
int operator<(const pqrLine &rhs) const
{
   if( this->x < rhs.x){ return 1; }
   else if ( this->x == rhs.x){
	   if (this->y < rhs.y){ return 1; }
	   else if (this->y == rhs.y){
		   if (this->z < rhs.z){ return 1; }
	   }
   }
   return 0;
}

  };
    

    list<pqrLine> pqrData;
  
  public:
    ///
    NumProcWritePotat (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcWritePotat();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcWritePotat (pde, flags);
    }
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "WritePotat";
    }
    
    void readpqrData();
    

  };





  NumProcWritePotat :: NumProcWritePotat (PDE & apde, const Flags & flags)
    : NumProc (apde), point(1)
  {
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""), 1); 
    lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", ""), 1);
    gfu  = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""), 0);
    gfv = pde.GetGridFunction (flags.GetStringFlag ("gridfunction2", ""), 1); 
    pqrfile = flags.GetStringFlag ("pqrfile","pqr");
    noOfAtoms = flags.GetNumFlag ("numberatoms", 1);
    potatfile = flags.GetStringFlag ("potatfile","potat");

    states = flags.GetNumFlag ("statenr", 1);
    nrOfValues = flags.GetNumFlag ("nrOfValues", noOfAtoms);
      
    //pqrData.SetSize(noOfAtoms);
      
    readpqrData();

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

  NumProcWritePotat :: ~NumProcWritePotat()
  {
    ;
  }
  
  void NumProcWritePotat::readpqrData()
  {
    ifstream inFile(pqrfile.c_str());
    
    list<pqrLine>::iterator it;
    while (inFile.good()){
       pqrLine t;

      inFile >> t.fieldName >> t.atomSerialNumber >> t.atomName >> t.residueName >> t.chainIdentifier >> t.residueSequenceNumber >> t.x >> t.y >> t.z >> t.charge >> t.radius >> t.segmentIdentifier;      
       if (t.fieldName == "ATOM"){
		pqrData.insert(pqrData.begin(), t);
       }
    
  }  
    //pqrData.sort();
    inFile.close();
  }



  void NumProcWritePotat :: PrintDoc (ostream & ost)

  {
    ost <<
      "\n\nNumproc WritePotat:\n" \
      "-----------------\n" \
      "(documentation not updated!)" \
      "Evaluates linearform or bilinearform, or pointvalues:\n"\
      "Required flags:\n" \
      "-gridfunction=<gfname>\n" \
      "    gridfunction to evaluate\n" \
      "\nOptional flags:\n" \
      "-bilinearform=<bfname>\n" \
      "    evaluates bilinear-form bfname(gfname,gfname2)\n"\
      "    needs second gridfunction <gfname2>\n" \
      "-gridfunction2=<gfname2>\n" \
      "    second gridfunciton for bilinear-form evaluation\n" \
      "-cachecomp=<n>"\
      "    for gridfunctions with cachesize > 1 use the given component"
      "-linearform=<lfname>\n" \
      "    evaluates linearform lfname(gfname)\n" \
      "-point=[x,y,z]\n" \
      "    evaluates diffop applied to gridfunction in point p\n" \
      "    diffop taken from first term in bilinear-form\n" \
      "-domains=[d1,...,dn]\n" \
      "    the point has to lie in one of the specified domains\n" \
      "-point2=[x2,y2,z2]\n" \
      "    evaluates diffop applied to gridfunction in line p1-p2, and writes to file\n" \
      "-applyd\n" \
      "    evaluates flux instead of derivatives\n" \
      "-integrateonplanes\n" \
      "    diffop applied to gridfunction is integrated on a sequence of rectangles bounded by a brick\n" \
      "       -point=... and -point2=... now set the minimal and maximal point of the brick\n" \
      "       -point3=..., -point4=...\n"\
      "            if these additional points are set, the bounds are not given by a brick but by\n" \
      "            a parallelepiped\n"\
      "       -variabledirection=1|2|3\n" \
      "            sets the direction of the normal vector of the rectangles (x|y|z)\n"\
      "       -n1=<number of points>\n" \
      "       -n2=<number of points>\n" \
      "       -n3=<number of points>\n" \
      "            set the number of integration points resp. the number of rectangles in the 3 directions\n" \
      "-text=blabla \n" \
      "    prints text \n" \
      " -resultvariabe=<varname> \n" \
      "    stores scalar, real results in variable <varname>\n" \
      "-filename=<fn> \n" \
      "    writes values in file \n" \
      "-outputprecision=<n> \n" \
      "    sets the number of digits for the output of the results\n";
  }




  void NumProcWritePotat :: Do(LocalHeap & lh)
  {
    float readOutShift = 0;

    double result = 0;
    ofstream ofile (filename.c_str());

    int old_cout_precision = cout.precision();
    int old_ofile_precision = ofile.precision();

    if(outputprecision > 0)
    {
    	cout.precision(outputprecision);
    	ofile.precision(outputprecision);
    }

    if (lff)
    {
    	if (!gfu)
	  throw Exception ("evaluate linear-form needs an argument -gridfunction=u");

    	cout << "<" << lff->GetName() << ", " << gfu->GetName() << "> = " << flush;
	
	
    	if (!lff->GetFESpace().IsComplex())
    	{
    		result = S_InnerProduct<double>(lff->GetVector(), gfu->GetVector());
    		cout << result << endl;
    	}
    	else
	  cout << S_InnerProduct<Complex>(lff->GetVector(), gfu->GetVector()) << endl;
	
    }
    else if (point.Size() >= 2)
    {
	  
      if (MyMPI_GetNTasks() == 1 || MyMPI_GetId() != 0){
	    
    	  const BilinearFormIntegrator & bfi = (bfa) ? *bfa->GetIntegrator(0) : *gfu->GetFESpace().GetIntegrator();

    	
    	ofstream potfile;

		
		  
    	cout << IM(1) << "looking for " << potatfile.c_str() << endl;

    	struct stat st;
    	if(stat(potatfile.c_str(),&st) != 0){
	  cout << IM(1)  << "file not found. creating and writing states" << endl;
	  potfile.open( potatfile.c_str());
    	} else {
	  potfile.open( potatfile.c_str(), ios::app );
    	}

    	if (states != -1)
	  potfile << states << endl;

    	cout << IM(1)  << "Writing noOfAtoms for next state: " << noOfAtoms << endl;
    	potfile << noOfAtoms << endl;

	cout << IM(1)  << "init searchtree" << endl;
	Vec<3> any_point(0,0,0);
	IntegrationPoint ip;
	ma.FindElementOfPoint(any_point,ip,true);



    	list<pqrLine>::iterator it;
    	int position = 1;
    	pqrData.sort();
     
    	for (it=pqrData.begin(); it!=pqrData.end(); ++it){
	  pqrLine line;
	  line = *it;
	  float realX = line.x;
	  point(0) = line.x + readOutShift;
	  point(1) = line.y;
	  point(2) = line.z;
		
	  FlatVector<double> pflux(bfi.DimFlux(), lh);

	  if (!gfu->GetFESpace().IsComplex())
	  {

	    CalcPointFlux (ma, *gfu, point, domains, pflux, bfi, applyd, lh, component);
	      
	    result = pflux(0);
	    // Writing [V] into potat!
	    potfile.precision(-1);
	    potfile << realX << endl;
	    potfile << point(1) << endl;
	    potfile << point(2) << endl;
	    potfile.precision(6);
	    potfile << (1.60217646e-19 * 6.0221415e23 * 1 * pflux(0))/1000*1/(0.008314510*300) << endl;
	    position++;
    	}
    }
    potfile.close();
   }
  }
    
    pde.GetVariable(variablename,true) = result;
    
    cout.precision(old_cout_precision);
    ofile.precision(old_ofile_precision);
    ofile.close();
  }

  void NumProcWritePotat :: PrintReport (ostream & ost)
  {
    ost << "NumProcWritePotat:" << endl;
  }


namespace WritePotat_cpp
{
  class Init
  { 
  public: 
    Init ();
  };
    
  Init::Init()
  {
    GetNumProcs().AddNumProc ("writepotat", NumProcWritePotat::Create, NumProcWritePotat::PrintDoc);
  }
    
  Init init;
}
  
