// implement a write out of potential map or potential map difference
 
#include <solve.hpp>
using namespace ngsolve;
#include <vector>
#include <deque>
#include <limits>

typedef std::numeric_limits< double > dbl;

class NumProcWriteDX : public NumProc
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
    bool showsteps;
    bool dx;

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
  float m_length;
  float m_ox, m_oy, m_oz;
  float m_h;
  string dxfile;
  float m_boxlength;
  
  public:
    ///
    NumProcWriteDX (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcWriteDX();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcWriteDX (pde, flags);
    }
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "writedx";
    }
    
  void readMolecule(float &m_length, float &m_ox, float &m_oy, float &m_oz);

  };


NumProcWriteDX :: NumProcWriteDX (PDE & apde, const Flags & flags)
: NumProc (apde), point(1)
{
	bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""), 1); 
    lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", ""), 1);
    gfu  = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""), 0);
    gfv = pde.GetGridFunction (flags.GetStringFlag ("gridfunction2", ""), 1); 
    pqrfile = flags.GetStringFlag ("pqrfile","pqr");
    dxfile = flags.GetStringFlag ("dxfile","output.dx");
    m_h = flags.GetNumFlag ("h",1);
    m_ox = flags.GetNumFlag ("ox",-1);
    m_oy = flags.GetNumFlag ("oy",-1);
    m_oz = flags.GetNumFlag ("oz",-1);
    m_boxlength = flags.GetNumFlag ("length",-1);

   
    showsteps = flags.GetDefineFlag("showsteps");


    dx = flags.GetDefineFlag ("dx");
    
    readMolecule(m_length, m_ox, m_oy, m_oz);
    cout << "DX properties:" << endl;
    cout << "molecule max length: " << m_length << endl;
    cout << "origin set to:   (" << m_ox << ", " << m_oy << ", " << m_oz << ")" << endl; 
    cout << "h = " << m_h << endl;

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

NumProcWriteDX :: ~NumProcWriteDX()
{
    ;
}
  
void NumProcWriteDX::readMolecule(float &m_length, float &m_ox, float &m_oy, float &m_oz)
{
  string line;
  ifstream in(pqrfile.c_str());

    if (!in) {
    	in.close();
    	cout << "Cannot open pqr file:" << pqrfile << endl;
    	exit(0);
    }

    float minx, miny, minz, maxx, maxy, maxz;
    minx = miny = minz = 99999999;
    maxx = maxy = maxz = -9999999;
    float ox, oy, oz;
    ox = oy = oz = 0;
    float max_vdW = -999;

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

		ox += t.x;
		oy += t.y;
		oz += t.z;

		if (t.x >= maxx)
		  maxx = t.x;
		if (t.x <= minx)
		  minx = t.x;

		if (t.y >= maxy)
		  maxy = t.y;
		if (t.y <= miny)
		  miny = t.y;

		if (t.z >= maxz)
		  maxz = t.z;
		if (t.z <= minz)
		  minz = t.z;

		// charge
		t.charge = atof(line.substr(54,6).c_str());
		   
		// radius
		t.radius = atof(line.substr(60,6).c_str());

		if (t.radius >= max_vdW)
		  max_vdW = t.radius;

		molecule.push_back(t);
		if (showsteps)
			cout << i << ": " << t.fieldName << " " << t.x << " " << t.y  << " " << t.z  << " " << t.charge  << " " << t.radius <<  endl;
    }    

    float lx, ly, lz;
    lx = maxx - minx;
    ly = maxy - miny;
    lz = maxz - minz;
    float length = 0;

    if (lx >= ly && lx >= lz)
      length = lx;
    else if (ly > lx && ly > lz)
      length = ly;
    else 
      length = lz;

    length += max_vdW*2+2;

    //    if (m_length == -1)
      m_length = length;

    cout << "#atoms: " << i << endl;
    cout << "molecule origin: " << ox/i << ", " << oy/i << ", " << oz/i << endl;

    if (m_ox == -1)
      m_ox = ox/i; //-(lx/2);
    if (m_oy == -1)
      m_oy = oy/i; //-(ly/2);
    if (m_oz == -1)
      m_oz = oz/i; //-(lz/2);

}



void NumProcWriteDX :: PrintDoc (ostream & ost)
{
	ost << "(documentation not updated!)";
}

void NumProcWriteDX :: Do(LocalHeap & lh)
{

  double temp_result = 0;
  double potential = 0;
  double V2kT_div_e = 1/2.58519972e-2;

 
  deque<float> result_left;

  const BilinearFormIntegrator & bfi = (bfa) ? *bfa->GetIntegrator(0) : *gfu->GetFESpace().GetIntegrator();
	
  cout.precision(dbl::digits10);
  cout << "init searchtree" << endl;
  Vec<3> any_point(0,0,0);
  IntegrationPoint ip;
  ma.FindElementOfPoint(any_point,ip,true);

  if (m_boxlength == -1)
    m_boxlength = m_length;

  if (m_boxlength < m_length)
    m_boxlength = m_length;

  const int n = ceil(m_boxlength/m_h)+1;

  cout << "xx=yy=zz: " << n << ", m_h: " << m_h << endl;
  cout << "boxlength = " << m_boxlength << ", molecule mod length: " << m_length << endl;
  FlatVector<double> pflux(bfi.DimFlux(), lh);
  for(unsigned int xx=0; xx<n; xx++) {
    
    cout << "\rGetting potentials at plane " << xx+1 << "/" << n << flush;
    for(unsigned int yy=0; yy<n; yy++) {
      for(unsigned int zz=0; zz<n; zz++) {
	point(0) = ((n-xx-1)*m_h)+m_ox-m_boxlength/2;
	point(1) = ((n-yy-1)*m_h)+m_oy-m_boxlength/2;
	point(2) = ((n-zz-1)*m_h)+m_oz-m_boxlength/2;
	



	//		cout << xx+1 << " -> looking at: " << point(0) << ", " << point(1) << ", " << point(2) << endl;

	CalcPointFlux (ma, *gfu, point, domains,
		       pflux, bfi, applyd, lh, component);
	      
	result_left.push_back(pflux(0)*V2kT_div_e);
      }
    }
  }
  cout << endl;
  if (MyMPI_GetNTasks() == 1 || MyMPI_GetId() != 0){
    if (dx){
      cout << "Write out dx to " << dxfile << " ... ";
      int i = 0;
      ofstream resultFile;
      resultFile.open(dxfile.c_str());
      resultFile << "# Generated with mFES" << endl;
      resultFile << "# Units: kT/e" << endl;
      resultFile << "object 1 class gridpositions counts " << n << " " << n << " " << n << endl;
      resultFile << "origin " << m_ox-(m_boxlength/2) << " " << m_oy-(m_boxlength/2) << " " << m_oz-(m_boxlength/2) << endl;
      resultFile << "delta " << m_h << " 0 0" << endl;
      resultFile << "delta 0 "<< m_h << " 0" << endl;
      resultFile << "delta 0 0 " << m_h << endl;
      resultFile << "object 2 class gridconnections counts " << n << " " << n << " " << n << endl;
      resultFile << "object 3 class array type double rank 0 items " << n*n*n << " data follows" << endl;
      
      for (int xx = 0; xx < n; xx++){
      cout << "\rWriting potentials at plane " << xx+1 << "/" << n << flush;
	for (int yy = 0; yy < n; yy++){
	  for (int zz = 0; zz < n; zz++){
	    i++;
	    resultFile << result_left.back() << " ";
	    result_left.pop_back();
	    if (i % 3 == 0)
	      resultFile << "\n";
	  }
	}
      }
      if (i%3 != 0)
	resultFile << "\n";

      resultFile << "attribute \"dep\" string \"positions\"" << endl;
      resultFile << "object \"regular positions regular connections\" class field" << endl;
      resultFile << "component \"positions\" value 1" << endl;
      resultFile << "component \"connections\" value 2" << endl;
      resultFile << "component \"data\" value 3" << endl;
      
      resultFile.close();
      cout << "..finished!" << endl;
    }
  }

}

void NumProcWriteDX :: PrintReport (ostream & ost)
{
	ost << "NumProcWriteDX:" << endl;
}

static RegisterNumProc<NumProcWriteDX> init_np_pts("writedx");
