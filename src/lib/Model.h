/** @file Model.h
 *  @brief Different models for meshing are set up within this file.
 *
 *  A surface or volume mesh is set up using the one- or two-cycle approach.
 *  Models are either the molecular surface, the ion exclusion layer of the 
 *  protein or the surfaces for the titratable groups (if using the two-cycle)
 *  approach.
 *
 *  @author Ilkay Sakalli
 */

//#include <solve.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/timer.hpp>

#include "LSMS.h"
#include "Voxel.h"
#include "VCG.h"


using namespace std;
using namespace vcg;

#include <meshing.hpp>  

namespace nglib {

#include <nglib.h>


  


  // Changes the bc number of a VOL mesh
  DLL_HEADER Ng_Mesh * Ng_SetProperties(Ng_Mesh * mesh, int surfnr, int bcnr, int domin, int domout)
  {
    netgen::Mesh * m = (netgen::Mesh*) mesh;
    int maxsteps = m->GetNSE();
    SurfaceElementIndex sei;
    for (sei = 0; sei < maxsteps; sei++)
      {
	//  if ((*m)[sei].GetIndex()){
	m->GetFaceDescriptor((*m)[sei].GetIndex ()).SetSurfNr (surfnr);
	m->GetFaceDescriptor((*m)[sei].GetIndex ()).SetBCProperty (bcnr);
	m->GetFaceDescriptor((*m)[sei].GetIndex ()).SetDomainIn (domin);
	m->GetFaceDescriptor((*m)[sei].GetIndex ()).SetDomainOut (domout);
	//  }
	 
      }
    return ( (Ng_Mesh*)m );
  }
 


  // Merge another mesh file into the currently loaded one
  DLL_HEADER Ng_Result Ng_Alt_MergeMesh( Ng_Mesh* mesh1, Ng_Mesh* mesh2)
  {

    Ng_Result status = NG_OK;
    Mesh * m = (Mesh*)mesh1;
    Mesh * m2 = (Mesh*)mesh2;

    if(!m)
      {
	status = NG_ERROR;
      }

    if(status == NG_OK)
      {
	const int num_pts = m->GetNP();
	int surfindex_offset = 0;

	int i;

	int oldnp = m->GetNP();
	int oldne = m->GetNSeg();

	for(SurfaceElementIndex si = 0; si < m->GetNSE(); si++)
	  for(int j=1; j<=(*m)[si].GetNP(); j++) (*m)[si].GeomInfoPi(j).trignum = -1;

	int max_surfnr = 0;
	for (i = 1; i <= m->GetNFD(); i++)
	  max_surfnr = netgen::max2 (max_surfnr, m->GetFaceDescriptor(i).SurfNr());
	max_surfnr++;

	if(max_surfnr < surfindex_offset) max_surfnr = surfindex_offset;

	         
	int maxsteps = m2->GetNSE();
	SurfaceElementIndex sei;
	for (sei = 0; sei < maxsteps; sei++)
	  {
	    int j;
	    int surfnr, bcp, domin, domout, nep = 3, faceind = 0;
	     
	    surfnr = m2->GetFaceDescriptor((*m2)[sei].GetIndex ()).SurfNr()+1;
	    bcp = m2->GetFaceDescriptor((*m2)[sei].GetIndex ()).BCProperty()+1;
	    domin = m2->GetFaceDescriptor((*m2)[sei].GetIndex ()).DomainIn()+1; 
	    domout = m2->GetFaceDescriptor((*m2)[sei].GetIndex ()).DomainOut()+1; 

	    for (j = 1; j <= m->GetNFD(); j++)
	      if (m->GetFaceDescriptor(j).SurfNr() == surfnr &&
		  m->GetFaceDescriptor(j).BCProperty() == bcp &&
		  m->GetFaceDescriptor(j).DomainIn() == domin &&
		  m->GetFaceDescriptor(j).DomainOut() == domout)
		faceind = j;

	    if (!faceind)
	      {
		faceind = m->AddFaceDescriptor (FaceDescriptor(surfnr, domin, domout, 0));
		if(m->GetDimension() == 2) bcp++;
		m->GetFaceDescriptor(faceind).SetBCProperty (bcp);
	      }

	    Element2d tri(nep);
	    tri.SetIndex(faceind);

	    Element2d sel = (*m2)[sei];
	                      
	    for (j = 1; j <= nep; j++)
	      {
		tri.PNum(j) = sel.PNum(j) + oldnp;
	      }


	    for (j = 1; j <= nep; j++)
	      {
		tri.GeomInfoPi(j).trignum = -1;
	      }

	    m->AddSurfaceElement (tri);
	  }
	        
	for (ElementIndex ei = 0; ei < m2->GetNE(); ei++)
	  {
	    Element el;
	      
	    int hi = (*m2)[ei].GetIndex() + 1;
	    if (hi == 0) hi = 1;
	    el.SetIndex(hi);
	      
	    int nep = (*m2)[ei].GetNP();
	    el.SetNP(nep);

	    for (int j = 0; j < nep; j++)
	      el[j] = (*m2)[ei][j] + oldnp;

	    m->AddVolumeElement (el);
	  }

	         
	PointIndex pi;
	maxsteps = m2->GetNP();
	for (pi = PointIndex::BASE; pi < maxsteps+PointIndex::BASE; pi++)
	  {
	    netgen::Point3d p;
	    p.X() = (*m2)[pi](0)/1.0;
	    p.Y() = (*m2)[pi](1)/1.0;
	    p.Z() = (*m2)[pi](2)/1.0;
	    m->AddPoint (p);
	     
	  }  
	         
	for (unsigned int i = 1; i <= m2->GetNSeg(); i++){
	  Segment & seg = m2->LineSegment (i);
	   
	  seg.surfnr1--;
	  seg.surfnr2--;
	  if(seg.surfnr1 >= 0)  seg.surfnr1 = seg.surfnr1 + max_surfnr;
	  if(seg.surfnr2 >= 0)  seg.surfnr2 = seg.surfnr2 + max_surfnr;
	  seg[0] = seg[0] +oldnp;
	  seg[1] = seg[1] +oldnp;
	  seg.edgenr = seg.edgenr + oldne;
	  seg.epgeominfo[1].edgenr = seg.epgeominfo[1].edgenr + oldne;
	  m->AddSegment (seg);
	}
	         

	m->CalcSurfacesOfNode ();

        /*
	m->topology -> Update();
	m->clusters -> Update();
        */
        m->UpdateTopology();
	m->SetNextMajorTimeStamp();

	if(m->GetNP() > num_pts)
	  {
	    status = NG_OK;
	  }
	else
	  {
	    status = NG_ERROR;
	  }
      }
    return status;

  }

}
using namespace nglib;


typedef boost::property_tree::ptree INI;

/**
 * @class Model
 *
 * @brief A model is a surface for the protein or titratable group depending on 
 *        the approach used. Also building up the whole volume models merging
 *        with the boundary is implemented here using NETGEN.
 *
 * A model is either a surface or a volume. A surface is generated for the molecular 
 * surface or the ion exclusion layer for the protein or a surface is generated
 * for every titratable group (also one for C- and N-terminus if selected). Whole
 * surface management, also smoothing and merging is performed using NETGEN file
 * structure conversion to VCG lib and vice versa. Exporting surfaces at every step
 * is also possible. First second order surface generation methods are commented out.
 * Meshing options are carefully handled using OPT files defined in the config INI 
 * files which is read in by mFES at starting point.
 *
 * @author Ilkay Sakalli
 *
 */ 

class Model {
public:
  /// True if the ion exclusion surface be computed. Is the case if ion concentration 
  /// ionc != 0.
  static bool exclusionSurface;

  /** @brief The ion exclusion layer is computed using an atom list, INI options 
   *         file and an output file (exclusion.vol)
   *
   *  The ion exclusion layer (IEL) is computed for the protein and titratable groups.
   *  The IEL is often defined as the van-der-Waals radius which is inflated by the
   *  radius of ions. mFES uses same definition but also rolls a probe sphere over the
   *  inflated van-der-Waals surface with a given radius (default: 0.5 Angstroems).
   *  This is done for robustness and it seems to be plausible by physical means.
   *
   *  @param atomList List of atoms be modeled.
   *  @param ini INI file which defines e.g. the probe radius used for the
   *                 ion exclusion layer.
   *  @param fname Name of the volume file which will be generated as a result.
   *               (default: exclusion.vol)
   *  @return Void.
   */
  void calcIonLayer(vector<Atom> &atomList, INI &ini, string fname = "exclusion.vol"){
    LSMS lmSurface;
    mMesh mSurface;
    int mask = 1;
    
    /// The file exclusion.stl has to be computed before using LSMS. If it is there
    /// use it, otherwise compute it using LSMS and export is as exclusion.stl.
    if (!boost::filesystem::exists( "exclusion.stl" )){
      cout << "Calculating ion exclusion layer " << fname << "..." << endl;
      boost::timer t;
      ofstream time;
      /// A file called "times" if created before to perform some timing and
      /// writing out to file
      time.open ("times");
      /// The ion exlcusion layer is computed here.
      if (lmSurface.calcMC(mSurface, atomList, ini, "exclusion")){	      
	/// The surface is cleaned which means removing duplicated points and vertices
	clean(mSurface);
	/// The surface is smoothed WITHOUT shrinkage using Taubin smoothing or other
	/// method defined in the INI config file.
	smooth(mSurface, ini);
	/// Resulting surface is exported as exclusion.stl to be read in later on
	/// or reused in a next computation.
	tri::io::ExporterSTL<mMesh>::Save(mSurface,"exclusion.stl",false);

	/// Timing is written out to times and also to console
	double currentTime = t.elapsed();
	time << "exclusion_surface " << currentTime << " s" << endl;
	cout << "mFES: ion exclusion layer " << fname << "surface calculation took " << currentTime << " seconds." <<endl;
    	time.close();
      }
    }
  }


  /** @brief A model (surface) is computed. This may be the molecular surface of 
   *         the protein, a titratable group or the cavity surface.
   *
   *  A model is computed using NETGEN. A new molecular surface is computed using the 
   *  a fine grained template with the help of the advancing front method. Here, a 
   *  model may be the molecular surface as a vol file for the protein, the titratable
   *  group or the cavity.
   *
   *  @param atomList List of atoms be modeled.
   *  @param ini INI file which defines e.g. the probe radius used for the
   *                 ion exclusion layer.
   *  @param fname Name of the volume file which will be generated as a result.
   *  @param cavity If the cavity should be computed, this flag is set to true.
   *  @return Void.
   */
  void calcModel(vector<Atom> &atomList, INI &ini, string fname = "", bool cavity = false){

    /// Template surface generated using LSMS
    LSMS lmSurface;
    /// Surface conversion to use with VCG lib
    mMesh mSurface;
	  
    /// Volume file name (VOL file)
    string volume_vol = ini.get<string>("model.volume_vol");
    /// Surface file name (STL file)
    string surface_stl = ini.get<string>("model.surface_stl");
    /// Generator to use for protien. Either standard (default; may be always used) 
    /// or voxelizer (used for small molecules, may more accurate but cannot handle 
    /// cavities)
    string generator = ini.get<string>("model.generator");
    /// Resolution to use for the surface generator (default: 512)
    int generatorResolution = atoi(ini.get<string>("model.grid_resolution").c_str());

    /// Generator to use for titratable groups
    string generatorResidue = ini.get_optional<string>("model.generator_residue").get_value_or("standard");
    /// Resolution to use for the surface generator (default: 256; because this is a
    /// smaller molecule and not a protein)
    int generatorModelResolution = atoi(ini.get_optional<string>("model.grid_residue_resolution").get_value_or("256").c_str());

    /// ionc is the ion concentration given in M
    boost::optional<string> ionc = ini.get_optional<string>("experiment.ionc");
    
    /// if ionc is set and it is not zero compute the ion exclusion surfaces for every
    /// model
    if(ionc){
      string ionc = ini.get<string>("experiment.ionc");
      if (atof(ionc.c_str()) != 0){
	exclusionSurface = true;
      }
    }
    

    int mask = 1;
    if (generator == "standard" && cavity && !boost::filesystem::exists( "cavity.vol" )){
      /// Cavity is calculated. A cavity can just be computed using the generator
      /// "standard" which is default behaviour
      cout << "Calculating cavities..." << endl;

      /// Perform some timing on the computation of the cavity
      boost::timer t;
      ofstream time;
      time.open ("times");
      
      /// First compute the cavity using LSMS and generate a template surface
      if (lmSurface.calcMC(mSurface, atomList, ini, "cavity")){
	/// Clean and smooth the surface smoothly
	clean(mSurface);
	smooth(mSurface, ini);
	
	/// Export the resulting surface as "cavity.stl" (default)
	if (surface_stl != "" )
	  tri::io::ExporterSTL<mMesh>::Save(mSurface,"cavity.stl",false);
	
	double currentTime = t.elapsed();
	time << "cavity_surface " << currentTime << " s" << endl;
	cout << "mFES: cavity model surface calculation took " << currentTime << " seconds." <<endl;
	
	t.restart();
	/// Here, the template cavity surface is regularized using the 
	/// advancing front method perfomed by NETGEN
	convert(mSurface, ini, "cavity.vol", "cavity");

	/// Timing is performed and written out to times and console
	currentTime = t.elapsed();
	time << "cavity_volume " << currentTime << " s" << endl;
	cout << "mFES: cavity model VOL meshing took " << currentTime << " seconds." <<endl;
      }
      time.close();
      
    } else if ( fname == "" && !boost::filesystem::exists( volume_vol ) ){
      // Protein model is calculated
      boost::timer t;
      ofstream time;
      time.open ("times");

      /// Test if the protein molecular surface was already be computed.
      if (boost::filesystem::exists( surface_stl )) {
	/// If the molecular surface of the protein already exists, import it.
	cout << surface_stl << " found. " << endl;
	tri::io::ImporterSTL<mMesh>::Open(mSurface, surface_stl.c_str(), mask);
      } else {
	/// The protein molecular surface was not computed yet.
	if (generator == "standard"){
	  /// If the generator is standard, compute the molecular surface using 
	  /// LSMS tSurface and convert it to mMesh the VCG surface data structure
	  LSMS tSurface;
	  mMesh tmSurface;
	  /// Here, the LSMS computation is performed using the whole atom list
	  /// and options defined in config INI.
	  tSurface.calcMC(tmSurface, atomList, ini);
	  /// The surface is cleaned and smoothed leaving the surface shape invariant
	  clean(tmSurface);
	  smooth(tmSurface, ini);
	  /// Define a standard output file name for the fine grained molecular surface
	  /// (default: protein.stl)
	  if (surface_stl == "" )
	    surface_stl = "protein.stl";
	  
	  /// Export STL
	  tri::io::ExporterSTL<mMesh>::Save(tmSurface,surface_stl.c_str(),false);
	  /// Import STL to VCG lib structure
	  /// This is done to be sure that the structure is loaded correctly
	  /// and is easier tom implement
	  tri::io::ImporterSTL<mMesh>::Open(mSurface, surface_stl.c_str(), mask);
	} else {
	  /// If the generator is Voxelizer
	  Voxel vSurface;
	  /// If the molecular surface does not yet exist, compute it using the 
	  /// voxelizer
	  if (!boost::filesystem::exists( "protein.stl" )){
	    /// Here, the voxelizer is used to compute the molecular surface.
	    /// and the result is written into protein.stl
	    vSurface.calcSurface(mSurface, atomList, ini, "protein.stl", generatorResolution);
	    
	    /// The surface is cleaned as before (removing duplicates, etc.)
	    clean(mSurface);
	    
	  } 
	  /// Afterwards the output (default: protein.stl) is read in again as VCG lib
	  /// molecular surface data structure.
	  tri::io::ImporterSTL<mMesh>::Open(mSurface, string("protein.stl").c_str(), mask);			     
	}
	
      }

      /// Timing is performed and written into times and output to console.
      double currentTime = t.elapsed();
      time << "protein_surface " << currentTime << " s" << endl;
      cout << "mFES: protein model surface calculation took " << currentTime << " seconds." <<endl;
	    
      t.restart();
      /// Here, the STL molecular surface is converted to a VOL file (NETGEN
      /// datastructure performing the advancing front method and using the 
      /// parameters defined in config INI and options files OPT.
      convert(mSurface, ini);
      currentTime = t.elapsed();
      time << "protein_volume " << currentTime << " s" << endl;
      cout << "mFES: protein model VOL meshing took " << currentTime << " seconds." <<endl;
      time.close();
      
    } else if (fname != ""){
      /// Molecular surface for titratable group with id nr will be calculated
      /// Optionally the ion exclusion layer is computed
      string exclusionstlFile = fname+"_exclusion.stl";
      string exclusionGroups = ini.get<string>("pka.explicit_models");

      /// If the ion exclusion layer surface for the specific titratable group
      /// is not yet computed, compute it after choosing the appropriate surface 
      /// generator (standard or voxelizer).
      if (!boost::filesystem::exists( exclusionstlFile.c_str() ) ){
	cout << "group model calculation " << fname << endl;
	
	if (generatorResidue == "standard"){
	  /// Here, the ion exclusion layer is computed using the mode
	  /// "residue_exclusion" in the LSMS wrapper method.
	  lmSurface.calcMC(mSurface, atomList, ini, "residue_exclusion");
	  clean(mSurface);
	  smooth(mSurface, ini);
	  /// After cleaning and smoothing the molecular surface leaving the 
	  /// surface invariant, the surface is exported as an STL file.
	  tri::io::ExporterSTL<mMesh>::Save(mSurface,exclusionstlFile.c_str(),false);
	} else {
	  /// If the generator is Voxelizer, compute it with this generator
	  /// and save the surface as STL.
	  Voxel vSurface;
	  vSurface.calcSurface(mSurface, atomList, ini, exclusionstlFile, generatorModelResolution, true);
	  clean(mSurface);
	  /// Import the surface created into a VCG lib surface datastructure
	  tri::io::ImporterSTL<mMesh>::Open(mSurface, exclusionstlFile.c_str(), mask);
	}
	
      }
      
      /// After the ion exclusion layer is computed, create the molecular surface of
      /// the titratable group.
      if ( !boost::filesystem::exists( fname )  ){
	string stlFile = fname+".stl";
	ofstream time;

	time.open ("times", ios::app);
	boost::timer t;
	// group model calculation
	if ( !boost::filesystem::exists( stlFile )){

	  /// If the generator is "standard" (default), compute the molecular surface
	  /// using the LSMS wrapper method.
	  if (generatorResidue == "standard"){
	    lmSurface.calcMC(mSurface, atomList, ini, "residue");
	    clean(mSurface);
	    smooth(mSurface, ini);
	    /// After cleaning and smoothing the molecular surface, export the result
	    /// as an STL file.
	    tri::io::ExporterSTL<mMesh>::Save(mSurface,stlFile.c_str(),false);
	  } else {
	    /// If the generator "Voxelizer" is chosen, compute the molecular surface
	    /// export it as an STL file and import it as an VCG lib surface.
	    Voxel vSurface;
	    vSurface.calcSurface(mSurface, atomList, ini, stlFile, generatorModelResolution);
	    clean(mSurface);
	    tri::io::ImporterSTL<mMesh>::Open(mSurface, stlFile.c_str(), mask);
	  }
	  
	  double currentTime = t.elapsed();
	  time << "model_surface " << currentTime << " s" << endl;
	  cout << "mFES: residue model surface calculation took " << currentTime << " seconds." <<endl;
	  
	  t.restart();
	  /// Here, the molecular surface is generated. Now, we need to regularize
	  /// the molecular surface using NETGEN by performing the advancing
	  /// front method. Options defined in the INI file are used.
	  convert(mSurface, ini, fname, "residue");
	  /// Elapsed time for timing purposes is written into times and
	  /// also to the console.
	  currentTime = t.elapsed();
	  time << "model_volume " << currentTime << " s" << endl;
	  cout << "mFES: residue model VOL meshing took " << currentTime << " seconds." <<endl;
	  time.close();
	}
	else {
	  /// The current molecular surface for the titratable group was computed
	  /// previously and saved as an STL file. So it just needs to be imported
	  /// as a VCG lib surface.
	  cout << "model surface found: " << stlFile << endl;
	  tri::io::ImporterSTL<mMesh>::Open(mSurface, stlFile.c_str(), mask);
	  t.restart();
	  /// ... and converted into a VOL file (NETGEN) regularizing the surface
	  /// using the advancing front method and the meshing options defined
	  /// in the INI options files.
	  convert(mSurface, ini, fname, "residue");
	  double currentTime = t.elapsed();
	  time << "model_volume " << currentTime << " s" << endl;
	  cout << "mFES: residue model VOL meshing took " << currentTime << " seconds." <<endl;
	  time.close();
	}
      }
    }    
  }


 private:

  /** @brief Nearly all available meshing options are printed out for debugging purposes.
   *
   *  Meshing options are provided by NETGEN. One may play with them using the 
   *  NETGEN GUI. Here we use options which are found out to be useful for our
   *  purposes. 
   *
   *  @param Ng_Meshing_Parameters Data structure to store meshing parameters
   *                               for NETGEN.
   *  @param string This is some text to get sure, which options are used for
   *                which model step.
   *  @return Void.
   */  
  void printMeshingOptions(Ng_Meshing_Parameters &mp, string prefix){
    cout << endl;
    cout << prefix << endl;
    cout << "========================" << endl;
    cout << "options.localh = " << mp.uselocalh << endl;
    cout << "options.meshsize = " << mp.maxh << endl;
    cout << "options.minmeshsize = " << mp.minh << endl;
    cout << "meshoptions.fineness = " << mp.fineness << endl;
    cout << "************** blockfill currently not available in netgen 6 ******* " << endl;
    /*
    cout << "meshoptions.blockfill = " << mp.blockfill << endl;
    cout << "meshoptions.filldist = " << mp.filldist << endl;
    */
    cout << "options.grading = " << mp.grading << endl;
    cout << "options.curvaturesafety = " << mp.elementspercurve << endl;
    cout << "options.segmentsperedge = " << mp.elementsperedge << endl;
    cout << "options.secondorder = " << mp.second_order << endl;
    cout << "options.quad = " << mp.quad_dominated << endl;
    cout << "stloptions.resthcloseedgeenable = " << mp.closeedgeenable << endl;
    cout << "stloptions.resthcloseedgefac = " << mp.closeedgefact << endl;
    cout << "options.optsteps2d = " << mp.optsteps_2d << endl;
    cout << "options.optsteps3d = " << mp.optsteps_3d << endl;
    cout << "options.inverttets = " << mp.invert_tets << endl;
    cout << "options.inverttrigs = " << mp.invert_trigs << endl;
    cout << "options.checkoverlap = " << mp.check_overlap << endl << endl;
  }

  /** @brief Meshing options are set up using an options file
   *
   *  An OPT file is defined for every meshing step. This is the meshing of 
   *  the surface of a protein or titratable and and this is the meshing of 
   *  the volume of a protein or titratable group. An options file has a 
   *  variable and a value. The value is assigned to the variable name which
   *  is one of the meshing parameters.
   *
   *  @param Ng_Meshing_Parameters Data structure to store meshing parameters
   *                               for NETGEN.
   *  @param string options file (OPT) to be read in and parsed so that it is 
   *                used for specific meshing step.
   *  @return Void.
   */  
  void setMeshingOptions(Ng_Meshing_Parameters &mp, string optFile){
    /// Open the OPT file
    ifstream in(optFile.c_str());
    if (!in) {
      in.close();
      cout << "Cannot open mesh options file:" << optFile << endl;
      exit(0);
    }
    string pqrline;
    string variable;
    string value;
    
    while( !in.eof() ) {
      getline(in, pqrline);
      
      istringstream ss(pqrline);
      /// For every variable and value pair
      ss >> variable >> value;
      
      /// If the variable has a specific name, fill the specific 
      /// Meshing property with the given value.
        /*
      if (variable == "meshoptions.blockfill"){
	mp.blockfill = atoi(value.c_str());
      }
      
      if (variable == "meshoptions.filldist"){
	mp.filldist = atof(value.c_str());
      }
      */
      cout << "************** blockfill currently not available in netgen 6 ******* " << endl;      
      
      if (variable == "options.localh"){
	/// Switch to enable / disable usage of local mesh size modifiers
	mp.uselocalh = atoi(value.c_str());
	continue;
      }
      
      if (variable == "options.meshsize"){
	/// Maximum global mesh size allowed
	mp.maxh = atof(value.c_str());
	continue;
      }
      
      if (variable == "options.minmeshsize"){
	/// Minimum global mesh size allowed
	mp.minh = atof(value.c_str());
	continue;
      }
      
      if (variable == "meshoptions.fineness"){
	/// Mesh density: 0...1 (0 => coarse; 1 => fine)
	mp.fineness = atof(value.c_str());
	continue;
      }
      
      if (variable == "options.grading"){
	/// Mesh grading: 0...1 (0 => uniform mesh; 1 => aggressive local grading)
	mp.grading = atof(value.c_str());
	continue;
      }
      
      
      if (variable == "options.curvaturesafety"){
	/// Elements to generate per curvature radius
	mp.elementspercurve = atoi(value.c_str());
	continue;
      }
      
      if (variable == "options.segmentsperedge"){
	/// Number of elements to generate per edge of the geometry
	mp.elementsperedge = atof(value.c_str());
	continue;
      }
      
      if (variable == "options.secondorder"){
	/// Generate second-order surface and volume elements
	mp.second_order = atoi(value.c_str());
	continue;
      }
      
      if (variable == "options.quad"){
	/// Creates a Quad-dominated mesh
	mp.quad_dominated = atoi(value.c_str());
	continue;
      }
      
//			if (variable == "options.meshsizefilename"){
//				//!< Optional external mesh size file
//				mp.meshsize_filename = value.c_str();
//				continue;
//			}

      if (variable == "stloptions.resthcloseedgeenable"){
	/// Enable / Disable mesh refinement at close edges
	mp.closeedgeenable = atoi(value.c_str());
	continue;
      }
      
      if (variable == "stloptions.resthcloseedgefac"){
	/// Factor to use for refinement at close edges (larger => finer)
	mp.closeedgefact = atof(value.c_str());
	continue;
      }
      
      if (variable == "options.optsteps2d"){
	/// Number of optimize steps to use for 2-D mesh optimization
	mp.optsteps_2d = atoi(value.c_str());
	continue;
      }
      
      if (variable == "options.optsteps3d"){
	/// Number of optimize steps to use for 3-D mesh optimization
	mp.optsteps_3d = atoi(value.c_str());
	continue;
      }
      
      if (variable == "options.inverttets"){
	/// Invert all the volume elements
	mp.invert_tets = atoi(value.c_str());
	continue;
      }
      
      if (variable == "options.inverttrigs"){
	/// Invert all the surface triangle elements
	mp.invert_trigs = atoi(value.c_str());
	continue;
      }
      
      if (variable == "options.checkoverlap"){
	/// Check for overlapping surfaces during Surface meshing
	mp.check_overlap = atoi(value.c_str());
	/// Check for overlapping surface elements before volume meshing
	mp.check_overlapping_boundary = atoi(value.c_str());
      }
      
      
    }
    // optimize by default is 1, because if number of
    // optsteps is 0, then no optimization made
    
    in.close();
  }
  
  /** @brief A fine grained mesh is converted into a NETGEN volume mesh 
   *         including a boundary surface.
   *
   *  A VCG mMesh structure is converted into an Ng_STL_Geometry. The geometry is 
   *  then regularized using parameters defined in OPT files using the advancing
   *  front algorithm implemented by Joachim Schoeberl. If a volume meshed cavity 
   *  surface is available it is merged with the molecular surface. The protein is
   *  volume meshed and merged with an ion exclusion layer if it is available. 
   *  If the protein is merged with the ion exclusion layer, the volume between the
   *  molecular surface and the exclusion layer is also meshed. In the end, the
   *  whole preceding model is merged with a boundary layer and volume meshed again.
   *  The whole molecular model is saved as defined in the config INI file.
   *
   *  @param mMesh A VCG lib surface structure.
   *  @param INI The options defined in the config INI files are read in
   *  @param fname Suffix of the filename of the resulting molecular model VOL file
   *  @param mode Defining which type of surface is going to be converted 
   *              (possible values: protein / cavity / residue).
   *  @return Integer 1 if everything went fine. Now, just an exit occurs if there
   *          there was a failure.
   */  
  int convert(mMesh &mSurface, INI &ini, string fname = "", string mode = "protein")
  {
    
    /// NETGEN stl geometry structure
    Ng_STL_Geometry *stl_geom = Ng_STL_NewGeometry();
    
    /// Initialise the Netgen Core library
    Ng_Init();
    
    /// Actually create the mesh structure
    ngVolume = Ng_NewMesh();
	  
    int np, ne;
    
    /// Convert VCG MyMesh to Ng_STL_Geometry
    
    double p1[3];
    double p2[3];
    double p3[3];
    double n[3];
    
    mMesh::FaceIterator   fi;
    
    Point3f p;
    int triangles = 0;
    for(fi=mSurface.face.begin(); fi!=mSurface.face.end(); ++fi) if( !(*fi).IsD() ){
	// For each triangle write the normal, the three coords and a short set to zero
	p.Import(vcg::NormalizedNormal(*fi));
	n[0] = p[0]; n[1] = p[1]; n[2] = p[2];
	
	p.Import((*fi).V(0)->P());
	p1[0] = p[0]; p1[1] = p[1]; p1[2] = p[2];
	p.Import((*fi).V(1)->P());
	p2[0] = p[0]; p2[1] = p[1]; p2[2] = p[2];
	p.Import((*fi).V(2)->P());
	p3[0] = p[0]; p3[1] = p[1]; p3[2] = p[2];
	triangles++;
	Ng_STL_AddTriangle(stl_geom, p1, p2, p3, n);
      }
    
    cout << triangles << " triangles in surface" << endl;
    if(!stl_geom) {
      cout << "Error reading VCG STL File" << endl;
      exit(1);
    }
    cout << "Successfully loaded VCG STL File" << endl;
    
    /// Set the Meshing Parameters to be used
    string debug = ini.get<string>("model.debug");
    string jobname = ini.get<string>("general.jobname");
    string meshMoleculeSurface, meshMoleculeVolume;
    if (mode == "protein" || mode == "cavity" ){
      meshMoleculeSurface = ini.get<string>("meshing.molecule_surface");
      meshMoleculeVolume  = ini.get<string>("meshing.molecule_volume");
    } else if (mode == "residue") {
      meshMoleculeSurface = ini.get<string>("meshing.residue_surface");
      meshMoleculeVolume  = ini.get<string>("meshing.residue_volume");
    }
    cout << "using: " << meshMoleculeSurface << endl;
    cout << "using: " << meshMoleculeVolume << endl;
    
    Ng_Meshing_Parameters mp;
    setMeshingOptions(mp, meshMoleculeSurface);
    mp.optsurfmeshenable = 1;
    mp.optvolmeshenable  = 1;
    
    printMeshingOptions(mp, "Meshing options for molecular surface");
    
    cout << "Initialise the STL Geometry structure...." << endl;
    ngSurface = Ng_STL_InitSTLGeometry(stl_geom);
    if(ngSurface != NG_OK) {
      cout << "Error Initialising the STL Geometry....Aborting!!" << endl;
      exit(1);
    }
    
    cout << "Start Edge Meshing...." << endl;
    ngSurface = Ng_STL_MakeEdges(stl_geom, ngVolume, &mp);
    if(ngSurface != NG_OK) {
      cout << "Error in Edge Meshing....Aborting!!" << endl;
      exit(1);
    }
    
    cout << "Start Surface Meshing...." << endl;
    ngSurface = Ng_STL_GenerateSurfaceMesh(stl_geom, ngVolume, &mp);
    if(ngSurface != NG_OK) {
      cout << "Error in Surface Meshing....Aborting!!" << endl;
      exit(1);
    }
    
    if (ini.get_optional<string>("meshing.second_order_surface").get_value_or("no") == "yes"){
      cout << "Generating second order geometry mesh" << endl;
      Ng_STL_Generate_SecondOrder(stl_geom, ngVolume);
    }
	  
    if (debug == "analyze" && mode == "protein"){
      stringstream fileName;
      fileName << jobname << "_surface_so.vol";
      Ng_SaveMesh(ngVolume,fileName.str().c_str());
    }
    
    if (mode == "protein" && boost::filesystem::exists( "cavity.vol" )){
      string cavity = "cavity.vol";
      Ng_Mesh* cSurface;
      cSurface = nglib::Ng_LoadMesh(cavity.c_str());
      
      Ng_SetProperties(ngVolume, 1, 1, 1, 0);
      Ng_SetProperties(cSurface, 1, 1, 1, 0);
      
      
      cout << "Merging Mesh with cavity....." << endl;
      ngSurface = Ng_Alt_MergeMesh( ngVolume, cSurface );
      if(ngSurface != NG_OK) {
	cout << "Error in cavity merging....Aborting!!" << endl;
	exit(1);
      }
      
    }
    
    // Example: Set up a global refinement if using a local h file
    // float globalH = 1e6;
    // cout << "setting global H to " << globalH << endl;
    // Ng_RestrictMeshSizeGlobal(ngVolume, globalH);
    
    setMeshingOptions(mp, meshMoleculeVolume);
    mp.optsurfmeshenable = 1;
    mp.optvolmeshenable  = 1;
    
    printMeshingOptions(mp, "Meshing options for molecular volume");
    
    cout << "Start Volume Meshing of molecule...." << endl;
    ngSurface = Ng_GenerateVolumeMesh (ngVolume, &mp);
    if(ngSurface != NG_OK) {
      cout << "Error in Volume Meshing....Aborting!!" << endl;
      exit(1);
    }

    // Example: Saving a volume file (VOL, NETGEN format) in current state
    // Ng_SaveMesh(ngVolume,"meshed.vol");


    // Example: Generate a second order surface mesh
    // if (ini.get<string>("meshing.second_order_surface") == "yes")
    //    Ng_STL_Generate_SecondOrder (stl_geom, ngVolume);

    /// If the ion exclusion layer has to be computed (ionc != 0)
    if (exclusionSurface){
      
      string exclusionVol = "exclusion.vol";
      string exclusionStl = "exclusion.stl";
      
      if (fname != "")
	exclusionVol = fname+"_exclusion.vol";
      if (fname != "")
	exclusionStl = fname+"_exclusion.stl";
      
      Ng_Mesh *exclusionVolume;
      exclusionVolume = Ng_NewMesh();
      setMeshingOptions(mp, meshMoleculeSurface);
      mp.optsurfmeshenable = 1;
      mp.optvolmeshenable  = 0;
      
      if (!boost::filesystem::exists( exclusionVol )){
	Ng_STL_Geometry *exclusion_geom = Ng_STL_NewGeometry();
	
	cout << "Loading exclusion " <<  exclusionStl << endl;
	exclusion_geom = Ng_STL_LoadGeometry( exclusionStl.c_str());
	if(!exclusion_geom)
	  {
	    cout << "Error reading in STL File: " << exclusionStl << endl;
	    return 1;
	  }
	cout << "Successfully loaded STL File: " << exclusionStl << endl;
	
	cout << "Initialise the STL Geometry structure for ion exclusion layer...." << endl;
	ngSurface = Ng_STL_InitSTLGeometry(exclusion_geom);
	if(ngSurface != NG_OK) {
	  cout << "Error Initialising the STL Geometry for ion exclusion layer....Aborting!!" << endl;
	  exit(1);
	}
	
	cout << "Start Edge Meshing...." << endl;
	ngSurface = Ng_STL_MakeEdges(exclusion_geom, exclusionVolume, &mp);
	if(ngSurface != NG_OK) {
	  cout << "Error in Edge Meshing ion exclusion layer....Aborting!!" << endl;
	  exit(1);
	}
	
	cout << "Start Surface Meshing...." << endl;
	ngSurface = Ng_STL_GenerateSurfaceMesh(exclusion_geom, exclusionVolume, &mp);
	if(ngSurface != NG_OK) {
	  cout << "Error in Surface Meshing ion exclusion layer....Aborting!!" << endl;
	  exit(1);
	}
	Ng_SaveMesh(exclusionVolume,exclusionVol.c_str());
	
      } else {
	exclusionVolume = nglib::Ng_LoadMesh(exclusionVol.c_str());
      }
      
      Ng_SetProperties(ngVolume, 1, 1, 1, 0);
      Ng_SetProperties(exclusionVolume, 1, 1, 1, 0);
      
      cout << "Merging Mesh with ion exclusion layer....." << endl;
      ngSurface = Ng_Alt_MergeMesh( exclusionVolume, ngVolume );
      if(ngSurface != NG_OK) {
	cout << "Error in exclusion layer merging....Aborting!!" << endl;
	exit(1);
      }
      	    	    
      mp.optsurfmeshenable = 1;
      mp.optvolmeshenable  = 1;
      
      cout << "Start Volume meshing of exclusion layer...." << endl;
      ngSurface = Ng_GenerateVolumeMesh (exclusionVolume, &mp);
      if(ngSurface != NG_OK) {
	cout << "Error in Volume Meshing of exclusion layer....Aborting!!" << endl;
	exit(1);
      }
      cout << "Meshing successfully completed....!!" << endl;

      Ng_Mesh* bSurface;
      
      cout << "Loading boundary settings ....." << endl;
      string boundary = ini.get<string>("model.boundary");
      bSurface = nglib::Ng_LoadMesh(boundary.c_str());
      
      Ng_SetProperties(bSurface, 1, 1, 1, 0);
      
      cout << "Merging Mesh with boundary....." << endl;
      ngSurface = Ng_Alt_MergeMesh( bSurface, exclusionVolume );
      if(ngSurface != NG_OK) {
	cout << "Error in Surface merging with ion exclusion....Aborting!!" << endl;
	exit(1);
      }
      
      cout << "Merging complete" << endl;
      
      string meshBoundaryVolume = ini.get<string>("meshing.boundary_volume");
      setMeshingOptions(mp, meshBoundaryVolume);
      
      mp.optsurfmeshenable = 1;
      mp.optvolmeshenable  = 1;
      printMeshingOptions(mp, "Meshing options for boundary volume");
      
      
      cout << "Start Volume meshing of whole model...." << endl;
      ngSurface = Ng_GenerateVolumeMesh (bSurface, &mp);
      if(ngSurface != NG_OK) {
	cout << "Error in Volume Meshing....Aborting!!" << endl;
	exit(1);
      }
      
      // volume mesh output
      np = Ng_GetNP(bSurface);
      cout << "Points: " << np << endl;
      
      ne = Ng_GetNE(bSurface);
      cout << "Elements: " << ne << endl;
      
      string volumeVol = ini.get<string>("model.volume_vol");
      
      if (fname != "")
	volumeVol = fname;
      
      if (volumeVol != ""){
	cout << "Saving Mesh in VOL Format...." << endl;
	Ng_SaveMesh(bSurface,volumeVol.c_str());   	
      }
    } else {

      /// If the ion exclusion layer is not computed (ionc = 0)
      Ng_Mesh* bSurface;
      
      if (mode != "cavity"){
	cout << "Loading boundary settings ....." << endl;
	string boundary = ini.get<string>("model.boundary");
	bSurface = nglib::Ng_LoadMesh(boundary.c_str());
	
	if (ini.get_optional<string>("meshing.second_order_surface").get_value_or("no") == "yes"){
	  cout << "Generating second order for volume mesh" << endl;
	  Ng_Generate_SecondOrder(bSurface);
	}
	
	Ng_SetProperties(ngVolume, 1, 1, 1, 0);
	Ng_SetProperties(bSurface, 1, 1, 1, 0);
	
	cout << "Merging Mesh with boundary....." << endl;
	ngSurface = Ng_Alt_MergeMesh( bSurface, ngVolume );
	if(ngSurface != NG_OK) {
	  cout << "Error in Surface merging....Aborting!!" << endl;
	  exit(1);
	}
	      
	cout << "Merging complete" << endl;
	
	
	string refineFile = ini.get_optional<string>("model.refine_file").get_value_or("");

	/// It is possible to use a file with points and h parameter to have control
	/// over density manually (not recommended yet)
	if (refineFile != ""){
	  // global refinement
	  // void Ng_RestrictMeshSizeGlobal (Ng_Mesh * mesh, double h);
	  cout << "Setting local refinement ....." << endl;
	  ifstream in(refineFile.c_str());
	  if (!in) {
	    in.close();
	    cout << "Cannot open refinement file:" << refineFile << endl;
	    exit(0);
	  }
	  
	  string currentLine;
	  double p[3]; double h;
	  unsigned int refPoints = 0;
	  while( !in.eof() ) {
	    getline(in, currentLine);
	    if (currentLine != ""){
	      istringstream ss(currentLine);
	      ss >> p[0] >> p[1] >> p[2] >> h;
	      Ng_RestrictMeshSizePoint (bSurface, p, h);
	      refPoints++;
	    }
	  }
	  cout << refPoints << " local refinement point(s) set." << endl;
	}

	/// Loading the boundary layer
	string meshBoundaryVolume = ini.get<string>("meshing.boundary_volume");
	setMeshingOptions(mp, meshBoundaryVolume);
	mp.optsurfmeshenable = 1;
	mp.optvolmeshenable  = 1;
	printMeshingOptions(mp, "Meshing options for boundary volume");
	

	/// Meshing the whole model
	cout << "Start Volume meshing of whole model...." << endl;
	ngSurface = Ng_GenerateVolumeMesh (bSurface, &mp);
	if(ngSurface != NG_OK) {
	  cout << "Error in Volume Meshing....Aborting!!" << endl;
	  exit(1);
	}
      } else {
	bSurface = ngVolume;
	Ng_SetProperties(ngVolume, 2, 2, 2, 1);
	
      }
      
      
      if (ini.get_optional<string>("meshing.second_order_surface").get_value_or("no") == "yes"){
	cout << "Last second order meshing" << endl;
	Ng_Generate_SecondOrder(bSurface);
      }
	  
      cout << "Meshing successfully completed....!!" << endl;
	  
      // volume mesh output
      np = Ng_GetNP(bSurface);
      cout << "Points: " << np << endl;
      
      ne = Ng_GetNE(bSurface);
      cout << "Elements: " << ne << endl;
      
      string volumeVol = ini.get<string>("model.volume_vol");
      
      /// Volume mesh is saved
      if (fname != "")
	volumeVol = fname;

      if (volumeVol != ""){
	cout << "Saving Mesh in VOL Format...." << endl;
	Ng_SaveMesh(bSurface,volumeVol.c_str());
      }
    }
    
    return 1;
    
  }


  /** @brief This method performs some cleaning of the mesh. Here, it removes
   *         duplicate vertices and removes unreferenced vertices (if any available)
   *
   *  Here, the surface may be cleaned or modified. Because we want to leave the 
   *  molecular surface invariant just a removing of duplicate vertices and a
   *  removing of unreferenced vertices is performed using methods from the VCG lib.
   *
   *  @param mMesh A VCG lib surface structure.
   *  @return Void.
   *
   */    
  void clean(mMesh &mSurface){
    // some cleaning to get rid of bad file formats like stl that duplicate vertexes..
    int dup = tri::Clean<mMesh>::RemoveDuplicateVertex(mSurface);
    int unref = tri::Clean<mMesh>::RemoveUnreferencedVertex(mSurface);
    printf("Removed %i duplicate and %i unreferenced vertices from mesh\n",dup,unref);
    tri::UpdateTopology<mMesh>::VertexFace(mSurface);
  }
  
  /** @brief A smoothing of the molecular surface is performed.
   *
   *  The aim of this method is to smooth the molecular surface, leaving the 
   *  surface invariant and without injecting self intersections. This is possible
   *  using the Taubin smoothing algorithm which may be performed many times.
   *  This step has to be done to boost up the advancing front algorithm implemented
   *  with NETGEN and to make the molecular surface meshing more robust.
   *  Other smoothing algorithms are possible but did not show better performance:
   *  Laplace, HC Laplace, Angle weigthed Laplace, Scale dependent Laplace.
   *
   *  @param mMesh A VCG lib surface structure.
   *  @param INI Options file to declare which smoothing to perform with how many
   *             iterations.
   *  @return Void.
   *
   */    
  void smooth(mMesh &mSurface, INI &ini){
    
    string line = ini.get<string>("model.smoothing");
    istringstream ss(line);
    string mode = "";
    unsigned int steps = 0;
    string progress = "*";
    
    while (ss >> mode >> steps){
      if (mode == "t"){
	cout << "taubin smoothing steps: " << steps << endl;
	for (unsigned int i = 0; i < steps; i++){
	  cout << "\r" << i+1 << "/" << steps << flush;
	  tri::UpdateNormal<mMesh>::PerFace(mSurface);
	  tri::Smooth<mMesh>::VertexCoordTaubin(mSurface,1,TAUBIN_LAMBDA,TAUBIN_MU);
	}
      } else if (mode == "lap"){
	cout << "laplace smoothing steps: " << steps << endl;
	for (unsigned int i = 0; i < steps; i++){
	  cout << "\r" << i+1 << "/" << steps << flush;
	  tri::UpdateNormal<mMesh>::PerFace(mSurface);
	  tri::Smooth<mMesh>::VertexCoordLaplacian(mSurface,1);
	}
      } else if (mode == "hc"){
	cout << "hc laplace smoothing steps: " << steps << endl;
	for (unsigned int i = 0; i < steps; i++){
	  cout << "\r" << i+1 << "/" << steps << flush;
	  tri::UpdateNormal<mMesh>::PerFace(mSurface);
	  tri::Smooth<mMesh>::VertexCoordLaplacianHC(mSurface,1);
	}
      } else if (mode == "aw"){
	cout << "angle weighted laplace smoothing steps: " << steps << endl;
	for (unsigned int i = 0; i < steps; i++){
	  cout << "\r" << i+1 << "/" << steps << flush;
	  tri::UpdateNormal<mMesh>::PerFace(mSurface);
	  tri::Smooth<mMesh>::VertexCoordLaplacianAngleWeighted(mSurface,1, CF_LAMBDA);
	}
      } else if (mode == "lapsd"){
	cout << "scale dependent laplace smoothing steps: " << steps << endl;
	for (unsigned int i = 0; i < steps; i++){
	  cout << "\r" << i+1 << "/" << steps << flush;
	  tri::UpdateNormal<mMesh>::PerFace(mSurface);
	  // tri::UpdateNormal<mMesh>::PerFaceNormalized(mSurface);
	  tri::Smooth<mMesh>::VertexCoordScaleDependentLaplacian_Fujiwara(mSurface,1, 0.0025);
	}	
      }
      cout << endl;
			
    }
    
    // Example code to test for self intersections, Abort here if true?
    //		vector<mMesh::FaceType *> SelfIntersectList;
    // tri::Clean<mMesh>::SelfIntersections(mSurface, SelfIntersectList);
    //	    int SelfIntersections = SelfIntersectList.size();
    //	cout << SelfIntersections << " intersections found." << endl;
    //		if (SelfIntersections > 0)
    //	cout << "This surface will probably not mesh! Keep care." << endl;
    
  }
  
  Ng_Result ngSurface;
  Ng_Mesh *ngVolume;
  
};

bool Model::exclusionSurface = false;
