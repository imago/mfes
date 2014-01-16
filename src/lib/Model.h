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


namespace nglib {
	#include <nglib.h>
}
using namespace nglib;


typedef boost::property_tree::ptree INI;

class Model {

public:
	void calcModel(vector<Atom> &atomList, INI &ini, string fname = "", bool cavity = false){

	  LSMS lmSurface;
	  mMesh mSurface;
	  
	  string volume_vol = ini.get<string>("model.volume_vol");
	  string surface_stl = ini.get<string>("model.surface_stl");
	  string generator = ini.get<string>("model.generator");
	  string generatorResidue = ini.get<string>("model.generator_residue");
	  int generatorResolution = atoi(ini.get<string>("model.grid_resolution").c_str());
	  int generatorModelResolution = atoi(ini.get<string>("model.grid_residue_resolution").c_str());


	  int mask = 1;

	  if (generator == "standard" && cavity && !boost::filesystem::exists( "cavity.vol" )){
	    // cavity is calculated
	    cout << "Calculating cavities..." << endl;
	    boost::timer t;
	    ofstream time;
	    time.open ("times");
	    
	    if (lmSurface.calcMC(mSurface, atomList, ini, "cavity")){
	      
	      clean(mSurface);
	      smooth(mSurface, ini);
	      
	      if (surface_stl != "" )
		tri::io::ExporterSTL<mMesh>::Save(mSurface,"cavity.stl",false);

	      double currentTime = t.elapsed();
	      time << "cavity_surface " << currentTime << " s" << endl;
	      cout << "mFES: cavity model surface calculation took " << currentTime << " seconds." <<endl;
	      
	      t.restart();
	      convert(mSurface, ini, "cavity.vol", "cavity");
	      currentTime = t.elapsed();
	      time << "cavity_volume " << currentTime << " s" << endl;
	      cout << "mFES: cavity model VOL meshing took " << currentTime << " seconds." <<endl;
	    }
	    time.close();
	    
	  } else if ( fname == "" && !boost::filesystem::exists( volume_vol ) ){
	    // protein model is calculated
	    boost::timer t;
	    ofstream time;
	    time.open ("times");
	    
	    if (boost::filesystem::exists( surface_stl )) {
	      cout << surface_stl << " found. " << endl;
	      tri::io::ImporterSTL<mMesh>::Open(mSurface, surface_stl.c_str(), mask);
	      // clean(mSurface);
	      // smooth(mSurface, ini);
	      
	    } else {
	      if (generator == "standard"){
		lmSurface.calcMC(mSurface, atomList, ini);
		clean(mSurface);
		smooth(mSurface, ini);
		if (surface_stl != "" )
		  tri::io::ExporterSTL<mMesh>::Save(mSurface,surface_stl.c_str(),false);
	      } else {
		Voxel vSurface;
		if (!boost::filesystem::exists( "protein.stl" )){
		  vSurface.calcSurface(mSurface, atomList, ini, "protein.stl", generatorResolution);
		  clean(mSurface);
		} 
		tri::io::ImporterSTL<mMesh>::Open(mSurface, string("protein.stl").c_str(), mask);			     
	      }
	      
	    }
	    
//			MeshModel* mesh = md.mm();
//			MeshModel::MM_FACEFACETOPO | MeshModel::MM_FACEFLAGBORDER);
//			tri::UpdateTopology<mMesh>::FaceFace(mSurface);
//			tri::UpdateTopology<mMesh>::VertexFace(mSurface);
//
//			tri::SmallComponent<mMesh>::Select(mSurface);
//			tri::SmallComponent<mMesh>::DeleteFaceVert(mSurface);
//
//			tri::UpdateTopology<mMesh>::FaceFace(mSurface);
//			tri::UpdateTopology<mMesh>::VertexFace(mSurface);
//
//
//			cout << "cc: " << tri::Clean<mMesh>::CountHoles(mSurface) << endl;
//			RemoveSmallConnectedComponentsSize(mSurface, 1);
//
//
//    		tri::io::ExporterSTL<mMesh>::Save(mSurface,"output.stl",false);
//    		exit(0);


	    double currentTime = t.elapsed();
	    time << "protein_surface " << currentTime << " s" << endl;
	    cout << "mFES: protein model surface calculation took " << currentTime << " seconds." <<endl;
	    
	    t.restart();
	    convert(mSurface, ini);
	    currentTime = t.elapsed();
	    time << "protein_volume " << currentTime << " s" << endl;
	    cout << "mFES: protein model VOL meshing took " << currentTime << " seconds." <<endl;
	    time.close();
	  }
	  else if (fname != ""){
	    // group model with id nr is calculated
	    if ( !boost::filesystem::exists( fname ) ){
	      string stlFile = fname+".stl";
	      ofstream time;
	      time.open ("times", ios::app);
	      boost::timer t;
	      if ( !boost::filesystem::exists( stlFile )){
		
		if (generatorResidue == "standard"){
		  lmSurface.calcMC(mSurface, atomList, ini);
		  clean(mSurface);
		  smooth(mSurface, ini);
		  tri::io::ExporterSTL<mMesh>::Save(mSurface,stlFile.c_str(),false);
		} else {
		  Voxel vSurface;
		  vSurface.calcSurface(mSurface, atomList, ini, stlFile, generatorModelResolution);
		  clean(mSurface);
		  tri::io::ImporterSTL<mMesh>::Open(mSurface, stlFile.c_str(), mask);
		}
		
		double currentTime = t.elapsed();
		time << "model_surface " << currentTime << " s" << endl;
		cout << "mFES: residue model surface calculation took " << currentTime << " seconds." <<endl;
		
		t.restart();
		convert(mSurface, ini, fname, "residue");
		currentTime = t.elapsed();
		time << "model_volume " << currentTime << " s" << endl;
		cout << "mFES: residue model VOL meshing took " << currentTime << " seconds." <<endl;
		time.close();
	      }
	      else {
		cout << "model surface found: " << stlFile << endl;
		tri::io::ImporterSTL<mMesh>::Open(mSurface, stlFile.c_str(), mask);
		t.restart();
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

	void printMeshingOptions(Ng_Meshing_Parameters &mp, string prefix){
		cout << endl;
		cout << prefix << endl;
		cout << "========================" << endl;
		cout << "options.localh = " << mp.uselocalh << endl;
		cout << "options.meshsize = " << mp.maxh << endl;
		cout << "options.minmeshsize = " << mp.minh << endl;
		cout << "meshoptions.fineness = " << mp.fineness << endl;
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
	void setMeshingOptions(Ng_Meshing_Parameters &mp, string optFile){
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
			ss >> variable >> value;

			if (variable == "options.localh"){
				//!< Switch to enable / disable usage of local mesh size modifiers
				mp.uselocalh = atoi(value.c_str());
				continue;
			}

			if (variable == "options.meshsize"){
				//!< Maximum global mesh size allowed
			        mp.maxh = atof(value.c_str());
			        continue;
			}

			if (variable == "options.minmeshsize"){
				//!< Minimum global mesh size allowed
				mp.minh = atof(value.c_str());
				continue;
			}

			if (variable == "meshoptions.fineness"){
				//!< Mesh density: 0...1 (0 => coarse; 1 => fine)
				mp.fineness = atof(value.c_str());
				continue;
			}

			if (variable == "options.grading"){
				 //!< Mesh grading: 0...1 (0 => uniform mesh; 1 => aggressive local grading)
				mp.grading = atof(value.c_str());
				continue;
			}


			if (variable == "options.curvaturesafety"){
				//!< Elements to generate per curvature radius
				mp.elementspercurve = atoi(value.c_str());
				continue;
			}

			if (variable == "options.segmentsperedge"){
				//!< Number of elements to generate per edge of the geometry
				mp.elementsperedge = atof(value.c_str());
				continue;
			}

			if (variable == "options.secondorder"){
				//!< Generate second-order surface and volume elements
				mp.second_order = atoi(value.c_str());
				continue;
			}

			if (variable == "options.quad"){
				//!< Creates a Quad-dominated mesh
				mp.quad_dominated = atoi(value.c_str());
				continue;
			}

//			if (variable == "options.meshsizefilename"){
//				//!< Optional external mesh size file
//				mp.meshsize_filename = value.c_str();
//				continue;
//			}

			if (variable == "stloptions.resthcloseedgeenable"){
				//!< Enable / Disable mesh refinement at close edges
				mp.closeedgeenable = atoi(value.c_str());
				continue;
			}

			if (variable == "stloptions.resthcloseedgefac"){
				//!< Factor to use for refinement at close edges (larger => finer)
				mp.closeedgefact = atof(value.c_str());
				continue;
			}

			if (variable == "options.optsteps2d"){
				//!< Number of optimize steps to use for 2-D mesh optimization
				mp.optsteps_2d = atoi(value.c_str());
				continue;
			}

			if (variable == "options.optsteps3d"){
				//!< Number of optimize steps to use for 3-D mesh optimization
				mp.optsteps_3d = atoi(value.c_str());
				continue;
			}

			if (variable == "options.inverttets"){
				//!< Invert all the volume elements
				mp.invert_tets = atoi(value.c_str());
				continue;
			}

			if (variable == "options.inverttrigs"){
				//!< Invert all the surface triangle elements
				mp.invert_trigs = atoi(value.c_str());
				continue;
			}

			if (variable == "options.checkoverlap"){
				//!< Check for overlapping surfaces during Surface meshing
				mp.check_overlap = atoi(value.c_str());
				//!< Check for overlapping surface elements before volume meshing
				mp.check_overlapping_boundary = atoi(value.c_str());
			}


		}			// optimize by default 1, because if number of
		// optsteps is 0, then no optimization made

		in.close();
	}

	int convert(mMesh &mSurface, INI &ini, string fname = "", string mode = "protein"){

//		using namespace nglib;

		Ng_STL_Geometry *stl_geom = Ng_STL_NewGeometry();

		// Initialise the Netgen Core library
		Ng_Init();

		// Actually create the mesh structure
		ngVolume = Ng_NewMesh();

		int np, ne;

		// Convert VCG MyMesh to Ng_STL_Geometry

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

		// Set the Meshing Parameters to be used
		string debug = ini.get<string>("model.debug");
		string meshMoleculeSurface, meshMoleculeVolume;
		if (mode == "protein" || mode == "cavity"){
			meshMoleculeSurface = ini.get<string>("meshing.molecule_surface");
			meshMoleculeVolume  = ini.get<string>("meshing.molecule_volume");
		} else if (mode == "residue") {
			meshMoleculeSurface = ini.get<string>("meshing.residue_surface");
			meshMoleculeVolume  = ini.get<string>("meshing.residue_volume");
		}
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



		if (debug == "analyze"){
			Ng_SaveMesh(ngVolume,"proteinSurface.vol");
		}

		if (mode == "protein" && boost::filesystem::exists( "cavity.vol" )){
			string cavity = "cavity.vol";
			Ng_Mesh* cSurface;
			cSurface = nglib::Ng_LoadMesh(cavity.c_str());

			Ng_SetProperties(ngVolume, 1, 1, 1, 0);
			Ng_SetProperties(cSurface, 1, 1, 1, 0);


			cout << "Merging Mesh with cavity....." << endl;
	     	ngSurface = Ng_MergeMesh( ngVolume, cSurface );
			if(ngSurface != NG_OK) {
				cout << "Error in cavity merging....Aborting!!" << endl;
				exit(1);
			}

		}
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
//		Ng_SaveMesh(ngVolume,"meshed.vol");

		if (ini.get<string>("meshing.second_order_surface") == "yes")
			Ng_STL_Generate_SecondOrder (stl_geom, ngVolume);



		if (debug == "analyze"){
			Ng_SaveMesh(ngVolume,"proteinVolume.vol");
		}

		Ng_Mesh* bSurface;

		if (mode != "cavity"){
			cout << "Loading boundary settings ....." << endl;
			string boundary = ini.get<string>("model.boundary");
			bSurface = nglib::Ng_LoadMesh(boundary.c_str());

			Ng_SetProperties(ngVolume, 1, 1, 1, 0);
			Ng_SetProperties(bSurface, 1, 1, 1, 0);

			cout << "Merging Mesh with boundary....." << endl;
			ngSurface = Ng_MergeMesh( bSurface, ngVolume );
			if(ngSurface != NG_OK) {
				cout << "Error in Surface merging....Aborting!!" << endl;
				exit(1);
			}


			string refineFile = ini.get<string>("model.refine_file");
			if (refineFile != ""){
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
					// global refinement
					// void Ng_RestrictMeshSizeGlobal (Ng_Mesh * mesh, double h);
				}
				cout << refPoints << " local refinement point(s) set." << endl;
			}

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
		} else {
				bSurface = ngVolume;
				Ng_SetProperties(ngVolume, 2, 2, 2, 1);

		}

		cout << "Meshing successfully completed....!!" << endl;

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

		return 1;

	}

	void clean(mMesh &mSurface){
		// some cleaning to get rid of bad file formats like stl that duplicate vertexes..
	    int dup = tri::Clean<mMesh>::RemoveDuplicateVertex(mSurface);
	    int unref = tri::Clean<mMesh>::RemoveUnreferencedVertex(mSurface);
	    printf("Removed %i duplicate and %i unreferenced vertices from mesh\n",dup,unref);
	    tri::UpdateTopology<mMesh>::VertexFace(mSurface);
	}

	void smooth(mMesh &mSurface, INI &ini){

		string line = ini.get<string>("model.smoothing");
		istringstream ss(line);
		string mode = "";
		unsigned int steps = 0;

		while (ss >> mode >> steps){
			if (mode == "t"){
			        cout << "taubin smoothing steps: " << steps << endl;
				for (unsigned int i = 0; i < steps; i++){
				        cout << "*";
					tri::UpdateNormal<mMesh>::PerFace(mSurface);
					tri::Smooth<mMesh>::VertexCoordTaubin(mSurface,1,TAUBIN_LAMBDA,TAUBIN_MU);
				}
			} else if (mode == "lap"){
                       	        cout << "laplace smoothing steps: " << steps << endl;
				for (unsigned int i = 0; i < steps; i++){
				        cout << "*";
				        tri::UpdateNormal<mMesh>::PerFace(mSurface);
					tri::Smooth<mMesh>::VertexCoordLaplacian(mSurface,1);
				}
			} else if (mode == "hc"){
  			        cout << "hc laplace smoothing steps: " << steps << endl;
				for (unsigned int i = 0; i < steps; i++){
                                        cout << "*";
					tri::UpdateNormal<mMesh>::PerFace(mSurface);
					tri::Smooth<mMesh>::VertexCoordLaplacianHC(mSurface,1);
				}
			} else if (mode == "aw"){
			        cout << "angle weighted laplace smoothing steps: " << steps << endl;
				for (unsigned int i = 0; i < steps; i++){
                                        cout << "*";
					tri::UpdateNormal<mMesh>::PerFace(mSurface);
					tri::Smooth<mMesh>::VertexCoordLaplacianAngleWeighted(mSurface,1, CF_LAMBDA);
				}
			}	else if (mode == "lapsd"){
			        cout << "scale dependent laplace smoothing steps: " << steps << endl;
				for (unsigned int i = 0; i < steps; i++){
				        cout << "*";
					tri::UpdateNormal<mMesh>::PerFace(mSurface);
//					tri::UpdateNormal<mMesh>::PerFaceNormalized(mSurface);
					tri::Smooth<mMesh>::VertexCoordScaleDependentLaplacian_Fujiwara(mSurface,1, 0.0025);
				}

			}
			cout << endl;

		}
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
