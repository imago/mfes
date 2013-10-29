//#include <solve.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "LSMS.h"
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
	void calcModel(vector<Atom> &atomList, INI &ini, string fname = ""){
	    LSMS lmSurface;
	    mMesh mSurface;

	    string volume_vol = ini.get<string>("model.volume_vol");
    	string surface_stl = ini.get<string>("model.surface_stl");

	    if ( fname == "" && !boost::filesystem::exists( volume_vol ) ){
	    // protein model is calculating
	    	lmSurface.calcMC(mSurface, atomList, ini);

	    	clean(mSurface);
	    	smooth(mSurface, ini);

	    	if (surface_stl != "" )
	    		tri::io::ExporterSTL<mMesh>::Save(mSurface,surface_stl.c_str(),false);

	    	convert(mSurface, ini);
	    } else if ( fname != "" ){
	    	// group model with id nr is calculating
	    	if ( !boost::filesystem::exists( fname ) ){
	    		    lmSurface.calcMC(mSurface, atomList, ini);

	    		    clean(mSurface);
	    		    smooth(mSurface, ini);

	    		    convert(mSurface, ini, fname);
	    	}
	    }
	}



private:

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

	void convert(mMesh &mSurface, INI &ini, string fname = ""){

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
		string meshOptions = ini.get<string>("model.options");
		Ng_Meshing_Parameters mp;
		setMeshingOptions(mp, meshOptions);

		mp.optsurfmeshenable = 1;
		mp.optvolmeshenable  = 1;


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

		cout << "Loading boundary settings ....." << endl;
		string boundary = ini.get<string>("model.boundary");
		Ng_Mesh* bSurface;
		bSurface = nglib::Ng_LoadMesh(boundary.c_str());

		Ng_SetBCProp(bSurface, 1);
		Ng_SetBCProp(ngVolume, 2);

		cout << "Merging Mesh with boundary....." << endl;
     	ngSurface = Ng_MergeMesh( bSurface, ngVolume );
		if(ngSurface != NG_OK) {
			cout << "Error in Surface merging....Aborting!!" << endl;
			exit(1);
		}

		mp.maxh = 1e6;

		cout << "Start Volume Meshing...." << endl;
		ngSurface = Ng_GenerateVolumeMesh (bSurface, &mp);
		if(ngSurface != NG_OK) {
			cout << "Error in Volume Meshing....Aborting!!" << endl;
		    exit(1);
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
				for (int i = 0; i < steps; i++){
					tri::UpdateNormal<mMesh>::PerFaceNormalized(mSurface);
					tri::Smooth<mMesh>::VertexCoordTaubin(mSurface,1,TAUBIN_LAMBDA,TAUBIN_MU);
				}
			} else if (mode == "lap"){
				for (int i = 0; i < steps; i++){
					tri::UpdateNormal<mMesh>::PerFaceNormalized(mSurface);
					tri::Smooth<mMesh>::VertexCoordLaplacian(mSurface,1);
				}
			} else if (mode == "hc"){
				for (int i = 0; i < steps; i++){
					tri::UpdateNormal<mMesh>::PerFaceNormalized(mSurface);
					tri::Smooth<mMesh>::VertexCoordLaplacianHC(mSurface,1);
				}
			} else if (mode == "aw"){
				for (int i = 0; i < steps; i++){
					tri::UpdateNormal<mMesh>::PerFaceNormalized(mSurface);
					tri::Smooth<mMesh>::VertexCoordLaplacianAngleWeighted(mSurface,1, CF_LAMBDA);
				}
			}
		}
		vector<mMesh::FaceType *> SelfIntersectList;
		tri::Clean<mMesh>::SelfIntersections(mSurface, SelfIntersectList);
	    int SelfIntersections = SelfIntersectList.size();
		cout << SelfIntersections << " intersections found." << endl;
		if (SelfIntersections > 0)
			cout << "This surface will probably not mesh! Keep care." << endl;

	}

	Ng_Result ngSurface;
	Ng_Mesh *ngVolume;

};
