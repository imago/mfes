#include <iostream>
#include <fstream>
#include <string>

#include "LSMS.h"
#include "VCG.h"

using namespace std;

typedef boost::property_tree::ptree INI;

class Model {

public:
	void calcSurface(vector<Atom> &atomList, INI &ini){
		mMesh mSurface;
	    LSMS lmSurface;
	    lmSurface.calcMC(mSurface, atomList, ini);

	    clean(mSurface);
	    smooth(mSurface, ini);

	    tri::io::ExporterSTL<mMesh>::Save(mSurface,"output.stl",false);
	    /*
       // convert moleculeSurface to NETGEN surface
	   …. netgenSurface
	   boundary = createBoundary(radius)
	   mergedSurface = merge(boundary, netgenSurface)
	   surface = mergedSurface
*/
	}


	void calcVolume(){
//	    loadOptions()
//	   NG_VOL newVolume = calcVolume
//	   volume = newVolume
	}

	void saveSTL(string fileName){
//	  surface.saveSTL(fileName+“.stl“)<- with netgen
	}

	void saveVOL(string fileName){
//	  volume.saveVOL(fileName+“.vol“) <- with netgen
	}

private:

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

	};

//  NG_STL surface;
//  NG_VOL volume;

};
