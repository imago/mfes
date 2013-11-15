#ifndef VCG_LIB

#define VCG_LIB

// VCG library
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/normal.h>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/smooth.h>

// input output
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_stl.h>


using namespace vcg;


class mFace;
class mVertex;

// Meshdeclaration
struct mUsedTypes : public UsedTypes<	Use<mVertex>::AsVertexType, Use<mFace>::AsFaceType>{};
class mVertex  : public Vertex< mUsedTypes, vertex::VFAdj, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
//class mVertex  : public Vertex< mUsedTypes, vertex::Coord3f, vertex::BitFlags, vertex::Normal3f, vertex::Mark,vertex::Color4b, vertex::Qualityf,vertex::VFAdj, vertex::Curvaturef >{};

class mFace    : public Face  < mUsedTypes, face::VFAdj, face::Normal3f, face::VertexRef, face::Mark, face::BitFlags, face::VertexRef > {};
//class mFace    : public Face  < mUsedTypes, face::VertexRef,face::BitFlags,face::Mark, face::Normal3f, vcg::face::VFAdj, vcg::face::FFAdj>{};

class mMesh    : public vcg::tri::TriMesh<vector<mVertex>, vector<mFace> > {};

#include <vcg/space/point3.h>
#include <vcg/complex/algorithms/clean.h>

namespace vcg
{
template<class MeshType>
std::pair<int,int>  RemoveSmallConnectedComponentsSize(MeshType &m, int maxCCSize)
{
  std::vector< std::pair<int, typename MeshType::FacePointer> > CCV;
      int TotalCC=vcg::tri::Clean<MeshType>::ConnectedComponents(m, CCV);
      cout << "# connected components: " << TotalCC << endl;
			int DeletedCC=0;

      tri::ConnectedIterator<MeshType> ci;
      for(unsigned int i=0;i<CCV.size();++i)
      {
        std::vector<typename MeshType::FacePointer> FPV;
        if(CCV[i].first<maxCCSize)
        {
					DeletedCC++;
          for(ci.start(m,CCV[i].second);!ci.completed();++ci)
            FPV.push_back(*ci);

          typename std::vector<typename MeshType::FacePointer>::iterator fpvi;
          for(fpvi=FPV.begin(); fpvi!=FPV.end(); ++fpvi)
						tri::Allocator<MeshType>::DeleteFace(m,(**fpvi));
        }
      }
			return std::pair<int,int>(TotalCC,DeletedCC);
}
}



namespace vcg {
namespace tri {

template<class _MeshType>
class SmallComponent
{

public:
	typedef _MeshType MeshType;
	typedef typename MeshType::VertexType     VertexType;
	typedef typename MeshType::VertexPointer  VertexPointer;
	typedef typename MeshType::VertexIterator VertexIterator;
	typedef typename MeshType::FaceType       FaceType;
	typedef typename MeshType::FacePointer    FacePointer;
	typedef typename MeshType::FaceIterator   FaceIterator;

	static int Select(MeshType &m, float nbFaceRatio = 0.1, bool nonClosedOnly = false)
	{
		assert(tri::HasFFAdjacency(m) && "The small component selection procedure requires face to face adjacency.");

		// the different components as a list of face pointer
		std::vector< std::vector<FacePointer> > components;

		for(uint faceSeed = 0; faceSeed<m.face.size(); )
		{
			// find the first not selected face
			bool foundSeed = false;
			while (faceSeed<m.face.size())
			{
				mFace& f = m.face[faceSeed];
				if (!f.IsS())
				{
					if (nonClosedOnly)
					{
						for (int k=0; k<3; ++k)
							if (face::IsBorder(f,k))
							{
								foundSeed = true;
								break;
							}
					}
					else
						foundSeed = true;
					if (foundSeed)
						break;
				}
				++faceSeed;
			}
			if (!foundSeed) // no more seed, stop
				break;

			// expand the region from this face...
			components.resize(components.size()+1);
			std::vector<FacePointer> activefaces;

			activefaces.push_back(&m.face[faceSeed]);

			while (!activefaces.empty())
			{
				FacePointer f = activefaces.back();
				activefaces.pop_back();
				if (f->IsS())
					continue;

				f->SetS();
				components.back().push_back(f);

				for (int k=0; k<3; ++k)
				{
					if (face::IsBorder(*f,k))
						continue;

					FacePointer of = f->FFp(k);
					if (!of->IsS()){
						activefaces.push_back(of);
					}
				}
			}
			++faceSeed;
		}
		cout << "update selection" << endl;
		UpdateSelection<MeshType>::FaceClear(m);

		// now the segmentation is done, let's compute the absolute face count threshold
		int total_selected = 0;
		int maxComponent = 0;
		for (uint i=0; i<components.size(); ++i)
		{
//			cout << "Component " << i << " -> " << components[i].size() << "\n";
			total_selected += components[i].size();
			maxComponent = std::max<int>(maxComponent,components[i].size());
		}
		int remaining = m.face.size() - total_selected;
		uint th = std::max(maxComponent,remaining) * nbFaceRatio;

		int selCount = 0;
		for (uint i=0; i<components.size(); ++i)
		{
			if (components[i].size()<th)
			{
				selCount += components[i].size();
				for (uint j=0; j<components[i].size(); ++j)
					components[i][j]->SetS();
			}
		}
		cout << selCount << endl;
		return selCount;
	}

	static void DeleteFaceVert(MeshType &m)
	{
		typename MeshType::FaceIterator fi;
		typename MeshType::VertexIterator vi;
		UpdateSelection<MeshType>::VertexClear(m);
		UpdateSelection<MeshType>::VertexFromFaceStrict(m);

		for(fi=m.face.begin();fi!=m.face.end();++fi)
			if (!(*fi).IsD() && (*fi).IsS() )
				tri::Allocator<mMesh>::DeleteFace(m,*fi);
		for(vi=m.vert.begin();vi!=m.vert.end();++vi)
			if (!(*vi).IsD() && (*vi).IsS() )
				tri::Allocator<mMesh>::DeleteVertex(m,*vi);
	}

}; // end class

}	// End namespace
}	// End namespace

#endif
