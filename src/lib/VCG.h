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
class mFace    : public Face  < mUsedTypes, face::VFAdj, face::Normal3f, face::VertexRef, face::Mark, face::BitFlags > {};
class mMesh    : public vcg::tri::TriMesh<vector<mVertex>, vector<mFace> > {};


#endif
