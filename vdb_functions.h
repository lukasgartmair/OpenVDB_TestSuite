#ifndef VDB_FUNCTIONS_H
#define VDB_FUNCTIONS_H

#include <openvdb/openvdb.h>
#include <openvdb/Grid.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/math/Coord.h>

typedef struct coord_main {
float x;
float y;
float z;
} coord;

namespace VDB{

// pure vdb functions
openvdb::FloatGrid::Ptr loadData();
openvdb::FloatGrid::Ptr createBlock(float radius, float value);
std::vector<openvdb::Vec3s> volumeToMeshVertices(openvdb::FloatGrid::Ptr grid, float isovalue, float adaptivity);

// helper functions to process the vdb code
int roundUp(float numToRound, float multiple);
coord getVoxelIndex(coord *vec, float voxsize);

// conversion functions for lpcvt
std::vector<std::vector<unsigned int> > splitQuadsToTriangles(std::vector<openvdb::Vec3s> points, std::vector<openvdb::Vec4I> quads);
std::vector<std::vector<unsigned int> > concatenateTriangleVectors(std::vector<openvdb::Vec3I> triangles, std::vector<std::vector<unsigned int> > triangles_from_splitted_quads);
std::vector<std::vector<float> > convertOpenVDBVectorToStandardVector(std::vector<openvdb::Vec3s> points);
std::vector<std::vector<unsigned int> > increaseTriangleVertexIndicesByN(std::vector<std::vector<unsigned int> > triangles, int N);
std::vector<std::vector<unsigned int> > decreaseTriangleVertexIndicesByN(std::vector<std::vector<unsigned int> > triangles, int N);

// vector stuff for normal calculation
float getLengthOfVector(std::vector<float> vec);
std::vector<float> normalizeVector(std::vector<float> vec);
std::vector<float> getCrossProduct(std::vector<float> vec1, std::vector<float> vec2);
std::vector<std::vector<float> > computeTriangleNormals(std::vector<std::vector<float> > points, std::vector<std::vector<unsigned int> > triangles);
std::vector<std::vector<float> > computeTriangleNormalsVDB(std::vector<openvdb::Vec3s> points, std::vector<std::vector<unsigned int> > triangles);
std::vector<std::vector<float> > computeVertexNormals(std::vector<std::vector<unsigned int> > triangles, std::vector<openvdb::Vec3s> points, std::vector<std::vector<float> > triangle_normals);

std::vector<float> computeTriangleAreas(std::vector<openvdb::Vec3s> points, std::vector<std::vector<unsigned int> > triangles);

// IO
void exportVDBMeshAsObj(std::vector<openvdb::Vec3s> points, std::vector<openvdb::Vec3I> triangles, std::vector<openvdb::Vec4I> quads);
void exportMeshAsVDB(openvdb::FloatGrid::Ptr grid);
void exportTriangleAreas(std::vector<unsigned int> triangle_areas);
void exportTriangleMeshAsObj(std::vector<openvdb::Vec3s> points, std::vector<std::vector<unsigned int> > triangles);

}

#endif
