#ifndef VDBTEST1_H
#define VDBTEST1_H

#include <openvdb/openvdb.h>
//#include <openvdb/math/Maps.cc>
#include <openvdb/Grid.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <stdio.h>
#include <stdlib.h>
#include <openvdb/tools/Composite.h>

typedef struct coord_main {
float x;
float y;
float z;
} coord;


// pure vdb functions
openvdb::FloatGrid::Ptr loadData();
openvdb::FloatGrid::Ptr createBlock(float radius, float value);
std::vector<openvdb::Vec3s> volumeToMeshVertices(openvdb::FloatGrid::Ptr grid, double isovalue, double adaptivity);

// helper functions to process the vdb code
int roundUp(float numToRound, float multiple);
coord GetVoxelIndex(coord *vec, float voxsize);

// conversion functions for lpcvt
std::vector<std::vector<float> > splitQuadsToTriangles(std::vector<openvdb::Vec3s> points, std::vector<openvdb::Vec4I> quads);
std::vector<std::vector<float> > concatenateTriangleVectors(std::vector<openvdb::Vec3I> triangles, std::vector<std::vector<float> > triangles_from_splitted_quads);
std::vector<std::vector<float> > convertOpenVDBVectorToStandardVector(std::vector<openvdb::Vec3s> points);
std::vector<std::vector<float> > IncreaseTriangleVertexIndicesByN(std::vector<std::vector<float> > triangles, int N);
std::vector<std::vector<float> > DecreaseTriangleVertexIndicesByN(std::vector<std::vector<float> > triangles, int N);

#endif
