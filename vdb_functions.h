#ifndef VDBTEST1_H
#define VDBTEST1_H

#include <openvdb/openvdb.h>
//#include <openvdb/math/Maps.cc>
#include <openvdb/Grid.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <stdio.h>
#include <stdlib.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/math/Coord.h>

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

// vector stuff for normal calculation

float GetLengthOfVector(std::vector<float> vec);
std::vector<float> NormalizeVector(std::vector<float> vec);
std::vector<float> GetCrossProduct(std::vector<float> vec1, std::vector<float> vec2);
std::vector<std::vector<float> > ComputeTriangleNormals(std::vector<std::vector<float> > points, std::vector<std::vector<float> > triangles);
std::vector<std::vector<float> > ComputeTriangleNormalsVDB(std::vector<openvdb::Vec3s> points, std::vector<std::vector<float> > triangles);
std::vector<std::vector<float> >  ComputeVertexNormals(std::vector<std::vector<float> > triangles, std::vector<openvdb::Vec3s> points, std::vector<std::vector<float> > triangle_normals);

// IO
void ExportMeshAsObj(std::string filename, std::vector<openvdb::Vec3s> points, std::vector<openvdb::Vec3I> triangles, std::vector<openvdb::Vec4I> quads);
void ExportMeshAsVDB(openvdb::FloatGrid::Ptr grid);


#endif
