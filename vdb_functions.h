#ifndef VDBTEST1_H
#define VDBTEST1_H

#include <openvdb/openvdb.h>
//#include <openvdb/math/Maps.cc>
#include <openvdb/Grid.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <stdio.h>
#include <stdlib.h>
#include <openvdb/tools/Composite.h>


openvdb::FloatGrid::Ptr loadData();
openvdb::FloatGrid::Ptr createBlock(float radius, float value);

#endif
