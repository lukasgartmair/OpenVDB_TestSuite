# include "vdb_functions.h"

// sudo g++ -L/usr/lib/x86_64-linux-gnu  vdb_test1.cpp -ltbb -lopenvdb -lHalf -o openvdb_test1

// home ubuntu compile: sudo g++ vdb_test1.cpp -lopenvdb -ltbb -lHalf



openvdb::FloatGrid::Ptr loadData()
{
	openvdb::initialize();
	// Create a VDB file object.
	openvdb::io::File file("bunny.vdb");
	// Open the file.  This reads the file header, but not any grids.
	file.open();
	// Loop over all grids in the file and retrieve a shared pointer
	// to the one named "LevelSetBlock".  (This can also be done
	// more simply by calling file.readGrid("LevelSetBlock").)
	openvdb::GridBase::Ptr baseGrid;
	baseGrid = file.readGrid("ls_bunny");
	file.close();
	// From the example above, "LevelSetBlock" is known to be a FloatGrid,
	// so cast the generic grid pointer to a FloatGrid pointer.
	openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

	return grid;
}

openvdb::FloatGrid::Ptr createBlock(float radius, float value)
{
	
	openvdb::FloatGrid::Ptr grid =
	openvdb::FloatGrid::create(/*background value=*/0);
	    // Get a voxel accessor.
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
	// Compute the signed distance from the surface of the Block of each
	// voxel within the bounding box and insert the value into the grid
	// if it is smaller in magnitude than the background value.
	openvdb::Coord ijk;
	int &i = ijk[0], &j = ijk[1], &k = ijk[2];
	for (i = -radius; i < radius; ++i) {
		for (j = -radius; j < radius; ++j) {
		    for (k = -radius; k < radius; ++k) {
			accessor.setValue(ijk, value);
		    }
		}
	}
	return grid;
}

std::vector<openvdb::Vec3s> volumeToMeshVertices(openvdb::FloatGrid::Ptr grid, double isovalue, double adaptivity)
{

openvdb::initialize();

// volume to mesh

  std::vector<openvdb::Vec3s> points;
  std::vector<openvdb::Vec3I> triangles;
  std::vector<openvdb::Vec4I> quads;
  
openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, triangles, quads, isovalue, adaptivity);

//std::cout << "points size" << " = " << points.size() << std::endl;

return points;

}

int roundUp(float numToRound, float multiple)
{
    if (multiple == 0)
        return numToRound;
    float remainder = fmod(abs(numToRound), multiple);
    if (remainder == 0)
        return numToRound;
    if (numToRound < 0)
        return -(abs(numToRound) - remainder);
    else
        return numToRound + multiple - remainder;
}

coord GetVoxelIndex(coord *vec, float voxsize)
{
    int cx = 0;
    int cy = 0;
    int cz = 0;

    cx = roundUp(vec->x, voxsize);
    cy = roundUp(vec->y, voxsize);
    cz = roundUp(vec->z, voxsize);

    coord voxel_index = {cx/voxsize, cy/voxsize , cz/voxsize};
    return voxel_index;
}


std::vector<std::vector<float> > ConvertOpenVDBVectorToStandardVector(std::vector<openvdb::Vec3s> points)
{
  
  int dimension = points.size();
  std::vector<std::vector<float> > standard_points(dimension, std::vector<float>(dimension));
  //http://stackoverflow.com/questions/21663256/how-to-initialize-a-vector-of-vectors
  
  
  for (int i=0;i<dimension;i++)
  {
  	standard_points[i][0] = points[i].x();
  	standard_points[i][1] = points[i].y();
  	standard_points[i][2] = points[i].z();
  }
	
  return standard_points;

}

/*
void writeVDBData(grid, points, triangles, quads, isovalue, adaptivity)
{

	FILE* f = fopen("testing_bunny.obj","wt");

	for(int i=0;i<points.size();i++) fprintf(f, "v %lf %lf %lf\n", points[i].x(), points[i].y(), points[i].z());
	for(int i=0;i<triangles.size();i++) fprintf(f, "f %d %d %d\n", triangles[i][2]+1, triangles[i][1]+1, triangles[i][0]+1);
	for(int i=0;i<quads.size();i++) fprintf(f, "f %d %d %d %d\n", quads[i][3]+1, quads[i][2]+1, quads[i][1]+1, quads[i][0]+1);

	fclose(f);

	openvdb::io::File file("vdb_test_grid_pos.vdb");
	openvdb::GridPtrVec grids;
	grids.push_back(grid);

	file.write(grids);
	file.close();

	// create pts file

	FILE* fpts = fopen("testing_bunny.pts","wt");

	for(int i=0;i<points.size();i++) fprintf(fpts, "v %.3f %.3f %.3f\n", points[i].x(), points[i].y(), points[i].z());

	fclose(fpts);

}
*/
