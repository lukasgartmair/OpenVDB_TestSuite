#include "vdb_functions.h"
#include <math.h>
#include <numeric> // for inner product
#include <vector>
#include <functional>

// sudo g++ -L/usr/lib/x86_64-linux-gnu  vdb_test1.cpp -ltbb -lopenvdb -lHalf -o openvdb_test1

// home ubuntu compile: sudo g++ vdb_test1.cpp -lopenvdb -ltbb -lHalf

// // create shared library of this with sudo g++ -std=c++11 -fPIC -shared vdb_functions.cpp -o /usr/lib/libvdb_functions.so

namespace VDB{

	const unsigned int xyzs = 3;
	const unsigned int vertices_per_triangle = 3;
	
openvdb::FloatGrid::Ptr createBlock(float radius, float value)
{
	
	openvdb::FloatGrid::Ptr grid =
	openvdb::FloatGrid::create(/*background value=*/0);
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
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

std::vector<openvdb::Vec3s> volumeToMeshVertices(openvdb::FloatGrid::Ptr grid, float isovalue, float adaptivity)
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

coord getVoxelIndex(coord *vec, float voxsize)
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


std::vector<std::vector<float> > convertOpenVDBVectorToStandardVector(std::vector<openvdb::Vec3s> points)
{
  	
	int dimension = points.size();
	std::vector<std::vector<float> > standard_points(dimension, std::vector<float>(xyzs));
	//http://stackoverflow.com/questions/21663256/how-to-initialize-a-vector-of-vectors

	for (int i=0;i<dimension;++i)
	{
		standard_points[i][0] = points[i].x();
		standard_points[i][1] = points[i].y();
		standard_points[i][2] = points[i].z();
	}

	return standard_points;

}

std::vector<std::vector<unsigned int> > splitQuadsToTriangles(std::vector<openvdb::Vec3s> points, std::vector<openvdb::Vec4I> quads)
{

	unsigned int number_of_splitted_triangles = 2*quads.size();
	std::vector<std::vector<unsigned int> > triangles_from_splitted_quads(number_of_splitted_triangles, std::vector<unsigned int>(xyzs));

	for(int ui=0;ui<quads.size();++ui)
	{
		openvdb::Vec3s A = points[quads[ui][0]];
		openvdb::Vec3s B = points[quads[ui][1]];
		openvdb::Vec3s C = points[quads[ui][2]];
		openvdb::Vec3s D = points[quads[ui][3]];
	
		float distanceAC = 0;
		float distanceBD = 0;
		distanceAC = sqrt((A.x()-C.x())*(A.x()-C.x()) + (A.y()-C.y())*(A.y()-C.y()) +  (A.z()-C.z())*(A.z()-C.z()));
		distanceBD = sqrt((B.x()-D.x())*(B.x()-D.x()) + (B.y()-D.y())*(B.y()-D.y()) +  (B.z()-D.z())*(B.z()-D.z()));
		
		if (distanceAC <= distanceBD)
		{	
			//tri 1 ABC
			triangles_from_splitted_quads[ui*2][0] = quads[ui][0];
			triangles_from_splitted_quads[ui*2][1] = quads[ui][1];
			triangles_from_splitted_quads[ui*2][2] = quads[ui][2];
			//tri2 ACD
			triangles_from_splitted_quads[(ui*2)+1][0] = quads[ui][0];
			triangles_from_splitted_quads[(ui*2)+1][1] = quads[ui][2];
			triangles_from_splitted_quads[(ui*2)+1][2] = quads[ui][3];	
		}
		if (distanceAC > distanceBD)
		{
			//tri 1 ABC
			triangles_from_splitted_quads[ui*2][0] = quads[ui][0];
			triangles_from_splitted_quads[ui*2][1] = quads[ui][1];
			triangles_from_splitted_quads[ui*2][2] = quads[ui][3];
			//tri2 ACD
			triangles_from_splitted_quads[(ui*2)+1][0] = quads[ui][3];
			triangles_from_splitted_quads[(ui*2)+1][1] = quads[ui][1];
			triangles_from_splitted_quads[(ui*2)+1][2] = quads[ui][2];
		}
	}
	
	return triangles_from_splitted_quads;
}

// combine normal and splitted triangles

std::vector<std::vector<unsigned int> > concatenateTriangleVectors(std::vector<openvdb::Vec3I> triangles, std::vector<std::vector<unsigned int> > triangles_from_splitted_quads)
{
	unsigned int number_of_total_triangles = triangles.size() + triangles_from_splitted_quads.size();
	std::vector<std::vector<unsigned int> > triangles_combined(number_of_total_triangles, std::vector<unsigned int>(xyzs));
	
	for (unsigned int i=0;i<triangles.size();++i)
	{
		for(unsigned int uj=0;uj<xyzs;++uj)
		{
			triangles_combined[i][uj] = triangles[i][uj];
		}
	}
	
	for (unsigned int i=triangles.size();i<number_of_total_triangles;++i)
	{
		int shifted_index = i - triangles.size();
		for(unsigned int uj=0;uj<xyzs;++uj)
		{
			triangles_combined[i][uj] = triangles_from_splitted_quads[shifted_index][uj];
		}
	}
	
	return triangles_combined;
}

std::vector<std::vector<unsigned int> > increaseTriangleVertexIndicesByN(std::vector<std::vector<unsigned int> > triangles, int N)
{
	std::vector<std::vector<unsigned int> > triangles_indices_increased(triangles.size(), std::vector<unsigned int>(xyzs));
	for (unsigned int i=0;i<triangles.size();++i)
	{
		triangles_indices_increased[i][0] = triangles[i][0] + N;
		triangles_indices_increased[i][1] = triangles[i][1] + N;
		triangles_indices_increased[i][2] = triangles[i][2] + N;
	}
	return triangles_indices_increased;
}

std::vector<std::vector<unsigned int> > decreaseTriangleVertexIndicesByN(std::vector<std::vector<unsigned int> > triangles, int N)
{
	std::vector<std::vector<unsigned int> > triangles_indices_decreased(triangles.size(), std::vector<unsigned int>(xyzs));
	for (unsigned int i=0;i<triangles.size();++i)
	{
		triangles_indices_decreased[i][0] = triangles[i][0] - N;
		triangles_indices_decreased[i][1] = triangles[i][1] - N;
		triangles_indices_decreased[i][2] = triangles[i][2] - N;
	}
	return triangles_indices_decreased;
}


float getLengthOfVector(std::vector<float> vec)
{
	float length = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	return length;
}

std::vector<float> normalizeVector(std::vector<float> vec)
{
	float length = getLengthOfVector(vec);
	
	std::vector<float> normal(xyzs);
	
	normal[0] = vec[0]/length;
	normal[1] = vec[1]/length;
	normal[2] = vec[2]/length;
	
	return normal;
}


std::vector<float>  getCrossProduct(std::vector<float> vec1, std::vector<float> vec2)
{
	//http://www.cplusplus.com/forum/general/77959/

	std::vector<float> crossproduct(xyzs);

        //Cross product formula 
	crossproduct[0] = (vec1[1] * vec2[2]) - (vec1[2] * vec2[1]);
	crossproduct[1] = (vec1[2] * vec2[0]) - (vec1[0] * vec2[2]);
	crossproduct[2] = (vec1[0] * vec2[1]) - (vec1[1] * vec2[0]);

	return crossproduct;
}

std::vector<std::vector<float> > computeTriangleNormalsVDB(std::vector<openvdb::Vec3s> points, std::vector<std::vector<unsigned int> > triangles)
{
	
	unsigned int number_of_triangles = triangles.size();
	std::vector<std::vector<float> > triangle_normals(number_of_triangles, std::vector<float>(vertices_per_triangle));
	
	for (unsigned int i=0;i<triangles.size();++i)
	{
		std::vector<float> vec1(xyzs);
		std::vector<float> vec2(xyzs);
		std::vector<float> normal(xyzs);
		std::vector<float> crossproduct;
		
		//conversion
		//https://www.opengl.org/wiki/Calculating_a_Surface_Normal
		// U = p2 - p1 and the vector V = p3 - p1
		vec1[0] = points[triangles[i][1]].x() - points[triangles[i][0]].x();
		vec1[1] = points[triangles[i][1]].y() - points[triangles[i][0]].y();
		vec1[2] = points[triangles[i][1]].z() - points[triangles[i][0]].z();
		vec2[0] = points[triangles[i][2]].x() - points[triangles[i][0]].x();
		vec2[1] = points[triangles[i][2]].y() - points[triangles[i][0]].y();
		vec2[2] = points[triangles[i][2]].z() - points[triangles[i][0]].z();
		
		// calculation
		crossproduct = getCrossProduct(vec1,vec2);
		
		normal = normalizeVector(crossproduct);
		
		// write
		triangle_normals[i] = normal;
	}
	return triangle_normals;

}

std::vector<std::vector<float> > computeTriangleNormals(std::vector<std::vector<float> > points, std::vector<std::vector<unsigned int> > triangles)
{
	unsigned int number_of_triangles = triangles.size();
	std::vector<std::vector<float> > triangle_normals(number_of_triangles, std::vector<float>(vertices_per_triangle));
	
	for (unsigned int i=0;i<triangles.size();++i)
	{
		std::vector<float> vec1(xyzs);
		std::vector<float> vec2(xyzs);
		std::vector<float> normal(xyzs);
		std::vector<float> crossproduct;
		
		//conversion
		//https://www.opengl.org/wiki/Calculating_a_Surface_Normal
		// U = p2 - p1 and the vector V = p3 - p1
		vec1[0] = points[triangles[i][1]][0] - points[triangles[i][0]][0];
		vec1[1] = points[triangles[i][1]][1] - points[triangles[i][0]][1];
		vec1[2] = points[triangles[i][1]][2] - points[triangles[i][0]][2];
		vec2[0] = points[triangles[i][2]][0] - points[triangles[i][0]][0];
		vec2[1] = points[triangles[i][2]][1] - points[triangles[i][0]][1];
		vec2[2] = points[triangles[i][2]][2] - points[triangles[i][0]][2];
		
		// calculation
		crossproduct = getCrossProduct(vec1,vec2);
		
		normal = normalizeVector(crossproduct);
		
		// write
		
		triangle_normals[i] = normal;
	}
	return triangle_normals;

}


std::vector<std::vector<float> >  computeVertexNormals(std::vector<std::vector<unsigned int> > triangles, std::vector<openvdb::Vec3s> points, std::vector<std::vector<float> > triangle_normals)
{
	std::vector<std::vector<float> > vertex_normals(points.size(), std::vector<float>(xyzs));
	//average the normals of all the faces that share a triangle vertex
	for (unsigned int i=0;i<triangles.size();++i) 
	{
		// http://stackoverflow.com/questions/16340931/calculating-vertex-normals-of-a-mesh
		//Think it otherway round: Iterate over the faces and add to the normal of the vertex.
		
		
		for (unsigned int j=0;j<xyzs;++j)
		{
			unsigned int current_vertex = triangles[i][j];
			for (unsigned int k=0;k<xyzs;++k)
			{
				vertex_normals[current_vertex][k] += triangle_normals[i][k];
			}
		}	
	}
	
	// https://www.opengl.org/discussion_boards/showthread.php/126927-Averaging-normals
	// AvarageNormal = (Normal1 + Normal2 + Normal3) / length(Normal1 + Normal2 + Normal3); 
	
	for (unsigned int i=0;i<vertex_normals.size();++i) 	
	{
		vertex_normals[i] = normalizeVector(vertex_normals[i]);
	}
	return vertex_normals;
}

std::vector<float> computeTriangleAreas(std::vector<openvdb::Vec3s> points, std::vector<std::vector<unsigned int> > triangles)
{
// http://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
//Say you have 3 points A,B,C. Find the angle between AB and AC using dot product (i.e. AB⋅AC=|AB||AC|cosθ) and then you can find the area of the triangle using 

	std::vector<float> triangle_areas(triangles.size());

	//1st create 2 common starting point vectors
	std::vector<float> P12(3);
	std::vector<float> P13(3);

	for (unsigned int i=0;i<triangles.size();++i) 
	{

		// P1(a1,a2,a3);   P2(b1,b2,b3);   P3(c1,c2,c3)
		// P1P2→=<b1−a1,b2−a2,b3−a3>=<x1,y1,z1>
		// P1P3→=<c1−a1,c2−a2,c3−a3>=<x2,y2,z2>

		P12[0] = points[triangles[i][1]].x() - points[triangles[i][0]].x();
		P12[1] = points[triangles[i][1]].y() - points[triangles[i][0]].y();
		P12[2] = points[triangles[i][1]].z() - points[triangles[i][0]].z();

		P13[0] = points[triangles[i][2]].x() - points[triangles[i][0]].x();
		P13[1] = points[triangles[i][2]].y() - points[triangles[i][0]].y();
		P13[2] = points[triangles[i][2]].z() - points[triangles[i][0]].z();

		std::vector<float> crossproduct(xyzs);

		crossproduct = getCrossProduct(P12, P13);
		
		float area_of_parallelogram = 0;
		area_of_parallelogram = getLengthOfVector(crossproduct);

		float triangle_area = 0.5 * area_of_parallelogram;

		triangle_areas[i] = triangle_area;

	}

	return triangle_areas;
}

void exportTriangleAreas(std::vector<float> triangle_areas)
{

	FILE* f = fopen("exported_triangle_areas.txt","wt");

	for(unsigned int i=0;i<triangle_areas.size();i++) fprintf(f, "%i %lf\n", i ,triangle_areas[i]);
	fclose(f);

}

void exportTriangleMeshAsObj(std::vector<openvdb::Vec3s> points, std::vector<std::vector<unsigned int> > triangles)
{

	FILE* f = fopen("3depict_triangle_isosurface_mesh.obj","wt");

	for(unsigned int i=0;i<points.size();++i) fprintf(f, "v %lf %lf %lf\n", points[i].x(), points[i].y(), points[i].z());
	for(unsigned int i=0;i<triangles.size();++i) fprintf(f, "f %i %i %i\n", (unsigned int)(triangles[i][2]+1), (unsigned int)(triangles[i][1]+1), (unsigned int)(triangles[i][0]+1));

	fclose(f);

}


void exportVDBMeshAsObj(std::vector<openvdb::Vec3s> points, std::vector<openvdb::Vec3I> triangles, std::vector<openvdb::Vec4I> quads)
{

	FILE* f = fopen("3depict_vdb_isosurface_mesh.obj","wt");

	for(unsigned int i=0;i<points.size();++i) fprintf(f, "v %lf %lf %lf\n", points[i].x(), points[i].y(), points[i].z());
	for(unsigned int i=0;i<triangles.size();++i) fprintf(f, "f %d %d %d\n", triangles[i][2]+1, triangles[i][1]+1, triangles[i][0]+1);
	for(unsigned int i=0;i<quads.size();++i) fprintf(f, "f %d %d %d %d\n", quads[i][3]+1, quads[i][2]+1, quads[i][1]+1, quads[i][0]+1);

	fclose(f);

}

void exportMeshAsVDB(openvdb::FloatGrid::Ptr grid)
{

	openvdb::io::File file("grid.vdb");
	openvdb::GridPtrVec grids;
	grids.push_back(grid);

	file.write(grids);
	file.close();

}




}











