# include "vdb_functions.h"
#include <math.h>

// sudo g++ -L/usr/lib/x86_64-linux-gnu  vdb_test1.cpp -ltbb -lopenvdb -lHalf -o openvdb_test1

// home ubuntu compile: sudo g++ vdb_test1.cpp -lopenvdb -ltbb -lHalf

// // create shared library of this with sudo g++ -std=c++11 -fPIC -shared vdb_functions.cpp -o /usr/lib/libvdb_functions.so


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


std::vector<std::vector<float> > convertOpenVDBVectorToStandardVector(std::vector<openvdb::Vec3s> points)
{
  	int xyzs = 3;
	int dimension = points.size();
	std::vector<std::vector<float> > standard_points(dimension, std::vector<float>(xyzs));
	//http://stackoverflow.com/questions/21663256/how-to-initialize-a-vector-of-vectors

	for (int i=0;i<dimension;i++)
	{
		standard_points[i][0] = points[i].x();
		standard_points[i][1] = points[i].y();
		standard_points[i][2] = points[i].z();
	}

	return standard_points;

}

std::vector<std::vector<float> > splitQuadsToTriangles(std::vector<openvdb::Vec3s> points, std::vector<openvdb::Vec4I> quads)
{

	int number_of_splitted_triangles = 2*quads.size();
	int xyzs = 3;
	std::vector<std::vector<float> > triangles_from_splitted_quads(number_of_splitted_triangles, std::vector<float>(xyzs));

	for(int ui=0;ui<quads.size();ui++)
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

std::vector<std::vector<float> > concatenateTriangleVectors(std::vector<openvdb::Vec3I> triangles, std::vector<std::vector<float> > triangles_from_splitted_quads)
{
	int xyzs = 3;
	int number_of_total_triangles = triangles.size() + triangles_from_splitted_quads.size();
	std::vector<std::vector<float> > triangles_combined(number_of_total_triangles, std::vector<float>(xyzs));
	
	for (int i=0;i<triangles.size();i++)
	{
		for(int uj=0;uj<xyzs;uj++)
		{
			triangles_combined[i][uj] = triangles[i][uj];
		}
	}
	
	for (int i=triangles.size();i<number_of_total_triangles;i++)
	{
		int shifted_index = i - triangles.size();
		for(unsigned int uj=0;uj<xyzs;uj++)
		{
			triangles_combined[i][uj] = triangles_from_splitted_quads[shifted_index][uj];
		}
	}
	
	return triangles_combined;
}

std::vector<std::vector<float> > IncreaseTriangleVertexIndicesByN(std::vector<std::vector<float> > triangles, int N)
{
	int xyzs = 3;
	std::vector<std::vector<float> > triangles_indices_increased(triangles.size(), std::vector<float>(xyzs));
	for (int i=0;i<triangles.size();i++)
	{
		triangles_indices_increased[i][0] = triangles[i][0] + N;
		triangles_indices_increased[i][1] = triangles[i][1] + N;
		triangles_indices_increased[i][2] = triangles[i][2] + N;
	}
	return triangles_indices_increased;
}

std::vector<std::vector<float> > DecreaseTriangleVertexIndicesByN(std::vector<std::vector<float> > triangles, int N)
{
	int xyzs = 3;
	std::vector<std::vector<float> > triangles_indices_decreased(triangles.size(), std::vector<float>(xyzs));
	for (int i=0;i<triangles.size();i++)
	{
		triangles_indices_decreased[i][0] = triangles[i][0] - N;
		triangles_indices_decreased[i][1] = triangles[i][1] - N;
		triangles_indices_decreased[i][2] = triangles[i][2] - N;
	}
	return triangles_indices_decreased;
}


float GetLengthOfVector(std::vector<float> vec)
{
	float length = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	return length;
}

std::vector<float> NormalizeVector(std::vector<float> vec)
{
	int xyzs = 3;
	float length = GetLengthOfVector(vec);
	
	std::vector<float> normal(xyzs);
	
	normal[0] = vec[0]/length;
	normal[1] = vec[1]/length;
	normal[2] = vec[2]/length;
	
	return normal;
}


std::vector<float>  GetCrossProduct(std::vector<float> vec1, std::vector<float> vec2)
{
	//http://www.cplusplus.com/forum/general/77959/

	int xyzs = 3;
	std::vector<float> crossproduct(xyzs);

        //Cross product formula 
	crossproduct[0] = (vec1[1] * vec2[2]) - (vec1[2] * vec2[1]);
	crossproduct[1] = (vec1[2] * vec2[0]) - (vec1[0] * vec2[2]);
	crossproduct[2] = (vec1[0] * vec2[1]) - (vec1[1] * vec2[0]);

	return crossproduct;
}

std::vector<std::vector<float> > ComputeTriangleNormalsVDB(std::vector<openvdb::Vec3s> points, std::vector<std::vector<float> > triangles)
{
	int vertices_per_triangle = 3;
	int number_of_triangles = triangles.size();
	std::vector<std::vector<float> > triangle_normals(number_of_triangles, std::vector<float>(vertices_per_triangle));
	
	for (int i=0;i<triangles.size();i++)
	{
		// just take both the first two verts of the triangle
		int xyzs = 3;
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
		crossproduct = GetCrossProduct(vec1,vec2);
		
		normal = NormalizeVector(crossproduct);
		
		// write
		triangle_normals[i] = normal;
	}
	return triangle_normals;

}

std::vector<std::vector<float> > ComputeTriangleNormals(std::vector<std::vector<float> > points, std::vector<std::vector<float> > triangles)
{
	int vertices_per_triangle = 3;
	int number_of_triangles = triangles.size();
	std::vector<std::vector<float> > triangle_normals(number_of_triangles, std::vector<float>(vertices_per_triangle));
	
	for (int i=0;i<triangles.size();i++)
	{
		int xyzs = 3;
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
		crossproduct = GetCrossProduct(vec1,vec2);
		
		normal = NormalizeVector(crossproduct);
		
		// write
		
		triangle_normals[i] = normal;
	}
	return triangle_normals;

}

std::vector<std::vector<float> >  ComputeVertexNormals(std::vector<std::vector<float> > triangles, std::vector<openvdb::Vec3s> points, std::vector<std::vector<float> > triangle_normals)
{
	std::vector<std::vector<float> > vertex_normals;
	//average the normals of all the faces that share a triangle vertex
	int xyzs = 0;
	for (int v=0;v<points.size();v++) 
	{
		
		// i dont want to run through all elements but is there a faster way?
		float x;
		float y; 
		float z;
		float shared_triangle_counter = 0;
		
		for (int j=0;j<triangles.size();j++)
		{
			if ((triangles[j][0] == v) || (triangles[j][1] == v) || (triangles[j][2] == v))
			{
				
				x += triangle_normals[j][0];
				y += triangle_normals[j][1];
				z += triangle_normals[j][2];
				shared_triangle_counter += 1;
			}
		} 
		 
		std::vector<float> vertex_normal = {0,0,0};

		vertex_normal[0] = x/shared_triangle_counter;
		vertex_normal[1] = y/shared_triangle_counter;
		vertex_normal[2] = z/shared_triangle_counter;
		
		std::vector<float> normalized_vertex_normal;
		
		normalized_vertex_normal = NormalizeVector(vertex_normal);
		
		vertex_normals.push_back(normalized_vertex_normal);
	}
	return vertex_normals;
}


