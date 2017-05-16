#ifndef VDB_TESTSUITE_H
#define VDB_TESTSUITE_H


#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>


#include "vdb_functions.h"

#include <algorithm>    // std::min

using namespace VDB;
 
class TestOpenVDB : public CppUnit::TestFixture {


public:
 
	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestOpenVDB");
 

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - Test the Test itsself",
				&TestOpenVDB::testOpenVDB_TestTheTest ));
		
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - Division of Data1 small by big",
				&TestOpenVDB::testOpenVDB_DivisionOfData1 ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - Division of Data2 big by small",
				&TestOpenVDB::testOpenVDB_DivisionOfData2 ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - Division of Data3 misc",
				&TestOpenVDB::testOpenVDB_DivisionOfData3 ));
		
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - Division of Data4 misc",
				&TestOpenVDB::testOpenVDB_DivisionOfData4 ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - Volume to mesh DMC",
				&TestOpenVDB::testOpenVDB_VolumeToMeshDMC ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - round up ",
				&TestOpenVDB::testOpenVDB_RoundUp ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - get voxel index ",
				&TestOpenVDB::testOpenVDB_GetVoxelIndex ));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - conversion openvdb vec3s to std::vec ",
				&TestOpenVDB::testOpenVDB_VectorConversion ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - split openvdb quads to triangles ",
				&TestOpenVDB::testOpenVDB_SplitQuadsToTriangles ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - concatenate triangle lists ",
				&TestOpenVDB::testOpenVDB_ConcatenateTriangles ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - increase triangle vertex indices by n ",
				&TestOpenVDB::testOpenVDB_IncreaseTrianglesVertexIndices ));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - decrease triangle vertex indices by n ",
				&TestOpenVDB::testOpenVDB_DecreaseTrianglesVertexIndices ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - compute vertex normals ",
				&TestOpenVDB::testOpenVDB_ComputeVertexNormals ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - compute triangle normals ",
				&TestOpenVDB::testOpenVDB_ComputeTriangleNormals ));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - bounding box ",
				&TestOpenVDB::testOpenVDB_BoundingBox));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - triangle areas ",
				&TestOpenVDB::testOpenVDB_TriangleAreas));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - export triangles to obj ",
				&TestOpenVDB::testOpenVDB_ExportTriangles));
	
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - active / inactive voxels ",
				&TestOpenVDB::testOpenVDB_ActiveInactiveVoxels));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - get voxel coords ",
				&TestOpenVDB::testOpenVDB_GetVoxelCoords));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - find coord in vector ",
				&TestOpenVDB::testOpenVDB_FindCoordInVector));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - signed distance field",
				&TestOpenVDB::testOpenVDB_SignedDistanceField));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - std::map creation",
				&TestOpenVDB::testOpenVDB_StdMap));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - sort",
				&TestOpenVDB::testOpenVDB_Sort));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - shell assignement",
				&TestOpenVDB::testOpenVDB_ShellAnalysis));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - rearrange the proximitiy ranges",
				&TestOpenVDB::testOpenVDB_RearrangeProximities));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - voxel intersections",
				&TestOpenVDB::testOpenVDB_VoxelIntersections));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - set proximity ranges",
				&TestOpenVDB::testOpenVDB_SetProximityRanges));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test - new shell check",
				&TestOpenVDB::testOpenVDB_ShellCheck));


		return suiteOfTests;
	}
 
	/// Setup method
	void setUp() {}
 
	/// Teardown method
	void tearDown() {}
 
protected:

	const unsigned int xyzs = 3;

	void testOpenVDB_TestTheTest() {
		int x = 1;
		int z = 2;
		int u = x + z;
		CPPUNIT_ASSERT_EQUAL(3, u);
	}
	
	void testOpenVDB_DivisionOfData1() {
	
		openvdb::initialize();
	
		// this test divides the small Block by the big one
		// that means there should be no division by zero because of background values.

		openvdb::FloatGrid::Ptr small_grid = openvdb::FloatGrid::create(/*background value=*/0);;
		openvdb::FloatGrid::Ptr big_grid = openvdb::FloatGrid::create(/*background value=*/0);;
		// radius, value
		small_grid = createBlock(10,1);
		big_grid = createBlock(20,1);
		
		// composite.h operations modify the first grid and leave the second grid emtpy
		// compute a = a / b
		openvdb::tools::compDiv(*small_grid, *big_grid);

		int number_of_nfs = 0;
		
		for (openvdb::FloatGrid::ValueOnIter iter = small_grid->beginValueOn(); iter; ++iter)
		{
		    if (std::isfinite(iter.getValue()) == false)
		{
			number_of_nfs += 1;		    
		}
		}

		//std::cout << "number of nfs" << " = " << number_of_nfs << std::endl;
		CPPUNIT_ASSERT_EQUAL(0, number_of_nfs);
		
	}
	
	void testOpenVDB_DivisionOfData2() {
	
		// this test divides the big Block by the small one
		// that means there should be division by zero
		// another thing occured here i did not see before:
		// instead of -nf from 0/0 division inf occurs from the division 
		// so i have to check both cases http://stackoverflow.com/questions/4095337/how-to-check-for-inf-and-or-nf-in-a-float-variable
		// with std::isfinite(x) --> if inf or -nf is finite results in false
		openvdb::initialize();
		openvdb::FloatGrid::Ptr small_grid = openvdb::FloatGrid::create(/*background value=*/0);
		openvdb::FloatGrid::Ptr big_grid = openvdb::FloatGrid::create(/*background value=*/0);
		// radius, value

		small_grid = createBlock(10,1);
		big_grid = createBlock(20,1);
		
		// composite.h operations modify the first grid and leave the second grid emtpy
		// compute a = a / b
		openvdb::tools::compDiv(*big_grid, *small_grid);

		// nf equals non finite
		int number_of_nfs = 0;
		
		for (openvdb::FloatGrid::ValueOnIter iter = big_grid->beginValueOn(); iter; ++iter)
		{
		    if (std::isfinite(iter.getValue()) == false)
		{
			number_of_nfs += 1;	    
		}
		}
		
		//std::cout << "number of nfs" << " = " << number_of_nfs << std::endl;
		float minVal = 0.0;
		float maxVal = 0.0;
		big_grid->evalMinMax(minVal,maxVal);
		//std::cout << " eval min max big grid" << " = " << minVal << " , " << maxVal << std::endl;
		//std::cout << " active voxel count big grid" << " = " << big_grid->activeVoxelCount() << std::endl;
		CPPUNIT_ASSERT( number_of_nfs > 0);

	}
	
	void testOpenVDB_DivisionOfData3() {
	
		// only one point is set to > 0 in order to obtain every active voxel division nf except 
		// grids (0,0,0) value is 0 divided by second_grids (0,0,0) value which is < 0
		openvdb::initialize();
		openvdb::FloatGrid::Ptr numerator_grid = openvdb::FloatGrid::create(/*background value=*/0);
		openvdb::FloatGrid::Ptr denominator_grid = openvdb::FloatGrid::create(/*background value=*/0);
		
		numerator_grid = createBlock(10,0);
		denominator_grid = createBlock(10,0);
		
		// set one single point with value > 0 in the 
		
		openvdb::FloatGrid::Accessor accessor = denominator_grid->getAccessor();
		openvdb::Coord ijk(0,0,0);
		int value = 1;
		accessor.setValue(ijk, value);

		// composite.h operations modify the first grid and leave the second grid emtpy
		// compute a = a / b
		openvdb::tools::compDiv(*numerator_grid, *denominator_grid);

		int number_of_nfs = 0;
		
		for (openvdb::FloatGrid::ValueOnIter iter = numerator_grid->beginValueOn(); iter; ++iter)
		{
		    if (std::isfinite(iter.getValue()) == false)
		{
			number_of_nfs += 1;	    
		}
		}

		//std::cout << "number of nfs" << " = " << number_of_nfs << std::endl; 
		int assert_number_of_nfs = numerator_grid->activeVoxelCount() - 1;
		CPPUNIT_ASSERT_EQUAL(assert_number_of_nfs,  number_of_nfs);
	}
	
	void testOpenVDB_DivisionOfData4() {
	
		// in this test the background is set to one which should produce not a single -nf
		openvdb::initialize();
		openvdb::FloatGrid::Ptr numerator_grid;
		openvdb::FloatGrid::Ptr denominator_grid;
		
		numerator_grid = createBlock(10,1);
		denominator_grid = createBlock(10,1);
		
		// set few points in grid
		
		openvdb::FloatGrid::Accessor accessor = numerator_grid->getAccessor();
		openvdb::Coord ijk(0,0,0);
		int value = 2;
		accessor.setValue(ijk, value);

		// composite.h operations modify the first grid and leave the second grid emtpy
		// compute a = a / b
		openvdb::tools::compDiv(*numerator_grid, *denominator_grid);
		
		int number_of_nfs = 0;
		
		for (openvdb::FloatGrid::ValueOnIter iter = numerator_grid->beginValueOn(); iter; ++iter)
		{
		    if (std::isfinite(iter.getValue()) == false)
		{
			number_of_nfs += 1;	    
		}
		}

		//std::cout << "number of nfs" << " = " << number_of_nfs << std::endl;
		CPPUNIT_ASSERT_EQUAL( 0, number_of_nfs);
	}
	
	
	void testOpenVDB_VolumeToMeshDMC()
	{
	
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(/*background value=*/0);
		// radius, value
		grid = createBlock(10,1);
		
		std::vector<openvdb::Vec3s> vertices;
		int vx = 0;
		// change the isovalue with no adaptivity
		vertices = volumeToMeshVertices(grid, 0.5, 0);
		vx = vertices.size();
		CPPUNIT_ASSERT_EQUAL( 2402, vx);
		vertices = volumeToMeshVertices(grid, 0.01, 0);
		vx = vertices.size();
		CPPUNIT_ASSERT_EQUAL(2402, vx);
		vx = vertices.size();
		vertices = volumeToMeshVertices(grid, 0.99, 0);
		CPPUNIT_ASSERT_EQUAL(2402, vx);
		
		// change the adaptivity with constant isovalue
		vertices = volumeToMeshVertices(grid, 0.5, 0.25);
		vx = vertices.size();
		CPPUNIT_ASSERT_EQUAL(566, vx);
		vertices = volumeToMeshVertices(grid, 0.5, 0.5);
		vx = vertices.size();
		CPPUNIT_ASSERT_EQUAL(422, vx);
		vertices = volumeToMeshVertices(grid, 0.5, 1);
		vx = vertices.size();
		CPPUNIT_ASSERT_EQUAL(56, vx);
		
		
		
		// got -nans as coordinates even if there were no -nans in the grid
		// this test should check how this is possible
		// i.e. not only to check the amount but the content!
		int non_finites_counter = 0;
		for (unsigned int i=0;i<vertices.size();i++)
		{
			for (unsigned int j=0;j<xyzs;j++)
			{
			    if (std::isfinite(vertices[i][j]) == false)
				{	
					non_finites_counter += 1;    
				}
			}
		}
		int assert_non_finites = 0;
		CPPUNIT_ASSERT_EQUAL(assert_non_finites, non_finites_counter);
	
		// got -nans as coordinates even if there were no -nans in the grid
		// this test should check how this is possible
		// i.e. not only to check the amount but the content!
		non_finites_counter = 0;
		
		vertices[0][0] = 0.0 / 0.0; // which results in -nan
 		
		for (unsigned int i=0;i<vertices.size();i++)
		{
			for (unsigned int j=0;j<xyzs;j++)
			{
			    if (std::isfinite(vertices[i][j]) == false)
				{	
					non_finites_counter += 1;    
				}
			}
		}
		
		assert_non_finites = 1;
		CPPUNIT_ASSERT_EQUAL(assert_non_finites, non_finites_counter);	
	
	}
	
	void testOpenVDB_RoundUp()
	{
	// check the helper function round which is used to determine the voxel indices
	// it has certainly to be discussed whether this behaviour is suitable to fill
	// the ion grids !
	
		int rounded_result = 0;
		rounded_result = roundUp(0.0, 5);
		//std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT_EQUAL(0, rounded_result);
		
		rounded_result = 0;
		rounded_result = roundUp(0.5, 5);
		//std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT_EQUAL(0, rounded_result);
		
		rounded_result = 0;
		rounded_result = roundUp(0.99, 5);
		//std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT_EQUAL(0, rounded_result);
		
		rounded_result = 0;
		rounded_result = roundUp(1.0, 5);
		//std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT_EQUAL(5, rounded_result);
		
		rounded_result = 0;
		rounded_result = roundUp(4.9, 5);
		//std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT_EQUAL(5, rounded_result);

		rounded_result = 0;
		rounded_result = roundUp(5.1, 5);
		//std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT_EQUAL(5, rounded_result);
		
		rounded_result = 0;
		rounded_result = roundUp(5.9, 5);
		//std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT_EQUAL(5, rounded_result);

		rounded_result = 0;
		rounded_result = roundUp(6.0, 5);
		//std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT_EQUAL(10, rounded_result);

	}
	
	
	void testOpenVDB_GetVoxelIndex()
	{
		// this function is of course highly dependent on its helper function RoundUp
	
		coord test_coord = {0,0,0};
		coord result_coord = {0,0,0};
		float voxelsize = 3;
		result_coord = getVoxelIndex(&test_coord, voxelsize);
		coord assert_coord = {0,0,0};
		CPPUNIT_ASSERT_EQUAL(result_coord.x, assert_coord.x);
		CPPUNIT_ASSERT_EQUAL(result_coord.y, assert_coord.y);
		CPPUNIT_ASSERT_EQUAL(result_coord.z, assert_coord.z);	
		
		test_coord = {12.5,-22,-3};
		result_coord = {0,0,0};
		voxelsize = 3;
		result_coord = getVoxelIndex(&test_coord, voxelsize);
		assert_coord = {4,-7,-1};
		//std::cout << "x" << " = " << result_coord.x << std::endl; 
		//std::cout << "y" << " = " << result_coord.y << std::endl; 
		//std::cout << "z" << " = " << result_coord.z << std::endl; 
		CPPUNIT_ASSERT_EQUAL(assert_coord.x, result_coord.x);
		CPPUNIT_ASSERT_EQUAL(assert_coord.y, result_coord.y);
		CPPUNIT_ASSERT_EQUAL(assert_coord.z, result_coord.z);	
		
		test_coord = {12.5,0.1,-3};
		result_coord = {0,0,0};
		voxelsize = 3;
		result_coord = getVoxelIndex(&test_coord, voxelsize);
		assert_coord = {4,0,-1};
		//std::cout << "x" << " = " << result_coord.x << std::endl; 
		//std::cout << "y" << " = " << result_coord.y << std::endl; 
		//std::cout << "z" << " = " << result_coord.z << std::endl; 
		CPPUNIT_ASSERT_EQUAL(assert_coord.x, result_coord.x);
		CPPUNIT_ASSERT_EQUAL(assert_coord.y, result_coord.y);
		CPPUNIT_ASSERT_EQUAL(assert_coord.z, result_coord.z);			
	}	
	
	
	void testOpenVDB_VectorConversion()
	{
		// set up some vertices
		
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(/*background value=*/0);
		grid = createBlock(10,1);	
		std::vector<openvdb::Vec3s> vertices;
		vertices = volumeToMeshVertices(grid, 0.5, 0);

		std::vector<std::vector<float> > standard_points(vertices.size(), std::vector<float>(xyzs));

		standard_points = convertOpenVDBVectorToStandardVector(vertices);
		
		for (unsigned int i=0;i<vertices.size();i++)
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL(vertices[i].x(), standard_points[i][0],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(vertices[i].y(), standard_points[i][1],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(vertices[i].z(), standard_points[i][2],0.01);
		}
	
	}
	
	void testOpenVDB_SplitQuadsToTriangles()
	{	

		// test if the sum is correct
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(2,1);	
		std::vector<openvdb::Vec3s> points;
		std::vector<openvdb::Vec3I> triangles;
		std::vector<openvdb::Vec4I> quads;

		float isovalue=0.5;
		float adaptivity=0;
		openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, triangles, quads, isovalue, adaptivity);

		std::vector<std::vector<unsigned int> > triangles_from_splitted_quads;

		//triangles_from_splitted_quads = splitQuadsToTriangles(points, quads);
		
		//std::cout << triangles_from_splitted_quads.size() << std::endl; 
		
		//CPPUNIT_ASSERT_EQUAL(quads.size()*2,triangles_from_splitted_quads.size());
		
		// faces outputs from openvdb volume to mesh starting from zero!
		unsigned int minimum  = 20;
		unsigned int minimum_tris = 20;
		for (unsigned int i=0;i<triangles.size();i++)
		{
			for (unsigned int j=0;j<3;j++)
			{
				if (triangles[i][j] < minimum)
				{
					minimum = triangles[i][j];
				}
			}
		}
		
		unsigned int minimum_quads = 20;
		for (unsigned int i=0;i<quads.size();i++)
		{
			for (unsigned int j=0;j<3;j++)
			{
				if (quads[i][j] < minimum_quads)
				{
					minimum_quads = quads[i][j];
				}
			}
		}
		
		if (minimum_tris < minimum_quads)
		{
			minimum = minimum_tris;
		}
		else
		{
			minimum = minimum_quads;
		}

		unsigned int assert_faces_minimum = 0;
		CPPUNIT_ASSERT_EQUAL(assert_faces_minimum, minimum);

	
		// simple plane 
		
		grid = createBlock(1,1);	
		isovalue=0.1;
		adaptivity=0;
		openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, quads, isovalue);
		points.resize(4);
		points[0][0] = -1.0;
		points[0][1] = 0.0;
		points[0][2] = 1.0;
		points[1][0] = 1.0;
		points[1][1] = 0.0;
		points[1][2] = 1.0;
		points[2][0] = -1.0;
		points[2][1] = 0.0;
		points[2][2] = -1.0;
		points[3][0] = 1.0;
		points[3][1] = 0.0;
		points[3][2] = -1.0;
		quads.resize(1);
		quads[0][0] = 0;
		quads[0][1] = 1;
		quads[0][2] = 2;
		quads[0][3] = 3;
		
		triangles_from_splitted_quads = splitQuadsToTriangles(points, quads);
		unsigned int tri_size = triangles_from_splitted_quads.size();
		unsigned int assert_size = 2;
		CPPUNIT_ASSERT_EQUAL(assert_size,tri_size);
		
		//check the contents
		std::vector<std::vector<unsigned int> > assert_triangles(tri_size, std::vector<unsigned int>(xyzs));
		
		assert_triangles[0][0] = 0; 
		assert_triangles[0][1] = 1; 
		assert_triangles[0][2] = 2; 
		assert_triangles[1][0] = 0; 
		assert_triangles[1][1] = 2; 
		assert_triangles[1][2] = 3; 	
		
		for (unsigned int i=0;i<tri_size;i++)
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][0], triangles_from_splitted_quads[i][0],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][1], triangles_from_splitted_quads[i][1],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][2], triangles_from_splitted_quads[i][2],0.01);
		}
	
	}
	
	void testOpenVDB_ConcatenateTriangles()
	{
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(2,1);	
		std::vector<openvdb::Vec3s> points;
		std::vector<openvdb::Vec3I> triangles;
		std::vector<openvdb::Vec4I> quads;

		float isovalue=0.5;
		float adaptivity=0.5;
		openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, triangles, quads, isovalue, adaptivity);

		std::vector<std::vector<unsigned int> > triangles_from_splitted_quads;
		triangles_from_splitted_quads = splitQuadsToTriangles(points, quads);
		
		std::vector<std::vector<unsigned int> > triangles_combined;
		
		triangles_combined = concatenateTriangleVectors(triangles, triangles_from_splitted_quads);
		
		//std::cout << triangles_combined[triangles.size()][0] << std::endl; 
		//std::cout << triangles_combined[triangles.size()][1] << std::endl; 
		//std::cout << triangles_combined[triangles.size()][2] << std::endl; 

		// check the amount
		unsigned int tri_size = triangles_combined.size();
		unsigned int assert_size = (quads.size()*2) + triangles.size();
		CPPUNIT_ASSERT_EQUAL(assert_size, tri_size);

		//check the contents
		
		for (unsigned int i=0;i<triangles.size();i++)
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles[i][0], triangles_combined[i][0],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles[i][1], triangles_combined[i][1],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles[i][2], triangles_combined[i][2],0.01);
		}
		
		for (unsigned int i=triangles.size();i<2*quads.size()+triangles.size();i++)
		{
			int shifted_index = i - triangles.size();
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles_from_splitted_quads[shifted_index][0], triangles_combined[i][0],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles_from_splitted_quads[shifted_index][1], triangles_combined[i][1],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles_from_splitted_quads[shifted_index][2], triangles_combined[i][2],0.01);
		}

	}
	void testOpenVDB_IncreaseTrianglesVertexIndices()
	{
		unsigned int tri_size = 2;
		std::vector<std::vector<unsigned int> > triangles(tri_size, std::vector<unsigned int>(xyzs));
	
		triangles[0][0] = 0; 
		triangles[0][1] = 1; 
		triangles[0][2] = 2; 
		triangles[1][0] = 0; 
		triangles[1][1] = 2; 
		triangles[1][2] = 3; 	

		std::vector<std::vector<unsigned int> > assert_triangles(tri_size, std::vector<unsigned int>(xyzs));
		
		assert_triangles[0][0] = 1; 
		assert_triangles[0][1] = 2; 
		assert_triangles[0][2] = 3; 
		assert_triangles[1][0] = 1; 
		assert_triangles[1][1] = 3; 
		assert_triangles[1][2] = 4; 
		
		std::vector<std::vector<unsigned int> > result_triangles(tri_size, std::vector<unsigned int>(xyzs));
		int N = 1;
		result_triangles = increaseTriangleVertexIndicesByN(triangles, N);
		
		//check the contents
		for (unsigned int i=0;i<tri_size;i++)
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][0], result_triangles[i][0],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][1], result_triangles[i][1],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][2], result_triangles[i][2],0.01);
		}
	}
	
	void testOpenVDB_DecreaseTrianglesVertexIndices()
	{
		const unsigned int tri_size = 2;
		std::vector<std::vector<unsigned int> > triangles(tri_size, std::vector<unsigned int>(xyzs));
	
		triangles[0][0] = 1; 
		triangles[0][1] = 2; 
		triangles[0][2] = 3; 
		triangles[1][0] = 1; 
		triangles[1][1] = 3; 
		triangles[1][2] = 4; 	

		std::vector<std::vector<unsigned int> > assert_triangles(tri_size, std::vector<unsigned int>(xyzs));
		
		assert_triangles[0][0] = 0; 
		assert_triangles[0][1] = 1; 
		assert_triangles[0][2] = 2; 
		assert_triangles[1][0] = 0; 
		assert_triangles[1][1] = 2; 
		assert_triangles[1][2] = 3; 
		
		std::vector<std::vector<unsigned int> > result_triangles(tri_size, std::vector<unsigned int>(xyzs));
		int N = 1;
		result_triangles = decreaseTriangleVertexIndicesByN(triangles, N);
		
		//check the contents
		for (unsigned int i=0;i<tri_size;i++)
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][0], result_triangles[i][0],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][1], result_triangles[i][1],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][2], result_triangles[i][2],0.01);
		}
	}
	
	void testOpenVDB_ComputeVertexNormals()
	{

		const unsigned int tri_size = 2;
		std::vector<std::vector<unsigned int> > triangles(tri_size, std::vector<unsigned int>(xyzs));
		
		triangles[0][0] = 0; 
		triangles[0][1] = 1; 
		triangles[0][2] = 2; 
		triangles[1][0] = 0; 
		triangles[1][1] = 2; 
		triangles[1][2] = 3; 
	
		// vdb version
		
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(1,1);	
		std::vector<openvdb::Vec3s> points;
		std::vector<openvdb::Vec4I> quads;

		float isovalue=0.5;
		float adaptivity=0;
		openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, quads, isovalue);
		points.resize(4);
		
		points[0][0] = 0;
		points[0][1] = 0;
		points[0][2] = 0;
		points[1][0] = 1.0;
		points[1][1] = 0;
		points[1][2] = 0;
		points[2][0] = 1.0;
		points[2][1] = 1.0;
		points[2][2] = 1.0;
		points[3][0] = 0.0;
		points[3][1] = 1.0;
		points[3][2] = 1.0;

		std::vector<std::vector<float> > triangle_normals;
		
		triangle_normals = computeTriangleNormalsVDB(points, triangles);
	
		std::vector<std::vector<float> > vertex_normals;
	
		vertex_normals = computeVertexNormals(triangles, points, triangle_normals);

	
		// check the amount
		const unsigned int size_normals = 4;
		CPPUNIT_ASSERT_DOUBLES_EQUAL(size_normals, vertex_normals.size(),0.01);
		// check the content
		/*
		std::cout << vertex_normals[0][0] << std::endl; 
		std::cout << vertex_normals[0][1] << std::endl; 
		std::cout << vertex_normals[0][2] << std::endl; 
		std::cout << vertex_normals[1][0] << std::endl; 
		std::cout << vertex_normals[1][1] << std::endl; 
		std::cout << vertex_normals[1][2] << std::endl; 
		std::cout << vertex_normals[2][0] << std::endl; 
		std::cout << vertex_normals[2][1] << std::endl; 
		std::cout << vertex_normals[2][2] << std::endl;
		std::cout << vertex_normals[3][0] << std::endl; 
		std::cout << vertex_normals[3][1] << std::endl; 
		std::cout << vertex_normals[3][2] << std::endl;  
		*/
	}


	void testOpenVDB_ComputeTriangleNormals()
	{
	
		const unsigned int tri_size = 2;
		std::vector<std::vector<float> > normals;
		
		// vdb version
		
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(1,1);	
		std::vector<openvdb::Vec3s> points;
		std::vector<openvdb::Vec4I> quads;

		float isovalue=0.5;
		float adaptivity=0;
		openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, quads, isovalue);
		
		points.resize(3);
		// both zs == 0 i.e. the triangle is parallel to xy plane
		// normal should be 001
		
		std::vector<std::vector<unsigned int> > triangle(1, std::vector<unsigned int>(xyzs));
		triangle[0][0] = 1; 
		triangle[0][1] = 2; 
		triangle[0][2] = 0; 
		
		points[0][0] = 0.0;
		points[0][1] = 0.0;
		points[0][2] = 0.0;
		points[1][0] = 1.0;
		points[1][1] = 0.0;
		points[1][2] = 0.0;
		points[2][0] = 1.0;
		points[2][1] = 1.0;
		points[2][2] = 0.0;
		
		normals = computeTriangleNormalsVDB(points, triangle);

		CPPUNIT_ASSERT_DOUBLES_EQUAL(0, normals[0][0],0.01);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(0, normals[0][1],0.01);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(1, normals[0][2],0.01);

		// std::vector version
		int points_size = 3;
		std::vector<std::vector<float> > std_points(points_size, std::vector<float>(xyzs));
		
		std_points[0][0] = 0.0;
		std_points[0][1] = 0.0;
		std_points[0][2] = 0.0;
		std_points[1][0] = 1.0;
		std_points[1][1] = 0.0;
		std_points[1][2] = 0.0;
		std_points[2][0] = 1.0;
		std_points[2][1] = 1.0;
		std_points[2][2] = 0.0;
		
		normals = computeTriangleNormals(std_points, triangle);

		CPPUNIT_ASSERT_DOUBLES_EQUAL(0, normals[0][0],0.01);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(0, normals[0][1],0.01);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(1, normals[0][2],0.01);

	}

	void testOpenVDB_BoundingBox()
	{

		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(15,1);	

		openvdb::math::CoordBBox bounding_box = openvdb::math::CoordBBox();
		bounding_box = grid->evalActiveVoxelBoundingBox();
/*
		std::cout << bounding_box.getStart() << std::endl;
            	std::cout << bounding_box.getStart()[0] << std::endl;

		std::cout << bounding_box.getEnd() << std::endl;
		std::cout << bounding_box.getCenter() << std::endl;
		std::cout << bounding_box.dim() << std::endl;
*/		
	}

	void testOpenVDB_TriangleAreas()
	{
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(2,1);	
		std::vector<openvdb::Vec3s> points;
		std::vector<openvdb::Vec3I> triangles;
		std::vector<openvdb::Vec4I> quads;

		float isovalue=0.5;
		float adaptivity=0;
		openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, triangles, quads, isovalue, adaptivity);


		std::vector<std::vector<unsigned int> > triangles_from_splitted_quads;
		triangles_from_splitted_quads = splitQuadsToTriangles(points, quads);
		
		std::vector<std::vector<unsigned int> > triangles_combined;
		
		triangles_combined = concatenateTriangleVectors(triangles, triangles_from_splitted_quads);


		std::vector<float> triangle_areas(triangles_combined.size());
		triangle_areas = computeTriangleAreas(points, triangles_combined);
		/*
		std::cout << triangle_areas[0] << std::endl; 
		std::cout << triangle_areas[0] << std::endl; 
		std::cout << triangle_areas[0] << std::endl; 
		*/
	}

	void testOpenVDB_ExportTriangles()
	{
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(2,1);	
		std::vector<openvdb::Vec3s> points;
		std::vector<openvdb::Vec3I> triangles;
		std::vector<openvdb::Vec4I> quads;

		float isovalue=0.5;
		float adaptivity=0;
		openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, triangles, quads, isovalue, adaptivity);

		std::vector<std::vector<unsigned int> > triangles_from_splitted_quads;
		triangles_from_splitted_quads = splitQuadsToTriangles(points, quads);
		
		std::vector<std::vector<unsigned int> > triangles_combined;
		
		triangles_combined = concatenateTriangleVectors(triangles, triangles_from_splitted_quads);

		
		exportTriangleMeshAsObj(points, triangles_combined);

	}



	void testOpenVDB_ActiveInactiveVoxels()
	{

		// only iterate the active voxels
		// so both the active and inactive voxels should have the value zero
		// but nevertheless different activation states - is that possible?
		// another possibility is to set the initial active value to one or something 
		// unequal to zero and just overwrite it the first time of writing 

		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(2,1);
		openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

		int active_voxels = grid->activeVoxelCount();
		int assert_voxel_number = 64;

		float minVal = 0.0;
		float maxVal = 0.0;
		grid->evalMinMax(minVal,maxVal);
		//std::cout << " eval min max denominator grid" << " = " << minVal << " , " << maxVal << std::endl;
		//std::cout << " active voxel count denominator grid" << " = " << grid->activeVoxelCount() << std::endl;

		for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter; ++iter)
		{   
				iter.setValue(0.0);
		}

		grid->evalMinMax(minVal,maxVal);

		// are now all values zero but nevertheless different voxel states present?
		CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_voxel_number, active_voxels,0.01);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(0, maxVal,0.01);

		// yes

		// further check how to get the current state of a particular voxel
	

		openvdb::Coord ijk(0,0,0);
		
	}



	void testOpenVDB_GetVoxelCoords()
	{

		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(2,1);
		openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

		std::vector<std::vector<float> > active_voxel_indices(grid->activeVoxelCount(), std::vector<float>(xyzs));
		openvdb::Coord hkl;

		unsigned int voxel_counter = 0;
		for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter; ++iter)
		{   
	    			iter.setValue(0.0);
			
				hkl = iter.getCoord();
				active_voxel_indices[voxel_counter][0] = hkl.x();
				active_voxel_indices[voxel_counter][1] = hkl.y();
				active_voxel_indices[voxel_counter][2] = hkl.z();
				voxel_counter += 1;
		}
		
		// is the voxel counter running properly?
		CPPUNIT_ASSERT_DOUBLES_EQUAL(grid->activeVoxelCount(), voxel_counter,0.01);
	}

	void testOpenVDB_FindCoordInVector()
	{

		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(2,1);
		openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

		std::vector<std::vector<float> > active_voxel_indices(grid->activeVoxelCount(), std::vector<float>(xyzs));
		openvdb::Coord hkl;

		unsigned int voxel_counter = 0;
		for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter; ++iter)
		{   
	    			iter.setValue(0.0);
			
				hkl = iter.getCoord();
				active_voxel_indices[voxel_counter][0] = hkl.x();
				active_voxel_indices[voxel_counter][1] = hkl.y();
				active_voxel_indices[voxel_counter][2] = hkl.z();
				voxel_counter += 1;
		}

		//initialize a (0,0,0) vector
		std::vector<float> vector_to_search(xyzs);

		vector_to_search[0] = 0;

		bool vector_found = false;

		if(std::find(active_voxel_indices.begin(), active_voxel_indices.end(), vector_to_search ) != active_voxel_indices.end()) 
		{
			/* v contains x */
			vector_found = true;
		} 
		else 
		{
			/* v does not contain x */
			vector_found = false;		
		}	

		CPPUNIT_ASSERT_DOUBLES_EQUAL(true, vector_found,0.01);	

	}

	void testOpenVDB_SignedDistanceField()
	{
		// extract mesh
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(2,1);	
		std::vector<openvdb::Vec3s> points;
		std::vector<openvdb::Vec3I> triangles;
		std::vector<openvdb::Vec4I> quads;

		float isovalue=0.5;
		float adaptivity=0;
		openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, triangles, quads, isovalue, adaptivity);

		float in_bandwidth = 2;
		float ex_bandwidth = 2;

		float voxelsize_levelset = 0.5;

		// this is the right way to call the sdf !!!
		// signed distance field
		openvdb::math::Transform::Ptr trans = openvdb::math::Transform::createLinearTransform(voxelsize_levelset);
		openvdb::FloatGrid::Ptr sdf = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(*trans, points, triangles, quads, ex_bandwidth, in_bandwidth);

		unsigned int active_voxels_meshgrid = grid->activeVoxelCount();
		unsigned int active_voxels_sdf = sdf->activeVoxelCount();

		//std::cout << " active_voxels_meshgrid" << " = " << active_voxels_meshgrid << std::endl;
		//std::cout << " active_voxels_sdf" << " = " << active_voxels_sdf << std::endl;

		bool narrowband = false;
		if (active_voxels_meshgrid < active_voxels_sdf)	
		{
			narrowband = true;
		}
		CPPUNIT_ASSERT_DOUBLES_EQUAL(true, narrowband,0.01);
	
	}

	void testOpenVDB_StdMap()
	{
		std::map<float, int> indicesMap;
		int increaser = 10;
		for (int i=0;i<3;i++)
		{
			indicesMap.insert(std::make_pair(increaser, i));
			increaser += 10;
		}

		int assert_key = 2;
		int ask_for_key = indicesMap[30];
		

		CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_key,ask_for_key,0.01);
	}


	void testOpenVDB_Sort()
	{
		std::vector<unsigned int> vec = {3,1,2};

		std::sort(vec.begin(), vec.end());

		CPPUNIT_ASSERT_DOUBLES_EQUAL(1,vec[0],0.01);
	}


	void testOpenVDB_ShellAnalysis()
	{
		std::vector<float> unique_distances = {-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5};
		std::vector<float> voxel_content = {1,1,1,100,100,100,100,100,10,10,10,10,10,10};

		// does the summarization of the shells work as expected?

		const unsigned int number_of_proximity_ranges = 3;

		std::vector<float> proximity_ranges_ends = {-0.5, 2, 5};

		std::vector<float> shell_content(number_of_proximity_ranges);

		unsigned int proximity_range_index = 0;
		for (unsigned int i=0;i<unique_distances.size();i++)
		{
			if (proximity_range_index < number_of_proximity_ranges)
			{
				if (unique_distances[i] > proximity_ranges_ends[proximity_range_index])
				{
					proximity_range_index += 1;				
				}
		
				shell_content[proximity_range_index] += voxel_content[i];	

			}	
		}

		CPPUNIT_ASSERT_DOUBLES_EQUAL(3,shell_content[0],0.01);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(500,shell_content[1],0.01);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(60,shell_content[2],0.01);

		// seems to work as expected
	}

	void testOpenVDB_RearrangeProximities()
	{
		
		const unsigned int number_of_proximity_ranges = 5;
		float shell_width = 0.5;

		std::vector<float> proximity_ranges_ends = {-1,-0.5,0,0.5,1};

		std::vector<float> proximity_ranges_plotting(number_of_proximity_ranges);

		for (unsigned int i=0;i<proximity_ranges_plotting.size();i++)
		{
			if (proximity_ranges_ends[i] <= 0)
			{
				proximity_ranges_plotting[i] = proximity_ranges_ends[i] - shell_width;
			}			
			else
			{
				proximity_ranges_plotting[i] = proximity_ranges_ends[i];		
			}			
		}		

		for (unsigned int i=0;i<proximity_ranges_plotting.size();i++)
		{		
			std::cout << " proximity_ranges_plotting [i] " << " = " << proximity_ranges_plotting[i] << std::endl;	
		}

		std::vector<float> assert_ranges_ends = {-1.5,-1,-0.5,0.5,1};

		for (unsigned int i=0;i<proximity_ranges_plotting.size();i++)
		{		
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_ranges_ends[i], proximity_ranges_plotting[i], 0.01);	
		}

	}


	void testOpenVDB_VoxelIntersections()
	{

		// extract mesh
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);
		grid = createBlock(2,1);	
		std::vector<openvdb::Vec3s> points;
		std::vector<openvdb::Vec3I> triangles;
		std::vector<openvdb::Vec4I> quads;

		float isovalue=0.5;
		float adaptivity=0;
		openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, triangles, quads, isovalue, adaptivity);

		float in_bandwidth = 2;
		float ex_bandwidth = 2;

		float voxelsize_levelset = 0.5;

		// signed distance field
		openvdb::math::Transform::Ptr trans = openvdb::math::Transform::createLinearTransform(voxelsize_levelset);
		openvdb::FloatGrid::Ptr sdf = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(*trans, points, triangles, quads, ex_bandwidth, in_bandwidth);


	}

	void testOpenVDB_SetProximityRanges()
	{
		// the narrowband has to be set to the lowest proximity range end + its other half shell
		// with a given shell with of 1 and a max distance of 2 i want to end up with
		// the proximity ends of -1.5,-0.5,0.5,1,5 -1.5, 2.5 
		// so the narrowband has to include the distance -2.5 to 2.5
		// the centers lie at -2,-1,0,1,2

		std::vector<float> proximity_ranges_ends;
		std::vector<float> proximity_ranges_centers;

		const unsigned int rescue_counter = 10000;
		float max_distance = 2;
		float shell_width = 1;
		float current_end = shell_width/2;
		// 1st entries
		proximity_ranges_ends.push_back(current_end);
		proximity_ranges_ends.push_back(-current_end);

		unsigned int counter = 0;
		while (current_end < max_distance)
		{	
			proximity_ranges_ends.push_back(current_end + shell_width);
			proximity_ranges_ends.push_back(-(current_end + shell_width));
			current_end += shell_width;
			counter += 1;
			if (counter > rescue_counter)
			{
				break;			
			}
		}

		// now sort the ends
		std::sort(proximity_ranges_ends.begin(), proximity_ranges_ends.end());


		counter = 0;
		float current_center = 0;
		// 1st entry
		proximity_ranges_centers.push_back(current_center);
		while (current_center < max_distance)
		{	
			proximity_ranges_centers.push_back(current_center + shell_width);
			proximity_ranges_centers.push_back(-(current_center + shell_width));
			current_center += shell_width;
			counter += 1;
			if (counter > rescue_counter)
			{
				break;			
			}
		}

		// now sort the centers
		std::sort(proximity_ranges_centers.begin(), proximity_ranges_centers.end());
	
		std::vector<float> assert_ranges_ends = {-2.5,-1.5,-0.5,0.5,1.5,2.5};
		std::vector<float> assert_ranges_centers = {-2,-1,0,1,2};

		for (unsigned int i=0;i<assert_ranges_ends.size();i++)
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_ranges_ends[i], proximity_ranges_ends[i], 0.01);
		}

		for (unsigned int i=0;i<assert_ranges_centers.size();i++)
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_ranges_centers[i], proximity_ranges_centers[i], 0.01);
		}
	}

	void testOpenVDB_ShellCheck()
	{ 

		std::vector<float> proximity_ranges_limits {-2.5,-1.5,-0.5,0.5,1.5,2.5};
		std::vector<float> assert_ranges_centers {-2,-1,0,1,2};
		std::vector<float> shell_content(assert_ranges_centers.size());

		std::vector<float> unique_distances {-4,-2.3,-2.0,-1.4,-1.1,-0.3,0.2,1.1,1.4,2.1,2.3,3.0};
		std::vector<float> voxel_content;
		for (unsigned int i=0;i<unique_distances.size();i++)
		{
			voxel_content.push_back(1);
		}

		unsigned int proximity_range_index = 0;
		for (unsigned int i=0;i<unique_distances.size();i++)
		{
			
			if ((unique_distances[i] >= proximity_ranges_limits[proximity_ranges_limits.size()-1]))
			{
				break;
			}
			
			
			if ((unique_distances[i] >= proximity_ranges_limits[proximity_range_index]))
			{

				if (unique_distances[i] > proximity_ranges_limits[proximity_range_index+1])
				{
					proximity_range_index += 1;				
				}
				shell_content[proximity_range_index] += voxel_content[i];	

			}	
		}
		for (unsigned int i=0;i<shell_content.size();i++)
		{
			std::cout << " shell_content [i] " << " = " << shell_content[i] << std::endl;
			CPPUNIT_ASSERT_DOUBLES_EQUAL(2,shell_content[i],0.01);
		}


	}

};

#endif
