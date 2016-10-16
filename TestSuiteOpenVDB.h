
#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include "vdb_functions.h"

#include <algorithm>    // std::min
 
class TestOpenVDB : public CppUnit::TestFixture {


public:
 
	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestOpenVDB");
 
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test0 - Test the Test itsself",
				&TestOpenVDB::testOpenVDB_TestTheTest ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test1 - Division of Data1 small by big",
				&TestOpenVDB::testOpenVDB_DivisionOfData1 ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test2 - Division of Data2 big by small",
				&TestOpenVDB::testOpenVDB_DivisionOfData2 ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test3 - Division of Data3 misc",
				&TestOpenVDB::testOpenVDB_DivisionOfData3 ));
		
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test4 - Division of Data4 misc",
				&TestOpenVDB::testOpenVDB_DivisionOfData4 ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test5 - Volume to mesh DMC",
				&TestOpenVDB::testOpenVDB_VolumeToMeshDMC ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test6 - round up ",
				&TestOpenVDB::testOpenVDB_RoundUp ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test7 - get voxel index ",
				&TestOpenVDB::testOpenVDB_GetVoxelIndex ));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test8 - conversion openvdb vec3s to std::vec ",
				&TestOpenVDB::testOpenVDB_VectorConversion ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test9 - split openvdb quads to triangles ",
				&TestOpenVDB::testOpenVDB_SplitQuadsToTriangles ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test10 - concatenate triangle lists ",
				&TestOpenVDB::testOpenVDB_ConcatenateTriangles ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test11 - increase triangle vertex indices by n ",
				&TestOpenVDB::testOpenVDB_IncreaseTrianglesVertexIndices ));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test12 - decrease triangle vertex indices by n ",
				&TestOpenVDB::testOpenVDB_DecreaseTrianglesVertexIndices ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test13 - compute vertex normals ",
				&TestOpenVDB::testOpenVDB_ComputeVertexNormals ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test14 - compute triangle normals ",
				&TestOpenVDB::testOpenVDB_ComputeTriangleNormals ));


		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test15 - bounding box ",
				&TestOpenVDB::testOpenVDB_BoundingBox));


		return suiteOfTests;
	}
 
	/// Setup method
	void setUp() {}
 
	/// Teardown method
	void tearDown() {}
 
protected:
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
		// so i have to check both cases http://stackoverflow.com/questions/4095337/how-to-check-for-inf-and-or-nf-in-a-double-variable
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
		int xyzs = 3;
		for (int i=0;i<vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
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
 		
		for (int i=0;i<vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
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
		result_coord = GetVoxelIndex(&test_coord, voxelsize);
		coord assert_coord = {0,0,0};
		CPPUNIT_ASSERT_EQUAL(result_coord.x, assert_coord.x);
		CPPUNIT_ASSERT_EQUAL(result_coord.y, assert_coord.y);
		CPPUNIT_ASSERT_EQUAL(result_coord.z, assert_coord.z);	
		
		test_coord = {12.5,-22,-3};
		result_coord = {0,0,0};
		voxelsize = 3;
		result_coord = GetVoxelIndex(&test_coord, voxelsize);
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
		result_coord = GetVoxelIndex(&test_coord, voxelsize);
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
		
	  	int xyzs = 3;
		std::vector<std::vector<float> > standard_points(vertices.size(), std::vector<float>(xyzs));

		standard_points = convertOpenVDBVectorToStandardVector(vertices);
		
		for (int i=0;i<vertices.size();i++)
		{
			CPPUNIT_ASSERT_EQUAL(vertices[i].x(), standard_points[i][0]);
			CPPUNIT_ASSERT_EQUAL(vertices[i].y(), standard_points[i][1]);
			CPPUNIT_ASSERT_EQUAL(vertices[i].z(), standard_points[i][2]);
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

		std::vector<std::vector<float> > triangles_from_splitted_quads;

		//triangles_from_splitted_quads = splitQuadsToTriangles(points, quads);
		
		//std::cout << triangles_from_splitted_quads.size() << std::endl; 
		
		//CPPUNIT_ASSERT_EQUAL(quads.size()*2,triangles_from_splitted_quads.size());
		
		// faces outputs from openvdb volume to mesh starting from zero!
		int minimum  = 20;
		int minimum_tris = 20;
		for (int i=0;i<triangles.size();i++)
		{
			for (int j=0;j<3;j++)
			{
				if (triangles[i][j] < minimum)
				{
					minimum = triangles[i][j];
				}
			}
		}
		
		int minimum_quads = 20;
		for (int i=0;i<quads.size();i++)
		{
			for (int j=0;j<3;j++)
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

		const int assert_faces_minimum = 0;
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
		int tri_size = triangles_from_splitted_quads.size();
		CPPUNIT_ASSERT_EQUAL(2,tri_size);
		
		//check the contents
		
		int xyzs = 3;
		std::vector<std::vector<float> > assert_triangles(tri_size, std::vector<float>(xyzs));
		
		assert_triangles[0][0] = 0; 
		assert_triangles[0][1] = 1; 
		assert_triangles[0][2] = 2; 
		assert_triangles[1][0] = 0; 
		assert_triangles[1][1] = 2; 
		assert_triangles[1][2] = 3; 	
		
		for (int i=0;i<tri_size;i++)
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

		std::vector<std::vector<float> > triangles_from_splitted_quads;
		triangles_from_splitted_quads = splitQuadsToTriangles(points, quads);
		
		std::vector<std::vector<float> > triangles_combined;
		
		triangles_combined = concatenateTriangleVectors(triangles, triangles_from_splitted_quads);
		
		//std::cout << triangles_combined[triangles.size()][0] << std::endl; 
		//std::cout << triangles_combined[triangles.size()][1] << std::endl; 
		//std::cout << triangles_combined[triangles.size()][2] << std::endl; 

		// check the amount
		int tri_size = triangles_combined.size();
		int assert_size = (quads.size()*2) + triangles.size();
		CPPUNIT_ASSERT_EQUAL(assert_size, tri_size);

		//check the contents
		
		for (int i=0;i<triangles.size();i++)
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles[i][0], triangles_combined[i][0],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles[i][1], triangles_combined[i][1],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles[i][2], triangles_combined[i][2],0.01);
		}
		
		for (int i=triangles.size();i<2*quads.size()+triangles.size();i++)
		{
			int shifted_index = i - triangles.size();
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles_from_splitted_quads[shifted_index][0], triangles_combined[i][0],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles_from_splitted_quads[shifted_index][1], triangles_combined[i][1],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(triangles_from_splitted_quads[shifted_index][2], triangles_combined[i][2],0.01);
		}

	}
	void testOpenVDB_IncreaseTrianglesVertexIndices()
	{
		int tri_size = 2;
		int xyzs = 3;
		std::vector<std::vector<float> > triangles(tri_size, std::vector<float>(xyzs));
	
		triangles[0][0] = 0; 
		triangles[0][1] = 1; 
		triangles[0][2] = 2; 
		triangles[1][0] = 0; 
		triangles[1][1] = 2; 
		triangles[1][2] = 3; 	

		std::vector<std::vector<float> > assert_triangles(tri_size, std::vector<float>(xyzs));
		
		assert_triangles[0][0] = 1; 
		assert_triangles[0][1] = 2; 
		assert_triangles[0][2] = 3; 
		assert_triangles[1][0] = 1; 
		assert_triangles[1][1] = 3; 
		assert_triangles[1][2] = 4; 
		
		std::vector<std::vector<float> > result_triangles(tri_size, std::vector<float>(xyzs));
		int N = 1;
		result_triangles = IncreaseTriangleVertexIndicesByN(triangles, N);
		
		//check the contents
		for (int i=0;i<tri_size;i++)
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][0], result_triangles[i][0],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][1], result_triangles[i][1],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][2], result_triangles[i][2],0.01);
		}
	}
	
	void testOpenVDB_DecreaseTrianglesVertexIndices()
	{
		int tri_size = 2;
		int xyzs = 3;
		std::vector<std::vector<float> > triangles(tri_size, std::vector<float>(xyzs));
	
		triangles[0][0] = 1; 
		triangles[0][1] = 2; 
		triangles[0][2] = 3; 
		triangles[1][0] = 1; 
		triangles[1][1] = 3; 
		triangles[1][2] = 4; 	

		std::vector<std::vector<float> > assert_triangles(tri_size, std::vector<float>(xyzs));
		
		assert_triangles[0][0] = 0; 
		assert_triangles[0][1] = 1; 
		assert_triangles[0][2] = 2; 
		assert_triangles[1][0] = 0; 
		assert_triangles[1][1] = 2; 
		assert_triangles[1][2] = 3; 
		
		std::vector<std::vector<float> > result_triangles(tri_size, std::vector<float>(xyzs));
		int N = 1;
		result_triangles = DecreaseTriangleVertexIndicesByN(triangles, N);
		
		//check the contents
		for (int i=0;i<tri_size;i++)
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][0], result_triangles[i][0],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][1], result_triangles[i][1],0.01);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_triangles[i][2], result_triangles[i][2],0.01);
		}
	}
	
	void testOpenVDB_ComputeVertexNormals()
	{

		int tri_size = 2;
		int xyzs = 3;
		std::vector<std::vector<float> > triangles(tri_size, std::vector<float>(xyzs));
		
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
		
		triangle_normals = ComputeTriangleNormalsVDB(points, triangles);
	
		std::vector<std::vector<float> > vertex_normals;
	
		vertex_normals = ComputeVertexNormals(triangles, points, triangle_normals);

	
		// check the amount
		int size_normals = 4;
		CPPUNIT_ASSERT_DOUBLES_EQUAL(size_normals, vertex_normals.size(),0.01);
		// check the content
		
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
		
	}


	void testOpenVDB_ComputeTriangleNormals()
	{
	
		int tri_size = 2;
		int xyzs = 3;
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
		
		std::vector<std::vector<float> > triangle(1, std::vector<float>(xyzs));
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
		
		normals = ComputeTriangleNormalsVDB(points, triangle);

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
		
		normals = ComputeTriangleNormals(std_points, triangle);

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

		std::cout << bounding_box.getStart() << std::endl;
            	std::cout << bounding_box.getStart()[0] << std::endl;

		std::cout << bounding_box.getEnd() << std::endl;
		std::cout << bounding_box.getCenter() << std::endl;
		std::cout << bounding_box.dim() << std::endl;
		

	}

















};
