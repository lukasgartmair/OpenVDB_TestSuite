
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

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test16 - triangle areas ",
				&TestOpenVDB::testOpenVDB_TriangleAreas));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test17 - export triangles to obj ",
				&TestOpenVDB::testOpenVDB_ExportTriangles));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test18 - active / inactive voxels ",
				&TestOpenVDB::testOpenVDB_ActiveInactiveVoxels));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test19 - get voxel coords ",
				&TestOpenVDB::testOpenVDB_GetVoxelCoords));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test20 - find coord in vector ",
				&TestOpenVDB::testOpenVDB_FindCoordInVector));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test21 - signed distance field",
				&TestOpenVDB::testOpenVDB_SignedDistanceField));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test22 - get active voxels mask",
				&TestOpenVDB::testOpenVDB_ActiveVoxelsMask));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test23 - sdf test",
				&TestOpenVDB::testOpenVDB_SDFTest));




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


		std::vector<std::vector<float> > triangles_from_splitted_quads;
		triangles_from_splitted_quads = splitQuadsToTriangles(points, quads);
		
		std::vector<std::vector<float> > triangles_combined;
		
		triangles_combined = concatenateTriangleVectors(triangles, triangles_from_splitted_quads);


		std::vector<float> triangle_areas(triangles_combined.size());
		triangle_areas = ComputeTriangleAreas(points, triangles_combined);
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

		std::vector<std::vector<float> > triangles_from_splitted_quads;
		triangles_from_splitted_quads = splitQuadsToTriangles(points, quads);
		
		std::vector<std::vector<float> > triangles_combined;
		
		triangles_combined = concatenateTriangleVectors(triangles, triangles_from_splitted_quads);

		
		ExportTriangleMeshAsObj(points, triangles_combined);

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

		const int xyzs = 3;

		std::vector<std::vector<float> > active_voxel_indices(grid->activeVoxelCount(), std::vector<float>(xyzs));
		openvdb::Coord hkl;

		int voxel_counter = 0;
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

		const int xyzs = 3;

		std::vector<std::vector<float> > active_voxel_indices(grid->activeVoxelCount(), std::vector<float>(xyzs));
		openvdb::Coord hkl;

		int voxel_counter = 0;
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

		// signed distance field
		openvdb::FloatGrid::Ptr sdf = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(openvdb::math::Transform(), points, triangles, quads, ex_bandwidth, in_bandwidth);

		int active_voxels_meshgrid = grid->activeVoxelCount();
		int active_voxels_sdf = sdf->activeVoxelCount();

		//std::cout << " active_voxels_meshgrid" << " = " << active_voxels_meshgrid << std::endl;
		//std::cout << " active_voxels_sdf" << " = " << active_voxels_sdf << std::endl;

		bool narrowband = false;
		if (active_voxels_meshgrid < active_voxels_sdf)	
		{
			narrowband = true;
		}
		CPPUNIT_ASSERT_DOUBLES_EQUAL(true, narrowband,0.01);
	
	}


	void testOpenVDB_ActiveVoxelsMask()
	{

		// only available in 3.2.0 -> see releasenotes
	}

	void testOpenVDB_SDFTest()
	{
		// extract mesh
		openvdb::initialize();

		// is the distance from sdf in voxel units or in real world units and
		// is the higher resoltution really passed, when the mesh extracted from the lower
		// resolution is passed in
		// and third how is the distance calculated

		float initial_voxelsize = 1.0;		

		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0);
		grid = createBlock(1,1);

		grid->setTransform(openvdb::math::Transform::createLinearTransform(initial_voxelsize));

		std::vector<openvdb::Vec3s> points;
		std::vector<openvdb::Vec3I> triangles;
		std::vector<openvdb::Vec4I> quads;

		float isovalue=0.5;
		float adaptivity=0;
		openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*grid, points, triangles, quads, isovalue, adaptivity);

		float voxelsize_levelset1 = 1.0;
		float voxelsize_levelset2 = 0.5;

		float in_bandwidth = 60;
		float ex_bandwidth = 60;

		openvdb::FloatGrid::Ptr sdf1 = openvdb::FloatGrid::create(0.0); 
		openvdb::FloatGrid::Ptr sdf2 = openvdb::FloatGrid::create(0.0); 

		// signed distance field
		sdf1 = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(openvdb::math::Transform(), points, triangles, quads, ex_bandwidth, in_bandwidth);
		sdf2 = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(openvdb::math::Transform(), points, triangles, quads, ex_bandwidth, in_bandwidth);

		sdf1->setTransform(openvdb::math::Transform::createLinearTransform(voxelsize_levelset1));
		sdf2->setTransform(openvdb::math::Transform::createLinearTransform(voxelsize_levelset2));

		std::cout << " active_voxels_sdf 1 " << " = " << sdf1->activeVoxelCount() << std::endl;
		std::cout << " active_voxels_sdf 2 " << " = " << sdf2->activeVoxelCount() << std::endl;

		float minVal = 0.0;
		float maxVal = 0.0;
		sdf1->evalMinMax(minVal,maxVal);
		std::cout << " eval min max sdf1" << " = " << minVal << " , " << maxVal << std::endl;
		sdf2->evalMinMax(minVal,maxVal);
		std::cout << " eval min max sdf2" << " = " << minVal << " , " << maxVal << std::endl;

		// are the distances just in voxel units or is the transform not passed correctly ?!
		// maybe try index world here
		
		openvdb::math::Transform::Ptr linearTransform1 =
		    openvdb::math::Transform::createLinearTransform(voxelsize_levelset1);

		openvdb::math::Transform::Ptr linearTransform2 =
		    openvdb::math::Transform::createLinearTransform(voxelsize_levelset2);

		openvdb::Coord ijk(6,6,6);
		openvdb::Vec3d worldSpacePoint1 = linearTransform1->indexToWorld(ijk);
		openvdb::Vec3d worldSpacePoint2 = linearTransform2->indexToWorld(ijk);

		std::cout << " worldSpacePoint1 " << " = " << worldSpacePoint1 << std::endl;
		std::cout << " worldSpacePoint2 " << " = " << worldSpacePoint2 << std::endl;
		
		// that's working as expected but what about the distances?
		// how can i get a distance in world space		

		// or do i just need to double the bandwidth when using e.g. half the initial voxel size for a higher resolution
		// in order to get the same spatial information
		// is the translation implied in the contribution function when the voxel size is passed ? 


		// passing two different treansforms on sdf and check whether the distances are still the same or not
		// the case i want to achieve with lower voxelsize in the sdf i want smaller discrete distance steps -
		// not just more voxels with the same distance the big voxel would have


		// now lets try what happens if the voxel distances in voxel units are converted to world space by 
		// bringing the voxelsize in

		// voxel distance of ijk and hkl
		openvdb::FloatGrid::Accessor sdf2_accessor = sdf2->getAccessor();
		openvdb::FloatGrid::Accessor sdf1_accessor = sdf1->getAccessor();
		std::cout << " sdf1 ijk " << " = " << sdf1_accessor.getValue(ijk) << std::endl;
		std::cout << " sdf2 ijk " << " = " << sdf2_accessor.getValue(ijk) << std::endl;

		openvdb::Coord hkl(5,5,5);
		std::cout << " sdf1 ijk " << " = " << sdf1_accessor.getValue(hkl) << std::endl;
		std::cout << " sdf2 ijk " << " = " << sdf2_accessor.getValue(hkl) << std::endl;

		//real distance of ijk and hkl
		// is multiplication with the voxelsize the right conversion?
		// to my understanding a voxel (5,5,5) with distance 10 voxelunits and a size of one should have a distance of 
		// 10 in the end. Whereas the voxel in the grid with the lower resolution (5,5,5) should have the same distance
		// in voxel units of course but with a voxel size of only 0.5 i.e. the half of the first one the real distance should also
		// be half the size of the first one. So 10 times 0.5 equals 5 and is half the distance of the first in real space?!

		std::cout << " sdf1 ijk real" << " = " << sdf1_accessor.getValue(ijk) * voxelsize_levelset1 << std::endl;
		std::cout << " sdf2 ijk real " << " = " << sdf2_accessor.getValue(ijk) * voxelsize_levelset2 << std::endl;

		std::cout << " sdf1 ijk real " << " = " << sdf1_accessor.getValue(hkl) * voxelsize_levelset1 << std::endl;
		std::cout << " sdf2 ijk real " << " = " << sdf2_accessor.getValue(hkl) * voxelsize_levelset2 << std::endl;

	}







































};
