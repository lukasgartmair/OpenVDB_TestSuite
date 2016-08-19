
#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include "vdb_functions.h"
 
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
		CPPUNIT_ASSERT( u == 3);
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

		std::cout << "number of nfs" << " = " << number_of_nfs << std::endl;
		CPPUNIT_ASSERT( number_of_nfs == 0);
		
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
		
		std::cout << "number of nfs" << " = " << number_of_nfs << std::endl;
		float minVal = 0.0;
		float maxVal = 0.0;
		big_grid->evalMinMax(minVal,maxVal);
		std::cout << " eval min max big grid" << " = " << minVal << " , " << maxVal << std::endl;
		std::cout << " active voxel count big grid" << " = " << big_grid->activeVoxelCount() << std::endl;
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

		std::cout << "number of nfs" << " = " << number_of_nfs << std::endl; 
		int assert_number_of_nfs = numerator_grid->activeVoxelCount() - 1;
		CPPUNIT_ASSERT( number_of_nfs == assert_number_of_nfs);
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

		std::cout << "number of nfs" << " = " << number_of_nfs << std::endl;
		CPPUNIT_ASSERT( number_of_nfs == 0);
	}
	
	
	void testOpenVDB_VolumeToMeshDMC()
	{
	
		openvdb::initialize();
		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(/*background value=*/0);
		// radius, value
		grid = createBlock(10,1);
		
		std::vector<openvdb::Vec3s> vertices;
		// change the isovalue with no adaptivity
		vertices = volumeToMeshVertices(grid, 0.5, 0);
		CPPUNIT_ASSERT(vertices.size() == 2402);
		
		vertices = volumeToMeshVertices(grid, 0.01, 0);
		CPPUNIT_ASSERT(vertices.size() == 2402);
		
		vertices = volumeToMeshVertices(grid, 0.99, 0);
		CPPUNIT_ASSERT(vertices.size() == 2402);
		
		// change the adaptivity with constant isovalue
		vertices = volumeToMeshVertices(grid, 0.5, 0.25);
		CPPUNIT_ASSERT(vertices.size() == 566);
		
		vertices = volumeToMeshVertices(grid, 0.5, 0.5);
		CPPUNIT_ASSERT(vertices.size() == 422);
		
		vertices = volumeToMeshVertices(grid, 0.5, 1);
		CPPUNIT_ASSERT(vertices.size() == 56);

	}
	
	void testOpenVDB_RoundUp()
	{
	// check the helper function round which is used to determine the voxel indices
	// it has certainly to be discussed whether this behaviour is suitable to fill
	// the ion grids !
	
		int rounded_result = 0;
		rounded_result = roundUp(0.0, 5);
		std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT(rounded_result == 0);
		
		
		rounded_result = 0;
		rounded_result = roundUp(0.5, 5);
		std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT(rounded_result == 0);
		
		rounded_result = 0;
		rounded_result = roundUp(0.99, 5);
		std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT(rounded_result == 0);
		
		rounded_result = 0;
		rounded_result = roundUp(1.0, 5);
		std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT(rounded_result == 5);
		
		rounded_result = 0;
		rounded_result = roundUp(4.9, 5);
		std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT(rounded_result == 5);

		rounded_result = 0;
		rounded_result = roundUp(5.1, 5);
		std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT(rounded_result == 5);
		
		rounded_result = 0;
		rounded_result = roundUp(5.9, 5);
		std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT(rounded_result == 5);

		rounded_result = 0;
		rounded_result = roundUp(6.0, 5);
		std::cout << "rounded result" << " = " << rounded_result << std::endl;
		CPPUNIT_ASSERT(rounded_result == 10);

		
	}
	

};
