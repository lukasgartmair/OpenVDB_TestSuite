
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
 
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test1 - Test the Test itsself",
				&TestOpenVDB::testOpenVDB_TestTheTest ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test2 - Division of Data small by big",
				&TestOpenVDB::testOpenVDB_DivisionOfData1 ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test3 - Division of Data big by small",
				&TestOpenVDB::testOpenVDB_DivisionOfData2 ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestOpenVDB>("Test4 - Division of Data misc",
				&TestOpenVDB::testOpenVDB_DivisionOfData3 ));

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
	
		// this test divides the small sphere by the big one
		// that means there should be no division by zero

		openvdb::FloatGrid::Ptr sphere_small;
		openvdb::FloatGrid::Ptr sphere_big;
		// radius, value
		float radius_small = 10;
		float radius_big = 20;
		
		float value_small = 1;
		float value_big = 1;
		sphere_small = createSphere(radius_small,value_small);
		sphere_big = createSphere(radius_big, value_big);
		
		// composite.h operations modify the first grid and leave the second grid emtpy
		// compute a = a / b
		openvdb::tools::compDiv(*sphere_small, *sphere_big);
		
		bool check_for_nan = false;
		int number_of_nans = 0;
		
		for (openvdb::FloatGrid::ValueAllIter iter = sphere_small->beginValueAll(); iter; ++iter)
		{
		    if (std::isnan(iter.getValue()) == true)
		{
			check_for_nan = true;
			number_of_nans += 1;		    
		}
		}

		//CPPUNIT_ASSERT( check_for_nan == false);
		std::cout << "number of nans" << " = " << number_of_nans << std::endl;
		CPPUNIT_ASSERT( number_of_nans == 0);
		
	}
	
	void testOpenVDB_DivisionOfData2() {
	
		// this test divides the big sphere by the small one
		// that means there should be division by zero

		openvdb::FloatGrid::Ptr sphere_small;
		openvdb::FloatGrid::Ptr sphere_big;
		// radius, value
		float radius_small = 10;
		float radius_big = 20;
		
		float value_small = 1;
		float value_big = 1;
		sphere_small = createSphere(radius_small,value_small);
		sphere_big = createSphere(radius_big, value_big);
		
		// composite.h operations modify the first grid and leave the second grid emtpy
		// compute a = a / b
		openvdb::tools::compDiv(*sphere_big, *sphere_small);
		
		bool check_for_nan = false;
		int number_of_nans = 0;
		
		for (openvdb::FloatGrid::ValueAllIter iter = sphere_big->beginValueAll(); iter; ++iter)
		{
		    if (std::isnan(iter.getValue()) == true)
		{
			check_for_nan = true;	
			number_of_nans += 1;	    
		}
		}

		//CPPUNIT_ASSERT( check_for_nan == true);
		std::cout << "number of nans" << " = " << number_of_nans << std::endl;
		CPPUNIT_ASSERT( number_of_nans > 0);

	}
	
	void testOpenVDB_DivisionOfData3() {

		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(/*background value=*/0);
		openvdb::FloatGrid::Ptr grid_of_zeros = openvdb::FloatGrid::create(/*background value=*/0);
		
		// set few points in grid
		
		openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
		openvdb::Coord ijk(0,0,0);
		int value = 0;
		accessor.setValue(ijk, value);

		// composite.h operations modify the first grid and leave the second grid emtpy
		// compute a = a / b
		openvdb::tools::compDiv(*grid, *grid_of_zeros);
		
		bool check_for_nan = false;
		int number_of_nans = 0;
		
		for (openvdb::FloatGrid::ValueAllIter iter = grid->beginValueAll(); iter; ++iter)
		{
		    if (std::isnan(iter.getValue()) == true)
		{
			check_for_nan = true;	
			number_of_nans += 1;	    
		}
		}

		//CPPUNIT_ASSERT( check_for_nan == true);
		std::cout << "number of nans" << " = " << number_of_nans << std::endl;
		CPPUNIT_ASSERT( number_of_nans == 1);

	}
	
	
	

};
