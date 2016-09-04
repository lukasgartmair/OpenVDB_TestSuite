#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include "TestSuiteOpenVDB.h"
#include "vdb_functions.h"


// compile: sudo g++ -I/usr/include/openvdb TestRunner.cpp -lcppunit -ltbb -o RunTests
//latest successful: sudo g++ -I/usr/include/openvdb ~/LpCVT/combinatorics/mesh.cpp TestRunner.cpp -lcppunit -ltbb -o RunTests
// sudo g++ -I/usr/include/openvdb vdb_test1.cpp TestRunner.cpp -lcppunit -lopenvdb -ltbb -lHalf -o RunTests
//sudo g++ -I/usr/include/openvdb vdb_functions.cpp TestRunner.cpp -lcppunit -lopenvdb -ltbb -lHalf -o RunTests

// sudo g++ -std=c++11 -I/usr/include/openvdb vdb_functions.cpp TestRunner.cpp -lcppunit -lopenvdb -ltbb -lHalf -o RunTests


using namespace std;
 
int main() {
	CppUnit::TextUi::TestRunner runner;
 
	cout << "Creating Test Suites:" << endl;
	runner.addTest(TestOpenVDB::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();
 
	return 0;
}
