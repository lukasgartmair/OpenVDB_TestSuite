#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include "vdb_testsuite.h"
#include "vdb_functions.h"

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
