/**
* @file
* @author Zhengliang Liu
* @brief Main function.
* @note .
*/
#include "General.h"
#include "Obj_Manager.h"


using namespace std;
ofstream* Log_function::logfile = nullptr;
Obj_Manager* Obj_Manager::pointer_me;

int main()
{
	Timer tmr;
	// Create application log file
	std::ofstream logfile;
	logfile.open("log", ios::out);
	Log_function::logfile = &logfile;
	if (!logfile.is_open()) {
		cout << "Can't open logfile" << endl;
	};

	Obj_Manager obj_manager;
	Obj_Manager::pointer_me = &obj_manager;

	double t0 = tmr.elapsed();

	obj_manager.initial();

	double t1 = tmr.elapsed();
	std::cout << "time for initialization: " << t1 - t0 << std::endl;

	obj_manager.time_marching_management();

	double t2 = tmr.elapsed();
	std::cout << "time for time marching: " << t2 - t1 << std::endl;

	obj_manager.output();

	logfile.close();
	return 0;
}