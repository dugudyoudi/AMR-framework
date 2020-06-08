/**
* @file
* @author Zhengliang Liu
* @brief Include stardard libraries, headers and functions.
* @note 1) include libraries; 2) define macros; 3) define inline funciotns.
*/

#ifndef GENERAL_H
#define GENERAL_H

// libraries
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cstdlib>
#include <cmath>
#include <bitset>
#include <direct.h>

typedef double D_real;
typedef int D_int;
typedef unsigned long long D_uint;

// detect system to choose approriate libraries for mikdir
#ifdef _WIN32   
#include <io.h> 
#include <direct.h>  
#include <Windows.h>
#elif __linux__ 
#include <unistd.h>  
#include <sys/types.h>  
#include <sys/stat.h> 
#endif

// timer
#ifdef __unix__
#include <ctime>
class Timer
{
public:
	Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }

	double elapsed() {
		clock_gettime(CLOCK_REALTIME, &end_);
		return end_.tv_sec - beg_.tv_sec +
			(end_.tv_nsec - beg_.tv_nsec) / 1000000000.;
	}

	void reset() { clock_gettime(CLOCK_REALTIME, &beg_); }

private:
	timespec beg_, end_;
};
#else
#include <chrono>
class Timer
{
public:
	Timer() : beg_(clock_::now()) {}
	void reset() { beg_ = clock_::now(); }
	double elapsed() const {
		return std::chrono::duration_cast<second_>
			(clock_::now() - beg_).count();
	}

private:
	typedef std::chrono::high_resolution_clock clock_;
	typedef std::chrono::duration<double, std::ratio<1> > second_;
	std::chrono::time_point<clock_> beg_;
};
#endif

// headers
#include "Constants.h"

// define expressions
#define SQ(x) ((x) * (x)) ///< Square

/**
* @brief The structure store functions to store log file
*/
struct Log_function
{
public:
	static std::ofstream* logfile;			///< pointer point to log file
};

/**
* @brief function to write normal information to the logfile.
* @param[in]  msg        information write to the logfile.
* @param[in]  logfile   pointer to the logfile address.
*/
inline void log_infor(const std::string &msg, std::ofstream *logfile)
{

	*logfile << "Infor: " << msg << std::endl;
	std::cout << "Infor:" << msg << std::endl;
	//*logfile << msg << std::endl;
};

/**
* @brief function to write warning information to the logfile.
* @param[in]  msg        information write to the logfile.
* @param[in]  logfile   pointer to the logfile address.
*/
inline  void log_warning(const std::string &msg, std::ofstream *logfile)
{

	*logfile << "Warning: " << msg << std::endl;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN);
	std::cout << "Warning:" << msg << std::endl;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
};

/**
* @brief function to write error information to the logfile.
* @param[in]  msg        information write to the logfile.
* @param[in]  logfile   pointer to the logfile address.
*/
inline  void log_error(const std::string &msg, std::ofstream *logfile)
{

	*logfile << "Error: " << msg << std::endl;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | FOREGROUND_RED);
	std::cout << "Error:" << msg << std::endl;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
	exit(0);
};

/**
* @brief function to calculate 2 power N.
* @param[in]  n        power exponent.
*/
inline D_uint two_power_n(D_uint n)
{
	D_uint two = 1;
	two  <<= n;
	return two;
};

#endif
