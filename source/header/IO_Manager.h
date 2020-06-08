/**
* @file
* @author Zhengliang Liu
* @brief This class used to manager IO classes
*.
*/
#ifndef IO_MANAGER_H
#define IO_MANAGER_H
#include "General.h"
//#include "Tecplot.h"
//#include "CGNS.h"
//#include "VTK.h"
#include "HDF5H.h"
class IO_Manager
{
public:
	static IO_Manager* pointer_me;    ///< static pointer point to the IO manager
	unsigned int method;              ///< format used to write flowfield data
	std::string outfile;              ///< name of the file to write flowfield data
	std::vector<unsigned int> vlevel; ///< refiement levels of blocks for flowfield output
public:
	void control();                   ///< the maind funtion of the IO manager
private:
//	Tecplot io_tecplot;
//	CGNS io_cgns;
//	VTK io_vtk;
	HDF5H io_hdf5;
};
#endif