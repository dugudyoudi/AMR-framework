/**
* @file
* @author Zhengliang Liu
* @brief mah.
* @note .
*/
#include "IO_Manager.h"
void IO_Manager::control()
{

	if (method == 1)
	{
		//io_tecplot.grid_vlevel = vlevel;
		//io_tecplot.write_manager();	
	}
	else if (method == 2)
	{
		//io_cgns.grid_vlevel = vlevel;
		//io_cgns.write_manager(outfile);
	}
	else if (method == 3)
	{
		//io_vtk.grid_vlevel = vlevel;
		//io_vtk.write_manager(outfile);
	}
	else if (method == 4)
	{
		io_hdf5.grid_vlevel = vlevel;
		io_hdf5.write_manager(outfile);
	}
}


