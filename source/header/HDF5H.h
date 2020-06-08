/**
* @file
* @author Zhengliang Liu
* @brief Define Tecplot Class.
*/
#ifndef HDF5_H
#define HDF5_H
#include "General.h"
#include "hdf5.h"
#include "Grid_Manager.h"
#include "Solid_Manager.h"

class HDF5H
{
	friend class IO_Manager;
private:
	std::vector<unsigned int> grid_vlevel;
	std::vector<D_uint> grid_npoints, grid_nelements;
	std::array<std::vector<D_uint>, C_overlap * 2 + 1> bouandry_fine2coarse_npoints, bouandry_fine2coarse_nelements;
	std::array <std::array<D_uint, 2>, C_DIMS> domain_boundary_npoints, domain_boundary_nelements;
private:
	int write_manager(std::string outfile);
//	void compute_points(unsigned int ilevel);
	void write_hdf5(std::string outfile);
	template <class T_grptr> void write_grid(T_grptr &grid_ptr, int level_index, hid_t &group0_id, std::string &ini_group);
	void write_xdmf(std::string outfile);

	void write_solid_polyline(const hid_t file_id);
	void write_xdmf_solid_ployline(std::ofstream &out, const std::string &outfile);

	void write_numerical_boundary(const hid_t file_id);
	void write_xdmf_numerical_boundary(std::ofstream &out, const std::string &outfile);

	void write_domain_boundary(const hid_t file_id);
	void write_xdmf_domain_boundary(std::ofstream &out, const std::string &outfile);

#if (C_DEBUG == 1)
	void write_output_points_update(const hid_t file_id);
	void write_xdm_output_points_update(std::ofstream &out, const std::string &outfile);
#endif
};

#endif