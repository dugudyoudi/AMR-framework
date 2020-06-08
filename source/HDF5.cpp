/**
* @file
* @author Zhengliang Liu
* @brief  Functions to write flow field in HDF5 format. XDMF is used to described the structure of the grid.
*/
#include "HDF5H.h"
#include "Morton_Assist.h"
bool flag_write_numerical_boundary = true;
bool flag_write_domain_boundary = true;

/**
* @brief      function to manage output in HDF5 and XDMF format.
* @param[in]  outfile      name for output file.
*/
int HDF5H::write_manager(std::string outfile)
{
	write_hdf5(outfile);
	write_xdmf(outfile);

	return 0;
}

/**
* @brief      function to write and call other functions to write output data in HDF5 file.
* @param[in]  outfile      name for output file.
*/
void HDF5H::write_hdf5(std::string outfile)
{
	// create file
	std::string outname;
	outname = outfile + ".h5";
	hid_t file_id, group0_id;
	file_id = H5Fcreate(outname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// write fuild
	std::string ini_group = "/Ini_grid";
	group0_id = H5Gcreate(file_id, ini_group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		
	for (int level_index = 0; level_index < grid_vlevel.size(); ++level_index)
	{
		if (grid_vlevel.at(level_index) == C_max_level)
		{
			write_grid(Grid_Manager::pointer_me->gr_inner, level_index, group0_id, ini_group);
		}
		else
		{
			write_grid(Grid_Manager::pointer_me->gr_NoIB.at(grid_vlevel.at(level_index)), level_index, group0_id, ini_group);
		}
	}
	H5Gclose(group0_id);

	// write node.flag to output file
	ini_group = "/Flag";
	group0_id = H5Gcreate(file_id, ini_group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	for (int level_index = 0; level_index < grid_vlevel.size(); ++level_index)
	{
		unsigned int ilevel = grid_vlevel.at(level_index);
		std::string block_group = "/block" + std::to_string(ilevel);
		hid_t dataspace_id, dataset_id;

		int *scalar = nullptr;
		scalar = new int[grid_npoints.at(level_index)];
		D_uint  icount = 0;
		if (grid_vlevel.at(level_index) == C_max_level)
		{
			for (D_mapIB::iterator iter = Grid_Manager::pointer_me->gr_inner.grid.begin(); iter != Grid_Manager::pointer_me->gr_inner.grid.end(); ++iter)
			{
				scalar[icount] = iter->second.flag;
				++icount;
			}
		}
		else
		{
			for (D_map::iterator iter = Grid_Manager::pointer_me->gr_NoIB.at(grid_vlevel.at(level_index)).grid.begin(); iter != Grid_Manager::pointer_me->gr_NoIB.at(grid_vlevel.at(level_index)).grid.end(); ++iter)
			{

				scalar[icount] = iter->second.flag;
				++icount;
			}
		}
		
		hid_t group_id = H5Gcreate(group0_id, (ini_group + block_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hsize_t dim[1] = { grid_npoints.at(level_index) };
		dataspace_id = H5Screate_simple(1, dim, NULL);
		dataset_id = H5Dcreate(group_id, (ini_group + block_group + "/flag").c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, scalar);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);

		delete[] scalar;
		H5Gclose(group_id);
	}
	H5Gclose(group0_id);

	// write solid
	if (Solid_Manager::pointer_me->numb_solids > 0)
	{
		write_solid_polyline(file_id);
	}

	// write numerical boundary
	if (flag_write_numerical_boundary)
	{
		write_numerical_boundary(file_id);
	}

	if (flag_write_domain_boundary)
	{
		write_domain_boundary(file_id);
	}

#if	(C_SOLID_BOUNDARY == 2)	
#if (C_DEBUG == 1)
	write_output_points_update(file_id);
#endif
#endif

	H5Fclose(file_id);

}

/**
* @brief      function to write infomation of coordinates and connectivities of nodes.
* @param[in]  grid_ptr    pointer point to the class storing grid information at ilevel.
* @param[in]  level_index    level of refinement at which grid information needs to be output.
* @param[in]   group0_id      group id for HDF output.
* @param[in]   ini_group      group name.
*/
template <class T_grptr>
void HDF5H::write_grid(T_grptr &grid_ptr, int level_index, hid_t &group0_id, std::string &ini_group)
{

	unsigned int ilevel = grid_vlevel.at(level_index);
	std::string block_group = "/block" + std::to_string(ilevel);
	hid_t group_id, dataspace_id, dataset_id;
	group_id = H5Gcreate(group0_id, (ini_group + block_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	D_uint icount_points = 0, icount_cells = 0;

	grid_npoints.push_back(static_cast <D_uint>(grid_ptr.grid.size()));

	D_mapint count;
	D_real *xyz = nullptr;
	xyz = new D_real[grid_npoints.at(level_index)*C_DIMS];

	D_morton ref_one = 0;
	hsize_t dims[2];
	dims[0] = grid_npoints.at(level_index);
	dims[1] = C_DIMS;

	ref_one.set(Morton_Assist::bit_otherlevel - ilevel * C_DIMS, true); // reference value which is 1 at the sepecified position for a given level

	// write coodinates of nodes
	for (auto iter = grid_ptr.grid.begin(); iter != grid_ptr.grid.end(); ++iter)
	{
		count.insert({ iter->first, icount_points });
		xyz[icount_points*C_DIMS] = 0;
		xyz[icount_points*C_DIMS + 1] = 0;
		D_morton key_current = iter->first;
#if (C_DIMS==2)
		D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
		D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
		Morton_Assist::pointer_me->compute_coordinate(iter->first, ilevel, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1]);
		xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
		xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
		D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);

		if ((grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()))
		{
			++icount_cells;
		}
#endif
#if (C_DIMS==3)
		xyz[icount_points*C_DIMS + 2] = 0;
		Morton_Assist::pointer_me->compute_coordinate(iter->first, ilevel, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1], xyz[icount_points*C_DIMS + 2]);
		xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
		xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
		xyz[icount_points * C_DIMS + 2] -= Solid_Manager::pointer_me->shape_offest_z0_grid;
		D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
		D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
		D_morton key_z1 = Morton_Assist::find_z1(key_current, ilevel);
		D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);
		D_morton key_xz1 = Morton_Assist::find_z1(key_x1, ilevel);
		D_morton key_yz1 = Morton_Assist::find_z1(key_y1, ilevel);
		D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, ilevel);

		bool bool_xyz = (grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_z1) != grid_ptr.grid.end())
			&& (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_xz1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_yz1) != grid_ptr.grid.end())
			&& (grid_ptr.grid.find(key_xyz1) != grid_ptr.grid.end());

		if (bool_xyz)
		{
			++icount_cells;
		}
#endif
		++icount_points;
	}
	grid_nelements.push_back(icount_cells);

	dataspace_id = H5Screate_simple(2, dims, NULL);
	dataset_id = H5Dcreate(group_id, (ini_group + block_group + "/xyz").c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);

	delete[] xyz;
	int *cell = nullptr;
	hsize_t dims2[2];
	dims2[0] = grid_nelements.at(level_index);
	dims2[1] = two_power_n(C_DIMS);
	cell = new int[static_cast<int> (dims2[0]) * static_cast<int> (dims2[1])];
	D_uint sum_numb = 0;

	ref_one.set(Morton_Assist::bit_otherlevel - ilevel * C_DIMS, true); // reference value which is 1 at the sepecified position for a given level

	// write connectivities of the nodes
	for (auto iter = grid_ptr.grid.begin(); iter != grid_ptr.grid.end(); ++iter)
	{
		D_morton key_current = iter->first;
#if (C_DIMS==2)
		D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
		D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
		D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);
		if ((grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()))
		{
			cell[sum_numb] = count[key_current];
			++sum_numb;
			cell[sum_numb] = count[key_x1];
			++sum_numb;
			cell[sum_numb] = count[key_xy1];
			++sum_numb;
			cell[sum_numb] = count[key_y1];
			++sum_numb;
		}
#endif
#if (C_DIMS==3)
		D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
		D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
		D_morton key_z1 = Morton_Assist::find_z1(key_current, ilevel);
		D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);
		D_morton key_xz1 = Morton_Assist::find_z1(key_x1, ilevel);
		D_morton key_yz1 = Morton_Assist::find_z1(key_y1, ilevel);
		D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, ilevel);

		bool bool_xyz = (grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_z1) != grid_ptr.grid.end())
			&& (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_xz1) != grid_ptr.grid.end()) && (grid_ptr.grid.find(key_yz1) != grid_ptr.grid.end())
			&& (grid_ptr.grid.find(key_xyz1) != grid_ptr.grid.end());

		if (bool_xyz)
		{
			D_morton key_xz1 = Morton_Assist::find_z1(key_x1, ilevel);
			D_morton key_yz1 = Morton_Assist::find_z1(key_y1, ilevel);

			cell[sum_numb] = count[key_current];
			++sum_numb;
			cell[sum_numb] = count[key_x1];
			++sum_numb;
			cell[sum_numb] = count[key_xy1];
			++sum_numb;
			cell[sum_numb] = count[key_y1];
			++sum_numb;
			cell[sum_numb] = count[key_z1];
			++sum_numb;
			cell[sum_numb] = count[key_xz1];
			++sum_numb;
			cell[sum_numb] = count[key_xyz1];
			++sum_numb;
			cell[sum_numb] = count[key_yz1];
			++sum_numb;
		}
#endif
	}


	dataspace_id = H5Screate_simple(2, dims2, NULL);
	dataset_id = H5Dcreate(group_id, (ini_group + block_group + "/cell").c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);

	delete[] cell;
	H5Gclose(group_id);

}

/**
* @brief      function to write numercial boundary in HDF5 file.
* @param[in]  file_id      id of the open file.
*/
void HDF5H::write_numerical_boundary(const hid_t file_id)
{
	std::string ini_group = "/Ini_numb";
	hid_t group0_id = H5Gcreate(file_id, ini_group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	for (int level_index = 0; level_index < grid_vlevel.size(); ++level_index)
	{
		Grid * grid_ptr;
		if (grid_vlevel.at(level_index) == C_max_level)
		{
			grid_ptr = &(Grid_Manager::pointer_me->gr_inner);
		}
		else
		{
			grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(grid_vlevel.at(level_index)));
		}

		unsigned int ilevel = grid_vlevel.at(level_index);
		std::string block_group = "/block" + std::to_string(ilevel);
		hid_t group_id, dataspace_id, dataset_id;
		group_id = H5Gcreate(group0_id, (ini_group + block_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		// fine2coarse

		//for (int iboundary = 0; iboundary < C_overlap * 2 + 1; ++iboundary)
		for (int iboundary = 1; iboundary < 3; ++iboundary)
		{
			hid_t group1_id = H5Gcreate(group_id, (ini_group + block_group + "/boundary" + std::to_string(iboundary)).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			D_uint npoints = static_cast<D_uint> (grid_ptr->fine2coarse.at(iboundary).size());
			D_uint icount_points = 0, icount_cells = 0;

			D_mapint count;
			D_real *xyz = nullptr;
			xyz = new D_real[npoints*C_DIMS];

			hsize_t dims[2];
			dims[0] = npoints;
			dims[1] = C_DIMS;

			// write coodinates of nodes
			for (D_mapint::iterator iter = grid_ptr->fine2coarse.at(iboundary).begin(); iter != grid_ptr->fine2coarse.at(iboundary).end(); ++iter)
			{
				count.insert({ iter->first, icount_points });
				xyz[icount_points*C_DIMS] = 0;
				xyz[icount_points*C_DIMS + 1] = 0;
				D_morton key_current = iter->first;
#if (C_DIMS==2)
				Morton_Assist::pointer_me->compute_coordinate(iter->first, ilevel, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1]);
				xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
				xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
#endif
#if (C_DIMS==3)
				xyz[icount_points*C_DIMS + 2] = 0;
				Morton_Assist::pointer_me->compute_coordinate(iter->first, ilevel, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1], xyz[icount_points*C_DIMS + 2]);
				xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
				xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
				xyz[icount_points * C_DIMS + 2] -= Solid_Manager::pointer_me->shape_offest_z0_grid;
				bool bool_xy, bool_xz, bool_yz;
				D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
				D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
				D_morton key_z1 = Morton_Assist::find_z1(key_current, ilevel);
				D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);
				D_morton key_xz1 = Morton_Assist::find_z1(key_x1, ilevel);
				D_morton key_yz1 = Morton_Assist::find_z1(key_y1, ilevel);

				bool_xy = (grid_ptr->fine2coarse.at(iboundary).find(key_x1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_y1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_xy1) != grid_ptr->fine2coarse.at(iboundary).end());
				bool_xz = (grid_ptr->fine2coarse.at(iboundary).find(key_x1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_z1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_xz1) != grid_ptr->fine2coarse.at(iboundary).end());
				bool_yz = (grid_ptr->fine2coarse.at(iboundary).find(key_y1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_z1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_yz1) != grid_ptr->fine2coarse.at(iboundary).end());

				if (bool_xy)
				{
					++icount_cells;
				}
				if (bool_xz)
				{
					++icount_cells;
				}
				if (bool_yz)
				{
					++icount_cells;
				}
#endif
				++icount_points;
			}



			D_uint ncells = icount_cells;

			bouandry_fine2coarse_nelements.at(iboundary).push_back(ncells);
			bouandry_fine2coarse_npoints.at(iboundary).push_back(npoints);

			dataspace_id = H5Screate_simple(2, dims, NULL);
			dataset_id = H5Dcreate(group1_id, (ini_group + block_group + "/boundary" + std::to_string(iboundary) + "/xyz").c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);


			delete[] xyz;
#if (C_DIMS==3)
			int *cell = nullptr;

			hsize_t dims2[2];
			dims2[0] = ncells;
			dims2[1] = 4;
			cell = new int[static_cast<int> (dims2[0]) * static_cast<int> (dims2[1])];

			int sum_numb = 0;
			for (D_mapint::iterator iter = grid_ptr->fine2coarse.at(iboundary).begin(); iter != grid_ptr->fine2coarse.at(iboundary).end(); ++iter)
			{
				D_morton key_current = iter->first;
				bool bool_xy, bool_xz, bool_yz;
				D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
				D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
				D_morton key_z1 = Morton_Assist::find_z1(key_current, ilevel);
				D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);
				D_morton key_xz1 = Morton_Assist::find_z1(key_x1, ilevel);
				D_morton key_yz1 = Morton_Assist::find_z1(key_y1, ilevel);

				bool_xy = (grid_ptr->fine2coarse.at(iboundary).find(key_x1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_y1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_xy1) != grid_ptr->fine2coarse.at(iboundary).end());
				bool_xz = (grid_ptr->fine2coarse.at(iboundary).find(key_x1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_z1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_xz1) != grid_ptr->fine2coarse.at(iboundary).end());
				bool_yz = (grid_ptr->fine2coarse.at(iboundary).find(key_y1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_z1) != grid_ptr->fine2coarse.at(iboundary).end()) && (grid_ptr->fine2coarse.at(iboundary).find(key_yz1) != grid_ptr->fine2coarse.at(iboundary).end());

				if (bool_xy)
				{
					cell[sum_numb] = count[key_current];
					++sum_numb;
					cell[sum_numb] = count[key_x1];
					++sum_numb;
					cell[sum_numb] = count[key_xy1];
					++sum_numb;
					cell[sum_numb] = count[key_y1];
					++sum_numb;
				}
				if (bool_xz)
				{
					cell[sum_numb] = count[key_current];
					++sum_numb;
					cell[sum_numb] = count[key_x1];
					++sum_numb;
					cell[sum_numb] = count[key_xz1];
					++sum_numb;
					cell[sum_numb] = count[key_z1];
					++sum_numb;
				}
				if (bool_yz)
				{
					cell[sum_numb] = count[key_current];
					++sum_numb;
					cell[sum_numb] = count[key_y1];
					++sum_numb;
					cell[sum_numb] = count[key_yz1];
					++sum_numb;
					cell[sum_numb] = count[key_z1];
					++sum_numb;
				}
			}

			dataspace_id = H5Screate_simple(2, dims2, NULL);

			dataset_id = H5Dcreate(group1_id, (ini_group + block_group + "/boundary" + std::to_string(iboundary) + "/cell").c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);

			delete[] cell;
			H5Gclose(group1_id);
#endif
		}
		H5Gclose(group_id);
	}
	H5Gclose(group0_id);
}


/**
* @brief      function to write solid points using polyline in HDF5 format.
* @param[in]  file_id      id of the open file.
*/
void HDF5H::write_solid_polyline(const hid_t file_id)
{
	std::string ini_group = "/Ini_solid";
	hid_t group0_id = H5Gcreate(file_id, ini_group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	for (unsigned int iter = 0; iter < Solid_Manager::pointer_me->numb_solids; ++iter)
	{
		D_uint npoints = Solid_Manager::pointer_me->shape_solids.at(iter).numb_nodes;
		D_real *xyz = nullptr;
		xyz = new D_real[npoints*C_DIMS];

		std::string shaps_dir = "/solid" + std::to_string(iter);
		std::cout << "cencter of solid: " << iter << " is x = " << Solid_Manager::pointer_me->shape_solids.at(iter).x0 - Solid_Manager::pointer_me->shape_solids.at(iter).shape_offest_x0 << ", y = " << Solid_Manager::pointer_me->shape_solids.at(iter).y0 - Solid_Manager::pointer_me->shape_solids.at(iter).shape_offest_y0 << ", z = " << Solid_Manager::pointer_me->shape_solids.at(iter).z0 - Solid_Manager::pointer_me->shape_solids.at(iter).shape_offest_z0 << std::endl;
		for (D_uint i = 0; i < npoints; ++i)
		{
			xyz[i*C_DIMS] = Solid_Manager::pointer_me->shape_solids.at(iter).node.at(i).x  - Solid_Manager::pointer_me->shape_solids.at(iter).shape_offest_x0;
			xyz[i*C_DIMS + 1] = Solid_Manager::pointer_me->shape_solids.at(iter).node.at(i).y - Solid_Manager::pointer_me->shape_solids.at(iter).shape_offest_y0;
#if (C_DIMS == 3)
			xyz[i*C_DIMS + 2] = Solid_Manager::pointer_me->shape_solids.at(iter).node.at(i).z - Solid_Manager::pointer_me->shape_solids.at(iter).shape_offest_z0;
#endif
		}

		hid_t group_id = H5Gcreate(group0_id, (ini_group + shaps_dir).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hsize_t dim[1] = { npoints*C_DIMS };
		hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
		hid_t dataset_id = H5Dcreate(group_id, (ini_group + shaps_dir + "/xyz").c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		delete[] xyz;

		int *cell = nullptr;
		dim[0] = { npoints };
		cell = new int[npoints];
		
		for (unsigned int i = 0; i < npoints; ++i)
		{
			cell[i] = i;
		}

		dataspace_id = H5Screate_simple(1, dim, NULL);
		dataset_id = H5Dcreate(group_id, (ini_group + shaps_dir + "/cell").c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);

		delete[] cell;
		H5Gclose(group_id);
	}
	H5Gclose(group0_id);

}



/**
* @brief      function to write boundries of the compuational domain (i.e. background mesh) in HDF5 format.
* @param[in]  file_id      id of the open file.
*/
void HDF5H::write_domain_boundary(const hid_t file_id)
{
	std::string ini_group = "/Domain_boundary";
	hid_t group0_id = H5Gcreate(file_id, ini_group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	Grid * grid_ptr;
	grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(0));
#if (C_DIMS==3)
	std::array<D_mapint, 2> bk_x, bk_y, bk_z;
#endif
	// write coodinates of nodes (x boundary)
	for (int iboundary = 0; iboundary < 2; ++iboundary)
	{
#if (C_DIMS==3)
		for (unsigned int ilevel_boundary = 0; ilevel_boundary < C_max_level + 1; ++ilevel_boundary)
		{
			for (D_mapNodePtr::iterator iter = Grid_Manager::pointer_me->bk_boundary_x.at(ilevel_boundary).at(iboundary).begin(); iter != Grid_Manager::pointer_me->bk_boundary_x.at(ilevel_boundary).at(iboundary).end(); ++iter)
			{
				if (bk_x.at(iboundary).find(iter->first) == bk_x.at(iboundary).end())
				{
					bk_x.at(iboundary).insert({ iter->first , ilevel_boundary });
				}
				else if (bk_x.at(iboundary)[iter->first] > ilevel_boundary)
				{
					bk_x.at(iboundary)[iter->first] = ilevel_boundary;
				}
			}
			for (D_mapNodePtr::iterator iter = Grid_Manager::pointer_me->bk_boundary_y.at(ilevel_boundary).at(iboundary).begin(); iter != Grid_Manager::pointer_me->bk_boundary_y.at(ilevel_boundary).at(iboundary).end(); ++iter)
			{
				if (bk_y.at(iboundary).find(iter->first) == bk_y.at(iboundary).end())
				{
					bk_y.at(iboundary).insert({ iter->first , ilevel_boundary });
				}
				else if (bk_y.at(iboundary)[iter->first] > ilevel_boundary)
				{
					bk_y.at(iboundary)[iter->first] = ilevel_boundary;
				}
			}
			for (D_mapNodePtr::iterator iter = Grid_Manager::pointer_me->bk_boundary_z.at(ilevel_boundary).at(iboundary).begin(); iter != Grid_Manager::pointer_me->bk_boundary_z.at(ilevel_boundary).at(iboundary).end(); ++iter)
			{
				if (bk_z.at(iboundary).find(iter->first) == bk_z.at(iboundary).end())
				{
					bk_z.at(iboundary).insert({ iter->first , ilevel_boundary });
				}
				else if (bk_z.at(iboundary)[iter->first] > ilevel_boundary)
				{
					bk_z.at(iboundary)[iter->first] = ilevel_boundary;
				}
			}
		}
#endif

		hid_t group1_id = H5Gcreate(group0_id, (ini_group + "/boundary_x" + std::to_string(iboundary)).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		D_uint npoints = 0;
#if (C_DIMS==2)
		for (unsigned int ilevel_boundary = 0; ilevel_boundary < C_max_level + 1; ++ilevel_boundary)
		{
			npoints += static_cast<D_uint> (Grid_Manager::pointer_me->bk_boundary_x.at(ilevel_boundary).at(iboundary).size());
		}
#endif
#if (C_DIMS==3)
		npoints = bk_x.at(iboundary).size();
#endif

		D_uint icount_points = 0, icount_cells = 0;

		D_mapint count;
		D_real *xyz = nullptr;
		xyz = new D_real[npoints*C_DIMS];

		hsize_t dims[2];
		dims[0] = npoints;
		dims[1] = C_DIMS;


#if (C_DIMS==2)
		for (unsigned int ilevel_boundary = 0; ilevel_boundary < C_max_level + 1; ++ilevel_boundary)
		{
			for (D_mapNodePtr::iterator  iter = Grid_Manager::pointer_me->bk_boundary_x.at(ilevel_boundary).at(iboundary).begin(); iter != Grid_Manager::pointer_me->bk_boundary_x.at(ilevel_boundary).at(iboundary).end(); ++iter)
			{
				count.insert({ iter->first, icount_points });
				xyz[icount_points*C_DIMS] = 0;
				xyz[icount_points*C_DIMS + 1] = 0;
				D_morton key_current = iter->first;
				Morton_Assist::pointer_me->compute_coordinate(iter->first, ilevel_boundary, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1]);
				xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
				xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
				++icount_points;
			}
			
	    }
#endif
#if (C_DIMS==3)
		for (D_mapint::iterator iter = bk_x.at(iboundary).begin(); iter != bk_x.at(iboundary).end(); ++iter)
		{
			count.insert({ iter->first, icount_points });
			xyz[icount_points*C_DIMS] = 0;
			xyz[icount_points*C_DIMS + 1] = 0;
			xyz[icount_points*C_DIMS + 2] = 0;
			D_morton key_current = iter->first;
			Morton_Assist::pointer_me->compute_coordinate(iter->first, iter->second, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1], xyz[icount_points*C_DIMS + 2]);
			xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
			xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
			xyz[icount_points * C_DIMS + 2] -= Solid_Manager::pointer_me->shape_offest_z0_grid;

			bool bool_yz;
			D_morton key_x1 = Morton_Assist::find_x1(key_current, iter->second);
			D_morton key_y1 = Morton_Assist::find_y1(key_current, iter->second);
			D_morton key_z1 = Morton_Assist::find_z1(key_current, iter->second);
			//D_morton key_xy1 = Morton_Assist::find_y1(key_x1, iter->second);
			//D_morton key_xz1 = Morton_Assist::find_z1(key_x1, iter->second);
			D_morton key_yz1 = Morton_Assist::find_z1(key_y1, iter->second);

			//bool_xy = (bk_x.at(iboundary).find(key_x1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_y1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_xy1) != bk_x.at(iboundary).end());
			//bool_xz = (bk_x.at(iboundary).find(key_x1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_z1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_xz1) != bk_x.at(iboundary).end());
			bool_yz = (bk_x.at(iboundary).find(key_y1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_z1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_yz1) != bk_x.at(iboundary).end());

			//if (bool_xy)
			//{
			//	++icount_cells;
			//}
			//if (bool_xz)
			//{
			//	++icount_cells;
			//}
			if (bool_yz)
			{
				++icount_cells;
			}


			++icount_points;
		}
		
#endif

		D_uint ncells = icount_cells;

		domain_boundary_nelements.at(0).at(iboundary) = ncells;
		domain_boundary_npoints.at(0).at(iboundary) = npoints;

		hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
		hid_t dataset_id = H5Dcreate(group1_id, (ini_group + "/boundary_x" + std::to_string(iboundary) + "/xyz").c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		delete[] xyz;

#if (C_DIMS==3)
		int *cell = nullptr;

		hsize_t dims2[2];
		dims2[0] = ncells;
		dims2[1] = 4;
		cell = new int[static_cast<int> (dims2[0]) * static_cast<int> (dims2[1])];
		int sum_numb = 0;
		for (D_mapint::iterator iter = bk_x.at(iboundary).begin(); iter != bk_x.at(iboundary).end(); ++iter)
		{
			D_morton key_current = iter->first;
			bool bool_yz;
			D_morton key_x1 = Morton_Assist::find_x1(key_current, iter->second);
			D_morton key_y1 = Morton_Assist::find_y1(key_current, iter->second);
			D_morton key_z1 = Morton_Assist::find_z1(key_current, iter->second);
			//D_morton key_xy1 = Morton_Assist::find_y1(key_x1, iter->second);
			//D_morton key_xz1 = Morton_Assist::find_z1(key_x1, iter->second);
			D_morton key_yz1 = Morton_Assist::find_z1(key_y1, iter->second);

			//bool_xy = (bk_x.at(iboundary).find(key_x1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_y1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_xy1) != bk_x.at(iboundary).end());
			//bool_xz = (bk_x.at(iboundary).find(key_x1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_z1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_xz1) != bk_x.at(iboundary).end());
			bool_yz = (bk_x.at(iboundary).find(key_y1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_z1) != bk_x.at(iboundary).end()) && (bk_x.at(iboundary).find(key_yz1) != bk_x.at(iboundary).end());

			//if (bool_xy)
			//{
			//	cell[sum_numb] = count[key_current];
			//	++sum_numb;
			//	cell[sum_numb] = count[key_x1];
			//	++sum_numb;
			//	cell[sum_numb] = count[key_xy1];
			//	++sum_numb;
			//	cell[sum_numb] = count[key_y1];
			//	++sum_numb;
			//}
			//if (bool_xz)
			//{
			//	cell[sum_numb] = count[key_current];
			//	++sum_numb;
			//	cell[sum_numb] = count[key_x1];
			//	++sum_numb;
			//	cell[sum_numb] = count[key_xy1];
			//	++sum_numb;
			//	cell[sum_numb] = count[key_y1];
			//	++sum_numb;
			//}
			if (bool_yz)
			{
				cell[sum_numb] = count[key_current];
				++sum_numb;
				cell[sum_numb] = count[key_y1];
				++sum_numb;
				cell[sum_numb] = count[key_yz1];
				++sum_numb;
				cell[sum_numb] = count[key_z1];
				++sum_numb;
			}
		}

		dataspace_id = H5Screate_simple(2, dims2, NULL);

		dataset_id = H5Dcreate(group1_id, (ini_group + "/boundary_x" + std::to_string(iboundary) + "/cell").c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);

		delete[] cell;
		H5Gclose(group1_id);
#endif

	}

	// write coodinates of nodes (y boundary)

	for (int iboundary = 0; iboundary < 2; ++iboundary)
	{
		hid_t group1_id = H5Gcreate(group0_id, (ini_group + "/boundary_y" + std::to_string(iboundary)).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		D_uint npoints = 0;
#if (C_DIMS==2)
		for (unsigned int ilevel_boundary = 0; ilevel_boundary < C_max_level + 1; ++ilevel_boundary)
		{
			npoints += static_cast<D_uint> (Grid_Manager::pointer_me->bk_boundary_y.at(ilevel_boundary).at(iboundary).size());
		}
#endif
#if (C_DIMS==3)
		npoints = bk_y.at(iboundary).size();
#endif		
		D_uint icount_points = 0, icount_cells = 0;

		D_mapint count;
		D_real *xyz = nullptr;
		xyz = new D_real[npoints*C_DIMS];

		hsize_t dims[2];
		dims[0] = npoints;
		dims[1] = C_DIMS;


#if (C_DIMS==2)
		for (unsigned int ilevel_boundary = 0; ilevel_boundary < C_max_level + 1; ++ilevel_boundary)
		{
			for (D_mapNodePtr::iterator iter = Grid_Manager::pointer_me->bk_boundary_y.at(ilevel_boundary).at(iboundary).begin(); iter != Grid_Manager::pointer_me->bk_boundary_y.at(ilevel_boundary).at(iboundary).end(); ++iter)
			{
				count.insert({ iter->first, icount_points });
				xyz[icount_points*C_DIMS] = 0;
				xyz[icount_points*C_DIMS + 1] = 0;
				D_morton key_current = iter->first;
				Morton_Assist::pointer_me->compute_coordinate(iter->first, ilevel_boundary, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1]);
				xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
				xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
				++icount_points;
			}		
		}
#endif
#if (C_DIMS==3)
		for (D_mapint::iterator iter = bk_y.at(iboundary).begin(); iter != bk_y.at(iboundary).end(); ++iter)
		{
			count.insert({ iter->first, icount_points });
			xyz[icount_points*C_DIMS] = 0;
			xyz[icount_points*C_DIMS + 1] = 0;
			xyz[icount_points*C_DIMS + 2] = 0;
			D_morton key_current = iter->first;
			Morton_Assist::pointer_me->compute_coordinate(iter->first, iter->second, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1], xyz[icount_points*C_DIMS + 2]);
			xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
			xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
			xyz[icount_points * C_DIMS + 2] -= Solid_Manager::pointer_me->shape_offest_z0_grid;

			bool bool_xz;
			D_morton key_x1 = Morton_Assist::find_x1(key_current, iter->second);
			D_morton key_y1 = Morton_Assist::find_y1(key_current, iter->second);
			D_morton key_z1 = Morton_Assist::find_z1(key_current, iter->second);
			D_morton key_xz1 = Morton_Assist::find_z1(key_x1, iter->second);

			bool_xz = (bk_y.at(iboundary).find(key_x1) != bk_y.at(iboundary).end()) && (bk_y.at(iboundary).find(key_z1) != bk_y.at(iboundary).end()) && (bk_y.at(iboundary).find(key_xz1) != bk_y.at(iboundary).end());

			if (bool_xz)
			{
				++icount_cells;
			}

			++icount_points;
		}
#endif		

		D_uint ncells = icount_cells;

		domain_boundary_nelements.at(1).at(iboundary) = ncells;
		domain_boundary_npoints.at(1).at(iboundary) = npoints;

		hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
		hid_t dataset_id = H5Dcreate(group1_id, (ini_group + "/boundary_y" + std::to_string(iboundary) + "/xyz").c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		delete[] xyz;

#if (C_DIMS==3)
		int *cell = nullptr;

		hsize_t dims2[2];
		dims2[0] = ncells;
		dims2[1] = 4;
		cell = new int[static_cast<int> (dims2[0]) * static_cast<int> (dims2[1])];

		int sum_numb = 0;
		for (D_mapint::iterator iter = bk_y.at(iboundary).begin(); iter != bk_y.at(iboundary).end(); ++iter)
		{
			D_morton key_current = iter->first;
			bool bool_xz;
			D_morton key_x1 = Morton_Assist::find_x1(key_current, iter->second);
			D_morton key_y1 = Morton_Assist::find_y1(key_current, iter->second);
			D_morton key_z1 = Morton_Assist::find_z1(key_current, iter->second);
			D_morton key_xz1 = Morton_Assist::find_z1(key_x1, iter->second);

			bool_xz = (bk_y.at(iboundary).find(key_x1) != bk_y.at(iboundary).end()) && (bk_y.at(iboundary).find(key_z1) != bk_y.at(iboundary).end()) && (bk_y.at(iboundary).find(key_xz1) != bk_y.at(iboundary).end());

			if (bool_xz)
			{
				cell[sum_numb] = count[key_current];
				++sum_numb;
				cell[sum_numb] = count[key_x1];
				++sum_numb;
				cell[sum_numb] = count[key_xz1];
				++sum_numb;
				cell[sum_numb] = count[key_z1];
				++sum_numb;
			}
		}

		dataspace_id = H5Screate_simple(2, dims2, NULL);

		dataset_id = H5Dcreate(group1_id, (ini_group + "/boundary_y" + std::to_string(iboundary) + "/cell").c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);

		delete[] cell;
		H5Gclose(group1_id);
#endif

	}
#if (C_DIMS==3)
	// write coodinates of nodes (z boundary)
	for (int iboundary = 0; iboundary < 2; ++iboundary)
	{
		hid_t group1_id = H5Gcreate(group0_id, (ini_group + "/boundary_z" + std::to_string(iboundary)).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		D_uint npoints = bk_z.at(iboundary).size();

		D_uint icount_points = 0, icount_cells = 0;

		D_mapint count;
		D_real *xyz = nullptr;
		xyz = new D_real[npoints*C_DIMS];

		hsize_t dims[2];
		dims[0] = npoints;
		dims[1] = C_DIMS;
		
		for (D_mapint::iterator iter = bk_z.at(iboundary).begin(); iter != bk_z.at(iboundary).end(); ++iter)
		{
			count.insert({ iter->first, icount_points });
			xyz[icount_points*C_DIMS] = 0;
			xyz[icount_points*C_DIMS + 1] = 0;
			xyz[icount_points*C_DIMS + 2] = 0;
			D_morton key_current = iter->first;
			Morton_Assist::pointer_me->compute_coordinate(iter->first, iter->second, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1], xyz[icount_points*C_DIMS + 2]);
			xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
			xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
			xyz[icount_points * C_DIMS + 2] -= Solid_Manager::pointer_me->shape_offest_z0_grid;

			bool bool_xy;
			D_morton key_x1 = Morton_Assist::find_x1(key_current, iter->second);
			D_morton key_y1 = Morton_Assist::find_y1(key_current, iter->second);
			D_morton key_z1 = Morton_Assist::find_z1(key_current, iter->second);
			D_morton key_xy1 = Morton_Assist::find_y1(key_x1, iter->second);

			bool_xy = (bk_z.at(iboundary).find(key_x1) != bk_z.at(iboundary).end()) && (bk_z.at(iboundary).find(key_y1) != bk_z.at(iboundary).end()) && (bk_z.at(iboundary).find(key_xy1) != bk_z.at(iboundary).end());

			if (bool_xy)
			{
				++icount_cells;
			}
			
			++icount_points;
		}
		

		D_uint ncells = icount_cells;

		domain_boundary_nelements.at(2).at(iboundary) = ncells;
		domain_boundary_npoints.at(2).at(iboundary) = npoints;

		hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
		hid_t dataset_id = H5Dcreate(group1_id, (ini_group + "/boundary_z" + std::to_string(iboundary) + "/xyz").c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		delete[] xyz;

		int *cell = nullptr;

		hsize_t dims2[2];
		dims2[0] = ncells;
		dims2[1] = 4;
		cell = new int[static_cast<int> (dims2[0]) * static_cast<int> (dims2[1])];

		int sum_numb = 0;
		for (D_mapint::iterator iter = bk_z.at(iboundary).begin(); iter != bk_z.at(iboundary).end(); ++iter)
		{
			D_morton key_current = iter->first;
			bool bool_xy;
			D_morton key_x1 = Morton_Assist::find_x1(key_current, iter->second);
			D_morton key_y1 = Morton_Assist::find_y1(key_current, iter->second);
			D_morton key_z1 = Morton_Assist::find_z1(key_current, iter->second);
			D_morton key_xy1 = Morton_Assist::find_y1(key_x1, iter->second);

			bool_xy = (bk_z.at(iboundary).find(key_x1) != bk_z.at(iboundary).end()) && (bk_z.at(iboundary).find(key_y1) != bk_z.at(iboundary).end()) && (bk_z.at(iboundary).find(key_xy1) != bk_z.at(iboundary).end());

			
			if (bool_xy)
			{
				cell[sum_numb] = count[key_current];
				++sum_numb;
				cell[sum_numb] = count[key_x1];
				++sum_numb;
				cell[sum_numb] = count[key_xy1];
				++sum_numb;
				cell[sum_numb] = count[key_y1];
				++sum_numb;
			}
		}

		dataspace_id = H5Screate_simple(2, dims2, NULL);

		dataset_id = H5Dcreate(group1_id, (ini_group + "/boundary_z" + std::to_string(iboundary) + "/cell").c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);

		delete[] cell;
		H5Gclose(group1_id);
	}
#endif

	H5Gclose(group0_id);
}

#if	(C_SOLID_BOUNDARY == 2)
#if (C_DEBUG == 1)
void HDF5H::write_output_points_update(const hid_t file_id)
{
	std::string ini_group = "/Points_update";
	hid_t group0_id = H5Gcreate(file_id, ini_group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	unsigned int icount_points = 0;

	D_uint npoints = Grid_Manager::pointer_me->output_points_update.size();
	D_real *xyz = nullptr;
	xyz = new D_real[npoints*C_DIMS];
	for (D_mapint::iterator iter = Grid_Manager::pointer_me->output_points_update.begin(); iter != Grid_Manager::pointer_me->output_points_update.end(); ++iter)
	{
#if(C_DIMS==2)
		Morton_Assist::pointer_me->compute_coordinate(iter->first, iter->second, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1]);
		xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
		xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
#endif
#if(C_DIMS==3)
		Morton_Assist::pointer_me->compute_coordinate(iter->first, iter->second, xyz[icount_points*C_DIMS], xyz[icount_points*C_DIMS + 1], xyz[icount_points*C_DIMS + 2]);
		xyz[icount_points * C_DIMS] -= Solid_Manager::pointer_me->shape_offest_x0_grid;
		xyz[icount_points * C_DIMS + 1] -= Solid_Manager::pointer_me->shape_offest_y0_grid;
		xyz[icount_points * C_DIMS + 2] -= Solid_Manager::pointer_me->shape_offest_z0_grid;
#endif
		++icount_points;
	}


	hsize_t dim[1] = { npoints*C_DIMS };
	hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
	hid_t dataset_id = H5Dcreate(group0_id, (ini_group + "/xyz").c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	delete[] xyz;

	H5Gclose(group0_id);
}
#endif
#endif

/**
* @brief      function to write infomation of fluid nodes in XDMF file.
* @param[in]  outfile      name for output file.
*/
void HDF5H::write_xdmf(std::string outfile)
{
	std::string outname;
	outname = outfile + ".xmf";
	std::string ini_group = "/Ini_grid";
	std::ofstream out;
	out.open(outname);
	std::string space2="  ";
	std::string space4 = space2 + space2;
	std::string space6 = space4 + space2;

	// Header
	out << "<?xml version=\"1.0\" ?>" << std::endl;
	out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
	out << "<Xdmf Version = \"2.0\">" << std::endl;
	out << "<Domain>" << std::endl;
	// Write fluid
	// block tree
	out << "<Grid Name=\"Ini_Blcok\" GridType=\"Tree\">" << std::endl;
	// different blocks
	for (int level_index = 0; level_index < grid_vlevel.size(); ++level_index)
	{
		unsigned int ilevel = grid_vlevel.at(level_index);
		std::string block_level = "/block" + std::to_string(ilevel);

		out << space2 << "<Grid Name=\"Refine " << grid_vlevel.at(level_index) <<"\" GridType=\"Uniform\">" << std::endl;

		//out << space4 << "<Topology TopologyType=\"Polyvertex\"" << std::endl;
		//out << space4 << "NodesPerElement = \"" << grid_npoints.at(level_index) << "\">" << std::endl;
#if (C_DIMS==2)
		out << space4 << "<Topology TopologyType=\"Quadrilateral\" NumberOfElements=\"" << grid_nelements.at(level_index) << "\">" << std::endl;
#endif

#if (C_DIMS==3)
		out << space4 << "<Topology TopologyType=\"Hexahedron\" NumberOfElements=\"" << grid_nelements.at(level_index) << "\">" << std::endl;
#endif
		out << space4 << "<DataItem Name=\"cell\"" << std::endl;
		out << space6 << "DataType=\"Int\"" << std::endl;
		out << space6 << "Dimensions=\"" << grid_nelements.at(level_index) << " " << two_power_n(C_DIMS) << "\" " << std::endl;
		out << space6 << "Format=\"HDF\">" << std::endl;
		out << space6 << outfile << ".h5:"<< ini_group << block_level <<"/cell" << std::endl;
		out << space4 << "</DataItem>" << std::endl;
		out << space4 << "</Topology>" << std::endl;
		
#if (C_DIMS==2)
		out << space4 << "<Geometry GeometryType=\"XY\">" << std::endl;
#endif
#if (C_DIMS==3)
		out << space4 << "<Geometry GeometryType=\"XYZ\">" << std::endl;
#endif

		out << space4 << "<DataItem Name=\"xyz\"" << std::endl;
		out << space6 << "Dimensions=\"" << grid_npoints.at(level_index) * C_DIMS << "\"" << std::endl;
		out << space6 << "DataType=\"Float\" Precision=\"8\"" << std::endl;
		out << space6 << "Format=\"HDF\">" << std::endl;
		out << space6 << outfile << ".h5:" << ini_group << block_level << "/xyz" << std::endl;
		out << space4 << "</DataItem>" << std::endl;

		out << space4 << "</Geometry>" << std::endl;

		out << space4 << "<Attribute Center = \"Node\" Name = \"Flag\" Type = \"Scalar\">" << std::endl;
		out << space4 << "<DataItem " << std::endl;
		out << space6 << "Dimensions=\"" << grid_npoints.at(level_index) << "\"" << std::endl;
		out << space6 << "DataType=\"Int\"" << std::endl;
		out << space6 << "Format=\"HDF\">" << std::endl;
		out << space6 << outfile << ".h5:" << "/Flag" << block_level << "/flag" << std::endl;
		out << space4 << "</DataItem>" << std::endl;
		out << space4 << "</Attribute>" << std::endl;

		out << space2 << "</Grid>" << std::endl;

	}
	out << "</Grid>" << std::endl;
	if (Solid_Manager::pointer_me->numb_solids > 0)
	{
		write_xdmf_solid_ployline(out, outfile);
	}

	if (flag_write_numerical_boundary)
	{
		write_xdmf_numerical_boundary(out, outfile);
	}

	if (flag_write_domain_boundary)
	{
		write_xdmf_domain_boundary(out, outfile);
	}

#if	(C_SOLID_BOUNDARY == 2)
#if (C_DEBUG == 1)
	write_xdm_output_points_update(out, outfile);
#endif
#endif

	out << "</Domain>" << std::endl;
	out << "</Xdmf>" << std::endl;
}

/**
* @brief      function to write solid using ployline in XDMF file.
* @param[in]  out      stream writing to XDMF file.
* @param[in]  outfile      name for output file.
*/
void HDF5H::write_xdmf_numerical_boundary(std::ofstream &out, const std::string &outfile)
{
	std::string ini_group = "/Ini_numb";
	std::string space2 = "  ";
	std::string space4 = space2 + space2;
	std::string space6 = space4 + space2;
#if (C_DIMS==2)
	int dims = 1;
#endif
#if (C_DIMS==3)
	int dims = 4;
#endif

	// block tree
	out << "<Grid Name=\"Ini_numb\" GridType=\"Tree\">" << std::endl;
	// different blocks
	for (int level_index = 0; level_index < grid_vlevel.size(); ++level_index)
	{
		unsigned int ilevel = grid_vlevel.at(level_index);
		//if (ilevel > 0)
		{
			std::string block_group = "/block" + std::to_string(ilevel);
			out << "<Grid Name=\"Refine_numb " << grid_vlevel.at(level_index) << "\" GridType=\"Tree\">" << std::endl;
			for (int iboundary = 1; iboundary < 3; ++iboundary)
			{
				out << space2 << "<Grid Name=\"lvl_" << ilevel << "F2C_ " << iboundary << "\" GridType=\"Uniform\">" << std::endl;

#if (C_DIMS==2)
				out << space4 << "<Topology TopologyType=\"Polyvertex\"" << std::endl;
				out << space4 << "NodesPerElement = \"" << bouandry_fine2coarse_npoints.at(iboundary).at(level_index) << "\">" << std::endl;
#endif
#if (C_DIMS==3)
				//out << space4 << "<Topology TopologyType=\"Polyvertex\"" << std::endl;
				//out << space4 << "NodesPerElement = \"" << bouandry_fine2coarse_npoints.at(iboundary).at(level_index) << "\">" << std::endl;

				out << space4 << "<Topology TopologyType=\"Quadrilateral\" NumberOfElements=\"" << bouandry_fine2coarse_nelements.at(iboundary).at(level_index) << "\">" << std::endl;
				out << space4 << "<DataItem Name=\"cell\"" << std::endl;
				out << space6 << "DataType=\"Int\"" << std::endl;
				out << space6 << "Dimensions=\"" << bouandry_fine2coarse_nelements.at(iboundary).at(level_index) << " " << dims << "\" " << std::endl;
				out << space6 << "Format=\"HDF\">" << std::endl;
				out << space6 << outfile << ".h5:" << ini_group << block_group << "/boundary" << iboundary << "/cell" << std::endl;
				out << space4 << "</DataItem>" << std::endl;
#endif
				out << space4 << "</Topology>" << std::endl;

#if (C_DIMS==2)
				out << space4 << "<Geometry GeometryType=\"XY\">" << std::endl;
#endif
#if (C_DIMS==3)
				out << space4 << "<Geometry GeometryType=\"XYZ\">" << std::endl;
#endif

				out << space4 << "<DataItem Name=\"xyz\"" << std::endl;
				out << space6 << "Dimensions=\"" << bouandry_fine2coarse_npoints.at(iboundary).at(level_index) * C_DIMS << "\"" << std::endl;
				out << space6 << "DataType=\"Float\" Precision=\"8\"" << std::endl;
				out << space6 << "Format=\"HDF\">" << std::endl;
				out << space6 << outfile << ".h5:" << ini_group << block_group << "/boundary" << iboundary << "/xyz" << std::endl;
				out << space4 << "</DataItem>" << std::endl;

				out << space4 << "</Geometry>" << std::endl;

				out << space2 << "</Grid>" << std::endl;
			}
			out << "</Grid>" << std::endl;
		}
	}
	out << "</Grid>" << std::endl;
}

/**
* @brief      function to write solid using ployline in XDMF file.
* @param[in]  out      stream writing to XDMF file.
* @param[in]  outfile      name for output file.
* @note       for both 2D and 3D cases, the nodes are connected sequentially
*/
void HDF5H::write_xdmf_solid_ployline(std::ofstream &out, const std::string &outfile)
{
	std::string ini_group = "/Ini_solid";
	std::string space2 = "  ";
	std::string space4 = space2 + space2;
	std::string space6 = space4 + space2;
	out << "<Grid Name=\"Ini_Solid\" GridType=\"Tree\">" << std::endl;
	for(unsigned int index_solid = 0; index_solid < Solid_Manager::pointer_me->numb_solids; ++index_solid)
	{
		D_uint npoints = Solid_Manager::pointer_me->shape_solids.at(index_solid).numb_nodes;
		std::string shaps_dir = "/solid" + std::to_string(index_solid);
		
		out << space2 << "<Grid Name=\"Solid " << index_solid << "\" GridType=\"Uniform\">" << std::endl;
		out << space4 << "<Topology TopologyType=\"Polyline\">" << std::endl;
		out << space4 << "<DataItem Name=\"connection\"" << std::endl;
		out << space6 << "DataType=\"Int\"" << std::endl;
		out << space6 << "Dimensions=\"" << 1 << " " << npoints << "\" " << std::endl;
		out << space6 << "Format=\"HDF\">" << std::endl;
	    out << space6 << outfile << ".h5:" << ini_group << shaps_dir << "/cell" << std::endl;
		out << space4 << "</DataItem>" << std::endl;
		out << space4 << "</Topology>" << std::endl;

#if (C_DIMS==2)
		out << space4 << "<Geometry GeometryType=\"XY\">" << std::endl;
#endif
#if (C_DIMS==3)
		out << space4 << "<Geometry GeometryType=\"XYZ\">" << std::endl;
#endif

		out << space4 << "<DataItem Name=\"xyz\"" << std::endl;
		out << space6 << "Dimensions=\"" << npoints * C_DIMS << "\"" << std::endl;
		out << space6 << "DataType=\"Float\" Precision=\"8\"" << std::endl;
		out << space6 << "Format=\"HDF\">" << std::endl;
		out << space6 << outfile << ".h5:" << ini_group << shaps_dir << "/xyz" << std::endl;
		out << space4 << "</DataItem>" << std::endl;

		out << space4 << "</Geometry>" << std::endl;

		out << space2 << "</Grid>" << std::endl;

	}
	out << "</Grid>" << std::endl;
}


/**
* @brief      function to write boundries of the compuational domain (i.e. background mesh) in XDMF file.
* @param[in]  out      stream writing to XDMF file.
* @param[in]  outfile      name for output file.
* @note       for 2D cases, the boundary is represented by separated vertices; for 3D cases, the boundary is represented by quadrilateral.
*/
void HDF5H::write_xdmf_domain_boundary(std::ofstream &out, const std::string &outfile)
{
	std::string ini_group = "/Domain_boundary";
	std::string space2 = "  ";
	std::string space4 = space2 + space2;
	std::string space6 = space4 + space2;
#if (C_DIMS==2)
	int dims = 1;
#endif
#if (C_DIMS==3)
	int dims = 4;
#endif
	// block tree
	out << "<Grid Name=\"Domain_boundary\" GridType=\"Tree\">" << std::endl;
	// boundary x
		for (int iboundary = 0; iboundary < 2; ++iboundary)
		{
			out << space2 << "<Grid Name=\"Domain_boundary_x " << iboundary << "\" GridType=\"Uniform\">" << std::endl;

#if (C_DIMS==2)
			out << space4 << "<Topology TopologyType=\"Polyvertex\"" << std::endl;
			out << space4 << "NodesPerElement = \"" << domain_boundary_npoints.at(0).at(iboundary) << "\">" << std::endl;
#endif
#if (C_DIMS==3)
			//out << space4 << "<Topology TopologyType=\"Polyvertex\"" << std::endl;
			//out << space4 << "NodesPerElement = \"" << domain_boundary_npoints.at(0).at(iboundary) << "\">" << std::endl;

			out << space4 << "<Topology TopologyType=\"Quadrilateral\" NumberOfElements=\"" << domain_boundary_nelements.at(0).at(iboundary) << "\">" << std::endl;
			out << space4 << "<DataItem Name=\"cell\"" << std::endl;
			out << space6 << "DataType=\"Int\"" << std::endl;
			out << space6 << "Dimensions=\"" << domain_boundary_nelements.at(0).at(iboundary) << " " << dims << "\" " << std::endl;
			out << space6 << "Format=\"HDF\">" << std::endl;
			out << space6 << outfile << ".h5:" << ini_group << "/boundary_x" << iboundary << "/cell" << std::endl;
			out << space4 << "</DataItem>" << std::endl;
#endif
			out << space4 << "</Topology>" << std::endl;

#if (C_DIMS==2)
			out << space4 << "<Geometry GeometryType=\"XY\">" << std::endl;
#endif
#if (C_DIMS==3)
			out << space4 << "<Geometry GeometryType=\"XYZ\">" << std::endl;
#endif

			out << space4 << "<DataItem Name=\"xyz\"" << std::endl;
			out << space6 << "Dimensions=\"" << domain_boundary_npoints.at(0).at(iboundary)  * C_DIMS << "\"" << std::endl;
			out << space6 << "DataType=\"Float\" Precision=\"8\"" << std::endl;
			out << space6 << "Format=\"HDF\">" << std::endl;
			out << space6 << outfile << ".h5:" << ini_group << "/boundary_x" << iboundary << "/xyz" << std::endl;
			out << space4 << "</DataItem>" << std::endl;

			out << space4 << "</Geometry>" << std::endl;

			out << space2 << "</Grid>" << std::endl;
		}

		// boundary y
		for (int iboundary = 0; iboundary < 2; ++iboundary)
		{
			out << space2 << "<Grid Name=\"Domain_boundary_y " << iboundary << "\" GridType=\"Uniform\">" << std::endl;

#if (C_DIMS==2)
			out << space4 << "<Topology TopologyType=\"Polyvertex\"" << std::endl;
			out << space4 << "NodesPerElement = \"" << domain_boundary_npoints.at(1).at(iboundary) << "\">" << std::endl;
#endif
#if (C_DIMS==3)
			out << space4 << "<Topology TopologyType=\"Quadrilateral\" NumberOfElements=\"" << domain_boundary_nelements.at(1).at(iboundary) << "\">" << std::endl;
			out << space4 << "<DataItem Name=\"cell\"" << std::endl;
			out << space6 << "DataType=\"Int\"" << std::endl;
			out << space6 << "Dimensions=\"" << domain_boundary_nelements.at(1).at(iboundary) << " " << dims << "\" " << std::endl;
			out << space6 << "Format=\"HDF\">" << std::endl;
			out << space6 << outfile << ".h5:" << ini_group << "/boundary_y" << iboundary << "/cell" << std::endl;
			out << space4 << "</DataItem>" << std::endl;
#endif
			out << space4 << "</Topology>" << std::endl;

#if (C_DIMS==2)
			out << space4 << "<Geometry GeometryType=\"XY\">" << std::endl;
#endif
#if (C_DIMS==3)
			out << space4 << "<Geometry GeometryType=\"XYZ\">" << std::endl;
#endif

			out << space4 << "<DataItem Name=\"xyz\"" << std::endl;
			out << space6 << "Dimensions=\"" << domain_boundary_npoints.at(1).at(iboundary)  * C_DIMS << "\"" << std::endl;
			out << space6 << "DataType=\"Float\" Precision=\"8\"" << std::endl;
			out << space6 << "Format=\"HDF\">" << std::endl;
			out << space6 << outfile << ".h5:" << ini_group << "/boundary_y" << iboundary << "/xyz" << std::endl;
			out << space4 << "</DataItem>" << std::endl;

			out << space4 << "</Geometry>" << std::endl;

			out << space2 << "</Grid>" << std::endl;
		}

#if (C_DIMS==3)
		// boundary z
		for (int iboundary = 0; iboundary < 2; ++iboundary)
		{
			out << space2 << "<Grid Name=\"Domain_boundary_z " << iboundary << "\" GridType=\"Uniform\">" << std::endl;

			//out << space4 << "<Topology TopologyType=\"Polyvertex\"" << std::endl;
           //out << space4 << "NodesPerElement = \"" << domain_boundary_npoints.at(2).at(iboundary) << "\">" << std::endl;

			out << space4 << "<Topology TopologyType=\"Quadrilateral\" NumberOfElements=\"" << domain_boundary_nelements.at(2).at(iboundary) << "\">" << std::endl;
			out << space4 << "<DataItem Name=\"cell\"" << std::endl;
			out << space6 << "DataType=\"Int\"" << std::endl;
			out << space6 << "Dimensions=\"" << domain_boundary_nelements.at(2).at(iboundary) << " " << dims << "\" " << std::endl;
			out << space6 << "Format=\"HDF\">" << std::endl;
			out << space6 << outfile << ".h5:" << ini_group<< "/boundary_z" << iboundary << "/cell" << std::endl;
			out << space4 << "</DataItem>" << std::endl;
			out << space4 << "</Topology>" << std::endl;

			out << space4 << "<Geometry GeometryType=\"XYZ\">" << std::endl;
			out << space4 << "<DataItem Name=\"xyz\"" << std::endl;
			out << space6 << "Dimensions=\"" << domain_boundary_npoints.at(2).at(iboundary)  * C_DIMS << "\"" << std::endl;
			out << space6 << "DataType=\"Float\" Precision=\"8\"" << std::endl;
			out << space6 << "Format=\"HDF\">" << std::endl;
			out << space6 << outfile << ".h5:" << ini_group << "/boundary_z" << iboundary << "/xyz" << std::endl;
			out << space4 << "</DataItem>" << std::endl;
			out << space4 << "</Geometry>" << std::endl;

			out << space2 << "</Grid>" << std::endl;
		}
#endif
		out << "</Grid>" << std::endl;
}

#if	(C_SOLID_BOUNDARY == 2)
#if(C_DEBUG==1)
void HDF5H::write_xdm_output_points_update(std::ofstream &out, const std::string &outfile)
{
	std::string ini_group = "/Points_update";
	std::string space2 = "  ";
	std::string space4 = space2 + space2;
	std::string space6 = space4 + space2;
	out << "<Grid Name=\"Points_update\" GridType=\"Tree\">" << std::endl;

		D_uint npoints = Grid_Manager::pointer_me->output_points_update.size();
		std::string shaps_dir = "/Points_update";

		out << space2 << "<Grid Name=\"Points_update "<< "\" GridType=\"Uniform\">" << std::endl;
		out << space4 << "<Topology TopologyType=\"Polyvertex\"" << std::endl;
		out << space4 << "NodesPerElement = \"" << npoints << "\">" << std::endl;

		out << space4 << "</Topology>" << std::endl;

#if (C_DIMS==2)
		out << space4 << "<Geometry GeometryType=\"XY\">" << std::endl;
#endif
#if (C_DIMS==3)
		out << space4 << "<Geometry GeometryType=\"XYZ\">" << std::endl;
#endif

		out << space4 << "<DataItem Name=\"xyz\"" << std::endl;
		out << space6 << "Dimensions=\"" << npoints * C_DIMS << "\"" << std::endl;
		out << space6 << "DataType=\"Float\" Precision=\"8\"" << std::endl;
		out << space6 << "Format=\"HDF\">" << std::endl;
		out << space6 << outfile << ".h5:" << ini_group << "/xyz" << std::endl;
		out << space4 << "</DataItem>" << std::endl;

		out << space4 << "</Geometry>" << std::endl;

		out << space2 << "</Grid>" << std::endl;

	
	out << "</Grid>" << std::endl;
}
#endif
#endif