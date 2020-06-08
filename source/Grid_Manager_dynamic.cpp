/**
* @file
* @author Zhengliang Liu
* @brief Mesh update.
* @note .
*/
#include "General.h"
#include "Grid_Manager.h"
#include "Solid_Manager.h"
#include "Morton_assist.h"
int iter_total = 100; // total number of iterations in while loop to avoid infinite itering
#define C_CHECK_BOUDNARY_METHOD 1  // if is 1, erase node when removing nodes, if is 2, erase node at the end of  constructing boundary for removing nodes

#if	(C_SOLID_BOUNDARY == 2)
/**
* @brief      function to initilize reference coordinates (e.g. xref, xref2) of the solid geomery.
* @param[in]  ishape    order of shape.
* @param[in]  ipoint    order of the solid node.
*/
void Grid_Manager::initial_ref_coordinate(unsigned int ishape, D_uint ipoint)
{
	Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).xref = floor(Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).x / gr_inner.dx + C_eps) * gr_inner.dx;
	Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).yref = floor(Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).y / gr_inner.dx + C_eps) * gr_inner.dx;
	Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).xref2 = floor(Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).x / gr_inner.dx / 2. + C_eps) * gr_inner.dx * 2;
	Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).yref2 = floor(Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).y / gr_inner.dx / 2. + C_eps) * gr_inner.dx * 2;
#if(C_DIMS==3)
	Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).zref = floor(Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).z / gr_inner.dx + C_eps) * gr_inner.dx;
	Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).zref2 = floor(Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).z / gr_inner.dx / 2. + C_eps) * gr_inner.dx * 2;
#endif
}

/**
* @brief      function to update information of nodes within one grid space of the solid boundary.
* @param[in]  ishape    order of shape.
* @param[out]  map_add_nodes    nodes newly added to those within one grid space of the solid boundary.
* @param[out]  map_remove_nodes    nodes are nolong within one grid space of the solid boundary.
*/
void Grid_Manager::update_nodes_near_solid(unsigned int ishape, D_mapint& map_add_nodes, D_mapint& map_remove_nodes)
{
	int flag0 = 0;
	unsigned int ilevel = C_max_level;
	D_mapint map_temp_remove = {}, map_temp_add = {};

	Infor_near_solid node_temp;
	node_temp.icount = 1;
	D_real x_dis, y_dis;
	D_uint xint, yint;
	D_morton key_temp, key_temp1, key_temp2;
#if(C_DIMS==3)
	D_real z_dis;
	D_uint zint;
	D_morton key_temp3;
#endif
	for (unsigned int ipoint = 0; ipoint < Solid_Manager::pointer_me->shape_solids.at(ishape).numb_nodes; ++ipoint)
	{
		// check if solid point move to another cell
		x_dis = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).x - Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).xref;
		y_dis = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).y - Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).yref;
		bool bool_xmove = ((x_dis >= gr_inner.dx) || (x_dis < 0)), bool_ymove = ((y_dis >= gr_inner.dx) || (y_dis < 0));
#if(C_DIMS==3)
		z_dis = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).z - Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).zref;
		bool bool_zmove = ((z_dis >= gr_inner.dx) || (z_dis < 0));
#endif


		if (bool_xmove && ((x_dis >= 2 * gr_inner.dx) || (x_dis < -gr_inner.dx)))
		{
			std::stringstream error;
			error << "ipoint = " << ipoint << " of ishape= " << ishape << " moves more than 1 grid space (dx) in x direcion" << std::endl;
			log_error(error.str(), Log_function::logfile);
		}
		if (bool_ymove && ((y_dis >= 2 * gr_inner.dx) || (y_dis < -gr_inner.dx)))
		{
			std::stringstream error;
			error << "ipoint = " << ipoint << " of ishape= " << ishape << " moves more than 1 grid space (dx) in y direcion" << std::endl;
			log_error(error.str(), Log_function::logfile);
		}
#if (C_DIMS==3)
		if (bool_zmove && ((z_dis >= 2 * gr_inner.dx) || (z_dis < -gr_inner.dx)))
		{
			std::stringstream error;
			error << "ipoint = " << ipoint << " of ishape= " << ishape << " moves more than 1 grid space (dx) in z direcion" << std::endl;
			log_error(error.str(), Log_function::logfile);
		}
#endif

#if (C_DIMS==2)
		if (bool_xmove || bool_ymove)
#endif
#if (C_DIMS==3)
			if (bool_xmove || bool_ymove || bool_zmove)
#endif
			{
				// find if nodes need to be removed
				xint = static_cast<D_uint>((Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).xref) / gr_inner.dx + C_eps);
				yint = static_cast<D_uint>((Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).yref) / gr_inner.dx + C_eps);
#if (C_DIMS==2)
				key_temp = static_cast <D_morton>(Morton_Assist::pointer_me->morton_encode(xint, yint));
#endif
#if (C_DIMS==3)
				zint = static_cast<D_uint>((Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).zref) / gr_inner.dx + C_eps);
				key_temp = static_cast <D_morton>(Morton_Assist::pointer_me->morton_encode(xint, yint, zint));
#endif	
				// (x0, y0 ,z0) i.e. xref, yref, zref
				if (map_temp_remove.find(key_temp) == map_temp_remove.end())
				{
					map_temp_remove.insert({ key_temp , flag0 });
				}
				--map_node_near_solid[key_temp].icount;
#if (C_DIMS==3)
				// (x0, y0 ,z1)
				key_temp3 = Morton_Assist::find_z1(key_temp, ilevel);
				if (map_temp_remove.find(key_temp3) == map_temp_remove.end())
				{
					map_temp_remove.insert({ key_temp3 , flag0 });
				}
				--map_node_near_solid[key_temp3].icount;
#endif

				// (x1, y0 ,z0)
				key_temp1 = Morton_Assist::find_x1(key_temp, ilevel);
				if (map_temp_remove.find(key_temp1) == map_temp_remove.end())
				{
					map_temp_remove.insert({ key_temp1 , flag0 });
				}
				--map_node_near_solid[key_temp1].icount;
#if (C_DIMS==3)
				// (x1, y0 ,z1)
				key_temp3 = Morton_Assist::find_z1(key_temp1, ilevel);
				if (map_temp_remove.find(key_temp3) == map_temp_remove.end())
				{
					map_temp_remove.insert({ key_temp3 , flag0 });
				}
				--map_node_near_solid[key_temp3].icount;
#endif

				// (x1, y1 ,z0)
				key_temp2 = Morton_Assist::find_y1(key_temp1, ilevel);
				if (map_temp_remove.find(key_temp2) == map_temp_remove.end())
				{
					map_temp_remove.insert({ key_temp2 , flag0 });
				}
				--map_node_near_solid[key_temp2].icount;
#if (C_DIMS==3)
				// (x1, y1 ,z1)
				key_temp3 = Morton_Assist::find_z1(key_temp2, ilevel);
				if (map_temp_remove.find(key_temp3) == map_temp_remove.end())
				{
					map_temp_remove.insert({ key_temp3 , flag0 });
				}
				--map_node_near_solid[key_temp3].icount;
#endif

				// (x0, y1 ,z0)
				key_temp1 = Morton_Assist::find_y1(key_temp, ilevel);
				if (map_temp_remove.find(key_temp1) == map_temp_remove.end())
				{
					map_temp_remove.insert({ key_temp1 , flag0 });
				}
				--map_node_near_solid[key_temp1].icount;
#if (C_DIMS==3)
				// (x0, y1 ,z1)
				key_temp3 = Morton_Assist::find_z1(key_temp1, ilevel);
				if (map_temp_remove.find(key_temp3) == map_temp_remove.end())
				{
					map_temp_remove.insert({ key_temp3 , flag0 });
				}
				--map_node_near_solid[key_temp3].icount;
#endif

				// find if nodes need to be added
				if (bool_xmove)
				{
					if (x_dis >= gr_inner.dx)
					{
						++xint;
					}
					else
					{
						--xint;
					}
				}
				if (bool_ymove)
				{
					if (y_dis >= gr_inner.dx)
					{
						++yint;
					}
					else
					{
						--yint;
					}
				}
				//xint = static_cast<D_uint>((Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).x) / gr_inner.dx + C_eps);
				//yint = static_cast<D_uint>((Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).y) / gr_inner.dx + C_eps);
#if (C_DIMS==2)
				key_temp = static_cast <D_morton>(Morton_Assist::pointer_me->morton_encode(xint, yint));
#endif
#if (C_DIMS==3)
				if (bool_zmove)
				{
					if (z_dis >= gr_inner.dx)
					{
						++zint;
					}
					else
					{
						--zint;
					}
				}
				//zint = static_cast<D_uint>((Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).z) / gr_inner.dx + C_eps);
				key_temp = static_cast <D_morton>(Morton_Assist::pointer_me->morton_encode(xint, yint, zint));
#endif	
				// (x0, y0 ,z0)
				if (map_node_near_solid.find(key_temp) == map_node_near_solid.end())
				{
					map_node_near_solid.insert({ key_temp, node_temp });
					map_add_nodes.insert({ key_temp , flag0 });
					map_temp_add.insert({ key_temp , flag0 });
				}
				else
				{
					++map_node_near_solid[key_temp].icount;
				}

#if (C_DIMS==3)
				// (x0, y0 ,z1)
				key_temp3 = Morton_Assist::find_z1(key_temp, ilevel);
				if (map_node_near_solid.find(key_temp3) == map_node_near_solid.end())
				{
					map_node_near_solid.insert({ key_temp3, node_temp });
					map_add_nodes.insert({ key_temp3 , flag0 });
					map_temp_add.insert({ key_temp3 , flag0 });
				}
				else
				{
					++map_node_near_solid[key_temp3].icount;
				}

#endif

				// (x1, y0 ,z0)
				key_temp1 = Morton_Assist::find_x1(key_temp, ilevel);
				if (map_node_near_solid.find(key_temp1) == map_node_near_solid.end())
				{
					map_node_near_solid.insert({ key_temp1, node_temp });
					map_add_nodes.insert({ key_temp1 , flag0 });
					map_temp_add.insert({ key_temp1 , flag0 });
				}
				else
				{
					++map_node_near_solid[key_temp1].icount;
				}

#if (C_DIMS==3)
				// (x1, y0 ,z1)
				key_temp3 = Morton_Assist::find_z1(key_temp1, ilevel);
				if (map_node_near_solid.find(key_temp3) == map_node_near_solid.end())
				{
					map_node_near_solid.insert({ key_temp3, node_temp });
					map_add_nodes.insert({ key_temp3 , flag0 });
					map_temp_add.insert({ key_temp3 , flag0 });
				}
				else
				{
					++map_node_near_solid[key_temp3].icount;
				}
#endif

				// (x1, y1 ,z0)
				key_temp2 = Morton_Assist::find_y1(key_temp1, ilevel);
				if (map_node_near_solid.find(key_temp2) == map_node_near_solid.end())
				{
					map_node_near_solid.insert({ key_temp2, node_temp });
					map_add_nodes.insert({ key_temp2 , flag0 });
					map_temp_add.insert({ key_temp2 , flag0 });
				}
				else
				{
					++map_node_near_solid[key_temp2].icount;
				}
#if (C_DIMS==3)
				// (x1, y1 ,z1)
				key_temp3 = Morton_Assist::find_z1(key_temp2, ilevel);
				if (map_node_near_solid.find(key_temp3) == map_node_near_solid.end())
				{
					map_node_near_solid.insert({ key_temp3, node_temp });
					map_add_nodes.insert({ key_temp3 , flag0 });
					map_temp_add.insert({ key_temp3 , flag0 });
				}
				else
				{
					++map_node_near_solid[key_temp3].icount;
				}
#endif

				// (x0, y1 ,z0)
				key_temp1 = Morton_Assist::find_y1(key_temp, ilevel);
				if (map_node_near_solid.find(key_temp1) == map_node_near_solid.end())
				{
					map_node_near_solid.insert({ key_temp1, node_temp });
					map_add_nodes.insert({ key_temp1 , flag0 });
					map_temp_add.insert({ key_temp1 , flag0 });
				}
				else
				{
					++map_node_near_solid[key_temp1].icount;
				}
#if (C_DIMS==3)
				// (x0, y1 ,z1)
				key_temp3 = Morton_Assist::find_z1(key_temp1, ilevel);
				if (map_node_near_solid.find(key_temp3) == map_node_near_solid.end())
				{
					map_node_near_solid.insert({ key_temp3, node_temp });
					map_add_nodes.insert({ key_temp3 , flag0 });
					map_temp_add.insert({ key_temp3 , flag0 });
				}
				else
				{
					++map_node_near_solid[key_temp3].icount;
				}
#endif
			}

		// reset xref, yref and zref
		if (bool_xmove)
		{
			if (x_dis >= gr_inner.dx)
			{
				Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).xref += gr_inner.dx;
			}
			else
			{
				Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).xref -= gr_inner.dx;
			}
		}
		if (bool_ymove)
		{
			if (y_dis >= gr_inner.dx)
			{
				Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).yref += gr_inner.dx;
			}
			else
			{
				Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).yref -= gr_inner.dx;
			}
		}
#if (C_DIMS==3)
		if (bool_zmove)
		{
			if (z_dis >= gr_inner.dx)
			{
				Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).zref += gr_inner.dx;
			}
			else
			{
				Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).zref -= gr_inner.dx;
			}
		}
#endif
		// if (bool_xmove) Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).xref = floor(Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).x / gr_inner.dx + C_eps) * gr_inner.dx;
		// if (bool_ymove) Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).yref = floor(Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).y / gr_inner.dx + C_eps) * gr_inner.dx;
		// if (bool_zmove) Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).zref = floor(Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).z / gr_inner.dx + C_eps) * gr_inner.dx;
	}


	// remove nodes
	bool bool_N_find_flag = true;
	D_mapint map_remove_remain = {};
	for (D_mapint::iterator iter = map_temp_remove.begin(); iter != map_temp_remove.end(); ++iter)
	{
		if (map_node_near_solid[iter->first].icount <= 0)
		{
			map_node_near_solid.erase(iter->first);
			if (map_remove_nodes.find(iter->first) == map_remove_nodes.end())
			{
				map_remove_nodes.insert({ iter->first, iter->second });
			}

			//set flag of current node as its neighbour
			bool_N_find_flag = true;
			key_temp = Morton_Assist::find_x0(iter->first, ilevel);
			if (gr_inner.grid[key_temp].flag != flag_near_solid)
			{
				gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
				bool_N_find_flag = false;
			}
			if (bool_N_find_flag)
			{
				key_temp = Morton_Assist::find_x1(iter->first, ilevel);
				if (gr_inner.grid[key_temp].flag != flag_near_solid)
				{
					gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
					bool_N_find_flag = false;
				}
			}
			if (bool_N_find_flag)
			{
				key_temp = Morton_Assist::find_y0(iter->first, ilevel);
				if (gr_inner.grid[key_temp].flag != flag_near_solid)
				{
					gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
					bool_N_find_flag = false;
				}
			}
			if (bool_N_find_flag)
			{
				key_temp = Morton_Assist::find_y1(iter->first, ilevel);
				if (gr_inner.grid[key_temp].flag != flag_near_solid)
				{
					gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
					bool_N_find_flag = false;
				}
			}
#if(C_DIMS == 3)
			if (bool_N_find_flag)
			{
				key_temp = Morton_Assist::find_z0(iter->first, ilevel);
				if (gr_inner.grid[key_temp].flag != flag_near_solid)
				{
					gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
					bool_N_find_flag = false;
				}
			}
			if (bool_N_find_flag)
			{
				key_temp = Morton_Assist::find_z1(iter->first, ilevel);
				if (gr_inner.grid[key_temp].flag != flag_near_solid)
				{
					gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
					bool_N_find_flag = false;
				}
			}
#endif
			if (bool_N_find_flag)
			{
				map_remove_remain.insert({ iter->first , flag0 });
			}
		}
	}

	// set flag of removed nodes
	D_mapint map_remove_remain_temp = {};
	D_mapint map_remove_flaged = {};
	for (D_mapint::iterator iter = map_remove_remain.begin(); iter != map_remove_remain.end(); ++iter)
	{
		bool_N_find_flag = true;
		key_temp = Morton_Assist::find_x0(iter->first, ilevel);
		if (gr_inner.grid[key_temp].flag != flag_near_solid)
		{
			gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
			map_remove_flaged.insert({ iter->first , flag0 });
			bool_N_find_flag = false;
		}
		if (bool_N_find_flag)
		{
			key_temp = Morton_Assist::find_x1(iter->first, ilevel);
			if (gr_inner.grid[key_temp].flag != flag_near_solid)
			{
				gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
				map_remove_flaged.insert({ iter->first , flag0 });
				bool_N_find_flag = false;
			}
		}
		if (bool_N_find_flag)
		{
			key_temp = Morton_Assist::find_y0(iter->first, ilevel);
			if (gr_inner.grid[key_temp].flag != flag_near_solid)
			{
				gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
				map_remove_flaged.insert({ iter->first , flag0 });
				bool_N_find_flag = false;
			}
		}
		if (bool_N_find_flag)
		{
			key_temp = Morton_Assist::find_y1(iter->first, ilevel);
			if (gr_inner.grid[key_temp].flag != flag_near_solid)
			{
				gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
				map_remove_flaged.insert({ iter->first , flag0 });
				bool_N_find_flag = false;
			}
		}
#if(C_DIMS == 3)
		if (bool_N_find_flag)
		{
			key_temp = Morton_Assist::find_z0(iter->first, ilevel);
			if (gr_inner.grid[key_temp].flag != flag_near_solid)
			{
				gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
				map_remove_flaged.insert({ iter->first , flag0 });
				bool_N_find_flag = false;
			}
		}
		if (bool_N_find_flag)
		{
			key_temp = Morton_Assist::find_z1(iter->first, ilevel);
			if (gr_inner.grid[key_temp].flag != flag_near_solid)
			{
				gr_inner.grid[iter->first].flag = gr_inner.grid[key_temp].flag;
				map_remove_flaged.insert({ iter->first , flag0 });
				bool_N_find_flag = false;
			}
		}
#endif
		if (bool_N_find_flag)
		{
			map_remove_remain_temp.insert({ iter->first , flag0 });
		}
	}

	int i = 0;
	while (!map_remove_flaged.empty() && i < iter_total)
	{
		D_mapint map_remove_flaged_temp = {};
		for (D_mapint::iterator iter = map_remove_flaged.begin(); iter != map_remove_flaged.end(); ++iter)
		{
			key_temp = Morton_Assist::find_x0(iter->first, ilevel);
			if (map_remove_remain_temp.find(key_temp) != map_remove_remain_temp.end())
			{
				gr_inner.grid[key_temp].flag = gr_inner.grid[iter->first].flag;
				map_remove_remain_temp.erase(key_temp);
				map_remove_flaged_temp.insert({ key_temp, flag0 });
			}
			key_temp = Morton_Assist::find_x1(iter->first, ilevel);
			if (map_remove_remain_temp.find(key_temp) != map_remove_remain_temp.end())
			{
				gr_inner.grid[key_temp].flag = gr_inner.grid[iter->first].flag;
				map_remove_remain_temp.erase(key_temp);
				map_remove_flaged_temp.insert({ key_temp, flag0 });
			}
			key_temp = Morton_Assist::find_y0(iter->first, ilevel);
			if (map_remove_remain_temp.find(key_temp) != map_remove_remain_temp.end())
			{
				gr_inner.grid[key_temp].flag = gr_inner.grid[iter->first].flag;
				map_remove_remain_temp.erase(key_temp);
				map_remove_flaged_temp.insert({ key_temp, flag0 });
			}
			key_temp = Morton_Assist::find_y1(iter->first, ilevel);
			if (map_remove_remain_temp.find(key_temp) != map_remove_remain_temp.end())
			{
				gr_inner.grid[key_temp].flag = gr_inner.grid[iter->first].flag;
				map_remove_remain_temp.erase(key_temp);
				map_remove_flaged_temp.insert({ key_temp, flag0 });
			}
#if (C_DIMS==3)
			key_temp = Morton_Assist::find_z0(iter->first, ilevel);
			if (map_remove_remain_temp.find(key_temp) != map_remove_remain_temp.end())
			{
				gr_inner.grid[key_temp].flag = gr_inner.grid[iter->first].flag;
				map_remove_remain_temp.erase(key_temp);
				map_remove_flaged_temp.insert({ key_temp, flag0 });
			}
			key_temp = Morton_Assist::find_z1(iter->first, ilevel);
			if (map_remove_remain_temp.find(key_temp) != map_remove_remain_temp.end())
			{
				gr_inner.grid[key_temp].flag = gr_inner.grid[iter->first].flag;
				map_remove_remain_temp.erase(key_temp);
				map_remove_flaged_temp.insert({ key_temp, flag0 });
			}
#endif
		}
		map_remove_flaged.clear();
		map_remove_flaged = map_remove_flaged_temp;
		++i;
	}

	for (D_mapint::iterator iter = map_remove_remain_temp.begin(); iter != map_remove_remain_temp.end(); ++iter)
	{
		gr_inner.grid[iter->first].flag = flag_ghost;
	}


	if (i == iter_total)
	{
		std::stringstream warning;
		warning << "process of flag removed nodes near the solid exceed the maximum number of iteration: " << iter_total << ", try to incease iter_total to solve this" << std::endl;
		log_warning(warning.str(), Log_function::logfile);
	}

	// set flag of added nodes as flag_near_solid
	for (D_mapint::iterator iter = map_temp_add.begin(); iter != map_temp_add.end(); ++iter)
	{
		gr_inner.grid[iter->first].flag = flag_near_solid;
	}
}

#if	(C_FSI_INTERFACE == 1)
/**
* @brief      function add and remove elements in map_node_IB.
* @param[in]  ilevel           refinement level.
* @param[in]  map_add_nodes_in    nodes used to search if there are nodes needed to be added at ilevel.
* @param[in]  map_remove_nodes_in    nodes used to search if there are nodes needed to be removed at ilevel.
*/
void Grid_Manager::update_map_node_IB(const unsigned int ilevel, D_mapint& map_add_nodes_in, D_mapint& map_remove_nodes_in)
{
	Infor_IB node_temp;
	node_temp.icount = 1;
	// search for nodes used for IBM
	if (C_extend_IB == 1)
	{
		for (D_mapint::iterator iter = map_add_nodes_in.begin(); iter != map_add_nodes_in.end(); ++iter)
		{
			if (map_node_IB.find(iter->first) == map_node_IB.end())
			{
				map_node_IB.insert({ iter->first, node_temp });
			}
		}
		for (D_mapint::iterator iter = map_remove_nodes_in.begin(); iter != map_remove_nodes_in.end(); ++iter)
		{
			if (map_node_IB.find(iter->first) != map_node_IB.end())
			{
				map_node_IB.erase(iter->first);
			}
		}
	}
	else if (C_extend_IB > 1)
	{
		unsigned int extend_x0 = C_extend_IB - 1;
		unsigned int extend_x1 = C_extend_IB - 1;
		unsigned int extend_y0 = C_extend_IB - 1;
		unsigned int extend_y1 = C_extend_IB - 1;
#if (C_DIMS == 3)
		unsigned int extend_z0 = C_extend_IB - 1;
		unsigned int extend_z1 = C_extend_IB - 1;
#endif

		unsigned int extend_temp_x0 = extend_x0, extend_temp_x1 = extend_x1, extend_temp_y0 = extend_y0, extend_temp_y1 = extend_y1;
#if (C_DIMS==3)
		unsigned int extend_temp_z0 = extend_z0, extend_temp_z1 = extend_z1;
#endif
		D_morton morton_temp_x0, morton_temp_x1, morton_temp_y0, morton_temp_y1, morton_temp, morton_temp0, morton_temp1;
#if (C_DIMS == 3)
		D_morton morton_temp_z0, morton_temp_z1, morton_temp2;
#endif
		// update map_node_IB due to adding
		for (D_mapint::iterator iter = map_add_nodes_in.begin(); iter != map_add_nodes_in.end(); ++iter)
		{
			D_real x, y;
#if (C_DIMS == 2)
			Morton_Assist::compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
			D_real z;
			Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
			extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
			extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif
			if ((x - extend_x0 * gr_inner.dx) < (C_dx * C_x0b_offset))
			{
				extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / gr_inner.dx + C_eps);
			}
			if ((x + extend_x1 * gr_inner.dx) > xb_domain)
			{
				extend_temp_x1 = static_cast<unsigned int>((xb_domain - x) / gr_inner.dx + C_eps);
			}
			if ((y - extend_y0 * gr_inner.dx) < (C_dx * C_y0b_offset))
			{
				extend_temp_y0 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / gr_inner.dx + C_eps);
			}
			if ((y + extend_y1 * gr_inner.dx) > yb_domain)
			{
				extend_temp_y1 = static_cast<unsigned int>((yb_domain - y) / gr_inner.dx + C_eps);
			}
#if (C_DIMS == 3)
			if ((z - extend_z0 * gr_inner.dx) < (C_dx * C_z0b_offset))
			{
				extend_temp_z0 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / gr_inner.dx + C_eps);
			}
			if ((z + extend_z1 * gr_inner.dx) > zb_domain)
			{
				extend_temp_z1 = static_cast<unsigned int>((zb_domain - z) / gr_inner.dx + C_eps);
			}
#endif
#endif
#if (C_DIMS == 2)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx - extend_temp_x0 + C_eps), static_cast<D_int>(y / gr_inner.dx - extend_temp_y0 + C_eps);
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx - extend_temp_x0 + C_eps), static_cast<D_int>(y / gr_inner.dx - extend_temp_y0 + C_eps), static_cast<D_int>(z / gr_inner.dx - extend_temp_z0 + C_eps));
			morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
			for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
			{
				morton_temp1 = morton_temp2;
#endif
				for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
				{
					morton_temp = morton_temp1;
					for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
					{
						if (map_node_IB.find(morton_temp) == map_node_IB.end())
						{
							map_node_IB.insert({ morton_temp, node_temp });
						}
						else
						{
							++map_node_IB[morton_temp].icount;
						}
						morton_temp = Morton_Assist::find_x1(morton_temp, ilevel);
					}
					morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel);
				}
#if (C_DIMS == 3)
				morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
			}
#endif		
		}

		// update map_node_IB due to removing
		for (D_mapint::iterator iter = map_remove_nodes_in.begin(); iter != map_remove_nodes_in.end(); ++iter)
		{
			D_real x, y;
#if (C_DIMS == 2)
			Morton_Assist::compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
			D_real z;
			Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
			extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
			extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif
			if ((x - extend_x0 * gr_inner.dx) < (C_dx * C_x0b_offset))
			{
				extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / gr_inner.dx + C_eps);
			}
			if ((x + extend_x1 * gr_inner.dx) > xb_domain)
			{
				extend_temp_x1 = static_cast<unsigned int>((xb_domain - x) / gr_inner.dx + C_eps);
			}
			if ((y - extend_y0 * gr_inner.dx) < (C_dx * C_y0b_offset))
			{
				extend_temp_y0 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / gr_inner.dx + C_eps);
			}
			if ((y + extend_y1 * gr_inner.dx) > yb_domain)
			{
				extend_temp_y1 = static_cast<unsigned int>((yb_domain - y) / gr_inner.dx + C_eps);
			}
#if (C_DIMS == 3)
			if ((z - extend_z0 * gr_inner.dx) < (C_dx * C_z0b_offset))
			{
				extend_temp_z0 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / gr_inner.dx + C_eps);
			}
			if ((z + extend_z1 * gr_inner.dx) > zb_domain)
			{
				extend_temp_z1 = static_cast<unsigned int>((zb_domain - z) / gr_inner.dx + C_eps);
			}
#endif
#endif
#if (C_DIMS == 2)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx - extend_temp_x0 + C_eps), static_cast<D_int>(y / gr_inner.dx - extend_temp_y0 + C_eps);
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx - extend_temp_x0 + C_eps), static_cast<D_int>(y / gr_inner.dx - extend_temp_y0 + C_eps), static_cast<D_int>(z / gr_inner.dx - extend_temp_z0 + C_eps));
			morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
			for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
			{
				morton_temp1 = morton_temp2;
#endif
				for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
				{
					morton_temp = morton_temp1;
					for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
					{
						if (map_node_IB.find(morton_temp) != map_node_IB.end())
						{
							--map_node_IB[morton_temp].icount;
							if (map_node_IB[morton_temp].icount <= 0)
							{
								map_node_IB.erase(morton_temp);
							}
						}
						morton_temp = Morton_Assist::find_x1(morton_temp, ilevel);
					}
					morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel);
				}
#if (C_DIMS == 3)
				morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
			}
#endif		
		}
	}
}
#endif

/**
* @brief      add and remove ghost nodes.
* @param[in]  ilevel           refinement level.
* @param[in]  map_add_nodes_in    nodes used to search if there are nodes needed to be added at ilevel.
* @param[in]  map_remove_nodes_in    nodes used to search if there are nodes needed to be removed at ilevel.
*/
void Grid_Manager::update_ghost_node(const unsigned int ilevel, D_mapint& map_add_nodes_in, D_mapint& map_remove_nodes_in)
{
	unsigned int extend_x0 = C_extend_ghost;
	unsigned int extend_x1 = C_extend_ghost;
	unsigned int extend_y0 = C_extend_ghost;
	unsigned int extend_y1 = C_extend_ghost;
#if (C_DIMS == 3)
	unsigned int extend_z0 = C_extend_ghost;
	unsigned int extend_z1 = C_extend_ghost;
#endif
	D_morton morton_temp, morton_temp0, morton_temp1;
#if (C_DIMS==3)
	D_morton morton_temp2;
#endif

	unsigned int extend_temp_x0 = extend_x0, extend_temp_x1 = extend_x1, extend_temp_y0 = extend_y0, extend_temp_y1 = extend_y1;
#if (C_DIMS==3)
	unsigned int extend_temp_z0 = extend_z0, extend_temp_z1 = extend_z1;
#endif
	Node_IB node_temp;
	node_temp.flag = flag_ghost;
	node_temp.icount_ghost = 1;
	for (D_mapint::iterator iter = map_add_nodes_in.begin(); iter != map_add_nodes_in.end(); ++iter)
	{
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
		extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
		extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif
		if ((x - extend_temp_x0 * gr_inner.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / gr_inner.dx + C_eps);
		}
		if ((x + extend_x1 * gr_inner.dx) > xb_domain)
		{
			extend_temp_x1 = static_cast<unsigned int>((xb_domain - x) / gr_inner.dx + C_eps);
		}
		if ((y - extend_y0 * gr_inner.dx) < (C_dx * C_y0b_offset))
		{
			extend_temp_y0 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / gr_inner.dx + C_eps);
		}
		if ((y + extend_y1 * gr_inner.dx) > yb_domain)
		{
			extend_temp_y1 = static_cast<unsigned int>((yb_domain - y) / gr_inner.dx + C_eps);
		}
#if (C_DIMS == 3)
		if ((z - extend_z0 * gr_inner.dx) < (C_dx * C_z0b_offset))
		{
			extend_temp_z0 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / gr_inner.dx + C_eps);
		}
		if ((z + extend_z1 * gr_inner.dx) > zb_domain)
		{
			extend_temp_z1 = static_cast<unsigned int>((zb_domain - z) / gr_inner.dx + C_eps);
		}
#endif
#endif
#if (C_DIMS == 2)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx - extend_temp_x0 + C_eps), static_cast<D_int>(y / gr_inner.dx - extend_temp_y0 + C_eps);
		morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx - extend_temp_x0 + C_eps), static_cast<D_int>(y / gr_inner.dx - extend_temp_y0 + C_eps), static_cast<D_int>(z / gr_inner.dx - extend_temp_z0 + C_eps));
		morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
		for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp1 = morton_temp2;
#endif
			for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
			{
				morton_temp = morton_temp1;
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (gr_inner.grid.find(morton_temp) == gr_inner.grid.end())
					{
						gr_inner.grid.insert({ morton_temp , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel)[0].find(morton_temp) == bk_boundary_x.at(ilevel)[0].end())
							{
								bk_boundary_x.at(ilevel)[0].insert({ morton_temp , &gr_inner.grid[morton_temp] });
							}
						}
						else if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel)[1].find(morton_temp) == bk_boundary_x.at(ilevel)[1].end())
							{
								bk_boundary_x.at(ilevel)[1].insert({ morton_temp , &gr_inner.grid[morton_temp] });
							}
						}
						if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel)[0].find(morton_temp) == bk_boundary_y.at(ilevel)[0].end())
							{
								bk_boundary_y.at(ilevel)[0].insert({ morton_temp , &gr_inner.grid[morton_temp] });
							}
						}
						else if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel)[1].find(morton_temp) == bk_boundary_y.at(ilevel)[1].end())
							{
								bk_boundary_y.at(ilevel)[1].insert({ morton_temp , &gr_inner.grid[morton_temp] });
							}
						}
#if(C_DIMS == 3)
						if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel)[0].find(morton_temp) == bk_boundary_z.at(ilevel)[0].end())
							{
								bk_boundary_z.at(ilevel)[0].insert({ morton_temp , &gr_inner.grid[morton_temp] });
							}
						}
						else if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel)[1].find(morton_temp) == bk_boundary_z.at(ilevel)[1].end())
							{
								bk_boundary_z.at(ilevel)[1].insert({ morton_temp , &gr_inner.grid[morton_temp] });
							}
						}
#endif
#endif
					}
					else
					{
						++gr_inner.grid[morton_temp].icount_ghost;
					}
					morton_temp = Morton_Assist::find_x1(morton_temp, ilevel);
				}
				morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel);
			}
#if (C_DIMS == 3)
			morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
		}
#endif		
	}

	for (D_mapint::iterator iter = map_remove_nodes_in.begin(); iter != map_remove_nodes_in.end(); ++iter)
	{
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
		extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
		extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif
		if ((x - extend_temp_x0 * gr_inner.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / gr_inner.dx + C_eps);
		}
		if ((x + extend_x1 * gr_inner.dx) > xb_domain)
		{
			extend_temp_x1 = static_cast<unsigned int>((xb_domain - x) / gr_inner.dx + C_eps);
		}
		if ((y - extend_y0 * gr_inner.dx) < (C_dx * C_y0b_offset))
		{
			extend_temp_y0 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / gr_inner.dx + C_eps);
		}
		if ((y + extend_y1 * gr_inner.dx) > yb_domain)
		{
			extend_temp_y1 = static_cast<unsigned int>((yb_domain - y) / gr_inner.dx + C_eps);
		}
#if (C_DIMS == 3)
		if ((z - extend_z0 * gr_inner.dx) < (C_dx * C_z0b_offset))
		{
			extend_temp_z0 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / gr_inner.dx + C_eps);
		}
		if ((z + extend_z1 * gr_inner.dx) > zb_domain)
		{
			extend_temp_z1 = static_cast<unsigned int>((zb_domain - z) / gr_inner.dx + C_eps);
		}
#endif
#endif
#if (C_DIMS == 2)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx - extend_temp_x0 + C_eps), static_cast<D_int>(y / gr_inner.dx - extend_temp_y0 + C_eps);
		morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx - extend_temp_x0 + C_eps), static_cast<D_int>(y / gr_inner.dx - extend_temp_y0 + C_eps), static_cast<D_int>(z / gr_inner.dx - extend_temp_z0 + C_eps));
		morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
		for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp1 = morton_temp2;
#endif
			for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
			{
				morton_temp = morton_temp1;
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (gr_inner.grid.find(morton_temp) != gr_inner.grid.end())
					{
						--gr_inner.grid[morton_temp].icount_ghost;
						if ((gr_inner.grid[morton_temp].icount_ghost <= 0) && (gr_inner.grid[morton_temp].flag == flag_ghost))
						{
							gr_inner.grid.erase(morton_temp);
#if (C_CHECK_MORTON_BOUNDARY == 1)
							if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
							{
								if (bk_boundary_x.at(ilevel)[0].find(morton_temp) != bk_boundary_x.at(ilevel)[0].end())
								{
									bk_boundary_x.at(ilevel)[0].erase(morton_temp);
								}
							}
							else if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
							{
								if (bk_boundary_x.at(ilevel)[1].find(morton_temp) != bk_boundary_x.at(ilevel)[1].end())
								{
									bk_boundary_x.at(ilevel)[1].erase(morton_temp);
								}
							}
							if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
							{
								if (bk_boundary_y.at(ilevel)[0].find(morton_temp) != bk_boundary_y.at(ilevel)[0].end())
								{
									bk_boundary_y.at(ilevel)[0].erase(morton_temp);
								}
							}
							else if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
							{
								if (bk_boundary_y.at(ilevel)[1].find(morton_temp) != bk_boundary_y.at(ilevel)[1].end())
								{
									bk_boundary_y.at(ilevel)[1].erase(morton_temp);
								}
							}
#if(C_DIMS == 3)
							if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
							{
								if (bk_boundary_z.at(ilevel)[0].find(morton_temp) != bk_boundary_z.at(ilevel)[0].end())
								{
									bk_boundary_z.at(ilevel)[0].erase(morton_temp);
								}
							}
							else if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
							{
								if (bk_boundary_z.at(ilevel)[1].find(morton_temp) != bk_boundary_z.at(ilevel)[1].end())
								{
									bk_boundary_z.at(ilevel)[1].erase(morton_temp);
								}
							}
#endif
#endif
						}
					}
					morton_temp = Morton_Assist::find_x1(morton_temp, ilevel);
				}
				morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel);
			}
#if (C_DIMS == 3)
			morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
		}
#endif		
	}
}

/**
* @brief      funciton to call functions with the same number of parameters but different data type, not sure if it's the only way to use tamplate to achieve what I want.
* @param[in]  ilevel           refinement level.
* @param[in]  map_add_nodes_in    nodes used to search if there are nodes needed to be added at ilevel.
* @param[in]  map_remove_nodes_in    nodes used to search if there are nodes needed to be removed at ilevel.
* @param[out]  map_add_nodes_out    nodes used to search if there are nodes needed to be added at (ilevel - 1).
* @param[out]  map_remove_nodes_out    nodes used to search if there are nodes needed to be removed at (ilevel - 1).
*/
D_mapint map_add2 = {}, map_remove2 = {};
void Grid_Manager::call_update_nodes(const unsigned int ilevel, D_mapint& map_add_nodes_in, D_mapint& map_remove_nodes_in, D_mapint& map_add_nodes_out, D_mapint& map_remove_nodes_out)
{

	if (ilevel == C_max_level)
	{
		Node_IB node_temp;

		//map_add_nodes_in.clear();
		//map_remove_nodes_in.clear();

		//D_mapint map_boundary = {}, map_add = {}, map_remove = {};
		//map_boundary = Grid_Manager::pointer_me->gr_NoIB.at(3).coarse2fine.at(1);
		//if (map_add_nodes_out.size() < 1)
		//{
		//	map_add2.clear();
		//}
		//if (map_remove_nodes_out.size() < 1)
		//{
		//	map_remove2.clear();
		//}

		update_nodes(ilevel, Grid_Manager::pointer_me->gr_inner, node_temp, map_add_nodes_in, map_remove_nodes_in, map_add_nodes_out, map_remove_nodes_out);

		// delete nodes added and removed in funciton update_nodes
		std::vector<D_morton> map_duplicate_add = {};
		for (D_mapint::iterator iter = map_add_nodes_out.begin(); iter != map_add_nodes_out.end(); ++iter)
		{
			if (map_remove_nodes_out.find(iter->first) != map_remove_nodes_out.end())
			{
				map_duplicate_add.push_back(iter->first);
				map_remove_nodes_out.erase(iter->first);
			}
		}

		for (std::vector<D_morton>::iterator iter = map_duplicate_add.begin(); iter != map_duplicate_add.end(); ++iter)
		{
			if (map_add_nodes_out.find(*iter) != map_add_nodes_out.end())
			{
				map_add_nodes_out.erase(*iter);
			}
		}


		//for (auto iter = map_boundary.begin(); iter != map_boundary.end(); ++iter)
		//{
		//	if (Grid_Manager::pointer_me->gr_NoIB.at(3).coarse2fine.at(1).find(iter->first) == Grid_Manager::pointer_me->gr_NoIB.at(3).coarse2fine.at(1).end())
		//	{
		//		map_remove.insert({ iter->first, 1 });
		//	}
		//}
		//for (auto iter = Grid_Manager::pointer_me->gr_NoIB.at(3).coarse2fine.at(1).begin(); iter != Grid_Manager::pointer_me->gr_NoIB.at(3).coarse2fine.at(1).end(); ++iter)
		//{
		//	if (map_boundary.find(iter->first) == map_boundary.end())
		//	{
		//		map_add.insert({ iter->first, 1 });
		//	}
		//}
		//map_add2.insert(map_add.begin(), map_add.end());
		//map_remove2.insert(map_remove.begin(), map_remove.end());
		//for (auto iter = map_add_nodes_out.begin(); iter != map_add_nodes_out.end(); ++iter)
		//{
		//	//if (map_add2.find(iter->first) == map_add2.end())
		//	if (map_remove_nodes_out.find(iter->first) != map_remove_nodes_out.end())
		//	{
		//		//output_points_update.insert({ iter->first, ilevel });
		//	}		
		//}

		//std::cout << "add size: " << map_add2.size() << " " << map_add_nodes_out.size() << std::endl;
		//std::cout << "remove size: " << map_remove2.size() << " " << map_remove_nodes_out.size() << std::endl;
	}
	else
	{
		Node node_temp;
		update_nodes(ilevel, Grid_Manager::pointer_me->gr_NoIB.at(ilevel), node_temp, map_add_nodes_in, map_remove_nodes_in, map_add_nodes_out, map_remove_nodes_out);

		// delete nodes added and removed in funciton update_nodes
		std::vector<D_morton> map_duplicate_add = {};
		for (D_mapint::iterator iter = map_add_nodes_out.begin(); iter != map_add_nodes_out.end(); ++iter)
		{
			if (map_remove_nodes_out.find(iter->first) != map_remove_nodes_out.end())
			{
				map_duplicate_add.push_back(iter->first);
				map_remove_nodes_out.erase(iter->first);
			}
		}

		for (std::vector<D_morton>::iterator iter = map_duplicate_add.begin(); iter != map_duplicate_add.end(); ++iter)
		{
			if (map_add_nodes_out.find(*iter) != map_add_nodes_out.end())
			{
				map_add_nodes_out.erase(*iter);
			}
		}
	}

}

/**
* @brief function to check nodes to be added or removed.
* @param[in]  ilevel                refinement level.
* @param[in]  grid_ptr              pointer point to the class storing grid information at ilevel.
* @param[in]  map_add_nodes_in      nodes used to search if there are nodes needed to be added at ilevel, i.e. nodes at (ilevel + 1) for (grid_NIB) and nodes near solid for (grid_IB)
* @param[in]  map_remove_nodes_in   nodes used to search if there are nodes needed to be removed at ilevel, i.e. nodes at (ilevel + 1) for (grid_NIB) and nodes near solid for (grid_IB)
* @param[out]  map_add_nodes_out     nodes used to search if there are nodes needed to be added at (ilevel - 1)
* @param[out]  map_remove_nodes_out  nodes used to search if there are nodes needed to be added at (ilevel - 1)
* @todo {to avoid search backward to check if nodes need to be removed after searching surrounding boundary, count number can be stored at the boundary. the problem is how to handle block splitting when the removed nodes is not on the numerical boundary}
*/
template <class T_grid, class T_node>
void Grid_Manager::update_nodes(const unsigned int ilevel, T_grid& grid_ptr, T_node node_temp, D_mapint& map_add_nodes_in, D_mapint& map_remove_nodes_in, D_mapint& map_add_nodes_out, D_mapint& map_remove_nodes_out)
{
	//output_points_update.clear();
	Grid_NIB* grid_ptr_coarse = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel - 1));

	Timer tmr;
	double t0 = tmr.elapsed();

	const unsigned int ilevel_1 = ilevel - 1;
	D_mapint map_add_out_temp = {};
	D_mapint map_remove_out_temp = {};
	D_morton morton_temp, morton_temp_ilevel_1;
	int flag_remove = -1;
	int flag_add = 1;
	D_morton morton_temp0, morton_temp1; // 0: faces; 1: edges; 2: corners
#if(C_DIMS==3)
	D_morton morton_temp2;
#endif
	D_morton morton_xyz = Morton_Assist::morton_xyz.at(ilevel);
	D_morton morton_xyz_1 = Morton_Assist::morton_xyz_flip.at(ilevel);
	// convert nodes in (map_remove_nodes_in) at (ilevel) to (map_remove_temp) at (ilevel - 1)
	D_mapint map_remove_temp = {};

	if (ilevel == C_max_level)
	{
		for (D_mapint::iterator iter = map_remove_nodes_in.begin(); iter != map_remove_nodes_in.end(); ++iter)
		{
			morton_temp = iter->first & morton_xyz_1;
			if (map_remove_temp.find(morton_temp) == map_remove_temp.end())
			{
				map_remove_temp.insert({ morton_temp, 1 });
			}
			else
			{
				++map_remove_temp[morton_temp];
			}
		}
	}
	else
	{
		for (D_mapint::iterator iter = map_remove_nodes_in.begin(); iter != map_remove_nodes_in.end(); ++iter)
		{
			morton_temp = iter->first & morton_xyz_1;
			if (map_remove_temp.find(morton_temp) == map_remove_temp.end())
			{
				map_remove_temp.insert({ morton_temp, 1 });
			}
			else
			{
				++map_remove_temp[morton_temp];
			}
		}
	}
	// convert nodes in (map_add_nodes_in) at (ilevel) to (map_add_temp) at (ilevel - 1)
	D_mapint map_add_temp = {};
	if (ilevel == C_max_level)
	{
		for (D_mapint::iterator iter = map_add_nodes_in.begin(); iter != map_add_nodes_in.end(); ++iter)
		{
			morton_temp = iter->first & morton_xyz_1;
			if (map_add_temp.find(morton_temp) == map_add_temp.end())
			{
				map_add_temp.insert({ morton_temp, 1 });
			}
			else
			{
				++map_add_temp[morton_temp];
			}
		}
	}
	else
	{
		for (D_mapint::iterator iter = map_add_nodes_in.begin(); iter != map_add_nodes_in.end(); ++iter)
		{
			morton_temp = iter->first & morton_xyz_1;
			if (map_add_temp.find(morton_temp) == map_add_temp.end())
			{
				map_add_temp.insert({ morton_temp, 1 });
			}
			else
			{
				++map_add_temp[morton_temp];
			}
		}
	}

	// delete duplicated nodes
	std::vector<D_morton> map_duplicate_add = {};
	for (D_mapint::iterator iter = map_add_temp.begin(); iter != map_add_temp.end(); ++iter)
	{
		if (map_remove_temp.find(iter->first) != map_remove_temp.end())
		{
			if (map_add_temp[iter->first] > map_remove_temp[iter->first])
			{
				map_add_temp[iter->first] -= map_remove_temp[iter->first];
				map_remove_temp.erase(iter->first);
			}
			else if (map_add_temp[iter->first] == map_remove_temp[iter->first])
			{
				map_duplicate_add.push_back(iter->first);
				map_remove_temp.erase(iter->first);
			}
			else
			{
				map_remove_temp[iter->first] -= map_add_temp[iter->first];
				map_duplicate_add.push_back(iter->first);
			}		
		}
	}

	for (std::vector<D_morton>::iterator iter = map_duplicate_add.begin(); iter != map_duplicate_add.end(); ++iter)
	{
		if (map_add_temp.find(*iter) != map_add_temp.end())
		{
			map_add_temp.erase(*iter);
		}
	}
	// find points on the searching boundary (vector_search_point_boundary[ilevel]) whose searching center is at nodes in (map_remove_nodes_in)
	D_mapint map_remove_surroundings = {}; // points, given by (vector_search_point_boundary[ilevel]), around the nodes removed in the previous step
	unsigned int extend_x0, extend_x1, extend_y0, extend_y1;
	D_morton morton_temp_x0, morton_temp_x1, morton_temp_y0, morton_temp_y1;
#if (C_DIMS == 3)
	D_morton morton_temp_z0, morton_temp_z1;
	unsigned int extend_z0, extend_z1;
#endif
	if (ilevel == C_max_level)
	{
		extend_x0 = C_extend_inner / 2 + C_extend_inner_x0 / 2 + 1;
		extend_x1 = C_extend_inner / 2 + C_extend_inner_x1 / 2 + 1;
		extend_y0 = C_extend_inner / 2 + C_extend_inner_y0 / 2 + 1;
		extend_y1 = C_extend_inner / 2 + C_extend_inner_y1 / 2 + 1;
#if (C_DIMS == 3)
		extend_z0 = C_extend_inner / 2 + C_extend_inner_z0 / 2 + 1;
		extend_z1 = C_extend_inner / 2 + C_extend_inner_z1 / 2 + 1;
#endif
	}
	else
	{
		extend_x0 = C_extend / 2 + C_extend_outer_x0 / 2 + 1;
		extend_x1 = C_extend / 2 + C_extend_outer_x1 / 2 + 1;
		extend_y0 = C_extend / 2 + C_extend_outer_y0 / 2 + 1;
		extend_y1 = C_extend / 2 + C_extend_outer_y1 / 2 + 1;
#if (C_DIMS == 3)
		extend_z0 = C_extend / 2 + C_extend_outer_z0 / 2 + 1;
		extend_z1 = C_extend / 2 + C_extend_outer_z1 / 2 + 1;
#endif
	}

	D_mapint map_add_surroundings = {};
	unsigned int extend_temp_x0_1 = extend_x0 - 1, extend_temp_x1_1 = extend_x1 - 1, extend_temp_y0_1 = extend_y0 - 1, extend_temp_y1_1 = extend_y1 - 1;
	unsigned int extend_temp_x0 = extend_x0, extend_temp_x1 = extend_x1, extend_temp_y0 = extend_y0, extend_temp_y1 = extend_y1;
#if (C_DIMS==3)
	unsigned int extend_temp_z0 = extend_z0, extend_temp_z1 = extend_z1;
	unsigned int extend_temp_z0_1 = extend_z0 - 1, extend_temp_z1_1 = extend_z1 - 1;
#endif
	unsigned int bit_move = Morton_Assist::bit_otherlevel - ilevel * C_DIMS;

#if (C_SEARCH_METHOD == 1)
	// find nodes on the boundary of a given region for removing
	for (D_mapint::iterator iter = map_add_temp.begin(); iter != map_add_temp.end(); ++iter)
	{
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
		extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
		extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif
		if ((x - extend_x0 * grid_ptr.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((x + extend_x1 * grid_ptr.dx) > xb_domain)
		{
			extend_temp_x1 = static_cast<unsigned int>((xb_domain - x) / grid_ptr.dx + C_eps);
		}
		if ((y - extend_y0 * grid_ptr.dx) < (C_dx * C_y0b_offset))
		{
			extend_temp_y0 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((y + extend_y1 * grid_ptr.dx) > yb_domain)
		{
			extend_temp_y1 = static_cast<unsigned int>((yb_domain - y) / grid_ptr.dx + C_eps);
		}
#if (C_DIMS == 3)
		if ((z - extend_z0 * grid_ptr.dx) < (C_dx * C_z0b_offset))
		{
			extend_temp_z0 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((z + extend_z1 * grid_ptr.dx) > zb_domain)
		{
			extend_temp_z1 = static_cast<unsigned int>((zb_domain - z) / grid_ptr.dx + C_eps);
		}
#endif
#endif
#if (C_DIMS == 3)
		// (x, y) plane, negative z
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp1 = morton_temp << bit_move;
		for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end())
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		}

		// (x, y) plane, positive z
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) + extend_temp_z1 * 2);
		morton_temp1 = morton_temp << bit_move;
		for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end())
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		}

		// (x, z) plane, negative y
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end())
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (x, z) plane, positive y
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) + extend_temp_y1 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end())
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (y, z) plane, negative x
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int iy = 1; iy < extend_temp_y0 + extend_temp_y1; ++iy)
			{
				if (grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end())
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				morton_temp = Morton_Assist::find_y1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (y, z) plane, positve x
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) + extend_temp_x1 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int iy = 1; iy < extend_temp_y0 + extend_temp_y1; ++iy)
			{
				if (grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end())
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				morton_temp = Morton_Assist::find_y1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}
#endif
	}

	// find nodes on the boundary of a given region for removing (not complete 2D case for searching on the boundary)
	for (D_mapint::iterator iter = map_remove_temp.begin(); iter != map_remove_temp.end(); ++iter)
	{
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel_1, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel_1, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
		extend_temp_x0 = extend_x0; extend_temp_x1 = extend_x1; extend_temp_y0 = extend_y0; extend_temp_y1 = extend_y1;
#if (C_DIMS==3)
		extend_temp_z0 = extend_z0; extend_temp_z1 = extend_z1;
#endif
		if ((x - extend_x0 * grid_ptr.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((x + extend_x1 * grid_ptr.dx) > xb_domain)
		{
			extend_temp_x1 = static_cast<unsigned int>((xb_domain - x) / grid_ptr.dx + C_eps);
		}
		if ((y - extend_y0 * grid_ptr.dx) < (C_dx * C_y0b_offset))
		{
			extend_temp_y0 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((y + extend_y1 * grid_ptr.dx) > yb_domain)
		{
			extend_temp_y1 = static_cast<unsigned int>((yb_domain - y) / grid_ptr.dx + C_eps);
		}
#if (C_DIMS == 3)
		if ((z - extend_z0 * grid_ptr.dx) < (C_dx * C_z0b_offset))
		{
			extend_temp_z0 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / grid_ptr.dx);
		}
		if ((z + extend_z1 * grid_ptr.dx) > zb_domain)
		{
			extend_temp_z1 = static_cast<unsigned int>((zb_domain - z) / grid_ptr.dx);
		}
#endif
#endif
#if (C_DIMS == 3)
		// (x, y) plane, negative z
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp1 = morton_temp << bit_move;
		for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid.find(morton_temp) != grid_ptr.grid.end())
				{
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		}

		// (x, y) plane, positive z
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) + extend_temp_z1 * 2);
		morton_temp1 = morton_temp << bit_move;
		for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid.find(morton_temp) != grid_ptr.grid.end())
				{
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		}

		// (x, z) plane, negative y
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid.find(morton_temp) != grid_ptr.grid.end())
				{
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (x, z) plane, positive y
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) + extend_temp_y1 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid.find(morton_temp) != grid_ptr.grid.end())
				{
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (y, z) plane, negative x
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int iy = 1; iy < extend_temp_y0 + extend_temp_y1; ++iy)
			{
				if (grid_ptr.grid.find(morton_temp) != grid_ptr.grid.end())
				{
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_y1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (y, z) plane, positve x
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) + extend_temp_x1 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int iy = 1; iy < extend_temp_y0 + extend_temp_y1; ++iy)
			{
				if (grid_ptr.grid.find(morton_temp) != grid_ptr.grid.end())
				{
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_y1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}
#endif
	}

	// check if nodes need to be removed at (ilevel - 1) and remove nodes
	extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
	extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif

	for (D_mapint::iterator iter = map_remove_surroundings.begin(); iter != map_remove_surroundings.end(); ++iter)
	{
		bool flag_need2remove = true;
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
		extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
		extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif
		if ((x - extend_x0 * grid_ptr.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x1 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((x + extend_x1 * grid_ptr.dx) > xb_domain)
		{
			extend_temp_x0 = static_cast<unsigned int>((xb_domain - x) / grid_ptr.dx);
		}
		if ((y - extend_y0 * grid_ptr.dx) < (C_dx * C_y0b_offset))
		{
			extend_temp_y1 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((y + extend_y1 * grid_ptr.dx) > yb_domain)
		{
			extend_temp_y0 = static_cast<unsigned int>((yb_domain - y) / grid_ptr.dx);
		}
#if (C_DIMS == 3)
		if ((z - extend_z0 * grid_ptr.dx) < (C_dx * C_z0b_offset))
		{
			extend_temp_z1 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((z + extend_z1 * grid_ptr.dx) > zb_domain)
		{
			extend_temp_z0 = static_cast<unsigned int>((zb_domain - z) / grid_ptr.dx);
		}
#endif
#endif
#if (C_DIMS == 2)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2);
		morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
#endif
#if (C_DIMS == 3)
		for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp1 = morton_temp2;
#endif
			for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
			{
				morton_temp = morton_temp1;
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (grid_ptr.grid.find(morton_temp) != grid_ptr.grid.end())
					{
						if (ilevel == C_max_level)
						{
							// when the node is within one grid space to the solid node
							if (grid_ptr.grid[morton_temp].flag == flag_near_solid)
							{
								flag_need2remove = false;
							}
						}
						else
						{
							if ((grid_ptr.coarse2fine.at(iboundary0_coarse).find(morton_temp) != grid_ptr.coarse2fine.at(iboundary0_coarse).end()) || (grid_ptr.coarse2fine.at(iboundary2_coarse).find(morton_temp) != grid_ptr.coarse2fine.at(iboundary2_coarse).end()))
							{
								flag_need2remove = false;
							}
						}
					}
					morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
				}
				morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
			}
#if (C_DIMS == 3)
			morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		}
#endif
		if (flag_need2remove && ((grid_ptr.grid[iter->first].flag == flag_iboundary[iboundary2]) || (grid_ptr.grid[iter->first].flag == flag_refine))) // delete nodes and its neighbours at (ilevel)
		{
			if (map_remove_out_temp.find(iter->first) == map_remove_out_temp.end())
			{
				map_remove_out_temp.insert({ iter->first , flag_remove });
			}
		}
	}
#endif

#if (C_SEARCH_METHOD == 2)
	Node_index node_index_temp;
	// compute icount of nodes on the boundary of a given region for adding
	for (D_mapint::iterator iter = map_add_temp.begin(); iter != map_add_temp.end(); ++iter)
	{
		int icount_temp = iter->second;
		node_index_temp.icount_refine = icount_temp;
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel_1, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel_1, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
		extend_temp_x0 = extend_x0; extend_temp_x1 = extend_x1; extend_temp_y0 = extend_y0; extend_temp_y1 = extend_y1;
#if (C_DIMS==3)
		extend_temp_z0 = extend_z0; extend_temp_z1 = extend_z1;
#endif
		if ((x - extend_x0 * grid_ptr.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((x + extend_x1 * grid_ptr.dx) > xb_domain)
		{
			extend_temp_x1 = static_cast<unsigned int>((xb_domain - x) / grid_ptr.dx + C_eps);
		}
		if ((y - extend_y0 * grid_ptr.dx) < (C_dx * C_y0b_offset))
		{
			extend_temp_y0 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((y + extend_y1 * grid_ptr.dx) > yb_domain)
		{
			extend_temp_y1 = static_cast<unsigned int>((yb_domain - y) / grid_ptr.dx + C_eps);
		}
#if (C_DIMS == 3)
		if ((z - extend_z0 * grid_ptr.dx) < (C_dx * C_z0b_offset))
		{
			extend_temp_z0 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((z + extend_z1 * grid_ptr.dx) > zb_domain)
		{
			extend_temp_z1 = static_cast<unsigned int>((zb_domain - z) / grid_ptr.dx + C_eps);
		}
#endif
#endif
#if (C_DIMS == 3)
		// (x, y) plane, negative z
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp1 = morton_temp << bit_move;
		for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if ((grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end()))
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				if ((grid_ptr.grid_index.find(morton_temp) == grid_ptr.grid_index.end()))
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
				}
				else
				{
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		}

		// (x, y) plane, positive z
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) + extend_temp_z1 * 2);
		morton_temp1 = morton_temp << bit_move;
		for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if ((grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end()))
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				if ((grid_ptr.grid_index.find(morton_temp) == grid_ptr.grid_index.end()))
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
				}
				else
				{
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		}

		// (x, z) plane, negative y
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if ((grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end()))
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				if ((grid_ptr.grid_index.find(morton_temp) == grid_ptr.grid_index.end()))
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
				}
				else
				{
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (x, z) plane, positive y
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) + extend_temp_y1 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if ((grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end()))
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				if ((grid_ptr.grid_index.find(morton_temp) == grid_ptr.grid_index.end()))
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
				}
				else
				{
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (y, z) plane, negative x
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int iy = 1; iy < extend_temp_y0 + extend_temp_y1; ++iy)
			{
				if ((grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end()))
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				if ((grid_ptr.grid_index.find(morton_temp) == grid_ptr.grid_index.end()))
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
				}
				else
				{
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				morton_temp = Morton_Assist::find_y1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (y, z) plane, positve x
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) + extend_temp_x1 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int iy = 1; iy < extend_temp_y0 + extend_temp_y1; ++iy)
			{
				if ((grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end()))
				{
					if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
					{
						map_add_surroundings.insert({ morton_temp, 1 });
					}
					//else
					//{
					//	++map_add_surroundings[morton_temp];
					//}
				}
				if ((grid_ptr.grid_index.find(morton_temp) == grid_ptr.grid_index.end()))
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
				}
				else
				{
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				morton_temp = Morton_Assist::find_y1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}
#endif
	}

	// compute icount of nodes on the boundary of a given region for removing
	for (D_mapint::iterator iter = map_remove_temp.begin(); iter != map_remove_temp.end(); ++iter)
	{
		int icount_temp = iter->second;
		node_index_temp.icount_refine = icount_temp;
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel_1, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel_1, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
		extend_temp_x0 = extend_x0; extend_temp_x1 = extend_x1; extend_temp_y0 = extend_y0; extend_temp_y1 = extend_y1;
#if (C_DIMS==3)
		extend_temp_z0 = extend_z0; extend_temp_z1 = extend_z1;
#endif
		if ((x - extend_x0 * grid_ptr.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((x + extend_x1 * grid_ptr.dx) > xb_domain)
		{
			extend_temp_x1 = static_cast<unsigned int>((xb_domain - x) / grid_ptr.dx + C_eps);
		}
		if ((y - extend_y0 * grid_ptr.dx) < (C_dx * C_y0b_offset))
		{
			extend_temp_y0 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((y + extend_y1 * grid_ptr.dx) > yb_domain)
		{
			extend_temp_y1 = static_cast<unsigned int>((yb_domain - y) / grid_ptr.dx + C_eps);
		}
#if (C_DIMS == 3)
		if ((z - extend_z0 * grid_ptr.dx) < (C_dx * C_z0b_offset))
		{
			extend_temp_z0 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((z + extend_z1 * grid_ptr.dx) > zb_domain)
		{
			extend_temp_z1 = static_cast<unsigned int>((zb_domain - z) / grid_ptr.dx + C_eps);
		}
#endif
#endif
#if (C_DIMS == 3)
		// (x, y) plane, negative z
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp1 = morton_temp << bit_move;
		for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid_index.find(morton_temp) != grid_ptr.grid_index.end())
				{
					grid_ptr.grid_index[morton_temp].icount_refine -= icount_temp;
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		}

		// (x, y) plane, positive z
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) + extend_temp_z1 * 2);
		morton_temp1 = morton_temp << bit_move;
		for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid_index.find(morton_temp) != grid_ptr.grid_index.end())
				{
					grid_ptr.grid_index[morton_temp].icount_refine -= icount_temp;
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		}

		// (x, z) plane, negative y
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid_index.find(morton_temp) != grid_ptr.grid_index.end())
				{
					grid_ptr.grid_index[morton_temp].icount_refine -= icount_temp;
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (x, z) plane, positive y
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) + extend_temp_y1 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
			{
				if (grid_ptr.grid_index.find(morton_temp) != grid_ptr.grid_index.end())
				{
					grid_ptr.grid_index[morton_temp].icount_refine -= icount_temp;
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (y, z) plane, negative x
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int iy = 1; iy < extend_temp_y0 + extend_temp_y1; ++iy)
			{
				if (grid_ptr.grid_index.find(morton_temp) != grid_ptr.grid_index.end())
				{
					grid_ptr.grid_index[morton_temp].icount_refine -= icount_temp;
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_y1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}

		// (y, z) plane, positve x
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) + extend_temp_x1 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
		morton_temp1 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
		for (unsigned int iz = 1; iz < extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp = morton_temp1;
			for (unsigned int iy = 1; iy < extend_temp_y0 + extend_temp_y1; ++iy)
			{
				if (grid_ptr.grid_index.find(morton_temp) != grid_ptr.grid_index.end())
				{
					grid_ptr.grid_index[morton_temp].icount_refine -= icount_temp;
					if (map_remove_surroundings.find(morton_temp) == map_remove_surroundings.end())
					{
						map_remove_surroundings.insert({ morton_temp, flag_refine });
					}
				}
				morton_temp = Morton_Assist::find_y1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}
#endif
	}

	D_mapint map_remove_surroundings_temp = {};
	for (D_mapint::iterator iter = map_remove_surroundings.begin(); iter != map_remove_surroundings.end(); ++iter)
	{
		
		if (grid_ptr.grid_index[iter->first].icount_refine <= 0)
		{
			if (grid_ptr.grid.find(iter->first) != grid_ptr.grid.end())
			{
				map_remove_surroundings_temp.insert({ iter->first, flag_remove });
			}		
			grid_ptr.grid_index.erase(iter->first);
		}
	}

	// check if nodes need to be removed at (ilevel - 1) and remove nodes
	extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
	extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif

	for (D_mapint::iterator iter = map_remove_surroundings_temp.begin(); iter != map_remove_surroundings_temp.end(); ++iter)
	{
		bool flag_need2remove = true;
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
		extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
		extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif
		if ((x - extend_x0 * grid_ptr.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x1 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((x + extend_x1 * grid_ptr.dx) > xb_domain)
		{
			extend_temp_x0 = static_cast<unsigned int>((xb_domain - x) / grid_ptr.dx + C_eps);
		}
		if ((y - extend_y0 * grid_ptr.dx) < (C_dx * C_y0b_offset))
		{
			extend_temp_y1 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((y + extend_y1 * grid_ptr.dx) > yb_domain)
		{
			extend_temp_y0 = static_cast<unsigned int>((yb_domain - y) / grid_ptr.dx + C_eps);
		}
#if (C_DIMS == 3)
		if ((z - extend_z0 * grid_ptr.dx) < (C_dx * C_z0b_offset))
		{
			extend_temp_z1 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((z + extend_z1 * grid_ptr.dx) > zb_domain)
		{
			extend_temp_z0 = static_cast<unsigned int>((zb_domain - z) / grid_ptr.dx + C_eps);
		}
#endif
#endif
#if (C_DIMS == 2)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2);
		morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
#endif
#if (C_DIMS == 3)
		for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp1 = morton_temp2;
#endif
			for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
			{
				morton_temp = morton_temp1;
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (grid_ptr.grid.find(morton_temp) != grid_ptr.grid.end())
					{
						if (ilevel == C_max_level)
						{
							// when the node is within one grid space to the solid node
							if (grid_ptr.grid[morton_temp].flag == flag_near_solid)
							{
								flag_need2remove = false;
							}
						}
						else
						{
							if ((grid_ptr.coarse2fine.at(iboundary0_coarse).find(morton_temp) != grid_ptr.coarse2fine.at(iboundary0_coarse).end()) || (grid_ptr.coarse2fine.at(iboundary2_coarse).find(morton_temp) != grid_ptr.coarse2fine.at(iboundary2_coarse).end()))
							{
								flag_need2remove = false;
							}
						}
					}
					morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
				}
				morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
			}
#if (C_DIMS == 3)
			morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		}
#endif
		if (flag_need2remove && ((grid_ptr.grid[iter->first].flag == flag_iboundary[iboundary2]) || (grid_ptr.grid[iter->first].flag == flag_refine))) // delete nodes and its neighbours at (ilevel)
		{
			if (map_remove_out_temp.find(iter->first) == map_remove_out_temp.end())
			{
				map_remove_out_temp.insert({ iter->first , flag_remove });
			}
		}
	}

	
#endif

#if(C_SEARCH_METHOD == 3)
	// compute icount of nodes within a given region
	Node_index node_index_temp;
	for (D_mapint::iterator iter = map_add_temp.begin(); iter != map_add_temp.end(); ++iter)
	{
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
		extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
		extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif
		if ((x - extend_x0 * grid_ptr.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((x + extend_x1 * grid_ptr.dx) > xb_domain)
		{
			extend_temp_x1 = static_cast<unsigned int>((xb_domain - x) / grid_ptr.dx + C_eps);
		}
		if ((y - extend_y0 * grid_ptr.dx) < (C_dx * C_y0b_offset))
		{
			extend_temp_y0 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((y + extend_y1 * grid_ptr.dx) > yb_domain)
		{
			extend_temp_y1 = static_cast<unsigned int>((yb_domain - y) / grid_ptr.dx + C_eps);
		}
#if (C_DIMS == 3)
		if ((z - extend_z0 * grid_ptr.dx) < (C_dx * C_z0b_offset))
		{
			extend_temp_z0 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((z + extend_z1 * grid_ptr.dx) > zb_domain)
		{
			extend_temp_z1 = static_cast<unsigned int>((zb_domain - z) / grid_ptr.dx + C_eps);
		}
#endif
#endif
#if (C_DIMS == 2)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2);
		morton_temp1 = morton_temp << bit_move;
#endif
#if (C_DIMS == 3)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
#endif
		int icount_temp = iter->second;
#if (C_DIMS == 3)
		for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp1 = morton_temp2;
#endif
			for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
			{
				morton_temp = morton_temp1;
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (grid_ptr.grid.find(morton_temp) == grid_ptr.grid.end())
					{
						if (map_add_surroundings.find(morton_temp) == map_add_surroundings.end())
						{
							map_add_surroundings.insert({ morton_temp, flag_refine });
						}
					}
					if (grid_ptr.grid_index.find(morton_temp) == grid_ptr.grid_index.end())
					{
						node_index_temp.icount_refine = icount_temp;
						grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
					}
					else
					{
						grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
					}
					morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
				}
				morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
			}
#if (C_DIMS == 3)
			morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		}
#endif		
	}

    for (D_mapint::iterator iter = map_remove_temp.begin(); iter != map_remove_temp.end(); ++iter)
	{
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
		extend_temp_x0 = extend_x1; extend_temp_x1 = extend_x0; extend_temp_y0 = extend_y1; extend_temp_y1 = extend_y0;
#if (C_DIMS==3)
		extend_temp_z0 = extend_z1; extend_temp_z1 = extend_z0;
#endif
		if ((x - extend_x0 * grid_ptr.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((x + extend_x1 * grid_ptr.dx) > xb_domain)
		{
			extend_temp_x1 = static_cast<unsigned int>((xb_domain - x) / grid_ptr.dx + C_eps);
		}
		if ((y - extend_y0 * grid_ptr.dx) < (C_dx * C_y0b_offset))
		{
			extend_temp_y0 = static_cast<unsigned int>((y - C_dx * C_y0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((y + extend_y1 * grid_ptr.dx) > yb_domain)
		{
			extend_temp_y1 = static_cast<unsigned int>((yb_domain - y) / grid_ptr.dx + C_eps);
		}
#if (C_DIMS == 3)
		if ((z - extend_z0 * grid_ptr.dx) < (C_dx * C_z0b_offset))
		{
			extend_temp_z0 = static_cast<unsigned int>((z - C_dx * C_z0b_offset) / grid_ptr.dx + C_eps);
		}
		if ((z + extend_z1 * grid_ptr.dx) > zb_domain)
		{
			extend_temp_z1 = static_cast<unsigned int>((zb_domain - z) / grid_ptr.dx + C_eps);
		}
#endif
#endif
#if (C_DIMS == 2)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2);
		morton_temp1 = morton_temp << bit_move;
#endif
#if (C_DIMS == 3)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / grid_ptr.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / grid_ptr.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / grid_ptr.dx + C_eps) - extend_temp_z0 * 2);
		morton_temp2 = morton_temp << bit_move;
#endif

		int icount_temp = iter->second;
#if (C_DIMS == 3)
		for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
		{
			morton_temp1 = morton_temp2;
#endif
			for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
			{
				morton_temp = morton_temp1;
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if ((grid_ptr.grid_index.find(morton_temp) != grid_ptr.grid_index.end()))
					{
						grid_ptr.grid_index[morton_temp].icount_refine -= icount_temp;
						if (grid_ptr.grid_index[morton_temp].icount_refine <= 0)
						{

							if (map_remove_out_temp.find(morton_temp) == map_remove_out_temp.end())
							{
								map_remove_out_temp.insert({ morton_temp , flag_remove });
							}
							grid_ptr.grid_index.erase(morton_temp);
						}
					}
					morton_temp = Morton_Assist::find_x1(morton_temp, ilevel_1);
				}
				morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel_1);
			}
#if (C_DIMS == 3)
			morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel_1);
		}
#endif		
	}

#endif
	
//	double t1 = tmr.elapsed();
//std::cout << "time for adding: " << t1 - t0 << std::endl;


#if (C_SEARCH_METHOD == 1)
	// reconstruct numerical boundaries due to removing nodes
    reconstruct_numerical_boundary_removing(ilevel, grid_ptr, node_temp, map_remove_out_temp, map_add_nodes_out, map_remove_nodes_out);

	// find if nodes need to be added
	for (D_mapint::iterator iter = map_add_surroundings.begin(); iter != map_add_surroundings.end(); ++iter)
	{
		bool flag_need2add = false;
		// (-x) 
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
			{
				flag_need2add = true;
			}
		}

		// (-x, -y) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (-x, -y, -z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-x, -y, +z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-x, -z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}
#endif

		// (-x, +y)
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (-x, +y, -z)
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-x, +y, +z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-x, +z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#endif
		// (+x) 
		if (!flag_need2add)
		{
			morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+x, -y) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (+x, -y, -z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+x, -y, +z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+x, -z)
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}


#endif
		// (+x, +y) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (+x, +y, -z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+x, +y, +z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+x, +z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#endif
		// (-y) 
		if (!flag_need2add)
		{
			morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (-y, -z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-y, +z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#endif
		// (+y) 
		if (!flag_need2add)
		{
			morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (+y, -z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+y, +z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-z) 
		if (!flag_need2add)
		{
			morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+z) 
		if (!flag_need2add)
		{
			morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#endif
		if (flag_need2add)
		{
			if (map_add_out_temp.find(iter->first) == map_add_out_temp.end())
			{
				map_add_out_temp.insert({ iter->first , iter->second });
			}
		}

	}


	// reconstruct numerical boundaries due to adding nodes
	reconstruct_numerical_boundary_adding(ilevel, grid_ptr, node_temp, map_add_out_temp, map_add_nodes_out, map_remove_nodes_out);

#else
   // find if nodes need to be added
for (D_mapint::iterator iter = map_add_surroundings.begin(); iter != map_add_surroundings.end(); ++iter)
{
	if (grid_ptr.grid_index.find(iter->first) != grid_ptr.grid_index.end())
	{
		bool flag_need2add = false;
		// (-x) 
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
			{
				flag_need2add = true;
			}
		}

		// (-x, -y) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}


#if(C_DIMS==3)
		// (-x, -y, -z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-x, -y, +z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-x, -z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}
#endif

		// (-x, +y)
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (-x, +y, -z)
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-x, +y, +z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-x, +z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#endif
		// (+x) 
		if (!flag_need2add)
		{
			morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+x, -y) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (+x, -y, -z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+x, -y, +z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+x, -z)
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}


#endif
		// (+x, +y) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (+x, +y, -z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+x, +y, +z) 
		if (!flag_need2add)
		{
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+x, +z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#endif
		// (-y) 
		if (!flag_need2add)
		{
			morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (-y, -z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-y, +z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#endif
		// (+y) 
		if (!flag_need2add)
		{
			morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#if(C_DIMS==3)
		// (+y, -z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+y, +z) 
		if (!flag_need2add)
		{
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (-z) 
		if (!flag_need2add)
		{
			morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

		// (+z) 
		if (!flag_need2add)
		{
			morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					flag_need2add = true;
				}
			}
		}

#endif
		if (flag_need2add)
		{
			if (map_add_out_temp.find(iter->first) == map_add_out_temp.end())
			{
				map_add_out_temp.insert({ iter->first , iter->second });
			}
		}

	}
}
	// reconstruct numerical boundaries due to adding nodes
	reconstruct_numerical_boundary_adding(ilevel, grid_ptr, node_temp, map_add_out_temp, map_add_nodes_out, map_remove_nodes_out);
	// reconstruct numerical boundaries due to removing nodes
	reconstruct_numerical_boundary_removing(ilevel, grid_ptr, node_temp, map_remove_out_temp, map_add_nodes_out, map_remove_nodes_out);
#endif


}

/**
* @brief function to update boundary information due to removing nodes.
* @param[in]  ilevel                refinement level.
* @param[in]  grid_ptr              pointer point to the class storing grid information at ilevel.
* @param[in]  map_add_nodes_in      nodes used to search if there are nodes needed to be added at ilevel, i.e. nodes at (ilevel + 1) for (grid_NIB) and nodes near solid for (grid_IB)
* @param[in]  map_remove_nodes_in   nodes used to search if there are nodes needed to be removed at ilevel, i.e. nodes at (ilevel + 1) for (grid_NIB) and nodes near solid for (grid_IB)
* @param[out]  map_add_nodes_out     nodes used to search if there are nodes needed to be added at (ilevel - 1)
* @param[out]  map_remove_nodes_out  nodes used to search if there are nodes needed to be added at (ilevel - 1)
* @todo {to avoid search backward to check if nodes need to be removed after searching surrounding boundary, count number can be stored at the boundary. the problem is how to handle block splitting when the removed nodes is not on the numerical boundary}
*/
template <class T_grid, class T_node>
void Grid_Manager::reconstruct_numerical_boundary_removing(const unsigned int ilevel, T_grid& grid_ptr, T_node node_temp, D_mapint& map_remove_temp, D_mapint& map_add_nodes_out, D_mapint& map_remove_nodes_out)
{
	D_morton morton_temp0, morton_temp1, morton_temp0_ilevel_1, morton_temp1_ilevel_1;
	unsigned int ilevel_1 = ilevel - 1;
	int flag_remove = -1;
	Grid_NIB* grid_ptr_coarse = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel_1));
#if(C_DIMS==3)
	D_morton  morton_temp2, morton_temp2_ilevel_1;
#endif
	D_mapint remove_temp = {};
	D_mapint map_remove_out_temp = {};

	// delete nodes in (map_remove_temp) and within one grid space at (ilevel) in map (grid_ptr.grid, grid_ptr.fine2coarse, grid_ptr_coarse->coarse2fine, and bk_boundary_?)
	for (D_mapint::iterator iter = map_remove_temp.begin(); iter != map_remove_temp.end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) != grid_ptr.grid.end())
		{
			if ((grid_ptr.grid[iter->first].flag == flag_iboundary[iboundary2]) || (grid_ptr.grid[iter->first].flag == flag_refine))
			{
				if (map_remove_out_temp.find(iter->first) == map_remove_out_temp.end())
				{
					map_remove_out_temp.insert({ iter->first , flag_remove });
				}
			}
			//(center)
			if (grid_ptr.grid[iter->first].flag == flag_iboundary[iboundary2])
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(iter->first);
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).erase(iter->first);
				if (map_remove_nodes_out.find(iter->first) == map_remove_nodes_out.end())
				{
					map_remove_nodes_out.insert({ iter->first , flag_refine });
				}
			}

			grid_ptr.grid.erase(iter->first);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(iter->first) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(iter->first);
				}
			}
			else if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(iter->first) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(iter->first);
				}
			}
			if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(iter->first) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(iter->first);
				}
			}
			else if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(iter->first) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(iter->first);
				}
			}
#if(C_DIMS == 3)
			if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(iter->first) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(iter->first);
				}
			}
			else if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(iter->first) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(iter->first);
				}
			}
#endif
#endif
			// (-x) 
			morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
					}
				}
				else if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
				}
				grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
					}
				}
				if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
					}
				}
#endif
#endif
			}
			// (-x, -y) 
			morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
#endif
			}
#if(C_DIMS==3)
			// (-x, -y, -z) 
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
					}
				}
				else if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
				}
				grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
					}
				}
#endif
			}
			// (-x, -y, +z) 
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
					}
				}
				else if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
				}
				grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
					}
				}
#endif
			}
			// (-x, -z) 
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
#endif
			}
#endif
			// (-x, +y) 
			morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
#endif
			}
#if(C_DIMS==3)
			// (-x, +y, -z) 
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
					}
				}
				else if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
				}
				grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
					}
				}
#endif
			}
			// (-x, +y, +z) 
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
					}
				}
				else if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
				}
				grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
					}
				}
#endif
			}
			// (-x, +z) 
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
#endif
			}
#endif
			// (+x) 
			morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
					}
				}
				else if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
				}
				grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
					}
				}
				if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
					}
				}
#endif
#endif
			}
			// (+x, -y) 
			morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
#endif
			}
#if(C_DIMS==3)
			// (+x, -y, -z) 
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
					}
				}
				else if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
				}
				grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
					}
				}
#endif
			}
			// (+x, -y, +z) 
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
					}
				}
				else if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
				}
				grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
					}
				}
#endif
			}
			// (+x, -z) 
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
#endif
			}
#endif
			// (+x, +y) 
			morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
#endif
			}
#if(C_DIMS==3)
			// (+x, +y, -z) 
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
					}
				}
				else if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
				}
				grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
					}
				}
#endif
			}
			// (+x, +y, +z) 
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
					}
				}
				else if (grid_ptr.grid[morton_temp2].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
				}
				grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
					}
				}
				if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
					}
				}
				else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
					}
				}
#endif
			}
			// (+x, +z) 
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
			}
#endif
			// (-y) 
			morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
					}
				}
				else if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
				}
				grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
					}
				}
				if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
					}
				}
#endif
#endif
			}
#if(C_DIMS==3)
			// (-y, -z) 
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
			}
			// (-y, +z) 
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
			}
#endif
			// (+y) 
			morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
					}
				}
				else if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
				}
				grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
					}
				}
				if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
					}
				}
#if(C_DIMS == 3)
				if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
					}
				}
#endif
#endif
			}
#if(C_DIMS==3)
			// (+y, -z) 
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}
				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
			}
			// (+y, +z) 
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
					}
				}
				else if (grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
				}

				grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
					}
				}
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
					}
				}
#endif
			}
			// (-z) 
			morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
					}
				}
				else if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
				}
				grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
					}
				}
				if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
					}
				}
				if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
					}
				}
#endif
			}
			// (+z) 
			morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
					}
				}
				else if (grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary1])
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
				}
				grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
				if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
					}
				}
				if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
					}
				}
				if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
					}
				}
				else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
					}
				}
#endif
			}
#endif
		}
	}

	Node node_temp_ilevel_1;
	node_temp_ilevel_1.flag = flag_refine;
#if (C_CHECK_MORTON_BOUNDARY == 1)
	bool bool_bx0, bool_bx1, bool_by0, bool_by1;
#if (C_DIMS==3)
	bool bool_bz0, bool_bz1;
#endif
#endif

	// find nodes near (map_remove_temp) and which will be new iboundary2 and iboudnary2_coarse
	for (D_mapint::iterator iter = map_remove_out_temp.begin(); iter != map_remove_out_temp.end(); ++iter)
	{
#if (C_CHECK_MORTON_BOUNDARY == 1)
		bool_bx0 = true;  bool_bx1 = true;  bool_by0 = true;  bool_by1 = true;
#if (C_DIMS==3)
		bool_bz0 = true; bool_bz1 = true;
#endif
		if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
		{
			bool_bx0 = false;
		}
		else if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
		{
			bool_bx1 = false;
		}
		if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
		{
			bool_by0 = false;
		}
		else if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
		{
			bool_by1 = false;
		}
#if (C_DIMS==3)
		if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
		{
			bool_bz0 = false;
		}
		else if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
		{
			bool_bz1 = false;
		}
#endif
#endif
		// (-x) 
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp0) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp0) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp0, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp0, flag_refine });
				if (map_add_nodes_out.find(morton_temp0) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp0 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp0);
			}
			if (remove_temp.find(morton_temp0) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp0,flag_remove });
			}
		}

		// (-x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}


#if(C_DIMS==3)
		// (-x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp2, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp2) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp2) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp2, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2, &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp2, flag_refine });
				if (map_add_nodes_out.find(morton_temp2) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp2 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp2);
			}
			if (remove_temp.find(morton_temp2) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp2,flag_remove });
			}
		}

		// (-x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp2, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp2) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp2) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp2, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2, &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp2, flag_refine });
				if (map_add_nodes_out.find(morton_temp2) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp2 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp2);
			}
			if (remove_temp.find(morton_temp2) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp2,flag_remove });
			}
		}

		// (-x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

#endif

		// (-x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

#if(C_DIMS==3)
		// (-x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp2, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp2) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp2) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp2, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2, &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp2, flag_refine });
				if (map_add_nodes_out.find(morton_temp2) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp2 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp2);
			}
			if (remove_temp.find(morton_temp2) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp2,flag_remove });
			}
		}

		// (-x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp2, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp2) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp2) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp2, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2, &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp2, flag_refine });
				if (map_add_nodes_out.find(morton_temp2) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp2 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp2);
			}
			if (remove_temp.find(morton_temp2) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp2,flag_remove });
			}
		}

		// (-x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

#endif
		// (+x) 
		morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp0) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp0) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp0, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp0, flag_refine });
				if (map_add_nodes_out.find(morton_temp0) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp0 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp0);
			}
			if (remove_temp.find(morton_temp0) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp0,flag_remove });
			}
		}

		// (+x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

#if(C_DIMS==3)
		// (+x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp2, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp2) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp2) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp2, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2, &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp2, flag_refine });
				if (map_add_nodes_out.find(morton_temp2) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp2 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp2);
			}
			if (remove_temp.find(morton_temp2) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp2,flag_remove });
			}
		}

		// (+x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp2, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp2) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp2) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp2, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2, &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp2, flag_refine });
				if (map_add_nodes_out.find(morton_temp2) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp2 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp2);
			}
			if (remove_temp.find(morton_temp2) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp2,flag_remove });
			}
		}

		// (+x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

#endif
		// (+x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

#if(C_DIMS==3)
		// (+x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp2, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp2) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp2) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp2, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2, &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp2, flag_refine });
				if (map_add_nodes_out.find(morton_temp2) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp2 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp2);
			}
			if (remove_temp.find(morton_temp2) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp2,flag_remove });
			}
		}

		// (+x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp2, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp2) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp2) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp2, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2, &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2 , &grid_ptr_coarse->grid[morton_temp2] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp2, flag_refine });
				if (map_add_nodes_out.find(morton_temp2) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp2 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp2);
			}
			if (remove_temp.find(morton_temp2) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp2,flag_remove });
			}
		}

		// (+x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

#endif
		// (-y) 
		morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp0) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp0) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp0, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp0, flag_refine });
				if (map_add_nodes_out.find(morton_temp0) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp0 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp0);
			}
			if (remove_temp.find(morton_temp0) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp0,flag_remove });
			}
		}

#if(C_DIMS==3)
		// (-y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

		// (-y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

#endif
		// (+y) 
		morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp0) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp0) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp0, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp0, flag_refine });
				if (map_add_nodes_out.find(morton_temp0) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp0 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp0);
			}
			if (remove_temp.find(morton_temp0) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp0,flag_remove });
			}
		}

#if(C_DIMS==3)
		// (+y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

		// (+y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp1, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp1) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp1) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp1, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1, &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1 , &grid_ptr_coarse->grid[morton_temp1] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp1, flag_refine });
				if (map_add_nodes_out.find(morton_temp1) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp1 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp1);
			}
			if (remove_temp.find(morton_temp1) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp1,flag_remove });
			}
		}

		// (-z) 
		morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp0) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp0) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp0, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp0, flag_refine });
				if (map_add_nodes_out.find(morton_temp0) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp0 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp0);
			}
			if (remove_temp.find(morton_temp0) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp0,flag_remove });
			}
		}

		// (+z) 
		morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(morton_temp0) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				if (grid_ptr_coarse->grid.find(morton_temp0) == grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.insert({ morton_temp0, node_temp_ilevel_1 });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0) == bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0 , &grid_ptr_coarse->grid[morton_temp0] });
						}
					}
#endif
				}
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ morton_temp0, flag_refine });
				if (map_add_nodes_out.find(morton_temp0) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ morton_temp0 , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(morton_temp0);
			}
			if (remove_temp.find(morton_temp0) == remove_temp.end())
			{
				remove_temp.insert({ morton_temp0,flag_remove });
			}
		}

#endif
	}

	bool bool_xy;
	bool bool_x_temp0, bool_y_temp0, bool_x_temp1, bool_y_temp1;
#if(C_DIMS == 3)
	bool  bool_yz, bool_xz;
	bool bool_z_temp0, bool_z_temp1;
#endif

#if(C_DIMS == 3)
	D_mapint iboundary2_temp = {};
#endif
	D_mapint map_cell_delete = {};
	// update the boundary information based on nodes at (ilevel) on iboundary2
	for (D_mapint::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
	{
		// (-x)
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_x0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
						grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
			}
			else
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp0_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp0_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp0_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#if (C_DIMS == 3)
						if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#endif
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp0_ilevel_1, flag_remove });
				}
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_remove });
					}
				}
			}
		}

		// (-x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#if (C_DIMS == 3)
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
#if (C_DIMS == 2)
				if (map_cell_delete.find(morton_temp1) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp1, flag_remove });
				}
#endif
#if (C_DIMS == 3)
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
#endif
			}

		}

#if(C_DIMS==3)
		// (-x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_remove });
				}
			}
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp2_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp2_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp2_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1, &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp2_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_remove });
				}
			}
		}

		// (-x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_remove });
				}
			}
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp2_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp2_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp2_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1, &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp2_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_remove });
				}
			}
		}


		// (-x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
			}
		}
#endif

		// (-x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#if (C_DIMS == 3)
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
#if (C_DIMS == 2)
				if (map_cell_delete.find(morton_temp1) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp1, flag_remove });
				}
#endif
#if (C_DIMS == 3)
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
#endif
			}
		}

#if(C_DIMS==3)
		// (-x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_remove });
				}
			}
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp2_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp2_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp2_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1, &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp2_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_remove });
				}
			}
		}

		// (-x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_remove });
				}
			}
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp2_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp2_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp2_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1, &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp2_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_remove });
				}
			}
		}

		// (-x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
			}
		}

#endif
		// (+x) 
		morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_x1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
						grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
			}
			else
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp0_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp0_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp0_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#if (C_DIMS == 3)
						if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#endif
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp0_ilevel_1, flag_remove });
				}
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_remove });
					}
				}
			}
		}

		// (+x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#if (C_DIMS == 3)
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
#if (C_DIMS == 2)
				if (map_cell_delete.find(morton_temp1) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp1, flag_remove });
				}
#endif
#if (C_DIMS == 3)
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
#endif
			}
		}

#if(C_DIMS==3)
		// (+x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_remove });
				}
			}
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp2_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp2_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp2_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1, &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp2_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_remove });
				}
			}
		}

		// (+x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_remove });
				}
			}
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp2_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp2_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp2_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1, &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp2_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_remove });
				}
			}
		}

		// (+x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
			}
		}

#endif
		// (+x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#if (C_DIMS == 3)
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
#if (C_DIMS == 2)
				if (map_cell_delete.find(morton_temp1) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp1, flag_remove });
				}
#endif
#if (C_DIMS == 3)
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
#endif
			}
		}

#if(C_DIMS==3)
		// (+x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_remove });
				}
			}
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp2_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp2_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp2_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1, &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp2_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_remove });
				}
			}
		}

		// (+x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_remove });
				}
			}
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp2_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp2_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp2_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp2_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1, &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
						else if ((morton_temp2_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp2_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp2_ilevel_1 , &grid_ptr_coarse->grid[morton_temp2_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp2_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_remove });
				}
			}
		}

		// (+x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
			}
		}

#endif
		// (-y) 
		morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_y0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
						grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
			}
			else
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp0_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp0_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp0_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#if (C_DIMS == 3)
						if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#endif
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp0_ilevel_1, flag_remove });
				}
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_remove });
					}
				}
			}
		}

#if(C_DIMS==3)
		// (-y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
			}
		}

		// (-y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
			}
		}

#endif
		// (+y) 
		morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_y1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{

			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
						grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
			}
			else
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp0_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp0_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp0_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#if (C_DIMS == 3)
						if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#endif
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp0_ilevel_1, flag_remove });
				}
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_remove });
					}
				}
			}
		}

#if(C_DIMS==3)
		// (+y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
			}
		}

		// (+y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag != flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_remove });
					}
				}
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp1_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp1_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp1_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp1_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1, &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
						else if ((morton_temp1_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp1_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp1_ilevel_1 , &grid_ptr_coarse->grid[morton_temp1_ilevel_1] });
							}
						}
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp1_ilevel_1, flag_remove });
				}
			}
			else
			{
				if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
				{
					if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
					{
						iboundary2_temp.insert({ morton_temp1, flag_remove });
					}
				}
			}
		}

		// (-z) 
		morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_z0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
						grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
			}
			else
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp0_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp0_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp0_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#if (C_DIMS == 3)
						if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#endif
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp0_ilevel_1, flag_remove });
				}
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_remove });
					}
				}
			}
		}

		// (+z) 
		morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_z1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
						grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_remove });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
			}
			else
			{
				if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(morton_temp0_ilevel_1) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
				{
					if (grid_ptr_coarse->grid.find(morton_temp0_ilevel_1) == grid_ptr_coarse->grid.end())
					{
						grid_ptr_coarse->grid.insert({ morton_temp0_ilevel_1 , node_temp_ilevel_1 });
						//grid_ptr_coarse->grid[morton_temp0_ilevel_1].flag = flag_inner_boundary;
#if (C_CHECK_MORTON_BOUNDARY == 1)
						if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
						{
							if (bk_boundary_x.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[0].end())
							{
								bk_boundary_x.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
						{
							if (bk_boundary_x.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_x.at(ilevel_1)[1].end())
							{
								bk_boundary_x.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
						{
							if (bk_boundary_y.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[0].end())
							{
								bk_boundary_y.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
						{
							if (bk_boundary_y.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_y.at(ilevel_1)[1].end())
							{
								bk_boundary_y.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#if (C_DIMS == 3)
						if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
						{
							if (bk_boundary_z.at(ilevel_1)[0].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[0].end())
							{
								bk_boundary_z.at(ilevel_1)[0].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
						else if ((morton_temp0_ilevel_1 & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
						{
							if (bk_boundary_z.at(ilevel_1)[1].find(morton_temp0_ilevel_1) == bk_boundary_z.at(ilevel_1)[1].end())
							{
								bk_boundary_z.at(ilevel_1)[1].insert({ morton_temp0_ilevel_1 , &grid_ptr_coarse->grid[morton_temp0_ilevel_1] });
							}
						}
#endif
#endif
					}
					grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ morton_temp0_ilevel_1, flag_remove });
				}
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_remove });
					}
				}
			}
		}

#endif
	}

	// delete isolate cell, i.e. all its faces (3D) or edges (2D) are on iboundary2 
	D_morton morton_ref0 = ~Morton_Assist::morton_xyz.at(ilevel);
	D_morton morton_corner;
	D_mapint map_iboudnary2_delete_temp = {};
	//	@todo use flag,i.e. x0x1y0y1z0z1 in previous for loop, to identify if flag of nodes is flag_iboundary[2]
	for (D_mapint::iterator iter = map_cell_delete.begin(); iter != map_cell_delete.end(); ++iter)
	{
		bool bool_need2delete = true;
		// (x0, y0, z0)
		morton_corner = iter->first & morton_ref0;
		if (grid_ptr.grid.find(morton_corner) == grid_ptr.grid.end() || grid_ptr.grid[morton_corner].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

		// (x1, y0, z0)
		morton_temp0_ilevel_1 = Morton_Assist::find_x1(morton_corner, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

		// (x0, y1, z0)
		morton_temp0_ilevel_1 = Morton_Assist::find_y1(morton_corner, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

		// (x1, y1, z0)
		morton_temp1_ilevel_1 = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

#if (C_DIMS==3)
		// (x0, y0, z1)
		morton_temp0_ilevel_1 = Morton_Assist::find_z1(morton_corner, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}
		// (x1, y0, z1)
		morton_temp1_ilevel_1 = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

		// (x0, y1, z1)
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

		// (x1, y1, z1)
		morton_temp2_ilevel_1 = Morton_Assist::find_x1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}
#endif

		if (bool_need2delete)
		{
			D_morton morton_temp0_inner, morton_temp1_inner;
#if(C_DIMS==3)
			D_morton morton_temp2_inner;
#endif
			// (x0, y0, z0)
			if (grid_ptr.grid.find(morton_corner) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_corner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_corner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_corner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_corner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
#if (C_DIMS == 3)
				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_corner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_corner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_corner) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_corner, flag_remove });
					}
				}
			}

			// (x1, y0, z0)
			morton_temp0_ilevel_1 = Morton_Assist::find_x1(morton_corner, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp0_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
#if (C_DIMS == 3)
				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp0_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp0_ilevel_1, flag_remove });
					}
				}
			}


			// (x0, y1, z0)
			morton_temp0_ilevel_1 = Morton_Assist::find_y1(morton_corner, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp0_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
#if (C_DIMS == 3)
				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp0_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp0_ilevel_1, flag_remove });
					}
				}
			}

			// (x1, y1, z0)
			morton_temp1_ilevel_1 = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp1_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
#if (C_DIMS == 3)
				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp1_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp1_ilevel_1, flag_remove });
					}
				}
			}
#if (C_DIMS==3)
			// (x0, y0, z1)
			morton_temp0_ilevel_1 = Morton_Assist::find_z1(morton_corner, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp0_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp0_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp0_ilevel_1, flag_remove });
					}
				}
			}

			// (x1, y0, z1)
			morton_temp1_ilevel_1 = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp1_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp1_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp1_ilevel_1, flag_remove });
					}
				}
			}


			// (x0, y1, z1)
			morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp1_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}


				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp1_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp1_ilevel_1, flag_remove });
					}
				}
			}


			// (x1, y1, z1)
			morton_temp2_ilevel_1 = Morton_Assist::find_x1(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp2_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp2_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp2_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp2_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp2_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp2_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp2_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp2_ilevel_1, flag_remove });
					}
				}
			}
#endif

		}
	}

	// delete boundary information based on nodes got in previous for loop
	for (D_mapint::iterator iter = map_iboudnary2_delete_temp.begin(); iter != map_iboudnary2_delete_temp.end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(iter->first);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(iter->first) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(iter->first);
				}
			}
			else if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(iter->first) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(iter->first);
				}
			}
			if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(iter->first) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(iter->first);
				}
			}
			else if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(iter->first) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(iter->first);
				}
			}
#if(C_DIMS == 3)
			if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(iter->first) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(iter->first);
				}
			}
			else if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(iter->first) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(iter->first);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(iter->first) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(iter->first);
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(iter->first) != grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).erase(iter->first);
				if (map_remove_nodes_out.find(iter->first) == map_remove_nodes_out.end())
				{
					map_remove_nodes_out.insert({ iter->first , flag_refine });
				}
			}
		}

		// (-x) 
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

		// (-x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#if(C_DIMS==3)
		// (-x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (-x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (-x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}
#endif

		// (-x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#if(C_DIMS==3)
		// (-x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (-x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (-x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#endif
		// (+x) 
		morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

		// (+x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#if(C_DIMS==3)
		// (+x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (+x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (+x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}
#endif

		// (+x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#if(C_DIMS==3)
		// (+x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (+x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (+x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#endif
		// (-y) 
		morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

#if(C_DIMS==3)
		// (-y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

		// (-y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#endif
		// (+y) 
		morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

#if(C_DIMS==3)
		// (+y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

		// (+y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

		// (-z) 
		morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

		// (+z) 
		morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

#endif

	}

#if (C_DIMS == 3)
	// update iboundary2 information on (xy, xz, or yz) planes
	for (D_mapint::iterator iter = iboundary2_temp.begin(); iter != iboundary2_temp.end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[iter->first].flag != flag_iboundary[iboundary1])
			{
				if (grid_ptr.fine2coarse.at(iboundary2).find(iter->first) == grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.grid[iter->first].flag = flag_iboundary[iboundary2];
					grid_ptr.fine2coarse.at(iboundary2).insert({ iter->first, flag_remove });
				}
			}
			else
			{
				morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
				bool_x_temp0 = grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2];
				if (bool_x_temp0)
				{
					morton_temp1 = Morton_Assist::find_x1(iter->first, ilevel);
					bool_x_temp1 = grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2];
				}

				morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
				bool_y_temp0 = grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2];
				if (bool_y_temp0)
				{
					morton_temp1 = Morton_Assist::find_y1(iter->first, ilevel);
					bool_y_temp1 = grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2];
				}

				morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
				bool_z_temp0 = grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2];
				if (bool_z_temp0)
				{
					morton_temp1 = Morton_Assist::find_z1(iter->first, ilevel);
					bool_z_temp1 = grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2];
				}
				// check if node is on any of the plane which is on iboudanry2
				bool_xy = bool_x_temp0 && bool_x_temp1 && bool_y_temp0 && bool_y_temp1;
				bool_xz = bool_x_temp0 && bool_x_temp1 && bool_z_temp0 && bool_z_temp1;
				bool_yz = bool_y_temp0 && bool_y_temp1 && bool_z_temp0 && bool_z_temp1;
				if (bool_xy || bool_xz || bool_yz)
				{
					grid_ptr.grid[iter->first].flag = flag_iboundary[iboundary2];
					grid_ptr.fine2coarse.at(iboundary2).insert({ iter->first, flag_remove });
					if (grid_ptr.fine2coarse.at(iboundary1).find(iter->first) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(iter->first);
					}
				}
			}
		}
	}
#endif

	// remove keys stored for domain boundary (bk_boundary)
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 2))
	std::vector<D_morton> boundary_temp;
	for (D_mapNodePtr::iterator iter = bk_boundary_x.at(ilevel).at(0).begin(); iter != bk_boundary_x.at(ilevel).at(0).end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) == grid_ptr.grid.end())
		{
			boundary_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = boundary_temp.begin(); iter != boundary_temp.end(); ++iter)
	{
		bk_boundary_x.at(ilevel).at(0).erase(*iter);
	}
	boundary_temp.clear();
	for (D_mapNodePtr::iterator iter = bk_boundary_x.at(ilevel).at(1).begin(); iter != bk_boundary_x.at(ilevel).at(1).end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) == grid_ptr.grid.end())
		{
			boundary_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = boundary_temp.begin(); iter != boundary_temp.end(); ++iter)
	{
		bk_boundary_x.at(ilevel).at(1).erase(*iter);
	}
	boundary_temp.clear();
	for (D_mapNodePtr::iterator iter = bk_boundary_y.at(ilevel).at(0).begin(); iter != bk_boundary_y.at(ilevel).at(0).end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) == grid_ptr.grid.end())
		{
			boundary_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = boundary_temp.begin(); iter != boundary_temp.end(); ++iter)
	{
		bk_boundary_y.at(ilevel).at(0).erase(*iter);
	}
	boundary_temp.clear();
	for (D_mapNodePtr::iterator iter = bk_boundary_y.at(ilevel).at(1).begin(); iter != bk_boundary_y.at(ilevel).at(1).end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) == grid_ptr.grid.end())
		{
			boundary_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = boundary_temp.begin(); iter != boundary_temp.end(); ++iter)
	{
		bk_boundary_y.at(ilevel).at(1).erase(*iter);
	}
	boundary_temp.clear();
#if (C_DIMS == 3)
	for (D_mapNodePtr::iterator iter = bk_boundary_z.at(ilevel).at(0).begin(); iter != bk_boundary_z.at(ilevel).at(0).end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) == grid_ptr.grid.end())
		{
			boundary_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = boundary_temp.begin(); iter != boundary_temp.end(); ++iter)
	{
		bk_boundary_z.at(ilevel).at(0).erase(*iter);
	}
	boundary_temp.clear();
	for (D_mapNodePtr::iterator iter = bk_boundary_z.at(ilevel).at(1).begin(); iter != bk_boundary_z.at(ilevel).at(1).end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) == grid_ptr.grid.end())
		{
			boundary_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = boundary_temp.begin(); iter != boundary_temp.end(); ++iter)
	{
		bk_boundary_z.at(ilevel).at(1).erase(*iter);
	}
	boundary_temp.clear();
#endif
#endif
}

/**
* @brief function to add nodes and update boundary information due to adding nodes.
* @param[in]  ilevel                refinement level.
* @param[in]  grid_ptr              pointer point to the class storing grid information at ilevel.
* @param[in]  node_temp             template for node at ilevel.
* @param[in]  map_add_temp          nodes need to be added
* @param[out]  map_add_nodes_out     nodes used to search if there are nodes needed to be added at (ilevel - 1)
* @param[out]  map_remove_nodes_out  nodes used to search if there are nodes needed to be added at (ilevel - 1)
*/
template <class T_grid, class T_node>
void Grid_Manager::reconstruct_numerical_boundary_adding(const unsigned int ilevel, T_grid& grid_ptr, T_node node_temp, D_mapint& map_add_temp, D_mapint& map_add_nodes_out, D_mapint& map_remove_nodes_out)
{
	unsigned int ilevel_1 = ilevel - 1;
	Grid_NIB* grid_ptr_coarse = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel - 1));

	D_morton morton_temp0_ilevel_1, morton_temp1_ilevel_1;
	D_morton morton_temp0, morton_temp1;
#if(C_DIMS == 3)
	D_morton morton_temp2_ilevel_1;
	D_morton morton_temp2;
#endif

	// add nodes to grid
	int flag_add = flag_refine;
	for (D_mapint::iterator iter = map_add_temp.begin(); iter != map_add_temp.end(); ++iter)
	{
		node_temp.flag = iter->second;
		grid_ptr.grid.insert(make_pair(iter->first, node_temp));
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
		{
			if (bk_boundary_x.at(ilevel)[0].find(iter->first) == bk_boundary_x.at(ilevel)[0].end())
			{
				bk_boundary_x.at(ilevel)[0].insert({ iter->first, &grid_ptr.grid[iter->first] });
			}
		}
		else if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
		{
			if (bk_boundary_x.at(ilevel)[1].find(iter->first) == bk_boundary_x.at(ilevel)[1].end())
			{
				bk_boundary_x.at(ilevel)[1].insert({ iter->first , &grid_ptr.grid[iter->first] });
			}
		}
		if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
		{
			if (bk_boundary_y.at(ilevel)[0].find(iter->first) == bk_boundary_y.at(ilevel)[0].end())
			{
				bk_boundary_y.at(ilevel)[0].insert({ iter->first , &grid_ptr.grid[iter->first] });
			}
		}
		else if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
		{
			if (bk_boundary_y.at(ilevel)[1].find(iter->first) == bk_boundary_y.at(ilevel)[1].end())
			{
				bk_boundary_y.at(ilevel)[1].insert({ iter->first, &grid_ptr.grid[iter->first] });
			}
		}
#if(C_DIMS == 3)
		if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
		{
			if (bk_boundary_z.at(ilevel)[0].find(iter->first) == bk_boundary_z.at(ilevel)[0].end())
			{
				bk_boundary_z.at(ilevel)[0].insert({ iter->first , &grid_ptr.grid[iter->first] });
			}
		}
		else if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
		{
			if (bk_boundary_z.at(ilevel)[1].find(iter->first) == bk_boundary_z.at(ilevel)[1].end())
			{
				bk_boundary_z.at(ilevel)[1].insert({ iter->first , &grid_ptr.grid[iter->first] });
			}
		}
#endif
#endif
	}

	// find nodes surround added nodes, and check if added nodes are on the numerical boundary, i.e. if there are nodes adjcent to the added nodes
#if (C_CHECK_MORTON_BOUNDARY == 1)
	bool bool_bx0, bool_bx1, bool_by0, bool_by1;
#if (C_DIMS==3)
	bool bool_bz0, bool_bz1;
#endif
#endif
	D_mapint map_surrounding_nodes = {};
	for (D_mapint::iterator iter = map_add_temp.begin(); iter != map_add_temp.end(); ++iter)
	{
		bool bool_merge = true;
#if (C_CHECK_MORTON_BOUNDARY == 1)
		bool_bx0 = true;  bool_bx1 = true;  bool_by0 = true;  bool_by1 = true;
#if (C_DIMS==3)
		bool_bz0 = true; bool_bz1 = true;
#endif
		if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
		{
			bool_bx0 = false;
		}
		else if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
		{
			bool_bx1 = false;
		}
		if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
		{
			bool_by0 = false;
		}
		else if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
		{
			bool_by1 = false;
		}
#if (C_DIMS==3)
		if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
		{
			bool_bz0 = false;
		}
		else if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
		{
			bool_bz1 = false;
		}
#endif
#endif
		// (-x) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_x0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp0_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp0_ilevel_1 , flag_refine });
				}
				morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if(C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
		// (-x, -y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#if(C_DIMS==3)
		// (-x, -y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by0 && bool_bz0)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp2_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp2_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, -y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by0 && bool_bz1)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp2_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp2_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_bz0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif

		// (-x, +y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (-x, +y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by1 && bool_bz0)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp2_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp2_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, +y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by1 && bool_bz1)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp2_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp2_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_bz1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif

		// (+x) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_x1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp0_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp0_ilevel_1 , flag_refine });
				}
				morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if(C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
					}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, -y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (+x, -y, -z)
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by0 && bool_bz0)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp2_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp2_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, -y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by0 && bool_bz1)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp2_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp2_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_bz0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#endif
		// (+x, +y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (+x, +y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by1 && bool_bz0)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp2_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp2_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, +y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by1 && bool_bz1)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp2_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp2_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, +z)
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_bz1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#endif
		// (-y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by0)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_y0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp0_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp0_ilevel_1 , flag_refine });
				}
				morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if(C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
					}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (-y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by0 && bool_bz0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by0 && bool_bz1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#endif
		// (+y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by1)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_y1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp0_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp0_ilevel_1 , flag_refine });
				}
				morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if(C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
					}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (+y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by1 && bool_bz0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by1 && bool_bz1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp1_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp1_ilevel_1 , flag_refine });
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bz0)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_z0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp0_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp0_ilevel_1 , flag_refine });
				}
				morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
					}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bz1)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_z1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end())
			{
				bool_merge = false;
			}
			else
			{
				if ((map_surrounding_nodes.find(morton_temp0_ilevel_1) == map_surrounding_nodes.end()))
				{
					map_surrounding_nodes.insert({ morton_temp0_ilevel_1 , flag_refine });
				}
				morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
					}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif
		if (!bool_merge)// added nodes are on iboundary2
		{
			grid_ptr.grid[iter->first].flag = flag_iboundary[iboundary2];
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).insert({ iter->first, flag_refine });
				if (map_add_nodes_out.find(iter->first) == map_add_nodes_out.end())
				{
					map_add_nodes_out.insert({ iter->first , flag_refine });
				}
			}
			grid_ptr.fine2coarse.at(iboundary2).insert({ iter->first, flag_refine });
		}
		else
		{
			grid_ptr.grid[iter->first].flag = flag_refine;
		}
	}

	// check if surrounding nodes are on the numerical boundary, including all the added nodes
	D_mapint map_not_boundary = {};
	D_mapint map_is_boundary = {};
	map_surrounding_nodes.insert(map_add_temp.begin(), map_add_temp.end());
	for (D_mapint::iterator iter = map_surrounding_nodes.begin(); iter != map_surrounding_nodes.end(); ++iter)
	{
		bool bool_not_bounadry = true;
#if (C_CHECK_MORTON_BOUNDARY == 1)
		bool_bx0 = true;  bool_bx1 = true;  bool_by0 = true;  bool_by1 = true;
#if (C_DIMS==3)
		bool_bz0 = true; bool_bz1 = true;
#endif
		if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
		{
			bool_bx0 = false;
		}
		else if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
		{
			bool_bx1 = false;
		}
		if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
		{
			bool_by0 = false;
		}
		else if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
		{
			bool_by1 = false;
		}
#if (C_DIMS==3)
		if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
		{
			bool_bz0 = false;
		}
		else if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
		{
			bool_bz1 = false;
		}
#endif
#endif
		// (-x) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_x0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end())
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, -y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (-x, -y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by0 && bool_bz0)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, -y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by0 && bool_bz1)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_bz0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif

		// (-x, +y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (-x, +y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by1 && bool_bz0)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, +y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by1 && bool_bz1)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_bz1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#endif
		// (+x) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_x1(iter->first, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, -y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (+x, -y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by0 && bool_bz0)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, -y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by0 && bool_bz1)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_bz0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#endif
		// (+x, +y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (+x, +y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by1 && bool_bz0)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, +y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by1 && bool_bz1)
		{
#endif
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_bz1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#endif
		// (-y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by0)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_y0(iter->first, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (-y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by0 && bool_bz0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by0 && bool_bz1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#endif
		// (+y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by1)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_y1(iter->first, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (+y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by1 && bool_bz0)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by1 && bool_bz1)
		{
#endif
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bz0)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_z0(iter->first, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bz1)
		{
#endif
			morton_temp0_ilevel_1 = Morton_Assist::find_z1(iter->first, ilevel_1);
			if (bool_not_bounadry && (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end()))
			{
				bool_not_bounadry = false;
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif
		if (bool_not_bounadry)
		{
			grid_ptr.grid[iter->first].flag = flag_refine;
			// delete from iboundary2_coarse and add to iboundary0_coarse 
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(iter->first) != grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).erase(iter->first);
				if (map_remove_nodes_out.find(iter->first) == map_remove_nodes_out.end())
				{
					map_remove_nodes_out.insert({ iter->first , flag_refine });
				}
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			// delete from iboundary2
			if (grid_ptr.fine2coarse.at(iboundary2).find(iter->first) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(iter->first);
			}
			map_not_boundary.insert({ iter->first, flag_refine });
		}
		else
		{
			grid_ptr.grid[iter->first].flag = flag_iboundary[iboundary2];
			map_is_boundary.insert({ iter->first, flag_iboundary[iboundary2] });
		}
	}

	// update boundary information based on surrounded nodes which are not on iboundary2 at (ilevel - 1)
	D_mapint map_boundary0 = {};
	D_mapint map_boundary1_erase = {};
	node_temp.flag = flag_refine;

	for (D_mapint::iterator iter = map_not_boundary.begin(); iter != map_not_boundary.end(); ++iter)
	{
#if (C_CHECK_MORTON_BOUNDARY == 1)
		bool_bx0 = true;  bool_bx1 = true;  bool_by0 = true;  bool_by1 = true;
#if (C_DIMS==3)
		bool_bz0 = true; bool_bz1 = true;
#endif
		if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
		{
			bool_bx0 = false;
		}
		else if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
		{
			bool_bx1 = false;
		}
		if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
		{
			bool_by0 = false;
		}
		else if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
		{
			bool_by1 = false;
		}
#if (C_DIMS==3)
		if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
		{
			bool_bz0 = false;
		}
		else if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
		{
			bool_bz1 = false;
		}
#endif
#endif
		// (-x) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0)
		{
#endif
			morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
			morton_temp0_ilevel_1 = Morton_Assist::find_x0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp0].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
					}
				}
				else
				{
					if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_refine;
					}
					if (map_boundary0.find(morton_temp0_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp0_ilevel_1, flag_refine });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, -y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by0)
		{
#endif
			morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
#if(C_DIMS==3)
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
#endif
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (-x, -y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by0 && bool_bz0)
		{
#endif
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp2) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp2 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp2].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp2_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp2_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
					}
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, -y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by0 && bool_bz1)
		{
#endif
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp2) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp2 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp2].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp2_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp2_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
					}
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_bz0)
		{
#endif
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif

		// (-x, +y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by1)
		{
#endif
			morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
#if(C_DIMS==3)
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
#endif
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (-x, +y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by1 && bool_bz0)
		{
#endif
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp2) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp2 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp2].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp2_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp2_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
					}
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, +y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_by1 && bool_bz1)
		{
#endif
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp2) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp2 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp2].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp2_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp2_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
					}
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-x, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx0 && bool_bz1)
		{
#endif
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif
		// (+x)
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1)
		{
#endif
			morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
			morton_temp0_ilevel_1 = Morton_Assist::find_x1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp0].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
					}
				}
				else
				{
					if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_refine;
					}
					if (map_boundary0.find(morton_temp0_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp0_ilevel_1, flag_refine });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, -y)
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by0)
		{
#endif
			morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
#if(C_DIMS==3)
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
#endif
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (+x, -y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by0 && bool_bz0)
		{
#endif
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp2) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp2 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp2].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp2_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp2_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
					}
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, -y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by0 && bool_bz1)
		{
#endif
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp2) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp2 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp2].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp2_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp2_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
					}
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_bz0)
		{
#endif
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif
		// (+x, +y)
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by1)
		{
#endif
			morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
#if(C_DIMS==3)
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
#endif
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (+x, +y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by1 && bool_bz0)
		{
#endif
			morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
			morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp2) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp2 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp2].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp2_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp2_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
					}
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, +y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_by1 && bool_bz1)
		{
#endif
			morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
			morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp2) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp2 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp2, &grid_ptr.grid[morton_temp2] });
						}
					}
					if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
					else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp2 , &grid_ptr.grid[morton_temp2] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp2].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp2_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp2_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
					}
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+x, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bx1 && bool_bz1)
		{
#endif
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif
		// (-y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by0)
		{
#endif
			morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
			morton_temp0_ilevel_1 = Morton_Assist::find_y0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp0].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
					}
				}
				else
				{
					if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_refine;
					}
					if (map_boundary0.find(morton_temp0_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp0_ilevel_1, flag_refine });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (-y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by0 && bool_bz0)
		{
#endif
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by0 && bool_bz1)
		{
#endif
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
				if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
					}
				}
				else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
					}
				}
				if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
					}
				}
				else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
					}
				}
				if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
					}
				}
				else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
					}
				}
#endif
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif
		// (+y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by1)
		{
#endif
			morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
			morton_temp0_ilevel_1 = Morton_Assist::find_y1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp0].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
					}
				}
				else
				{
					if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_refine;
					}
					if (map_boundary0.find(morton_temp0_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp0_ilevel_1, flag_refine });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

#if(C_DIMS==3)
		// (+y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by1 && bool_bz0)
		{
#endif
			morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_by1 && bool_bz1)
		{
#endif
			morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
			morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp1) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp1 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, &grid_ptr.grid[morton_temp1] });
						}
					}
					if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
					else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp1 , &grid_ptr.grid[morton_temp1] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp1].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
					}
				}
				else
				{
					if (map_boundary0.find(morton_temp1_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp1_ilevel_1, flag_add });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (-z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bz0)
		{
#endif
			morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
			morton_temp0_ilevel_1 = Morton_Assist::find_z0(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp0].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
					}
				}
				else
				{
					if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_refine;
					}
					if (map_boundary0.find(morton_temp0_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp0_ilevel_1, flag_refine });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif

		// (+z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
		if (bool_bz1)
		{
#endif
			morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
			morton_temp0_ilevel_1 = Morton_Assist::find_z1(iter->first, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				if (grid_ptr.grid.find(morton_temp0) == grid_ptr.grid.end())
				{
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
				}
				else
				{
					grid_ptr.grid[morton_temp0].flag = flag_refine;
				}
				if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
				{
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
						grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
					}
				}
				else
				{
					if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
					{
						grid_ptr.grid[morton_temp0].flag = flag_refine;
					}
					if (map_boundary0.find(morton_temp0_ilevel_1) == map_boundary0.end())
					{
						map_boundary0.insert({ morton_temp0_ilevel_1, flag_refine });
					}
					if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
					{
						grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
					}
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
				}
			}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		}
#endif
#endif
	}

	// delete nodes from iboundary0_coarse
	for (D_mapint::iterator iter = map_boundary0.begin(); iter != map_boundary0.end(); ++iter)
	{
		bool bool_not_iboundary0_coarse = true;
		//(-x)
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_x0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
				}
			}
			bool_not_iboundary0_coarse = false;
		}

		// (-x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}


#if(C_DIMS==3)
		// (-x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
				}
			}
		}

		// (-x, -y, +z)
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
				}
			}
		}

		// (-x, -z)
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}
#endif

		// (-x, +y)
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}


#if(C_DIMS==3)
		// (-x, +y, -z)
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
				}
			}
		}

		// (-x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
				}
			}
		}

		// (-x, +z)
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}

#endif
		// (+x)
		morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_x1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
				}
			}
		}

		// (+x, -y)
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}

#if(C_DIMS==3)
		// (+x, -y, -z)
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
				}
			}
		}

		// (+x, -y, +z)
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
				}
			}
		}

		// (+x, -z)
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}


#endif
		// (+x, +y)
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}

#if(C_DIMS==3)
		// (+x, +y, -z)
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
				}
			}
		}

		// (+x, +y, +z)
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp2].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp2, flag_add });
				}
			}
		}

		// (+x, +z)
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}

#endif
		// (-y)
		morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_y0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
				}
			}
		}

#if(C_DIMS==3)
		// (-y, -z)
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}

		// (-y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}

#endif
		// (+y)
		morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_y1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
				}
			}
		}

#if(C_DIMS==3)
		// (+y, -z)
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}

		// (+y, +z)
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp1].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp1, flag_add });
				}
			}
		}

		// (-z) 
		morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_z0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
				}
			}
		}

		// (+z)
		morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_z1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
			bool_not_iboundary0_coarse = false;
			if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
			{
				grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary1];
				if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).insert({ morton_temp0, flag_add });
				}
			}
		}

#endif
		if (bool_not_iboundary0_coarse)
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(iter->first);
				if (grid_ptr_coarse->grid.find(iter->first) != grid_ptr_coarse->grid.end())
				{
					grid_ptr_coarse->grid.erase(iter->first);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel_1)[0].find(iter->first) != bk_boundary_x.at(ilevel_1)[0].end())
						{
							bk_boundary_x.at(ilevel_1)[0].erase(iter->first);
						}
					}
					else if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel_1)[1].find(iter->first) != bk_boundary_x.at(ilevel_1)[1].end())
						{
							bk_boundary_x.at(ilevel_1)[1].erase(iter->first);
						}
					}
					if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel_1)[0].find(iter->first) != bk_boundary_y.at(ilevel_1)[0].end())
						{
							bk_boundary_y.at(ilevel_1)[0].erase(iter->first);
						}
					}
					else if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel_1)[1].find(iter->first) != bk_boundary_y.at(ilevel_1)[1].end())
						{
							bk_boundary_y.at(ilevel_1)[1].erase(iter->first);
						}
					}
#if(C_DIMS == 3)
					if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel_1)[0].find(iter->first) != bk_boundary_z.at(ilevel_1)[0].end())
						{
							bk_boundary_z.at(ilevel_1)[0].erase(iter->first);
						}
					}
					else if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel_1).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel_1)[1].find(iter->first) != bk_boundary_z.at(ilevel_1)[1].end())
						{
							bk_boundary_z.at(ilevel_1)[1].erase(iter->first);
						}
					}
#endif
#endif
				}

			}
		}
	}

#if(C_DIMS==3)
	D_mapint iboundary2_temp = {};
#endif
	// update boundary information based on nodes on iboundary2 at (ilevel - 1)
	D_mapint map_cell_delete = {};
	for (D_mapint::iterator iter = map_is_boundary.begin(); iter != map_is_boundary.end(); ++iter)
	{
		// (-x)
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_x0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
				}
				else
				{
					node_temp.flag = flag_iboundary[iboundary2];
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_add });
				}
			}
		}

		// (-x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{		
#if (C_DIMS == 2)
				if (map_cell_delete.find(morton_temp1) == map_cell_delete.end())
				{
						map_cell_delete.insert({ morton_temp1, flag_add });
				}
#endif
#if(C_DIMS==3)
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
#endif	
			}
		}

#if(C_DIMS==3)
		// (-x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_add });
				}
			}
		}

		// (-x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_add });
				}
			}
		}

		// (-x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
			}
		}
#endif

		// (-x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
#if (C_DIMS == 2)
				if (map_cell_delete.find(morton_temp1) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp1, flag_add });
				}
#endif
#if(C_DIMS==3)
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
#endif
			}
		}
#if(C_DIMS==3)
		// (-x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_add });
				}
			}
		}

		// (-x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_add });
				}
			}
		}

		// (-x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
			}
		}

#endif
		// (+x)
		morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_x1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
				}
				else
				{
					node_temp.flag = flag_iboundary[iboundary2];
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_iboundary[iboundary2] });
				}
			}
		}


		// (+x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
#if (C_DIMS == 2)
				if (map_cell_delete.find(morton_temp1) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp1, flag_add });
				}
#endif
#if(C_DIMS==3)
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
#endif
			}
		}
#if (C_DIMS == 3)
		// (+x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_add });
				}
			}
		}

		// (+x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_add });
				}
			}
		}
#endif
		// (+x, +y)
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
#if (C_DIMS == 2)
				if (map_cell_delete.find(morton_temp1) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp1, flag_add });
				}
#endif
#if(C_DIMS==3)
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
#endif
			}
		}
#if(C_DIMS==3)
		// (+x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_add });
				}
			}
		}

		// (+x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (map_cell_delete.find(morton_temp2) == map_cell_delete.end())
				{
					map_cell_delete.insert({ morton_temp2, flag_add });
				}
			}
		}

		// (+x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
			}
		}

		// (+x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
			}
		}
#endif
		// (-y) 
		morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_y0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
				}
				else
				{
					node_temp.flag = flag_iboundary[iboundary2];
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_iboundary[iboundary2] });
				}
			}
		}

#if(C_DIMS==3)
		// (-y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
			}
		}

		// (-y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
			}
		}

#endif
		// (+y) 
		morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_y1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
				}
				else
				{
					node_temp.flag = flag_iboundary[iboundary2];
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
#if (C_DIMS == 3)
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
#endif
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_iboundary[iboundary2] });
				}
			}
		}

#if(C_DIMS==3)
		// (+y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
			}
		}

		// (+y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (iboundary2_temp.find(morton_temp1) == iboundary2_temp.end())
				{
					iboundary2_temp.insert({ morton_temp1, flag_add });
				}
			}
		}

		// (-z) 
		morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_z0(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
				}
				else
				{
					node_temp.flag = flag_iboundary[iboundary2];
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_iboundary[iboundary2] });
				}
			}
		}

		// (+z) 
		morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_z1(iter->first, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
			{
				if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
				{
					grid_ptr.grid[morton_temp0].flag = flag_iboundary[iboundary2];
				}
				else
				{
					node_temp.flag = flag_iboundary[iboundary2];
					grid_ptr.grid.insert({ morton_temp0 , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
					if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
					{
						if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) == bk_boundary_x.at(ilevel)[0].end())
						{
							bk_boundary_x.at(ilevel)[0].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
					{
						if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) == bk_boundary_x.at(ilevel)[1].end())
						{
							bk_boundary_x.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
					{
						if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) == bk_boundary_y.at(ilevel)[0].end())
						{
							bk_boundary_y.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
					{
						if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) == bk_boundary_y.at(ilevel)[1].end())
						{
							bk_boundary_y.at(ilevel)[1].insert({ morton_temp0, &grid_ptr.grid[morton_temp0] });
						}
					}
					if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
					{
						if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) == bk_boundary_z.at(ilevel)[0].end())
						{
							bk_boundary_z.at(ilevel)[0].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
					else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
					{
						if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) == bk_boundary_z.at(ilevel)[1].end())
						{
							bk_boundary_z.at(ilevel)[1].insert({ morton_temp0 , &grid_ptr.grid[morton_temp0] });
						}
					}
#endif
				}
				if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) == grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.fine2coarse.at(iboundary2).insert({ morton_temp0, flag_iboundary[iboundary2] });
				}
			}
		}

#endif

	}

	// delete isolate cell, i.e. all its faces (3D) or edges (2D) are on iboundary2 
	D_morton morton_ref0 = ~Morton_Assist::morton_xyz.at(ilevel);
	D_morton morton_corner;
	D_mapint map_iboudnary2_delete_temp = {};
	//	@todo use flag,i.e. x0x1y0y1z0z1 in previous for loop, to identify if flag of nodes is flag_iboundary[2]
	for (D_mapint::iterator iter = map_cell_delete.begin(); iter != map_cell_delete.end(); ++iter)
	{
		bool bool_need2delete = true;
		// (x0, y0, z0)
		morton_corner = iter->first & morton_ref0;
		if (grid_ptr.grid.find(morton_corner) == grid_ptr.grid.end() || grid_ptr.grid[morton_corner].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

		// (x1, y0, z0)
		morton_temp0_ilevel_1 = Morton_Assist::find_x1(morton_corner, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

		// (x0, y1, z0)
		morton_temp0_ilevel_1 = Morton_Assist::find_y1(morton_corner, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

		// (x1, y1, z0)
		morton_temp1_ilevel_1 = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

#if (C_DIMS==3)
		// (x0, y0, z1)
		morton_temp0_ilevel_1 = Morton_Assist::find_z1(morton_corner, ilevel_1);
		if (grid_ptr.grid.find(morton_temp0_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp0_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}
		// (x1, y0, z1)
		morton_temp1_ilevel_1 = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

		// (x0, y1, z1)
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp1_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp1_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}

		// (x1, y1, z1)
		morton_temp2_ilevel_1 = Morton_Assist::find_x1(morton_temp1_ilevel_1, ilevel_1);
		if (grid_ptr.grid.find(morton_temp2_ilevel_1) == grid_ptr.grid.end() || grid_ptr.grid[morton_temp2_ilevel_1].flag == flag_iboundary[iboundary2])
		{
		}
		else
		{
			bool_need2delete = false;
		}
#endif

		if (bool_need2delete)
		{
			D_morton morton_temp0_inner, morton_temp1_inner;
#if(C_DIMS==3)
			D_morton morton_temp2_inner;
#endif
			// (x0, y0, z0)
			if (grid_ptr.grid.find(morton_corner) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_corner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_corner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_corner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_corner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
#if (C_DIMS == 3)
				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_corner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_corner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_corner) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_corner, flag_add });
					}
				}
			}

			// (x1, y0, z0)
			morton_temp0_ilevel_1 = Morton_Assist::find_x1(morton_corner, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp0_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
#if (C_DIMS == 3)
				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp0_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp0_ilevel_1, flag_add });
					}
				}
			}


			// (x0, y1, z0)
			morton_temp0_ilevel_1 = Morton_Assist::find_y1(morton_corner, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp0_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
#if (C_DIMS == 3)
				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp0_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp0_ilevel_1, flag_add });
					}
				}
			}

			// (x1, y1, z0)
			morton_temp1_ilevel_1 = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp1_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#endif
				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

#if(C_DIMS==3)
				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#if (C_DIMS == 3)
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
#if (C_DIMS == 3)
				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
#endif
				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp1_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp1_ilevel_1, flag_add });
					}
				}
			}
#if (C_DIMS==3)
			// (x0, y0, z1)
			morton_temp0_ilevel_1 = Morton_Assist::find_z1(morton_corner, ilevel_1);
			if (grid_ptr.grid.find(morton_temp0_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp0_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp0_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp0_ilevel_1, flag_add });
					}
				}
			}

			// (x1, y0, z1)
			morton_temp1_ilevel_1 = Morton_Assist::find_x1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp1_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp1_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp1_ilevel_1, flag_add });
					}
				}
			}


			// (x0, y1, z1)
			morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp1_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp1_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}


				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp1_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp1_ilevel_1, flag_add });
					}
				}
			}


			// (x1, y1, z1)
			morton_temp2_ilevel_1 = Morton_Assist::find_x1(morton_temp1_ilevel_1, ilevel_1);
			if (grid_ptr.grid.find(morton_temp2_ilevel_1) != grid_ptr.grid.end())
			{
				bool bool_need2delete_inner = true;
				// (-x) 
				morton_temp0_inner = Morton_Assist::find_x0(morton_temp2_ilevel_1, ilevel_1);
				if (grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x) 
				morton_temp0_inner = Morton_Assist::find_x1(morton_temp2_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y) 
				morton_temp1_inner = Morton_Assist::find_y0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, -y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y) 
				morton_temp1_inner = Morton_Assist::find_y1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, -z) 
				morton_temp2_inner = Morton_Assist::find_z0(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+x, +y, +z) 
				morton_temp2_inner = Morton_Assist::find_z1(morton_temp1_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp2_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp2_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-y) 
				morton_temp0_inner = Morton_Assist::find_y0(morton_temp2_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (+y) 
				morton_temp0_inner = Morton_Assist::find_y1(morton_temp2_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z0(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				morton_temp1_inner = Morton_Assist::find_z1(morton_temp0_inner, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp1_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				// (-z) 
				morton_temp0_inner = Morton_Assist::find_z0(morton_temp2_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}
				// (+z) 
				morton_temp0_inner = Morton_Assist::find_z1(morton_temp2_ilevel_1, ilevel_1);
				if (bool_need2delete_inner && grid_ptr.grid.find(morton_temp0_inner) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0_inner].flag != flag_iboundary[iboundary2])
				{
					bool_need2delete_inner = false;
				}

				if (bool_need2delete_inner)
				{
					if (map_iboudnary2_delete_temp.find(morton_temp2_ilevel_1) == map_iboudnary2_delete_temp.end())
					{
						map_iboudnary2_delete_temp.insert({ morton_temp2_ilevel_1, flag_add });
					}
				}
			}
#endif

		}
	}

	// delete boundary information based on nodes got in previous for loop
	for (D_mapint::iterator iter = map_iboudnary2_delete_temp.begin(); iter != map_iboudnary2_delete_temp.end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(iter->first);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(iter->first) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(iter->first);
				}
			}
			else if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(iter->first) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(iter->first);
				}
			}
			if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(iter->first) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(iter->first);
				}
			}
			else if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(iter->first) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(iter->first);
				}
			}
#if(C_DIMS == 3)
			if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(iter->first) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(iter->first);
				}
			}
			else if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(iter->first) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(iter->first);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(iter->first) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(iter->first);
			}
			if (grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).find(iter->first) != grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).erase(iter->first);
				if (map_remove_nodes_out.find(iter->first) == map_remove_nodes_out.end())
				{
					map_remove_nodes_out.insert({ iter->first , flag_refine });
				}
			}
		}

		// (-x) 
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

		// (-x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#if(C_DIMS==3)
		// (-x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (-x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (-x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}
#endif

		// (-x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#if(C_DIMS==3)
		// (-x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (-x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (-x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#endif
		// (+x) 
		morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

		// (+x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#if(C_DIMS==3)
		// (+x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (+x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (+x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}
#endif

		// (+x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#if(C_DIMS==3)
		// (+x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (+x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (grid_ptr.grid.find(morton_temp2) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp2);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp2) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp2) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp2) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp2) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp2);
				}
			}
			if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp2) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp2);
				}
			}
			else if ((morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp2) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp2);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp2);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp2) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp2);
			}
		}

		// (+x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#endif
		// (-y) 
		morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

#if(C_DIMS==3)
		// (-y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

		// (-y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

#endif
		// (+y) 
		morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

#if(C_DIMS==3)
		// (+y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

		// (+y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp1);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp1) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp1) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp1);
				}
			}
			if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp1) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp1) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp1);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp1) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp1);
				}
			}
			else if ((morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp1) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp1);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp1);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp1) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp1);
			}
		}

		// (-z) 
		morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

		// (+z) 
		morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
		if (grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end())
		{
			grid_ptr.grid.erase(morton_temp0);
#if ((C_CHECK_MORTON_BOUNDARY == 1) && (C_CHECK_BOUDNARY_METHOD == 1))
			if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (bk_boundary_x.at(ilevel)[0].find(morton_temp0) != bk_boundary_x.at(ilevel)[0].end())
				{
					bk_boundary_x.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (bk_boundary_x.at(ilevel)[1].find(morton_temp0) != bk_boundary_x.at(ilevel)[1].end())
				{
					bk_boundary_x.at(ilevel)[1].erase(morton_temp0);
				}
			}
			if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (bk_boundary_y.at(ilevel)[0].find(morton_temp0) != bk_boundary_y.at(ilevel)[0].end())
				{
					bk_boundary_y.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (bk_boundary_y.at(ilevel)[1].find(morton_temp0) != bk_boundary_y.at(ilevel)[1].end())
				{
					bk_boundary_y.at(ilevel)[1].erase(morton_temp0);
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (bk_boundary_z.at(ilevel)[0].find(morton_temp0) != bk_boundary_z.at(ilevel)[0].end())
				{
					bk_boundary_z.at(ilevel)[0].erase(morton_temp0);
				}
			}
			else if ((morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (bk_boundary_z.at(ilevel)[1].find(morton_temp0) != bk_boundary_z.at(ilevel)[1].end())
				{
					bk_boundary_z.at(ilevel)[1].erase(morton_temp0);
				}
			}
#endif
#endif
			if (grid_ptr.fine2coarse.at(iboundary2).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary2).end())
			{
				grid_ptr.fine2coarse.at(iboundary2).erase(morton_temp0);
			}
			if (grid_ptr.fine2coarse.at(iboundary1).find(morton_temp0) != grid_ptr.fine2coarse.at(iboundary1).end())
			{
				grid_ptr.fine2coarse.at(iboundary1).erase(morton_temp0);
			}
		}

#endif
	}

#if (C_DIMS == 3)
	// update iboundary2 information on (xy, xz, or yz) planes
	node_temp.flag = flag_iboundary[iboundary2];
	bool bool_xy, bool_yz, bool_xz;
	bool bool_x_temp0, bool_y_temp0, bool_z_temp0, bool_x_temp1, bool_y_temp1, bool_z_temp1;
	for (D_mapint::iterator iter = iboundary2_temp.begin(); iter != iboundary2_temp.end(); ++iter)
	{
		if (grid_ptr.grid.find(iter->first) != grid_ptr.grid.end())
		{
			if (grid_ptr.grid[iter->first].flag != flag_iboundary[iboundary1])
			{
				if (grid_ptr.fine2coarse.at(iboundary2).find(iter->first) == grid_ptr.fine2coarse.at(iboundary2).end())
				{
					grid_ptr.grid[iter->first].flag = flag_iboundary[iboundary2];
					grid_ptr.fine2coarse.at(iboundary2).insert({ iter->first, flag_add });
				}
			}
			else
			{
				morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
				bool_x_temp0 = grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2];
				if (bool_x_temp0)
				{
					morton_temp1 = Morton_Assist::find_x1(iter->first, ilevel);
					bool_x_temp1 = grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2];
				}

				morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
				bool_y_temp0 = grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2];
				if (bool_y_temp0)
				{
					morton_temp1 = Morton_Assist::find_y1(iter->first, ilevel);
					bool_y_temp1 = grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2];
				}

				morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
				bool_z_temp0 = grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2];
				if (bool_z_temp0)
				{
					morton_temp1 = Morton_Assist::find_z1(iter->first, ilevel);
					bool_z_temp1 = grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2];
				}
				// check if node is on any of the plane which is on iboudanry2
				bool_xy = bool_x_temp0 && bool_x_temp1 && bool_y_temp0 && bool_y_temp1;
				bool_xz = bool_x_temp0 && bool_x_temp1 && bool_z_temp0 && bool_z_temp1;
				bool_yz = bool_y_temp0 && bool_y_temp1 && bool_z_temp0 && bool_z_temp1;
				if (bool_xy || bool_xz || bool_yz)
				{

					grid_ptr.grid[iter->first].flag = flag_iboundary[iboundary2];
					grid_ptr.fine2coarse.at(iboundary2).insert({ iter->first, flag_add });
				}
			}
		}
		else
		{
			morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
			bool_x_temp0 = grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2];
			if (bool_x_temp0)
			{
				morton_temp1 = Morton_Assist::find_x1(iter->first, ilevel);
				bool_x_temp1 = grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2];
			}

			morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
			bool_y_temp0 = grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2];
			if (bool_y_temp0)
			{
				morton_temp1 = Morton_Assist::find_y1(iter->first, ilevel);
				bool_y_temp1 = grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2];
			}

			morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
			bool_z_temp0 = grid_ptr.grid.find(morton_temp0) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp0].flag == flag_iboundary[iboundary2];
			if (bool_z_temp0)
			{
				morton_temp1 = Morton_Assist::find_z1(iter->first, ilevel);
				bool_z_temp1 = grid_ptr.grid.find(morton_temp1) != grid_ptr.grid.end() && grid_ptr.grid[morton_temp1].flag == flag_iboundary[iboundary2];
			}
			// check if node is on any of the plane which is on iboudanry2
			bool_xy = bool_x_temp0 && bool_x_temp1 && bool_y_temp0 && bool_y_temp1;
			bool_xz = bool_x_temp0 && bool_x_temp1 && bool_z_temp0 && bool_z_temp1;
			bool_yz = bool_y_temp0 && bool_y_temp1 && bool_z_temp0 && bool_z_temp1;
			if (bool_xy || bool_xz || bool_yz)
			{
				grid_ptr.grid.insert({ iter->first , node_temp });
#if (C_CHECK_MORTON_BOUNDARY == 1)
				if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
				{
					if (bk_boundary_x.at(ilevel)[0].find(iter->first) == bk_boundary_x.at(ilevel)[0].end())
					{
						bk_boundary_x.at(ilevel)[0].insert({ iter->first, &grid_ptr.grid[iter->first] });
					}
				}
				else if ((iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
				{
					if (bk_boundary_x.at(ilevel)[1].find(iter->first) == bk_boundary_x.at(ilevel)[1].end())
					{
						bk_boundary_x.at(ilevel)[1].insert({ iter->first , &grid_ptr.grid[iter->first] });
					}
				}
				if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
				{
					if (bk_boundary_y.at(ilevel)[0].find(iter->first) == bk_boundary_y.at(ilevel)[0].end())
					{
						bk_boundary_y.at(ilevel)[0].insert({ iter->first , &grid_ptr.grid[iter->first] });
					}
				}
				else if ((iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
				{
					if (bk_boundary_y.at(ilevel)[1].find(iter->first) == bk_boundary_y.at(ilevel)[1].end())
					{
						bk_boundary_y.at(ilevel)[1].insert({ iter->first, &grid_ptr.grid[iter->first] });
					}
				}
				if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
				{
					if (bk_boundary_z.at(ilevel)[0].find(iter->first) == bk_boundary_z.at(ilevel)[0].end())
					{
						bk_boundary_z.at(ilevel)[0].insert({ iter->first , &grid_ptr.grid[iter->first] });
					}
				}
				else if ((iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
				{
					if (bk_boundary_z.at(ilevel)[1].find(iter->first) == bk_boundary_z.at(ilevel)[1].end())
					{
						bk_boundary_z.at(ilevel)[1].insert({ iter->first , &grid_ptr.grid[iter->first] });
					}
				}
#endif
				grid_ptr.fine2coarse.at(iboundary2).insert({ iter->first, flag_add });
				if (grid_ptr.fine2coarse.at(iboundary1).find(iter->first) != grid_ptr.fine2coarse.at(iboundary1).end())
				{
					grid_ptr.fine2coarse.at(iboundary1).erase(iter->first);
				}
			}
		}
	}
#endif
}

#endif