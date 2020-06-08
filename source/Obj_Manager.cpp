/**
* @file
* @author Zhengliang Liu
* @brief Main function.
* @note .
*/
#include "Obj_Manager.h"
Grid_Manager* Grid_Manager::pointer_me;
IO_Manager* IO_Manager::pointer_me;
Solid_Manager* Solid_Manager::pointer_me;

// declarations
Grid_Manager gr_manager;
Solid_Manager soild_manager;
IO_Manager io_manager;



void Obj_Manager::initial()
{

	// Initialize solid geometry
	std::vector<Ini_Shape> ini_shape(1); // number of solids

	ini_shape.at(0).shape_type = geofile;
	ini_shape.at(0).bool_moving = true;
	ini_shape.at(0).length.push_back(C_xb / 2);
    ini_shape.at(0).length.push_back(C_yb / 2);
	ini_shape.at(0).length.push_back(C_zb / 2);

	//ini_shape.at(0).shape_type = circle;
	//ini_shape.at(0).bool_moving = true;
	//ini_shape.at(0).numb_nodes = 2000;
	//ini_shape.at(0).length.push_back(C_xb / 2);
	//ini_shape.at(0).length.push_back(C_yb / 2);
	//ini_shape.at(0).length.push_back(1);


	//ini_shape.at(1).shape_type = line_fillx;
	//ini_shape.at(1).x0 = 20.;
	//ini_shape.at(1).y0 = 20.;
	//ini_shape.at(1).numb_nodes = 200;
	//ini_shape.at(1).length.push_back(-1.);
	//ini_shape.at(1).length.push_back(20.);

	Solid_Manager::pointer_me = &soild_manager;
	soild_manager.initial(ini_shape);

	std::vector <Ini_Shape>().swap(ini_shape);

	// Generate and initialize mesh
	Grid_Manager::pointer_me = &gr_manager;
	gr_manager.initial();

	// dermine run order of each refinement level
	std::array<unsigned int, C_max_level + 1> accumulate_t = {};
	if (C_max_level == 1)
	{
		run_order.push_back(0);
		run_order.push_back(C_max_level);
		run_order.push_back(C_max_level);
	}
	else
	{
		unsigned int ilevel = 1;
		accumulate_t[0] = two_power_n(C_max_level);
		run_order.push_back(0);
		while (ilevel > 0)
		{
			if (accumulate_t[ilevel- 1] == accumulate_t[ilevel])
			{
				--ilevel;
			}
			else
			{
				run_order.push_back(ilevel);
				accumulate_t[ilevel] += two_power_n(C_max_level - ilevel);
				
				if (ilevel + 1 < C_max_level)
				{
					++ilevel;
				}
				else
				{
					run_order.push_back(ilevel + 1);
					run_order.push_back(ilevel + 1);
					accumulate_t[ilevel + 1] += two_power_n(C_max_level - ilevel);
					++ilevel;
				}
			}
		}
	}
}

void Obj_Manager::time_marching_management()
{
	D_real sum_t = 0.;
	std::array<D_mapint, C_max_level + 1> map_add_nodes, map_remove_nodes;
	while (sum_t < 3 * C_dx)
	{
		time_marching(sum_t, map_add_nodes, map_remove_nodes);
		sum_t += C_dx;
	}

	
}

void Obj_Manager::time_marching(D_real sum_t, std::array<D_mapint, C_max_level + 1>  &map_add_nodes, std::array<D_mapint, C_max_level + 1>  &map_remove_nodes)
{
	D_real dt = C_dx / static_cast<D_real> (two_power_n(C_max_level));
	std::array<unsigned int, C_max_level + 1> accumulate_t{};
	std::array<unsigned int, C_max_level + 1> flag_time_step{}; // record number of time step at ilevel

	for (std::vector<unsigned int>::iterator iter = run_order.begin(); iter != run_order.end(); ++iter)
	{
		Timer tmr;
		double t0 = tmr.elapsed();

		unsigned int ilevel = *iter;
		++flag_time_step[ilevel];

		D_real current_t;

		accumulate_t[ilevel] += two_power_n(C_max_level - ilevel);
		current_t = sum_t + dt * static_cast<D_real> (accumulate_t[ilevel]);

		if (ilevel == C_max_level)
		{
#if (C_SOLID_BOUNDARY == 2)
			// update shape information
			unsigned int numb_solids = Solid_Manager::pointer_me->numb_solids;
			for (unsigned int ishape = 0; ishape < numb_solids; ++ishape)
			{
				if (Solid_Manager::pointer_me->shape_solids.at(ishape).bool_moving)
				{
					Solid_Manager::pointer_me->renew(ishape, current_t);				
					Grid_Manager::pointer_me->update_nodes_near_solid(ishape, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel));
					Grid_Manager::pointer_me->update_map_node_IB(ilevel, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel));

					if (Solid_Manager::pointer_me->shape_solids.at(ishape).bool_enclosed)
					{
						Grid_Manager::pointer_me->update_ghost_node(ilevel, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel));
					}
				}

			}
			
			if ((flag_time_step[ilevel] % 2) == 0)
			{
				Grid_Manager::pointer_me->call_update_nodes(ilevel, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel), map_add_nodes.at(ilevel - 1), map_remove_nodes.at(ilevel - 1));
				map_add_nodes[ilevel].clear();
				map_remove_nodes[ilevel].clear();
			}

#endif			

		}
		else if (ilevel == 0)
		{

		}
		else
		{
#if (C_SOLID_BOUNDARY == 2)
			if ((flag_time_step[ilevel] % 2) == 0)
			{
				Grid_Manager::pointer_me->call_update_nodes(ilevel, map_add_nodes[ilevel], map_remove_nodes[ilevel], map_add_nodes[ilevel - 1], map_remove_nodes[ilevel - 1]);
				map_add_nodes[ilevel].clear();
				map_remove_nodes[ilevel].clear();
			}
#endif	
		}

		//double t1 = tmr.elapsed();
		//double t2 = tmr.elapsed();
		//std::cout << ilevel<< ": " << t1 - t0 << ", " << t2 - t1 << std::endl;
	}
	
}

void Obj_Manager::output()
{
	// Write flowfield
	IO_Manager::pointer_me = &io_manager;
	io_manager.method = 4;
	io_manager.outfile = "test";
	//std::vector<unsigned int> out_vlevel = { C_max_level};
	//std::vector<unsigned int> out_vlevel = { 2, 1, 0 };
	std::vector<unsigned int> out_vlevel;
	for (unsigned ilevel = 0; ilevel < C_max_level + 1; ++ilevel)
	{
		out_vlevel.push_back(ilevel);
	}
	io_manager.vlevel = out_vlevel;

	io_manager.control();
}