/**
* @file
* @author Zhengliang Liu
* @brief Mesh generation.
* @note .
*/
#include "General.h"
#include "Grid_Manager.h"
#include "Solid_Manager.h"
#include "Morton_assist.h"
//////////////////////////////////////
//              C_BIT               //
//  {---------------------------}   //
//  xxxxxxxxxxxxxx xxxxxxxxxxxxxx   //
//  {------------} {------------}   //
//  bit_background bit_otherlevel   //
///////////////////////////////
Morton_Assist* Morton_Assist::pointer_me;
std::array<D_morton, C_max_level + 1> Morton_Assist::ref_one;
std::array<std::array<D_morton, 2>, C_max_level + 1> Morton_Assist::xmorton_ref;
std::array<std::array<D_morton, 2>, C_max_level + 1> Morton_Assist::ymorton_ref;
#if (C_DIMS == 3)
std::array<std::array<D_morton, 2>, C_max_level + 1> Morton_Assist::zmorton_ref;
#endif
std::array<D_morton, C_max_level + 1> Morton_Assist::morton_xyz;
std::array<D_morton, C_max_level + 1> Morton_Assist::morton_xyz_flip;
D_morton Morton_Assist::mortonx_min, Morton_Assist::mortonx_max;
D_morton Morton_Assist::mortony_min, Morton_Assist::mortony_max;
#if (C_DIMS == 3)
D_morton Morton_Assist::mortonz_min, Morton_Assist::mortonz_max;
#endif
int flag0 = 0;
int flag_ghost_temp = two_power_n(0);

/**
* @brief function to generate initial mesh.
*/
void Grid_Manager::initial()
{
	Morton_Assist mor_assist;
	Morton_Assist::pointer_me = &mor_assist;
	Morton_Assist::pointer_me->morton_initial();
	// compuate refence morton code xmorton_ref, ymorton_ref and zmorton_ref
	for (unsigned int ilevel = 0; ilevel < C_max_level + 1; ++ilevel)
	{
		Morton_Assist::compute_ref(ilevel);
	}

	D_uint nix = static_cast <D_uint> (C_xb / C_dx + C_eps) + C_x0b_offset;
	D_uint niy = static_cast <D_uint> (C_yb / C_dx + C_eps) + C_y0b_offset;
#if (C_DIMS==2)
	Morton_Assist::mortonx_min = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(C_x0b_offset, 0));
	Morton_Assist::mortony_min = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode( 0, C_y0b_offset));
#endif
#if (C_DIMS==3)
	Morton_Assist::mortonx_min = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(C_x0b_offset, 0, 0));
	Morton_Assist::mortony_min = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(0, C_y0b_offset, 0));
	Morton_Assist::mortonz_min = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(0, 0, C_z0b_offset));
#endif
	Morton_Assist::mortonx_min = Morton_Assist::mortonx_min << Morton_Assist::bit_otherlevel;
	Morton_Assist::mortony_min = Morton_Assist::mortony_min << Morton_Assist::bit_otherlevel;

	xb_domain = C_dx * static_cast <D_real> (nix);
	if (C_xb > xb_domain + C_eps)
	{
		std::stringstream infor;
		infor << "set the maximum domain boundary in x direction as xb_domain = " << xb_domain << " instead of C_xb = " << C_xb << " in settings." << std::endl;
		log_infor(infor.str(), Log_function::logfile);
	}

	yb_domain = C_dx * static_cast <D_real> (niy);
	if (C_yb > yb_domain + C_eps)
	{
		std::stringstream infor;
		infor << "set the maximum domain boundary in y direction as yb_domain = " << yb_domain << " instead of C_yb = " << C_yb << " in settings." << std::endl;
		log_infor(infor.str(), Log_function::logfile);

	}
#if (C_DIMS==2)
	Morton_Assist::mortonx_max = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(nix, 0));
	Morton_Assist::mortony_max = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(0, niy));
#endif
#if (C_DIMS == 3)
	Morton_Assist::mortonx_max = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(nix, 0, 0));
	Morton_Assist::mortony_max = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(0, niy, 0));
#endif
	Morton_Assist::mortonx_max = Morton_Assist::mortonx_max << Morton_Assist::bit_otherlevel;
	Morton_Assist::mortony_max = Morton_Assist::mortony_max << Morton_Assist::bit_otherlevel;

#if (C_DIMS == 3)

	D_uint niz = static_cast <D_uint> (C_zb / C_dx + C_eps) + C_z0b_offset;
	zb_domain = C_dx * static_cast <D_real> (niz);
	if (C_zb > zb_domain + C_eps)
	{
		std::stringstream infor;
		infor << "set the maximum domain boundary in z direction as zb_domain = " << zb_domain << " instead of C_zb = " << C_zb << " in settings." << std::endl;
		log_infor(infor.str(), Log_function::logfile);
	}
	
	Morton_Assist::mortonz_min = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(0, 0, C_z0b_offset));
	Morton_Assist::mortonz_min = Morton_Assist::mortonz_min << Morton_Assist::bit_otherlevel;
	Morton_Assist::mortonz_max = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(0, 0, niz));
	Morton_Assist::mortonz_max = Morton_Assist::mortonz_max << Morton_Assist::bit_otherlevel;
#endif

	generate_inner();
	

	std::cout << "ilevel = " << C_max_level << ", No. of nodes:" << gr_inner.grid.size() << std::endl;
	for (unsigned int ilevel = C_max_level - 1; ilevel > 0; --ilevel)
	{
		generate_intermediate(ilevel);
		std::cout << "ilevel = " << ilevel << ", No. of nodes:" << gr_NoIB[ilevel].grid.size() << std::endl;
	}
	
	generate_background();

	std::cout << "ilevel = " << 0 << ", No. of nodes:" << gr_NoIB[0].grid.size() << std::endl;

	identify_domain_boundary();


#if (C_SEARCH_METHOD != 1)
	for (unsigned int ilevel = 1; ilevel <= C_max_level; ++ilevel)
	{
		if (ilevel == C_max_level)
		{
			Node_IB node_temp;
			initial_icount_refine(ilevel, Grid_Manager::pointer_me->gr_inner);
		}
		else
		{		
			Node node_temp;
			initial_icount_refine(ilevel, Grid_Manager::pointer_me->gr_NoIB.at(ilevel));
		}
	}
#endif

	gr_NoIB[0].fine2coarse.at(2).clear();
	gr_NoIB[0].fine2coarse.at(2) = gr_NoIB[0].coarse2fine.at(1);
	gr_NoIB[0].fine2coarse.at(1).clear();
	gr_NoIB[0].fine2coarse.at(1) = gr_NoIB[0].coarse2fine.at(0);

}

/**
* @brief function to generate nodes at the fineset level, i.e. ilevel = C_max_level.
*/
void Grid_Manager::generate_inner()
{
	unsigned int ilevel = C_max_level ;
	gr_inner.dx = C_dx / static_cast<D_real> (two_power_n(C_max_level));

	// find fluid nodes near the solid
	D_mapint2 refine_nodes = {}; ///< nodes need to refine. the first array element is the flag for identifying node features, the second is the number of solid points replis on this node
	D_mapint2 cells_level1 = {}; // cells centers at (ilevel - 1)
	std::array<D_mapint, 2> boundary_x_temp, boundary_y_temp, boundary_z_temp; // temporary array to store domain boundaries at ilevel, only work if C_CHECK_MORTON_BOUNDARY = 1
	for (unsigned int ishape = 0; ishape < Solid_Manager::pointer_me->numb_solids; ++ishape)
	{
		D_mapint nearest_center = {};
		std::array<D_real, C_DIMS> xyz;
#if (C_CHECK_MORTON_BOUNDARY == 1)
		std::array<D_real, C_DIMS> temp_xyz; // the latest xyz which does not exceed the computational domain
#endif
		for (D_uint ipoint = 0; ipoint < Solid_Manager::pointer_me->shape_solids.at(ishape).numb_nodes; ++ipoint)
		{
			xyz.at(0) = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).x;
			xyz.at(1) = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).y;
#if(C_DIMS == 3)
			xyz.at(2) = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).z;
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
			bool exceed_domain = false;
			if (xyz.at(0) > xb_domain + C_eps)
			{
				std::stringstream warning;
				warning << "coordinate x = " << xyz.at(0) << " of point " << ipoint << " in shape "<< ishape << " exceeds the compuational domain where xb_domain  = " << xb_domain << std::endl;
				log_warning(warning.str(), Log_function::logfile);
				exceed_domain = true;
			}
			if (xyz.at(1) > yb_domain + C_eps)
			{
				std::stringstream warning;
				warning << "coordinate y = " << xyz.at(1) << " of point " << ipoint << " in shape " << ishape << " exceeds the compuational domain where yb_domain  = " << yb_domain << std::endl;
				log_warning(warning.str(), Log_function::logfile);
				exceed_domain = true;
			}
#if(C_DIMS == 2)
			if (exceed_domain)
			{
				xyz = temp_xyz;
			}
			else
			{
				temp_xyz = xyz;
			}
#endif
#if(C_DIMS == 3)
			if (xyz.at(2) > zb_domain  + C_eps)
			{
				std::stringstream warning;
				warning << "coordinate z = " << xyz.at(2) << "of point " << ipoint << " in shape " << ishape << " exceeds the compuational domain where zb_domain  = " << zb_domain  << std::endl;
				log_warning(warning.str(), Log_function::logfile);
				exceed_domain = true;
			}
			if (exceed_domain)
			{
				xyz = temp_xyz;
			}
			else
			{
				temp_xyz = xyz;
			}
#endif
#endif
			find_nodes_near_solid(ilevel, xyz, refine_nodes);

#if	(C_SOLID_BOUNDARY == 2)
			initial_ref_coordinate(ishape, ipoint);
#endif
		}
		
		
		initial_map_node_IB(ilevel);

		D_morton morton_xyz = Morton_Assist::morton_xyz.at(ilevel);
		D_morton morton_xyz_1 = morton_xyz;
		morton_xyz_1.flip();
		D_morton morton_temp;

		for (D_map_Infor_near_solid::iterator iter = map_node_near_solid.begin(); iter != map_node_near_solid.end(); ++iter)
		{
			morton_temp = iter->first & morton_xyz_1;
			if (nearest_center.find(morton_temp) == nearest_center.end())
			{
				nearest_center.insert({ morton_temp, flag_refine });
			}
		}

		// mark nodes within C_extend_ghost grip space from nodes in (map_near_solid)
		if (Solid_Manager::pointer_me->shape_solids.at(ishape).bool_enclosed == true)
		{
			mark_ghost_points(ilevel, refine_nodes);
		}

		// find cells at (ilevel - 1) near the solid
		D_mapint nodes_level1 = nearest_center;

		Timer tmr;
		search_nodes_near_solid(ilevel, nearest_center, nodes_level1, refine_nodes, boundary_x_temp, boundary_y_temp, boundary_z_temp);
		for (D_mapint::iterator iter = nodes_level1.begin(); iter != nodes_level1.end(); ++iter)
		{
			if (cells_level1.find(iter->first) == cells_level1.end())
			{
				cells_level1.insert(make_pair(iter->first, std::array<int, 2>{iter->second, 0}));
			}
		}
	}

	// find a node inside each encolsed geometry and do flood fill
	for (unsigned int ishape = 0; ishape < Solid_Manager::pointer_me->numb_solids; ++ishape)
	{
		if (Solid_Manager::pointer_me->shape_solids.at(ishape).bool_enclosed == true)
		{
			bool_any_enclosed = true;
			// find a node inside the solid boundary
			D_real x_index = (Solid_Manager::pointer_me->shape_solids.at(ishape).x0) / gr_inner.dx;
			D_real y_index = (Solid_Manager::pointer_me->shape_solids.at(ishape).y0) / gr_inner.dx;
			D_uint xint = static_cast<D_uint>(x_index);
			D_uint yint = static_cast<D_uint>(y_index);
			D_uint xint_max = static_cast<D_uint>(xb_domain / gr_inner.dx + C_eps);
			D_uint yint_max = static_cast<D_uint>(yb_domain / gr_inner.dx + C_eps);
			D_morton morton_center;
			D_morton morton_temp = 0, keyin = 0;
			bool bool_key_Nfound = true;
			D_uint icount = 0;
#if (C_DIMS==2)
			morton_center = Morton_Assist::pointer_me->morton_encode(xint, yint);
#endif
#if (C_DIMS==3)
			D_real z_index = (Solid_Manager::pointer_me->shape_solids.at(ishape).z0) / gr_inner.dx;
			D_uint zint = static_cast<D_uint>(z_index);
			D_uint zint_max = static_cast<D_uint>(zb_domain / gr_inner.dx + C_eps);
			morton_center = Morton_Assist::pointer_me->morton_encode(xint, yint, zint);
#endif
			// search in -x direction from the solid center untill meet the first point in refine_node
			morton_temp = morton_center;
			icount = 0;
			while (refine_nodes.find(morton_temp) == refine_nodes.end())
			{
				morton_temp = Morton_Assist::find_x0(morton_temp, ilevel);
				if ((refine_nodes.find(morton_temp) != refine_nodes.end()) && (refine_nodes[morton_temp].at(0) == flag_refine))
				//if (refine_nodes.find(morton_temp) != refine_nodes.end())
				{
					keyin = morton_temp;
					bool_key_Nfound = false;
					break;
				}
				++icount;
				if (icount > xint)
				{
					break;
				}
			}


			// search in -y direction from the solid center untill meet the first point in refine_node
			if (bool_key_Nfound)
			{
				morton_temp = morton_center;
				icount = 0;
				while (refine_nodes.find(morton_temp) == refine_nodes.end())
				{
					morton_temp = Morton_Assist::find_y0(morton_temp, ilevel);
					if ((refine_nodes.find(morton_temp) != refine_nodes.end()) && (refine_nodes[morton_temp].at(0) == flag_refine))
					//if (refine_nodes.find(morton_temp) != refine_nodes.end())
					{
						keyin = morton_temp;
						bool_key_Nfound = false;
						break;
					}
					++icount;
					if (icount > yint)
					{
						break;
					}
				}
			}

			if (bool_key_Nfound)
			{
				std::stringstream error;
				error << "can't find a node inside the solid boundary" << std::endl;
				log_error(error.str(), Log_function::logfile);
			}
			else
			{
				std::stringstream infor;
				D_real x, y;
#if(C_DIMS == 2)
				Morton_Assist::compute_coordinate(keyin, ilevel, x, y);
				infor << "seed node inside the solid boundary for flood fill: x = " << x << ", y=" << y << std::endl;
#endif
#if(C_DIMS == 3)
				D_real z;
				Morton_Assist::compute_coordinate(keyin, ilevel, x, y, z);
				infor << "seed node inside the solid boundary for flood fill: x = " << x << ", y=" << y << ", z=" << z << std::endl;
#endif
				log_infor(infor.str(), Log_function::logfile);
			}

			// mark fluid nodes inside the solid 
			flood_fill_inner(ilevel, keyin, cells_level1, refine_nodes, boundary_x_temp, boundary_y_temp, boundary_z_temp);
		}
	}

	// find the outer layer
	D_mapint map_outer = {};
	search_outer_boundary(ilevel, cells_level1, cells_level1, refine_nodes, map_outer, boundary_x_temp, boundary_y_temp, boundary_z_temp);


	// add boundary nodes to refinenodes
	for (D_mapint::iterator iter = gr_inner.fine2coarse.at(iboundary1).begin(); iter != gr_inner.fine2coarse.at(iboundary1).end(); ++iter)
	{
		if (refine_nodes.find(iter->first) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(iter->first, std::array<int, 2>{flag0, 0}));
		}
		else
		{
			refine_nodes[iter->first].at(0) = flag0;
		}
	}

	// delete nodes on overlaping layers
	check_merge_nodes(ilevel, refine_nodes, map_outer);
	gr_NoIB[3].coarse2fine.at(iboundary2_coarse) = map_outer;

	//// find iboundary2
	search_inter_boundary(ilevel, iboundary1, iboundary2, iboundary2_coarse, map_outer, refine_nodes, boundary_x_temp, boundary_y_temp, boundary_z_temp);

	// add boundary nodes to refinenodes
	for (int iboundary = 1; iboundary < C_overlap * 2 + 1; ++iboundary)
	{
		for (D_mapint::iterator iter = gr_inner.fine2coarse.at(iboundary).begin(); iter != gr_inner.fine2coarse.at(iboundary).end(); ++iter)
		{
			if (refine_nodes.find(iter->first) == refine_nodes.end())
			{
				refine_nodes.insert(make_pair(iter->first, std::array<int, 2>{flag_iboundary[iboundary], 0}));
			}
			else
			{
				refine_nodes[iter->first].at(0) = flag_iboundary[iboundary];
			}
		}
	}

	// generate nodes
	for (D_mapint2::iterator iter = refine_nodes.begin(); iter != refine_nodes.end(); ++iter)
	{

		Node_IB temp;
		 temp.flag = iter->second.at(0);
		 temp.icount_ghost = iter->second.at(1);
		 gr_inner.grid.insert({ iter->first, temp });
	}

}


/**
* @brief      function to search nodes (C_max_level - 1) nearest to the solid node.
* @param[in]  ilevel    refinement level. 
* @param[in]  xyz       coordinate of the solid node.
* @param[out] nearest_center      nodes (C_max_level - 1) nearest to the solid node.
* @param[out] refine_nodes        nodes (C_max_level) near the solid.
* @note Here param[in]:ilevel is C_max_level.
*/
void Grid_Manager::find_nodes_near_solid(unsigned int ilevel, std::array<D_real, C_DIMS> xyz,D_mapint2 &refine_nodes)
{
	Grid * grid_ptr;
	if (ilevel == C_max_level)
	{
		grid_ptr = &gr_inner;
	}
	else
	{
		grid_ptr = &gr_NoIB.at(ilevel);
	}

	D_real x_index, y_index;
	D_uint xint, yint;
#if (C_DIMS==3)
	D_real  z_index;
	D_uint zint;
#endif

	// fluid node nearest to the solid node
	x_index = xyz.at(0) / grid_ptr->dx;
	y_index = xyz.at(1) / grid_ptr->dx;

#if (C_CHECK_MORTON_BOUNDARY == 1)
	if (xyz.at(0) >= xb_domain)
	{
		// to ensure the solid points located in a cell in the computational domain
		x_index = xb_domain / grid_ptr->dx - C_eps;
		//std::stringstream warning;
		//warning << "the x cooridnat solid node is " << xyz.at(0) << ": x boundary of the background mesh is" << C_xb << std::endl;
		//log_warning(warning.str(), Log_function::logfile);
}
	if (xyz.at(1) >= yb_domain )
	{
		y_index = yb_domain  / grid_ptr->dx - C_eps;
		//std::stringstream warning;
		//warning << "the y cooridnat solid node is " << xyz.at(1) << ": y boundary of the background mesh is" << yb_domain  << std::endl;
		//log_warning(warning.str(), Log_function::logfile);
	}
#if(C_DIMS == 3)
	if (xyz.at(2) >= zb_domain )
	{
		z_index = zb_domain  / grid_ptr->dx - C_eps;
		//std::stringstream warning;
		//warning << "the z cooridnat solid node is " << xyz.at(0) << ": z boundary of the background mesh is" << zb_domain  << std::endl;
		//log_warning(warning.str(), Log_function::logfile);
	}
#endif
	if (xyz.at(0) < 0)
	{
		std::stringstream error;
		error << "the x cooridnate of solid node is less than 0" << std::endl;
		log_error(error.str(), Log_function::logfile);
	}
	if (xyz.at(1) < 0)
	{
		std::stringstream error;
		error << "the y cooridnate of solid node is less than 0" << std::endl;
		log_error(error.str(), Log_function::logfile);
	}
#if(C_DIMS == 3)
	if (xyz.at(2) < 0)
	{
		std::stringstream error;
		error << "the z cooridnate of solid node is less than 0" << std::endl;
		log_error(error.str(), Log_function::logfile);
	}
#endif
#endif
	xint = static_cast<D_uint>(x_index + C_eps);
	yint = static_cast<D_uint>(y_index + C_eps);
#if (C_DIMS==3)
	z_index = xyz.at(2) / gr_inner.dx;
	zint = static_cast<D_uint>(z_index + C_eps);
#endif

	std::vector<D_morton> vector_key(two_power_n(C_DIMS));
#if (C_DIMS == 2)
	vector_key.at(0) = Morton_Assist::pointer_me->morton_encode(xint, yint);
#endif
#if (C_DIMS == 3)
	vector_key.at(0) = Morton_Assist::pointer_me->morton_encode(xint, yint, zint);
#endif
	D_real x, y, z;
	Morton_Assist::compute_coordinate(vector_key.at(0), ilevel, x, y, z);
	vector_key.at(1) = Morton_Assist::find_x1(vector_key.at(0), ilevel);
	vector_key.at(2) = Morton_Assist::find_y1(vector_key.at(1), ilevel);
	vector_key.at(3) = Morton_Assist::find_y1(vector_key.at(0), ilevel);
#if (C_DIMS == 3)
	vector_key.at(4) = Morton_Assist::find_z1(vector_key.at(0), ilevel);
	vector_key.at(5) = Morton_Assist::find_z1(vector_key.at(1), ilevel);
	vector_key.at(6) = Morton_Assist::find_z1(vector_key.at(2), ilevel);
	vector_key.at(7) = Morton_Assist::find_z1(vector_key.at(3), ilevel);
#endif

	for (int ikey = 0; ikey < vector_key.size(); ++ikey)
	{
		if (refine_nodes.find(vector_key.at(ikey)) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(vector_key.at(ikey), std::array<int, 2>{flag_near_solid, 0}));
		}
		else
		{
			refine_nodes[vector_key.at(ikey)].at(0) = flag_near_solid;
		}

		Infor_near_solid node_temp;
		if (map_node_near_solid.find(vector_key.at(ikey)) == map_node_near_solid.end())
		{
			node_temp.icount = 1;
			map_node_near_solid.insert({ vector_key.at(ikey), node_temp });
		}
		else
		{
			++map_node_near_solid[vector_key.at(ikey)].icount;
		}
	}

}

#if	(C_FSI_INTERFACE == 1)
void Grid_Manager::initial_map_node_IB(unsigned int ilevel)
{
	Infor_IB node_temp;
	node_temp.icount = 1;
	// search for nodes used for IBM
	if (C_extend_IB == 1)
	{
		for (D_map_Infor_near_solid::iterator iter = map_node_near_solid.begin(); iter != map_node_near_solid.end(); ++iter)
		{
			if (map_node_IB.find(iter->first) == map_node_IB.end())
			{
				map_node_IB.insert({ iter->first, node_temp });
			}
		}
	}
	else if (C_extend_IB >1)
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

		// compute icount of nodes within a given region
		for (D_map_Infor_near_solid::iterator iter = map_node_near_solid.begin(); iter != map_node_near_solid.end(); ++iter)
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
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0;
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0);
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
	}
}
#endif

/**
* @brief      mark nodes within C_extend_ghost grid space extended from solid points
* @param[in]  ilevel    refinement level.
* @param[out] refine_nodes        nodes (C_max_level) near the solid.
* @note Here param[in]:ilevel is C_max_level.
*/
void Grid_Manager::mark_ghost_points(unsigned int ilevel, D_mapint2& refine_nodes)
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
	if (C_extend_inner <= C_extend_ghost)
	{
		std::stringstream error;
		error << "C_extend_inner should be greater than C_extend_ghost" << std::endl;
		log_error(error.str(), Log_function::logfile);
	}

	unsigned int extend_temp_x0 = extend_x0, extend_temp_x1 = extend_x1, extend_temp_y0 = extend_y0, extend_temp_y1 = extend_y1;
#if (C_DIMS==3)
	unsigned int extend_temp_z0= extend_z0, extend_temp_z1 = extend_z1;
#endif

	for (D_map_Infor_near_solid::iterator iter = map_node_near_solid.begin(); iter != map_node_near_solid.end(); ++iter)
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
		// find out if there are nodes on the domain boundary
		if (extend_temp_x0 != extend_temp_x0)
		{
#if (C_DIMS == 2)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0);
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0);
			morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
			for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
			{
				morton_temp1 = morton_temp2;
#endif
				for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
				{
					if (bk_boundary_x.at(ilevel).at(0).find(morton_temp1) == bk_boundary_x.at(ilevel).at(0).end())
					{
						bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel);
				}
#if (C_DIMS == 3)
				morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
			}
		}
#endif
		if (extend_temp_x1 != extend_temp_x1)
		{
#if (C_DIMS == 2)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) + extend_temp_x1, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0);
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) + extend_temp_x1, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0);
			morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
			for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
			{
				morton_temp1 = morton_temp2;
#endif
				for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
				{
					if (bk_boundary_x.at(ilevel).at(1).find(morton_temp1) == bk_boundary_x.at(ilevel).at(1).end())
					{
						bk_boundary_x.at(ilevel)[1].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel);
				}
#if (C_DIMS == 3)
				morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
			}
		}
#endif
		if (extend_temp_y0 != extend_temp_y0)
		{
#if (C_DIMS == 2)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0);
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0);
			morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
			for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
			{
				morton_temp1 = morton_temp2;
#endif
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (bk_boundary_y.at(ilevel).at(0).find(morton_temp1) == bk_boundary_y.at(ilevel).at(0).end())
					{
						bk_boundary_y.at(ilevel)[0].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_x1(morton_temp1, ilevel);
				}
#if (C_DIMS == 3)
				morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
			}
	}
#endif
		if (extend_temp_y1 != extend_temp_y1)
		{
#if (C_DIMS == 2)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) + extend_temp_y1);
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) + extend_temp_y1, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0);
			morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
			for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
			{
				morton_temp1 = morton_temp2;
#endif
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (bk_boundary_y.at(ilevel).at(1).find(morton_temp1) == bk_boundary_y.at(ilevel).at(1).end())
					{
						bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_x1(morton_temp1, ilevel);
				}
#if (C_DIMS == 3)
				morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
			}
	}
#endif
#if (C_DIMS == 3)
		if (extend_temp_z0 != extend_temp_z0)
		{
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0);
			morton_temp2 = morton_temp;
			for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
			{
				morton_temp1 = morton_temp2;
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (bk_boundary_z.at(ilevel).at(0).find(morton_temp1) == bk_boundary_z.at(ilevel).at(0).end())
					{
						bk_boundary_z.at(ilevel)[0].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_x1(morton_temp1, ilevel);
				}
				morton_temp2 = Morton_Assist::find_y1(morton_temp2, ilevel);
			}
		}
		if (extend_temp_z1 != extend_temp_z1)
		{
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0, static_cast<D_int>(z / gr_inner.dx + C_eps) + extend_temp_z1);
			morton_temp2 = morton_temp;
			for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
			{
				morton_temp1 = morton_temp2;
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (bk_boundary_z.at(ilevel).at(1).find(morton_temp1) == bk_boundary_z.at(ilevel).at(1).end())
					{
						bk_boundary_z.at(ilevel)[1].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_x1(morton_temp1, ilevel);
				}
				morton_temp2 = Morton_Assist::find_y1(morton_temp2, ilevel);
			}
		}
#endif
#endif
		// seach nodes in the seraching area
#if (C_DIMS == 2)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0);
		morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0);
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
					if (refine_nodes.find(morton_temp) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp, std::array<int, 2>{flag_ghost_temp, 1}));
					}
					else
					{
						++refine_nodes[morton_temp].at(1);
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
* @brief      function to search nodes (C_max_level - 1) within a predefined distance (C_extend_inner/2) from solid nodes.
* @param[in]  ilevel    refinement level.
* @param[in]  nearest_center      nodes (C_max_level - 1) nearest to the solid node.
* @param[out] nodes_level1        nodes (C_max_level - 1) within the predefined distance from solid nodes.
* @param[out] refine_nodes        nodes (C_max_level) near the solid.
* @note Here param[in]:ilevel is C_max_level.
*/
void Grid_Manager::search_nodes_near_solid(unsigned int ilevel, D_mapint &nearest_center, D_mapint &nodes_level1, D_mapint2 &refine_nodes, std::array<D_mapint, 2> &boundary_x_temp, std::array<D_mapint, 2> &boundary_y_temp, std::array<D_mapint, 2> &boundary_z_temp)
{
	Grid * grid_ptr;
	if (ilevel == C_max_level)
	{
		grid_ptr = &gr_inner;
	}
	else
	{
		grid_ptr = &gr_NoIB.at(ilevel);
	}

	unsigned int extend_x0 = C_extend_inner / 2 + C_extend_inner_x0 / 2;
	unsigned int extend_x1 = C_extend_inner / 2 + C_extend_inner_x1 / 2;
	unsigned int extend_y0 = C_extend_inner / 2 + C_extend_inner_y0 / 2;
	unsigned int extend_y1 = C_extend_inner / 2 + C_extend_inner_y1 / 2;
#if (C_DIMS == 3)
	unsigned int extend_z0 = C_extend_inner / 2 + C_extend_inner_z0 / 2;
	unsigned int extend_z1 = C_extend_inner / 2 + C_extend_inner_z1 / 2;
#endif

	nodes_level1 = nearest_center;
	const unsigned int ilevel_1 = ilevel - 1;
	D_morton morton_temp, morton_temp0, morton_temp1;
#if (C_DIMS==3)
	D_morton morton_temp2;
#endif
	unsigned int extend_temp_x0 = extend_x0, extend_temp_x1 = extend_x1, extend_temp_y0 = extend_y0, extend_temp_y1 = extend_y1;
#if (C_DIMS==3)
	unsigned int extend_temp_z0 = extend_z0, extend_temp_z1 = extend_z1;
#endif
	// nodes at ilevel near the solid node
	for (D_mapint::iterator iter = nearest_center.begin(); iter != nearest_center.end(); ++iter)
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
		// find out if there are nodes on the domain boundary
		if (extend_temp_x0 != extend_temp_x0)
		{
#if (C_DIMS == 2)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0 * 2);
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0 * 2);
			morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
			for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
			{
				morton_temp1 = morton_temp2;
#endif
				for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
				{
					if (bk_boundary_x.at(ilevel).at(0).find(morton_temp1) == bk_boundary_x.at(ilevel).at(0).end())
					{
						bk_boundary_x.at(ilevel)[0].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel);
				}
#if (C_DIMS == 3)
				morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
			}
	}
#endif
		if (extend_temp_x1 != extend_temp_x1)
		{
#if (C_DIMS == 2)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) + extend_temp_x1 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0 * 2;
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) + extend_temp_x1 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0 * 2);
			morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
			for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
			{
				morton_temp1 = morton_temp2;
#endif
				for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
				{
					if (bk_boundary_x.at(ilevel).at(1).find(morton_temp1) == bk_boundary_x.at(ilevel).at(1).end())
					{
						bk_boundary_x.at(ilevel)[1].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel);
				}
#if (C_DIMS == 3)
				morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
			}
		}
#endif
		if (extend_temp_y0 != extend_temp_y0)
		{
#if (C_DIMS == 2)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0 * 2);
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0 * 2);
			morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
			for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
			{
				morton_temp1 = morton_temp2;
#endif
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (bk_boundary_y.at(ilevel).at(0).find(morton_temp1) == bk_boundary_y.at(ilevel).at(0).end())
					{
						bk_boundary_y.at(ilevel)[0].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_x1(morton_temp1, ilevel);
				}
#if (C_DIMS == 3)
				morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
			}
		}
#endif
		if (extend_temp_y1 != extend_temp_y1)
		{
#if (C_DIMS == 2)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) + extend_temp_y1 * 2);
			morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) + extend_temp_y1 * 2, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0 * 2);
			morton_temp2 = morton_temp;
#endif
#if (C_DIMS == 3)
			for (unsigned int iz = 0; iz <= extend_temp_z0 + extend_temp_z1; ++iz)
			{
				morton_temp1 = morton_temp2;
#endif
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (bk_boundary_y.at(ilevel).at(1).find(morton_temp1) == bk_boundary_y.at(ilevel).at(1).end())
					{
						bk_boundary_y.at(ilevel)[1].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_x1(morton_temp1, ilevel);
				}
#if (C_DIMS == 3)
				morton_temp2 = Morton_Assist::find_z1(morton_temp2, ilevel);
			}
		}
#endif
#if (C_DIMS == 3)
		if (extend_temp_z0 != extend_temp_z0)
		{
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0 * 2);
			morton_temp2 = morton_temp;
			for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
			{
				morton_temp1 = morton_temp2;
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (bk_boundary_z.at(ilevel).at(0).find(morton_temp1) == bk_boundary_z.at(ilevel).at(0).end())
					{
						bk_boundary_z.at(ilevel)[0].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_x1(morton_temp1, ilevel);
				}
				morton_temp2 = Morton_Assist::find_y1(morton_temp2, ilevel);
			}
		}
		if (extend_temp_z1 != extend_temp_z1)
		{
			morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / gr_inner.dx + C_eps) + extend_temp_z1 * 2);
			morton_temp2 = morton_temp;
			for (unsigned int iy = 0; iy <= extend_temp_y0 + extend_temp_y1; ++iy)
			{
				morton_temp1 = morton_temp2;
				for (unsigned int ix = 0; ix <= extend_temp_x0 + extend_temp_x1; ++ix)
				{
					if (bk_boundary_z.at(ilevel).at(1).find(morton_temp1) == bk_boundary_z.at(ilevel).at(1).end())
					{
						bk_boundary_z.at(ilevel)[1].insert({ morton_temp1, nullptr });
					}
					morton_temp1 = Morton_Assist::find_x1(morton_temp1, ilevel);
				}
				morton_temp2 = Morton_Assist::find_y1(morton_temp2, ilevel);
			}
		}
#endif
#endif
#if (C_DIMS == 2)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0 * 2;
		morton_temp1 = morton_temp;
#endif
#if (C_DIMS == 3)
		morton_temp = Morton_Assist::pointer_me->morton_encode(static_cast<D_int>(x / gr_inner.dx + C_eps) - extend_temp_x0 * 2, static_cast<D_int>(y / gr_inner.dx + C_eps) - extend_temp_y0 * 2, static_cast<D_int>(z / gr_inner.dx + C_eps) - extend_temp_z0 * 2);
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
					if (nodes_level1.find(morton_temp) == nodes_level1.end())
					{
						nodes_level1.insert(make_pair(morton_temp, flag_refine));
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

	D_real bx0_temp = Solid_Manager::pointer_me->shape_offest_x0_grid + C_eps;
	D_real by0_temp = Solid_Manager::pointer_me->shape_offest_y0_grid + C_eps;
#if(C_DIMS == 3)
	D_real bz0_temp = Solid_Manager::pointer_me->shape_offest_z0_grid + C_eps;
#endif
	// generate nodes at ilevel in refine_nodes based on those at ilevel - 1 in nodes_level1 
	for (D_mapint::iterator iter = nodes_level1.begin(); iter != nodes_level1.end(); ++iter)
	{
		D_morton morton_temp0, morton_temp1; // 0: faces; 1: edges; 2: corners
		D_morton morton_tempx, morton_tempy;
#if(C_DIMS==3)
		D_morton morton_temp2, morton_tempz;
#endif
		// cell center
		D_morton morton_temp_center = iter->first;
		if (refine_nodes.find(morton_temp_center) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp_center, std::array<int, 2> { flag_refine, 0 }));
		}

#if (C_CHECK_MORTON_BOUNDARY == 1)
#if (C_DIMS == 2)
		D_real x, y;
		Morton_Assist::pointer_me->compute_coordinate(morton_temp_center, ilevel, x, y);
#endif
#if (C_DIMS == 3)
		D_real x, y, z;
		Morton_Assist::pointer_me->compute_coordinate(morton_temp_center, ilevel, x, y, z);
#endif
		if (x < bx0_temp)
		{
			morton_temp_center = Morton_Assist::find_x1(morton_temp_center, ilevel);
		}
#if (C_SPEED_UP == 0)
		else if ((xb_domain - x) < C_eps)
		{
			morton_temp_center = Morton_Assist::find_x0(morton_temp_center, ilevel);
		}
#endif
		if (y < by0_temp)
		{
			morton_temp_center = Morton_Assist::find_y1(morton_temp_center, ilevel);
		}
#if (C_SPEED_UP == 0)
		else if ((yb_domain - y) < C_eps)
		{
			morton_temp_center = Morton_Assist::find_y0(morton_temp_center, ilevel);
		}
#endif
#if (C_DIMS == 3)
		if (z < by0_temp)
		{
			morton_temp_center = Morton_Assist::find_z1(morton_temp_center, ilevel);
	    }
#if (C_SPEED_UP == 0)
		else if ((zb_domain - z) < C_eps)
		{
			morton_temp_center = Morton_Assist::find_z0(morton_temp_center, ilevel);
		}
#endif
#endif
#endif
		// (+x) node in the cell
		morton_temp0 = Morton_Assist::find_x1(morton_temp_center, ilevel);
		if (refine_nodes.find(morton_temp0) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
		}

		// (+x, +y) node in the cell
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}

#if(C_DIMS==3)
		// (+x, +y, +z) node in the cell
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (refine_nodes.find(morton_temp2) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
		}
#if(C_SPEED_UP==0)
		// (+x, +y, - z) node in the cell
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (refine_nodes.find(morton_temp2) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}
#endif
#endif
#if(C_SPEED_UP==0)
		// (+x, -y) node in the cell
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}
#if(C_DIMS==3)
		// (+x, -y, -z) node in the cell
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (refine_nodes.find(morton_temp2) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
		}
		// (+x, -y, +z) node in the cell
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (refine_nodes.find(morton_temp2) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
		}
#endif
#endif
#if(C_DIMS==3)
		// (+x, +z) node in the cell
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}
#if(C_SPEED_UP==0)
		// (+x, -z) node in the cell
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}
#endif
#endif

#if(C_SPEED_UP==0)
		// (-x) node in the cell
		morton_temp0 = Morton_Assist::find_x0(morton_temp_center, ilevel);
		if (refine_nodes.find(morton_temp0) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
		}
		// (-x, -y) node in the cell
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}
#if(C_DIMS==3)
		// (-x, -y, -z) node in the cell
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (refine_nodes.find(morton_temp2) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
		}
		// (-x, -y, +z) node in the cell
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (refine_nodes.find(morton_temp2) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
		}
		// (-x, -z) node in the cell
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}
		// (-x, +z) node in the cell
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}

#endif
		// (-x, +y) node in the cell
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}
#if(C_DIMS==3)
		// (-x, +y, -z) node in the cell
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (refine_nodes.find(morton_temp2) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
		}
		// (+x, +y, +z) node in the cell
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (refine_nodes.find(morton_temp2) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
		}
#endif
#endif

		// (+y) node in the cell
		morton_temp0 = Morton_Assist::find_y1(morton_temp_center, ilevel);
		if (refine_nodes.find(morton_temp0) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
		}
#if(C_DIMS==3)
		// (+y, +z) node in the cell
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0}));
		}
#if(C_SPEED_UP==0)
		// (+y, -z) node in the cell
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}
#endif
#endif

#if(C_SPEED_UP==0)
		// (-y) node in the cell
		morton_temp0 = Morton_Assist::find_y0(morton_temp_center, ilevel);
		if (refine_nodes.find(morton_temp0) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
		}
#if(C_DIMS==3)
		// (-y, -z) node in the cell
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}
		// (-y, +z) node in the cell
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (refine_nodes.find(morton_temp1) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
		}
#endif
#endif
#if(C_DIMS==3)
		// (+z) node in the cell
		morton_temp0 = Morton_Assist::find_z1(morton_temp_center, ilevel);
		if (refine_nodes.find(morton_temp0) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
		}

#if(C_SPEED_UP==0)
		// (-z) node in the cell
		morton_temp0 = Morton_Assist::find_z1(morton_temp_center, ilevel);
		if (refine_nodes.find(morton_temp0) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
		}
#endif
#endif

	}
}


/**
* @brief      function to delete fluid nodes inside the solid.
* @param[in]  ilevel    refinement level.
* @param[in]  key_in    a morton key representing a node inside in the enclosed solid.
* @param[in]  nodes_level1        nodes (C_max_level - 1) within the predefined distance from solid nodes.
* @param[out] refine_nodes        nodes (C_max_level) near the solid.
* @note This function will be called only if shape.enclosed == true.
*/
void Grid_Manager::flood_fill_inner(unsigned int ilevel, D_morton &key_in, D_mapint2 &nodes_level1, D_mapint2 &refine_nodes, std::array<D_mapint, 2> &boundary_x_temp, std::array<D_mapint, 2> &boundary_y_temp, std::array<D_mapint, 2> &boundary_z_temp)
{
	D_morton morton_xyz = Morton_Assist::morton_xyz.at(ilevel);
	std::vector<D_morton> stk;
	stk.push_back(key_in);
	while (!stk.empty())
	{
		D_morton seed = stk.back();
		D_morton morton_temp;
		stk.pop_back();

		morton_temp = seed;

		if ((refine_nodes.find(morton_temp) == refine_nodes.end()) || (refine_nodes[morton_temp].at(0) == flag_near_solid) || (refine_nodes[morton_temp].at(0) == flag_ghost))
		{
			//if (refine_nodes.find(morton_temp) != refine_nodes.end())
			//{
			//	//if (refine_nodes[morton_temp].at(0) == flag_near_solid)
			//	//{
			//	//	refine_nodes[morton_temp].at(0) = flag3;
			//	//}
			//}
			//continue;

		}
		else
		{
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if ((morton_temp&Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
			{
				if (boundary_x_temp[0].find(morton_temp) == boundary_x_temp[0].end())
				{
					boundary_x_temp[0].insert({ morton_temp, flag_refine });
				}
			}
			else if ((morton_temp&Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
			{
				if (boundary_x_temp[1].find(morton_temp) == boundary_x_temp[1].end())
				{
					boundary_x_temp[1].insert({ morton_temp, flag_refine });
				}
			}
			if ((morton_temp&Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
			{
				if (boundary_y_temp[0].find(morton_temp) == boundary_y_temp[0].end())
				{
					boundary_y_temp[0].insert({ morton_temp, flag_refine });
				}
			}
			else if ((morton_temp&Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
			{
				if (boundary_y_temp[1].find(morton_temp) == boundary_y_temp[1].end())
				{
					boundary_y_temp[1].insert({ morton_temp, flag_refine });
				}
			}
#if(C_DIMS == 3)
			if ((morton_temp&Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
			{
				if (boundary_z_temp[0].find(morton_temp) == boundary_z_temp[0].end())
				{
					boundary_z_temp[0].insert({ morton_temp, flag_refine });
				}
			}
			else if ((morton_temp&Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
			{
				if (boundary_z_temp[1].find(morton_temp) == boundary_z_temp[1].end())
				{
					boundary_z_temp[1].insert({ morton_temp, flag_refine });
				}
			}
#endif
#endif
			if (refine_nodes[morton_temp].at(0) == flag_ghost_temp)
			{
				refine_nodes[morton_temp].at(0) = flag_ghost;
			}
			else
			{
				refine_nodes.erase(morton_temp);
			}

			if ((morton_temp & morton_xyz)==0)
			{
				if (nodes_level1.find(morton_temp) != nodes_level1.end()) // colour the node at ilevel - 1 inside the solid boundary
				{
					nodes_level1[morton_temp].at(0) = flag_inner;
				}
			}

			morton_temp = Morton_Assist::find_x0(seed, ilevel);
			stk.push_back(morton_temp);
			morton_temp = Morton_Assist::find_x1(seed, ilevel);
			stk.push_back(morton_temp);
			morton_temp = Morton_Assist::find_y0(seed, ilevel);
			stk.push_back(morton_temp);
			morton_temp = Morton_Assist::find_y1(seed, ilevel);
			stk.push_back(morton_temp);
#if (C_DIMS==3)
			morton_temp = Morton_Assist::find_z0(seed, ilevel);
			stk.push_back(morton_temp);
			morton_temp = Morton_Assist::find_z1(seed, ilevel);
			stk.push_back(morton_temp);
#endif
		}
	}
}


/**
* @brief      function to search the innermost numerical boundary.
* @param[in]  ilevel    refinement level.
* @param[in]  iter_nodes               nodes used for searching.
* @param[in]  exist_nodes              exsiting nodes.
* @param[out] refine_nodes_temp        mesh where nodes on the numerical boudnary are added.
* @param[out] map_outer                nodes on the the outer layer.
* @param[out] boundary_x_temp          tempory array to store domain boundaries when x = xmin and x = xmax (only works when C_CHECK_MORTON_BOUNDARY = 1).
* @param[out] boundary_y_temp          tempory array to store domain boundaries when y = ymin and y = ymax (only works when C_CHECK_MORTON_BOUNDARY = 1).
* @param[out] boundary_z_temp          tempory array to store domain boundaries when z = xmin and z = zmax (only works when C_CHECK_MORTON_BOUNDARY = 1).
* @note  Nodes on the numerical boundaries will be used to transfer informtion between different refinement level.
*/
void Grid_Manager::search_outer_boundary(unsigned int ilevel, D_mapint2 &iter_nodes, D_mapint2 &exist_nodes, D_mapint2 &refine_nodes_temp, D_mapint &map_outer, std::array<D_mapint, 2> &boundary_x_temp, std::array<D_mapint, 2> &boundary_y_temp, std::array<D_mapint, 2> &boundary_z_temp)
{
	unsigned int ilevel_1 = ilevel - 1;
	Grid * grid_ptr;
	if (ilevel == C_max_level)
	{
		grid_ptr = &(Grid_Manager::pointer_me->gr_inner);
	}
	else
	{
		grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel));
	}

	Grid *grid_ptr_coarse = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel - 1));

	D_morton morton_temp0, morton_temp1;
	D_morton morton_temp0_ilevel, morton_temp1_ilevel;
	D_morton morton_tempx, morton_tempy;
#if (C_DIMS==3)
	D_morton morton_temp2, morton_temp2_ilevel;
	D_morton morton_tempz;
#endif
	for (D_mapint2::iterator iter = iter_nodes.begin(); iter != iter_nodes.end(); ++iter)
	{
		if (iter->second.at(0) != flag_inner)
		{
#if (C_CHECK_MORTON_BOUNDARY == 1)
			morton_tempx = iter->first & Morton_Assist::xmorton_ref.at(ilevel).at(1);
			morton_tempy = iter->first & Morton_Assist::ymorton_ref.at(ilevel).at(1);
#if (C_DIMS == 3)
			morton_tempz = iter->first & Morton_Assist::zmorton_ref.at(ilevel).at(1);
#endif
			// check if the current cell is not on the boundary of the compuational domain
			bool x_positive = true, x_negative = true, y_positive = true, y_negative = true;
			if ( morton_tempx == Morton_Assist::mortonx_min)
			{
				x_negative = false;
				if (boundary_x_temp[0].find(iter->first) == boundary_x_temp[0].end())
				{
					boundary_x_temp[0].insert({ iter->first, flag_refine });
				}
			}

			else if (morton_tempx == Morton_Assist::mortonx_max)
			{
				x_positive = false;
				if (boundary_x_temp[1].find(iter->first) == boundary_x_temp[1].end())
				{
					boundary_x_temp[1].insert({ iter->first, flag_refine });
				}
			}

			if (morton_tempy == Morton_Assist::mortony_min)
			{
				y_negative = false;
				if (boundary_y_temp[0].find(iter->first) == boundary_y_temp[0].end())
				{
					boundary_y_temp[0].insert({ iter->first, flag_refine });
				}
			}
			else if (morton_tempy == Morton_Assist::mortony_max )
			{
				y_positive = false;
				if (boundary_y_temp[1].find(iter->first) == boundary_y_temp[1].end())
				{
					boundary_y_temp[1].insert({ iter->first, flag_refine });
				}
			}

#if (C_DIMS == 3)
			bool z_positive = true, z_negative = true;
			if (morton_tempz == Morton_Assist::mortonz_min)
			{
				z_negative = false;
				if (boundary_z_temp[0].find(iter->first) == boundary_z_temp[0].end())
				{
					boundary_z_temp[0].insert({ iter->first, flag_refine });
				}
			}
			else if (morton_tempz == Morton_Assist::mortonz_max )
			{
				z_positive = false;
				if (boundary_z_temp[1].find(iter->first) == boundary_z_temp[1].end())
				{
					boundary_z_temp[1].insert({ iter->first, flag_refine });
				}
			}
#endif
#endif
			// check if cell at the (-x) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative)
			{
#endif
				morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel_1);
				morton_temp0_ilevel = Morton_Assist::find_x0(iter->first, ilevel);
				if (exist_nodes.find(morton_temp0) == exist_nodes.end())
				{
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
#if (C_DIMS == 3)
					morton_tempz = morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
#endif
					if (morton_tempx == Morton_Assist::mortonx_min) // if node (morton_temp0) is on the boundary of the computatioanl domain
					{
						if (boundary_x_temp[0].find(morton_temp0) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp0, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp0) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0, std::array<int, 2> { flag4, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp0_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp0) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp0) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#if(C_DIMS == 3)
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp0) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp0) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#endif

#endif
					// node at the outer interface of the boundary layer at ilevel_1
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_x1 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_x1;
					}

					// node at the middle interface of the boundary layer at ilevel
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}

				}
#if (C_CHECK_MORTON_BOUNDARY == 1)				
			}
#endif
			// check if cell at the (-x, -y) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel_1);
				morton_temp1_ilevel = Morton_Assist::find_y0(morton_temp0_ilevel, ilevel);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
#if (C_DIMS == 3)
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
#endif
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp1) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag4, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp1) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag6, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#if(C_DIMS == 3)
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp1) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp1) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}

				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (-x, -y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_negative && z_negative)
			{
#endif
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
				if (exist_nodes.find(morton_temp2) == exist_nodes.end())
				{
					morton_temp2_ilevel = Morton_Assist::find_z0(morton_temp1_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp2) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag4, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp2) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag6, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp2) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag8, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#endif
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-x, -y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_negative && z_positive)
			{
#endif
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
				if (exist_nodes.find(morton_temp2) == exist_nodes.end())
				{
					morton_temp2_ilevel = Morton_Assist::find_z1(morton_temp1_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp2) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag4, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp2) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag6, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp2) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag9, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#endif
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-x, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && z_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
					morton_temp1_ilevel = Morton_Assist::find_z0(morton_temp0_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp1) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag4, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp1) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag8, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}

					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp1) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp1) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp1, flag_refine });
						}
					}

#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{					
						map_outer.insert({ morton_temp1, 0 });
					}

					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel,flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-x, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && z_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
					morton_temp1_ilevel = Morton_Assist::find_z1(morton_temp0_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp1) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag4, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp1) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag9, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}

					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp1) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp1) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{					
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// check if cell at the (-x, +y) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel_1);
				morton_temp1_ilevel = Morton_Assist::find_y1(morton_temp0_ilevel, ilevel);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
#if (C_DIMS == 3)
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
#endif
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp1) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag4, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp1) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag7, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#if(C_DIMS == 3)
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp1) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp1) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (-x, +y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_positive && z_negative)
			{
#endif
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
				if (exist_nodes.find(morton_temp2) == exist_nodes.end())
				{
					morton_temp2_ilevel = Morton_Assist::find_z0(morton_temp1_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp2) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag4, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp2) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag7, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp2) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag8, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#endif
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-x, +y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_positive && z_positive)
			{
#endif
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
				if (exist_nodes.find(morton_temp2) == exist_nodes.end())
				{
					morton_temp2_ilevel = Morton_Assist::find_z1(morton_temp1_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp2) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag4, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp2) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag7, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp2) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag9, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#endif
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif

			// check if cell at the (+x) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive)
			{
#endif
				morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel_1);
				morton_temp0_ilevel = Morton_Assist::find_x1(iter->first, ilevel);
				if (exist_nodes.find(morton_temp0) == exist_nodes.end())
				{
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
#if (C_DIMS == 3)
					morton_tempz = morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
#endif
					if (morton_tempx == Morton_Assist::mortonx_max) // if node (morton_temp0) is on the boundary of the computatioanl domain
					{
						if (boundary_x_temp[1].find(morton_temp0) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp0, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp0) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0, std::array<int, 2> { flag5, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp0_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp0) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp0) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#if(C_DIMS == 3)
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp0) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp0) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#endif
#endif
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_x0 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_x0;
					}

					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first,flag_iboundary[0] });
						}
					}

				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+x, -y) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel_1);
				morton_temp1_ilevel = Morton_Assist::find_y0(morton_temp0_ilevel, ilevel);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
#if (C_DIMS == 3)
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
#endif
					if (morton_tempx == Morton_Assist::mortonx_max)
					{
						if (boundary_x_temp[1].find(morton_temp1) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag5, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp1) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag6, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#if(C_DIMS == 3)
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp1) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp1) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
#endif

					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}


					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (+x, -y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_negative && z_negative)
			{
#endif
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
				if (exist_nodes.find(morton_temp2) == exist_nodes.end())
				{
					morton_temp2_ilevel = Morton_Assist::find_z0(morton_temp1_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_max)
					{
						if (boundary_x_temp[1].find(morton_temp2) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag5, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp2) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag6, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp2) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag8, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#endif
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+x, -y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_negative && z_positive)
			{
#endif
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
				if (exist_nodes.find(morton_temp2) == exist_nodes.end())
				{
					morton_temp2_ilevel = Morton_Assist::find_z1(morton_temp1_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_max)
					{
						if (boundary_x_temp[1].find(morton_temp2) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag5, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp2) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag6, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp2) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag9, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#endif
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// check if cell at the (+x, +y) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel_1);
				morton_temp1_ilevel = Morton_Assist::find_y1(morton_temp0_ilevel, ilevel);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
#if (C_DIMS == 3)
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
#endif
					if (morton_tempx == Morton_Assist::mortonx_max)
					{
						if (boundary_x_temp[1].find(morton_temp1) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag5, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp1) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag7, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#if(C_DIMS == 3)
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp1) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp1) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (+x, +y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_positive & z_negative)
			{
#endif
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
				if (exist_nodes.find(morton_temp2) == exist_nodes.end())
				{
					morton_temp2_ilevel = Morton_Assist::find_z0(morton_temp1_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp2) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag5, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp2) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag7, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp2) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag8, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#endif
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+x, +y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_positive & z_positive)
			{
#endif
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
				if (exist_nodes.find(morton_temp2) == exist_nodes.end())
				{
					morton_temp2_ilevel = Morton_Assist::find_z1(morton_temp1_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp2 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp2 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp2 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp2) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag5, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp2) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag7, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp2) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp2, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp2) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2, std::array<int, 2> { flag9, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp2_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp2_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
#endif
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+x, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive & z_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
					morton_temp1_ilevel = Morton_Assist::find_z0(morton_temp0_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp1) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag5, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp1) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag8, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp1) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp1) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+x, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive & z_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
					morton_temp1_ilevel = Morton_Assist::find_z1(morton_temp0_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp1) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag5, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp1) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag9, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp1) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp1) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif

			// check if cell at the (-y) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_negative)
			{
#endif
				morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel_1);
				morton_temp0_ilevel = Morton_Assist::find_y0(iter->first, ilevel);
				if (exist_nodes.find(morton_temp0) == exist_nodes.end())
				{
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
#if (C_DIMS == 3)
					morton_tempz = morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
#endif
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp0) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp0, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp0) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0, std::array<int, 2> { flag6, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp0_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp0) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp0) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#if(C_DIMS == 3)
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp0) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp0) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#endif
#endif
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_y1 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_y1;
					}

					

					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}

				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (-y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_negative && z_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
					morton_temp1_ilevel = Morton_Assist::find_z0(morton_temp0_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp1) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag6, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp1) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag8, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp1) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp1) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_negative && z_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
					morton_temp1_ilevel = Morton_Assist::find_z1(morton_temp0_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp1) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag6, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp1) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag9, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp1) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp1) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// check if cell at the (+y) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_positive)
			{
#endif
				morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel_1);
				morton_temp0_ilevel = Morton_Assist::find_y1(iter->first, ilevel);
				if (exist_nodes.find(morton_temp0) == exist_nodes.end())
				{
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
#if (C_DIMS == 3)
					morton_tempz = morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
#endif
					if (morton_tempy == Morton_Assist::mortony_max ) // if node (morton_temp0) is on the boundary of the computatioanl domain
					{
						if (boundary_y_temp[1].find(morton_temp0) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp0, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp0) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0, std::array<int, 2> { flag7, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp0_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp0) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp0) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#if(C_DIMS == 3)
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp0) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp0) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#endif
#endif
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_y0 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_y0;
					}

					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (+y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_positive && z_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
					morton_temp1_ilevel = Morton_Assist::find_z0(morton_temp0_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp1) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag7, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_min)
					{
						if (boundary_z_temp[0].find(morton_temp1) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag8, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp1) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp1) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_positive && z_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
				if (exist_nodes.find(morton_temp1) == exist_nodes.end())
				{
					morton_temp1_ilevel = Morton_Assist::find_z1(morton_temp0_ilevel, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp1 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp1 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp1 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp1) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag7, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempz == Morton_Assist::mortonz_max )
					{
						if (boundary_z_temp[1].find(morton_temp1) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp1, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp1) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1, std::array<int, 2> { flag9, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp1_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp1_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp1) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp1, flag_refine });
						}
					}
					else if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp1) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp1, flag_refine });
						}
					}
#endif
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// check if cell at the (-z) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (z_negative)
			{
#endif
				morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel_1);
				if (exist_nodes.find(morton_temp0) == exist_nodes.end())
				{
					morton_temp0_ilevel = Morton_Assist::find_z0(iter->first, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempz == Morton_Assist::mortonz_min) // if node (morton_temp0) is on the boundary of the computatioanl domain
					{
						if (boundary_z_temp[0].find(morton_temp0) == boundary_z_temp[0].end())
						{
							boundary_z_temp[0].insert({ morton_temp0, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp0) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0, std::array<int, 2> { flag8, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp0_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp0) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp0) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp0) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp0) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#endif
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_z1 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_z1;
					}

					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+z) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (z_positive)
			{
#endif
				morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel_1);
				if (exist_nodes.find(morton_temp0) == exist_nodes.end())
				{
					morton_temp0_ilevel = Morton_Assist::find_z1(iter->first, ilevel);
#if (C_CHECK_MORTON_BOUNDARY == 1)
					morton_tempx = morton_temp0 & Morton_Assist::xmorton_ref.at(ilevel).at(1);
					morton_tempy = morton_temp0 & Morton_Assist::ymorton_ref.at(ilevel).at(1);
					morton_tempz = morton_temp0 & Morton_Assist::zmorton_ref.at(ilevel).at(1);
					if (morton_tempz == Morton_Assist::mortonz_max ) // if node (morton_temp0) is on the boundary of the computatioanl domain
					{
						if (boundary_z_temp[1].find(morton_temp0) == boundary_z_temp[1].end())
						{
							boundary_z_temp[1].insert({ morton_temp0, flag_refine });
						}
						if (refine_nodes_temp.find(morton_temp0) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0, std::array<int, 2> { flag9, 0 }));
						}
						if (refine_nodes_temp.find(morton_temp0_ilevel) == refine_nodes_temp.end())
						{
							refine_nodes_temp.insert(make_pair(morton_temp0_ilevel, std::array<int, 2> { flag_refine, 0 }));
						}
					}
					if (morton_tempx == Morton_Assist::mortonx_min)
					{
						if (boundary_x_temp[0].find(morton_temp0) == boundary_x_temp[0].end())
						{
							boundary_x_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempx == Morton_Assist::mortonx_max )
					{
						if (boundary_x_temp[1].find(morton_temp0) == boundary_x_temp[1].end())
						{
							boundary_x_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#if(C_DIMS == 3)
					if (morton_tempy == Morton_Assist::mortony_min)
					{
						if (boundary_y_temp[0].find(morton_temp0) == boundary_y_temp[0].end())
						{
							boundary_y_temp[0].insert({ morton_temp0, flag_refine });
						}
					}
					else if (morton_tempy == Morton_Assist::mortony_max )
					{
						if (boundary_y_temp[1].find(morton_temp0) == boundary_y_temp[1].end())
						{
							boundary_y_temp[1].insert({ morton_temp0, flag_refine });
						}
					}
#endif
#endif
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_z0 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_z0;
					}

					if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0_ilevel) == grid_ptr->fine2coarse.at(iboundary1).end())
					{
						grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0_ilevel, flag0 });
					}

					if (ilevel == C_max_level)
					{
						if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
						{
							grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_iboundary[0] });
						}
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
		    }
#endif
#endif
		}
	}

}


/**
* @brief      function to find the boudnary does not overlap with the coarser boundaries.
* @param[in]  ilevel    refinement level.
* @param[in]  iboundary1                order of numerical boundary at ilevel (inner).
* @param[in]  iboundary2                order of numerical boundary at ilevel (outer).
* @param[in]  iboundary2_coarse         order of numerical boundary at ilevel - 1.
* @param[out] refine_nodes        mesh where nodes on the numerical boudnary are added.
* @param[out] map_outer                nodes on the the outer layer.
* @param[out] boundary_x_temp          tempory array to store domain boundaries when x = xmin and x = xmax (only works when C_CHECK_MORTON_BOUNDARY = 1).
* @param[out] boundary_y_temp          tempory array to store domain boundaries when y = ymin and y = ymax (only works when C_CHECK_MORTON_BOUNDARY = 1).
* @param[out] boundary_z_temp          tempory array to store domain boundaries when z = xmin and z = zmax (only works when C_CHECK_MORTON_BOUNDARY = 1).
*/
template<typename  T> void Grid_Manager::search_inter_boundary(unsigned int ilevel, unsigned int iboundary1, unsigned int iboundary2, unsigned int iboundary2_coarse, D_mapint map_outer, T refine_nodes, std::array<D_mapint, 2> &boundary_x_temp, std::array<D_mapint, 2> &boundary_y_temp, std::array<D_mapint, 2> &boundary_z_temp)
{
	unsigned int ilevel_1 = ilevel - 1;
	Grid * grid_ptr;
	if (ilevel == C_max_level)
	{
		grid_ptr = &(Grid_Manager::pointer_me->gr_inner);
	}
	else
	{
		grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel));
	}
	Grid *grid_ptr_coarse = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel - 1));

	// nodes on iboundary2_coarse (ilevel - 1)
	grid_ptr_coarse->coarse2fine.at(iboundary2_coarse)= map_outer;
	// nodes on iboundary2 (ilevel)
	grid_ptr->fine2coarse.at(iboundary2) = map_outer;
	// search for nodes on iboundary2 (ilevel) which are not coincident with nodes at coarser refinement level (ilevel - 1)
/**
*@todo {searching iboundary2 could be faster if using map_outer, i.e. to check if there's a node in map_outer in the searching direction.}
*/
	D_morton morton_temp, morton_temp_ilevel_1;
	D_morton morton_temp_bool1, morton_temp_bool2;
#if(C_DIMS==3)
	D_morton morton_temp1;
#endif
	// find iboundary2 on the edges
	for (D_mapint::iterator iter = grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).begin(); iter != grid_ptr_coarse->coarse2fine.at(iboundary2_coarse).end(); ++iter)
	{
		// (-x)
		morton_temp_ilevel_1 = Morton_Assist::find_x0(iter->first, ilevel_1);
		if (map_outer.find(morton_temp_ilevel_1) != map_outer.end())
		{
			morton_temp = Morton_Assist::find_x0(iter->first, ilevel);
			if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary2).end())
			{
				grid_ptr->fine2coarse.at(iboundary2).insert({ morton_temp, flag0 });
//#if (C_CHECK_MORTON_BOUNDARY == 1)
//				if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
//				{
//					boundary_y_temp[0].insert({ morton_temp, flag_refine });
//				}
//				else if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
//				{
//					boundary_y_temp[1].insert({ morton_temp, flag_refine });
//				}
//#if (C_DIMS==3)
//				if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
//				{
//					boundary_z_temp[0].insert({ morton_temp, flag_refine });
//				}
//				else if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
//				{
//					boundary_z_temp[1].insert({ morton_temp, flag_refine });
//				}
//#endif
//#endif
			}
		}
		// (-y)
		morton_temp_ilevel_1 = Morton_Assist::find_y0(iter->first, ilevel_1);
		if (map_outer.find(morton_temp_ilevel_1) != map_outer.end())
		{
			morton_temp = Morton_Assist::find_y0(iter->first, ilevel);
			if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary2).end())
			{
				grid_ptr->fine2coarse.at(iboundary2).insert({ morton_temp, flag0 });
//#if (C_CHECK_MORTON_BOUNDARY == 1)
//				if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
//				{
//					boundary_x_temp[0].insert({ morton_temp, flag_refine });
//				}
//				else if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
//				{
//					boundary_x_temp[1].insert({ morton_temp, flag_refine });
//				}
//#if (C_DIMS==3)
//				if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
//				{
//					boundary_z_temp[0].insert({ morton_temp, flag_refine });
//				}
//				else if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
//				{
//					boundary_z_temp[1].insert({ morton_temp, flag_refine });
//				}
//#endif
//#endif
			}
		}
#if (C_DIMS == 3)
		// (-z)
		morton_temp_ilevel_1 = Morton_Assist::find_z0(iter->first, ilevel_1);
		if (map_outer.find(morton_temp_ilevel_1) != map_outer.end())
		{
			morton_temp = Morton_Assist::find_z0(iter->first, ilevel);
			if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary2).end())
			{
				grid_ptr->fine2coarse.at(iboundary2).insert({ morton_temp, flag0 });
//#if (C_CHECK_MORTON_BOUNDARY == 1)
//				if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
//				{
//					boundary_x_temp[0].insert({ morton_temp, flag_refine });
//				}
//				else if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
//				{
//					boundary_x_temp[1].insert({ morton_temp, flag_refine });
//				}
//				if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
//				{
//					boundary_y_temp[0].insert({ morton_temp, flag_refine });
//				}
//				else if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
//				{
//					boundary_y_temp[1].insert({ morton_temp, flag_refine });
//				}
//#endif
			}
		}
#endif
	}
#if (C_DIMS == 3)
	// find iboundary2 cell centers
	D_morton morton_xyz = Morton_Assist::morton_xyz.at(ilevel);
	for (D_mapint::iterator iter = grid_ptr->fine2coarse.at(iboundary1).begin(); iter != grid_ptr->fine2coarse.at(iboundary1).end(); ++iter)
	{
		if ((iter->first & morton_xyz) == morton_xyz)
		{
			// (-x)
			morton_temp = Morton_Assist::find_x0(iter->first, ilevel);
			morton_temp_ilevel_1 = Morton_Assist::find_x0(morton_temp, ilevel);
			if (refine_nodes.find(morton_temp_ilevel_1) == refine_nodes.end())
			{
				if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary1).end())
				{
					if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary2).end())
					{
						grid_ptr->fine2coarse.at(iboundary2).insert({ morton_temp, flag0 });
//#if (C_CHECK_MORTON_BOUNDARY == 1)
//						if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
//						{
//							boundary_y_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
//						{
//							boundary_y_temp[1].insert({ morton_temp, flag_refine });
//						}
//#if (C_DIMS==3)
//						if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
//						{
//							boundary_z_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
//						{
//							boundary_z_temp[1].insert({ morton_temp, flag_refine });
//						}
//#endif
//#endif
					}
				}
			}
			// (+x)
			morton_temp = Morton_Assist::find_x1(iter->first, ilevel);
			morton_temp_ilevel_1 = Morton_Assist::find_x1(morton_temp, ilevel);
			if (refine_nodes.find(morton_temp_ilevel_1) == refine_nodes.end())
			{
				if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary1).end())
				{
					if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary2).end())
					{
						grid_ptr->fine2coarse.at(iboundary2).insert({ morton_temp, flag0 });
//#if (C_CHECK_MORTON_BOUNDARY == 1)
//						if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
//						{
//							boundary_y_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
//						{
//							boundary_y_temp[1].insert({ morton_temp, flag_refine });
//						}
//#if (C_DIMS==3)
//						if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
//						{
//							boundary_z_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
//						{
//							boundary_z_temp[1].insert({ morton_temp, flag_refine });
//						}
//#endif
//#endif
					}
				}
			}
			// (-y)
			morton_temp = Morton_Assist::find_y0(iter->first, ilevel);
			morton_temp_ilevel_1 = Morton_Assist::find_y0(morton_temp, ilevel);
			if (refine_nodes.find(morton_temp_ilevel_1) == refine_nodes.end())
			{
				if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary1).end())
				{
					if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary2).end())
					{
						grid_ptr->fine2coarse.at(iboundary2).insert({ morton_temp, flag0 });
//#if (C_CHECK_MORTON_BOUNDARY == 1)
//						if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
//						{
//							boundary_x_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
//						{
//							boundary_x_temp[1].insert({ morton_temp, flag_refine });
//						}
//#if (C_DIMS==3)
//						if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
//						{
//							boundary_z_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
//						{
//							boundary_z_temp[1].insert({ morton_temp, flag_refine });
//						}
//#endif
//#endif
					}
				}
			}
			// (+y)
			morton_temp = Morton_Assist::find_y1(iter->first, ilevel);
			morton_temp_ilevel_1 = Morton_Assist::find_y1(morton_temp, ilevel);
			if (refine_nodes.find(morton_temp_ilevel_1) == refine_nodes.end())
			{
				if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary1).end())
				{
					if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary2).end())
					{
						grid_ptr->fine2coarse.at(iboundary2).insert({ morton_temp, flag0 });
//#if (C_CHECK_MORTON_BOUNDARY == 1)
//						if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
//						{
//							boundary_x_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
//						{
//							boundary_x_temp[1].insert({ morton_temp, flag_refine });
//						}
//#if (C_DIMS==3)
//						if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min)
//						{
//							boundary_z_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max)
//						{
//							boundary_z_temp[1].insert({ morton_temp, flag_refine });
//						}
//#endif
//#endif
					}
				}
			}
			// (-z)
			morton_temp = Morton_Assist::find_z0(iter->first, ilevel);
			morton_temp_ilevel_1 = Morton_Assist::find_z0(morton_temp, ilevel);
			if (refine_nodes.find(morton_temp_ilevel_1) == refine_nodes.end())
			{
				if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary1).end())
				{
					if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary2).end())
					{
						grid_ptr->fine2coarse.at(iboundary2).insert({ morton_temp, flag0 });
//#if (C_CHECK_MORTON_BOUNDARY == 1)
//						if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
//						{
//							boundary_x_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
//						{
//							boundary_x_temp[1].insert({ morton_temp, flag_refine });
//						}
//						if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
//						{
//							boundary_y_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
//						{
//							boundary_y_temp[1].insert({ morton_temp, flag_refine });
//						}
//#endif
					}
				}
			}
			// (+z)
			morton_temp = Morton_Assist::find_z1(iter->first, ilevel);
			morton_temp_ilevel_1 = Morton_Assist::find_z1(morton_temp, ilevel);
			if (refine_nodes.find(morton_temp_ilevel_1) == refine_nodes.end())
			{
				if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary1).end())
				{
					if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp) == grid_ptr->fine2coarse.at(iboundary2).end())
					{
						grid_ptr->fine2coarse.at(iboundary2).insert({ morton_temp, flag0 });
//#if (C_CHECK_MORTON_BOUNDARY == 1)
//						if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min)
//						{
//							boundary_x_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max)
//						{
//							boundary_x_temp[1].insert({ morton_temp, flag_refine });
//						}
//						if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min)
//						{
//							boundary_y_temp[0].insert({ morton_temp, flag_refine });
//						}
//						else if ((morton_temp & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max)
//						{
//							boundary_y_temp[1].insert({ morton_temp, flag_refine });
//						}
//#endif
					}
				}
			}
		}
	}
#endif

	// search for nodes on domain boundaries (ilevel) which are not coincident with nodes at coarser refinement level (ilevel - 1)
#if (C_CHECK_MORTON_BOUNDARY == 1)
	D_morton morton_temp_ilevel;
#if(C_DIMS==3)
	D_morton morton_temp1_ilevel;
#endif
	// x boundary
	std::array<D_mapint, 2> bk_boundary_x_map;
	bk_boundary_x_map = boundary_x_temp;

	// boundary (x = xmin)
	for (D_mapint::iterator iter = boundary_x_temp[0].begin(); iter != boundary_x_temp[0].end(); ++iter)
	{
		// delete node on numercial boundaries
		//if (grid_ptr->fine2coarse.at(iboundary2).find(iter->first) != grid_ptr->fine2coarse.at(iboundary2).end())
		//{
		//	grid_ptr->fine2coarse.at(iboundary2).erase(iter->first);
		//}
		if (iter->second != flag0)
		{
			morton_temp = Morton_Assist::find_y0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_y0(iter->first, ilevel);
			if (boundary_x_temp[0].find(morton_temp) != boundary_x_temp[0].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_x_map[0].find(morton_temp_ilevel) == bk_boundary_x_map[0].end())
				{
					bk_boundary_x_map[0].insert({ morton_temp_ilevel, flag_refine });
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
#if(C_DIMS==3)
			morton_temp = Morton_Assist::find_z0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_z0(iter->first, ilevel);
			if (boundary_x_temp[0].find(morton_temp) != boundary_x_temp[0].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_x_map[0].find(morton_temp_ilevel) == bk_boundary_x_map[0].end())
				{
					bk_boundary_x_map[0].insert({ morton_temp_ilevel, flag_refine });
				}
				if (boundary_x_temp[0].find(Morton_Assist::find_y0(morton_temp, ilevel_1)) != boundary_x_temp[0].end())
				{
					morton_temp1_ilevel = Morton_Assist::find_y0(morton_temp_ilevel, ilevel);
					if (bk_boundary_x_map[0].find(morton_temp1_ilevel) == bk_boundary_x_map[0].end())
					{
						bk_boundary_x_map[0].insert({ morton_temp1_ilevel, flag_refine });
					}
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
#endif
		}
	}
	for (D_mapint::iterator iter = bk_boundary_x_map[0].begin(); iter != bk_boundary_x_map[0].end(); ++iter)
	{
		bk_boundary_x.at(ilevel)[0].insert({ iter->first, nullptr });
	}
	// boundary (x = xmax)
	for (D_mapint::iterator iter = boundary_x_temp[1].begin(); iter != boundary_x_temp[1].end(); ++iter)
	{
		// delete node on numercial boundaries
		//if (grid_ptr->fine2coarse.at(iboundary2).find(iter->first) != grid_ptr->fine2coarse.at(iboundary2).end())
		//{
		//	grid_ptr->fine2coarse.at(iboundary2).erase(iter->first);
		//}
		if (iter->second != flag0)
		{
			morton_temp = Morton_Assist::find_y0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_y0(iter->first, ilevel);
			if (boundary_x_temp[1].find(morton_temp) != boundary_x_temp[1].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_x_map[1].find(morton_temp_ilevel) == bk_boundary_x_map[1].end())
				{
					bk_boundary_x_map[1].insert({ morton_temp_ilevel, flag_refine });
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
#if(C_DIMS==3)
			morton_temp = Morton_Assist::find_z0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_z0(iter->first, ilevel);
			if (boundary_x_temp[1].find(morton_temp) != boundary_x_temp[1].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_x_map[1].find(morton_temp_ilevel) == bk_boundary_x_map[1].end())
				{
					bk_boundary_x_map[1].insert({ morton_temp_ilevel, flag_refine });
				}
				if (boundary_x_temp[1].find(Morton_Assist::find_y0(morton_temp, ilevel_1)) != boundary_x_temp[1].end())
				{
					morton_temp1_ilevel = Morton_Assist::find_y0(morton_temp_ilevel, ilevel);
					if (bk_boundary_x_map[1].find(morton_temp1_ilevel) == bk_boundary_x_map[1].end())
					{
						bk_boundary_x_map[1].insert({ morton_temp1_ilevel, flag_refine });
					}
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
#endif
		}
	}
	for (D_mapint::iterator iter = bk_boundary_x_map[1].begin(); iter != bk_boundary_x_map[1].end(); ++iter)
	{
		bk_boundary_x.at(ilevel)[1].insert({ iter->first, nullptr });
	}

	// y boundary
	std::array<D_mapint, 2> bk_boundary_y_map;
	bk_boundary_y_map = boundary_y_temp;
	// boundary (y = ymin)
	for (D_mapint::iterator iter = boundary_y_temp[0].begin(); iter != boundary_y_temp[0].end(); ++iter)
	{
		// delete node on numercial boundaries
		//if (grid_ptr->fine2coarse.at(iboundary2).find(iter->first) != grid_ptr->fine2coarse.at(iboundary2).end())
		//{
		//	grid_ptr->fine2coarse.at(iboundary2).erase(iter->first);
		//}
		if (iter->second != flag0)
		{
			morton_temp = Morton_Assist::find_x0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_x0(iter->first, ilevel);
			if (boundary_y_temp[0].find(morton_temp) != boundary_y_temp[0].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_y_map[0].find(morton_temp_ilevel) == bk_boundary_y_map[0].end())
				{
					bk_boundary_y_map[0].insert({ morton_temp_ilevel, flag_refine });
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
#if(C_DIMS==3)
			morton_temp = Morton_Assist::find_z0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_z0(iter->first, ilevel);
			if (boundary_y_temp[0].find(morton_temp) != boundary_y_temp[0].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_y_map[0].find(morton_temp_ilevel) == bk_boundary_y_map[0].end())
				{
					bk_boundary_y_map[0].insert({ morton_temp_ilevel, flag_refine });
				}
				if (boundary_y_temp[0].find(Morton_Assist::find_x0(morton_temp, ilevel_1)) != boundary_y_temp[0].end())
				{
					morton_temp1_ilevel = Morton_Assist::find_x0(morton_temp_ilevel, ilevel);
					if (bk_boundary_y_map[0].find(morton_temp1_ilevel) == bk_boundary_y_map[0].end())
					{
						bk_boundary_y_map[0].insert({ morton_temp1_ilevel, flag_refine });
					}
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
#endif
		}
	}
	for (D_mapint::iterator iter = bk_boundary_y_map[0].begin(); iter != bk_boundary_y_map[0].end(); ++iter)
	{
		bk_boundary_y.at(ilevel)[0].insert({ iter->first, nullptr });
	}
	// boundary (y = ymax)
	for (D_mapint::iterator iter = boundary_y_temp[1].begin(); iter != boundary_y_temp[1].end(); ++iter)
	{
		// delete node on numercial boundaries
		//if (grid_ptr->fine2coarse.at(iboundary2).find(iter->first) != grid_ptr->fine2coarse.at(iboundary2).end())
		//{
		//	grid_ptr->fine2coarse.at(iboundary2).erase(iter->first);
		//}
		if (iter->second != flag0)
		{
			morton_temp = Morton_Assist::find_x0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_x0(iter->first, ilevel);
			if (boundary_y_temp[1].find(morton_temp) != boundary_y_temp[1].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_y_map[1].find(morton_temp_ilevel) == bk_boundary_y_map[1].end())
				{
					bk_boundary_y_map[1].insert({ morton_temp_ilevel, flag_refine });
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
#if(C_DIMS==3)
			morton_temp = Morton_Assist::find_z0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_z0(iter->first, ilevel);
			if (boundary_y_temp[1].find(morton_temp) != boundary_y_temp[1].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_y_map[1].find(morton_temp_ilevel) == bk_boundary_y_map[1].end())
				{
					bk_boundary_y_map[1].insert({ morton_temp_ilevel, flag_refine });
				}
				if (boundary_y_temp[1].find(Morton_Assist::find_x0(morton_temp, ilevel_1)) != boundary_y_temp[1].end())
				{
					morton_temp1_ilevel = Morton_Assist::find_x0(morton_temp_ilevel, ilevel);
					if (bk_boundary_y_map[1].find(morton_temp1_ilevel) == bk_boundary_y_map[1].end())
					{
						bk_boundary_y_map[1].insert({ morton_temp1_ilevel, flag_refine });
					}
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
#endif
		}
	}
	for (D_mapint::iterator iter = bk_boundary_y_map[1].begin(); iter != bk_boundary_y_map[1].end(); ++iter)
	{
		bk_boundary_y.at(ilevel)[1].insert({ iter->first, nullptr });
	}

#if(C_DIMS == 3)
	// z boundary
	std::array<D_mapint, 2> bk_boundary_z_map;
	bk_boundary_z_map = boundary_z_temp;
	// boundary (z = zmin)
	for (D_mapint::iterator iter = boundary_z_temp[0].begin(); iter != boundary_z_temp[0].end(); ++iter)
	{
		// delete node on numercial boundaries
		//if (grid_ptr->fine2coarse.at(iboundary2).find(iter->first) != grid_ptr->fine2coarse.at(iboundary2).end())
		//{
		//	grid_ptr->fine2coarse.at(iboundary2).erase(iter->first);
		//}
		//if (iter->second != flag0)
		{
			morton_temp = Morton_Assist::find_x0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_x0(iter->first, ilevel);
			if (boundary_z_temp[0].find(morton_temp) != boundary_z_temp[0].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_z_map[0].find(morton_temp_ilevel) == bk_boundary_z_map[0].end())
				{
					bk_boundary_z_map[0].insert({ morton_temp_ilevel, flag_refine });
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
			morton_temp = Morton_Assist::find_y0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_y0(iter->first, ilevel);
			if (boundary_z_temp[0].find(morton_temp) != boundary_z_temp[0].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_z_map[0].find(morton_temp_ilevel) == bk_boundary_z_map[0].end())
				{
					bk_boundary_z_map[0].insert({ morton_temp_ilevel, flag_refine });
				}
				if (boundary_z_temp[0].find(Morton_Assist::find_x0(morton_temp, ilevel_1)) != boundary_z_temp[0].end())
				{
					morton_temp1_ilevel = Morton_Assist::find_x0(morton_temp_ilevel, ilevel);
					if (bk_boundary_z_map[0].find(morton_temp1_ilevel) == bk_boundary_z_map[0].end())
					{
						bk_boundary_z_map[0].insert({ morton_temp1_ilevel, flag_refine });
					}
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
		}
	}
	for (D_mapint::iterator iter = bk_boundary_z_map[0].begin(); iter != bk_boundary_z_map[0].end(); ++iter)
	{
		bk_boundary_z.at(ilevel)[0].insert({ iter->first, nullptr });
	}
	// boundary (z = zmax)
	for (D_mapint::iterator iter = boundary_z_temp[1].begin(); iter != boundary_z_temp[1].end(); ++iter)
	{
		// delete node on numercial boundaries
		//if (grid_ptr->fine2coarse.at(iboundary2).find(iter->first) != grid_ptr->fine2coarse.at(iboundary2).end())
		//{
		//	grid_ptr->fine2coarse.at(iboundary2).erase(iter->first);
		//}
		//if (iter->second != flag0)
		{
			morton_temp = Morton_Assist::find_x0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_x0(iter->first, ilevel);
			if (boundary_z_temp[1].find(morton_temp) != boundary_z_temp[1].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_z_map[1].find(morton_temp_ilevel) == bk_boundary_z_map[1].end())
				{
					bk_boundary_z_map[1].insert({ morton_temp_ilevel, flag_refine });
				}
				// delete node on numercial boundaries
				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
			morton_temp = Morton_Assist::find_y0(iter->first, ilevel_1);
			morton_temp_ilevel = Morton_Assist::find_y0(iter->first, ilevel);
			if (boundary_z_temp[1].find(morton_temp) != boundary_z_temp[1].end())
			{
				// insert node on the domain boundary to domain boundary
				if (bk_boundary_z_map[1].find(morton_temp_ilevel) == bk_boundary_z_map[1].end())
				{
					bk_boundary_z_map[1].insert({ morton_temp_ilevel, flag_refine });
				}
				if (boundary_z_temp[1].find(Morton_Assist::find_x0(morton_temp, ilevel_1)) != boundary_z_temp[1].end())
				{
					morton_temp1_ilevel = Morton_Assist::find_x0(morton_temp_ilevel, ilevel);
					if (bk_boundary_z_map[1].find(morton_temp1_ilevel) == bk_boundary_z_map[1].end())
					{
						bk_boundary_z_map[1].insert({ morton_temp1_ilevel, flag_refine });
					}
				}
				// delete node on numercial boundaries

				//if (grid_ptr->fine2coarse.at(iboundary2).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary2).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary2).erase(morton_temp_ilevel);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp);
				//}
				//if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp_ilevel) != grid_ptr->fine2coarse.at(iboundary1).end())
				//{
				//	grid_ptr->fine2coarse.at(iboundary1).erase(morton_temp_ilevel);
				//}
			}
		}
	}
	for (D_mapint::iterator iter = bk_boundary_z_map[1].begin(); iter != bk_boundary_z_map[1].end(); ++iter)
	{
		bk_boundary_z.at(ilevel)[1].insert({ iter->first, nullptr });
	}
#endif
#endif

}


/**
* @brief      function to generate mesh at intermediate refinement level (from C_max_level - 1 to 1).
* @param[in]  ilevel    refinement level.
* @note       Mesh generation is based on the numerical boudanry (coarse2fine) obained in the previous step (ilevel + 1). The mesh is expanded one layer by one layer.
*/
void Grid_Manager::generate_intermediate(unsigned int ilevel)
{
	Grid * grid_ptr;	
	grid_ptr = &gr_NoIB.at(ilevel);

	grid_ptr->dx = C_dx / static_cast<D_real> (two_power_n(ilevel));

	Grid *grid_ptr_coarse = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel - 1));
	D_morton morton_xyz = Morton_Assist::morton_xyz.at(ilevel);

	D_mapint2 refine_nodes;
	for (D_mapint::iterator iter = grid_ptr->coarse2fine.at(iboundary0_coarse).begin(); iter != grid_ptr->coarse2fine.at(iboundary0_coarse).end(); ++iter)
	{
		if (refine_nodes.find(iter->first) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(iter->first, std::array<int, 2> { flag_refine, 0 }));
		}
	}
	for (D_mapint::iterator iter = grid_ptr->coarse2fine.at(iboundary2_coarse).begin(); iter != grid_ptr->coarse2fine.at(iboundary2_coarse).end(); ++iter)
	{
		if (refine_nodes.find(iter->first) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(iter->first, std::array<int, 2> { flag_refine, 0 }));
		}
	}

    D_mapint map_outer = grid_ptr->coarse2fine.at(iboundary2_coarse);
	std::array<D_mapint, 2> boundary_x_temp, boundary_y_temp, boundary_z_temp; // temporary array to store domain boundaries at ilevel, only work if C_CHECK_MORTON_BOUNDARY = 1
#if (C_CHECK_MORTON_BOUNDARY == 1)
	unsigned int ilevel_p1 = ilevel + 1;
	D_morton morton_xyz_p1 = Morton_Assist::morton_xyz.at(ilevel_p1);
	// add overlapping nodes to the domain boundary
	for (unsigned int iboundary = 0; iboundary < 2; ++iboundary)
	{
		for (D_mapNodePtr::iterator iter = bk_boundary_x.at(ilevel_p1).at(iboundary).begin(); iter != bk_boundary_x.at(ilevel_p1).at(iboundary).end(); ++iter)
		{
			if (((iter->first&morton_xyz_p1) == 0) && (grid_ptr->coarse2fine.at(0).find(iter->first) != grid_ptr->coarse2fine.at(0).end()))
			{
				boundary_x_temp[iboundary].insert({ iter->first, flag_refine });
			}
		}
		for (D_mapNodePtr::iterator iter = bk_boundary_y.at(ilevel_p1).at(iboundary).begin(); iter != bk_boundary_y.at(ilevel_p1).at(iboundary).end(); ++iter)
		{
			if (((iter->first &morton_xyz_p1) == 0) && (grid_ptr->coarse2fine.at(0).find(iter->first) != grid_ptr->coarse2fine.at(0).end()))
			{
				boundary_y_temp[iboundary].insert({ iter->first, flag_refine });
			}
		}
#if(C_DIMS == 3)
		for (D_mapNodePtr::iterator iter = bk_boundary_z.at(ilevel_p1).at(iboundary).begin(); iter != bk_boundary_z.at(ilevel_p1).at(iboundary).end(); ++iter)
		{
			if (((iter->first &morton_xyz_p1) == 0) && (grid_ptr->coarse2fine.at(0).find(iter->first) != grid_ptr->coarse2fine.at(0).end()))
			{
				boundary_z_temp[iboundary].insert({ iter->first, flag_refine });
			}
		}
#endif
	}
#endif


	D_mapint map_boundary_temp = {};
	D_morton morton_temp_center;
	D_morton morton_temp0, morton_temp1; // 0: faces; 1: edges; 2: corners
#if(C_DIMS==3)
	D_morton morton_temp2;
#endif

	unsigned int extend_temp = C_extend;
#if (C_SOLID_BOUNDARY == 2) // since searching added and removed nodes are based on (ilevel - 1) refinement
	                        //, to avoid the inital grid exceeds the searching region, reset number of extended nodes
	extend_temp -= 2;
#endif

	// identify in which direction mesh is expanded by the most number of layers
	unsigned int C_extend_max = C_extend_outer_x0;
	if (C_extend_outer_x1 > C_extend_max) C_extend_max = C_extend_outer_x1;
	if (C_extend_outer_y0 > C_extend_max) C_extend_max = C_extend_outer_y0;
	if (C_extend_outer_y1 > C_extend_max) C_extend_max = C_extend_outer_y1;
#if (C_DIMS == 3)
	if (C_extend_outer_z0 > C_extend_max) C_extend_max = C_extend_outer_z0;
	if (C_extend_outer_z1 > C_extend_max) C_extend_max = C_extend_outer_z1;
#endif

	// expand mesh layer based on previous layer
	bool bool_extend_x0 = true, bool_extend_x1 = true, bool_extend_y0 = true, bool_extend_y1 = true;
#if (C_DIMS == 3)
	bool bool_extend_z0 = true, bool_extend_z1 = true;
#endif

	D_real bx0_temp = Solid_Manager::pointer_me->shape_offest_x0_grid + C_eps;
	D_real by0_temp = Solid_Manager::pointer_me->shape_offest_y0_grid + C_eps;
#if(C_DIMS == 3)
	D_real bz0_temp = Solid_Manager::pointer_me->shape_offest_z0_grid + C_eps;
#endif
	for (unsigned int iextend = 1; iextend < (extend_temp + C_extend_max); ++iextend)
	{
		if (iextend > extend_temp)
		{
			unsigned int iextend_outer = iextend - extend_temp;
			if (C_extend_outer_x0 < iextend_outer) bool_extend_x0 = false;
			if (C_extend_outer_x1 < iextend_outer) bool_extend_x1 = false;
			if (C_extend_outer_y0 < iextend_outer) bool_extend_y0 = false;
			if (C_extend_outer_y1 < iextend_outer) bool_extend_y1 = false;
#if (C_DIMS == 3)
			if (C_extend_outer_z0 < iextend_outer) bool_extend_z0 = false;
			if (C_extend_outer_z1 < iextend_outer) bool_extend_z1 = false;
#endif
		}
		D_mapint map_current = map_outer;
		map_outer.clear();
		for (D_mapint::iterator iter = map_current.begin(); iter != map_current.end(); ++iter)
		{
#if (C_CHECK_MORTON_BOUNDARY == 1)
#if (C_DIMS == 2)
			D_real x, y;
			Morton_Assist::pointer_me->compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
			D_real x, y, z;
			Morton_Assist::pointer_me->compute_coordinate(iter->first, ilevel, x, y, z);
#endif
			// check if the current node is not at boundary of the computational domain
			bool x_positive = true, x_negative = true, y_positive = true, y_negative = true;
			if (x < bx0_temp)
			{
				x_negative = false;
				boundary_x_temp[0].insert({ iter->first, flag_refine });
			}
			else if (x + C_eps > xb_domain)
			{
				x_positive = false;
				boundary_x_temp[1].insert({ iter->first, flag_refine });
			}


			if (y < by0_temp)
			{
				y_negative = false;
				boundary_y_temp[0].insert({ iter->first, flag_refine });
			}
			else if (y + C_eps > yb_domain)
			{
				y_positive = false;
				boundary_y_temp[1].insert({ iter->first, flag_refine });
			}

#if (C_DIMS == 3)
			bool z_positive = true, z_negative = true;
			if (z < bz0_temp)
			{
				z_negative = false;
				boundary_z_temp[0].insert({ iter->first, flag_refine });
			}
			else if (z + C_eps > zb_domain)
			{
				z_positive = false;
				boundary_z_temp[1].insert({ iter->first, flag_refine });
			}
#endif
#endif
			// edge or corners, or nodes relavent to merged boundary in the previous step (merge_nodes)

			   // check if cell at the (-x) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative)
			{
#endif
				morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
				if (bool_extend_x0 && refine_nodes.find(morton_temp0) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp0) == map_outer.end())
					{

						map_outer.insert({ morton_temp0, bit_x1 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_x1;
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-x, -y) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
				if (bool_extend_x0 && bool_extend_y0 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (-x, -y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_negative && z_negative)
			{
#endif
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
				if (bool_extend_x0 && bool_extend_y0 && bool_extend_z0 && refine_nodes.find(morton_temp2) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-x, -y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_negative && z_positive)
			{
#endif
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
				if (bool_extend_x0 && bool_extend_y0 && bool_extend_z1 && refine_nodes.find(morton_temp2) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-x, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && z_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
				if (bool_extend_x0 && bool_extend_z0 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}

				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-x, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && z_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
				if (bool_extend_x0 && bool_extend_z1 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// check if cell at the (-x, +y) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
				if (bool_extend_x0 && bool_extend_y1 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (-x, +y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_positive && z_negative)
			{
#endif
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
				if (bool_extend_x0 && bool_extend_y1 && bool_extend_z0 && refine_nodes.find(morton_temp2) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-x, +y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_positive && z_positive)
			{
#endif
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
				if (bool_extend_x0 && bool_extend_y1 && bool_extend_z1 && refine_nodes.find(morton_temp2) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif

			// check if cell at the (+x) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive)
			{
#endif
				morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
				if (bool_extend_x1 && refine_nodes.find(morton_temp0) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_x0 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_x0;
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+x, -y) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
				if (bool_extend_x1 && bool_extend_y0 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (+x, -y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_negative && z_negative)
			{
#endif
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
				if (bool_extend_x1 && bool_extend_y0 && bool_extend_z0 && refine_nodes.find(morton_temp2) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+x, -y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_negative && z_positive)
			{
#endif
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
				if (bool_extend_x1 && bool_extend_y0 && bool_extend_z1 && refine_nodes.find(morton_temp2) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// check if cell at the (+x, +y) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
				if (bool_extend_x1 && bool_extend_y1 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (+x, +y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_positive & z_negative)
			{
#endif
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
				if (bool_extend_x1 && bool_extend_y1 && bool_extend_z0 && refine_nodes.find(morton_temp2) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+x, +y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_positive & z_positive)
			{
#endif
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
				if (bool_extend_x1 && bool_extend_y1 && bool_extend_z1 && refine_nodes.find(morton_temp2) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp2) == map_outer.end())
					{
						map_outer.insert({ morton_temp2, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+x, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive & z_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
				if (bool_extend_x1 && bool_extend_z0 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+x, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive & z_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
				if (bool_extend_x1 && bool_extend_z1 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif

			// check if cell at the (-y) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_negative)
			{
#endif
				morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
				if (bool_extend_y0 && refine_nodes.find(morton_temp0) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_y1 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_y1;
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (-y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_negative && z_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
				if (bool_extend_y0 && bool_extend_z0 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (-y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_negative && z_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
				if (bool_extend_y0 && bool_extend_z1 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// check if cell at the (+y) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_positive)
			{
#endif
				morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
				if (bool_extend_y1 && refine_nodes.find(morton_temp0) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_y0 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_y0;
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_DIMS==3)
			// check if cell at the (+y, -z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_positive && z_negative)
			{
#endif
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
				if (bool_extend_y1 && bool_extend_z0 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+y, +z) exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_positive && z_positive)
			{
#endif
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
				if (bool_extend_y1 && bool_extend_z1 && refine_nodes.find(morton_temp1) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp1) == map_outer.end())
					{
						map_outer.insert({ morton_temp1, 0 });
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// check if cell at the (-z) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (z_negative)
			{
#endif
				morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
				if (bool_extend_z0 && refine_nodes.find(morton_temp0) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_z1 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_z1;
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// check if cell at the (+z) face exist
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (z_positive)
			{
#endif
				morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
				if (bool_extend_z1 && refine_nodes.find(morton_temp0) == refine_nodes.end())
				{
					if (map_outer.find(morton_temp0) == map_outer.end())
					{
						map_outer.insert({ morton_temp0, bit_z0 });
					}
					else
					{
						map_outer[morton_temp0] = map_outer[morton_temp0] | bit_z0;
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
		}


		for (D_mapint::iterator iter = map_outer.begin(); iter != map_outer.end(); ++iter)
		{
			if (refine_nodes.find(iter->first) == refine_nodes.end())
			{
				refine_nodes.insert(make_pair(iter->first, std::array<int, 2> { flag_refine, 0 }));
			}
		}

		bool bool_merge; // identify if the node have two adjacent nodes in one direction, e.g. -x and +x
		std::vector<D_morton> merge_nodes;
		// find nodes are not the boudnary of the layer
		for (D_mapint::iterator iter = map_outer.begin(); iter != map_outer.end(); ++iter)
		{
			bool_merge = ((iter->second & bit_x01) == bit_x01) | ((iter->second & bit_y01) == bit_y01);
#if (C_DIMS==3)
			bool_merge = bool_merge | ((iter->second & bit_z01) == bit_z01);
#endif
			if (bool_merge)
			{
				merge_nodes.push_back(iter->first);
			}

		}

		// store nodes in the last two layers to search for iboundary0_coarse
		if (iextend > (extend_temp - 3)) // the last two layers
		{
			map_boundary_temp.insert(map_outer.begin(), map_outer.end());
		}

		// delete nodes which are not the boundary of the layer
		for (std::vector<D_morton>::iterator iter = merge_nodes.begin(); iter != merge_nodes.end(); ++iter)
		{
			map_outer.erase(*iter);
		}

	}


	unsigned int ilevel_1 = ilevel - 1;
	D_morton morton_temp0_ilevel_1, morton_temp1_ilevel_1;
#if (C_DIMS==3)
	D_morton morton_temp2_ilevel_1;
#endif
	for (D_mapint::iterator iter = map_boundary_temp.begin(); iter != map_boundary_temp.end(); ++iter)
	{
		if ((iter->first&morton_xyz) == 0)
		{
#if (C_CHECK_MORTON_BOUNDARY == 1)
#if (C_DIMS == 2)
			D_real x, y;
			Morton_Assist::pointer_me->compute_coordinate(iter->first, ilevel, x, y);
#endif
#if (C_DIMS == 3)
			D_real x, y, z;
			Morton_Assist::pointer_me->compute_coordinate(iter->first, ilevel, x, y, z);
#endif
			// check if the current node is not at boundary of the computational domain
			bool x_positive = true, x_negative = true, y_positive = true, y_negative = true;
			if (x < bx0_temp)
			{
				x_negative = false;
				boundary_x_temp[0].insert({ iter->first, flag_refine });
			}
			else if (x + C_eps > xb_domain)
			{
				x_positive = false;
				boundary_x_temp[1].insert({ iter->first, flag_refine });
			}


			if (y < by0_temp)
			{
				y_negative = false;
				boundary_y_temp[0].insert({ iter->first, flag_refine });
			}
			else if (y + C_eps > yb_domain)
			{
				y_positive = false;
				boundary_y_temp[1].insert({ iter->first, flag_refine });
			}

#if (C_DIMS == 3)
			bool z_positive = true, z_negative = true;
			if (z < by0_temp)
			{
				z_negative = false;
				boundary_z_temp[0].insert({ iter->first, flag_refine });
			}
			else if (z + C_eps > zb_domain)
			{
				z_positive = false;
				boundary_z_temp[1].insert({ iter->first, flag_refine });
			}
#endif
#endif
			bool bool_N_boundary = true;
			// (-x) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative)
			{
#endif
				morton_temp0_ilevel_1 = Morton_Assist::find_x0(iter->first, ilevel_1);
				morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
				if (refine_nodes.find(morton_temp0_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp0) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (-x, -y)
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_negative)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

#if(C_DIMS==3)
			// (-x, -y, -z)
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_negative && z_negative)
			{
#endif
				morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
				if (refine_nodes.find(morton_temp2_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp2) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (-x, -y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_negative && z_positive)
			{
#endif
				morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
				if (refine_nodes.find(morton_temp2_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp2) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (-x, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && z_negative)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// (-x, +y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_positive)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

#if(C_DIMS==3)
			// (-x, +y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_positive && z_negative)
			{
#endif
				morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
				if (refine_nodes.find(morton_temp2_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp2) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (-x, +y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && y_positive && z_positive)
			{
#endif
				morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
				if (refine_nodes.find(morton_temp2_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp2) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (-x, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_negative && z_positive)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// (+x) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive)
			{
#endif
				morton_temp0_ilevel_1 = Morton_Assist::find_x1(iter->first, ilevel_1);
				morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
				if (refine_nodes.find(morton_temp0_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp0) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (+x, -y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_negative)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
				if ( refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

#if(C_DIMS==3)
			// (+x, -y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_negative && z_negative)
			{
#endif
				morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
				if (refine_nodes.find(morton_temp2_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp2) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (+x, -y, +z)
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_negative && z_positive)
			{
#endif
				morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
				if (refine_nodes.find(morton_temp2_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp2) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (+x, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && z_negative)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// (+x, +y) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_positive)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

#if(C_DIMS==3)
			// (+x, +y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_positive && z_negative)
			{
#endif
				morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
				morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
				if (refine_nodes.find(morton_temp2_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp2) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			// (+x, +y, +z)
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && y_positive && z_positive)
			{
#endif
				morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
				morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
				if (refine_nodes.find(morton_temp2_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp2) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (+x, +z)
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (x_positive && z_positive)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// (-y)
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_negative)
			{
#endif
				morton_temp0_ilevel_1 = Morton_Assist::find_y0(iter->first, ilevel_1);
				morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
				if (refine_nodes.find(morton_temp0_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp0) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

#if(C_DIMS==3)
			// (-y, -z)
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_negative && z_negative)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (-y, +z)
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_negative && z_positive)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			// (+y)
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_positive)
			{
#endif
				morton_temp0_ilevel_1 = Morton_Assist::find_y1(iter->first, ilevel_1);
				morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
				if (refine_nodes.find(morton_temp0_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp0) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

#if(C_DIMS==3)
			// (+y, -z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_positive && z_negative)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (+y, +z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (y_positive && z_positive)
			{
#endif
				morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
				morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
				if (refine_nodes.find(morton_temp1_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp1) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (-z)
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (z_negative)
			{
#endif
				morton_temp0_ilevel_1 = Morton_Assist::find_z0(iter->first, ilevel_1);
				morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
				if (refine_nodes.find(morton_temp0_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp0) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif

			// (+z) 
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (z_positive)
			{
#endif
				morton_temp0_ilevel_1 = Morton_Assist::find_z1(iter->first, ilevel_1);
				morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
				if (refine_nodes.find(morton_temp0_ilevel_1) == refine_nodes.end())
				{
					bool_N_boundary = false;
				}
				else
				{
					if (refine_nodes.find(morton_temp0) == refine_nodes.end())
					{
						refine_nodes.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
					}
				}
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#endif
			if (!bool_N_boundary)
			{
				grid_ptr_coarse->coarse2fine.at(0).insert({ iter->first, flag_iboundary[0] });
			}
		}
	}

	// find the outer layer
	D_mapint2 refine_nodes_temp = {};
	D_mapint2 boundary_temp = {};
	for (D_mapint::iterator iter = grid_ptr_coarse->coarse2fine.at(0).begin(); iter != grid_ptr_coarse->coarse2fine.at(0).end(); ++iter)
	{
		boundary_temp.insert(make_pair(iter->first, std::array<int, 2> { flag_refine, 0 }));
	}

	map_outer.clear();
	
	search_outer_boundary(ilevel, boundary_temp, refine_nodes, refine_nodes_temp, map_outer, boundary_x_temp, boundary_y_temp, boundary_z_temp);
	
	for (D_mapint2::iterator iter = refine_nodes_temp.begin(); iter != refine_nodes_temp.end(); ++iter)
	{
		if (refine_nodes.find(iter->first) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(iter->first, std::array<int, 2> { iter->second.at(0), 0 }));
		}
	}
	
	// add boundary nodes to refinenodes
	for (D_mapint::iterator iter = grid_ptr->fine2coarse.at(iboundary1).begin(); iter != grid_ptr->fine2coarse.at(iboundary1).end(); ++iter)
	{
		if (refine_nodes.find(iter->first) == refine_nodes.end())
		{
			refine_nodes.insert(make_pair(iter->first, std::array<int, 2> { flag0, 0 }));
		}
		else
		{
			refine_nodes[iter->first].at(0) = flag0;
		}
	}

	// delete nodes on overlaping layers
	D_mapint merge_refine = {};
	check_merge_nodes(ilevel, refine_nodes, map_outer);

	// find iboundary2
	search_inter_boundary(ilevel, iboundary1, iboundary2, iboundary2_coarse, map_outer, refine_nodes, boundary_x_temp, boundary_y_temp, boundary_z_temp);


	// add boundary nodes to refinenodes
	for (int iboundary = 1; iboundary < C_overlap * 2 + 1; ++iboundary)
	{
		for (D_mapint::iterator iter = grid_ptr->fine2coarse.at(iboundary).begin(); iter != grid_ptr->fine2coarse.at(iboundary).end(); ++iter)
		{
			if (refine_nodes.find(iter->first) == refine_nodes.end())
			{
				refine_nodes.insert(make_pair(iter->first, std::array<int, 2> { flag0, 0 }));
			}
		}
	}

	// generate nodes
	for (D_mapint2::iterator iter = refine_nodes.begin(); iter != refine_nodes.end(); ++iter)
	{
		Node temp;
		temp.flag = flag_refine;
		gr_NoIB.at(ilevel).grid.insert({ iter->first, temp });
	}

	// set flag at numerical boundaries as flag_iboundary
	for (int iboundary = 1; iboundary < C_overlap * 2 + 1; ++iboundary)
	{
		for (D_mapint::iterator iter = grid_ptr->fine2coarse.at(iboundary).begin(); iter != grid_ptr->fine2coarse.at(iboundary).end(); ++iter)
		{
			gr_NoIB.at(ilevel).grid[iter->first].flag = flag_iboundary[iboundary];
		}

	}

}

/**
* @brief      function to check nodes need to be merged.
* @param[in]  ilevel              refinement level.
* @param[in]  nodes_exist         exsiting nodes.
* @param[out] map_outer           nodes at the outer layer.
*/

void Grid_Manager::check_merge_nodes(unsigned int ilevel, D_mapint2 &nodes_exist, D_mapint &map_outer)
{
	unsigned int ilevel_1 = ilevel - 1;
	Grid * grid_ptr;
	if (ilevel == C_max_level)
	{
		grid_ptr = &(Grid_Manager::pointer_me->gr_inner);
	}
	else
	{
		grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel));
	}

	Grid *grid_ptr_coarse = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel - 1));

	// add map_outer to nodes exist
	for (D_mapint::iterator iter = map_outer.begin(); iter != map_outer.end(); ++iter)
	{
		if (nodes_exist.find(iter->first) == nodes_exist.end())
		{
			nodes_exist.insert(make_pair(iter->first, std::array<int, 2> { flag_refine, 0 }));
		}
	}

	D_morton morton_temp0, morton_temp1, morton_temp0_ilevel_1, morton_temp1_ilevel_1;
#if (C_DIMS==3)
	D_morton morton_temp2, morton_temp2_ilevel_1;
#endif
	// find nodes in the megred layer
	D_mapint merge_nodes;
	for (D_mapint::iterator iter = map_outer.begin(); iter != map_outer.end(); ++iter)
	{
		bool bool_merge  = true;
		// (-x) 
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel_1);
		if (nodes_exist.find(morton_temp0) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (-x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}

#if(C_DIMS==3)
		// (-x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp2) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (-x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp2) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (-x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}
#endif

		// (-x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}

#if(C_DIMS==3)
		// (-x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp2) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (-x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp2) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (-x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}

#endif
		// (+x) 
		morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp0) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (+x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}

#if(C_DIMS==3)
		// (+x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp2) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (+x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp2) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (+x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}


#endif
		// (+x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}

#if(C_DIMS==3)
		// (+x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp2) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (+x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp2) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (+x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}

#endif
		// (-y) 
		morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp0) == nodes_exist.end())
		{
			bool_merge = false;
		}

#if(C_DIMS==3)
		// (-y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (-y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}

#endif
		// (+y) 
		morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp0) == nodes_exist.end())
		{
			bool_merge = false;
		}

#if(C_DIMS==3)
		// (+y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (+y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp1) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (-z) 
		morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp0) == nodes_exist.end())
		{
			bool_merge = false;
		}

		// (+z) 
		morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel_1);
		if (bool_merge&&nodes_exist.find(morton_temp0) == nodes_exist.end())
		{
			bool_merge = false;
		}

		if (bool_merge)
		{
			if (merge_nodes.find(iter->first) == merge_nodes.end())
			{
				merge_nodes.insert({ iter->first ,flag_refine });
			}
		}

#endif
	}

	// delete merged nodes from boundary2 and set flag as flag_refine
	for (D_mapint::iterator iter = merge_nodes.begin(); iter != merge_nodes.end(); ++iter)
	{
		if (map_outer.find(iter->first) != map_outer.end())
		{
			map_outer.erase(iter->first);
		}
		if (nodes_exist.find(iter->first) != nodes_exist.end())
		{
			nodes_exist[iter->first].at(0) = flag_refine;
		}
	}

	// find nodes surround the merge_nodes
	D_mapint surround_nodes_ilevel = {};
	D_mapint surround_nodes_ilevel_1 = {};
	for (D_mapint::iterator iter = merge_nodes.begin(); iter != merge_nodes.end(); ++iter)
	{
		// (-x)
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_x0(iter->first, ilevel_1);
		if (map_outer.find(morton_temp0_ilevel_1) == map_outer.end())
		{
			if (surround_nodes_ilevel_1.find(morton_temp0_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp0_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp0) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp0, flag0 });
			}
			if (nodes_exist.find(morton_temp0) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp0].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0 , flag0 });
			}
		}

		// (-x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

#if(C_DIMS==3)
		// (-x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp2_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp2_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp2_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp2) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp2, flag0 });
			}
			if (nodes_exist.find(morton_temp2) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp2].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2 , flag0 });
			}
		}

		// (-x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp2_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp2_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp2_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp2) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp2, flag0 });
			}
			if (nodes_exist.find(morton_temp2) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp2].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2 , flag0 });
			}
		}

		// (-x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}
#endif

		// (-x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

#if(C_DIMS==3)
		// (-x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp2_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp2_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp2_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp2) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp2, flag0 });
			}
			if (nodes_exist.find(morton_temp2) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp2].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2 , flag0 });
			}
		}

		// (-x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp2_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp2_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp2_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp2) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp2, flag0 });
			}
			if (nodes_exist.find(morton_temp2) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp2].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2 , flag0 });
			}
		}

		// (-x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

#endif
		// (+x) 
		morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_x1(iter->first, ilevel_1);
		if (map_outer.find(morton_temp0_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp0_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp0_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp0) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp0, flag0 });
			}
			if (nodes_exist.find(morton_temp0) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp0].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0 , flag0 });
			}
		}

		// (+x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

#if(C_DIMS==3)
		// (+x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp2_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp2_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp2_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp2) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp2, flag0 });
			}
			if (nodes_exist.find(morton_temp2) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp2].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2 , flag0 });
			}
		}

		// (+x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp2_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp2_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp2_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp2) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp2, flag0 });
			}
			if (nodes_exist.find(morton_temp2) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp2].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2 , flag0 });
			}
		}

		// (+x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

#endif
		// (+x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

#if(C_DIMS==3)
		// (+x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp2_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp2_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp2_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp2) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp2, flag0 });
			}
			if (nodes_exist.find(morton_temp2) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp2].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2 , flag0 });
			}
		}

		// (+x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp2_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp2_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp2_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp2) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp2, flag0 });
			}
			if (nodes_exist.find(morton_temp2) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp2, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp2].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp2) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp2 , flag0 });
			}
		}

		// (+x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

#endif
		// (-y) 
		morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_y0(iter->first, ilevel_1);
		if (map_outer.find(morton_temp0_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp0_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp0_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp0) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp0, flag0 });
			}
			if (nodes_exist.find(morton_temp0) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp0].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0 , flag0 });
			}
		}

#if(C_DIMS==3)
		// (-y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

		// (-y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

#endif
		// (+y) 
		morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_y1(iter->first, ilevel_1);
		if (map_outer.find(morton_temp0_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp0_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp0_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp0) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp0, flag0 });
			}
			if (nodes_exist.find(morton_temp0) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp0].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0 , flag0 });
			}
		}

#if(C_DIMS==3)
		// (+y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

		// (+y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (map_outer.find(morton_temp1_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp1_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp1_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp1) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp1, flag0 });
			}
			if (nodes_exist.find(morton_temp1) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp1, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp1].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp1) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp1 , flag0 });
			}
		}

		// (-z) 
		morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_z0(iter->first, ilevel_1);
		if (map_outer.find(morton_temp0_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp0_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp0_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp0) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp0, flag0 });
			}
			if (nodes_exist.find(morton_temp0) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp0].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0 , flag0 });
			}
		}

		// (+z) 
		morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
		morton_temp0_ilevel_1 = Morton_Assist::find_z1(iter->first, ilevel_1);
		if (map_outer.find(morton_temp0_ilevel_1) == map_outer.end())
		{

			if (surround_nodes_ilevel_1.find(morton_temp0_ilevel_1) == surround_nodes_ilevel_1.end())
			{
				surround_nodes_ilevel_1.insert({ morton_temp0_ilevel_1, flag_refine });
			}
			if (surround_nodes_ilevel.find(morton_temp0) == surround_nodes_ilevel.end())
			{
				surround_nodes_ilevel.insert({ morton_temp0, flag0 });
			}
			if (nodes_exist.find(morton_temp0) == nodes_exist.end())
			{
				nodes_exist.insert(make_pair(morton_temp0, std::array<int, 2> { flag_refine, 0 }));
			}
			else
			{
				nodes_exist[morton_temp0].at(0) = flag_refine;
			}
		}
		else
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) == grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).insert({ iter->first, flag_refine });
			}
			if (grid_ptr->fine2coarse.at(iboundary1).find(morton_temp0) == grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).insert({ morton_temp0 , flag0 });
			}
		}
#endif
	}

	// delete nodes which are not on the boundary in coarse2fine
	
	for (D_mapint::iterator iter = surround_nodes_ilevel_1.begin(); iter != surround_nodes_ilevel_1.end(); ++iter)
	{
		bool bool_N_boundary = true;
		// (-x)
		morton_temp0_ilevel_1 = Morton_Assist::find_x0(iter->first, ilevel_1);
		if (map_outer.find(morton_temp0_ilevel_1) != map_outer.end())
		{			
			bool_N_boundary = false;
		}

		// (-x, -y) 
		morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);		
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (-x, -y, -z) 
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp2_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-x, -y, +z) 
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp2_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-x, -z) 
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}
#endif

		// (-x, +y) 
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (-x, +y, -z) 
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp2_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-x, +y, +z) 
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp2_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-x, +z) 
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#endif
		// (+x) 
		morton_temp0_ilevel_1 = Morton_Assist::find_x1(iter->first, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp0_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+x, -y) 
		morton_temp1_ilevel_1 = Morton_Assist::find_y0(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (+x, -y, -z) 
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp2_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+x, -y, +z) 
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp2_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+x, -z) 
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}


#endif
		// (+x, +y) 
		morton_temp1_ilevel_1 = Morton_Assist::find_y1(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (+x, +y, -z) 
		morton_temp2_ilevel_1 = Morton_Assist::find_z0(morton_temp1_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp2_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+x, +y, +z) 
		morton_temp2_ilevel_1 = Morton_Assist::find_z1(morton_temp1_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp2_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+x, +z) 
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#endif
		// (-y) 
		morton_temp0_ilevel_1 = Morton_Assist::find_y0(iter->first, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp0_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (-y, -z) 
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-y, +z) 
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#endif
		// (+y) 
		morton_temp0_ilevel_1 = Morton_Assist::find_y1(iter->first, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp0_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (+y, -z) 
		morton_temp1_ilevel_1 = Morton_Assist::find_z0(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+y, +z) 
		morton_temp1_ilevel_1 = Morton_Assist::find_z1(morton_temp0_ilevel_1, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp1_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-z) 
		morton_temp0_ilevel_1 = Morton_Assist::find_z0(iter->first, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp0_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+z) 
		morton_temp0_ilevel_1 = Morton_Assist::find_z1(iter->first, ilevel_1);
		if (bool_N_boundary && map_outer.find(morton_temp0_ilevel_1) != map_outer.end())
		{
			bool_N_boundary = false;
		}
#endif
		if (bool_N_boundary)
		{
			if (grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).find(iter->first) != grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).end())
			{
				grid_ptr_coarse->coarse2fine.at(iboundary0_coarse).erase(iter->first);
			}
		}
	}

	for (D_mapint::iterator iter = surround_nodes_ilevel.begin(); iter != surround_nodes_ilevel.end(); ++iter)
	{
		bool bool_N_boundary = true;
		// (-x)
		morton_temp0 = Morton_Assist::find_x0(iter->first, ilevel);
		if (map_outer.find(morton_temp0) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (-x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp2) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp2) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}
#endif

		// (-x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (-x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp2) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp2) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#endif
		// (+x) 
		morton_temp0 = Morton_Assist::find_x1(iter->first, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp0) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+x, -y) 
		morton_temp1 = Morton_Assist::find_y0(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (+x, -y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp2) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+x, -y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp2) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+x, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}


#endif
		// (+x, +y) 
		morton_temp1 = Morton_Assist::find_y1(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (+x, +y, -z) 
		morton_temp2 = Morton_Assist::find_z0(morton_temp1, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp2) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+x, +y, +z) 
		morton_temp2 = Morton_Assist::find_z1(morton_temp1, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp2) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+x, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#endif
		// (-y) 
		morton_temp0 = Morton_Assist::find_y0(iter->first, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp0) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (-y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#endif
		// (+y) 
		morton_temp0 = Morton_Assist::find_y1(iter->first, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp0) != map_outer.end())
		{
			bool_N_boundary = false;
		}

#if(C_DIMS==3)
		// (+y, -z) 
		morton_temp1 = Morton_Assist::find_z0(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+y, +z) 
		morton_temp1 = Morton_Assist::find_z1(morton_temp0, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp1) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (-z) 
		morton_temp0 = Morton_Assist::find_z0(iter->first, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp0) != map_outer.end())
		{
			bool_N_boundary = false;
		}

		// (+z) 
		morton_temp0 = Morton_Assist::find_z1(iter->first, ilevel);
		if (bool_N_boundary && map_outer.find(morton_temp0) != map_outer.end())
		{
			bool_N_boundary = false;
		}
#endif
		if (bool_N_boundary)
		{
			if (grid_ptr->fine2coarse.at(iboundary1).find(iter->first) != grid_ptr->fine2coarse.at(iboundary1).end())
			{
				grid_ptr->fine2coarse.at(iboundary1).erase(iter->first);
			}
		}
	}
	
}

/**
* @brief function to generate backgound mesh.
*/
void Grid_Manager::generate_background()
{

	int ilevel = 0;
	Grid_NIB * grid_ptr;
	grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel));
	grid_ptr->dx = C_dx;
	nx = static_cast <D_uint> (C_xb / C_dx) + 1 + C_x0b_offset;
	ny = static_cast <D_uint> (C_yb / C_dx) + 1 + C_y0b_offset;

	D_uint nix = nx - 1;
	D_uint niy = ny - 1;

#if (C_DIMS==2)
	D_uint nmax = two_power_n(Morton_Assist::bit_background - 1);
	if ((nix > two_power_n(Morton_Assist::bit_background - 1)) || (niy > two_power_n(Morton_Assist::bit_background - 1)))
	{
		std::stringstream error;
		error << "Index of x = " << nix << " or y = " << niy << " exceeds the limit of " << D_uint nmax << std::endl;
		log_error(error.str(), Log_function::logfile);
	}
	gr_NoIB.at(ilevel).grid.reserve(nx*ny);
#endif
#if (C_DIMS==3)
	nz = static_cast <D_uint> (C_zb / C_dx) + 1 + C_z0b_offset;
	D_uint niz = nz - 1;
	D_uint nmax = two_power_n(Morton_Assist::bit_background - 2);
	if ((nix > nmax) || (niy > nmax) || (niz > nmax))
	{
		std::stringstream error;
		error << "Index of x = " << nix << " or y = " << niy << " or z = " << niz << " exceeds the limit of " << nmax << std::endl;
		log_error(error.str(), Log_function::logfile);
	}
#if (C_MAP_TYPE == 1)
	grid_ptr->grid.reserve(nx*ny*nz);
#endif
#endif

	D_uint  n_max = nx - C_x0b_offset;
	if (n_max < ny) n_max = ny - C_y0b_offset;
#if(C_DIMS == 3)
	if (n_max < nz) n_max = nz - C_z0b_offset;
	unsigned int z0_offset_quotient = C_z0b_offset;
#endif
	D_uint n_max_temp = n_max;
	unsigned int cnt = 0;
	while (n_max_temp >>= 1) cnt++;
	if ((sizeof(D_uint) * CHAR_BIT) < (cnt * C_DIMS + Morton_Assist::bit_otherlevel))
	{
		std::stringstream error;
		error << "size of D_uint: "<< (sizeof(D_uint) * CHAR_BIT) <<" bits is smaller than that required in Morton encoding:  " << (cnt * C_DIMS + Morton_Assist::bit_otherlevel)<< " bits" << std::endl;
		log_error(error.str(), Log_function::logfile);
	}

	D_uint n_max_power = static_cast<D_uint>(pow(two_power_n(cnt), C_DIMS));
	D_uint int_key = 0;
	D_morton  key = 0;
	bool bool_x0 = true, bool_x1 = true, bool_y0 = true, bool_y1 = true;
	bool bool_z0 = true, bool_z1 = true;
	D_uint mortonx_temp, mortony_temp;
#if (C_DIMS == 2)
	D_uint mortonx_max_temp = Morton_Assist::pointer_me->morton_encode(nx - C_x0b_offset - 1, 0);
	D_uint mortony_max_temp = Morton_Assist::pointer_me->morton_encode(0, ny - C_y0b_offset - 1);
#endif
#if (C_DIMS == 3)
	D_uint mortonz_temp;
	D_uint mortonx_max_temp = (Morton_Assist::pointer_me->morton_encode(nx -C_x0b_offset - 1, 0, 0)).to_ullong();
	D_uint mortony_max_temp = (Morton_Assist::pointer_me->morton_encode(0, ny - C_y0b_offset - 1, 0)).to_ullong();
	D_uint mortonz_max_temp = (Morton_Assist::pointer_me->morton_encode(0, 0, nz - C_z0b_offset - 1)).to_ullong();
#endif

	for (D_uint  i_n = 0 ; i_n < n_max_power; ++i_n)
	{
		key = static_cast <D_morton>(int_key);
		mortonx_temp = (key & Morton_Assist::xmorton_ref.at(C_max_level).at(1)).to_ullong();
		mortony_temp = (key & Morton_Assist::ymorton_ref.at(C_max_level).at(1)).to_ullong();
#if (C_DIMS == 3)
		mortonz_temp = (key & Morton_Assist::zmorton_ref.at(C_max_level).at(1)).to_ullong();
#endif
		if (mortonx_temp < mortonx_max_temp)
		{
			bool_x0 = true;
			bool_x1 = true;
		}
		else if (mortonx_temp == mortonx_max_temp)
		{
			bool_x0 = true;
			bool_x1 = false;
		}
		else
		{
			bool_x0 = false;
			bool_x1 = false;
		}
		if (mortony_temp < mortony_max_temp)
		{
			bool_y0 = true;
			bool_y1 = true;
		}
		else if (mortony_temp == mortony_max_temp)
		{
			bool_y0 = true;
			bool_y1 = false;
		}
		else
		{
			bool_y0 = false;
			bool_y1 = false;
		}
#if (C_DIMS == 3)
		if (mortonz_temp < mortonz_max_temp)
		{
			bool_z0 = true;
			bool_z1 = true;
		}
		else if (mortonz_temp == mortonz_max_temp)
		{
			bool_z0 = true;
			bool_z1 = false;
		}
		else
		{
			bool_z0 = false;
			bool_z1 = false;
		}
#endif
		// (x0, y0, z0)
		if (bool_x0 && bool_y0 && bool_z0)
		{
			key = static_cast <D_morton>(int_key);
			for (unsigned int iter = 0; iter < C_x0b_offset; ++iter)
			{
				key = Morton_Assist::find_x1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_y0b_offset; ++iter)
			{
				key = Morton_Assist::find_y1(key, C_max_level);
			}
#if (C_DIMS == 3)
			for (unsigned int iter = 0; iter < C_z0b_offset; ++iter)
			{
				key = Morton_Assist::find_z1(key, C_max_level);
			}
#endif
			key = key << Morton_Assist::bit_otherlevel;
			Node temp;
			temp.flag = flag_refine;
			gr_NoIB.at(ilevel).grid.insert({ key, temp });
		}
		++int_key;

		// (x1, y0, z0)
		if (bool_x1 && bool_y0 && bool_z0)
		{
			key = static_cast <D_morton>(int_key);
			for (unsigned int iter = 0; iter < C_x0b_offset; ++iter)
			{
				key = Morton_Assist::find_x1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_y0b_offset; ++iter)
			{
				key = Morton_Assist::find_y1(key, C_max_level);
			}
#if (C_DIMS == 3)
			for (unsigned int iter = 0; iter < C_z0b_offset; ++iter)
			{
				key = Morton_Assist::find_z1(key, C_max_level);
			}
#endif
			key = key << Morton_Assist::bit_otherlevel;
			Node temp;
			temp.flag = flag_refine;
			gr_NoIB.at(ilevel).grid.insert({ key, temp });
		}
		++int_key;

		// (x0, y1, z0)
		if (bool_x0 && bool_y1 && bool_z0)
		{
			key = static_cast <D_morton>(int_key);
			for (unsigned int iter = 0; iter < C_x0b_offset; ++iter)
			{
				key = Morton_Assist::find_x1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_y0b_offset; ++iter)
			{
				key = Morton_Assist::find_y1(key, C_max_level);
			}
#if (C_DIMS == 3)
			for (unsigned int iter = 0; iter < C_z0b_offset; ++iter)
			{
				key = Morton_Assist::find_z1(key, C_max_level);
			}
#endif
			key = key << Morton_Assist::bit_otherlevel;
			Node temp;
			temp.flag = flag_refine;
			gr_NoIB.at(ilevel).grid.insert({ key, temp });
		}
		++int_key;

		// (x1, y1, z0)
		if (bool_x1 && bool_y1 && bool_z0)
		{
			key = static_cast <D_morton>(int_key);
			for (unsigned int iter = 0; iter < C_x0b_offset; ++iter)
			{
				key = Morton_Assist::find_x1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_y0b_offset; ++iter)
			{
				key = Morton_Assist::find_y1(key, C_max_level);
			}
#if (C_DIMS == 3)
			for (unsigned int iter = 0; iter < C_z0b_offset; ++iter)
			{
				key = Morton_Assist::find_z1(key, C_max_level);
			}
#endif
			key = key << Morton_Assist::bit_otherlevel;
			Node temp;
			temp.flag = flag_refine;
			gr_NoIB.at(ilevel).grid.insert({ key, temp });
		}
		++int_key;

#if(C_DIMS == 3)
		// (x0, y0, z1)
		if (bool_x0 && bool_y0 && bool_z1)
		{
			key = static_cast <D_morton>(int_key);
			for (unsigned int iter = 0; iter < C_x0b_offset; ++iter)
			{
				key = Morton_Assist::find_x1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_y0b_offset; ++iter)
			{
				key = Morton_Assist::find_y1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_z0b_offset; ++iter)
			{
				key = Morton_Assist::find_z1(key, C_max_level);
			}
			key = key << Morton_Assist::bit_otherlevel;
			Node temp;
			temp.flag = flag_refine;
			gr_NoIB.at(ilevel).grid.insert({ key, temp });
		}
		++int_key;

		// (x1, y0, z1)
		if (bool_x1 && bool_y0 && bool_z1)
		{
			key = static_cast <D_morton>(int_key);
			for (unsigned int iter = 0; iter < C_x0b_offset; ++iter)
			{
				key = Morton_Assist::find_x1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_y0b_offset; ++iter)
			{
				key = Morton_Assist::find_y1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_z0b_offset; ++iter)
			{
				key = Morton_Assist::find_z1(key, C_max_level);
			}
			key = key << Morton_Assist::bit_otherlevel;
			Node temp;
			temp.flag = flag_refine;
			gr_NoIB.at(ilevel).grid.insert({ key, temp });
		}
		++int_key;

		// (x0, y1, z1)
		if (bool_x0 && bool_y1 && bool_z1)
		{
			key = static_cast <D_morton>(int_key);
			for (unsigned int iter = 0; iter < C_x0b_offset; ++iter)
			{
				key = Morton_Assist::find_x1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_y0b_offset; ++iter)
			{
				key = Morton_Assist::find_y1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_z0b_offset; ++iter)
			{
				key = Morton_Assist::find_z1(key, C_max_level);
			}
			key = key << Morton_Assist::bit_otherlevel;
			Node temp;
			temp.flag = flag_refine;
			gr_NoIB.at(ilevel).grid.insert({ key, temp });
		}
		++int_key;

		// (x1, y1, z1)
		if (bool_x1 && bool_y1 && bool_z1)
		{
			key = static_cast <D_morton>(int_key);
			for (unsigned int iter = 0; iter < C_x0b_offset; ++iter)
			{
				key = Morton_Assist::find_x1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_y0b_offset; ++iter)
			{
				key = Morton_Assist::find_y1(key, C_max_level);
			}
			for (unsigned int iter = 0; iter < C_z0b_offset; ++iter)
			{
				key = Morton_Assist::find_z1(key, C_max_level);
			}
			key = key << Morton_Assist::bit_otherlevel;
			Node temp;
			temp.flag = flag_refine;
			gr_NoIB.at(ilevel).grid.insert({ key, temp });
		}
		++int_key;
#endif
	}

//#if(C_DIMS == 3)
//	for (D_uint iz = C_z0b_offset; iz < nz; ++iz)
//	{
//#endif
//		// Compute morton code for background nodes
//		for (D_uint iy = C_y0b_offset; iy < ny; ++iy)
//		{
//			for (D_uint ix = C_x0b_offset; ix < nx; ++ix)
//			{
//#if (C_DIMS==2)
//				D_morton key = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(ix, iy));
//#endif
//#if (C_DIMS==3)
//				D_morton key = static_cast <D_morton> (Morton_Assist::pointer_me->morton_encode(ix, iy, iz));
//#endif
//				key = key << Morton_Assist::bit_otherlevel;
//				Node temp;
//				temp.flag = flag_refine;
//				gr_NoIB.at(ilevel).grid.insert({ key, temp });
//			}
//		}
//#if (C_DIMS==3)
//	}
//#endif

	// set flag of nodes (ilevel = 0) on the numerical boundaries
	for (D_mapint::iterator iter = grid_ptr->coarse2fine.at(iboundary0_coarse).begin(); iter != grid_ptr->coarse2fine.at(iboundary0_coarse).end(); ++iter)
	{
		gr_NoIB.at(ilevel).grid[iter->first].flag = flag_iboundary[0];
	}

	// delete background nodes which are refined at higher level (i.e. ilevel>0) or are not on the numerical boundary of the background mesh
	for (unsigned int ishape = 0; ishape < Solid_Manager::pointer_me->numb_solids; ++ishape)
	{
		D_morton key_in;
//		D_uint numb_half = Solid_Manager::pointer_me->shape_solids.at(ishape).numb_nodes / 2;
//		D_real x = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(numb_half).x, y = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(numb_half).y;
//#if(C_DIMS == 3)
//		D_real z = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(numb_half).z;
//#endif
		D_real x = Solid_Manager::pointer_me->shape_solids.at(ishape).x0, y = Solid_Manager::pointer_me->shape_solids.at(ishape).y0;
#if(C_DIMS == 3)
		D_real z = Solid_Manager::pointer_me->shape_solids.at(ishape).z0;
#endif
		D_real x_index = x / grid_ptr->dx;
		D_real y_index = y / grid_ptr->dx;
		D_uint xint = static_cast<D_uint>(x_index);
		D_uint yint = static_cast<D_uint>(y_index);
#if (C_DIMS==2)
		key_in = static_cast <D_morton>(Morton_Assist::pointer_me->morton_encode(xint, yint));
#endif
#if (C_DIMS==3)
		D_real z_index = z / grid_ptr->dx;
		D_uint zint = static_cast<D_uint>(z_index);
		key_in = static_cast <D_morton>(Morton_Assist::pointer_me->morton_encode(xint, yint, zint));
#endif
		key_in = key_in << Morton_Assist::bit_otherlevel;
		flood_fill_background(ilevel, key_in);
	}

}


/**
* @brief function to delete background nodes which are refined at higher level (i.e. ilevel>0) or are not on the numerical boundary of the background mesh.
* @param[in]  ilevel    refinement level.
* @param[in]  key_in    morton code representing the node at refinement level > 0.
* @note Here param[in]:ilevel = 0. param[in]:key_in can't be a node on the numerical boundary between the block at refinement level = 0 and the block at refinement level= 1.
* @todo {using scan line flood fill method maybe faster, but encounter some problems in 3D cases, i.e., some nodes are not deleted}
*/
void Grid_Manager::flood_fill_background(unsigned int ilevel, D_morton &key_in)
{
	Grid * grid_ptr;
	grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel));

	D_morton morton_xyz = Morton_Assist::morton_xyz.at(ilevel);
	std::vector<D_morton> stk;
	stk.push_back(key_in);
	while (!stk.empty())
	{
		D_morton seed = stk.back();
		D_morton morton_temp;
		stk.pop_back();
		morton_temp = seed;
		if ((gr_NoIB.at(ilevel).grid.find(morton_temp) == gr_NoIB.at(ilevel).grid.end()) || (gr_NoIB.at(ilevel).grid[morton_temp].flag == flag_iboundary[0]))
		{
			continue;
		}
		else
		{
			gr_NoIB.at(ilevel).grid.erase(morton_temp);

			morton_temp = Morton_Assist::find_x0(seed, ilevel);
			stk.push_back(morton_temp);
			morton_temp = Morton_Assist::find_x1(seed, ilevel);
			stk.push_back(morton_temp);
			morton_temp = Morton_Assist::find_y0(seed, ilevel);
			stk.push_back(morton_temp);
			morton_temp = Morton_Assist::find_y1(seed, ilevel);
			stk.push_back(morton_temp);
#if (C_DIMS==3)
			morton_temp = Morton_Assist::find_z0(seed, ilevel);
			stk.push_back(morton_temp);
			morton_temp = Morton_Assist::find_z1(seed, ilevel);
			stk.push_back(morton_temp);
#endif
		}
	}
}

#if (C_SEARCH_METHOD != 1)
/**
* @brief function to intialize icount_refine of nodes of which digitals of the morton corresponding to x, y, z at current refinment level are 0.
* @param[in]  ilevel    refinement level.
* @param[in]  grid_pt    grid information at ilevel.
*/
template <class T_grid>
void Grid_Manager::initial_icount_refine(unsigned int ilevel, T_grid& grid_ptr)
{

	unsigned int extend_x0, extend_x1, extend_y0, extend_y1;
#if (C_DIMS == 3)
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

	unsigned int ilevel_1 = ilevel - 1;

	D_morton morton_xyz = Morton_Assist::morton_xyz.at(ilevel);
	D_morton morton_xyz_1 = Morton_Assist::morton_xyz_flip.at(ilevel);
	D_mapint map_node_ilevel_1 = {};
	D_morton morton_temp;
	
	if (ilevel == C_max_level)
	{
		for (D_map_Infor_near_solid::iterator iter = map_node_near_solid.begin(); iter != map_node_near_solid.end(); ++iter)
		{
			morton_temp = iter->first & morton_xyz_1;
			if (map_node_ilevel_1.find(morton_temp) == map_node_ilevel_1.end())
			{
				map_node_ilevel_1.insert({ morton_temp, 1 });
			}
			else
			{
				++map_node_ilevel_1[morton_temp];
			}
		}
	}
	else
	{
		for (D_mapint::iterator iter = grid_ptr.coarse2fine.at(iboundary2_coarse).begin(); iter != grid_ptr.coarse2fine.at(iboundary2_coarse).end(); ++iter)
		{
			morton_temp = iter->first & morton_xyz_1;			
			if (map_node_ilevel_1.find(morton_temp) == map_node_ilevel_1.end())
			{
				map_node_ilevel_1.insert({ morton_temp, 1 });
			}
			else
			{
				++map_node_ilevel_1[morton_temp];
			}
		}
	}


	D_morton morton_temp0, morton_temp1;
#if (C_DIMS==3)
	D_morton morton_temp2;
#endif
	unsigned int extend_temp_x0_1 = extend_x0 - 1, extend_temp_x1_1 = extend_x1 - 1, extend_temp_y0_1 = extend_y0 - 1, extend_temp_y1_1 = extend_y1 - 1;
	unsigned int extend_temp_x0 = extend_x0, extend_temp_x1 = extend_x1, extend_temp_y0 = extend_y0, extend_temp_y1 = extend_y1;
	D_morton morton_temp_x0, morton_temp_x1, morton_temp_y0, morton_temp_y1;
#if (C_DIMS==3)
	unsigned int extend_temp_z0 = extend_z0, extend_temp_z1 = extend_z1;
	unsigned int extend_temp_z0_1 = extend_z0 - 1, extend_temp_z1_1 = extend_z1 - 1;
	D_morton morton_temp_z0, morton_temp_z1;
#endif
	Timer tmr;
	double t0 = tmr.elapsed();
	Node_index node_index_temp;
	unsigned int bit_move = Morton_Assist::bit_otherlevel - ilevel * C_DIMS;
#if(C_SEARCH_METHOD == 2)
	// compute icount of nodes on the boundary of a given region
	for (D_mapint::iterator iter = map_node_ilevel_1.begin(); iter != map_node_ilevel_1.end(); ++iter)
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
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				else
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
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
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				else
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
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
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				else
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
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
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				else
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
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
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				else
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
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
					grid_ptr.grid_index[morton_temp].icount_refine += icount_temp;
				}
				else
				{
					grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
				}
				morton_temp = Morton_Assist::find_y1(morton_temp, ilevel_1);
			}
			morton_temp1 = Morton_Assist::find_z1(morton_temp1, ilevel_1);
		}
#endif
	}
#endif
#if(C_SEARCH_METHOD == 3)
	// compute icount of nodes within a given region

	// in this for loop, searching begins at the corner of the searching region and in x, y, and z direction
	for (D_mapint::iterator iter = map_node_ilevel_1.begin(); iter != map_node_ilevel_1.end(); ++iter)
	{
#if (C_CHECK_MORTON_BOUNDARY == 1)
		D_real x, y;
#if (C_DIMS == 2)
		Morton_Assist::compute_coordinate(iter->first, ilevel_1, x, y);
#endif
#if (C_DIMS == 3)
		D_real z;
		Morton_Assist::compute_coordinate(iter->first, ilevel_1, x, y, z);
#endif
		if ((x - extend_x0 * grid_ptr.dx) < (C_dx * C_x0b_offset))
		{
			extend_temp_x0 = static_cast<unsigned int>((x - C_dx * C_x0b_offset)/ grid_ptr.dx + C_eps);
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
					if (grid_ptr.grid_index.find(morton_temp) != grid_ptr.grid_index.end())
					{
						grid_ptr.grid_index[morton_temp].icount_refine += iter->second;
					}
					else
					{
						node_index_temp.icount_refine = iter->second;
						grid_ptr.grid_index.insert({ morton_temp, node_index_temp });
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

	std::cout << "size of grid_index at " << ilevel << " level:" << grid_ptr.grid_index.size() << std::endl;
	double t1 = tmr.elapsed();
	std::cout << ilevel << ", " << t1 - t0 << std::endl;
}
#endif

/**
* @brief function to find physical boundaries of the computational domain.
*/
void Grid_Manager::identify_domain_boundary()
{
	unsigned int ilevel = 0;
	Grid * grid_ptr;
	grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(ilevel));

#if (C_CHECK_MORTON_BOUNDARY == 1)
	// bounadries not in the background mesh
	// Mesh has GridIB member
	std::vector<D_morton> remove_temp;
	unsigned int ilevel_boundary_max = C_max_level;
	// x boundary
	for (D_mapNodePtr::iterator iter = bk_boundary_x.at(ilevel_boundary_max).at(0).begin(); iter != bk_boundary_x.at(ilevel_boundary_max).at(0).end(); ++iter)
	{
		if (gr_inner.grid.find(iter->first) != gr_inner.grid.end())
		{
			iter->second = &gr_inner.grid[iter->first];
		}
		else
		{
			remove_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
	{
		bk_boundary_x.at(ilevel_boundary_max).at(0).erase(*iter);
	}
	remove_temp.clear();
	for (D_mapNodePtr::iterator iter = bk_boundary_x.at(ilevel_boundary_max).at(1).begin(); iter != bk_boundary_x.at(ilevel_boundary_max).at(1).end(); ++iter)
	{

		if (gr_inner.grid.find(iter->first) != gr_inner.grid.end())
		{
			iter->second = &gr_inner.grid[iter->first];
		}
		else
		{
			remove_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
	{
		bk_boundary_x.at(ilevel_boundary_max).at(1).erase(*iter);
	}
	remove_temp.clear();
	// y boundary
	for (D_mapNodePtr::iterator iter = bk_boundary_y.at(ilevel_boundary_max).at(0).begin(); iter != bk_boundary_y.at(ilevel_boundary_max).at(0).end(); ++iter)
	{

		if (gr_inner.grid.find(iter->first) != gr_inner.grid.end())
		{
			iter->second = &gr_inner.grid[iter->first];
		}
		else
		{
			remove_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
	{
		bk_boundary_y.at(ilevel_boundary_max).at(0).erase(*iter);
	}
	remove_temp.clear();
	for (D_mapNodePtr::iterator iter = bk_boundary_y.at(ilevel_boundary_max).at(1).begin(); iter != bk_boundary_y.at(ilevel_boundary_max).at(1).end(); ++iter)
	{

		if (gr_inner.grid.find(iter->first) != gr_inner.grid.end())
		{
			iter->second = &gr_inner.grid[iter->first];
		}
		else
		{
			remove_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
	{
		bk_boundary_y.at(ilevel_boundary_max).at(1).erase(*iter);
	}
	remove_temp.clear();
#if (C_DIMS==3)
	// z boundary
	for (D_mapNodePtr::iterator iter = bk_boundary_z.at(ilevel_boundary_max).at(0).begin(); iter != bk_boundary_z.at(ilevel_boundary_max).at(0).end(); ++iter)
	{

		if (gr_inner.grid.find(iter->first) != gr_inner.grid.end())
		{
			iter->second = &gr_inner.grid[iter->first];
		}
		else
		{
			remove_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
	{
		bk_boundary_z.at(ilevel_boundary_max).at(0).erase(*iter);
	}
	remove_temp.clear();
	for (D_mapNodePtr::iterator iter = bk_boundary_z.at(ilevel_boundary_max).at(1).begin(); iter != bk_boundary_z.at(ilevel_boundary_max).at(1).end(); ++iter)
	{

		if (gr_inner.grid.find(iter->first) != gr_inner.grid.end())
		{
			iter->second = &gr_inner.grid[iter->first];
		}
		else
		{
			remove_temp.push_back(iter->first);
		}
	}
	for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
	{
		bk_boundary_z.at(ilevel_boundary_max).at(1).erase(*iter);
	}
	remove_temp.clear();
#endif
	// Mesh has with Grid member
	for (unsigned int ilevel_boundary = 1; ilevel_boundary < C_max_level; ++ilevel_boundary)
	{
		// x boundary
		for (D_mapNodePtr::iterator iter = bk_boundary_x.at(ilevel_boundary).at(0).begin(); iter != bk_boundary_x.at(ilevel_boundary).at(0).end(); ++iter)
		{

			if (gr_NoIB.at(ilevel_boundary).grid.find(iter->first) != gr_NoIB.at(ilevel_boundary).grid.end())
			{
				iter->second = &gr_NoIB.at(ilevel_boundary).grid[iter->first];
			}
			else
			{
				remove_temp.push_back(iter->first);
			}
		}
		for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
		{
			bk_boundary_x.at(ilevel_boundary).at(0).erase(*iter);
		}
		remove_temp.clear();
		for (D_mapNodePtr::iterator iter = bk_boundary_x.at(ilevel_boundary).at(1).begin(); iter != bk_boundary_x.at(ilevel_boundary).at(1).end(); ++iter)
		{

			if (gr_NoIB.at(ilevel_boundary).grid.find(iter->first) != gr_NoIB.at(ilevel_boundary).grid.end())
			{
				iter->second = &gr_NoIB.at(ilevel_boundary).grid[iter->first];
			}
			else
			{
				remove_temp.push_back(iter->first);
			}
		}
		for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
		{
			bk_boundary_x.at(ilevel_boundary).at(1).erase(*iter);
		}
		remove_temp.clear();
		// y boundary
		for (D_mapNodePtr::iterator iter = bk_boundary_y.at(ilevel_boundary).at(0).begin(); iter != bk_boundary_y.at(ilevel_boundary).at(0).end(); ++iter)
		{

			if (gr_NoIB.at(ilevel_boundary).grid.find(iter->first) != gr_NoIB.at(ilevel_boundary).grid.end())
			{
				iter->second = &gr_NoIB.at(ilevel_boundary).grid[iter->first];
			}
			else
			{
				remove_temp.push_back(iter->first);
			}
		}
		for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
		{
			bk_boundary_y.at(ilevel_boundary).at(0).erase(*iter);
		}
		remove_temp.clear();
		for (D_mapNodePtr::iterator iter = bk_boundary_y.at(ilevel_boundary).at(1).begin(); iter != bk_boundary_y.at(ilevel_boundary).at(1).end(); ++iter)
		{

			if (gr_NoIB.at(ilevel_boundary).grid.find(iter->first) != gr_NoIB.at(ilevel_boundary).grid.end())
			{
				iter->second = &gr_NoIB.at(ilevel_boundary).grid[iter->first];
			}
			else
			{
				remove_temp.push_back(iter->first);
			}
		}
		for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
		{
			bk_boundary_y.at(ilevel_boundary).at(1).erase(*iter);
		}
		remove_temp.clear();
#if (C_DIMS == 3)
		// z boundary
		for (D_mapNodePtr::iterator iter = bk_boundary_z.at(ilevel_boundary).at(0).begin(); iter != bk_boundary_z.at(ilevel_boundary).at(0).end(); ++iter)
		{

			if (gr_NoIB.at(ilevel_boundary).grid.find(iter->first) != gr_NoIB.at(ilevel_boundary).grid.end())
			{
				iter->second = &gr_NoIB.at(ilevel_boundary).grid[iter->first];
			}
			else
			{
				remove_temp.push_back(iter->first);
			}
		}
		for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
		{
			bk_boundary_z.at(ilevel_boundary).at(0).erase(*iter);
		}
		remove_temp.clear();
		for (D_mapNodePtr::iterator iter = bk_boundary_z.at(ilevel_boundary).at(1).begin(); iter != bk_boundary_z.at(ilevel_boundary).at(1).end(); ++iter)
		{

			if (gr_NoIB.at(ilevel_boundary).grid.find(iter->first) != gr_NoIB.at(ilevel_boundary).grid.end())
			{
				iter->second = &gr_NoIB.at(ilevel_boundary).grid[iter->first];
			}
			else
			{
				remove_temp.push_back(iter->first);
			}
		}
		for (std::vector<D_morton>::iterator iter = remove_temp.begin(); iter != remove_temp.end(); ++iter)
		{
			bk_boundary_z.at(ilevel_boundary).at(1).erase(*iter);
		}
		remove_temp.clear();
#endif
	}
#endif


	// domain boundary in the background mesh
	D_morton morton_temp0, morton_temp1;
	// x boundary
	morton_temp0 = Morton_Assist::mortonx_min;
	morton_temp1 = Morton_Assist::mortonx_max;
#if(C_DIMS == 3)
	D_morton morton_temp_pre0, morton_temp_pre1;
	morton_temp_pre0 = Morton_Assist::mortonx_min;
	morton_temp_pre1 = Morton_Assist::mortonx_max;
#endif
#if(C_DIMS == 3)
	for (unsigned int iz = 0; iz < nz; ++iz)
	{
		morton_temp0 = morton_temp_pre0;
		morton_temp1 = morton_temp_pre1;
#endif
		for (unsigned int iy = 0; iy < ny; ++iy)
		{
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (gr_NoIB.at(ilevel).grid.find(morton_temp0) != gr_NoIB.at(ilevel).grid.end())
			{
#endif		
				bk_boundary_x.at(0).at(0).insert({ morton_temp0 , &gr_NoIB.at(ilevel).grid[morton_temp0] });
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (gr_NoIB.at(ilevel).grid.find(morton_temp1) != gr_NoIB.at(ilevel).grid.end())
			{
#endif
				bk_boundary_x.at(0).at(1).insert({ morton_temp1 , &gr_NoIB.at(ilevel).grid[morton_temp1] });
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			morton_temp0 = Morton_Assist::find_y1(morton_temp0, ilevel);
			morton_temp1 = Morton_Assist::find_y1(morton_temp1, ilevel);
		}
#if(C_DIMS == 3)
		morton_temp_pre0 = Morton_Assist::find_z1(morton_temp_pre0, ilevel);
		morton_temp_pre1 = Morton_Assist::find_z1(morton_temp_pre1, ilevel);
	}
#endif
	// y boundary
	morton_temp0 = Morton_Assist::mortony_min;
	morton_temp1 = Morton_Assist::mortony_max;
#if(C_DIMS == 3)
	morton_temp_pre0 = Morton_Assist::mortony_min;
	morton_temp_pre1 = Morton_Assist::mortony_max;
#endif
#if(C_DIMS == 3)
	for (unsigned int iz = 0; iz < nz; ++iz)
	{
		morton_temp0 = morton_temp_pre0;
		morton_temp1 = morton_temp_pre1;
#endif
		for (unsigned int ix = 0; ix < nx; ++ix)
		{
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (gr_NoIB.at(ilevel).grid.find(morton_temp0) != gr_NoIB.at(ilevel).grid.end())
			{
#endif		
				bk_boundary_y.at(0).at(0).insert({ morton_temp0 , &gr_NoIB.at(ilevel).grid[morton_temp0] });
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (gr_NoIB.at(ilevel).grid.find(morton_temp1) != gr_NoIB.at(ilevel).grid.end())
			{
#endif		
				bk_boundary_y.at(0).at(1).insert({ morton_temp1 , &gr_NoIB.at(ilevel).grid[morton_temp1] });
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			morton_temp0 = Morton_Assist::find_x1(morton_temp0, ilevel);
			morton_temp1 = Morton_Assist::find_x1(morton_temp1, ilevel);
		}
#if(C_DIMS == 3)
		morton_temp_pre0 = Morton_Assist::find_z1(morton_temp_pre0, ilevel);
		morton_temp_pre1 = Morton_Assist::find_z1(morton_temp_pre1, ilevel);
	}
#endif


#if(C_DIMS == 3)
	// z boundary
	morton_temp0 = Morton_Assist::mortonz_min;
	morton_temp1 = Morton_Assist::mortonz_max;
	morton_temp_pre0 = Morton_Assist::mortonz_min;
	morton_temp_pre1 = Morton_Assist::mortonz_max;
	for (unsigned int iy = 0; iy < ny; ++iy)
	{
		morton_temp0 = morton_temp_pre0;
		morton_temp1 = morton_temp_pre1;
		for (unsigned int ix = 0; ix < nx; ++ix)
		{
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (gr_NoIB.at(ilevel).grid.find(morton_temp0) != gr_NoIB.at(ilevel).grid.end())
			{
#endif		
				bk_boundary_z.at(0).at(0).insert({ morton_temp0 , &gr_NoIB.at(ilevel).grid[morton_temp0] });
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
#if (C_CHECK_MORTON_BOUNDARY == 1)
			if (gr_NoIB.at(ilevel).grid.find(morton_temp1) != gr_NoIB.at(ilevel).grid.end())
			{
#endif		
				bk_boundary_z.at(0).at(1).insert({ morton_temp1 , &gr_NoIB.at(ilevel).grid[morton_temp1] });
#if (C_CHECK_MORTON_BOUNDARY == 1)
			}
#endif
			morton_temp0 = Morton_Assist::find_x1(morton_temp0, ilevel);
			morton_temp1 = Morton_Assist::find_x1(morton_temp1, ilevel);
		}
		morton_temp_pre0 = Morton_Assist::find_y1(morton_temp_pre0, ilevel);
		morton_temp_pre1 = Morton_Assist::find_y1(morton_temp_pre1, ilevel);
	}

#endif
}