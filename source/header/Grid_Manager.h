/**
* @file
* @author Zhengliang Liu
* @brief Define grid manager class.
* @note .
*/
#ifndef GRID_MANAGER_H
#define GRID_MANAGER_H
#include "General.h"
#include "Grid_Class.h"
#include "Morton_assist.h"
/**
* @brief This class used to manager all grid related information.
*/
class Grid_Manager
{
	friend class Morton_Assist;
	friend class Obj_Manager;
	friend class LBM_Manager;
public:
	static Grid_Manager* pointer_me;         ///< pointer points to the class itself
	D_real xb_domain, yb_domain;
	D_uint nx, ny; ///< maximum number of nodes in the x and y direction of the background mesh
	std::array< std::array<D_mapNodePtr, 2>, C_max_level + 1> bk_boundary_x, bk_boundary_y; ///< domain boundary, i.e. bk_boundary_x[0] is the boudnary where x=0
#if (C_DIMS == 3)
	D_uint nz;
	D_real zb_domain;
	std::array< std::array<D_mapNodePtr, 2>, C_max_level + 1>  bk_boundary_z;
#endif
	std::array<Grid_NIB, C_max_level> gr_NoIB;
	Grid_IB gr_inner;

	D_map_Infor_near_solid map_node_near_solid;
#if (C_FSI_INTERFACE == 1)
	D_map_Infor_IB map_node_IB;
#endif

private:
	bool bool_any_enclosed = false;
	const std::array<int, 3> flag_iboundary = {-1, -2, -3}; // nodes on the numerical boundary;
	const int flag_refine = two_power_n(1); // nodes need to refine;
	const int flag_near_solid = 15;//two_power_n(11); // nodes within one grid space to the solid node
	const int flag_ghost = two_power_n(2); // nodes within C_extend_ghost grid space to the solid node
	const int flag_inner = two_power_n(3); //  nodes inside the solid
#if (C_CHECK_MORTON_BOUNDARY == 1)
	const int flag4 = two_power_n(4); //  nodes on the negative x boundary of the computational domain
	const int flag5 = two_power_n(5); //  nodes on the positive x boundary of the computational domain
	const int flag6 = two_power_n(6); //  nodes on the negative y boundary of the computational domain
	const int flag7 = two_power_n(7); //  nodes on the positive y boundary of the computational domain
#if (C_DIMS == 3)
	const int flag8 = two_power_n(8); //  nodes on the negative z boundary of the computational domain
	const int flag9 = two_power_n(9); //  nodes on the positive z boundary of the computational domain
#endif
#endif

// the last six bit of the int value in map_outer is used to flag if there's a node on the inner layer adjacent to the current node in the x, y or z direction
// the bit represents (z1 z0 y1 y0 x1 x0)
	const int bit_x0 = two_power_n(0);
	const int bit_x1 = two_power_n(1);
	const int bit_y0 = two_power_n(2);
	const int bit_y1 = two_power_n(3);
#if(C_DIMS==3)
	const int bit_z0 = two_power_n(4);
	const int bit_z1 = two_power_n(5);
	const int bit_z01 = bit_z0 | bit_z1;
#endif 
	const int bit_x01 = bit_x0 | bit_x1;
	const int bit_y01 = bit_y0 | bit_y1;

	unsigned int ioverlap = 1;  // the order of the overlapped boudary layer from inner to outer at the coarse level (ilevel_1)
	unsigned int iboundary0 = (ioverlap - 1) * 2; // the inner interface of the current boudnary layer at the fine level (ilevel)
	unsigned int iboundary1 = (ioverlap - 1) * 2 + 1; // the middle interface of the current boudnary layer at the fine level (ilevel)
	unsigned int iboundary2 = ioverlap * 2; // the outer interface of the current boudnary layer at the fine level (ilevel)
	unsigned int iboundary0_coarse = ioverlap - 1; // the interface at the coarse level (ilevel - 1) corresponding to iboundary0 
	unsigned int iboundary2_coarse = ioverlap; // the interface at the coarse level (ilevel - 1) corresponding to iboundary2 

private:
	void initial();

	void generate_inner();
	void find_nodes_near_solid(unsigned int ilevel, std::array<D_real, C_DIMS> xyz, D_mapint2 &refine_nodes);
#if	(C_FSI_INTERFACE == 1)
	void initial_map_node_IB(unsigned int ilevel);
#endif
	void mark_ghost_points(unsigned int ilevel, D_mapint2& refine_nodes);
	void search_nodes_near_solid(unsigned int ilevel, D_mapint &nearest_center, D_mapint &nodes_level1, D_mapint2  &refine_nodes, std::array<D_mapint, 2> &boundary_x_temp, std::array<D_mapint, 2> &boundary_y_temp, std::array<D_mapint, 2> &boundary_z_temp);
	void flood_fill_inner(unsigned int ilevel, D_morton &key_in, D_mapint2 &nodes_level1, D_mapint2  &refine_nodes, std::array<D_mapint, 2> &boundary_x_temp, std::array<D_mapint, 2> &boundary_y_temp, std::array<D_mapint, 2> &boundary_z_temp);

	void search_outer_boundary(unsigned int ilevel, D_mapint2 &iter_nodes, D_mapint2 &exist_nodes, D_mapint2 &refine_nodes_temp, D_mapint &map_outer, std::array<D_mapint, 2> &boundary_x_temp, std::array<D_mapint, 2> &boundary_y_temp, std::array<D_mapint, 2> &boundary_z_temp);
//#if(C_DIMS==3)
//	void check_isolate_nodes(unsigned int ilevel);
//#endif
	void check_merge_nodes(unsigned int ilevel, D_mapint2 &node_exist, D_mapint &map_outer);
	template<typename  T> void search_inter_boundary(unsigned int ilevel, unsigned int iboundary1, unsigned int iboundary2, unsigned int iboundary2_coarse, D_mapint map_outer, T refine_nodes, std::array<D_mapint, 2> &boundary_x_temp, std::array<D_mapint, 2> &boundary_y_temp, std::array<D_mapint, 2> &boundary_z_temp);

	void generate_intermediate(unsigned int ilevel);

	void generate_background();
	void flood_fill_background(unsigned int ilevel, D_morton &key_in);
	void identify_domain_boundary();

#if (C_SEARCH_METHOD != 1)
	template <class T_grid>
	void initial_icount_refine(unsigned int ilevel, T_grid& grid_ptr);
#endif

#if	(C_SOLID_BOUNDARY == 2)
	void initial_ref_coordinate(unsigned int ishape, D_uint ipoint);
	void update_nodes_near_solid(unsigned int ishape, D_mapint &map_add_nodes, D_mapint &map_remove_nodes);
	void call_update_nodes(const unsigned int ilevel, D_mapint &map_add_nodes_in, D_mapint &map_remove_nodes_in, D_mapint &map_add_nodes_out, D_mapint &map_remove_nodes_out);
#if	(C_FSI_INTERFACE == 1)
	void update_map_node_IB(const unsigned int ilevel, D_mapint &map_add_nodes_in, D_mapint &map_remove_nodes_in);
#endif
	void update_ghost_node(const unsigned int ilevel, D_mapint &map_add_nodes_in, D_mapint& map_remove_nodes_in);
	template <class T_grid, class T_node>
	void update_nodes(const unsigned int ilevel, T_grid &grid_ptr, T_node node_temp, D_mapint &map_add_nodes_in,  D_mapint &map_remove_nodes_in, D_mapint &map_add_nodes_out, D_mapint &map_remove_nodes_out);

	template <class T_grid, class T_node>
	void reconstruct_numerical_boundary_removing(const unsigned int ilevel, T_grid &grid_ptr, T_node node_temp, D_mapint &map_remove_temp, D_mapint &map_add_nodes_out, D_mapint &map_remove_nodes_out);
	//void reconstruct_numerical_boundary_removing(const unsigned int ilevel, Grid_NIB &grid_ptr, Node node_temp, D_mapint &map_remove_temp, D_mapint &map_remove_nodes_out);

	template <class T_grid, class T_node>
	void reconstruct_numerical_boundary_adding(const unsigned int ilevel, T_grid &grid_ptr, T_node node_temp, D_mapint &map_add_temp, D_mapint &map_add_nodes_out, D_mapint &map_remove_nodes_out);
	//void reconstruct_numerical_boundary_adding(const unsigned int ilevel, Grid_NIB &grid_ptr, Node node_temp, D_mapint &map_add_temp, D_mapint &map_add_nodes_out);
#if (C_DEBUG == 1)
public:
	D_mapint output_points_update;
#endif
#endif

};
#endif