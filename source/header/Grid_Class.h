/**
* @file
* @author Zhengliang Liu
* @brief Define grid class.
*/

#ifndef GRID_CLASS_H
#define GRID_CLASS_H
#include "General.h"

#if (C_MAP_TYPE == 1)
#include <unordered_map>
#elif(C_MAP_TYPE == 2)
#include <map>
#endif

/**
* @brief This structure will be used to store node information used to add and remove node
*/
struct Node_index
{
#if (C_SOLID_BOUNDARY==2)
	D_int icount_refine = 0;
#endif
};

/**
* @brief This structure will be used to store node information related to LBM 
*/
struct Node
{
public:
	// flag
	int flag = 0;                       ///< flag to categrary nodes
	// macroscopic and mesoscopic variables
//	std::array<D_real, C_Q> f;		    ///< Distribution functions
//	std::array<D_real, C_Q> f_collid;   ///< Distribution functions after collision
//	D_real rho = 1;			                ///< Macropscopic density
//	D_real rho_aver = 1;	                ///< Time-avaraged density
//	D_real u = 0, v = 0;			            ///< Macropscopic velocities
//	D_real u_aver = 0, v_aver = 0;	            ///< Time-avaraged velocities
//#if (C_DIMS==3)
//	D_real w = 0;			                ///< Macropscopic velocities
//	D_real w_aver = 0;	                    ///< Time-avaraged velocities
//#endif
public:
	void initial(unsigned int ilevel);
	Node operator+(const Node& node_R)
	{
		Node node_return;
		node_return.flag = this->flag + node_R.flag;
		return node_return;
	}
	Node& operator+=(const Node& node_R)
	{
		this->flag += node_R.flag;
		return *this;
	}
	Node operator*(D_real real_R)
	{
		Node node_return;
		node_return.flag = this->flag * real_R;
		return node_return;
	}
	Node operator-(const Node& node_R)
	{
		Node node_return;
		node_return.flag = this->flag - node_R.flag;
		return node_return;
	}
};

class Node_IB : public Node
{
public:
	D_int icount_ghost = 0;                       ///< number of solid points witchin a given distance at (ilevel - 1)
};

/**
* @brief This structure will be used to store node information near the solid
*/
struct Infor_near_solid
{
public:
	D_int icount;  ///< number of solid points witchin a given distance at ilevel
};

#if (C_FSI_INTERFACE == 1)
/**
* @brief This structure will be used to store node information related to immersed boundary method
*/
struct Infor_IB
{
public:
	D_int icount;  ///< number of solid points witchin a given distance at ilevel
	D_real fx, fy;  ///< body force
#if (C_DIMS == 3)
	D_real fz;
#endif
};

#endif

/**
* @brief The basic class represents mesh based on STL container
*/
#if (C_MAP_TYPE == 1)
template <typename value>
using D_map_define = std::unordered_map<D_morton, value>;
#elif(C_MAP_TYPE == 2)

class compare
{
public:
	bool operator()(const D_morton &a, const D_morton &b) const
	{
		return a.to_ullong() < b.to_ullong();
	}
};
template <typename value>
using D_map_define = std::map<D_morton, value, compare>;
#endif

using D_map_Infor_near_solid = D_map_define<Infor_near_solid>;
#if (C_FSI_INTERFACE == 1)
using D_map_Infor_IB = D_map_define<Infor_IB>;
#endif
using D_map_index = D_map_define<Node_index>; // Using unordered map to store mesh
using D_map = D_map_define<Node>; // Using unordered map to store mesh
using D_mapIB = D_map_define<Node_IB>; // Using unordered map to store mesh
using D_mapint = D_map_define<int>;
using D_mapint2 = D_map_define<std::array<int, 2>>;
using D_mapNodePtr = D_map_define<Node*>;


//struct Morton_Ptr
//{
//public:
//	Node* ptr = nullptr;
//	D_morton morton;
//};

class Grid
{
	friend class Grid_Manager;
	friend class Morton_Assist;
private:
	D_map_index grid_index;
protected:
	D_real dx;                       ///< grid space in the background mesh
//	array<D_real, C_DIMS> center;    ///< Center of the block
public:
	std::array <D_mapint, C_overlap * 2 + 1> fine2coarse;         ///< interface of the boundary layer overlapping coarser mesh
	std::array <D_mapint, C_overlap + 1> coarse2fine;         ///< interface of the boundary layer overlapping finner mesh
};

class Grid_NIB : public Grid
{
public:
	D_map grid;                ///< Map to store morton code of the grid
};

class Grid_IB : public Grid
{
public:
	D_mapIB grid;                ///< Map to store morton code of the grid
};



///**
//* @brief This class represents back ground mesh (refine_level is 0)
//*/
//class Grid_Background : public Grid
//{
//public:
//	D_map grid;                ///< Map to store morton code of the grid
//	D_int nx;                  ///< maximum number of grid points in the x direction  
//	D_int ny;                  ///< maximum number of grid points in the y direction 
//#if (C_DIMS==3)
//	D_int nz;                  ///< maximum number of grid points in the z direction 
//#endif
//};

#endif