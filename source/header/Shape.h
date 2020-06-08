/**
* @file
* @author Zhengliang Liu
* @brief Define objective manager class.
*/

#ifndef SHAPE_H
#define SHAPE_H
#include "General.h"
struct Solid_Node
{
public:
	D_real x, y;
#if (C_FSI_INTERFACE == 1)
	D_real xref, yref; ///< reference coordinate at previous time step, to check if the solid point move to adjacent cells at C_max_level
	D_real xref2, yref2; ///< reference coordinate at previous time step, to check if the solid point move to adjacent cells at (C_max_level - 1)
	D_real area;     ///< area near the node, used to calculate IB force
	D_real u, v;    /// velocity of the solid point
#endif
#if (C_DIMS==3)
	D_real z;
#if (C_FSI_INTERFACE == 1)
	D_real zref, zref2, w;
#endif
#endif
};

/**
* @brief This enumerate will be used to classify different solid shapes.
*/
enum Shape_enum
{
	circle = 1, // 2D
	line_fillx = 2,  //2D and 3D
	channel = 3, // 3D
	geofile = 4, // read from point cloud data (2D and 3D)
};


/**
* @brief This structure used to store intial inforamtion of solid.
*/
struct Ini_Shape
{
	Shape_enum shape_type;
public:
	bool bool_moving = false;
	bool bool_enclosed = false;    ///< when set as true, it's a enclosed shape
	D_real x0 = 0, y0 = 0;    ///< center of each solid
	D_uint numb_nodes = 0; ///< number of solid nodes
#if (C_DIMS==3)
	D_real z0;
#endif
	std::vector<D_real> length; ///< used to store character length of the solid. I.e. cycle: length[0] is the radius
};


/**
* @brief This class used to store all inforamtion of solid nodes.
*/
class Shape : public Ini_Shape
{
	friend class Solid_Manager;
public:
	std::vector<Solid_Node> node;
	D_real shape_offest_x0 = 0, shape_offest_y0 = 0;
#if (C_DIMS==3)
	D_real shape_offest_z0 = 0;
#endif
private:
	Shape& operator=(const Ini_Shape &c1);

	void cycle(std::vector<Solid_Node> &node, D_real t);
	void line_fillx(D_real a, D_real b, std::vector<Solid_Node> &node, D_real t);
	void channel(D_real a, D_real b, D_real c, D_real d, D_real radius, std::vector<Solid_Node> &node, D_real t);
	void geofile(D_real x0, D_real y0, D_real z0, std::vector<Solid_Node> &node, D_real t);
	void geofile(std::vector<Solid_Node> &node, D_real t);
	void sphere(D_real a, D_real b, D_real c, D_real d, D_real radius, std::vector<Solid_Node> &node, D_real t);
};
#endif
