/**
* @file
* @author Zhengliang Liu
* @brief Define solid manager class.
* @note .
*/

#ifndef SOLID_MANAGER_H
#define SOLID_MANAGER_H
#include "Shape.h"
class Solid_Manager
{
	friend class Obj_Manager;
public:
	static Solid_Manager* pointer_me;
	unsigned int numb_solids;           ///< total number of geometries
	std::vector<Shape> shape_solids;
	D_real shape_offest_x0_grid = 0, shape_offest_y0_grid = 0;
#if (C_DIMS==3)
	D_real shape_offest_z0_grid = 0;
#endif
private:
	void initial(const std::vector<Ini_Shape> &ini_shape);
	void renew(unsigned int ishape, D_real t);
};
#endif
