/**
* @file
* @author Zhengliang Liu
* @brief Main function.
* @note  If shape.enclosed is true, to use flood fill to delte nodes inside the solid boundary, the line connecting the first solid node and the center of the solid should be inside the solid boundary
*/
#include "Shape.h"
#include <algorithm>
Shape& Shape::operator=(const Ini_Shape &c2)
{
	this->shape_type = c2.shape_type;
	this->bool_moving = c2.bool_moving;
	this->x0 = c2.x0;
	this->y0 = c2.y0;
#if (C_DIMS==3)
	this->z0 = c2.z0;
#endif
	this->numb_nodes = c2.numb_nodes;
	this->length = c2.length;
	return *this;
}

#if (C_DIMS == 2)
/**
* @brief      function to generate 2D circle.
* @param[in]  radius     radius of the circle.
* @param[out] node  points on the circle.
* @param[in]  t     time.
*/
void Shape::cycle(std::vector<Solid_Node> &node, D_real t)
{
	bool_enclosed = true;
	x0 = length.at(0);
	y0 = length.at(1);
	D_real radius = length.at(2);

	if ((2 * C_pi*radius / static_cast<D_real>(numb_nodes)) > (static_cast<D_real>(C_dx) / static_cast<D_real>(two_power_n(C_max_level))))
	{
		std::stringstream warning;
		warning << "the distance between two solid points of the cycle is greater than the grid space at the finest refinement level" << std::endl;
		log_warning(warning.str(), Log_function::logfile);
	}

	D_real radius_final = 2*radius - radius * sin(2. * C_pi*t/C_dx/80);

	for (D_uint i = 0; i < numb_nodes; ++i)
	{
		node.at(i).x = radius_final * cos(2. * C_pi * (static_cast<D_real>(i) / static_cast<D_real>(numb_nodes))) + x0 + shape_offest_x0;
		node.at(i).y = radius_final * sin(2. * C_pi * (static_cast<D_real>(i) / static_cast<D_real>(numb_nodes))) + y0 + shape_offest_y0;
#if (C_FSI_INTERFACE == 1)
		node.at(i).u = 0;
		node.at(i).v = 0;
		node.at(i).area = 2 * C_pi*radius_final / static_cast<D_real>(numb_nodes);
#endif
	}

}
#endif

/**
* @brief      function to generate 2D line.
* @param[in]  a     constant for y = ax + b.
* @param[in]  b     constant for y = ax + b.
* @param[out] node  points on the line.
* @param[in]  t     time.
* @note       it is assumed that the range of x is from 0 to C_xb. When in 3D, the z is set as C_zb / 2.
*/
void Shape::line_fillx(D_real a, D_real b, std::vector<Solid_Node> &node, D_real t)
{
	bool_enclosed = false;

	if ((C_xb / static_cast<D_real>(numb_nodes)) > (static_cast<D_real>(C_dx) / static_cast<D_real>(two_power_n(C_max_level))))
	{
		std::stringstream warning;
		warning << "the distance between two solid points of the line is greater than the grid space at the finest refinement level" << std::endl;
		log_warning(warning.str(), Log_function::logfile);
	}


    D_real dx = C_xb / static_cast<D_real>(numb_nodes);
	node.at(0).x = 0. + shape_offest_x0;
	node.at(0).y = b + shape_offest_y0;

#if (C_FSI_INTERFACE == 1)
	D_real arc = dx * sqrt(SQ(a) + 1);
	node.at(0).u = 0;
	node.at(0).v = 0;
	node.at(0).area = arc;
#endif

	D_uint end_i = 0;
	for (D_uint i = 1; i < numb_nodes; ++i)
	{
		node.at(i).x = node.at(i - 1).x + dx + shape_offest_x0;
		node.at(i).y = a * node.at(i).x + b + shape_offest_y0;
		if ((node.at(i).y) < 0)
		{
			node.at(i).x = -b/ a + shape_offest_x0;
			node.at(i).y = 0 + shape_offest_y0;
			end_i = i + 1;
			break;
		}
#if (C_FSI_INTERFACE == 1)
		node.at(i).u = 0;
		node.at(i).v = 0;
		node.at(i).area = arc;
#endif
	}

	if (end_i > 0)
	{
		for (D_uint i = end_i; i < numb_nodes; ++i)
		{
			node.at(i).x = node.at(0).x + shape_offest_x0;
			node.at(i).y = node.at(0).y + shape_offest_y0;
#if (C_FSI_INTERFACE == 1)
			node.at(i).u = 0;
			node.at(i).v = 0;
			node.at(i).area = arc;
#endif
		}
	}

#if(C_DIMS == 3)
	for (D_uint i = 0; i < numb_nodes; ++i)
	{
		node.at(i).z = (C_zb + shape_offest_z0) / 2;
	}
#endif

	if (node.at(numb_nodes - 1).y > (C_yb + shape_offest_y0))
	{
		std::stringstream warning;
		warning << "coordinate y = "<< node.at(numb_nodes - 1).y<< " of the last node exceeds the top compuational domain C_yb = "<<C_yb << std::endl;
		log_warning(warning.str(), Log_function::logfile);
	}

	x0 = node.at(numb_nodes / 2).x;
	y0 = node.at(numb_nodes / 2).y;

#if(C_DIMS==3)
	z0 = node.at(numb_nodes / 2).z;
#endif

}


#if (C_DIMS == 3)
/**
* @brief      function to generate 3D channel.
* @param[in]  a     constant for y = a*x + b, centerline of the channel.
* @param[in]  b     constant for y = a*x + b.
* @param[in]  c     constant for z = c*x + d.
* @param[in]  d     constant for z = c*x + d.
* @param[out] node  points on the line.
* @param[in]  t     time.
* @note       it is assumed that the range of x is from 0 to C_xb. Both y and z rely on x. Points are distributed with equal distance.
*/
void Shape::channel(D_real a, D_real b, D_real c, D_real d, D_real radius, std::vector<Solid_Node> &node, D_real t)
{

	bool_enclosed = true;

	D_real circ = 2 * C_pi * radius;  // circumference
	D_real dx = (circ + sqrt(SQ(circ) + 4 * static_cast<D_real>(numb_nodes)*C_xb *circ)) / 2 / static_cast<D_real>(numb_nodes);
	if ((static_cast<D_real>(C_dx) / static_cast<D_real>(two_power_n(C_max_level))) < dx)
	{
		std::stringstream warning;
		warning << "the distance (" << dx << ") between two adjacent solid points of the channel is greater than the grid space(" << (static_cast<D_real>(C_dx) / static_cast<D_real>(two_power_n(C_max_level))) << ") at the finest refinement level" << std::endl;
		log_warning(warning.str(), Log_function::logfile);
	}

	D_uint nx = static_cast<D_uint>(C_xb / dx + C_eps) + 1;
	D_uint ncirc = static_cast<D_uint>(circ / dx + C_eps);
	D_real xtemp0, ytemp0, ytemp1;
	D_real ztemp0, ztemp1;
	D_real theta;
	D_uint icount = 0;
	
	for (unsigned int ix = 0; ix < nx; ++ix)
	{
		xtemp0 = static_cast<D_real>(ix) * dx;
		ytemp0 = a * xtemp0 + b;
		ztemp0 = c * xtemp0 + d;
		for (unsigned int icirc = 0; icirc < ncirc; ++icirc)
		{
			theta = 2 * C_pi*static_cast<D_real>(icirc)/(ncirc + 1);
			ytemp1 = radius * cos(2 * theta);
			ztemp1 = radius * sin(2 * theta);
			node.at(icount).x = xtemp0 + shape_offest_x0;
			node.at(icount).y = ytemp0 + ytemp1 + shape_offest_y0;
			node.at(icount).z = ztemp0 + ztemp1 + shape_offest_z0;
#if (C_FSI_INTERFACE == 1)
			node.at(icount).u = 0;
			node.at(icount).v = 0;
			node.at(icount).area = dx;
#endif
			++icount;
		}
	}
	x0 = C_xb / 2;
	y0 = a * x0 + b;
	z0 = c * x0 + d;

	if (icount < numb_nodes)
	{
		node.resize(icount);
		std::stringstream warning;
		warning << "the number of nodes is not a multiplier of number of segments in x direction, set numb_nodes = "<< numb_nodes <<" as "<<icount << std::endl;
		log_warning(warning.str(), Log_function::logfile);
		numb_nodes = icount;
	}

}


#endif

/**
* @brief      function to read geometry data (cloud point).
* @param[in]  x0     geometry center, must inside the geometry for foll fill method.
* @param[in]  y0     geometry center.
* @param[in]  z0     geometry center.
* @param[out] node   infomration of solid points.
* @param[in]  t     time.
* @note       it is assumed that the range of x is from 0 to C_xb. Both y and z rely on x. Points are distributed with equal distance.
*/
void Shape::geofile(D_real xc, D_real yc, D_real zc, std::vector<Solid_Node> &node, D_real t)
{
	//bool_enclosed = true;
	D_real x_offset = shape_offest_x0, y_offset = shape_offest_y0;
#if (C_DIMS == 3)
	D_real z_offset = shape_offest_z0;
#endif

	std::istringstream istr;
	std::string s;
	std::ifstream file_in("./stl/airplane.txt", std::ios::in);
	if (!file_in.is_open())
	{
		std::stringstream error;
		char* buffer;
		if ((buffer = _getcwd(NULL, 0)) == NULL)
		{
			perror("getcwd error");
		}
		else
		{
			printf("%s\n", buffer);
			free(buffer);
		}
		error << "Can't open geometry file " << std::endl;
		log_error(error.str(), Log_function::logfile);
	}

	// check number of vertices and resise vector (nodes)
	numb_nodes = std::count(std::istreambuf_iterator<char>(file_in), std::istreambuf_iterator<char>(), '\n');;
	node.resize(numb_nodes);
	file_in.seekg(0, std::ios::beg);
	//file_in.seekg(0, std::ios::end);
	//std::streampos fp = file_in.tellg();
	//if (int(fp) == 0)
	//{
	//	std::cout << "test" << std::endl;
	//}

	D_real x_min = 0, y_min = 0;
#if (C_DIMS == 3)
	D_real z_min = 0;
#endif
	for (unsigned int i = 0; i < numb_nodes; ++i)
	{
		std::getline(file_in, s); istr.str(s);

		istr >> node.at(i).x >> node.at(i).y;
#if (C_DIMS == 3)
		istr >> node.at(i).z;
#endif
		if (node.at(i).x < x_min)
		{
			x_min = node.at(i).x;
		}

		if (node.at(i).y < y_min)
		{
			y_min = node.at(i).y;
		}
#if (C_DIMS == 3)
		if (node.at(i).z < z_min)
		{
			z_min = node.at(i).z;
		}
#endif
		istr.clear();
	}

	x_offset += xc;
	y_offset += yc;
#if (C_DIMS == 3)
	z_offset += zc;
#endif

	for (unsigned int i = 0; i < numb_nodes; ++i)
	{
		node.at(i).x += x_offset;
		node.at(i).y += y_offset;
#if (C_DIMS == 3)
		node.at(i).z += z_offset;
#endif
	}

	x0 = x_offset;
	y0 = y_offset;
#if (C_DIMS == 3)
	z0 = z_offset;
#endif

	file_in.close();

	//std::cout <<"center of the geometry: " << xc << ", " << yc << ", " << zc << std::endl;
		//std::cout <<"center of the geometry: " << xc << ", " << yc << ", " << zc << std::endl;
}

/**
* @brief      function to update information of solid nodes.
* @param[in] node  solid points.
* @param[in]  t     time.
* @note       it is assumed that the range of x is from 0 to C_xb. Both y and z rely on x. Points are distributed with equal distance.
*/
void Shape::geofile(std::vector<Solid_Node> &node, D_real t)
{
	x0+= C_dx / two_power_n(C_max_level) * 0.999;
	for (unsigned int i = 0; i < numb_nodes; ++i)
	{
		node.at(i).x += C_dx / two_power_n(C_max_level) * 0.999;
		node.at(i).y += 0;
#if (C_DIMS == 3)
		node.at(i).z += 0;
#endif
	}
}