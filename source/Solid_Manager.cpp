/**
* @file
* @author Zhengliang Liu
* @brief  Function to manage solid.
* @note .
*/
#include "General.h"
#include "Solid_Manager.h"

/**
* @brief      function to intialize solid boundary.
* @param[in]  ini_shap    initial information of the solid shapes, see struct Ini_Shape.
*/
void Solid_Manager::initial(const std::vector<Ini_Shape> &ini_shap)
{
	numb_solids = ini_shap.size();
#if (C_CHECK_MORTON_BOUNDARY==1)
	shape_offest_x0_grid = C_dx * C_x0b_offset;
	shape_offest_y0_grid = C_dx * C_y0b_offset;
#if (C_DIMS==3)
	shape_offest_z0_grid = C_dx * C_z0b_offset;
#endif
#endif
	for (unsigned int i = 0; i < numb_solids; ++i)
	{
		Shape shape_temp;
		shape_temp = ini_shap.at(i);
		shape_solids.push_back(shape_temp);
		shape_solids.at(i) = ini_shap.at(i);
		shape_solids.at(i).shape_offest_x0 = shape_offest_x0_grid;
		shape_solids.at(i).shape_offest_y0 = shape_offest_y0_grid;
#if (C_DIMS==3)
		shape_solids.at(i).shape_offest_z0 = shape_offest_z0_grid;
#endif
		shape_solids.at(i).node = std::vector<Solid_Node>(shape_solids.at(i).numb_nodes);
		switch (ini_shap.at(i).shape_type)
		{
		case circle:
#if (C_DIMS == 2)		
			if (shape_temp.length.size() != 3)
			{
				std::stringstream error;
				error << "the number of input character lengths should be 1 for circle type, input is: (radius)" << std::endl;
				log_error(error.str(), Log_function::logfile);
			}
			else
			{
				shape_solids.at(i).cycle(shape_solids.at(i).node, 0);
			}
#endif				
			break;
		case	line_fillx:
			if (shape_temp.length.size() != 2)
			{
				std::stringstream error;
				error << "the number of input character lengths should be 2 for line_fillx type, inputs are: (a) and (b) " << std::endl;
				log_error(error.str(), Log_function::logfile);
			}
			else
			{
				shape_solids.at(i).line_fillx(shape_temp.length.at(0), shape_temp.length.at(1), shape_solids.at(i).node, 0);
			}
			break;
#if(C_DIMS == 3)
		case channel:
			if (shape_temp.length.size() != 5)
			{
				std::stringstream error;
				error << "the number of input character lengths should be 5 for line_fillx type, inputs are: (a), (b), (c), (d)  and radius" << std::endl;
				log_error(error.str(), Log_function::logfile);
			}
			else
			{
				shape_solids.at(i).channel(shape_temp.length.at(0), shape_temp.length.at(1), shape_temp.length.at(2), shape_temp.length.at(3), shape_temp.length.at(4), shape_solids.at(i).node, 0);
			}
			break;
#endif
		case geofile:
			shape_solids.at(i).geofile(shape_temp.length.at(0), shape_temp.length.at(1), shape_temp.length.at(2), shape_solids.at(i).node, 0);
			break;
		default:
			std::stringstream error;
			error << "input solid type does not match any of built in types" << std::endl;
			log_error(error.str(), Log_function::logfile);
			break;
		}

#if (C_FSI_INTERFACE == 1)
		for (unsigned int inode = 0; inode < shape_solids.at(i).numb_nodes; ++inode)
		{
			shape_solids.at(i).node.at(inode).xref = shape_solids.at(i).node.at(inode).x;
			shape_solids.at(i).node.at(inode).yref = shape_solids.at(i).node.at(inode).y;
#if(C_DIMS == 3)
			shape_solids.at(i).node.at(inode).zref = shape_solids.at(i).node.at(inode).z;
#endif
		}
#endif
		//if (ini_shap.at(i).shape_type == circle)
		//{
		//	
		//}
		//else if (ini_shap.at(i).shape_type == line_fillx)
		//{
		//	shape_solids.at(i).line_fillx(shape_temp.length.at(0), shape_temp.length.at(1), shape_solids.at(i).node, 0);
		//}

		
	}
}

/**
* @brief      function to update solid boundary.
* @param[in]  ishape    numbering of the solid boundary.
* @param[in]  t         time.
*/
void Solid_Manager::renew(unsigned int ishape, D_real t)
{	
	switch (shape_solids.at(ishape).shape_type)
	{
	case circle:
#if (C_DIMS == 2)
		shape_solids.at(ishape).cycle(shape_solids.at(ishape).node, t);
#endif
		break;
	case geofile:
		shape_solids.at(ishape).geofile(shape_solids.at(ishape).node, t);
		break;
	default:
		break;
	}
	
}

