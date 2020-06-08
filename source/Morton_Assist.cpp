/**
* @file
* @author Zhengliang Liu
* @brief Functions related to morton code.
* @note .
*/
// (x, y, z) is coded as zyx
#include "General.h"
#include "Grid_Manager.h"
#include "Morton_assist.h"
unsigned int Morton_Assist::bit_background;
unsigned int Morton_Assist::bit_otherlevel;
D_morton Morton_Assist::xbk_ref0, Morton_Assist::xbk_ref1, Morton_Assist::ybk_ref0, Morton_Assist::ybk_ref1;
#if (C_DIMS==3)
D_morton Morton_Assist::zbk_ref0, Morton_Assist::zbk_ref1;
#endif

/**
* @brief function to compute reference values used in Morton code calculation.
*/
void Morton_Assist::morton_initial()
{
	Morton_Assist::bit_background = C_BIT - C_max_level * C_DIMS;
	Morton_Assist::bit_otherlevel = C_max_level * C_DIMS;
	//bk_refine = 0;
	//for (unsigned int i = 0; i < Morton_Assist::bit_otherlevel; ++i)
	//{
	//	bk_refine = bk_refine.set(i, true);
	//}
	xbk_ref1 = 0;
	ybk_ref1 = 0;
#if (C_DIMS==3)
	zbk_ref1 = 0;
#endif
	for (unsigned int i = 0; i < Morton_Assist::bit_background / C_DIMS; ++i)
	{
		unsigned int pos = i * C_DIMS + Morton_Assist::bit_otherlevel;
		xbk_ref1.set(pos, true);
		ybk_ref1.set(pos + 1, true);
#if (C_DIMS==3)
		zbk_ref1.set(pos + 2, true);
#endif
	}
	xbk_ref0 = ~xbk_ref1;
	ybk_ref0 = ~ybk_ref1;
#if (C_DIMS==3)
	zbk_ref0 = ~zbk_ref1;
#endif
	for (unsigned int i = 0; i < Morton_Assist::bit_otherlevel; ++i)
	{
		xbk_ref0.set(i, false);
		ybk_ref0.set(i, false);
#if (C_DIMS==3)
		zbk_ref0.set(i, false);
#endif
	}
}

/**
* @brief function compute reference mortorn code.
* @param[in]  ilevel   refinement level.
*/
void Morton_Assist::compute_ref(unsigned int ilevel)
{
	morton_xyz.at(ilevel) = 0;
	morton_xyz.at(ilevel).set(Morton_Assist::bit_otherlevel - ilevel * C_DIMS, true);
	morton_xyz.at(ilevel).set(Morton_Assist::bit_otherlevel - ilevel * C_DIMS + static_cast<long long>(1), true);
#if (C_DIMS==3)
	morton_xyz.at(ilevel).set(Morton_Assist::bit_otherlevel - ilevel * C_DIMS + static_cast<long long>(2), true);
#endif
	morton_xyz_flip.at(ilevel) = morton_xyz.at(ilevel);
	morton_xyz_flip.at(ilevel).flip();

	ref_one.at(ilevel) = 0;
	ref_one.at(ilevel).set(Morton_Assist::bit_otherlevel - ilevel * C_DIMS, true); // reference value which is 1 at the sepecified position for a given level

	unsigned int pos;
	xmorton_ref.at(ilevel).at(1) = xbk_ref1;
	ymorton_ref.at(ilevel).at(1) = ybk_ref1;
#if (C_DIMS==3)
	zmorton_ref.at(ilevel).at(1) = zbk_ref1;
#endif
	for (unsigned int i = 0; i < ilevel + 1; ++i)
	{
		pos = Morton_Assist::bit_otherlevel - i * C_DIMS;
		xmorton_ref.at(ilevel).at(1).set(pos, true);
		pos = Morton_Assist::bit_otherlevel - i * C_DIMS + 1;
		ymorton_ref.at(ilevel).at(1).set(pos, true);
#if (C_DIMS==3)
		pos = Morton_Assist::bit_otherlevel - i * C_DIMS + 2;
		zmorton_ref.at(ilevel).at(1).set(pos, true);
#endif
	}
	xmorton_ref.at(ilevel).at(0) = ~xmorton_ref.at(ilevel).at(1);
	ymorton_ref.at(ilevel).at(0) = ~ymorton_ref.at(ilevel).at(1);
#if (C_DIMS==3)
	zmorton_ref.at(ilevel).at(0) = ~zmorton_ref.at(ilevel).at(1);
#endif
	for (unsigned int i = 0; i < Morton_Assist::bit_otherlevel - ilevel * C_DIMS; ++i)
	{
		xmorton_ref.at(ilevel).at(0).set(i, false);
		ymorton_ref.at(ilevel).at(0).set(i, false);
#if (C_DIMS==3)
		zmorton_ref.at(ilevel).at(0).set(i, false);
#endif
	}
}


#if (C_DIMS==2)
/**
* @brief function for Morton edcoding.
* @param[in]  ix   index in x direction.
* @param[in]  iy   index in y direction.
*/
D_morton Morton_Assist::morton_encode(D_uint ix, D_uint iy)
{
	D_morton key_int = 0;
	for (unsigned int i = 0; i < (C_BIT / C_DIMS); ++i) {
		key_int |= ((ix & ((static_cast <D_uint> (1)) << i)) << i) | ((iy & ((static_cast <D_uint> (1)) << i)) << (i + 1));
	}
	return key_int;
}

/**
* @brief function to compute the coordinates from Morton code.
* @param[in]  key      Morton code of current node.
* @param[in]  ilevel   refinement level.
* @param[out] x        x coordinate.
* @param[out] y        y coordinate.
*/
void Morton_Assist::compute_coordinate(D_morton key, unsigned int ilevel, D_real &x, D_real &y)
{
	D_real x_temp = 0.;
	D_real y_temp = 0.;
	for (unsigned int i = 1; i < ilevel + 1; ++i)
	{
		if (key.test(Morton_Assist::bit_otherlevel - i * C_DIMS))
		{
			x_temp += C_dx / two_power_n(i);
		}
		if (key.test(Morton_Assist::bit_otherlevel - i * C_DIMS + 1))
		{
			y_temp += C_dx / two_power_n(i);
		}

	}
	for (unsigned int i = 0; i < Morton_Assist::bit_background / C_DIMS; ++i)
	{
		if (key.test(Morton_Assist::bit_otherlevel + i * C_DIMS))
		{
			x_temp += C_dx * two_power_n(i);
		}
		if (key.test(Morton_Assist::bit_otherlevel + i * C_DIMS + 1))
		{
			y_temp += C_dx * two_power_n(i);
		}
	}
	x = x_temp;
	y = y_temp;
};
#endif
#if (C_DIMS==3)
D_morton Morton_Assist::morton_encode(D_uint ix, D_uint iy, D_uint iz)
{
	D_morton key_int = 0;
	for (unsigned int i = 0; i < (C_BIT / C_DIMS); ++i)
	{
		key_int |= ((ix & ((static_cast <D_uint> (1)) << i)) << (2 * i)) | ((iy & ((static_cast <D_uint> (1)) << i)) << (2 * i + 1)) | ((iz & ((static_cast <D_uint> (1)) << i)) << (2 * i + 2));
	}
	return key_int;
}

/**
* @brief function to compute the coordinates from Morton code.
* @param[in]  key      Morton code of current node.
* @param[in]  ilevel   refinement level.
* @param[out] x        x coordinate.
* @param[out] y        y coordinate.
* @param[out] z        z coordinate.
*/
void Morton_Assist::compute_coordinate(D_morton key, unsigned int ilevel, D_real &x, D_real &y, D_real &z)
{
	D_real x_temp = 0.;
	D_real y_temp = 0.;
	D_real z_temp = 0.;
	for (unsigned int i = 1; i < ilevel + 1; ++i)
	{
		if (key.test(Morton_Assist::bit_otherlevel - i * C_DIMS))
		{
			x_temp += C_dx / two_power_n(i);
		}
		if (key.test(Morton_Assist::bit_otherlevel - i * C_DIMS + 1))
		{
			y_temp += C_dx / two_power_n(i);
		}
		if (key.test(Morton_Assist::bit_otherlevel - i * C_DIMS + 2))
		{
			z_temp += C_dx / two_power_n(i);
		}
	}
	for (unsigned int i = 0; i < Morton_Assist::bit_background / C_DIMS; ++i)
	{
		if (key.test(Morton_Assist::bit_otherlevel + i * C_DIMS))
		{
			x_temp += C_dx * two_power_n(i);
		}
		if (key.test(Morton_Assist::bit_otherlevel + i * C_DIMS + 1))
		{
			y_temp += C_dx * two_power_n(i);
		}
		if (key.test(Morton_Assist::bit_otherlevel + i * C_DIMS + 2))
		{
			z_temp += C_dx * two_power_n(i);
		}
	}
	x = x_temp;
	y = y_temp;
	z = z_temp;
};
#endif