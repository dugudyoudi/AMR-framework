/**
* @file
* @author Zhengliang Liu
* @brief Constants used in the solver.
* @note define constants include: 1) lattice paramters; 2) phycical parameters; 3) numerical parameters; and 4) other parameters.
*/
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include "General.h"
#define C_DEBUG 1     // when set as 1, (a) add additional output to check grid points added and removed
#define C_CHECK_MORTON_BOUNDARY 1  // when set as 1, will check if the node exceeds boudnaries of the backgound block when generate mesh. Must be set as 1 when solid boudnary
                            // near the computatinal domain, such as channel flow
#define C_SPEED_UP 1 // when set as 1, the project will speed up by negleting some processess. This may result in some problems. These processes need further improvement
// 1) The process to add all nodes in cells on the innermost boundary. Since all nodes in cells on the boundary close to the innermost boundary will be added,
//    generally neglect this process is OK if the boundary is smooth. In Grid_Manager::search_nodes_near_solid(...). It will influence nodes inside the solid boundary. Seems time saving is insignificant.

#define C_MAP_TYPE 1 // when set as 1, use unordered_map, else, use map
 
#define C_SEARCH_METHOD 2 // method to identify if nodes need to be removed, when set as 3, use icount_refine in a regio;
                                                                          // when set as 2, use two way search and use icout_refine on sreaching boundary
                                                                          // when set as 1, use two way search
                                                                          // 1 and 2 does not impletement 2D, and does not check numerical boundary when find nodes on the searching boundary

// Lattice parameters
#define C_DIMS 3                      ///< Number of dimensions
#define C_Q 9                          ///< Number of discrete velocities
const D_real C_cs = sqrt(3);           ///< Lattice sound speed
const D_real C_cssq = 3.0;             ///< Square lattice sound speed

// Grid informaion
#define C_FSI_INTERFACE 1               /// if C_SOLID_BOUNDARY == 1, using IBM, an addtional map will be used to store number of points near the node;
#define C_SOLID_BOUNDARY 2              ///< Choose method to implement solid boundary condition, 1 static boundary, 2 moving boundary;
                                        //  if C_SOLID_BOUNDARY == 2: update nodes at each finest time step
const D_real C_dx = 256.;// / 1000 / 10;//256 for airplane;            ///< Grid space of the background mesh
const D_real C_xb = 20000;// / 10000;// 20000 for airplane;                ///< Background boundary distance in x direction
const D_real C_yb = 20000;// / 10000;//./1000;                ///< Background boundary distance in y direction
#if (C_DIMS==3)
const D_real C_zb = 20000;// / 10000;//./1000;    ///< Background boundary distance in z direction
#endif

const int C_x0b_offset = 1;    ///< offset to avoid Morton_assit::find_x0 exceeding the boundary limit, the offset distance is (x0b_offset * C_dx)
const int C_y0b_offset = 1;
#if (C_DIMS == 3)
const int C_z0b_offset = 1;
#endif

 // grid refinement
#define C_BIT 64                               ///< Bits to store morton code 
const unsigned int C_max_level = 4;            ///< Maximum refinement level, the background mesh is level 0, others are 1, 2, 3, ..., C_maxlevel
const unsigned int C_overlap = 1;              ///< Number of overlapping grids (corresponding to coarser grid at ilevel - 1)
const unsigned int C_extend_inner = 8;         ///< Number of nodes extended from solid points in the finest block
const unsigned int C_extend_inner_x0 = 0;      ///< additional extension in -x direction (inner block)
const unsigned int C_extend_inner_x1 = 0;      ///< additional extension in +x direction (inner block)
const unsigned int C_extend_inner_y0 = 0;      ///< additional extension in -y direction (inner block)
const unsigned int C_extend_inner_y1 = 0;      ///< additional extension in +y direction (inner block)
#if (C_DIMS==3)
const unsigned int C_extend_inner_z0 = 0;      ///< additional extension in -z direction (inner block)
const unsigned int C_extend_inner_z1 = 0;      ///< additional extension in +z direction (inner block)
#endif
const unsigned int C_extend = 8;                ///< Number of nodes extended from the inner block for outer blocks (iblock < C_max_level)
const unsigned int C_extend_outer_x0 = 0;      ///< additional extension in -x direction (outer block)
const unsigned int C_extend_outer_x1 = 0;      ///< additional extension in +x direction (outer block)
const unsigned int C_extend_outer_y0 = 0;      ///< additional extension in -y direction (outer block)
const unsigned int C_extend_outer_y1 = 0;      ///< additional extension in +y direction (outer block)
#if (C_DIMS==3)
const unsigned int C_extend_outer_z0 = 0;      ///< additional extension in -z direction (outer block)
const unsigned int C_extend_outer_z1 = 0;      ///< additional extension in +z direction (outer block)
#endif
const unsigned int C_extend_ghost = 3;      ///< Number of nodes extended from solid points inside the solid

#if (C_FSI_INTERFACE == 1)
const unsigned int C_extend_IB = 2;      ///< Number of nodes extended from solid points for IBM
#endif

using D_morton = std::bitset<C_BIT>;         // Bitset to store morton code, due to the method used to find neighbours in Morton_Assist.h, the bit is limited to unsigned long long

// Physical parameters
 // reference parameters to calculate Reynolds number and for initialization
const D_real C_u0 = 1.0;		       ///< Reference and initial velocity
const D_real C_rho0 = 1.0;      ///< Reference and initial density	
const D_real C_l0 = 1.0;		       ///< Reference length
const D_real C_mu = 0.0001;		       ///< Reference viscosoity
 // initial parameters
const double C_v0 = 0.0;		       ///< Reference and initial velocity
#if (C_DIMS==3)
const D_real C_w0 = 0.0;		       ///< Reference and initial velocity	
#endif
 // other paramters
const D_real C_gravity = 0.;		   ///< Force density due to gravity

// Numerical parameters
const D_real C_eps = 1e-10;            ///< Small number for comparing floats to zero
const D_real C_pi = 4.*atan(1.);       ///< Pi
#endif