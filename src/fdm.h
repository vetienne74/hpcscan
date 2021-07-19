
//-------------------------------------------------------------------------------------------------------
// Definition of FD operators:
// - FD_D2_OX_N1 = 2nd derivative along N1 axis with stencil order X (for grids 1D, 2D and 3D)
// - FD_D2_OX_N2 = 2nd derivative along N2 axis with stencil order X (for grids 2D and 3D)
// - FD_D2_OX_N3 = 2nd derivative along N3 axis with stencil order X (for grids 3D)
//
// All operators apply on 1d array (grid->3d_grid)
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_FDM_H_
#define HPCSCAN_FDM_H_

#include <cmath>
#include <vector>

#include "constant.h"

namespace hpcscan {

//=============================================== 1st differential ==========================================================

// NOTE: 1st differentials are not yet used in this version, for future use cases

// -------------------------------------%%%%%%%%%%% D1 FD space O8 %%%%%%%%%%%%%%-------------------------------------
const Myfloat FD_D1_O8_A1   =  1225./1024. ;
const Myfloat FD_D1_O8_A2   = -245./3072. ;
const Myfloat FD_D1_O8_A3   =  49./5120. ;
const Myfloat FD_D1_O8_A4   = -5./7168. ;
const Myint   FD_D1_O8_NOP  = 12 ;

#define FD_D1_O8_N1(U, i1, i2, i3, inv_d1, inv_d2, inv_d3, n1, n2, n3) \
		((FD_D1_O8_A1         * (U[i1   + i2*n1 + i3*n2*n1] - U[i1-1 + i2*n1 + i3*n2*n1])  \
				+ FD_D1_O8_A2 * (U[i1+1 + i2*n1 + i3*n2*n1] - U[i1-2 + i2*n1 + i3*n2*n1])  \
				+ FD_D1_O8_A3 * (U[i1+2 + i2*n1 + i3*n2*n1] - U[i1-3 + i2*n1 + i3*n2*n1])  \
				+ FD_D1_O8_A4 * (U[i1+3 + i2*n1 + i3*n2*n1] - U[i1-4 + i2*n1 + i3*n2*n1])) \
				* inv_d1)

#define FD_D1_O8_N2(U, i1, i2, i3, inv_d1, inv_d2, inv_d3, n1, n2, n3) \
		((FD_D1_O8_A1         * (U[i1 + (i2)  *n1 + i3*n2*n1] - U[i1 + (i2-1)*n1 + i3*n2*n1])  \
				+ FD_D1_O8_A2 * (U[i1 + (i2+1)*n1 + i3*n2*n1] - U[i1 + (i2-2)*n1 + i3*n2*n1])  \
				+ FD_D1_O8_A3 * (U[i1 + (i2+2)*n1 + i3*n2*n1] - U[i1 + (i2-3)*n1 + i3*n2*n1])  \
				+ FD_D1_O8_A4 * (U[i1 + (i2+3)*n1 + i3*n2*n1] - U[i1 + (i2-4)*n1 + i3*n2*n1])) \
				* inv_d2)

#define FD_D1_O8_N3(U, i1, i2, i3, inv_d1, inv_d2, inv_d3, n1, n2, n3) \
		((FD_D1_O8_A1         * (U[i1 + i2*n1 + (i3)  *n2*n1] - U[i1 + i2*n1 + (i3-1)*n2*n1])  \
				+ FD_D1_O8_A2 * (U[i1 + i2*n1 + (i3+1)*n2*n1] - U[i1 + i2*n1 + (i3-2)*n2*n1])  \
				+ FD_D1_O8_A3 * (U[i1 + i2*n1 + (i3+2)*n2*n1] - U[i1 + i2*n1 + (i3-3)*n2*n1])  \
				+ FD_D1_O8_A4 * (U[i1 + i2*n1 + (i3+3)*n2*n1] - U[i1 + i2*n1 + (i3-4)*n2*n1])) \
				* inv_d3)

//-=============================================== 2nd differential ==========================================================

//-------------------------------------%%%%%%%%%%% D2 FD space O2 %%%%%%%%%%%%%%------------------------------------
const Myfloat FD_D2_O2_A0   = -2.0 ;
const Myfloat FD_D2_O2_A1   =  1.0 ;
const Myint   FD_D2_O2_NOP  = 4 ;

#define FD_D2_O2_N1(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		((FD_D2_O2_A0 * U[i1 + i2*n1 + i3*n2*n1] \
				+ U[i1+1 + i2*n1 + i3*n2*n1] + U[i1-1 + i2*n1 + i3*n2*n1]) * inv2_d1)

#define FD_D2_O2_N2(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		((FD_D2_O2_A0 * U[i1 + i2*n1 + i3*n2*n1] \
				+ U[i1 + (i2+1)*n1 + i3*n2*n1] + U[i1 + (i2-1)*n1 + i3*n2*n1]) * inv2_d2)

#define FD_D2_O2_N3(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		((FD_D2_O2_A0 * U[i1 + i2*n1 + i3*n2*n1] \
				+ U[i1 + i2*n1 + (i3+1)*n2*n1] + U[i1 + i2*n1 + (i3-1)*n2*n1]) * inv2_d3)

// -------------------------------------%%%%%%%%%%%% D2 FD space O4 %%%%%%%%%%%%%-------------------------------------
const Myfloat FD_D2_O4_A0   = -5./2. ;
const Myfloat FD_D2_O4_A1   =  4./3. ;
const Myfloat FD_D2_O4_A2   = -1./12. ;
const Myint   FD_D2_O4_NOP  = 8 ;

#define FD_D2_O4_N1(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O4_A0 *  U[i1   + i2*n1 + i3*n2*n1] \
				+ FD_D2_O4_A1 * (U[i1+1 + i2*n1 + i3*n2*n1] + U[i1-1 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O4_A2 * (U[i1+2 + i2*n1 + i3*n2*n1] + U[i1-2 + i2*n1 + i3*n2*n1])) \
				* inv2_d1)

#define FD_D2_O4_N2(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O4_A0 *  U[i1 + i2    *n1 + i3*n2*n1] \
				+ FD_D2_O4_A1 * (U[i1 + (i2+1)*n1 + i3*n2*n1] + U[i1 + (i2-1)*n1 + i3*n2*n1])  \
				+ FD_D2_O4_A2 * (U[i1 + (i2+2)*n1 + i3*n2*n1] + U[i1 + (i2-2)*n1 + i3*n2*n1])) \
				* inv2_d2)

#define FD_D2_O4_N3(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O4_A0 *  U[i1 + i2*n1 + i3    *n2*n1] \
				+ FD_D2_O4_A1 * (U[i1 + i2*n1 + (i3+1)*n2*n1] + U[i1 + i2*n1 + (i3-1)*n2*n1])  \
				+ FD_D2_O4_A2 * (U[i1 + i2*n1 + (i3+2)*n2*n1] + U[i1 + i2*n1 + (i3-2)*n2*n1])) \
				* inv2_d3)

// -------------------------------------%%%%%%%%%%%% D2 FD space O6 %%%%%%%%%%%%%-------------------------------------
const Myfloat FD_D2_O6_A0   = -49./18. ;
const Myfloat FD_D2_O6_A1   =  3./2. ;
const Myfloat FD_D2_O6_A2   = -3./20. ;
const Myfloat FD_D2_O6_A3   =  1./90. ;
const Myint   FD_D2_O6_NOP  = 11 ;

#define FD_D2_O6_N1(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O6_A0 *  U[i1   + i2*n1 + i3*n2*n1] \
				+ FD_D2_O6_A1 * (U[i1+1 + i2*n1 + i3*n2*n1] + U[i1-1 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O6_A2 * (U[i1+2 + i2*n1 + i3*n2*n1] + U[i1-2 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O6_A3 * (U[i1+3 + i2*n1 + i3*n2*n1] + U[i1-3 + i2*n1 + i3*n2*n1])) \
				* inv2_d1)

#define FD_D2_O6_N2(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O6_A0 *  U[i1 + i2    *n1 + i3*n2*n1] \
				+ FD_D2_O6_A1 * (U[i1 + (i2+1)*n1 + i3*n2*n1] + U[i1 + (i2-1)*n1 + i3*n2*n1])  \
				+ FD_D2_O6_A2 * (U[i1 + (i2+2)*n1 + i3*n2*n1] + U[i1 + (i2-2)*n1 + i3*n2*n1])  \
				+ FD_D2_O6_A3 * (U[i1 + (i2+3)*n1 + i3*n2*n1] + U[i1 + (i2-3)*n1 + i3*n2*n1])) \
				* inv2_d2)

#define FD_D2_O6_N3(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O6_A0 *  U[i1 + i2*n1 + i3    *n2*n1] \
				+ FD_D2_O6_A1 * (U[i1 + i2*n1 + (i3+1)*n2*n1] + U[i1 + i2*n1 + (i3-1)*n2*n1])  \
				+ FD_D2_O6_A2 * (U[i1 + i2*n1 + (i3+2)*n2*n1] + U[i1 + i2*n1 + (i3-2)*n2*n1])  \
				+ FD_D2_O6_A3 * (U[i1 + i2*n1 + (i3+3)*n2*n1] + U[i1 + i2*n1 + (i3-3)*n2*n1])) \
				* inv2_d3)

// -------------------------------------%%%%%%%%%%%% D2 FD space O8 %%%%%%%%%%%%-------------------------------------
const Myfloat FD_D2_O8_A0   = -205./72. ;
const Myfloat FD_D2_O8_A1   =  8./5. ;
const Myfloat FD_D2_O8_A2   = -1/5. ;
const Myfloat FD_D2_O8_A3   =  8./315. ;
const Myfloat FD_D2_O8_A4   = -1/560. ;
const Myint   FD_D2_O8_NOP  = 14 ;

#define FD_D2_O8_N1(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O8_A0 *  U[i1   + i2*n1 + i3*n2*n1] \
				+ FD_D2_O8_A1 * (U[i1+1 + i2*n1 + i3*n2*n1] + U[i1-1 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O8_A2 * (U[i1+2 + i2*n1 + i3*n2*n1] + U[i1-2 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O8_A3 * (U[i1+3 + i2*n1 + i3*n2*n1] + U[i1-3 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O8_A4 * (U[i1+4 + i2*n1 + i3*n2*n1] + U[i1-4 + i2*n1 + i3*n2*n1])) \
				* inv2_d1)

#define FD_D2_O8_N2(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O8_A0 *  U[i1 + i2    *n1 + i3*n2*n1] \
				+ FD_D2_O8_A1 * (U[i1 + (i2+1)*n1 + i3*n2*n1] + U[i1 + (i2-1)*n1 + i3*n2*n1])  \
				+ FD_D2_O8_A2 * (U[i1 + (i2+2)*n1 + i3*n2*n1] + U[i1 + (i2-2)*n1 + i3*n2*n1])  \
				+ FD_D2_O8_A3 * (U[i1 + (i2+3)*n1 + i3*n2*n1] + U[i1 + (i2-3)*n1 + i3*n2*n1])  \
				+ FD_D2_O8_A4 * (U[i1 + (i2+4)*n1 + i3*n2*n1] + U[i1 + (i2-4)*n1 + i3*n2*n1])) \
				* inv2_d2)

#define FD_D2_O8_N3(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O8_A0 *  U[i1 + i2*n1 + i3    *n2*n1] \
				+ FD_D2_O8_A1 * (U[i1 + i2*n1 + (i3+1)*n2*n1] + U[i1 + i2*n1 + (i3-1)*n2*n1])  \
				+ FD_D2_O8_A2 * (U[i1 + i2*n1 + (i3+2)*n2*n1] + U[i1 + i2*n1 + (i3-2)*n2*n1])  \
				+ FD_D2_O8_A3 * (U[i1 + i2*n1 + (i3+3)*n2*n1] + U[i1 + i2*n1 + (i3-3)*n2*n1])  \
				+ FD_D2_O8_A4 * (U[i1 + i2*n1 + (i3+4)*n2*n1] + U[i1 + i2*n1 + (i3-4)*n2*n1])) \
				* inv2_d3)

// -------------------------------------%%%%%%%%%%%% D2 FD space 10 %%%%%%%%%%%%-------------------------------------
const Myfloat FD_D2_O10_A0   = -5269./1800. ;
const Myfloat FD_D2_O10_A1   =  5./3. ;
const Myfloat FD_D2_O10_A2   = -5./21. ;
const Myfloat FD_D2_O10_A3   =  5./126. ;
const Myfloat FD_D2_O10_A4   = -5./1008. ;
const Myfloat FD_D2_O10_A5   =  1./3150. ;
const Myint   FD_D2_O10_NOP  = 17 ;

#define FD_D2_O10_N1(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O10_A0 *  U[i1   + i2*n1 + i3*n2*n1] \
				+ FD_D2_O10_A1 * (U[i1+1 + i2*n1 + i3*n2*n1] + U[i1-1 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O10_A2 * (U[i1+2 + i2*n1 + i3*n2*n1] + U[i1-2 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O10_A3 * (U[i1+3 + i2*n1 + i3*n2*n1] + U[i1-3 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O10_A4 * (U[i1+4 + i2*n1 + i3*n2*n1] + U[i1-4 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O10_A5 * (U[i1+5 + i2*n1 + i3*n2*n1] + U[i1-5 + i2*n1 + i3*n2*n1])) \
				* inv2_d1)

#define FD_D2_O10_N2(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O10_A0 *  U[i1 + i2    *n1 + i3*n2*n1] \
				+ FD_D2_O10_A1 * (U[i1 + (i2+1)*n1 + i3*n2*n1] + U[i1 + (i2-1)*n1 + i3*n2*n1])  \
				+ FD_D2_O10_A2 * (U[i1 + (i2+2)*n1 + i3*n2*n1] + U[i1 + (i2-2)*n1 + i3*n2*n1])  \
				+ FD_D2_O10_A3 * (U[i1 + (i2+3)*n1 + i3*n2*n1] + U[i1 + (i2-3)*n1 + i3*n2*n1])  \
				+ FD_D2_O10_A4 * (U[i1 + (i2+4)*n1 + i3*n2*n1] + U[i1 + (i2-4)*n1 + i3*n2*n1])  \
				+ FD_D2_O10_A5 * (U[i1 + (i2+5)*n1 + i3*n2*n1] + U[i1 + (i2-5)*n1 + i3*n2*n1])) \
				* inv2_d2)

#define FD_D2_O10_N3(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O10_A0 *  U[i1 + i2*n1 + i3    *n2*n1] \
				+ FD_D2_O10_A1 * (U[i1 + i2*n1 + (i3+1)*n2*n1] + U[i1 + i2*n1 + (i3-1)*n2*n1])  \
				+ FD_D2_O10_A2 * (U[i1 + i2*n1 + (i3+2)*n2*n1] + U[i1 + i2*n1 + (i3-2)*n2*n1])  \
				+ FD_D2_O10_A3 * (U[i1 + i2*n1 + (i3+3)*n2*n1] + U[i1 + i2*n1 + (i3-3)*n2*n1])  \
				+ FD_D2_O10_A4 * (U[i1 + i2*n1 + (i3+4)*n2*n1] + U[i1 + i2*n1 + (i3-4)*n2*n1])  \
				+ FD_D2_O10_A5 * (U[i1 + i2*n1 + (i3+5)*n2*n1] + U[i1 + i2*n1 + (i3-5)*n2*n1])) \
				* inv2_d3)

// -------------------------------------%%%%%%%%%%%% D2 FD space 12 %%%%%%%%%%%%-------------------------------------
const Myfloat FD_D2_O12_A0   = -2598./871. ;
const Myfloat FD_D2_O12_A1   =  12./7.;
const Myfloat FD_D2_O12_A2   = -15./56. ;
const Myfloat FD_D2_O12_A3   =  10./189.;
const Myfloat FD_D2_O12_A4   = -1./112. ;
const Myfloat FD_D2_O12_A5   =  2./1925.;
const Myfloat FD_D2_O12_A6   = -1./16632. ;
const Myint   FD_D2_O12_NOP  = 20 ;

#define FD_D2_O12_N1(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O12_A0 *  U[i1   + i2*n1 + i3*n2*n1] \
				+ FD_D2_O12_A1 * (U[i1+1 + i2*n1 + i3*n2*n1] + U[i1-1 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O12_A2 * (U[i1+2 + i2*n1 + i3*n2*n1] + U[i1-2 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O12_A3 * (U[i1+3 + i2*n1 + i3*n2*n1] + U[i1-3 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O12_A4 * (U[i1+4 + i2*n1 + i3*n2*n1] + U[i1-4 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O12_A5 * (U[i1+5 + i2*n1 + i3*n2*n1] + U[i1-5 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O12_A6 * (U[i1+6 + i2*n1 + i3*n2*n1] + U[i1-6 + i2*n1 + i3*n2*n1])) \
				* inv2_d1)

#define FD_D2_O12_N2(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O12_A0 *  U[i1 + i2    *n1 + i3*n2*n1] \
				+ FD_D2_O12_A1 * (U[i1 + (i2+1)*n1 + i3*n2*n1] + U[i1 + (i2-1)*n1 + i3*n2*n1])  \
				+ FD_D2_O12_A2 * (U[i1 + (i2+2)*n1 + i3*n2*n1] + U[i1 + (i2-2)*n1 + i3*n2*n1])  \
				+ FD_D2_O12_A3 * (U[i1 + (i2+3)*n1 + i3*n2*n1] + U[i1 + (i2-3)*n1 + i3*n2*n1])  \
				+ FD_D2_O12_A4 * (U[i1 + (i2+4)*n1 + i3*n2*n1] + U[i1 + (i2-4)*n1 + i3*n2*n1])  \
				+ FD_D2_O12_A5 * (U[i1 + (i2+5)*n1 + i3*n2*n1] + U[i1 + (i2-5)*n1 + i3*n2*n1])  \
				+ FD_D2_O12_A6 * (U[i1 + (i2+6)*n1 + i3*n2*n1] + U[i1 + (i2-6)*n1 + i3*n2*n1])) \
				* inv2_d2)

#define FD_D2_O12_N3(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O12_A0 *  U[i1 + i2*n1 + i3    *n2*n1] \
				+ FD_D2_O12_A1 * (U[i1 + i2*n1 + (i3+1)*n2*n1] + U[i1 + i2*n1 + (i3-1)*n2*n1])  \
				+ FD_D2_O12_A2 * (U[i1 + i2*n1 + (i3+2)*n2*n1] + U[i1 + i2*n1 + (i3-2)*n2*n1])  \
				+ FD_D2_O12_A3 * (U[i1 + i2*n1 + (i3+3)*n2*n1] + U[i1 + i2*n1 + (i3-3)*n2*n1])  \
				+ FD_D2_O12_A4 * (U[i1 + i2*n1 + (i3+4)*n2*n1] + U[i1 + i2*n1 + (i3-4)*n2*n1])  \
				+ FD_D2_O12_A5 * (U[i1 + i2*n1 + (i3+5)*n2*n1] + U[i1 + i2*n1 + (i3-5)*n2*n1])  \
				+ FD_D2_O12_A6 * (U[i1 + i2*n1 + (i3+6)*n2*n1] + U[i1 + i2*n1 + (i3-6)*n2*n1])) \
				* inv2_d3)

// -------------------------------------%%%%%%%%%%%% D2 FD space 14 %%%%%%%%%%%%-------------------------------------
const Myfloat FD_D2_O14_A0   = -266681./88200. ;
const Myfloat FD_D2_O14_A1   =  7./4.;
const Myfloat FD_D2_O14_A2   = -7./24. ;
const Myfloat FD_D2_O14_A3   =  7./108.;
const Myfloat FD_D2_O14_A4   = -7./528. ;
const Myfloat FD_D2_O14_A5   =  7./3300.;
const Myfloat FD_D2_O14_A6   = -7./30888. ;
const Myfloat FD_D2_O14_A7   =  1./84084. ;
const Myint   FD_D2_O14_NOP  = 23 ;

#define FD_D2_O14_N1(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O14_A0 *  U[i1   + i2*n1 + i3*n2*n1] \
				+ FD_D2_O14_A1 * (U[i1+1 + i2*n1 + i3*n2*n1] + U[i1-1 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A2 * (U[i1+2 + i2*n1 + i3*n2*n1] + U[i1-2 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A3 * (U[i1+3 + i2*n1 + i3*n2*n1] + U[i1-3 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A4 * (U[i1+4 + i2*n1 + i3*n2*n1] + U[i1-4 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A5 * (U[i1+5 + i2*n1 + i3*n2*n1] + U[i1-5 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A6 * (U[i1+6 + i2*n1 + i3*n2*n1] + U[i1-6 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A7 * (U[i1+7 + i2*n1 + i3*n2*n1] + U[i1-7 + i2*n1 + i3*n2*n1])) \
				* inv2_d1)

#define FD_D2_O14_N2(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O14_A0 *  U[i1 + i2    *n1 + i3*n2*n1] \
				+ FD_D2_O14_A1 * (U[i1 + (i2+1)*n1 + i3*n2*n1] + U[i1 + (i2-1)*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A2 * (U[i1 + (i2+2)*n1 + i3*n2*n1] + U[i1 + (i2-2)*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A3 * (U[i1 + (i2+3)*n1 + i3*n2*n1] + U[i1 + (i2-3)*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A4 * (U[i1 + (i2+4)*n1 + i3*n2*n1] + U[i1 + (i2-4)*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A5 * (U[i1 + (i2+5)*n1 + i3*n2*n1] + U[i1 + (i2-5)*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A6 * (U[i1 + (i2+6)*n1 + i3*n2*n1] + U[i1 + (i2-6)*n1 + i3*n2*n1])  \
				+ FD_D2_O14_A7 * (U[i1 + (i2+7)*n1 + i3*n2*n1] + U[i1 + (i2-7)*n1 + i3*n2*n1])) \
				* inv2_d2)

#define FD_D2_O14_N3(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O14_A0 *  U[i1 + i2*n1 + i3    *n2*n1] \
				+ FD_D2_O14_A1 * (U[i1 + i2*n1 + (i3+1)*n2*n1] + U[i1 + i2*n1 + (i3-1)*n2*n1])  \
				+ FD_D2_O14_A2 * (U[i1 + i2*n1 + (i3+2)*n2*n1] + U[i1 + i2*n1 + (i3-2)*n2*n1])  \
				+ FD_D2_O14_A3 * (U[i1 + i2*n1 + (i3+3)*n2*n1] + U[i1 + i2*n1 + (i3-3)*n2*n1])  \
				+ FD_D2_O14_A4 * (U[i1 + i2*n1 + (i3+4)*n2*n1] + U[i1 + i2*n1 + (i3-4)*n2*n1])  \
				+ FD_D2_O14_A5 * (U[i1 + i2*n1 + (i3+5)*n2*n1] + U[i1 + i2*n1 + (i3-5)*n2*n1])  \
				+ FD_D2_O14_A6 * (U[i1 + i2*n1 + (i3+6)*n2*n1] + U[i1 + i2*n1 + (i3-6)*n2*n1])  \
				+ FD_D2_O14_A7 * (U[i1 + i2*n1 + (i3+7)*n2*n1] + U[i1 + i2*n1 + (i3-7)*n2*n1])) \
				* inv2_d3)

// -------------------------------------%%%%%%%%%%%% D2 FD space 16 %%%%%%%%%%%%-------------------------------------
const Myfloat FD_D2_O16_A0   = -1671./547. ;
const Myfloat FD_D2_O16_A1   =  16./9. ;
const Myfloat FD_D2_O16_A2   = -14./45. ;
const Myfloat FD_D2_O16_A3   =  112./1485.;
const Myfloat FD_D2_O16_A4   = -7./396. ;
const Myfloat FD_D2_O16_A5   =  65./18673. ;
const Myfloat FD_D2_O16_A6   = -2./3861.;
const Myfloat FD_D2_O16_A7   =  5./98536.;
const Myfloat FD_D2_O16_A8   = -1./411840. ;
const Myint   FD_D2_O16_NOP  = 26 ;

#define FD_D2_O16_N1(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O16_A0 *  U[i1   + i2*n1 + i3*n2*n1] \
				+ FD_D2_O16_A1 * (U[i1+1 + i2*n1 + i3*n2*n1] + U[i1-1 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A2 * (U[i1+2 + i2*n1 + i3*n2*n1] + U[i1-2 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A3 * (U[i1+3 + i2*n1 + i3*n2*n1] + U[i1-3 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A4 * (U[i1+4 + i2*n1 + i3*n2*n1] + U[i1-4 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A5 * (U[i1+5 + i2*n1 + i3*n2*n1] + U[i1-5 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A6 * (U[i1+6 + i2*n1 + i3*n2*n1] + U[i1-6 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A7 * (U[i1+7 + i2*n1 + i3*n2*n1] + U[i1-7 + i2*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A8 * (U[i1+8 + i2*n1 + i3*n2*n1] + U[i1-8 + i2*n1 + i3*n2*n1])) \
				* inv2_d1)

#define FD_D2_O16_N2(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O16_A0 *  U[i1 + i2    *n1 + i3*n2*n1] \
				+ FD_D2_O16_A1 * (U[i1 + (i2+1)*n1 + i3*n2*n1] + U[i1 + (i2-1)*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A2 * (U[i1 + (i2+2)*n1 + i3*n2*n1] + U[i1 + (i2-2)*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A3 * (U[i1 + (i2+3)*n1 + i3*n2*n1] + U[i1 + (i2-3)*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A4 * (U[i1 + (i2+4)*n1 + i3*n2*n1] + U[i1 + (i2-4)*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A5 * (U[i1 + (i2+5)*n1 + i3*n2*n1] + U[i1 + (i2-5)*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A6 * (U[i1 + (i2+6)*n1 + i3*n2*n1] + U[i1 + (i2-6)*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A7 * (U[i1 + (i2+7)*n1 + i3*n2*n1] + U[i1 + (i2-7)*n1 + i3*n2*n1])  \
				+ FD_D2_O16_A8 * (U[i1 + (i2+8)*n1 + i3*n2*n1] + U[i1 + (i2-8)*n1 + i3*n2*n1])) \
				* inv2_d2)

#define FD_D2_O16_N3(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
		(        (FD_D2_O16_A0 *  U[i1 + i2*n1 + i3    *n2*n1] \
				+ FD_D2_O16_A1 * (U[i1 + i2*n1 + (i3+1)*n2*n1] + U[i1 + i2*n1 + (i3-1)*n2*n1])  \
				+ FD_D2_O16_A2 * (U[i1 + i2*n1 + (i3+2)*n2*n1] + U[i1 + i2*n1 + (i3-2)*n2*n1])  \
				+ FD_D2_O16_A3 * (U[i1 + i2*n1 + (i3+3)*n2*n1] + U[i1 + i2*n1 + (i3-3)*n2*n1])  \
				+ FD_D2_O16_A4 * (U[i1 + i2*n1 + (i3+4)*n2*n1] + U[i1 + i2*n1 + (i3-4)*n2*n1])  \
				+ FD_D2_O16_A5 * (U[i1 + i2*n1 + (i3+5)*n2*n1] + U[i1 + i2*n1 + (i3-5)*n2*n1])  \
				+ FD_D2_O16_A6 * (U[i1 + i2*n1 + (i3+6)*n2*n1] + U[i1 + i2*n1 + (i3-6)*n2*n1])  \
				+ FD_D2_O16_A7 * (U[i1 + i2*n1 + (i3+7)*n2*n1] + U[i1 + i2*n1 + (i3-7)*n2*n1])  \
				+ FD_D2_O16_A8 * (U[i1 + i2*n1 + (i3+8)*n2*n1] + U[i1 + i2*n1 + (i3-8)*n2*n1])) \
				* inv2_d3)

//-------------------------------------------------------------------------------------------------------
// return vector of FD coefficients
//
static vector<Myfloat> getFD_D2coefVector(Myint fdOrder)
{
	vector<Myfloat> FD_coef ;

	if (fdOrder == 2)
	{
		FD_coef.push_back(FD_D2_O2_A0) ;
		FD_coef.push_back(FD_D2_O2_A1) ;
	}
	else if (fdOrder == 4)
	{
		FD_coef.push_back(FD_D2_O4_A0) ;
		FD_coef.push_back(FD_D2_O4_A1) ;
		FD_coef.push_back(FD_D2_O4_A2) ;
	}
	else if (fdOrder == 6)
	{
		FD_coef.push_back(FD_D2_O6_A0) ;
		FD_coef.push_back(FD_D2_O6_A1) ;
		FD_coef.push_back(FD_D2_O6_A2) ;
		FD_coef.push_back(FD_D2_O6_A3) ;
	}
	else if (fdOrder == 8)
	{
		FD_coef.push_back(FD_D2_O8_A0) ;
		FD_coef.push_back(FD_D2_O8_A1) ;
		FD_coef.push_back(FD_D2_O8_A2) ;
		FD_coef.push_back(FD_D2_O8_A3) ;
		FD_coef.push_back(FD_D2_O8_A4) ;
	}
	else if (fdOrder == 10)
	{
		FD_coef.push_back(FD_D2_O10_A0) ;
		FD_coef.push_back(FD_D2_O10_A1) ;
		FD_coef.push_back(FD_D2_O10_A2) ;
		FD_coef.push_back(FD_D2_O10_A3) ;
		FD_coef.push_back(FD_D2_O10_A4) ;
		FD_coef.push_back(FD_D2_O10_A5) ;
	}
	else if (fdOrder == 12)
	{
		FD_coef.push_back(FD_D2_O12_A0) ;
		FD_coef.push_back(FD_D2_O12_A1) ;
		FD_coef.push_back(FD_D2_O12_A2) ;
		FD_coef.push_back(FD_D2_O12_A3) ;
		FD_coef.push_back(FD_D2_O12_A4) ;
		FD_coef.push_back(FD_D2_O12_A5) ;
		FD_coef.push_back(FD_D2_O12_A6) ;
	}
	else if (fdOrder == 14)
	{
		FD_coef.push_back(FD_D2_O14_A0) ;
		FD_coef.push_back(FD_D2_O14_A1) ;
		FD_coef.push_back(FD_D2_O14_A2) ;
		FD_coef.push_back(FD_D2_O14_A3) ;
		FD_coef.push_back(FD_D2_O14_A4) ;
		FD_coef.push_back(FD_D2_O14_A5) ;
		FD_coef.push_back(FD_D2_O14_A6) ;
		FD_coef.push_back(FD_D2_O14_A7) ;
	}
	else if (fdOrder == 16)
	{
		FD_coef.push_back(FD_D2_O16_A0) ;
		FD_coef.push_back(FD_D2_O16_A1) ;
		FD_coef.push_back(FD_D2_O16_A2) ;
		FD_coef.push_back(FD_D2_O16_A3) ;
		FD_coef.push_back(FD_D2_O16_A4) ;
		FD_coef.push_back(FD_D2_O16_A5) ;
		FD_coef.push_back(FD_D2_O16_A6) ;
		FD_coef.push_back(FD_D2_O16_A7) ;
		FD_coef.push_back(FD_D2_O16_A8) ;
	}

	return FD_coef ;
}

//-------------------------------------------------------------------------------------------------------
// compute sum of abs(FD coefficients)
// this is used to compute the CFL for the propagator
//
static Myfloat64 getSumAbsFD_D2Coef(Myint fdOrder)
{
	Myfloat64 sumCoef = UNSPECIFIED ;

	auto FD_coef = getFD_D2coefVector(fdOrder) ;

	for (Myint ii = 0; ii < FD_coef.size(); ii++)
	{
		if (ii == 0) sumCoef = fabs(FD_coef[ii]) ;
		else sumCoef += TWO * fabs(FD_coef[ii]) ;
	}

	return(sumCoef) ;
}

//-------------------------------------------------------------------------------------------------------
// compute sum FD coefficients
// this is used for validation purpose only (in testCase_Util)
// this sum is expected to be equal to 0
//
static Myfloat64 getSumFD_D2Coef(Myint fdOrder)
{
	Myfloat64 sumCoef = UNSPECIFIED ;

	auto FD_coef = getFD_D2coefVector(fdOrder) ;

	for (Myint ii = 0; ii < FD_coef.size(); ii++)
	{
		if (ii == 0) sumCoef = FD_coef[ii] ;
		else sumCoef += TWO * FD_coef[ii] ;
	}

	return(sumCoef) ;
}

//-------------------------------------------------------------------------------------------------------
// get halo width associated with a given FD order
//
static Myint getFD_D2haloWidth(Myint fdOrder)
{
	auto FD_coef = getFD_D2coefVector(fdOrder) ;
	Myint haloWidth = FD_coef.size() - 1 ;
	return haloWidth  ;
}

//-------------------------------------------------------------------------------------------------------
// get number of mathematical operations for a given FD order
// for derivative computations along one axis
//
static Myint getFD_D2nMathOp(Myint fdOrder)
{
	Myint nMathOp = UNSPECIFIED ;

	if (fdOrder == 2)
	{
		nMathOp = FD_D2_O2_NOP ;
	}
	else if (fdOrder == 4)
	{
		nMathOp = FD_D2_O4_NOP ;
	}
	else if (fdOrder == 6)
	{
		nMathOp = FD_D2_O6_NOP ;
	}
	else if (fdOrder == 8)
	{
		nMathOp = FD_D2_O8_NOP ;
	}
	else if (fdOrder == 10)
	{
		nMathOp = FD_D2_O10_NOP ;
	}
	else if (fdOrder == 12)
	{
		nMathOp = FD_D2_O12_NOP ;
	}
	else if (fdOrder == 14)
	{
		nMathOp = FD_D2_O14_NOP ;
	}
	else if (fdOrder == 16)
	{
		nMathOp = FD_D2_O16_NOP ;
	}

	return(nMathOp) ;
}


} // namespace hpcscan

#endif
