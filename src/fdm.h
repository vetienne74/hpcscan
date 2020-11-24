
// FD operators
// All operators apply on 1d array (grid->3d_grid)

#ifndef HPCSCAN_FDM_H_
#define HPCSCAN_FDM_H_

#include <cmath>

namespace hpcscan {

//=============================================== 1st differential ==========================================================

// -------------------------------------%%%%%%%%%%% D1 FD space O8 %%%%%%%%%%%%%%-------------------------------------
const Myint   FD_D1_O8_HSTL = 4 ;
const Myfloat FD_D1_O8_A1   = 1225./1024. ;
const Myfloat FD_D1_O8_A2   = -245./3072. ;
const Myfloat FD_D1_O8_A3   = 49./5120. ;
const Myfloat FD_D1_O8_A4   = -5./7168. ;
const Myfloat FD_D1_O8_CFL  = 0.777 ;
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
const Myint   FD_D2_O2_HSTL = 1 ;
const Myfloat FD_D2_O2_A0   = -2.0 ;
const Myfloat FD_D2_O2_A1   = 1.0 ;
const Myfloat FD_D2_O2_CFL  = 1.0 / sqrt(3.0) ;
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
const Myint   FD_D2_O4_HSTL = 2 ; //half stencil
const Myfloat FD_D2_O4_A0   =  -5./2. ;
const Myfloat FD_D2_O4_A1   =  4./3. ;
const Myfloat FD_D2_O4_A2   =  -1./12. ;
const Myfloat FD_D2_O4_CFL  = 0.857 ;
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


// -------------------------------------%%%%%%%%%%%% D2 FD space O8 %%%%%%%%%%%%-------------------------------------
const Myint   FD_D2_O8_HSTL = 4 ;
const Myfloat FD_D2_O8_A0   = -205./72. ;
const Myfloat FD_D2_O8_A1   = 8./5. ;
const Myfloat FD_D2_O8_A2   = -1/5. ;
const Myfloat FD_D2_O8_A3   = 8./315. ;
const Myfloat FD_D2_O8_A4   = -1/560. ;
const Myfloat FD_D2_O8_CFL  = FD_D1_O8_CFL ;
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

// -------------------------------------%%%%%%%%%%%% D2 FD space 12 %%%%%%%%%%%%-------------------------------------
const Myint   FD_D2_O12_HSTL = 6 ;
const Myfloat FD_D2_O12_A0   = -2598./871. ;
const Myfloat FD_D2_O12_A1   =  12./7.;
const Myfloat FD_D2_O12_A2   = -15./56. ;
const Myfloat FD_D2_O12_A3   =  10./189.;
const Myfloat FD_D2_O12_A4   = -1./112. ;
const Myfloat FD_D2_O12_A5   =  2./1925.;
const Myfloat FD_D2_O12_A6   = -1./16632. ;
const Myfloat FD_D2_O12_CFL  = 0.746 ;
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


// -------------------------------------%%%%%%%%%%%% D2 FD space 16 %%%%%%%%%%%%-------------------------------------
const Myint   FD_D2_O16_HSTL = 8 ;
const Myfloat FD_D2_O16_A0   = -1671./547. ;
const Myfloat FD_D2_O16_A1   = 16./9. ;
const Myfloat FD_D2_O16_A2   = -14./45. ;
const Myfloat FD_D2_O16_A3   =  112./1485.;
const Myfloat FD_D2_O16_A4   = -7./396. ;
const Myfloat FD_D2_O16_A5   = 65./18673. ;
const Myfloat FD_D2_O16_A6   =  -2./3861.;
const Myfloat FD_D2_O16_A7   =  5./98536.;
const Myfloat FD_D2_O16_A8   = -1./411840. ;
const Myfloat FD_D2_O16_CFL  = 0.729 ;
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

static Myfloat64 getSumAbsFD_D2Coef(Myint fdOrder)
{
	Myfloat64 sumCoef = 0.0 ;
	if (fdOrder == 2)
	{
		sumCoef = fabs(FD_D2_O2_A0) + 2*fabs(FD_D2_O2_A1) ;
	}
	else if (fdOrder == 4)
	{
		sumCoef = fabs(FD_D2_O4_A0) + 2*fabs(FD_D2_O4_A1) + 2*fabs(FD_D2_O4_A2) ;
	}
	else if (fdOrder == 8)
	{
		sumCoef = fabs(FD_D2_O8_A0) + 2*fabs(FD_D2_O8_A1) + 2*fabs(FD_D2_O8_A2)
		+ 2*fabs(FD_D2_O8_A3) + 2*fabs(FD_D2_O8_A4);
	}
	else if (fdOrder == 12)
	{
		sumCoef = fabs(FD_D2_O12_A0) + 2*fabs(FD_D2_O12_A1) + 2*fabs(FD_D2_O12_A2)
		+ 2*fabs(FD_D2_O12_A3) + 2*fabs(FD_D2_O12_A4);
	}
	else if (fdOrder == 16)
	{
		sumCoef = fabs(FD_D2_O16_A0) + 2*fabs(FD_D2_O16_A1) + 2*fabs(FD_D2_O16_A2)
		+ 2*fabs(FD_D2_O16_A3) + 2*fabs(FD_D2_O16_A4)
		+ 2*fabs(FD_D2_O16_A5) + 2*fabs(FD_D2_O16_A6)
		+ 2*fabs(FD_D2_O16_A7) + 2*fabs(FD_D2_O16_A8) ;
	}

	return(sumCoef) ;
}

} // namespace hpcscan

#endif
