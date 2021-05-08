#ifndef HPCSCAN_TYPE_DEF_H_
#define HPCSCAN_TYPE_DEF_H_

using namespace std;

namespace hpcscan {

//--------------------------------------------------
// 
//          T Y P E   D E F I N I T I O N S
//
//--------------------------------------------------

// Basic types
#ifdef _DOUBLE_PRECISION_
typedef double Myfloat ;
#define MPI_MYFLOAT MPI_DOUBLE
#else
typedef float Myfloat ;
#define MPI_MYFLOAT MPI_FLOAT
#endif
#define MPI_MYFLOAT64 MPI_DOUBLE

typedef int        Myint ;
typedef int        Myint32 ;
typedef long int   Myint64 ;
typedef float      Myfloat32 ;
typedef double     Myfloat64 ;

// Return codes used by all methods
enum Rtn_code {RTN_CODE_OK=0, RTN_CODE_KO=-1, RTN_CODE_EXIT=-2} ;

// Debug level
enum Debug_level {NO_DEBUG=-1, LIGHT_DEBUG, MID_DEBUG, FULL_DEBUG} ;

// For display in output report
enum Display_type {MASTER=0, ALL=1} ;

// Spatial axis
enum Axis_type {AXIS1, AXIS2, AXIS3} ;

// Spatial dimension
enum Dim_type {DIM1=1, DIM2=2, DIM3=3} ;

// FD type
enum FD_type {FD_D1, FD_D2} ;

// Function type
enum Func_type {FUNC_NONE, FUNC_SINE, FUNC_COSINE, FUNC_CONST, FUNC_LINEAR, FUNC_DIRAC, FUNC_RANDOM} ;

// Point type
enum Point_type{INNER_POINTS, ALL_POINTS, MIDDLE_POINT,
	I1HALO1, I1HALO2, I2HALO1, I2HALO2, I3HALO1, I3HALO2,
	I1INNERHALO1, I1INNERHALO2, I2INNERHALO1, I2INNERHALO2, I3INNERHALO1, I3INNERHALO2} ;

// Grid type
enum Grid_type{GRID_LOCAL, GRID_GLOBAL} ;

// Special values
enum {UNSPECIFIED = -1} ;

// MPI comm. mode
enum MPI_comm_mode_type {MPI_COMM_MODE_SENDRECV} ;

// Boundary condition
enum BoundCond_type {NO_BOUND_COND, BOUND_COND_ANTI_MIRROR} ;

// Propagator initialization type
enum PropaInit_type {EIGEN_MODE} ;

} // namespace hpcscan

#endif
