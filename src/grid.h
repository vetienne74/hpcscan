
#ifndef HPCSCAN_GRID_H_
#define HPCSCAN_GRID_H_

#include <string>

#include "mpi.h"

#include "type_def.h"

namespace hpcscan {

class Grid
{
public:

	// constructor
	Grid(Grid_type) ;

	// constructor
	Grid(Grid_type, Dim_type, Myint64, Myint64, Myint64) ;

	// destructor
	~Grid(void) ;

	// initialise grid index and MPI data structure
	virtual void initializeGrid(void) ;

	// print info
	virtual void info(void) ;

	// write grid on disk
	void write(string) ;

	// L1 error between this grid and another
	virtual Myfloat L1Err(Point_type pointType, const Grid& gridIn) const ;

	// collective (all MPI process) L1 error between this grid and another
	Myfloat allProcL1Err(Point_type, const Grid&) const ;

	// Max error between this grid and another (point wise)
	virtual Myfloat maxErr(Point_type, const Grid&) const ;

	// get grid indexes
	Rtn_code getGridIndex(Point_type, Myint64* i1Start, Myint64* i1End,
			Myint64* i2Start, Myint64* i2End, Myint64* i3Start, Myint64* i3End) const ;

	// get min and max in grid
	virtual Myfloat getMin(Point_type) ;
	virtual Myfloat getMax(Point_type) ;

	// get sum of abs grid points
	virtual Myfloat getSumAbs(Point_type) const ;

	// get sum of abs diff of grid points
	virtual Myfloat getSumAbsDiff(Point_type, const Grid&) const;

	// get number of neighbours
	Myint getNumberOfNeighbour() ;

	// get number of grid points according to point type
	Myint64 getNumberOfGridPoint(Grid_type gridType, Point_type pointType) ;

	// get number of grid points involved in exchange halos
	Myint64 getNumberOfGridPointCommHalo(Grid_type gridType) ;

	// get number of grid points involved in boundary condition
	Myint64 getNumberOfGridPointBoundaryCondition(BoundCond_type boundCondType) ;

	// get min coordinates
	Myfloat getMinCoord(Axis_type) ;

	// get max coordinates
	Myfloat getMaxCoord(Axis_type) ;

	// get min spatial sampling
	Myfloat getMinSpaceSampling(void) ;

	// get neighbour proc
	Myint getNeighbourProc(Point_type) ;

	// check grids have same size
	bool sameSize(const Grid&) const ;

	// fill grid with constant value
	virtual void fill(Point_type pointType, Myfloat val) ;

	// fill grid with predefined functions
	virtual void fill(Point_type pType, Func_type t1,  Func_type t2, Func_type t3,
			Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp) ;

	// define unit grid
	void defineUnitGrid() ;

	// exchange halos with MPI
	Rtn_code exchangeHalos(MPI_comm_mode_type) ;

	// exchange one halo with MPI
	Rtn_code exchangeHalo(MPI_comm_mode_type, Point_type pointType) ;

	// apply boundary condition
	virtual Rtn_code applyBoundaryCondition(BoundCond_type boundCondType) ;

	// compute FD_D2 along N1
	virtual Rtn_code FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N2
	virtual Rtn_code FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N3
	virtual Rtn_code FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_LAPLACIAN
	virtual Rtn_code FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// update pressure
	virtual Rtn_code updatePressure(Point_type pType, const Grid& prcGrid,
			const Grid& coefGrid, const Grid& laplaGrid) ;

	// compute pressure with FD
	virtual Rtn_code computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder) ;

	// get number of FLOP per point for operator FD_D2
	Myint getFlopPerPtFD_D2(Myint fdOrder) ;

	// get number of FLOP per point for operator LAPLACIAN
	Myint getFlopPerPtFD_LAPLACIAN(Myint fdOrder) ;

	// get number of points in the stencil for operator FD_D2
	Myint getPtPerStencilFD_D2(Myint fdOrder) ;

	// get number of points in the stencil for operator LAPLACIAN
	Myint getPtPerStencilFD_LAPLACIAN(Myint fdOrder) ;

	// Space dimension
	Dim_type dim ;

	// Grid type
	Grid_type gridType ;

	// 3D grid array
	Myfloat * grid_3d ;   // this will be allocated by Grid::initializeGrid on the CPU (we should remove it soon)
	Myfloat * d_grid_3d ; // this will be a pointer to device (GPU) memory


	// Grid size (local)
	Myint n1, n2, n3 ;

	// grid sampling
	Myfloat64 d1, d2, d3 ;

	// halo width (half of FD stencil)
	Myint haloWidth ;

	// grid mode
	string gridMode ;

protected:

	// array size (points)
	Myint64 npoint ;

	// inner points
	Myint64 n1Inner, n2Inner, n3Inner ;

	// grid indexes (local)
	Myint i1Halo1Start , i1Halo1End , i2Halo1Start , i2Halo1End , i3Halo1Start , i3Halo1End ;
	Myint i1Halo2Start , i1Halo2End , i2Halo2Start , i2Halo2End , i3Halo2Start , i3Halo2End ;
	Myint i1Layer1Start, i1Layer1End, i2Layer1Start, i2Layer1End, i3Layer1Start, i3Layer1End ;
	Myint i1Layer2Start, i1Layer2End, i2Layer2Start, i2Layer2End, i3Layer2Start, i3Layer2End ;
	Myint i1InnerStart , i1InnerEnd , i2InnerStart , i2InnerEnd , i3InnerStart , i3InnerEnd ;
	Myint i1PadStart   , i1PadEnd   , i2PadStart   , i2PadEnd   , i3PadStart   , i3PadEnd ;

	// origin of local grid
	Myfloat Orig1, Orig2, Orig3 ;

	// grid offset (with respect to global grid)
	Myint i1OffsetGlob, i2OffsetGlob, i3OffsetGlob ;

	// proc id of neighbours (for subdomain decomposiiton)
	Myint i1ProcIdStart, i1ProcIdEnd, i2ProcIdStart, i2ProcIdEnd, i3ProcIdStart, i3ProcIdEnd ;

	// MPI data type used to exchanges halos
	MPI_Datatype i1HaloDataType, i2HaloDataType, i3HaloDataType ;

	// grid padding
	virtual void padGridn1(void) ;
	virtual void padGridn2(void) ;
	virtual void padGridn3(void) ;

} ;

} // namespace hpcscan

#endif
