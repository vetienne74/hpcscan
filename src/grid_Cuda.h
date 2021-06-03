
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode CUDA
// Derived class from Grid
// CUDA implementation (target GPU)
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_GRID_CUDA_H_
#define HPCSCAN_GRID_CUDA_H_

#include <string>

#include "mpi.h"

#include "grid.h"
#include "type_def.h"

namespace hpcscan {

class Grid_Cuda : public Grid
{
public:

	// constructor
	Grid_Cuda(Grid_type) ;

	// constructor
	Grid_Cuda(Grid_type, Dim_type, Myint64, Myint64, Myint64) ;

	// destructor
	~Grid_Cuda(void) ;

	// compute pressure with FD
	virtual Rtn_code computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder) ;

	// initialise grid index and MPI data structure
	virtual Rtn_code initializeGrid(void) ;

	// print info
	virtual void info(void) ;

	// write grid on disk
	virtual void write(string) ;

	// compute FD_D2 along N1
	virtual Rtn_code FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N2
	virtual Rtn_code FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N3
	virtual Rtn_code FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_LAPLACIAN
	virtual Rtn_code FD_LAPLACIAN(Point_type pType, const Grid&, Myint fdOrder) ;

	// fill grid with constant value
	virtual void fill(Point_type pointType, Myfloat val) ;

	// fill grid with predefined functions
	virtual void fill(Point_type pType, Func_type t1,  Func_type t2, Func_type t3,
			Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp) ;

	// get min and max in grid
	virtual Myfloat getMin(Point_type) ;
	virtual Myfloat getMax(Point_type) ;

	// fill array with constant value
	virtual void fillArray(Myfloat val) ;

	// copy array
	virtual void copyArray(const Grid& gridIn) ;

	// add array
	virtual void addArray(const Grid& gridIn1, const Grid& gridIn2) ;

	// multiply array
	virtual void multiplyArray(const Grid& gridIn1, const Grid& gridIn2) ;

	// add and update array
	virtual void addUpdateArray(const Grid& gridIn) ;

	// L1 error between this grid and another
	virtual Myfloat L1Err(Point_type pointType, const Grid& gridIn) const ;

	// collective (all MPI process) L1 error between this grid and another
	virtual Myfloat allProcL1Err(Point_type, const Grid&) const ;

	// Max error between this grid and another (point wise)
	virtual Myfloat maxErr(Point_type, const Grid&) const ;

	// update pressure
	virtual Rtn_code updatePressure(Point_type pType, const Grid& prcGrid,
			const Grid& coefGrid, const Grid& laplaGrid) ;

	// exchange one halo with MPI
	virtual Rtn_code exchangeHalo(MPI_comm_mode_type, Point_type pointType) ;

	// copy Grid device to host
	void copyGridDeviceToHost(Point_type) ;

	// copy Grid host to device
	void copyGridHostToDevice(Point_type) ;

	// apply boundary condition
	virtual Rtn_code applyBoundaryCondition(BoundCond_type boundCondType) ;

	// get sum of abs grid points
	virtual Myfloat getSumAbs(Point_type) const ;

	// get sum of abs diff of grid points
	virtual Myfloat getSumAbsDiff(Point_type, const Grid&) const;

	// send grid with MPI_Send
	virtual Rtn_code sendWithMPI(Myint64 nGridPoint, Myint procDestId) ;

	// receive grid with MPI_Recv
	virtual Rtn_code recvWithMPI(Myint64 nGridPoint, Myint procSrcId) ;

	// send and receive grid with MPI_Sendrecv
	virtual Rtn_code sendRecvWithMPI(const Grid& gridDest, Myint idSend, Myint idRecv, Myint64 nGridPoint) ;

protected:

	// pointer to device (GPU) memory
	// equivalent to grid_3d on the CPU side
	Myfloat * d_grid_3d ;

	// pointer to device (GPU) memory
	Myfloat * d_help_3d ;
	Myfloat * d_help_3d_2 ;

	// GPU block size (number of threads per block)
	Myint gpuBlkSize ;

	// GPU grid size (number of blocks per grid)
	Myint gpuGridSize ;

} ;

} // namespace hpcscan

#endif
