
//-------------------------------------------------------------------------------------------------------
// * 2nd order acoustic wave equation with constant density = 1
//   d2P /dt2 = Vp^2 Laplacian(P) + S
//
//   where
//   P is pressure wavefield (x,t)
//   Vp is the P-wave velocity (x)
//   S is the source term (x.t)
//
// * Discretisation with FD
//   2nd order in time
//   variable order in space from 2 to 16
//
//   [P(n+1) - 2 P(n) + P(n-1)] / (dt^2) = Vp^2 * Laplacian(P(n)) + S(n)
//
// * Implementation with 2 arrays for P and one array for physical coef.
//
//   pn = coef * LAP(pc) + dt^2 src + 2 pc - pn
//
//   with coef = vp^2 * dt^2
//   note at each time step pn and pc are swapped
//
// * Eigen mode
//   In a unit cube, with free surface boundary condition at all egdes,
//   there exist an analytical solution:
//
//   In 1D: p = sin (pi x1 N1) . sin (pi gamma t)
//          with gamma = N1
//
//   In 2D: p = sin (pi x1 N1) . sin (pi x2 N2) . sin (pi Gamma t)
//          with gamma = sqrt(N1^2 + N2^2)
//
//   In 3D: p = sin (pi x1 N1) . sin (pi x2 N2) . sin (pi x3 N3) . sin (pi Gamma t)
//          with gamma = sqrt(N1^2 + N2^2 + N3^2)
//
//   N1, N2, N3 are the number of modes (integer >=1) along x1,x2 and x2 respectively
//-------------------------------------------------------------------------------------------------------

#include "propagator_Ac2.h"

#include <algorithm>

#include "config.h"
#include "constant.h"
#include "fdm.h"
#include "grid_Factory.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Propagator_Ac2::Propagator_Ac2(void)
		{
	printDebug(FULL_DEBUG, "In Propagator_Ac2::Propagator_Ac2") ;

	propaWaveEq     = "Acoustic Wave Eq. 2nd order" ;
	propaKernelType = "Standard implementation" ;

	fdOrder     = Config::Instance()->fdOrder ;
	snapInc     = 0 ;
	tmax        = 0.0 ;
	dt          = 0.0 ;
	nt          = 0 ;
	stableDt    = 0.0 ;
	maxFreq     = 0.0 ;
	minVelocity = 0.0 ;
	maxVelocity = 0.0 ;

	coefGrid    = nullptr ;
	prnGrid     = nullptr ;
	prcGrid     = nullptr ;

	printDebug(FULL_DEBUG, "Out Propagator_Ac2::Propagator_Ac2") ;
		}

//-------------------------------------------------------------------------------------------------------

Rtn_code Propagator_Ac2::initialize(PropaInit_type propaInitType)
{

	printDebug(FULL_DEBUG, "In Propagator_Ac2::initialize") ;

	// allocate grids
	const string gridMode = Config::Instance()->testMode ;

	coefGrid = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	coefGrid->initializeGrid();

	prnGrid = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	prnGrid->initializeGrid();

	prcGrid = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	prcGrid->initializeGrid();

	if (propaInitType == EIGEN_MODE)
	{
		// define unit grids
		coefGrid->defineUnitGrid() ;
		prnGrid->defineUnitGrid() ;
		prcGrid->defineUnitGrid() ;

		//----------------
		// Set Vp = 1 m/s
		//----------------

		const Myfloat Vp = 1.0 ;
		minVelocity = Vp ;
		maxVelocity = Vp ;

		//--------
		// set dt
		//--------

		// compute stable dt = CFL * h / vmax
		stableDt = getCFL() * coefGrid->d1 / maxVelocity ;

		// if dt specified in config used that one
		if (Config::Instance()->dt > 0.0)
		{
			dt = Config::Instance()->dt ;
		}
		else
		{
			// apply ratio CFL
			dt = stableDt * Config::Instance()->ratioCFL ;
		}

		//-------------
		// set snapInc
		//-------------

		// if snapDt specified in config used that one
		if (Config::Instance()->snapDt > 0.0)
		{
			// make sure snapDt is a multiple of dt
			Myint incTmp = ceil(Config::Instance()->snapDt / dt) ;
			dt = Config::Instance()->snapDt / incTmp ;
			snapInc = Config::Instance()->snapDt / dt ;

			if ((Config::Instance()->snapDt - (snapInc * dt)) > dt/2.0) snapInc++ ;
		}
		else
		{
			snapInc = Config::Instance()->snapInc ;
		}

		//-----------------
		// set tmax and nt
		//-----------------

		if (Config::Instance()->tmax != UNSPECIFIED)
		{
			tmax = Config::Instance()->tmax ;
			nt   = (tmax / dt) + 1 ;
		}
		else
		{
			nt   = Config::Instance()->nt ;
			tmax = (nt-1) * dt ;
		}

		//---------------------
		// initialize coefGrid
		//---------------------

		coefGrid->fill(INNER_POINTS, Vp*Vp*dt*dt) ;

		//--------------------------------
		// initialize prnGrid and prcGrid
		//--------------------------------

		// initialize prn at t = -2 dt
		Myfloat timeSec = -2.0 * dt ;
		if (initializeGrid(*prnGrid, propaInitType, timeSec) != RTN_CODE_OK)
		{
			printError("In Propagator_Ac2::initialize, initializeGrid not Ok") ;
			return(RTN_CODE_KO) ;
		}

		// initialize prc at t = -dt
		timeSec = -dt ;
		if (initializeGrid(*prcGrid, propaInitType, timeSec) != RTN_CODE_OK)
		{
			printError("In Propagator_Ac2::initialize, initializeGrid not Ok") ;
			return(RTN_CODE_KO) ;
		}

	}
	else
	{
		printError("In Propagator_Ac2::initialize, propaInitType not supported", propaInitType) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "Out Propagator_Ac2::initialize") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Propagator_Ac2::info(void)
{
	printDebug(FULL_DEBUG, "IN Propagator_Ac2::info");

	print_blank() ;
	printInfo(MASTER, " * Propagator parameters *") ;
	printInfo(MASTER, " Wave equation\t", propaWaveEq) ;
	printInfo(MASTER, " Kernel type\t", propaKernelType) ;

	printInfo(MASTER, " FD space order\t", fdOrder) ;
	printInfo(MASTER, " CFL \t\t", getCFL()) ;

	printInfo(MASTER, " Stable time step (s)", stableDt) ;
	printInfo(MASTER, " No. time steps\t", nt) ;
	printInfo(MASTER, " Time step (s)\t", dt) ;
	printInfo(MASTER, " Ratio Eff. / Stable dt", dt / stableDt) ;

	printInfo(MASTER, " Snapshot inc. (sec.)", (snapInc-1)*dt) ;
	printInfo(MASTER, " Snapshot inc. (steps)", snapInc) ;

	printInfo(MASTER, "") ;
	printInfo(MASTER, " Max. frequency (Hz)", maxFreq) ;
	printInfo(MASTER, " Min. period T (s)", 1.0/maxFreq) ;
	printInfo(MASTER, " Max. time (s)\t", tmax) ;
	printInfo(MASTER, " No. of periods\t", tmax*maxFreq) ;
	printInfo(MASTER, " No. of steps / period", (1.0/maxFreq) / dt) ;

	printInfo(MASTER, "") ;
	printInfo(MASTER, " Min. velocity (m/s)", minVelocity) ;
	printInfo(MASTER, " Max. velocity (m/s)", maxVelocity) ;
	Myfloat lambdaMin = minVelocity/maxFreq ;
	printInfo(MASTER, " Min. lambda (m)", lambdaMin) ;
	printInfo(MASTER, " Min. grid spacing (m)", prnGrid->getMinSpaceSampling()) ;

	Myfloat size1 = prnGrid->getMaxCoord(AXIS1) - prnGrid->getMinCoord(AXIS1) ;
	printInfo(MASTER, " No. of lambda axis1", size1/lambdaMin) ;
	if (prnGrid->dim >= DIM2)
	{
		Myfloat size2 = prnGrid->getMaxCoord(AXIS2) - prnGrid->getMinCoord(AXIS2) ;
		printInfo(MASTER, " No. of lambda axis2", size2/lambdaMin) ;
	}
	if (prnGrid->dim >= DIM3)
	{
		Myfloat size3 = prnGrid->getMaxCoord(AXIS3) - prnGrid->getMinCoord(AXIS3) ;
		printInfo(MASTER, " No. of lambda axis3", size3/lambdaMin) ;
	}

	printInfo(MASTER, " No. grid pts / lambda", (minVelocity/maxFreq)/prnGrid->getMinSpaceSampling()) ;

	printDebug(FULL_DEBUG, "OUT Propagator_Ac2::info");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Propagator_Ac2::computeWavefieldNextTimeStep(Grid& prnGridIn, Grid& prcGridIn)
{

	printDebug(FULL_DEBUG, "In Propagator_Ac2::computeWavefieldNextTimeStep") ;

	// exchange halos between wavefield grid
	prcGridIn.exchangeHalos(MPI_COMM_MODE_SENDRECV) ;

	// apply boundary condition
	prcGridIn.applyBoundaryCondition(BOUND_COND_ANTI_MIRROR) ;

	// compute with FD kernel
	computePressureWithFD(prnGridIn, prcGridIn) ;

	printDebug(FULL_DEBUG, "Out Propagator_Ac2::computeWavefieldNextTimeStep") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Propagator_Ac2::computePressureWithFD(Grid& prnGridIn, Grid& prcGridIn)
{

	printDebug(FULL_DEBUG, "In Propagator_Ac2::computePressureWithFD") ;

	Rtn_code rtnCode = prnGridIn.computePressureWithFD(prcGridIn, *coefGrid, fdOrder) ;

	printDebug(FULL_DEBUG, "Out Propagator_Ac2::computePressureWithFD") ;
	return(rtnCode) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Propagator_Ac2::getCFL(void)
{
	printDebug(FULL_DEBUG, "In Propagator_Ac2::getCFL") ;

	// compute dt from Lines 1999
	// v dt / h < sqrt (a1 / a2)
	// a1 = sum of FD coef for time discretization
	// a2 = sum of FD coef for space discretization

	// get a1
	Myfloat64 a1 = getSumAbsFD_D2Coef(2) ;
	printDebug(MID_DEBUG, " a1 \t\t", a1) ;

	// get a2
	Myfloat64 a2 = getSumAbsFD_D2Coef(fdOrder) ;
	printDebug(MID_DEBUG, " a2 \t\t", a2) ;
	if (coefGrid == nullptr)
	{
		printInfo(MASTER, "In Propagator_Ac2::getCFL, coefGrid == nullptr") ;
		return(-1.0) ;
	}
	a2 *= coefGrid->dim ;
	printDebug(MID_DEBUG, " a2 * dim \t\t", a2) ;

	// CFL = sqrt(a1 / a2)
	Myfloat CFL = sqrt(a1/a2) ;

	printDebug(FULL_DEBUG, "Out Propagator_Ac2::getCFL") ;
	return(CFL) ;
}

Rtn_code Propagator_Ac2::initializeGrid(Grid& gridIn, PropaInit_type propaInitType, Myfloat timeSec)
{
	printDebug(FULL_DEBUG, "In Propagator_Ac2::initializeGrid") ;

	if (propaInitType == EIGEN_MODE)
	{

		const Myint nMode1 = ceil(Config::Instance()->param1) ;
		const Myint nMode2 = ceil(Config::Instance()->param2) ;
		const Myint nMode3 = ceil(Config::Instance()->param3) ;

		Myfloat64 val1=0.0, val2=0.0, val3=0.0, val4=0.0, gamma=0.0 ;

		//===========================================================
		// 1D CASE
		//===========================================================
		if (Config::Instance()->dim == DIM1)
		{
			//printInfo(MASTER, " # modes along axis1", nMode1) ;
			if (nMode1 < 1) {
				printError("nMode1 < 1") ;
				return(RTN_CODE_KO) ;
			}

			// parameters
			val1  = PI * nMode1 ;
			val2  = 1.0 ;
			val3  = 1.0 ;
			gamma = nMode1 ;
		}

		//===========================================================
		// 2D CASE
		//===========================================================
		else if (Config::Instance()->dim == DIM2)
		{

			//printInfo(MASTER, " # modes along axis1", nMode1) ;
			if (nMode1 < 1) {
				printError("nMode1 < 1") ;
				return(RTN_CODE_KO) ;
			}
			//printInfo(MASTER, " # modes along axis2", nMode2) ;
			if (nMode2 < 1) {
				printError("nMode2 < 1") ;
				return(RTN_CODE_KO) ;
			}

			// initialize parameters
			val1  = PI * nMode1 ;
			val2  = PI * nMode2 ;
			val3  = 1.0 ;
			gamma = sqrt(nMode1*nMode1 + nMode2*nMode2) ;
		}

		//===========================================================
		// 3D CASE
		//===========================================================
		else if (Config::Instance()->dim == DIM3)
		{

			//printInfo(MASTER, " # modes along axis1", nMode1) ;
			if (nMode1 < 1) {
				printError("nMode1 < 1") ;
				return(RTN_CODE_KO) ;
			}
			//printInfo(MASTER, " # modes along axis2", nMode2) ;
			if (nMode2 < 1) {
				printError("nMode2 < 1") ;
				return(RTN_CODE_KO) ;
			}
			//printInfo(MASTER, " # modes along axis3", nMode3) ;
			if (nMode3 < 1) {
				printError("nMode3 < 1") ;
				return(RTN_CODE_KO) ;
			}

			// initialize parameters
			val1  = PI * nMode1 ;
			val2  = PI * nMode2 ;
			val3  = PI * nMode3 ;
			gamma = sqrt(nMode1*nMode1 + nMode2*nMode2 + nMode3*nMode3) ;
		}

		// initialize at time = timeSec
		val4  = sin(PI* timeSec * gamma) ;
		gridIn.fill(ALL_POINTS, FUNC_SINE, FUNC_SINE, FUNC_SINE, val1, val2, val3, val4) ;

		// set max frequency
		// TO DO should be done once
		maxFreq = gamma / 2.0 ;

	}
	else
	{
		printError("In Propagator_Ac2::initializeGrid, propaInitType not supported", propaInitType) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "Out Propagator_Ac2::initializeGrid") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan
