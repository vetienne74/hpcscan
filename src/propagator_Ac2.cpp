
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

namespace hpcscan
{

	//-------------------------------------------------------------------------------------------------------

	Propagator_Ac2::Propagator_Ac2(void)
	{
		printDebug(FULL_DEBUG, "In Propagator_Ac2::Propagator_Ac2");

		propaWaveEq = "Acoustic Wave Eq. 2nd order";
		propaKernelType = "Standard implementation";

		boundCondType = BOUND_COND_ANTI_MIRROR;
		fdOrder = Config::Instance()->fdOrder;
		snapInc = 0;
		tmax = 0.0;
		dt = 0.0;
		nt = 0;
		stableDt = 0.0;
		freqMax = 0.0;
		minVelocity = 0.0;
		maxVelocity = 0.0;

		coefGrid = nullptr;
		prnGrid = nullptr;
		prcGrid = nullptr;

		printDebug(FULL_DEBUG, "Out Propagator_Ac2::Propagator_Ac2");
	}

	//-------------------------------------------------------------------------------------------------------

	Rtn_code Propagator_Ac2::initialize(PropaInit_type propaInitType)
	{

		printDebug(FULL_DEBUG, "In Propagator_Ac2::initialize");

		//----------------
		// allocate grids
		//----------------
		const string gridMode = Config::Instance()->testMode;

		coefGrid = Grid_Factory::create(gridMode, GRID_LOCAL);
		coefGrid->initializeGrid();

		prnGrid = Grid_Factory::create(gridMode, GRID_LOCAL);
		prnGrid->initializeGrid();

		prcGrid = Grid_Factory::create(gridMode, GRID_LOCAL);
		prcGrid->initializeGrid();

		//--------------------
		// boundary condition
		//--------------------
		if (Config::Instance()->boundary.compare("FreeSurf") == 0)
		{
			boundCondType = BOUND_COND_ANTI_MIRROR;
		}
		else if (Config::Instance()->boundary.compare("None") == 0)
		{
			boundCondType = NO_BOUND_COND;
		}
		else
		{
			printError("In Propagator_Ac2::initialize, invalid boundary", Config::Instance()->boundary);
			return (RTN_CODE_KO);
		}

		//------------------
		// initialize grids
		//------------------

		if (propaInitType == EIGEN_MODE)
		{
			// define unit grids
			coefGrid->defineUnitGrid();
			prnGrid->defineUnitGrid();
			prcGrid->defineUnitGrid();

			//----------------
			// Set Vp = 1 m/s
			//----------------

			const Myfloat Vp = 1.0;
			minVelocity = Vp;
			maxVelocity = Vp;
		}
		else if (propaInitType == MODELING)
		{
			//-------------------------
			// read Vp model from disk
			// temporary in coefGrid
			//-------------------------

			if (Config::Instance()->modelVpFile.compare("UNSPECIFIED") == 0)
			{
				printError("In Propagator_Ac2::initialize, modelVpFile has not been specified");
				return (RTN_CODE_KO);
			}

			coef_read_total_time = 0.0;

			double t0_read_time = MPI_Wtime();
			coefGrid->read(INNER_POINTS, Config::Instance()->modelVpFile);
			double t1_read_time = MPI_Wtime();
			coef_read_total_time = t1_read_time - t0_read_time;
			i_o_read_time = coefGrid->i_o_pread_time;

			minVelocity = coefGrid->allProcGetMin(INNER_POINTS);
			maxVelocity = coefGrid->allProcGetMax(INNER_POINTS);
		}
		else
		{
			printError("In Propagator_Ac2::initialize, propaInitType not supported", propaInitType);
			return (RTN_CODE_KO);
		}

		// check Vp model
		if (minVelocity <= 0.0)
		{
			printError("In Propagator_Ac2::initialize, minVelocity <= 0.0");
			return (RTN_CODE_KO);
		}

		//--------
		// set dt
		//--------

		// compute stable dt = CFL * min(h) / vmax
		stableDt = getCFL() * coefGrid->getMinSpaceSampling() / maxVelocity;

		// if dt specified in config used that one
		if (Config::Instance()->dt > 0.0)
		{
			dt = Config::Instance()->dt;
		}
		else
		{
			// apply ratio CFL
			dt = stableDt * Config::Instance()->ratioCFL;
		}

		//-------------
		// set snapInc
		//-------------

		// if snapDt specified in config used that one
		if (Config::Instance()->snapDt > 0.0)
		{
			// make sure snapDt is a multiple of dt
			Myint incTmp = ceil(Config::Instance()->snapDt / dt);
			dt = Config::Instance()->snapDt / incTmp;
			snapInc = Config::Instance()->snapDt / dt;

			if ((Config::Instance()->snapDt - (snapInc * dt)) > dt / 2.0)
				snapInc++;
		}
		else
		{
			snapInc = Config::Instance()->snapInc;
		}

		//-----------------
		// set tmax and nt
		//-----------------

		if (Config::Instance()->tmax != UNSPECIFIED)
		{
			tmax = Config::Instance()->tmax;
			nt = (tmax / dt) + 1;
		}
		else
		{
			nt = Config::Instance()->nt;
			tmax = (nt - 1) * dt;
		}

		//---------------------
		// initialize coefGrid
		//---------------------
		if (propaInitType == EIGEN_MODE)
		{
			coefGrid->fill(INNER_POINTS, minVelocity * minVelocity * dt * dt);
		}
		else if (propaInitType == MODELING)
		{
			Myfloat *const w = coefGrid->grid_3d;
		
			Myint64 i1Start,i1End,i2Start,i2End,i3Start,i3End;
			Myint n1, n2 ,n3;
			n1 = coefGrid->n1 ;
			n2 = coefGrid->n2 ;
			n3 = coefGrid->n3 ;

			coefGrid->getGridIndex(INNER_POINTS,&i1Start,&i1End,&i2Start,&i2End,&i3Start,&i3End);

#pragma omp parallel for collapse(2)
			for(Myint64 i3 = i3Start ; i3 <= i3End ; i3++)
			{
				for(Myint64 i2 = i2Start ; i2 <= i2End ; i2++)
				{
					for(Myint64 i1 = i1Start ; i1 <= i1End ; i1++)
					{
						Myint64 index = i1 + i2 * n1 + i3 * n2 * n1;
						w[index] = w[index] * w[index] * dt * dt;
					}
				}
			}
		}
		
		coefGrid->write("coefGrid") ;
		i_o_write_coefgrid_time = coefGrid->grid_write_time;

		//--------------------------------
		// initialize prnGrid and prcGrid
		//--------------------------------

		if (propaInitType == EIGEN_MODE)
		{
			// initialize prn at t = -2 dt
			Myfloat timeSec = -2.0 * dt;
			if (initializeGrid(*prnGrid, propaInitType, timeSec) != RTN_CODE_OK)
			{
				printError("In Propagator_Ac2::initialize, initializeGrid not Ok");
				return (RTN_CODE_KO);
			}

			// initialize prc at t = -dt
			timeSec = -dt;
			if (initializeGrid(*prcGrid, propaInitType, timeSec) != RTN_CODE_OK)
			{
				printError("In Propagator_Ac2::initialize, initializeGrid not Ok");
				return (RTN_CODE_KO);
			}
		}
		else if (propaInitType == MODELING)
		{
			// initialize prn and prc grids to 0
			prnGrid->fill(ALL_POINTS, 0.0);
			prcGrid->fill(ALL_POINTS, 0.0);

			// set maximum frequency
			freqMax = Config::Instance()->freqMax ;
		}

		// initialize time measure at 0
		halo_comm_time = 0.0;
		compute_pressure_time = 0.0;
		boundary_condition_time = 0.0;

		printDebug(FULL_DEBUG, "Out Propagator_Ac2::initialize");
		return (RTN_CODE_OK);
	}

	//-------------------------------------------------------------------------------------------------------

	void Propagator_Ac2::info(void)
	{
		printDebug(FULL_DEBUG, "IN Propagator_Ac2::info");

		print_blank();
		printInfo(MASTER, " * Propagator parameters *");
		printInfo(MASTER, " Wave equation\t", propaWaveEq);
		printInfo(MASTER, " Kernel type\t", propaKernelType);
		if (boundCondType == NO_BOUND_COND)
			printInfo(MASTER, " Boundary type\t", "None");
		else if (boundCondType == BOUND_COND_ANTI_MIRROR)
			printInfo(MASTER, " Boundary type\t", "Free Surface");

		printInfo(MASTER, " FD space order\t", fdOrder);
		printInfo(MASTER, " CFL \t\t", getCFL());

		printInfo(MASTER, " Stable time step (s)", stableDt);
		printInfo(MASTER, " No. time steps\t", nt);
		printInfo(MASTER, " Time step (s)\t", dt);
		printInfo(MASTER, " Ratio Eff. / Stable dt", dt / stableDt);

		printInfo(MASTER, " Snapshot inc. (sec.)", (snapInc - 1) * dt);
		printInfo(MASTER, " Snapshot inc. (steps)", snapInc);

		printInfo(MASTER, "");
		printInfo(MASTER, " Max. frequency (Hz)", freqMax);
		printInfo(MASTER, " Min. period T (s)", 1.0 / freqMax);
		printInfo(MASTER, " Max. time (s)\t", tmax);
		printInfo(MASTER, " No. of periods\t", tmax * freqMax);
		printInfo(MASTER, " No. of steps / period", (1.0 / freqMax) / dt);

		printInfo(MASTER, "");
		printInfo(MASTER, " Min. velocity (m/s)", minVelocity);
		printInfo(MASTER, " Max. velocity (m/s)", maxVelocity);
		Myfloat lambdaMin = minVelocity / freqMax;
		printInfo(MASTER, " Min. lambda (m)", lambdaMin);
		printInfo(MASTER, " Min. grid spacing (m)", prnGrid->getMinSpaceSampling());

		Myfloat size1 = prnGrid->getMaxCoord(AXIS1) - prnGrid->getMinCoord(AXIS1);
		printInfo(MASTER, " No. of lambda axis1", size1 / lambdaMin);
		if (prnGrid->dim >= DIM2)
		{
			Myfloat size2 = prnGrid->getMaxCoord(AXIS2) - prnGrid->getMinCoord(AXIS2);
			printInfo(MASTER, " No. of lambda axis2", size2 / lambdaMin);
		}
		if (prnGrid->dim >= DIM3)
		{
			Myfloat size3 = prnGrid->getMaxCoord(AXIS3) - prnGrid->getMinCoord(AXIS3);
			printInfo(MASTER, " No. of lambda axis3", size3 / lambdaMin);
		}

		printInfo(MASTER, " No. grid pts / lambda", (minVelocity / freqMax) / prnGrid->getMinSpaceSampling());

		printDebug(FULL_DEBUG, "OUT Propagator_Ac2::info");
	}

	//-------------------------------------------------------------------------------------------------------

	Rtn_code Propagator_Ac2::computeWavefieldNextTimeStep(Grid &prnGridIn, Grid &prcGridIn)
	{

		printDebug(FULL_DEBUG, "In Propagator_Ac2::computeWavefieldNextTimeStep");

		double t0_exchange = MPI_Wtime ();
		// exchange halos between wavefield grid
		prcGridIn.myExchangeAll_halos(); // faster
		// prcGridIn.Grid::exchangeHalos(MPI_COMM_MODE_SENDRECV);
		// prcGridIn.exchangeHalos(MPI_COMM_MODE_SENDRECV); // overload in grid NEC
		double t1_exchange = MPI_Wtime ();
		halo_comm_time += t1_exchange - t0_exchange;

		// apply boundary 
		double t0_boundary = MPI_Wtime();
		prcGridIn.applyBoundaryCondition(boundCondType);
		double t1_boundary = MPI_Wtime();
		boundary_condition_time += t1_boundary - t0_boundary ;

		// compute with FD kernel
		double t0_pressure = MPI_Wtime();
		computePressureWithFD(prnGridIn, prcGridIn);
		double t1_pressure = MPI_Wtime();
		compute_pressure_time += t1_pressure - t0_pressure ;

		printDebug(FULL_DEBUG, "Out Propagator_Ac2::computeWavefieldNextTimeStep");
		return (RTN_CODE_OK);
	}

	//-------------------------------------------------------------------------------------------------------

	Rtn_code Propagator_Ac2::computePressureWithFD(Grid &prnGridIn, Grid &prcGridIn)
	{

		printDebug(FULL_DEBUG, "In Propagator_Ac2::computePressureWithFD");

		Rtn_code rtnCode = prnGridIn.computePressureWithFD(prcGridIn, *coefGrid, fdOrder);

		printDebug(FULL_DEBUG, "Out Propagator_Ac2::computePressureWithFD");
		return (rtnCode);
	}

	//-------------------------------------------------------------------------------------------------------

	Myfloat Propagator_Ac2::getCFL(void)
	{
		printDebug(FULL_DEBUG, "In Propagator_Ac2::getCFL");

		// compute dt from Lines 1999
		// v dt / h < sqrt (a1 / a2)
		// a1 = sum of FD coef for time discretization
		// a2 = sum of FD coef for space discretization

		// get a1
		Myfloat64 a1 = getSumAbsFD_D2Coef(2);
		printDebug(MID_DEBUG, " a1 \t\t", a1);

		// get a2
		Myfloat64 a2 = getSumAbsFD_D2Coef(fdOrder);
		printDebug(MID_DEBUG, " a2 \t\t", a2);
		if (coefGrid == nullptr)
		{
			printInfo(MASTER, "In Propagator_Ac2::getCFL, coefGrid == nullptr");
			return (-1.0);
		}
		a2 *= coefGrid->dim;
		printDebug(MID_DEBUG, " a2 * dim \t\t", a2);

		// CFL = sqrt(a1 / a2)
		Myfloat CFL = sqrt(a1 / a2);

		printDebug(FULL_DEBUG, "Out Propagator_Ac2::getCFL");
		return (CFL);
	}

	Rtn_code Propagator_Ac2::initializeGrid(Grid &gridIn, PropaInit_type propaInitType, Myfloat timeSec)
	{
		printDebug(FULL_DEBUG, "In Propagator_Ac2::initializeGrid");

		if (propaInitType == EIGEN_MODE)
		{

			const Myint nMode1 = ceil(Config::Instance()->param1);
			const Myint nMode2 = ceil(Config::Instance()->param2);
			const Myint nMode3 = ceil(Config::Instance()->param3);

			Myfloat64 val1 = 0.0, val2 = 0.0, val3 = 0.0, val4 = 0.0, gamma = 0.0;

			//===========================================================
			// 1D CASE
			//===========================================================
			if (Config::Instance()->dim == DIM1)
			{
				// printInfo(MASTER, " # modes along axis1", nMode1) ;
				if (nMode1 < 1)
				{
					printError("nMode1 < 1");
					return (RTN_CODE_KO);
				}

				// parameters
				val1 = PI * nMode1;
				val2 = 1.0;
				val3 = 1.0;
				gamma = nMode1;
			}

			//===========================================================
			// 2D CASE
			//===========================================================
			else if (Config::Instance()->dim == DIM2)
			{

				// printInfo(MASTER, " # modes along axis1", nMode1) ;
				if (nMode1 < 1)
				{
					printError("nMode1 < 1");
					return (RTN_CODE_KO);
				}
				// printInfo(MASTER, " # modes along axis2", nMode2) ;
				if (nMode2 < 1)
				{
					printError("nMode2 < 1");
					return (RTN_CODE_KO);
				}

				// initialize parameters
				val1 = PI * nMode1;
				val2 = PI * nMode2;
				val3 = 1.0;
				gamma = sqrt(nMode1 * nMode1 + nMode2 * nMode2);
			}

			//===========================================================
			// 3D CASE
			//===========================================================
			else if (Config::Instance()->dim == DIM3)
			{

				// printInfo(MASTER, " # modes along axis1", nMode1) ;
				if (nMode1 < 1)
				{
					printError("nMode1 < 1");
					return (RTN_CODE_KO);
				}
				// printInfo(MASTER, " # modes along axis2", nMode2) ;
				if (nMode2 < 1)
				{
					printError("nMode2 < 1");
					return (RTN_CODE_KO);
				}
				// printInfo(MASTER, " # modes along axis3", nMode3) ;
				if (nMode3 < 1)
				{
					printError("nMode3 < 1");
					return (RTN_CODE_KO);
				}

				// initialize parameters
				val1 = PI * nMode1;
				val2 = PI * nMode2;
				val3 = PI * nMode3;
				gamma = sqrt(nMode1 * nMode1 + nMode2 * nMode2 + nMode3 * nMode3);
			}

			// initialize at time = timeSec
			val4 = sin(PI * timeSec * gamma);
			gridIn.fill(ALL_POINTS, FUNC_SINE, FUNC_SINE, FUNC_SINE, val1, val2, val3, val4);

			// set max frequency
			// TO DO should be done once
			freqMax = gamma / 2.0;
		}
		else
		{
			printError("In Propagator_Ac2::initializeGrid, propaInitType not supported", propaInitType);
			return (RTN_CODE_KO);
		}

		printDebug(FULL_DEBUG, "Out Propagator_Ac2::initializeGrid");
		return (RTN_CODE_OK);
	}

} // namespace hpcscan
