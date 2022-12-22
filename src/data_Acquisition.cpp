//-------------------------------------------------------------------------------------------------------
// This class handles the data acquisition:
// - Sources and receivers position
// - Source excitation (apply wavelet excitation at sources locations)
// - Storage of traces (data recorded at receivers locations)
//-------------------------------------------------------------------------------------------------------

#include "data_Acquisition.h"

#include <cassert>
#include <cmath>
#include <cstring> // for C memcpy
#include <fcntl.h> // for C open file operations
#include <unistd.h> // for C read/write file operations

#include "config.h"
#include "output_report.h"
#include "global.h" 

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

DataAcquisition::DataAcquisition(void)
{
    printDebug(MID_DEBUG, "IN DataAcquisition::DataAcquisition");

    nSrc = 0;
    nRec = 0;
    nSampleInTraceLoc = 0;
    total_nSampleInTrace = 0;
    nRecLoc=0;

    srcIdx = NULL;
    trace = NULL;

    printDebug(MID_DEBUG, "OUT DataAcquisition::DataAcquisition");
}

//-------------------------------------------------------------------------------------------------------

DataAcquisition::~DataAcquisition(void)
{
	printDebug(MID_DEBUG, "IN DataAcquisition::~DataAcquisition");

	if(nSrc != 0) delete[] srcIdx ;
    if(nRec != 0) delete[] tab_recIdxLoc ;
    if(nRec != 0) delete[] tab_idxglob ;

    delete[] trace;

	printDebug(MID_DEBUG, "OUT DataAcquisition::~DataAcquisition");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code DataAcquisition::initialize(DataAcqui_type acquiType, Grid& coefGrid, Myint nt)
{
	printDebug(MID_DEBUG, "In DataAcquisition::initialize") ;

    print_blank() ;
    printInfo(MASTER, " * Initialize Data Acquisition *");

    this->nt=nt;
    // get sources and receivers positions
    if (acquiType == ACQUI_BUILTIN)
    {
        printInfo(MASTER, " Acquisition type", "BUILTIN");

        // get indexes of inner points in local grid
        Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End;
        coefGrid.getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
        
        // get indexes of inner points in global grid
        Myint64 i1StartGlobal, i1EndGlobal, i2StartGlobal, i2EndGlobal, i3StartGlobal, i3EndGlobal;
        coefGrid.getGridIndexGlobal(INNER_POINTS, &i1StartGlobal, &i1EndGlobal, &i2StartGlobal, &i2EndGlobal, &i3StartGlobal, &i3EndGlobal);

        // get local grid size
        Myint64 n1 = coefGrid.n1;
        Myint64 n2 = coefGrid.n2;
        Myint64 n3 = coefGrid.n3;

        // get offset of local grid local
        Myint i1OffsetGlobInner = coefGrid.getOffsetGlobal(AXIS1);
        Myint i2OffsetGlobInner = coefGrid.getOffsetGlobal(AXIS2);
        Myint i3OffsetGlobInner = coefGrid.getOffsetGlobal(AXIS3);

        // Define source positions
        // Set one source located one point below top (n1)
        {           
            // For 1D: only one point
            // For 2D: in the middle of n2
            // For 3D: in the middle of n2 x n3

            Myint64 i1_glob = 1;
            Myint64 i2_glob = (i2EndGlobal - i2StartGlobal + 1) / 2;
            Myint64 i3_glob = (i3EndGlobal - i3StartGlobal + 1) / 2;

            if(coefGrid.isInMyDomain(i1_glob, i2_glob, i3_glob))
            {
                nSrc = 1;
                srcIdx = new Myint64[nSrc];
                
                Myint64 i1 = i1Start + i1_glob - i1OffsetGlobInner;
                Myint64 i2 = i2Start + i2_glob - i2OffsetGlobInner;
                Myint64 i3 = i3Start + i3_glob - i3OffsetGlobInner;

                srcIdx[0] = i1 + i2 * n1 + i3 * n2 * n1 ;
                
            }
            else
            {
                nSrc = 0;
            }

            Myint nsrc_glob;
            MPI_Reduce(&nSrc, &nsrc_glob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            printInfo(MASTER, " Number of sources", nsrc_glob);
        }

        // Define receiver positions
        // Set one plane a receivers located one point below top (n1)

        {
            // count the number of receivers globally
            Myint iRec = 0;

            // count how many reciever there is in this process
            Myint iRec_Loc = 0;

            if (coefGrid.dim == DIM1)
            {
                Myint64 i1_global = 1;
                Myint64 i2_global = 0;
                Myint64 i3_global = 0;

                if (coefGrid.isInMyDomain(i1_global, i2_global, i3_global))
                {
                    // For 1D: only one point
                    nRecLoc = 1;
                    tab_idxglob = new Myint[nRecLoc];
                    tab_recIdxLoc = new Myint64[nRecLoc];

                    Myint64 i1 = i1Start + i1_global - i1OffsetGlobInner;
                    Myint64 i2 = i2Start + i2_global - i2OffsetGlobInner;
                    Myint64 i3 = i3Start + i3_global - i3OffsetGlobInner;

                    printDebug(FULL_DEBUG, " i1", i1);
                    printDebug(FULL_DEBUG, " i2", i2);
                    printDebug(FULL_DEBUG, " i3", i3);

                    tab_recIdxLoc[iRec_Loc] = i1 + i2 * n1 + i3 * n2 * n1;
                    tab_idxglob[iRec_Loc] = iRec;
                }
            }
            else if (coefGrid.dim == DIM2)
            {
                // For 2D: covering n2
                //         one rec for every points for n2
                nRec = i2EndGlobal - i2StartGlobal + 1;

                // count local receptors
                for (Myint64 iRec2 = i2StartGlobal; iRec2 <= i2EndGlobal; iRec2++)
                {
                    Myint64 i1_Global = 1 ;
                    Myint64 i2_Global = iRec2 - i2StartGlobal;
                    Myint64 i3_Global = 0;
                    if(coefGrid.isInMyDomain(i1_Global, i2_Global, i3_Global))
                    {
                        nRecLoc++;
                    }
                }
                tab_recIdxLoc = new Myint64[nRecLoc];
                tab_idxglob = new Myint[nRecLoc];

                for (Myint64 iRec2 = i2StartGlobal; iRec2 <= i2EndGlobal; iRec2++)
                {   
                    //define postion of receptor in global grid
                    Myint64 i1_global = 1;
                    Myint64 i2_global = iRec2 - i2StartGlobal;
                    Myint64 i3_global = 0;
                   
                    if(coefGrid.isInMyDomain(i1_global,i2_global,i3_global))
                    {
                        Myint64 i1 = i1Start + i1_global - i1OffsetGlobInner;
                        Myint64 i2 = i2Start + i2_global - i2OffsetGlobInner;
                        Myint64 i3 = i3Start + i3_global - i3OffsetGlobInner;

                        printDebug(FULL_DEBUG, " i1", i1);
                        printDebug(FULL_DEBUG, " i2", i2);
                        printDebug(FULL_DEBUG, " i3", i3);
                        
                        tab_recIdxLoc[iRec_Loc] = i1 + i2 * n1 + i3 * n2 * n1;
                        tab_idxglob[iRec_Loc] = iRec;

                        iRec_Loc++;
                    }
                    iRec++;
                }
            }
            else if (coefGrid.dim == DIM3)
            {
                // For 3D: covering n2 x n3
                //         one rec for every points for n2
                //         one rec out of 10 points for n3
                Myint64 n2Rec = i2EndGlobal - i2StartGlobal + 1;
                Myint64 n3Rec = (i3EndGlobal - i3StartGlobal) / 10 + 1;
                nRec = n2Rec * n3Rec;

                for (Myint64 iRec3 = i3StartGlobal; iRec3 <= i3EndGlobal; iRec3=iRec3+10)
                {
                    for (Myint64 iRec2 = i2StartGlobal; iRec2 <= i2EndGlobal; iRec2++)
                    {
                        Myint64 i1_Global = 1 ;
                        Myint64 i2_Global = iRec2 - i2StartGlobal;
                        Myint64 i3_Global = iRec3 - i3StartGlobal;
                        if(coefGrid.isInMyDomain(i1_Global, i2_Global, i3_Global))
                        {
                          nRecLoc++;
                        }
                    }
                }
                
                tab_recIdxLoc = new Myint64[(nRecLoc)];
                tab_idxglob = new Myint[nRecLoc];

                for (Myint64 iRec3 = i3StartGlobal; iRec3 <= i3EndGlobal; iRec3=iRec3+10)
                {
                    for (Myint64 iRec2 = i2StartGlobal; iRec2 <= i2EndGlobal; iRec2++)
                    {
                        Myint64 i1_Global = 1 ;
                        Myint64 i2_Global = iRec2 - i2StartGlobal;
                        Myint64 i3_Global = iRec3 - i3StartGlobal;
                        if(coefGrid.isInMyDomain(i1_Global, i2_Global, i3_Global))
                        {
                            Myint64 i1 = i1Start + i1_Global - i1OffsetGlobInner;
                            Myint64 i2 = i2Start + i2_Global - i2OffsetGlobInner;
                            Myint64 i3 = i3Start + i3_Global - i3OffsetGlobInner;

                            printDebug(FULL_DEBUG, " i1", i1);
                            printDebug(FULL_DEBUG, " i2", i2);
                            printDebug(FULL_DEBUG, " i3", i3);

                            tab_recIdxLoc[iRec_Loc] = i1 + i2 * n1 + i3 * n2 * n1;
                            tab_idxglob[iRec_Loc] = iRec;
                            iRec_Loc++;
                            
                        }
                        iRec++;
                    }
                }
            }

            printInfo(MASTER, " Number of receivers", nRec);
            printInfo(ALL, " Number of receivers locally", nRecLoc);

             // allocate the array to store the traces
            // but only the points that contain values (no 0)
            nSampleInTraceLoc = nRecLoc * nt;

            total_nSampleInTrace = 0;
            MPI_Reduce(&nSampleInTraceLoc,&total_nSampleInTrace,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);

            printInfo(MASTER, " Trace file size (MB)", total_nSampleInTrace * sizeof(Myfloat) / 1e6);
            trace = new Myfloat[nSampleInTraceLoc];
            for(int i = 0 ; i < nSampleInTraceLoc; i++)
                trace[i]=0.0;
                
        }
    }
    else
    {
        printError("In DataAcquisition::initialize, unsupported acquiType", acquiType) ;
        return (RTN_CODE_KO);
    }    

    // allocate grid to store collected traces

    printDebug(MID_DEBUG, "Out DataAcquisition::initialize") ;
    return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code DataAcquisition::appliSourceTerm(Grid &grid, Modeling_type modType, Myint it, Myfloat dt)
{
    printDebug(MID_DEBUG, "In DataAcquisition::appliSourceTerm");

    if (modType == FORWARD)
    {
        if (nSrc != 0)
        {
            // compute amplitude of Ricker wavelet
            Myfloat rickerFreq = Config::Instance()->freqMax / 2.5; // dominant frequency of Ricker
            Myfloat da = M_PI * rickerFreq;
            Myfloat t00 = 1.5 * sqrt(6.) / da;
            Myfloat tt = it * dt;
            Myfloat aa = M_PI * rickerFreq * (tt - t00);
            Myfloat a2 = pow(aa, 2.);
            Myfloat srcTerm = (1. - 2. * a2) * exp(-a2);

            // apply source term
            grid.grid_3d[srcIdx[0]] = srcTerm;
        }
    }
    else if (modType == ADJOINT)
    {
        // TO DO
    }
    else
    {
        printError("In DataAcquisition::appliSourceTerm, unsupported modType", modType) ;
        return (RTN_CODE_KO);
    }

    printDebug(MID_DEBUG, "Out DataAcquisition::appliSourceTerm");
    return (RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code DataAcquisition::recordTrace(Grid &grid, Modeling_type modType, Myint it, Myfloat dt)
{
    printDebug(MID_DEBUG, "In DataAcquisition::recordTrace");

    double t0 = MPI_Wtime();
    if (modType == FORWARD)
    {
        // loop on all receivers
        for (Myint64 iRec = 0; iRec < nRecLoc; iRec++)
        {
            Myint64 traceIdx = it + nt * iRec ;
            trace[traceIdx] = grid.grid_3d[tab_recIdxLoc[iRec]];
        }
        
    }
    else if (modType == ADJOINT)
    {
        // TO DO
    }
    else
    {
        printError("In DataAcquisition::recordTrace, unsupported modType", modType) ;
        return (RTN_CODE_KO);
    }
    double t1 = MPI_Wtime();
    record_trace_time = t1 - t0;
    printDebug(MID_DEBUG, "Out DataAcquisition::recordTrace");
    return (RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code DataAcquisition::writeTrace(string fileName)
{
    printDebug(MID_DEBUG, "In DataAcquisition::writeTrace");

    double t0 = MPI_Wtime() ;

    // prepare write data
    Myint out_file;
    Myint64 nb_wrote;
    string fileName_bin = fileName + ".trace.bin";
    printInfo(MASTER, "* write on disk\t", fileName_bin.c_str());
    printDebug(LIGHT_DEBUG, "* write on disk\t", fileName_bin.c_str());

    //gather nSampleInTraceLoc to get the size of trace from each procs
    Myint64 tab_nSampleInTrace[nMpiProc];

    double t0_Gather_comm_time = MPI_Wtime();
    MPI_Gather(&nSampleInTraceLoc, 1, MPI_LONG_LONG, tab_nSampleInTrace, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    double t1_Gather_comm_time = MPI_Wtime();
    trace_gather_time += t1_Gather_comm_time-t0_Gather_comm_time;

    // gather nRecLoc from each other procs on proc 0
    // it gives us tab_nRec_procs
    Myint tab_nRec_procs [nMpiProc];
    t0_Gather_comm_time = MPI_Wtime();
    MPI_Gather(&nRecLoc, 1, MPI_INT, tab_nRec_procs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    t1_Gather_comm_time = MPI_Wtime();
    trace_gather_time += t1_Gather_comm_time - t0_Gather_comm_time;


    //prepare arrays for next gatherv
    Myint tab_deplacement [nMpiProc];
    tab_deplacement[0]=0;
    for(int i=1;i<nMpiProc;i++)
        tab_deplacement[i] = tab_deplacement[i-1] + tab_nRec_procs[i-1];

    Myint AllIdxrecglob [nRec];
    double t0_Gatherv_time = MPI_Wtime();
    MPI_Gatherv(tab_idxglob,nRecLoc,MPI_INT, AllIdxrecglob, tab_nRec_procs, tab_deplacement,MPI_INT, 0, MPI_COMM_WORLD);
    double t1_Gatherv_time = MPI_Wtime();
    trace_gather_time += t1_Gatherv_time - t0_Gatherv_time;

    //prepare arrays for next gatherv
    int tab_recv_count [nMpiProc] ; 
    for (int i=0; i < nMpiProc;i++)
        tab_recv_count[i]= tab_nSampleInTrace[i];

    int tab_displacement [nMpiProc];
    tab_displacement[0]=0;
    for(int i=1;i<nMpiProc;i++)
        tab_displacement[i] = tab_displacement[i-1] + tab_nSampleInTrace[i-1];
        

    //alloc array to store all traces
    Myfloat* all_trace_temp = new Myfloat[total_nSampleInTrace];
    for(int i =0;i<total_nSampleInTrace;i++)
        all_trace_temp[i]=0;

    t0_Gatherv_time = MPI_Wtime();
    MPI_Gatherv(trace, nSampleInTraceLoc, MPI_FLOAT, all_trace_temp, tab_recv_count, tab_displacement, MPI_FLOAT, 0, MPI_COMM_WORLD);
    t1_Gatherv_time = MPI_Wtime();
    trace_gather_time += t1_Gatherv_time - t0_Gatherv_time;
    
    // traces are gathered but in wrong order
    if(myMpiRank == 0){

        Myfloat all_traces[total_nSampleInTrace];
        for(int i=0;i<total_nSampleInTrace;i++)
            all_traces[i] = 0.0;

        // copy data at correct location
        Myint n_wrote=0;
        for(Myint i = 0 ; i< nRec; i++){
            memcpy(&all_traces[AllIdxrecglob[i]*nt], &all_trace_temp[n_wrote], nt * sizeof(float));
            n_wrote += nt;
        }

        double t0_trace_write = MPI_Wtime();

        out_file = open(fileName_bin.c_str(), O_WRONLY | O_APPEND | O_CREAT, S_IRUSR | S_IWUSR );
        if(out_file == -1){
            fprintf(stderr, "Error is %s",strerror( errno ));
            printError("Error while opening file on write trace NEC ") ;
            return (RTN_CODE_KO);
        }
        // write data
        Myint64 total_to_write = total_nSampleInTrace * sizeof(Myfloat);
        Myint nb_points_wrote = 0;
        do{
            nb_wrote=::write(out_file, &all_traces[nb_points_wrote],  (total_to_write > BLOCK_WRITE_SIZE ? BLOCK_WRITE_SIZE : total_to_write));
            if(nb_wrote < 0){
                printError("Error while writing trace NEC ");
                return (RTN_CODE_KO);
            }
            total_to_write -= nb_wrote;
            nb_points_wrote += nb_wrote / sizeof(Myfloat);
        }while(total_to_write > 0);


        close(out_file);

        // write info
        string fileName_info = fileName + ".trace.info";
        printDebug(LIGHT_DEBUG, "* write on disk\t", fileName_info.c_str());
        ofstream out_file2(fileName_info);
        out_file2 << total_nSampleInTrace / nRec << "\n";
        out_file2 << nRec << "\n";
        out_file2 << 1 << "\n";
        out_file2.close();

        double t1_trace_write = MPI_Wtime();
        trace_write_time = t1_trace_write - t0_trace_write;
    }

    delete[] all_trace_temp ;
    
    double t1 = MPI_Wtime();
    total_timeWriteTrace = t1 - t0 ;
    printInfo(MASTER, "* write time (s)", total_timeWriteTrace);
    printInfo(MASTER, "* write bwth (GB/s)", total_nSampleInTrace * sizeof(Myfloat)/total_timeWriteTrace/1e9);


    printDebug(MID_DEBUG, "Out DataAcquisition::writeTrace");
    return (RTN_CODE_OK);
}

} // namespace hpcscan
