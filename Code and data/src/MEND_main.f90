#define USE_MPI
    PROGRAM MEND_main
    !MEND: Microbial-ENzyme Decomposition model-CNP
    !C  Author: GANGSHENG WANG @ ORNL & ZEHAO LV @ WHU
    !C  Environmental Sciences Division
    !C  Wuhan university & Oak Ridge National Laboratory
    !C  Oak Ridge, TN 37831-6301
    !C  EMAIL: WANGG@ORNL.GOV  lvzehao@whu.edu.cn
    !C  March, 2014
    !C  Updated: November 03, 2023
#ifdef USE_MPI                 
    USE mpi
#endif                                

    USE MOD_MEND_TYPE
    USE MOD_MEND, ONLY: fMEND_OBJ, sOUT_tscale, sOUT_ALL_tscale, sCOFI_PAR
    USE MOD_MEND, ONLY: sSOBOL_VAR_Read,sSOBOL_VAR_OBJ, subMEND_Files_Close,TEST
    USE MOD_OPT_TYPE
    USE MOD_OPT, ONLY: SCE, sWRITE_OPT_out_head
    USE MOD_MCMC, ONLY: MCMC
    USE MOD_USRFS, ONLY: iRandSeedGen, sort1, indexx, gasdev,gasdev0,fMSE, Array_2dto1d, Array_1dto2d
    USE MOD_USRFS, ONLY: Sec2HMS, nMonsbwDates,nDaysbwDates
    USE MOD_String, only: StrCompress

    IMPLICIT NONE

    !!--------------------------------
    TYPE(sSCE_PAR)  :: sSCE
    TYPE(sMEND_PAR) :: sPAR
    TYPE(sMEND_OUT) :: sOUT
    TYPE(sMEND_INI) :: sINI
    !!--------------------------------
       
    INTEGER i, j, id,nRun, eof ,nrow, ncol !!, iDay
    INTEGER iFpar0,iFpar, iFvar, iFobj

    REAL(KIND=8) fObj !!,fObj_cr(2)  !!upper and lower bound for COFI
    REAL(KIND=8), ALLOCATABLE :: xx(:), bestPar(:,:), bestObj(:)
    integer, ALLOCATABLE :: iwk(:),iwk0(:)

    CHARACTER(LEN = 1000) sRead 
    CHARACTER(LEN = 200) format100,format101,format510,format521

    INTEGER t_start,t_end,t_rate,t_elapse  !!time_start, time_end, time_elapsed
    INTEGER, DIMENSION(3):: tHMS  !!CALL subroutine "Sec2HMS" in "WGSfunc.f"
    INTEGER time_end(8)

    CHARACTER(LEN = 200) sFile,sFile_in, sFile_out,sFile_out0
    INTEGER nRow_skip, nday, nVAR, tstep,flag_avg, nPAR
    REAL(KIND=8), ALLOCATABLE:: dSIM_t01d(:)
    REAL(KIND=8), ALLOCATABLE:: dSIM_t0(:,:)
    REAL(KIND=8), ALLOCATABLE:: dSIM_t(:,:)
    REAL(KIND=8), ALLOCATABLE:: Sobol_obj(:), Sobol_obj_slave(:)

    INTEGER nSA, nSA0 !! number of parameter sets for Sobol Sensitivity Analysis
    REAL(KIND=8) rRead
    REAL(KIND=8), ALLOCATABLE:: xpar(:,:), xpar0(:,:)
    REAL(KIND=8), ALLOCATABLE:: xpar1d(:)
    REAL(KIND=8), ALLOCATABLE:: yvar(:)

    !! MPI Variables:
    integer pid, root_process, num_procs, ierr
#ifdef USE_MPI
    integer status(MPI_STATUS_SIZE)
    integer islave
    INTEGER, PARAMETER :: send_tag(4)   = (/2001, 2002, 2003, 2004/)
    INTEGER, PARAMETER :: return_tag(4) = (/2021, 2022, 2023, 2024/)
#endif


    !!---------------------------------------------------------
    !! TEST: BEGIN
    !!    CALL TEST()
    !! TEST: END
    !!---------------------------------------------------------
    root_process = 0

#ifdef USE_MPI
    call MPI_INIT (ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD, pid, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

    sINI%pid    = pid
#else
    pid         = 0
    sINI%pid    = 0
    num_procs   = 1
#endif
    sINI%nprocs = num_procs
    open(unit=1010,file='hello.txt')
    if(pid == 0) then
        print*, "*---------------------------------------------*"
        print*, "* Microbial-ENzyme Decomposition (MEND) Model *"
        print*, "*--------------C-N-P Coupled MEND-------------*"
        print*, "*----MPI Parallel Version;November,03,2023----*"
        print*, "*-wang.gangsheng@gmail.com,lvzehao@whu.edu.cn-*"
        print*, "*---------------------------------------------*"

    end if

    !    write (*,*) '>>MEND RUN BEGINS:'


    CALL system_clock(t_start,t_rate)          
    
    sINI%nOutStep = 24    ![h], output interval, 24 h = daily

    CALL MENDIN(sSCE,sINI) !read Model parameters (initial value, lower and upper bounds), input data

    nPAR = sSCE%nPar
    Allocate(xx(nPAR))
    xx = sSCE%a
    
    !Select Model Run (0: run model; 1: optimization; 2/5: COFI uncertainty; 3: Sobol Sensitivity; 4: MCMC)
    SELECT CASE (sINI%iModel)                            
        !!==============================================================================================
        CASE (1) !Model Calibration/Parameter Optimization
            if (sSCE%nRun .gt. 0) then
                nRun = min(sSCE%nRun, 200)
            else
                nRun = 1
            end if
        
            ALLOCATE(yvar(sINI%nVARopt))  !!fobj for each variable

            if(pid.eq.root_process) then

                write(sSCE%iFout_end,*)"Parameter Optimization Results from Multiple Runs:"
                sRead = "    OBJ-1:"
                write(format510,*)"(/,",sSCE%nPar,"(a16),","' |  CRITERION'",",a10, I10)"
                write(format521,*)"(",nPar,"(f15.8,','),","' | '",",f10.4,",sINI%nVARopt,"f10.4)"
                write(format100,*)"(/a10,",nPar,"(a15))"

                write(sSCE%iFout_end,format510)sSCE%parName,sRead,sINI%nVARopt

                write(sSCE%iFout_ini,*) ">>MODEL CALIBRATION/OPTIMIZATION..."

                write(sSCE%iFout_ini,'(a10,I5)')"nRun= ", nRun
                write(sSCE%iFout_ini,'(a10,I5)')"nPar= ", sSCE%nPar
                write(sSCE%iFout_ini,'(a10,I5)')"nOpt= ", sSCE%nOpt
                write(format510,*)"(a12,",sSCE%nOpt,"(I3))"
                write(sSCE%iFout_ini,format510)"PAR_Opt[i]= ",sSCE%iOpt
                write(sSCE%iFout_ini,*)
                write(sSCE%iFout_ini,*)"Output File 1: model setup and log: ", trim(sINI%dirout)//"OPT_ini.out"
                write(sSCE%iFout_ini,*)"Output File 2: optimal parameters : ", trim(sINI%dirout)//"OPT_end.out"
                write(sSCE%iFout_ini,*)"Output File 3: all parameters     : ", trim(sINI%dirout)//"OPT_all.out"

                ALLOCATE(iwk(nRun))
                ALLOCATE(bestObj(nRun))
                ALLOCATE(bestPar(nRun, sSCE%nPar))
            end if

            !! head for OPT_out.ini
            CALL sWRITE_OPT_out_head(sSCE)                

            if(num_procs .eq. 1) then
                do i = 1, nRun
                    CALL iRandSeedGen(sSCE%iseed)  !!see MOD_USRFS             
                    !!            sSCE%iseed = -1919936880  !! for DEBUG: see MOD_USRFS.F90, very small variance<0 in fSTDDEV (caused by identical xx), bug fixed on 11/16/2018
                    !!            sSCE%iseed = -1869941160  !!for DEBUG: Active Microbial Biomass (MBA <0), see MOD_MEND.F90 Line 1097, bug fixed on 11/20/2018
                    CALL srand(sSCE%iseed)  !!reinitialize the random number generator

                    CALL SCE(sSCE, sPAR, sINI, sOUT)                   
                    if (sINI%iError < 0) then
                        write(*,*)"Mass balance error: sINI%iError = ", sINI%iError       
                        stop
                    end if
                    bestObj(i) = sSCE % bestObj
                    bestPar(i,:) = sSCE % bestPar
                    write(sSCE%iFout_end,format521) bestPar(i,:),bestObj(i),sINI%rOBJ
                end do
            else !! num_procs > 1
#ifdef USE_MPI
                if(pid.eq.root_process) then

                    do j = 1, nRun
                        call MPI_RECV(islave, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                            return_tag(1), MPI_COMM_WORLD, status, ierr)  !!MPI_ANY_TAG
                        call MPI_RECV(yvar, sINI%nVARopt, MPI_REAL8, MPI_ANY_SOURCE, &
                            return_tag(2), MPI_COMM_WORLD, status, ierr)
                        call MPI_RECV(fObj, 1, MPI_REAL8, MPI_ANY_SOURCE, &
                            return_tag(3), MPI_COMM_WORLD, status, ierr)
                        call MPI_RECV(xx, nPAR, MPI_REAL8, MPI_ANY_SOURCE, &
                            return_tag(4), MPI_COMM_WORLD, status, ierr)

                        bestObj(islave) = fObj
                        bestPar(islave,1:nPAR) = xx(1:nPAR)
                        write(sSCE%iFout_end,format521) bestPar(islave,:),bestObj(islave),yvar(:)
                    end do
                else !!slave process
                    do i = pid, nRun, (num_procs-1)  !!root process is NOT used for computing
                        sSCE%iRun = i
                        CALL iRandSeedGen(sSCE%iseed)  !!see MOD_USRFS
                        sSCE%iseed = pid*sSCE%iseed

                        CALL srand(sSCE%iseed)  !!reinitialize the random number generator

                        CALL SCE(sSCE, sPAR, sINI, sOUT)
                        fObj = sSCE % bestObj
                        xx(1:nPAR) = sSCE % bestPar
                        yvar = sINI%rOBJ

                        call MPI_SEND(i, 1, MPI_INTEGER, root_process, return_tag(1), MPI_COMM_WORLD, ierr)
                        call MPI_SEND(yvar, sINI%nVARopt, MPI_REAL8, root_process, return_tag(2), MPI_COMM_WORLD, ierr)
                        call MPI_SEND(fObj, 1, MPI_REAL8, root_process, &
                            return_tag(3), MPI_COMM_WORLD, ierr)
                        call MPI_SEND(xx, nPAR, MPI_REAL8, root_process, &
                            return_tag(4), MPI_COMM_WORLD, ierr)
                    end do !!i = pid, nSA, num_procs

                end if !!(pid.eq.0)
#endif
            end if !!(num_procs .eq. 1)

            if(pid.eq.root_process) then
                CALL indexx(nRun, bestObj, iwk)     !rank best objective function value (OBF)   
                xx = bestPar(iwk(1),:)              !pick parameter values resulting in best OBF

            
                DEALLOCATE(bestPar)
                DEALLOCATE(bestObj)
                DEALLOCATE(iwk)
            end if
        
            DEALLOCATE(yvar)
        !!==============================================================================================
        CASE (2) !read parameter from SCE samples to select PAR for COFI
            if(pid.eq.0) then
                write(sSCE%iFout_ini,*)">>Parameter Uncertainty Quantification by COFI, ", &
                    "Multiple PAR Values are Provided in <COFIpar.dat>"

                sFile_in   = trim(sINI%dirinp)//'/COFIpar.dat'
                sFile_out0 = trim(sINI%dirout)//'COFIpar.out'

                write(sSCE%iFout_ini,*)"Input  File: All Parameters from optimization: ", sFile_in
                write(sSCE%iFout_ini,*)"Output File: COFI feasible parameter sets: ", sFile_out0

                CALL sCOFI_PAR(nPAR, sSCE%parName,sFile_in,sFile_out0,nSA)     
            end if
        !!==============================================================================================
        CASE (5) !read parameter from SCE samples to compute statistics of response variables
            nVAR = sINI%nOBS_tot
            !        print*,pid, sINI%nOBS_tot, nVAR

            if(pid.eq.0) then
                write(sSCE%iFout_ini,*)">>Parameter & Variable Uncertainty Quantification by COFI, ", &
                    "Multiple PAR Values are Provided in <COFIpar.dat>"
                !            iFpar0   = 301            !File to store parameter values for sensitivity/uncertainty analysis
                iFpar    = 302            !PAR file with fObj < fObj_cr
                iFvar    = 303            !File to store response variable outputs from sensitivity/uncertainty analysis
            
                sFile_in   = trim(sINI%dirinp)//'/COFIpar.dat'
                sFile_out0 = trim(sINI%dirout)//'COFIpar.out'
                sFile_out  = trim(sINI%dirout)//'COFIvar.out'
                write(sSCE%iFout_ini,*)"Input  File  : All Parameters from optimization: ", sFile_in
                write(sSCE%iFout_ini,*)"Output File 1: COFI feasible parameter sets: ", sFile_out0
                write(sSCE%iFout_ini,*)"Output File 2: COFI simulated variables: ", sFile_out

                CALL sCOFI_PAR(nPAR,sSCE%parName,sFile_in,sFile_out0, nSA)

                open(unit = iFvar, file = trim(sFile_out) , status = 'unknown')
                write(iFvar,'(4A10,I5,A)')"ID", "fobj0", "fobj", "var[1:",nVAR,"]"
                !            write(format101, *) "(", sINI%nOBS_tot, "E15.3)"
                write(format101, *) "(I10,2f10.4,", sINI%nOBS_tot, "E15.3)"
            end if !! pid.eq.0
#ifdef USE_MPI
            call MPI_BCAST (nSA, 1, MPI_INTEGER, root_process, MPI_COMM_WORLD, ierr)
#endif
            !        write(*,*)pid,nSA,nPAR
            ALLOCATE(xpar(nSA,nPAR))
            ALLOCATE(xpar1d(nSA*nPAR))
            ALLOCATE(bestObj(nSA))
            ALLOCATE(yvar(nVAR))

            if (pid.eq.0) then
                open(unit = iFpar, file = trim(sFile_out0) , status = 'old')
                read(iFpar,*)sRead   !!1st line
                read(iFpar,*)sRead   !!2nd line: parameter names
                do i = 1, nSA
                    read(iFpar, *)bestObj(i),xpar(i, 1:nPAR)
                end do

                CALL Array_2dto1d(nSA, nPAR, xpar, xpar1d)
            !                write(*,*)xpar1d(((nSA-1)*nPAR+1):(nSA*nPAR))
            end if !! pid.eq.0
#ifdef USE_MPI
            call MPI_BCAST (xpar1d, nSA*nPAR, MPI_REAL8, root_process, MPI_COMM_WORLD, ierr)
#endif
            if(num_procs .eq. 1) then
                do i = 1, nSA   !!root process is used for computing
                    xx(1:nPAR) = xpar1d(((i-1)*nPAR+1):(i*nPAR))
                    !                write(*,*)pid,i,xx
                    fObj = fMEND_OBJ(xx, sPAR, sINI, sOUT)           
                    write(iFvar, format101)i,bestObj(i),fobj,sINI%dSIM_opt(:,2)
                    write(*,'(A5,I6,A1,2f20.4)')"fObj[",i,"]",bestObj(i),fObj
                end do !!i = 1, nSA
            else !! num_procs > 1
#ifdef USE_MPI
                !            write(*,*)"pid=",pid,xpar1d(((nSA-1)*nPAR+1):(nSA*nPAR))
                if(pid.eq.0) then
                    do j = 1, nSA
                        call MPI_RECV(islave, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                            return_tag(1), MPI_COMM_WORLD, status, ierr)  !!MPI_ANY_TAG
                        call MPI_RECV(fObj, 1, MPI_REAL8, MPI_ANY_SOURCE, &
                            return_tag(2), MPI_COMM_WORLD, status, ierr)
                        call MPI_RECV(yvar, nVAR, MPI_REAL8, MPI_ANY_SOURCE, &
                            return_tag(3), MPI_COMM_WORLD, status, ierr)

                        write(iFvar, format101)islave,bestObj(islave),fObj,yvar
                        write(*,'(A5,I6,A1,2f20.4)')"fObj[",islave,"]",bestObj(islave),fObj
                    end do
                else !!slave process
                    do i = pid, nSA, (num_procs-1)  !!root process is NOT used for computing
                        xx(1:nPAR) = xpar1d(((i-1)*nPAR+1):(i*nPAR))
                        fObj = fMEND_OBJ(xx, sPAR, sINI, sOUT)
                        yvar = sINI%dSIM_opt(1:sINI%nOBS_tot,2)

                        !                    print*,pid,i,fObj,yvar

                        call MPI_SEND(i, 1, MPI_INTEGER, root_process, return_tag(1), MPI_COMM_WORLD, ierr)
                        call MPI_SEND(fObj, 1, MPI_REAL8, root_process, return_tag(2), MPI_COMM_WORLD, ierr)
                        call MPI_SEND(yvar, nVAR, MPI_REAL8, root_process, &
                            return_tag(3), MPI_COMM_WORLD, ierr)
                    end do !!i = pid, nSA, num_procs

                end if !!(pid.eq.0)
#endif
            end if !!(num_procs .eq. 1)


            DEALLOCATE(xpar1d)
            DEALLOCATE(xpar)
            DEALLOCATE(yvar)
            DEALLOCATE(bestObj)

            if (pid.eq.0) close(iFvar)
        !!==============================================================================================
        CASE (3) !Sobol Sensitivity: Calculate Objective Function Values; Parameters are Generated by JAVA MOEA Framework
            iFpar    = 301            !File to store parameter values for sensitivity/uncertainty analysis      
            iFobj    = 302            !File to store fObj

            nVAR=const_nPOOL*3 + const_nPOOL_MN
            nRow_skip=3;            !!_VAR_hour.out
            tstep=1
            flag_avg=1

            nday = nDaysbwDates(sINI%sDate_beg_sim,sINI%sDate_end_sim)
            nSA0 = sINI%iSA_range(2) - sINI%iSA_range(1) + 1 !! number of parameter sets for SA

            write(format100,*)"(I10,",nVAR,"(e20.6))"
            write(format101,*)"(I10,",nPAR,"(e20.6))"


            !        ALLOCATE(xpar1d(nSA*nPAR))

            ALLOCATE(dSIM_t0(nday,nVAR))
            ALLOCATE(dSIM_t01d(nday*nVAR))
            ALLOCATE(dSIM_t(nday,nVAR))

            ALLOCATE(Sobol_obj(nVAR))
            ALLOCATE(Sobol_obj_slave(nVAR))


            sFile_in = sINI%sFile_VAR_hour !!_VAR_hour.out
            sFile_out = sFile_in(1:(len(trim(sFile_in))-8))//"day.out"
            sFile_out0 = trim(sFile_out)//'0'

            !---------------------------
            if(pid.eq.0) then
                ALLOCATE(iwk0(nSA0))
                ALLOCATE(xpar0(nSA0,nPAR))

                write(sSCE%iFout_ini,*)">>Sobol Sensitivity Analysis, Multiple PAR Values", &
                    " Generated by MOEA Framework are Provided in <SOBOLpar.txt>"             
        
                open(unit = iFpar, file = trim(sINI%dirinp)//'/SOBOLpar.txt' , status = 'old')

                read(iFpar, '(A)')sRead  !!1st line: "Median/Mean of Parameters:"
                read(iFpar, '(A)')sRead  !!2nd line Par names
                read(iFpar, *)id,xx(1:nPAR)  !!3rd line: Median/Mean of Parameters    

                CALL sSOBOL_VAR_OBJ(0,xx,sPAR, sINI, sOUT, sFile_in,sFile_out, &
                    nRow_skip, tstep,flag_avg, &
                    nPAR, nVAR, nday,dSIM_t0,Sobol_obj)            

                CALL Array_2dto1d(nday, nVAR, dSIM_t0, dSIM_t01d)

                write(*,'(A5,I5,A)')"pid=",pid,"; Sobol Sensitivity: Completed 1st step with parameter median/mean values"

                write(sRead,*)trim(sINI%dirout),'SOBOLobj_',sINI%iSA_range(1),'-',sINI%iSA_range(2),'.out'    
                sFile = StrCompress(sRead)
                write(sSCE%iFout_ini,*)"Sobol Objective Function Value output file: ",sFile
                open(unit = iFobj,file = trim(sFile) , status = 'unknown')

                write(iFobj,"(A10,200A20)")"ID",(trim(sINI%Name_POOL(i))//"_C",i=1,const_nPOOL), &
                    (trim(sINI%Name_MNPOOL(i)),i=1,const_nPOOL_MN), &
                    (trim(sINI%Name_POOL(i))//"_N",i=1,const_nPOOL), &
                    (trim(sINI%Name_POOL(i))//"_CN",i=1,const_nPOOL)

                read(iFpar, '(A)')sRead  !!4th line: "Parameter Samples Generated"

                if(sINI%iSA_range(3) .eq. 0) then  !!continuous id with the first id = 1
                    do i=1, sINI%iSA_range(1) - 1  !! skip these lines
                        read(iFpar, *)rRead
                    end do
                end if

                j = 0
                DO i=sINI%iSA_range(1),sINI%iSA_range(2)
                    read(iFpar, *, iostat = eof)id, xx(1:nPAR)

                    IF (eof < 0) THEN !end of file has reached
                        print*,'eof=',eof
                        EXIT
                    ELSE IF (eof > 0) THEN !input error
                        print*,">>>Skip this Line!"
                    ELSE
                        j = j + 1
                        iwk0(j) = id
                        xpar0(j,1:nPAR) = xx(1:nPAR)
                    END IF
                END DO
                close(iFpar)                

                nSA = j  !!actual nSA <= sINI%iSA_range(2) - sINI%iSA_range(1) + 1
                write(*,'(A5,I5)')"nSA = ",nSA
            end if !! if(pid.eq.0)

#ifdef USE_MPI
            call MPI_BCAST(nSA, 1, MPI_INTEGER, root_process, MPI_COMM_WORLD, ierr)
#endif

            ALLOCATE(xpar1d(nSA*nPAR))

            if(pid.eq.0) then
                ALLOCATE(iwk(nSA))
                ALLOCATE(xpar(nSA,nPAR))

                iwk(1:nSA) = iwk0(1:nSA)
                xpar(1:nSA,1:nPAR) = xpar0(1:nSA,1:nPAR)
                CALL Array_2dto1d(nSA, nPAR, xpar, xpar1d)

                DEALLOCATE(iwk0)
                DEALLOCATE(xpar0)
            end if
            !---------------------------
#ifdef USE_MPI
            call MPI_BCAST(xpar1d, nSA*nPAR, MPI_REAL8, root_process, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(dSIM_t01d, nday*nVAR, MPI_REAL8, root_process, MPI_COMM_WORLD, ierr)
#endif

            if(num_procs .eq. 1) then
                do i = 1, nSA   !!root process is used for computing
                    xx(1:nPAR) = xpar1d(((i-1)*nPAR+1):(i*nPAR))
                    !                write(*,*)pid,i,xx
                    !                write(*,*)"dSIM_t0=",pid,dSIM_t0(nday,:)
                    CALL sSOBOL_VAR_OBJ(1,xx,sPAR, sINI, sOUT, sFile_in,sFile_out, &
                        nRow_skip, tstep,flag_avg, &
                        nPAR, nVAR, nday,dSIM_t0,Sobol_obj)                 

                    write(*,    format100)iwk(i),Sobol_obj
                    write(iFobj,format100)iwk(i),Sobol_obj
                end do !!i = pid, nSA, num_procs
            else !! num_procs > 1
#ifdef USE_MPI
                !        do i = (pid +1), nSA, num_procs  !!root process is used for computing
                if(pid.eq.0) then
                    do j = 1, nSA
                        call MPI_RECV(islave, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                            return_tag(1), MPI_COMM_WORLD, status, ierr)
                        call MPI_RECV(Sobol_obj_slave, nVAR, MPI_REAL8, MPI_ANY_SOURCE, &
                            return_tag(2), MPI_COMM_WORLD, status, ierr)

                        write(*,    format100)iwk(islave),Sobol_obj_slave
                        write(iFobj,format100)iwk(islave),Sobol_obj_slave
                    end do
                else !!slave process
                    CALL Array_1dto2d(dSIM_t01d, nday, nVAR, dSIM_t0)
                    do i = pid, nSA, (num_procs-1)  !!root process is NOT used for computing
                        xx(1:nPAR) = xpar1d(((i-1)*nPAR+1):(i*nPAR))
                        !                    write(*,*)pid,i,xx
                        !                    write(*,*)"dSIM_t0=",pid,dSIM_t0(nday,:)
                        CALL sSOBOL_VAR_OBJ(1,xx,sPAR, sINI, sOUT, sFile_in,sFile_out, &
                            nRow_skip, tstep,flag_avg, &
                            nPAR, nVAR, nday,dSIM_t0,Sobol_obj)

                        call MPI_SEND(i, 1, MPI_INTEGER, root_process, return_tag(1), MPI_COMM_WORLD, ierr)
                        call MPI_SEND(Sobol_obj, nVAR, MPI_REAL8, root_process, return_tag(2), MPI_COMM_WORLD, ierr)
                    end do !!i = pid, nSA, num_procs

                end if !!(pid.eq.0)
#endif
            end if !!(num_procs .eq. 1)


            DEALLOCATE(dSIM_t01d)
            DEALLOCATE(dSIM_t0)
            DEALLOCATE(dSIM_t)
            DEALLOCATE(xpar1d)
            DEALLOCATE(Sobol_obj)
            DEALLOCATE(Sobol_obj_slave)

            if(pid.eq.0) then
                DEALLOCATE(iwk)
                DEALLOCATE(xpar)
                close(iFobj)
            end if
        !!==============================================================================================
        CASE (6) !Sobol Sensitivity: Find unfinished Parameter sets
            if(pid.eq.0) then
                iFpar0   = 301            !File to store parameter values for sensitivity/uncertainty analysis
                iFobj    = 302            !File to store fObj
                iFpar    = 303

                nVAR=const_nPOOL*3 + const_nPOOL_MN
                nRow_skip=3;                !!_VAR_hour.out
                tstep=1
                flag_avg=1

                nday = nDaysbwDates(sINI%sDate_beg_sim,sINI%sDate_end_sim)
                nSA0 = sINI%iSA_range(2) - sINI%iSA_range(1) + 1 !! number of parameter sets for SA

                !        write(format101,*)"(2I10,",nVAR,"(e20.6))"
                write(format100,*)"(I15,",nPAR,"(e15.4))"

                ALLOCATE(xpar(nSA0,nPAR))
                ALLOCATE(iwk(nSA0))

                !---------------------------

                write(sSCE%iFout_ini,*)">>Sobol Sensitivity Analysis: Find unFinished PAR Sets in Previous Run"
                sFile_in = trim(sINI%dirinp)//'/SOBOLpar.txt'
                sFile_out = sFile_in(1:len(trim(sFile_in))-4)//"_Finished.txt"

                open(unit = iFpar, file = trim(sFile_in) , status = 'old')
                open(unit = iFpar0,file = trim(sFile_out) , status = 'unknown')

                read(iFpar, '(A)')sRead  !!1st line: "Median/Mean of Parameters:"
                write(iFpar0,*)sRead
                read(iFpar, '(A)')sRead  !!2nd line Par names
                write(iFpar0,'(A)')sRead  !!2nd line: head: par names
                read(iFpar, *) id,xx(1:nPAR)  !!3rd line: Median/Mean of Parameters
                write(iFpar0,format100)id,xx(1:nPAR)

                read(iFpar, '(A)')sRead  !!4th line: "Parameter Samples Generated"
                write(iFpar0,*)sRead

                do i=1, sINI%iSA_range(1) - 1  !! skip these lines
                    read(iFpar, *)rRead
                end do

                j = 0
                DO i=sINI%iSA_range(1),sINI%iSA_range(2)
                    read(iFpar, *, iostat = eof)id, xx(1:nPAR)  !!the id in this file are continuous numbers

                    IF (eof < 0) THEN !end of file has reached
                        print*,'eof=',eof
                        EXIT
                    ELSE IF (eof > 0) THEN !input error
                        print*,">>>Skip this Line!"
                    ELSE
                        j = j + 1
                        xpar(j,1:nPAR) = xx(1:nPAR)
                    END IF
                END DO
                close(iFpar)

                write(sRead,*)trim(sINI%dirinp),trim(sINI%SITE),'_SOBOLobj_',&
                    sINI%iSA_range(1),'-',sINI%iSA_range(2),'.out'
                sFile = StrCompress(sRead)
                write(sSCE%iFout_ini,*)"Input File 1: Original Sobol Parameter Sets: ",sFile_in
                write(sSCE%iFout_ini,*)"Input File 2: Sobol Objective function values that have been calculated: ",sFile
                write(sSCE%iFout_ini,*)"Output File: Sobol Parameter Sets that have NOT been used", &
                    " to calculate Objective function values: ", sFile_out
                !            open(iFobj,*)"Check Sobol Sensitivity output file: ",sFile
                open(unit = iFobj,file = trim(sFile) , status = 'old')
                read(iFobj,*)sRead  !! 1st line, head

                nrow = 0  !!number of lines of Obj
                DO
                    read(iFobj, *, iostat = eof) j, sRead

                    IF (eof < 0) THEN !end of file has reached
                        EXIT
                    ELSE IF (eof > 0) THEN !input error
                        print*,">>>Skip this Line!"
                    ELSE
                        nrow = nrow+1
                        iwk(nrow) = j
                    END IF
                END DO
                close(iFobj)

                CALL sort1(nrow,iwk(1:nrow))

                DO j=2,nrow
                    ierr = iwk(j) - iwk(j-1)
                    if(ierr > 1) then
                        do i=1, ierr - 1
                            id = iwk(j-1) + i
                            write(*,*)"Missing ID = ",id
                            write(iFpar0,format100)id, xpar(id-sINI%iSA_range(1)+1,1:nPAR)
                        end do
                    end if

                END DO

                close(iFpar0)

                DEALLOCATE(xpar)
                DEALLOCATE(iwk)

            end if !! if(pid.eq.0)
        !---------------------------

        !!==============================================================================================
        CASE (4) ! MCMC
            !! todo: Parallel for multiple chains
            CALL iRandSeedGen(sSCE%iseed)  !!see MOD_USRFS                       
            CALL srand(sSCE%iseed)  !!reinitialize the random number generator    
            write(sSCE%iFout_ini,*)"Output File: ", trim(sINI%dirout)//"MCMC.out"
            write (*, '(A20,I5,A15,I15)') 'MCMC Run Number = ', i, '; Random Seed = ', sSCE%iseed
        
            CALL MCMC(sSCE, sPAR, sINI, sOUT)     
            xx = sSCE % bestPar
    !!==============================================================================================
    END SELECT !!CASE (sINI%iModel)                  
    
    !RUN MEND MODEL USING (i)PARAMETER VALUES IN LAST LINE OF "MEND_namelist.nml" or (ii)"BEST" PARAMETER VALUES FROM OPTIMIZATION
    IF(sINI % iModel.lt.2 .or.sINI % iModel.eq.4) then        
        if (pid.eq.0) then
            write (*,*) '>>FINAL RUN <fMEND_OBJ> with GIVEN or BEST PAR...'
            ncol = 11
            nrow = sSCE%nPar/ncol
            if (nrow.lt.1) then
                write(*,'(20a12)')sSCE%parName
                write(*,'(20f12.6)')xx
            else
                do i=1,nrow
                    write(*,'(20a12)')sSCE%parName(((i-1)*ncol+1):(i*ncol))
                    write(*,'(20f12.6)')xx(((i-1)*ncol+1):(i*ncol))
                end do

                if(sSCE%nPar.gt.(nrow*ncol)) then
                    write(*,'(20a12)')sSCE%parName((nrow*ncol+1):sSCE%nPar)
                    write(*,'(20f12.6)')xx((nrow*ncol+1):sSCE%nPar)
                end if
            end if

            write(*,'(/,4A)')"Input-Date_Period = ",sINI%sDate_beg_all," -- ",sINI%sDate_end_all
            write(*,'(4A)')"Simulation_Period = ",sINI%sDate_beg_sim," -- ",sINI%sDate_end_sim
            write(*,'(A,I5)')"nMon= ", nMonsbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)

            if(sINI%Carbon_only) then
                write(*,*)">>>MEND is RUNNING, PLEASE BE PATIENT: Carbon-Only"
            else                                                                         
                write(*,*)">>>MEND is RUNNING, PLEASE BE PATIENT: C-N-P Coupled"
            end if
          
            sINI % iModel = 0 !RUN MEND MODEL USING (i)PARAMETER VALUES IN LAST LINE OF "MEND_namelist.nml" or (ii)"BEST" PARAMETER VALUES FROM OPTIMIZATION
            fObj = fMEND_OBJ(xx, sPAR, sINI, sOUT)    

            if (sINI%iError < 0) then
                write(*,*)"Mass balance error: sINI%iError = ", sINI%iError
                stop
            end if

            !!write SIM vs. OBS           
            Print*,">>>Write output for SIM vs. OBS: "
            Print*,"ATTENTION: It does NOT matter if there are void rOBJ values due to NO OBS data available."
            do i=1,sINI%nOBS_tot     
                write(sINI%iFout_SIM_obs,'(2i5,i15,3e20.6)')i,int(sINI%dOBS_opt(i,3)),int(sINI%dOBS_opt(i,1)),&
                    sINI%dOBS_opt(i,2),sINI%dSIM_opt(i,2),sINI%dSIM_opt(i,3)   !!obs_mean,sim_mean,sim_sd
            end do

            write(sINI%iFout_SIM_obs,'(/,a)')"OBJ-FUNCTIONS & PARAMETERS:"    
            write(sINI%iFout_SIM_obs,'(a20,f10.4)')"fOBJ (best=0) = ", fObj   
            write(sINI%iFout_SIM_obs,'(a20,4x,20a10)')"Observations: ",sINI%Name_VARopt(1:sINI%nVARopt) 
            write(sINI%iFout_SIM_obs,'(a20,20f10.4)')"fOBJ[i] = ",sINI%rOBJ   
            write(sINI%iFout_SIM_obs,'(a20,20f10.4)')"fOBJ_weight[i] = ",sINI%rOBJw 

            write(format510,*)"(/,",sSCE%nPar,"(a16))"
            write(sINI%iFout_SIM_obs,format510)sSCE%parName

            write(format521,*)"(",sSCE%nPar,"(f15.8,','))"
            !        write(format521,*)"(",sSCE%nPar,"(f30.20,','))"
            write(sINI%iFout_SIM_obs,format521) xx

            !!close all output files
            CALL subMEND_Files_Close(sINI)                    

            !!Write to Screen:
            write(*,'(a20,f10.4)')"fOBJ (best=0) = ", fObj         
            write(*,'(a20,4x,20a10)')"Observations: ",sINI%Name_VARopt(1:sINI%nVARopt)  
            write(*,'(a20,20f10.4)')"fOBJ[i] = ",sINI%rOBJ          
            write(*,'(a20,20f10.4)')"fOBJ_weight[i] = ",sINI%rOBJw  
            
            !!Convert OUTPUTS from HOURLY to DAILY & MONTHLY for all STATE VARIABLEs | FLUXes | RATEs
            CALL sOUT_tscale(sINI,sINI%sDate_beg_sim,sINI%sDate_end_sim)    
        end if !!if (pid.eq.0) then
        
    END IF !!IF(sINI % iModel.lt.2)
    
    DEALLOCATE(xx)  
    DEALLOCATE(sINI%STP)
    DEALLOCATE(sINI%SWC)
    DEALLOCATE(sINI%SWP)
    DEALLOCATE(sINI%SpH)
    DEALLOCATE(sINI%SIN)
    
    IF(.NOT.sINI%Carbon_only) THEN
        DEALLOCATE(sINI%SIN_NH4)
        DEALLOCATE(sINI%SIN_NO3)
    END IF
    
    if (pid.eq.0) then
        CALL system_clock(t_end) !!timer(t_end)
        t_elapse = (t_end - t_start)/(t_rate)
        CALL Sec2HMS(t_elapse,tHMS)
        CALL DATE_AND_TIME(values = time_end) ! Get the current time
       
        write(*,'(a16,3(I3,a8))')">>Elapsed Time = ",tHMS(1),"Hours", tHMS(2),"Minutes",tHMS(3),"Seconds"
        write (*,'(a25,I4,a1,I2,a1,I2,I3,a1,I2,a1,I2)') ">>MEND RUN COMPLETED at: ", &
            time_end(1),'-',time_end(2),'-',time_end(3),time_end(5),':',time_end(6),':',time_end(7)

        write(sSCE%iFout_ini,'(a16,3(I3,a8))')">>Elapsed Time = ",tHMS(1),"Hours", tHMS(2),"Minutes",tHMS(3),"Seconds"
        write (sSCE%iFout_ini,'(a25,I4,a1,I2,a1,I2,I3,a1,I2,a1,I2)') ">>MEND RUN COMPLETED at: ", &
            time_end(1),'-',time_end(2),'-',time_end(3),time_end(5),':',time_end(6),':',time_end(7)
    end if !!if (pid.eq.0)


#ifdef USE_MPI
    CALL MPI_FILE_CLOSE(sSCE%iFout_all, ierr)
    CALL MPI_TYPE_FREE(sSCE%mpi_type_record, ierr)
#else
    close(sSCE%iFout_all)
#endif

    close(sSCE%iFout_end)
    close(sSCE%iFout_ini)

#ifdef USE_MPI
    call MPI_FINALIZE (ierr)
#endif

END PROGRAM MEND_main
!!==================================================================================================================================    

