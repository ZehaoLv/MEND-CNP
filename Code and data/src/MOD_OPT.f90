MODULE MOD_OPT
    ! File:   MOD_OPT.F90
    ! Author: GANGSHENG WANG @ ORNL & ZEHAO LV @ WHU
    ! Updated: May 5, 2015
    ! Updated: Updated: November 03, 2023
#ifdef USE_MPI
    USE mpi
#endif
    USE MOD_OPT_TYPE
    USE MOD_MEND
    USE MOD_USRFS, ONLY: sort,sort1,indexx,gasdev,selectINT

    IMPLICIT NONE
    PRIVATE :: cce
    PRIVATE :: getpnt
    PRIVATE :: parstt
    PRIVATE :: comp
    PRIVATE :: chkcst
    
    PUBLIC :: SCE
    PUBLIC :: sWRITE_OPT_out_head

CONTAINS
    
    
    SUBROUTINE SCE(sSCE, sPAR, sINI, sOUT)

        !    (a,bl,bu,nopt,maxn,kstop,pcento,iseed,ngs,npg,nps,nspl,mings,iniflg,iprint)
        !c$debug
        !c
        !c  MODIFIED by GANGSHENG WANG
        !C  Environmental Sciences Division
        !C  Oak Ridge National Laboratory
        !C  Oak Ridge, TN 37831-6301
        !C  March, 2013
        !C  Updated: April 28, 2015

        !c  SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
        !c     -- Version 2.1
        !c
        !c  by QINGYUN DUAN
        !c  DEPARTMENT OF HYDROLOGY & WATER RESOURCES
        !c  UNIVERSITY OF ARIZONA, TUCSON, AZ 85721
        !c  (602) 621-9360, email: duan@hwr.arizona.edu
        !c
        !c  WRITTEN IN OCTOBER 1990.
        !c  REVISED IN AUGUST 1991
        !c  REVISED IN APRIL 1992
        !c
        !c  STATEMENT BY AUTHOR:
        !c  --------------------
        !c
        !c     This general purpose global optimization program is developed at
        !c     the Department of Hydrology & Water Resources of the University
        !c     of Arizona.  Further information regarding the SCE-UA method can
        !c     be obtained from Dr. Q. Duan, Dr. S. Sorooshian or Dr. V.K. Gupta
        !c     at the address and phone number listed above.  We request all
        !c     users of this program make proper reference to the paper entitled
        !c     'Effective and Efficient Global Optimization for Conceptual
        !c     Rainfall-runoff Models' by Duan, Q., S. Sorooshian, and V.K. Gupta,
        !c     Water Resources Research, Vol 28(4), pp.1015-1031, 1992.
        !c
        !c
        !c  LIST OF INPUT ARGUEMENT VARIABLES
        !c
        !c     a(.) = initial parameter set
        !c     bl(.) = lower bound on parameters
        !c     bu(.) = upper bound on parameters
        !c     nopt = number of parameters to be optimized
        !c
        !c
        !c  LIST OF SCE ALGORITHMIC CONTROL PARAMETERS:
        !c
        !c     ngs = number of complexes in the initial population
        !c     npg = number of points in each complex
        !c     npt = total number of points in initial population (npt=ngs*npg)
        !c     nps = number of points in a sub-complex
        !c     nspl = number of evolution steps allowed for each complex before
        !c         complex shuffling
        !c     mings = minimum number of complexes required, if the number of
        !c         complexes is allowed to reduce as the optimization proceeds
        !c     iseed = initial random seed
        !c     iniflg = flag on whether to include the initial point in population
        !c         = 0, not included
        !c         = 1, included
        !c     iprint = flag for controlling print-out after each shuffling loop
        !c         = 0, print information on the best point of the population
        !c         = 1, print information on every point of the population
        !c
        !c
        !c  CONVERGENCE CHECK PARAMETERS
        !c
        !c     maxn = max no. of trials allowed before optimization is terminated
        !c     kstop = number of shuffling loops in which the criterion value must
        !c         chang by the given percentage before optimization is terminated
        !c     pcento = percentage by which the criterion value must change in
        !c         given number of shuffling loops
        !c     ipcnvg = flag indicating whether parameter convergence is reached
        !c         (i.e., check if gnrng is less than 0.001)
        !c         = 0, parameter convergence not satisfied
        !c         = 1, parameter convergence satisfied
        !c
        !c
        !c  LIST OF LOCAL VARIABLES
        !c     x(.,.) = coordinates of points in the population
        !c     xf(.) = function values of x(.,.)
        !c     xx(.) = coordinates of a single point in x
        !c     cx(.,.) = coordinates of points in a complex
        !c     cf(.) = function values of cx(.,.)
        !c     s(.,.) = coordinates of points in the current simplex
        !c     sf(.) = function values of s(.,.)
        !c     bestx(.) = best point at current shuffling loop
        !c     bestf = function value of bestx(.)
        !c     worstx(.) = worst point at current shuffling loop
        !c     worstf = function value of worstx(.)
        !c     xnstd(.) = standard deviation of parameters in the population
        !c     gnrng = normalized geometric mean of parameter ranges
        !c     lcs(.) = indices locating position of s(.,.) in x(.,.)
        !c     bound(.) = bound on ith variable being optimized
        !c     ngs1 = number of complexes in current population
        !c     ngs2 = number of complexes in last population
        !c     iseed1 = current random seed
        !c     criter(.) = vector containing the best criterion values of the last
        !c         10 shuffling loops
        !c
        !      implicit none!REAL*8 (a-h,o-z)

        !c  ARRAYS FROM THE INPUT DATA
        !!ARGUMENTS:
        TYPE(sSCE_PAR) ,intent(inout):: sSCE
        TYPE(sMEND_PAR),intent(inout):: sPAR
        TYPE(sMEND_INI),intent(inout):: sINI
        TYPE(sMEND_OUT),intent(inout):: sOUT

        !!LOCAL VARIABLES:
        INTEGER i, j, k, k1, k2, l, iCALL, ibound
        INTEGER igs, ngs1, ngs2, nloop, loop, ipcnvg, iseed1
        INTEGER nPar, nOpt, maxn, kstop, iseed
        INTEGER ngs, npg, nps, nspl, mings, npt, npt1, iniflg, iprint
        REAL(KIND=8) fobj, fa, bestf, worstf
        REAL(KIND=8) denomi, timeou, gnrng, pcento
      
        INTEGER iOpt(sSCE%nOpt)
        REAL(KIND=8) a(sSCE%nPar),bl(sSCE%nPar),bu(sSCE%nPar)
        REAL(KIND=8) xx(sSCE%nPar),bestx(sSCE%nPar),worstx(sSCE%nPar)
        REAL(KIND=8) xnstd(sSCE%nPar),bound(sSCE%nPar),criter(20),unit(sSCE%nPar)
      
        REAL(KIND=8) x(sSCE%npt,sSCE%nPar),xf(sSCE%npt)
        REAL(KIND=8) cx(sSCE%npg,sSCE%nPar),cf(sSCE%npg)
        REAL(KIND=8) s(sSCE%nps,sSCE%nPar),sf(sSCE%nps)
        INTEGER lcs(sSCE%nps) !!lcs(50)
      
        CHARACTER*256 format510,format520,format521, format610,format630,format650,format660,format001
        CHARACTER(len=sSCE%OPT_all_line_length + 1):: sline

#ifdef USE_MPI
        INTEGER ierr
#endif
        !      write (*,*) '>>ENTER SUBROUTINE <SCE>'
  
        a     = sSCE%a
        bl    = sSCE%bl
        bu    = sSCE%bu
        nPar  = sSCE%nPar
        nOpt  = sSCE%nOpt
        iOpt  = sSCE%iOpt
        maxn  = sSCE%maxn
        kstop = sSCE%kstop
        iseed = sSCE%iseed
        ngs   = sSCE%ngs
        npg   = sSCE%npg
        npt   = sSCE%npt
        nps   = sSCE%nps
        nspl  = sSCE%nspl
        pcento = sSCE%pcento
        mings  = sSCE%mings
        iniflg = sSCE%iniflg
        iprint = sSCE%iprint
      
      
      
        !      print*, sSCE%parName
            
        write(format510,*)"(/,' CRITERION',",nPar,"(a15),/1x,",&  !(6x,a4)
            &                   (nPar+1)*10,"(1h-))"
        write(format520,*)"(f10.4,",nPar,"f15.8)"

        !      write(format521,*)"(",nPar,"f15.8,","' | '",",f10.4,",sINI%nVARopt,"f10.4)"
        write(format521,*)"(",nPar,"(f15.8,','),","' | '",",f10.4,",sINI%nVARopt,"f10.4)"
        write(format610,*)"(/,1x,'LOOP',1x,'TRIALS',1x,'COMPLXS',2x,",&
            &       "'BESTF',3x,'WORSTF',3x,'PAR-RNG',1x,",nPar,"(a15))"
        write(format630,*)"(i5,1x,i5,3x,i5,3g10.3,",nPar,"(f15.8))"
        !      write(format660,*)"(15x,g10.3,20x,",nPar,"(f15.8))"

        !      write(format660,*)"(g10.3,",nPar,"(f15.8))"
        !      write(format650,*)"(/a10,",nPar,"(a15))"
        write(format660,*)"(I10,g10.3,",nPar,"(f15.8))"
        write(format650,*)"(2a10,",nPar,"(a15))"

        !c  INITIALIZE VARIABLES
      
        nloop = 0
        loop = 0
        igs = 0

        !c  INITIALIZE RANDOM SEED TO A NEGATIVE INTEGER

        iseed1 = iseed !-abs(iseed)

        !c  COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPUALTION
        !      npt = ngs * npg
        ngs1 = ngs
        npt1 = npt

        !      write(sSCE%iFout_ini,400)
        !!------------------------------------------------------------------------
        write (*,'(3(A,I5))') '*  Process = ',sINI%pid,' : Optimization Run = ', &
            sSCE%iRun, ' : Evolution Loop = ',nloop
      
        write(sline, '(A,I5,A,I5,A,I20)')'Process = ',sINI%pid,'; SCE Run Number = ', &
            sSCE%iRun, '; Random Seed = ', sSCE%iseed
        sline((sSCE%OPT_all_line_length+1): (sSCE%OPT_all_line_length+1)) = lineEnd
#ifdef USE_MPI
        CALL MPI_File_write_shared(sSCE%iFout_all, sline, 1, sSCE%mpi_type_record, MPI_STATUS_IGNORE, ierr)
#else
        write(sSCE%iFout_all,*)sline(1:(len(trim(sline))-1))
#endif
        !!------------------------------------------------------------------------

        !c  COMPUTE THE BOUND FOR PARAMETERS BEING OPTIMIZED
        do j = 1, nPar
            bound(j) = bu(j) - bl(j)
            unit(j) = 1.0
            xx(j) = a(j)              !wgs
        end do


        !!------------------------------------------------------------------------
        !c Check if the Initial Value of Parameter is within the defined range
        l = 0
        do j = 1, sSCE%nPar
            CALL chkcst(1,a(j),bl(j),bu(j),ibound)
            if (ibound .ge. 1) then
                l = l +1
                write(*,*)"Constrains are Violated for Initial Value of Parameter: ",j,sSCE%parName(j),a(j), bl(j), bu(j)
            end if
        end do
        if (l.ge.1) STOP
        !!------------------------------------------------------------------------
        !c  COMPUTE THE FUNCTION VALUE OF THE INITIAL POINT
        fa = fMEND_OBJ(a, sPAR, sINI, sOUT)
        if (sINI%iError < 0) return

        !c  PRINT THE INITIAL POINT AND ITS CRITERION VALUE
        !wgs:10/11/2019
        !      write(sSCE%iFout_ini,500)
        !      write(sSCE%iFout_ini,format510) sSCE%parName!(xname(j),j=1,nPar)
        !      write(sSCE%iFout_ini,format520) fa,a !(a(j),j=1,nPar)
      
        !c  GENERATE AN INITIAL SET OF npt1 POINTS IN THE PARAMETER SPACE
        !c  IF iniflg IS EQUAL TO 1, SET x(1,.) TO INITIAL POINT a(.)
      
        if (iniflg .eq. 1) then
            x(1,1:nPAR) = a(1:nPAR)
            xf(1) = fa
                                                        
        !c  ELSE, GENERATE A POINT RANDOMLY AND SET IT EQUAL TO x(1,.)
        else
            CALL getpnt(nPar,nOpt,iOpt,1,xx,bl,bu,unit,bl)          
            x(1,1:nPar) = xx(1:nPar)
            xf(1) = fMEND_OBJ(xx, sPAR, sINI, sOUT)                  
            if (sINI%iError < 0) return
        end if
        iCALL = 1
        if (iCALL .ge. maxn) go to 9000

        !c  GENERATE npt1-1 RANDOM POINTS DISTRIBUTED UNIFORMLY IN THE PARAMETER
        !c  SPACE, AND COMPUTE THE CORRESPONDING FUNCTION VALUES                
        do i = 2, npt1
            CALL getpnt(nPar,nOpt,iOpt,1,xx,bl,bu,unit,bl)           
            x(i,1:nPAR) = xx(1:nPAR)

            xf(i) = fMEND_OBJ(xx, sPAR, sINI, sOUT) !functn(nopt,xx)    
            if (sINI%iError < 0) return
            iCALL = iCALL + 1
            if (iCALL .ge. maxn) then
                npt1 = i
                go to 45
            end if
        end do

        !c  ARRANGE THE POINTS IN ORDER OF INCREASING FUNCTION VALUE
        !      print*, npt1
45      CALL sort(npt1,nPar,x,xf)                             

        !c  RECORD THE BEST AND WORST POINTS
        do j = 1, nPar
            bestx(j) = x(1,j)
            worstx(j) = x(npt1,j)
        end do
        bestf = xf(1)
        worstf = xf(npt1)

        !c  COMPUTE THE PARAMETER RANGE FOR THE INITIAL POPULATION
        CALL parstt(nPar,nOpt,iOpt,npt1,x,xnstd,bound,gnrng,ipcnvg)      

        !c  PRINT THE RESULTS FOR THE INITIAL POPULATION
        !wgs:10/11/2019
        !      write(sSCE%iFout_ini,format610) sSCE%parName
        !      write(sSCE%iFout_ini,format630) nloop,iCALL,ngs1,bestf,worstf,gnrng,&
        !     &               (bestx(j),j=1,nPar)

        if (iprint .eq. 1) then         
            !        write(sSCE%iFout_all,650) nloop
            do i = 1, npt1
                !          write(sSCE%iFout_all,format660) xf(i),x(i,1:nPar)
                !!------------------------------------------------------------------------
                write(sline, format660)sINI%pid,xf(i),x(i,1:nPar)
                sline(len(sline):len(sline)) = lineEnd

#ifdef USE_MPI
                CALL MPI_File_write_shared(sSCE%iFout_all, sline, 1, sSCE%mpi_type_record, MPI_STATUS_IGNORE, ierr)
#else
                write(format001,*)"(a",len(trim(sline))-1,")"
                write(sSCE%iFout_all,format001)sline(1:(len(trim(sline))-1))
#endif
            !!------------------------------------------------------------------------
            end do

        end if

        if (iCALL .ge. maxn) go to 9000
        if (ipcnvg .eq. 1) go to 9200

    !c  BEGIN THE MAIN LOOP ----------------
1000 continue
     nloop = nloop + 1

     write (*,'(3(A,I5))') '*  Process = ',sINI%pid,' : Optimization Run = ', &
         sSCE%iRun, ' : Evolution Loop = ',nloop

     !!================================================================
     !c  BEGIN LOOP ON COMPLEXES                      
     do 3000 igs = 1, ngs1                     

         !c  ASSIGN POINTS INTO COMPLEXES         
         do k1 = 1, npg
             k2 = (k1-1) * ngs1 + igs
             cx(k1,1:nPAR) = x(k2,1:nPAR)
             cf(k1) = xf(k2)
         end do

         !c  BEGIN INNER LOOP - RANDOM SELECTION OF SUB-COMPLEXES ---------------
         do 2000 loop = 1, nspl

             !c  CHOOSE A SUB-COMPLEX (nps points) ACCORDING TO A LINEAR
             !c  PROBABILITY DISTRIBUTION

             CALL selectINT(npg, nps, lcs)              
             !c  ARRANGE THE SUB-COMPLEX IN ORDER OF INCEASING FUNCTION VALUE
             CALL sort1(nps,lcs)

             !c  CREATE THE SUB-COMPLEX ARRAYS
             !   85 do k = 1, nps
             do k = 1, nps
                 s(k,1:nPAR) = cx(lcs(k),1:nPAR)
                 sf(k) = cf(lcs(k))
             end do
             !      write(*,'(/,100I10)')lcs(1:nps)
             !      write(*,'(100F10.4)')sf(1:nps)
             !c  USE THE SUB-COMPLEX TO GENERATE NEW POINT(S)
             CALL cce(sSCE, sPAR, sINI, sOUT, s,sf,xnstd,iCALL)     
             if (sINI%iError < 0) return

             !c  IF THE SUB-COMPLEX IS ACCEPTED, REPLACE THE NEW SUB-COMPLEX
             !c  INTO THE COMPLEX
             do k = 1, nps
                 cx(lcs(k),1:nPAR) = s(k,1:nPAR)
                 cf(lcs(k)) = sf(k)
             end do

             !c  SORT THE POINTS
             CALL sort(npg,nPar,cx,cf)

             !c  IF MAXIMUM NUMBER OF RUNS EXCEEDED, BREAK OUT OF THE LOOP
             if (iCALL .ge. maxn) go to 2222

         !c  END OF INNER LOOP ------------
2000     continue
2222 continue

     !c  REPLACE THE NEW COMPLEX INTO ORIGINAL ARRAY x(.,.)
     do k1 = 1, npg                                            
         k2 = (k1-1) * ngs1 + igs
         x(k2,1:nPAR) = cx(k1,1:nPAR)
         xf(k2) = cf(k1)
     end do
     if (iCALL .ge. maxn) go to 3333

 !c  END LOOP ON COMPLEXES
3000 continue
     !!================================================================

     !c  RE-SORT THE POINTS
     ! print*, npt1
3333 CALL sort(npt1,nPar,x,xf)                     

     !c  RECORD THE BEST AND WORST POINTS
     do j = 1, nPar
         bestx(j) = x(1,j)
         worstx(j) = x(npt1,j)
     end do
     bestf = xf(1)
     worstf = xf(npt1)

     !c  TEST THE POPULATION FOR PARAMETER CONVERGENCE
     CALL parstt(nPar,nOpt,iOpt,npt1,x,xnstd,bound,gnrng,ipcnvg)  
                                                                  
     !c  PRINT THE RESULTS FOR CURRENT POPULATION
     !      if (mod(nloop,5) .ne. 0) go to 501
     !      write(sSCE%iFout_all,format610) sSCE%parName!(xname(j),j=1,nopt)
     !  501 continue
     !wgs:10/11/2019
     !      write(sSCE%iFout_ini,format630) nloop,iCALL,ngs1,bestf,worstf,gnrng,&
     !     &               (bestx(j),j=1,nPar)

     if (iprint .eq. 1) then
         !        write(sSCE%iFout_all,650) nloop
         do i = 1, npt1
             !          write(sSCE%iFout_all,format660) xf(i),x(i,1:nPar)!!(x(i,j),j=1,nPar)
             !!------------------------------------------------------------------------
             write(sline, format660)sINI%pid, xf(i),x(i,1:nPar)
             sline(len(sline):len(sline)) = lineEnd

#ifdef USE_MPI
             CALL MPI_File_write_shared(sSCE%iFout_all, sline, 1, sSCE%mpi_type_record, MPI_STATUS_IGNORE, ierr)
#else
             
             write(format001,*)"(a",len(trim(sline))-1,")"
             write(sSCE%iFout_all,format001)sline(1:(len(trim(sline))-1))
#endif
         !!------------------------------------------------------------------------
         end do
     end if

     !c  TEST IF MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED
     if (iCALL .ge. maxn) go to 9000

     !c  COMPUTE THE COUNT ON SUCCESSIVE LOOPS W/O FUNCTION IMPROVEMENT
     criter(20) = bestf
     if (nloop .lt. (kstop+1)) go to 132
     denomi = dabs(criter(20-kstop) + criter(20)) / 2.
     timeou = dabs(criter(20-kstop) - criter(20)) / denomi
     if (timeou .lt. pcento) go to 9100
132 continue
    do l = 1, 19
        criter(l) = criter(l+1)
    end do

    !c  IF POPULATION IS CONVERGED INTO A SUFFICIENTLY SMALL SPACE
    if (ipcnvg .eq. 1) go to 9200

    !c  NONE OF THE STOPPING CRITERIA IS SATISFIED, CONTINUE SEARCH

    !c  CHECK FOR COMPLEX NUMBER REDUCTION
    if (ngs1 .gt. mings) then
        ngs2 = ngs1
        ngs1 = ngs1 - 1
        npt1 = ngs1 * npg
        !!wgs[7/17/2015]: correction: change 'nOpt' to 'npt1'
        !        CALL comp(nPar,nOpt,ngs1,ngs2,npg,x,xf,cx,cf)
        CALL comp(nPar,npt1,ngs1,ngs2,npg,x,xf,cx,cf)
    ! subroutine comp(n,npt,ngs1,ngs2,npg,a,af,b,bf)
    end if

    !c  END OF MAIN LOOP -----------
    go to 1000

!c  SEARCH TERMINATED                                  
9000 continue
     write(*,800) maxn,loop,igs,nloop
     go to 9999
9100 continue
     write(*,810) pcento*100.,kstop
     go to 9999
9200 write(*,820) gnrng*100.
9999 continue

     !c  PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
     !wgs:10/11/2019
     !      write(sSCE%iFout_ini,830)
     !      write(sSCE%iFout_ini,format510) sSCE%parName
     !      write(sSCE%iFout_ini,format520) bestf,bestx
      
     fObj = fMEND_OBJ(bestx, sPAR, sINI, sOUT)           
     if (sINI%iError < 0) return                          
     !      write(sSCE%iFout_end,format521) bestx,fobj,sINI%rOBJ
     write(*,format521) bestx,fobj,sINI%rOBJ
      
     sSCE%bestObj = bestf
     sSCE%bestPar = bestx
                              
     !      write (*,*) '>>EXIT SUBROUTINE <SCEUA>'
     !c  END OF SUBROUTINE SCEUA
     return
     !  400 format(//,2x,50(1h=),/,2x,'ENTER THE SHUFFLED COMPLEX EVOLUTION',&
     !     &       ' GLOBAL SEARCH',/,2x,50(1h=))
     !  500 format(//,'*** PRINT THE INITIAL POINT AND ITS CRITERION ',&
     !     &       'VALUE ***')
     !  510 format(/,' CRITERION',12(6x,a4),/1x,60(1h-))
     !  520 format(g10.3,12f10.3)
     !  530 format(10x,12(6x,a4))
     !  540 format(10x,12f10.3)
     !  600 format(//,1x,'*** PRINT THE RESULTS OF THE SCE SEARCH ***')
     !  610 format(/,1x,'LOOP',1x,'TRIALS',1x,'COMPLXS',2x,'BEST F',3x,
     !     &       'WORST F',3x,'PAR RNG',1x,8(6x,a4))
     !  620 format(49x,8(6x,a4))
     !  630 format(i5,1x,i5,3x,i5,3g10.3,8(f10.3))
     !  640 format(49x,8(f10.3))
     !  650 format(/,1x,'POPULATION AT LOOP ',i3,/,1x,22(1h-))
     !  660 format(15x,g10.3,20x,8(f10.3))
800  format(//,1x,'*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE',&
         &       ' LIMIT ON THE MAXIMUM',/,5x,'NUMBER OF TRIALS ',i5,&
         &       ' EXCEEDED.  SEARCH WAS STOPPED AT',/,5x,'SUB-COMPLEX ',&
         &       i3,' OF COMPLEX ',i3,' IN SHUFFLING LOOP ',i3,' ***')
810  format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE CRITERION',&
         &       ' VALUE HAS NOT CHANGED ',/,5x,f5.2,' PERCENT IN',i3,&
         &       ' SHUFFLING LOOPS ***')
820  format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE POPULATION',&
         &       ' HAS CONVERGED INTO ',/,4x,f5.2,' PERCENT OF THE',&
         &       ' FEASIBLE SPACE ***')
 !  830 format(//,'*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS',&
 !     &       ' CRITERION VALUE ***')
     
 END !!SUBROUTINE SCE

      
      
 !!------------------------------------------------------------------
 !SUBROUTINE cce(sSCE, sPAR, sINI, sOUT, s,sf,xnstd,iCALL,iseed)
 SUBROUTINE cce(sSCE, sPAR, sINI, sOUT, s,sf,xnstd,iCALL)
     !(nPar,nOpt,iOpt,nps,s,sf,bl,bu,xnstd,iCALL,maxn,iseed)
     !USE MOD_MEND
     !USE MOD_OPT
     !c$debug
     !c
     !c  ALGORITHM GENERATE A NEW POINT(S) FROM A SUB-COMPLEX
     !c
     !c  SUB-COMPLEX VARIABLES
     !      implicit REAL*8 (a-h,o-z)

     !      dimension s(50,16),sf(50),bu(16),bl(16),xnstd(16)

     !c  LIST OF LOCAL VARIABLES
     !c    sb(.) = the best point of the simplex
     !c    sw(.) = the worst point of the simplex
     !c    w2(.) = the second worst point of the simplex
     !c    fw = function value of the worst point
     !c    ce(.) = the centroid of the simplex excluding wo
     !c    snew(.) = new point generated from the simplex
     !c    iviol = flag indicating if constraints are violated
     !c          = 1 , yes
     !c          = 0 , no
     !!ARGUMENTS:
     TYPE(sSCE_PAR) ,intent(inout):: sSCE
     TYPE(sMEND_PAR),intent(inout):: sPAR
     TYPE(sMEND_INI),intent(inout):: sINI
     TYPE(sMEND_OUT),intent(inout):: sOUT

     !      parameter (c1=0.8,c2=0.4)
     INTEGER iCALL !!, iseed
     !      REAL(KIND=8):: s(50,sSCE%nPar),sf(50),xnstd(sSCE%nPar)
     REAL(KIND=8):: s(:,:),sf(:),xnstd(:)
      
     !!LOCAL VARIABLES:
     INTEGER :: iOpt(sSCE%nOpt)
     REAL(KIND=8) :: bu(sSCE%nPar),bl(sSCE%nPar)
     REAL(KIND=8) sw(sSCE%nPar),sb(sSCE%nPar),ce(sSCE%nPar),snew(sSCE%nPar)
     INTEGER nPar,nOpt,nps,maxn,n,m
     REAL(KIND=8) alpha,beta,fw,fnew
     INTEGER i,j,k,ibound

     !c  EQUIVALENCE OF VARIABLES FOR READABILTY OF CODE
     nPar = sSCE%nPar
     nOpt = sSCE%nOpt
     iOpt = sSCE%iOpt
     nps = sSCE%nps
     maxn = sSCE%maxn
     bl = sSCE%bl
     bu = sSCE%bu
      
     n = nps
     m = nOpt
     alpha = 1.0d0
     beta = 0.5d0
      
     do i = 1, nPar
         snew(i) = s(1, i)
         sb(i) = s(1, i)
     end do
     !c  IDENTIFY THE WORST POINT wo OF THE SUB-COMPLEX s
     !c  COMPUTE THE CENTROID ce OF THE REMAINING POINTS
     !c  COMPUTE step, THE VECTOR BETWEEN wo AND ce
     !c  IDENTIFY THE WORST FUNCTION VALUE fw
     do k = 1, m
         j = iOpt(k)
         sb(j) = s(1,j)
         sw(j) = s(n,j)
         ce(j) = 0.0
         do i = 1, n-1
             ce(j) = ce(j) + s(i,j)
         end do
         ce(j) = ce(j)/dble(n-1)
     end do
     fw = sf(n)

     !c  COMPUTE THE NEW POINT snew

     !c  FIRST TRY A REFLECTION STEP
     do k = 1, m
         j = iOpt(k)
         snew(j) = ce(j) + alpha * (ce(j) - sw(j))
     end do

     !c  CHECK IF snew SATISFIES ALL CONSTRAINTS
     CALL chkcst(nPar,snew,bl,bu,ibound)


     !c  snew IS OUTSIDE THE BOUND,
     !c  CHOOSE A POINT AT RANDOM WITHIN FEASIBLE REGION ACCORDING TO
     !c  A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
     !c  AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
     if (ibound .ge. 1) CALL getpnt(nPar,nOpt,iOpt,2,snew,bl,bu,xnstd,sb)

     !c  COMPUTE THE FUNCTION VALUE AT snew
     fnew = fMEND_OBJ(snew, sPAR, sINI, sOUT) !functn(nopt,snew)
     if (sINI%iError < 0) return
     iCALL = iCALL + 1
     !c
     !c  COMPARE fnew WITH THE WORST FUNCTION VALUE fw
     !c
     !c  fnew IS LESS THAN fw, ACCEPT THE NEW POINT snew AND RETURN
     if (fnew .le. fw) go to 2000
     if (iCALL .ge. maxn) go to 3000

     !c  fnew IS GREATER THAN fw, SO TRY A CONTRACTION STEP
     do k = 1, m
         j = iOpt(k)
         snew(j) = ce(j) - beta * (ce(j) - sw(j))
     end do

     !c  COMPUTE THE FUNCTION VALUE OF THE CONTRACTED POINT
      
     fnew = fMEND_OBJ(snew, sPAR, sINI, sOUT) !functn(nopt,snew)
     if (sINI%iError < 0) return
     iCALL = iCALL + 1
     !c  COMPARE fnew TO THE WORST VALUE fw
     !c  IF fnew IS LESS THAN OR EQUAL TO fw, THEN ACCEPT THE POINT AND RETURN
     if (fnew .le. fw) go to 2000
     if (iCALL .ge. maxn) go to 3000


     !c  IF BOTH REFLECTION AND CONTRACTION FAIL, CHOOSE ANOTHER POINT
     !c  ACCORDING TO A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
     !c  AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
     CALL getpnt(nPar,nOpt,iOpt,2,snew,bl,bu,xnstd,sb)

     !c  COMPUTE THE FUNCTION VALUE AT THE RANDOM POINT
     fnew = fMEND_OBJ(snew, sPAR, sINI, sOUT) !functn(nopt,snew)
     if (sINI%iError < 0) return
     iCALL = iCALL + 1


 !c  REPLACE THE WORST POINT BY THE NEW POINT
2000 continue

     s(n,1:nPAR) = snew(1:nPAR)
     sf(n) = fnew
3000 continue

     !c  END OF SUBROUTINE CCE
     return
 END !!SUBROUTINE CCE

 !!------------------------------------------------------------------
 !SUBROUTINE getpnt(nPar,nOpt,iOpt,idist,iseed,x,bl,bu,std,xi)
 SUBROUTINE getpnt(nPar,nOpt,iOpt,idist,x,bl,bu,std,xi)
     !c$debug
     !c
     !c     This subroutine generates a new point within feasible region
     !c
     !c     x(.) = new point
     !c     xi(.) = focal point
     !c     bl(.) = lower bound
     !c     bu(.) = upper bound
     !c     std(.) = standard deviation of probability distribution
     !c     idist = probability flag
     !c           = 1 - uniform distribution
     !c           = 2 - Gaussian distribution
    
     !!ARGUMENTS:
     INTEGER nPar,nOpt,idist !!,iseed
     INTEGER iOpt(nOpt)
     REAL(KIND=8) x(nPar),bl(nPar),bu(nPar),std(nPar),xi(nPar)
      
     !!LOCAL VARIABLES:
     INTEGER k,j,ibound
     REAL(KIND=8) randx

     1 do k=1, nOpt
         j = iOpt(k)
2        if (idist .eq. 1)call random_number(randx)!randx = rand() !ran1(iseed)
         if (idist .eq. 2) then
             !            randx = gasdev(iseed)
             randx = gasdev()
         !            print*,">>>Gaussian=",randx,iseed
         end if
         x(j) = xi(j) + std(j) * randx * (bu(j) - bl(j))

         !c     Check explicit constraints
        
         CALL chkcst(1,x(j),bl(j),bu(j),ibound)
         if (ibound .ge. 1) go to 2
     end do

     !c     Check implicit constraints
      
     CALL chkcst(nPar,x,bl,bu,ibound)
     if (ibound .ge. 1) go to 1

     return
 END

 !!------------------------------------------------------------------
 SUBROUTINE parstt(nPar,nOpt,iOpt,npt,x,xnstd,bound,gnrng,ipcnvg)
     !c$debug
     !c
     !c  SUBROUTINE CHECKING FOR PARAMETER CONVERGENCE
     !!ARGUMENTS:
     INTEGER nPar,nOpt,npt,ipcnvg
     INTEGER:: iOpt(nOpt)
     REAL(KIND=8):: gnrng
     !      REAL(KIND=8):: x(2000,nPar),xnstd(nPar),bound(nPar)
     REAL(KIND=8):: x(:,:),xnstd(nPar),bound(nPar)
      
     !!LOCAL VARIABLES:
     INTEGER i,j,k
     REAL(KIND=8), parameter:: delta = 1.0d-20,peps=1.0d-3
     REAL(KIND=8):: xmean(nPar),xmax(nPar),xmin(nPar)
     REAL(KIND=8) gsum,xsum1,xsum2

     !c  COMPUTE MAXIMUM, MINIMUM AND STANDARD DEVIATION OF PARAMETER VALUES
     gsum = 0.d0
     do j = 1, nOpt
         k = iOpt(j)
         xmax(k) = -1.0d+20
         xmin(k) = 1.0d+20
         xsum1 = 0.d0
         xsum2 = 0.d0
         do i = 1, npt
             xmax(k) = dmax1(x(i,k), xmax(k))
             xmin(k) = dmin1(x(i,k), xmin(k))
             xsum1 = xsum1 + x(i,k)
             xsum2 = xsum2 + x(i,k)*x(i,k)
         end do
         xmean(k) = xsum1 / dble(npt)
         xnstd(k) = (xsum2 / dble(npt) - xmean(k)*xmean(k))
         if (xnstd(k) .le. delta) xnstd(k) = delta
         xnstd(k) = dsqrt(xnstd(k))
         xnstd(k) = xnstd(k) / bound(k)
         gsum = gsum + dlog( delta + (xmax(k)-xmin(k))/bound(k) )
     end do
     gnrng = dexp(gsum/dble(nOpt))

     !c  CHECK IF NORMALIZED STANDARD DEVIATION OF PARAMETER IS <= eps
     ipcnvg = 0
     if (gnrng .le. peps) then
         ipcnvg = 1
     end if

     !c  END OF SUBROUTINE PARSTT
     return
 END

      
      
 !!------------------------------------------------------------------
 SUBROUTINE comp(n,npt,ngs1,ngs2,npg,a,af,b,bf)
     !c$debug
     !c
     !c
     !c  THIS SUBROUTINE REDUCE INPUT MATRIX a(n,ngs2*npg) TO MATRIX
     !c  b(n,ngs1*npg) AND VECTOR af(ngs2*npg) TO VECTOR bf(ngs1*npg)
     !      implicit REAL*8 (a-h,o-z)
     !!ARGUMENTS:
     !      REAL(KIND=8):: a(2000,n),af(2000),b(2000,n),bf(2000)
     REAL(KIND=8):: a(:,:),af(:),b(:,:),bf(:)
     INTEGER n,npt,ngs1,ngs2,npg
    
     !!LOCAL VARIABLES:
     INTEGER igs,ipg,k1,k2,j,i
     !      dimension a(2000,16),af(2000),b(2000,16),bf(2000)
     do igs=1, ngs1
         do ipg=1, npg
             k1=(ipg-1)*ngs2 + igs
             k2=(ipg-1)*ngs1 + igs
             do i=1, n
                 b(k2,i) = a(k1,i)
             end do
             bf(k2) = af(k1)
         end do
     end do

     do j=1, npt
         do i=1, n
             a(j,i) = b(j,i)
         end do
         af(j) = bf(j)
     end do

     !c  END OF SUBROUTINE COMP
     return
 END
      
 !!------------------------------------------------------------------
 SUBROUTINE chkcst(nPar,x,bl,bu,ibound)
     !c
     !c     This subroutine check if the trial point satisfies all
     !c     constraints.
     !c
     !c     ibound - violation indicator
     !c            = -1 initial value
     !c            = 0  no violation
     !c            = 1  violation
     !c     nopt = number of optimizing variables
     !c     ii = the ii'th variable of the arrays x, bl, and bu
     !c
     !!ARGUMENTS:
     INTEGER nPar,ibound
     REAL(KIND=8) x(nPar),bl(nPar),bu(nPar)
      
     !!LOCAL VARIABLES:
     INTEGER ii

     ibound = -1

     !c     Check if explicit constraints are violated

     do ii=1, nPar
         if (x(ii) .lt. bl(ii) .or. x(ii) .gt. bu(ii)) then
!!            write(*,*)"Constrains are violated for Parameter: ",ii,x(ii), bl(ii), bu(ii)
            go to 10
         end if
     end do
     if (nPar .eq. 1) go to 9

     !c     Check if implicit constraints are violated
     !c     (no implicit constraints for this function)
     !c
     !c     No constraints are violated
     !c
9    ibound = 0
     return

     !c     At least one of the constraints are violated
     
10   ibound = 1
   
     !      do k = 1, nPar
     !          print*, k, x(k), bl(k), bu(k)
     !      end do
     !      print*, ibound
    
     return
 END

 !!-----------------------------------------------------------------
 SUBROUTINE sWRITE_OPT_out_head(sSCE)
     !! head for OPT_out.ini
     IMPLICIT NONE
     !!ARGUMENTS:
     TYPE(sSCE_PAR), intent(in):: sSCE
#ifdef USE_MPI
     INTEGER ierr
#endif
     CHARACTER(len=256) format650,format001
     CHARACTER(len=sSCE%OPT_all_line_length + 1):: sline

     write(format650,*)"(2a10,",sSCE%nPar,"(a15))"
     write(sline, format650)'Process','fOBJ',sSCE%parName(1:sSCE%nPar)
     write(sSCE%iFout_all, format650)'Process','fOBJ',sSCE%parName(1:sSCE%nPar)
     sline(len(sline):len(sline)) = lineEnd
#ifdef USE_MPI
     CALL MPI_File_write_shared(sSCE%iFout_all, sline, 1, &
         sSCE%mpi_type_record, MPI_STATUS_IGNORE, ierr)
#else
     !! delete the last character "carriage return" to avoid a blank line after each record
     write(format001,*)"(a",len(trim(sline))-1,")"
     write(sSCE%iFout_all,format001)sline(1:(len(trim(sline))-1))
#endif
     RETURN
 END
 END MODULE MOD_OPT
!!c==============================================================


