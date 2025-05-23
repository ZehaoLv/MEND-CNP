MODULE MOD_MEND_TYPE    
    ! File:   STRUCT_MEND.F90
    ! Author: GANGSHENG WANG @ ORNL & ZEHAO LV @ WHU
    ! Updated: July 15, 2015
    ! Updated: Updated: November 03, 2023
    
    ! ``STRUCTURE /name/ ... END STRUCTURE'' becomes
    ! ``TYPE name ... END TYPE''
    ! ``RECORD /name/ variable'' becomes ``TYPE(name) variable''
    INTEGER, PARAMETER:: const_nPOOL        = 29                    !#### of state variables in sMEND_POOL         
    INTEGER, PARAMETER:: const_nPOOL_MN     = 22                    !#### of state variables in sMEND_POOL_MN      
    INTEGER, PARAMETER:: const_nFLUX        = 47                    !#### of flux variables in sMEND_FLUX          
    INTEGER, PARAMETER:: const_nFLUX_MN     = 56                    !#### of flux variables in sMEND_FLUX_MN       
    INTEGER, PARAMETER:: const_nPAR         = 86                    !#### of parameters in sMEND_PAR               
    INTEGER, PARAMETER:: const_nPAR0        = 75                    !#### of parameters for calibration            
    INTEGER, PARAMETER:: const_nVAR0        = 44                    !# of variables for calibration              
    INTEGER, PARAMETER:: const_nRATE        = 46                    !# of derived rates in '*_RATE_*.out', 19
    
    INTEGER, PARAMETER:: const_nPOM = 2                             !# of particular organic carbon pools
    INTEGER, PARAMETER:: const_nISO = 2                             !# of isotopes, e.g., C12, C14, C13
    INTEGER, PARAMETER:: const_nENZNm = 6                           !# of enzymes for N fixation (1), nitrification (1), denitrification (4)
    INTEGER, PARAMETER:: const_nDenit = 4                           !# of enzymes for denitrification
    
    INTEGER, PARAMETER:: const_nHourInYear = 8784                   !! = 366 days * 24 hours
    
    REAL(KIND=8), PARAMETER:: const_R = 8.314d0                          ![J/mol/K],universal gas constant
    REAL(KIND=8), PARAMETER:: const_tmp_C2K = 273.15d0                   !conversion of degree C to degree K
    REAL(KIND=8), PARAMETER:: const_Rstd(const_nISO - 1) = (/1d-12/)     !C14/C12
    REAL(KIND=8), PARAMETER:: const_cm2MPa = 98d-6                       !1cm water column 
    
    REAL(KIND=8), PARAMETER:: const_CN_Cellulose = 500.d0                !C:N ratio for Cellulose
    REAL(KIND=8), PARAMETER:: const_CP_Cellulose = 1000.d0                !#### C:P ratio for Cellulose
    
    REAL(KIND=8), PARAMETER:: const_Tref         = 20.d0                 ![degree c],reference temperature
    REAL(KIND=8), PARAMETER:: const_SWPmin       = -3.0d4                ![MPa], lowest SWP
    REAL(KIND=8), PARAMETER:: const_SWPFC       = -0.033d0                ![MPa], SWP at field capacity
    REAL(KIND=8), PARAMETER:: const_FillValue    = -999d0                !filled value
    
    REAL(KIND=8), PARAMETER:: const_Pa0 = 101.3d0                         !kPa, atmospheric pressure at sea level
    REAL(KIND=8), PARAMETER:: const_DiffAir0_NOx = 0.05148  ![m2 h-1],NO, N2O, N2 diffusivity at NTP (0 C, 101.3kPa) 
    
    character(len=4), PARAMETER:: const_obj_name(8) = (/"NSEC","NSEn","MARE","MARt", "MARn", "CORR", "CORl", "AVGr"/)
    character, PARAMETER :: cBackspace = char(13)  !!carriage return
    character, PARAMETER :: lineEnd = Achar(10)   !!?? carriage return??
    ! All the "lines" in the file will be this length
    
      
    !-----------------------------------------------------------------------------!
    !Adsorption and Desorption of DOC                                         
    TYPE sSORP_PAR
        REAL(KIND=8) Qmax                        ![mg C/g soil], adsorption capacity
        REAL(KIND=8) Kads                        ![mg C/mg C/h], specific adsorption rate
        REAL(KIND=8) Kdes                        ![mg C/g soil/h], desorption rate
    END TYPE sSORP_PAR
    
    TYPE sSORP_INP
        REAL(KIND=8) adsorbate                   !e.g., DOC
        REAL(KIND=8) adsorbent                   !e.g., QOC (MOC)
    END TYPE sSORP_INP
    
    TYPE sSORP_OUT
        REAL(KIND=8) ads                         ![mg C/g soil/h], adsorption flux
        REAL(KIND=8) des                         ![mg C/g soil/h], desorption flux
        REAL(KIND=8) ads_net                     ![mg C/g soil/h], net adsorption flux = ads - des
    END TYPE sSORP_OUT
    !-----------------------------------------------------------------------------!
    !Michaelis-Menten Kinetics
    TYPE sMM_PAR
        REAL(KIND=8) vm                         ![mg C/mg C/h], specific enzyme activity
        REAL(KIND=8) km                         ![mg C/g soil], half-saturation constant
    END TYPE sMM_PAR
     
    TYPE sMM_INP
        REAL(KIND=8) substrate                  ![mg C/g soil], substrate concentration
        REAL(KIND=8) enzyme                     ![mg C/g soil], enzyme concentration
    END TYPE sMM_INP
    !-----------------------------------------------------------------------------!
    !Diffusivity                                

    
    TYPE sDIFFUSIVITY_PAR
        REAL(KIND=8) Poro                        !porosity
        REAL(KIND=8) SWCFC                       !field capacity
        REAL(KIND=8) DCa                         !gas diffusion coefficient in air
    END TYPE sDIFFUSIVITY_PAR
    
    TYPE sDIFFUSIVITY_INP
        REAL(KIND=8) SWC                         !volumetric soil water content
    END TYPE sDIFFUSIVITY_INP
    
    TYPE sDIFFUSIVITY_OUT
        REAL(KIND=8) DCratio                     !=DCs/DCa
        REAL(KIND=8) DCs                         !gas diffusion coefficient in soil 
        REAL(KIND=8) AFP                         !air-filled porosity
        REAL(KIND=8) PoroIn                      !internal porosity
        REAL(KIND=8) PoroEx                      !external porosity
        REAL(KIND=8) AFPIn                       !internal air-filled porosity
        REAL(KIND=8) AFPEx                       !/external air-filled porosity
        REAL(KIND=8) x(3)
        REAL(KIND=8) y(3)
    END TYPE sDIFFUSIVITY_OUT
    
    TYPE sDIFFUSIVITY
        TYPE(sDIFFUSIVITY_PAR) PAR
        TYPE(sDIFFUSIVITY_INP) INP
        TYPE(sDIFFUSIVITY_OUT) OUT
    END TYPE sDIFFUSIVITY
    !-----------------------------------------------------------------------------!
    !=============================================================================!
    !Microbial-Enzyme-mediated Nitrification-Denitrification-Decomposition (MEND)
    !-----------------------------------------------------------------------------!
    TYPE sMEND_POOL
        REAL(KIND=8) TM                         ![mg C/g soil], TOTAL MASS, CPOOL%TM = TOM, NPOOL%TM = TOM + Nmine
        REAL(KIND=8) TOM                        ![mg C/g soil], total OM = SOM + DOM + MB + ENZ
        REAL(KIND=8) SOM                        ![mg C/g soil], total SOM = POM + MOM + QOM
        REAL(KIND=8) POM(const_nPOM)            ![mg C/g soil],Particulate Organic Carbon, size of {} determined by nPOM
        REAL(KIND=8) POMT                        ![mg C/g soil], POMT = sum(POM)
        REAL(KIND=8) MOMT                        ![mg C/g soil], MOMT = MOM + QOM
        REAL(KIND=8) MOM                        ![mg C/g soil],Mineral Associate Organic Carbon
        REAL(KIND=8) QOM                        ![mg C/g soil],adsorbed phase of DOC
        REAL(KIND=8) DOM                        ![mg C/g soil],Liquid Dissolved Organic matter
        REAL(KIND=8) DOMS                       ![mg C/g soil],Solid Dissolved Organic matter
        REAL(KIND=8) DOMT                       ![mg C/g soil],Total Dissolved Organic matter, TDOM = SDOM + LDOM
        REAL(KIND=8) MB                         ![mg C/g soil],Microbial Biomass Carbon
        REAL(KIND=8) MBA                        ![mg C/g soil],Active Microbial Biomass Carbon
        REAL(KIND=8) MBD                        ![mg C/g soil],Dormant Microbial Biomass Carbon
        REAL(KIND=8) ENZ                        ![mg C/g soil],TOTAL ENZYME BIOMASS = ENZD + ENZNmt
        REAL(KIND=8) ENZD                       ![mg C/g soil],TOTAL ENZYME BIOMASS = SUM(ENZP) + ENZM
        REAL(KIND=8) ENZP(const_nPOM)           ![mg C/g soil],ENZyme for POM
        REAL(KIND=8) ENZM                       ![mg C/g soil],Enzyme for MAOC
        REAL(KIND=8) ENZNmt                     ![mg C/g soil],Total ENZymes for mineral-N = sum(ENZNmine(1:const_nENZNmine) )
        REAL(KIND=8) ENZNm(const_nENZNm)        ![mg C/g soil],ENZymes for mineral-N
        REAL(KIND=8) PTASE                      !####  [mg C/g soil],ENZymes for mineral-P  
        REAL(KIND=8) TM_err                     ![-],default = 0, mass balance check error = (POOL_end - POOL_beg)-(TOTinp - TOTout)
    END TYPE sMEND_POOL
    !-----------------------------------------------------------------------------!
    TYPE sMEND_POOL_MN
        REAL(KIND=8) CO2                        ![mg C/g soil],cumulative CO2 in the soil
        REAL(KIND=8) Nmine                      ![mg N/g soil],
        REAL(KIND=8) Nmine_Free                 ![mg N/g soil], all mineral N except NH4ads
        REAL(KIND=8) Nmine_Solid                ![mg N/g soil], NH4tot + NO3 + NO2
        REAL(KIND=8) NH4tot                     ![mg N/g soil], total NH4
        REAL(KIND=8) NH4ads                     ![mg N/g soil], adsorbed NH4
        REAL(KIND=8) NH4                        ![mg N/g soil], free NH4
        REAL(KIND=8) NOx                        ![mg N/g soil], NOx = NO3 + NO2 + NO + N2O
        REAL(KIND=8) NO32                        ![mg N/g soil], NO23 = NO3 + NO2
        REAL(KIND=8) NO3                        ![mg N/g soil]
        REAL(KIND=8) NO2                        ![mg N/g soil]
        REAL(KIND=8) NO                         ![mg N/g soil]
        REAL(KIND=8) N2O                        ![mg N/g soil]
        REAL(KIND=8) N2                         ![mg N/g soil]
        REAL(KIND=8) NGas                       ![mg N/g soil], NGas=NO + N2O + N2
        REAL(KIND=8) DIP                        !#### [mg P/g soil]
        REAL(KIND=8) LIP                        !#### [mg P/g soil]
        REAL(KIND=8) SP                         !#### [mg P/g soil], SP=DIP + LIP
        REAL(KIND=8) PIP                         !#### [mg P/g soil]
        REAL(KIND=8) SIP                        !#### [mg P/g soil]
        REAL(KIND=8) OIP                        !#### [mg P/g soil]
        REAL(KIND=8) P_dep                      !#### 
    END TYPE sMEND_POOL_MN
    !-----------------------------------------------------------------------------!
    TYPE sMEND_FLUX
        REAL(KIND=8) TOTout                     ![mg POM/g soil/h], total outputs
        REAL(KIND=8) TOTinp                     ![mg POM/g soil/h], total inputs
        REAL(KIND=8) POMadd(const_nPOM)         ![mg POM/g soil/h],inputs to POM
        REAL(KIND=8) DOMadd                     ![mg DOC/g soil/h], inputs to DOC
         
        REAL(KIND=8) POMdec(const_nPOM)         ![mg C/g soil/h],decomposition of POM
        REAL(KIND=8) POMdec_to_DOM(const_nPOM)  ![mg C/g soil/h],decomposition of POM allocated to DOC
        REAL(KIND=8) POMdec_to_MOM(const_nPOM)  ![mg C/g soil/h],decomposition of POM allocated to MOC
        REAL(KIND=8) MOMdec                     ![mg C/g soil/h],decomposition of MOM
        REAL(KIND=8) MOM_to_DOM                 ![mg C/g soil/h],decomposition of MOC
        REAL(KIND=8) QOM_to_DOM                 ![mg C/g soil/h], desorption flux
        REAL(KIND=8) DOM_to_QOM                 ![mg C/g soil/h], adsorption flux
        REAL(KIND=8) DOM_to_QOM_net             ![mg C/g soil/h], net adsorption flux = DOC_to_QOC - QOC_to_DOC
        REAL(KIND=8) QOM_to_MOM                 ![mg C/g soil/h], unprotected to protected MOM
        REAL(KIND=8) MOM_to_QOM                 ![mg C/g soil/h], protected to unprotected
        REAL(KIND=8) MOM_to_QOM_net             ![mg C/g soil/h], net protected to unprotected = MOM_to_QOM - QOM_to_MOM
        REAL(KIND=8) DOM_to_MBA                 ![mg C/g soil/h], uptake of DOC by MBC
        REAL(KIND=8) MBA_mortality              ![mg C/g soil/h], microbial mortality
        REAL(KIND=8) MBA_to_ENZP(const_nPOM)    ![mg C/g soil/h], enzyme-POM production
        REAL(KIND=8) MBA_to_ENZM                ![mg C/g soil/h], enzyme-MOC production
        REAL(KIND=8) MBA_PM                     ![mg C/g soil/h], microbial mortality + enzyme production
        REAL(KIND=8) ENZP_to_DOM(const_nPOM)    ![mg C/g soil/h], enzyme-POM turnover
        REAL(KIND=8) ENZM_to_DOM                ![mg C/g soil/h], enzyme-MOC turnover
        REAL(KIND=8) MBA_to_DOM                 ![mg C/g soil/h], turnover of MBC allocated to DOC
        REAL(KIND=8) MBA_to_POM(const_nPOM)     ![mg C/g soil/h], turnover of MBC allocated to POM (POM_O and POM_H)
        REAL(KIND=8) MBA_to_MBD                 ![mg C/g soil/h], dormancy
        REAL(KIND=8) MBD_to_MBA                 ![mg C/g soil/h], reactivation
        REAL(KIND=8) MBA_to_ENZNm(const_nENZNm)    ![mg C/g soil/h], N-enzyme production
        REAL(KIND=8) ENZNm_to_DOM(const_nENZNm)    ![mg C/g soil/h], N-enzyme production
        REAL(KIND=8) MBA_to_PTASE                 !#### [mg C/g soil/h]
        REAL(KIND=8) PTASE_to_DOM                 !#### [mg C/g soil/h]
    END TYPE sMEND_FLUX
    !-----------------------------------------------------------------------------!
    TYPE sMEND_FLUX_MN
        REAL(KIND=8) CO2_soil                   ![mg C/g soil/h], soil respiration rate = CO2_root + CO2_gmo
        REAL(KIND=8) CO2_root                   ![mg C/g soil/h], root/autotrophic respiration rate
        REAL(KIND=8) CO2_gmo                    ![mg C/g soil/h], growth + maintenance + overflow respiration rate
        REAL(KIND=8) CO2_gm                     ![mg C/g soil/h], growth + maintenance respiration rate
        REAL(KIND=8) CO2_growth                 ![mg C/g soil/h], growth respiration rate
        REAL(KIND=8) CO2_maintn                 ![mg C/g soil/h], maintenance respiration rate
        REAL(KIND=8) CO2_ovflow_N                 ![mg C/g soil/h], overflow respiration rate     
        REAL(KIND=8) CO2_ovflow_P                 ![mg C/g soil/h], overflow respiration rate     
        REAL(KIND=8) CO2_maintn_MBA             ![mg C/g soil/h], maintenance respiration rate of active microbes
        REAL(KIND=8) CO2_maintn_MBD             ![mg C/g soil/h], maintenance respiration rate of dormant microbes
        REAL(KIND=8) CO2_ovflow_MBA_N             ![mg C/g soil/h], overflow respiration rate of active microbes  
        REAL(KIND=8) CO2_ovflow_MBD_N             ![mg C/g soil/h], overflow respiration rate of dormant microbes 
        REAL(KIND=8) CO2_ovflow_MBA_P             ![mg C/g soil/h], overflow respiration rate of active microbes  
        REAL(KIND=8) CO2_ovflow_MBD_P             ![mg C/g soil/h], overflow respiration rate of dormant microbes 
        
        REAL(KIND=8) Nmine_dep                  ![mg N/g soil/h], total mineral N deposition rate = NH4_dep + NO3_dep
        REAL(KIND=8) NH4_dep                    ![mg N/g soil/h], NH4 deposition rate
        REAL(KIND=8) NO3_dep                    ![mg N/g soil/h], NO3 deposition rate
        REAL(KIND=8) P_dep                      !####
         
        REAL(KIND=8) Nmn_net                    !NET N mineralization: = Nmn - Nim
        REAL(KIND=8) Nmn                        !N mineralization: MB = Nmn_MBA + Nmn_MBD
        REAL(KIND=8) Nmn_MBA                    !N mineralization: MBA
        REAL(KIND=8) Nmn_MBD                    !N mineralization: MBD
        REAL(KIND=8) Nim                        !N immobilization = Nim_NH4 + Nim_NO3: MBA
        REAL(KIND=8) Nim_NH4                    !NH4 immobilization: MBA
        REAL(KIND=8) Nim_NO3                    !NO3 immobilization: MBA
        
        REAL(KIND=8) Pmn_net                    !#### NET P mineralization: = Pmn - Pim  [mg P/g soil/h]
        REAL(KIND=8) Pmn                        !#### P mineralization: MB = Pmn_MBA + Pmn_MBD [mg P/g soil/h]
        REAL(KIND=8) Pmn_MBA                    !#### P mineralization: MBA  [mg P/g soil/h]
        REAL(KIND=8) Pmn_MBD                    !#### P mineralization: MBD  [mg P/g soil/h]
        REAL(KIND=8) Pim                        !#### P immobilization       [mg P/g soil/h]
        
        REAL(KIND=8) NFix                       !N fixation flux
        REAL(KIND=8) Nitrif                     !nitrification flux
		REAL(KIND=8) Nitrifnet                     !net nitrification flux
        REAL(KIND=8) Denit(4)                !denitrification fluxes
        REAL(KIND=8) NO_efflux                  !efflux from soil to air
        REAL(KIND=8) N2O_efflux
        REAL(KIND=8) N2_efflux
        
        REAL(KIND=8) Nim_VG                        !N immobilization = Nim_NH4 + Nim_NO3 by VG: Vegetation/Plant
        REAL(KIND=8) Nim_NH4_VG                    !NH4 immobilization: MBA
        REAL(KIND=8) Nim_NO3_VG                    !NO3 immobilization: MBA
        
        REAL(KIND=8) Pim_VG                        !#### P immobilization 
        REAL(KIND=8) Pim_LIP                       !#### P immobilization 
        
        REAL(KIND=8) NO32_Leaching               !NO3 + NO2 leaching
        REAL(KIND=8) NO3_Leaching
        REAL(KIND=8) NO2_Leaching
        REAL(KIND=8) PO4_Leaching
        
        REAL(KIND=8) POMdec_to_LIP(2)                 !#### 
        REAL(KIND=8) MOMdec_to_LIP                    !#### 
        
        REAL(KIND=8) weath_flux                    !####   inorganic P CONVERTION 
        REAL(KIND=8) occlu_flux                    !####   
        REAL(KIND=8) adsorp_flux                   !####   
        REAL(KIND=8) desorp_flux                   !####   
    END TYPE sMEND_FLUX_MN
    !-----------------------------------------------------------------------------!
    !     TYPE sMEND_ADD
    !         REAL(KIND=8) POMadd(const_nPOM)         ![mg POM/g soil/h],inputs to POM
    !         REAL(KIND=8) DOMadd                     ![mg DOC/g soil/h], inputs to DOC
    !         REAL(KIND=8) NH4_dep                    ![mg N/g soil/h], NH4 deposition rate
    !         REAL(KIND=8) NO3_dep                    ![mg N/g soil/h], NO3 deposition rate
    !     END TYPE sMEND_ADD
    !-----------------------------------------------------------------------------!
    TYPE sMEND_PAR
        !         INTEGER nPar
        !         INTEGER iKinetics                 !decomposition kinetics: 0-Michaelis-Menten, 1-First Order, 2-Second Order
        REAL(KIND=8) fRa                        !Root/Autotrophic respiration CO2_root = fRa*GPP
        REAL(KIND=8) fINP                       !actual litter input = fINP*SIN
        REAL(KIND=8) VdPOM(const_nPOM)          ![mg POM/mg ENZP/h], maximum reaction rate for conversion of POM by ENZP
        REAL(KIND=8) VdMOM                      ![mg MOC/mg ENZMAOC/h], maximum reaction rate for conversion of MAOC by ENZMAOC
        REAL(KIND=8) KsPOM(const_nPOM)          ![mg POM/cm3], half-saturation constant for conversion of POM by ENZP
        REAL(KIND=8) KsMOM                      ![mg MOC/cm3], half-saturation constant for conversion of MAOC by ENZMAOC
        REAL(KIND=8) Qmax                       ![mg C/g soil], adsorption capacity
        REAL(KIND=8) Kba                        ![mg C/g soil/h], binding affinity
        REAL(KIND=8) Kdes                       ![mg DOC/h],desorption rate constant
        REAL(KIND=8) Kads                       ![mg DOC/mg DOC/h], adsorption rate constant = Kba*Kdes
        REAL(KIND=8) Kp2u                       ![1/h], rate from protected MOM to unprotected MOM (QOM)
        REAL(KIND=8) Ku2p                       ![1/h], rate from unprotected MOM (QOM) to protected MOM
        REAL(KIND=8) rENZP(const_nPOM)          ![1/h],turnover rate of ENZP
        REAL(KIND=8) rENZM                      ![1/h], turnover rate of ENZMAOC
        REAL(KIND=8) pENZP                      ![mg ENZP/mg MBC/h], production rate of ENZP
        REAL(KIND=8) pENZM                      ![mg ENZM/mg MBC/h], production rate of ENZMAOC
        REAL(KIND=8) frMB2POMO                  ![-], fraction of dead MMB to Oxidative POM
        REAL(KIND=8) frPOM2DOM                  ![-], fraction of decomposed POM allocated to DOC
        REAL(KIND=8) frMB2DOM                   ![-], fraction of dead MBC allocated to DOC
        REAL(KIND=8) Vg                         ![mg DOC/mg MBC/h], maximum uptake rate of DOC by MB
        REAL(KIND=8) alpha                      ![-], = Vm/(Vg+Vm)
        REAL(KIND=8) Vm                         ![1/h], specific microbial maintenance rate
        REAL(KIND=8) KsDOM                      ![mg DOC/cm3], half-saturation constant for uptake of DOC by MB
        REAL(KIND=8) Yg                         ![-], carbon use efficiency in uptake of DOC by MB
        REAL(KIND=8) Q10                        ![-], exponential in SWP scalar (fSWP0) for microbial mortality
        REAL(KIND=8) gamma                      ![-], rMORT = gama*Vm
        REAL(KIND=8) rMORT                      ![1/h], microbial mortality
        REAL(KIND=8) beta                       ![-], VmD = beta*Vm = dormant maintenance rate
        REAL(KIND=8) VmD                        ![1/h], VmD = beta*Vm, specific microbial maintenance rate for dormant microbes
        REAL(KIND=8) VmA2D                      ![1/h], dormancy rate
        REAL(KIND=8) VmD2A                      ![1/h], resuscitation rate
        REAL(KIND=8) SWP_A2D                    ![MPa], negative value, threshold SWP for microbial dormancy (e.g., -0.4 MPa)
        REAL(KIND=8) tau                        ![-], SWP_D2A = tau*SWP_A2D
        REAL(KIND=8) SWP_D2A                    ![MPa], negative value, threshold SWP for microbial resuscitation (e.g., 1/4*SWP_A2D)
        REAL(KIND=8) wdorm                      ![-], exponential in SWP scalar for A2D and D2A
         
        REAL(KIND=8) VNup_MB                    ![1/h], Mineral N uptake rate by MBA
        REAL(KIND=8) VNup_VG                    ![1/h], Mineral NH4 uptake rate by PLANT
        REAL(KIND=8) rNleach                    ![1/h], additional coefficient for NO3 leaching
        REAL(KIND=8) rPleach                    ![1/h], additional coefficient for PO4 leaching
		REAL(KIND=8) rN2O                        !#### [-], sacle factor of N2O
        REAL(KIND=8) bNup_VG                    ![-], exponential to scale GPP effect on plant N uptake, (GPP/GPPref)^bNup_VG
        
        REAL(KIND=8) VdP_amb                    !#### [1/h], Mineral P uptake rate by MBA
        REAL(KIND=8) VdP_plant                  !#### [1/h], Mineral P uptake rate by plant
        REAL(KIND=8) VdP_surf                   !#### [1/h], Mineral N uptake rate by labile P  
        
        REAL(KIND=8) KsNH4_MB                   ![mg N/cm3], half-saturation constant for NH4 uptake by MBA
        REAL(KIND=8) KsNO3_MB                   ![mg N/cm3], half-saturation constant for NO3 uptake by MBA
        
        REAL(KIND=8) KsP_amb                    !#### [mg P/cm3], half-saturation constant for P uptake by MBA
        REAL(KIND=8) KsP_plant                  !#### [mg P/cm3], half-saturation constant for P uptake by plant (VG)
        REAL(KIND=8) KsP_lab                    !#### [mg P/cm3], half-saturation constant for P uptake by Labile
        
        REAL(KIND=8) YgN                        ![-],Nitrogen Use Index (NUI), the fraction of N uptake being incorporated into microbial biomass
        REAL(KIND=8) wNUI                       ![-],exponent to compute NUI, i.e., YgN
        
        REAL(KIND=8) YgP                        !#### [-],P Use Index (NUI), the fraction of P uptake being incorporated into microbial biomass
        REAL(KIND=8) wPUI                       !#### [-],exponent to compute NUI, i.e., YgP 
        
        REAL(KIND=8) Qmax_NH4                   ![mg N/cm3], NH4 adsorption capacity
        REAL(KIND=8) Kba_NH4                    ![1/(mg N/cm3)], NH4 adsorption binding affinity
         
        REAL(KIND=8) KsNH4_VG                   ![mg N/cm3], half-saturation constant for NH4 uptake by plant (VG)
        REAL(KIND=8) KsNO3_VG                   ![mg N/cm3], half-saturation constant for NO3 uptake by plant (VG)
        
        REAL(KIND=8) fpENZN                     ![mg ENZ/mg MBC/h], production rate of ENZs for mineral N
        REAL(KIND=8) VNif                       ![1/h], N fixation rate
        REAL(KIND=8) VNit                       ![1/h], nitrification rate
        REAL(KIND=8) VDenit(const_nDenit)       ![1/h], denitrification rate
        
        REAL(KIND=8) fpPTASE                     !#### [mg ENZ/mg MBC/h], production rate of ENZs for mineral P
        REAL(KIND=8) VdP(3)                     !#### [1/h], decomposition P rate
         
        REAL(KIND=8) KsNif                       ![mg N/cm3], half-saturation constant
        REAL(KIND=8) KsNit                       ![mg N/cm3], half-saturation constant
        REAL(KIND=8) KsDenit(const_nDenit)       ![mg N/cm3], half-saturation constant
        
        REAL(KIND=8) KsP(3)                     !#### [mg P/cm3], half-saturation constant
        
        
        REAL(KIND=8)micr_bio2enz                !#### biomass to enzyme
        REAL(KIND=8)enz_plant_p                 !#### transport enzyme
        REAL(KIND=8)weath_rate                  !#### weathering rate
        REAL(KIND=8)occlu_rate                  !#### occlusion rate
        REAL(KIND=8)p_adsorp_rate               !#### adsorption rate
        REAL(KIND=8)p_desorp_rate               !#### deadsorption rate
        
    END TYPE sMEND_PAR
    !-----------------------------------------------------------------------------!
    TYPE sMEND_INP
        REAL(KIND=8) GPP                                    ![mgC/cm3/h],GPP
        REAL(KIND=8) SIN                                    ![mgC/cm3/h],SOC input, e.g., litter
        REAL(KIND=8) SIN3                                   ![mgC/cm3/h],type3 SOC input, e.g., roots
        REAL(KIND=8) tmp                                    ![C],temperature
        REAL(KIND=8) pH                                     ![-],pH value
        REAL(KIND=8) SWC                                    ![m3/m3], soil water content
        REAL(KIND=8) SWP                                    ![MPa],soil water potential [-]
        REAL(KIND=8) dt                                     ![h], time-step, dt = 0.5 h for CLM

        !         TYPE(sMEND_CADD) CADD                          !external carbon input
        !         TYPE(sMEND_CADD) CADDI(const_nISO)             !external carbon input (isotopes)
        TYPE(sMEND_POOL)     ::  CPOOL                   !organic C pools
        TYPE(sMEND_POOL)     ::  NPOOL                   !organic N pools
        TYPE(sMEND_POOL)     ::  PPOOL                   !#### organic P pools
        TYPE(sMEND_POOL)     ::  CN                      !C:N Ratio in OM pools
        TYPE(sMEND_POOL)     ::  CN0                         !Initial C:N Ratio in OM pools
        TYPE(sMEND_POOL)     ::  CP                      !#### C:P Ratio in OM pools
        TYPE(sMEND_POOL)     ::  CP0                     !#### Initial C:P Ratio in OM pools
        TYPE(sMEND_POOL_MN)  ::  MNPOOL                  !mineral pools
    !         TYPE(sMEND_CPOOL) CPOOLI(const_nISO)           !isotopic carbon pools
    !         TYPE(sMEND_CPOOL) CPOOLIFR(const_nISO)         !isotopic fraction in carbon pools
    !         TYPE(sMEND_CPOOL) CPOOLI_SIG(const_nISO - 1)   ![��],isotopic signature in carbon pools
    END TYPE sMEND_INP
    !-----------------------------------------------------------------------------!
    TYPE sMEND_OUT
        TYPE(sMEND_POOL)     ::  CPOOL                   !organic C pools
        TYPE(sMEND_POOL)     ::  NPOOL                   !organic N pools
        TYPE(sMEND_POOL)     ::  PPOOL                   !#### organic P pools
        TYPE(sMEND_POOL)     ::  CN                      !C:N Ratio in OM pools
        TYPE(sMEND_POOL)     ::  CP                      !#### C:P Ratio in OM pools
        TYPE(sMEND_POOL_MN)  ::  MNPOOL                  !mineral pools
        !         TYPE(sMEND_CPOOL) CPOOLI(const_nISO)           !isotopic carbon pools
        !         TYPE(sMEND_CPOOL) CPOOLIFR(const_nISO)         !isotopic fraction in carbon pools
        !         TYPE(sMEND_CPOOL) CPOOLI_SIG(const_nISO - 1)   ![��],isotopic signature in carbon pools
        TYPE(sMEND_FLUX)     ::  CFLUX                   !CARBON FLUXES
        TYPE(sMEND_FLUX)     ::  NFLUX                   !N FLUXES
        TYPE(sMEND_FLUX)     ::  PFLUX                   !#### P FLUXES
        TYPE(sMEND_FLUX_MN)  ::  MNFLUX                  !mineralization/immobilization fluxes
         
        REAL(KIND=8) CO2_gm_iso                             !C14 or C13 flux
    !          REAL(KIND=8) TOCinp                                 !total external C input
    !          REAL(KIND=8) TOCout                                 !total output: respiration
    !          REAL(KIND=8) TOCbeg                                 !total mass at the beginning
    !          REAL(KIND=8) TOCend                                 !total mass at the end
    !          REAL(KIND=8) RE                                     ![-],default = 0, mass balance check = (POOL_end - POOL_beg)-(TOTinp - TOTout)
    END TYPE sMEND_OUT
      
    !-----------------------------------------------------------------------------!
    TYPE sMEND_INI
        TYPE(sMEND_INP) :: sINP                        !INPUT
        INTEGER pid
        INTEGER nprocs
        CHARACTER(len=10) SITE                         !!Site Name, used for prefix of output files
        CHARACTER(len=3) BIOME                         !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
        CHARACTER(len=3) SOM                           !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
        CHARACTER(len=100) dirinp                       !!input dir
        CHARACTER(len=100) dirout                       !!output dir
        character(len=8) sDate_beg_all                 !!"yyyymmdd": available input data: begin date for optimizaton
        character(len=8) sDate_end_all                 !!"yyyymmdd": available input data: end date for optimizaton
        character(len=8) sDate_beg_sim                 !!"yyyymmdd": simulation: begin date
        character(len=8) sDate_end_sim                 !!"yyyymmdd": simulation: end date
        character(len=8) sDate_beg_inp2                !!"yyyymmdd": begin date for available constant input data
        character(len=8) sDate_end_inp2                !!"yyyymmdd": end date for available constant input data
        CHARACTER(len=100) SOIL_INI_file                !soil initialization file
        CHARACTER(len=20),ALLOCATABLE:: VARfile(:)     !VAR data file
        CHARACTER(len=4),ALLOCATABLE:: VARobj(:)       !VAR objective function type (NSEC or MARE)
        LOGICAL Carbon_only                            !.true. for Carbon_only model, .false. for C-N coupled model
        INTEGER iError                                 !iError = -11 denotes 1st C balance error greater than tolerance; -12: 1st N balance error; -21 & -22: 2nd C & N balance check
        INTEGER iModel                                 !0: run model; 1: run model optimization; 2: sensitivity
        INTEGER iSA_range(3)                           !iSA_range(1:2) difine the range of rows to compute obj for iModel=3; iSA_range(3) =0 or 1 denotes continuous or discontinuous id
        INTEGER nPar                                   !# of parameters to be optimized
        INTEGER iKinetics                              !decomposition kinetics: 0-Michaelis-Menten, 1-First Order, 2-Second Order
        INTEGER iENZN_Allocate                         !Nitrogen Enzyme Allocation: 0-weighting factor = N pool, 1-weighting factor = N-Pool/Half-saturation-constant
        INTEGER iHR                                    !HR calculation method, fraction of [0: potential uptake, 1: actual uptake]
        INTEGER iTmp_Func                              !sINI%iTmp_Func: temperature response function [0: Arrhenius Equation, 1: Q10 Equation]
        INTEGER nVARopt                                !# of observed response variables for optimization
        CHARACTER(len=20):: Name_VARopt(const_nVAR0)   !Names of Variables selected for comparison/optimization
        INTEGER nOBS_tot                               !total # of observations for those output variables, = sum(VARopt_int(:, 2) )
        INTEGER nYear                                  !simulation years (no-use)
        INTEGER nHour                                  !available input data period length
        INTEGER nHour_sim                              !simulation period length
        INTEGER nOutStep                               !output time interval (step)
        INTEGER iFout_SIM_obs                          !file unit for response variables matching observations
        INTEGER iFout_SIM_day                          !file unit for response variables, mean daily output
        INTEGER iFout_SIM_mon                          !file unit for response variables, mean monthly output
        INTEGER iFout_VAR_hour                          !file unit for all response variables, user-defined output time-step
        INTEGER iFout_FLX_hour                          !file unit for all flux variables, user-defined output time-step
        INTEGER iFout_RATE_hour                        !file unit for derived rates, e.g., 1st-order decomposition rate, active fraction, phi
        INTEGER iFout_ITW_hour                         !file unit for Input (e.g., litter), Temperature, Water content & potential
        INTEGER iFout_PAR_hour                         !file unit for hourly parameters modified by Temperature, Water potential, and other factors
        REAL(KIND=8) r0                                     !initial active-microbe fraction
        REAL(KIND=8) LCI0                                   !initial lignocellulose index = Lignin/(Lignin+Cellulose)!file unit for output
        REAL(KIND=8) soilDepth                              ![cm],soil depth
        REAL(KIND=8) porosity                               ![cm3/cm3],soil porosity = 1 - Bulk_Density/2.65
        REAL(KIND=8) fDOM                                   ![-] fDOM = DOM/DOMT, DOMT = DOM + DOMS
        REAL(KIND=8) SWCFC                                  ![cm3/cm3], soil water content at field capacity
        REAL(KIND=8) Ksat                                   ![cm/h], saturated hydraulic conductivity
        REAL(KIND=8) Lambda                                 ![-], Slope of logarithmic tension-moisture curve = 1/B
        REAL(KIND=8) fNO3_Leaching                          ![-], fraction of NO3 leaching
        REAL(KIND=8) fPO4_Leaching                          ![-], fraction of PO4 leaching
        REAL(KIND=8) VNup_VG_GPPscaler                       ![-], fGPP scaler for plant N uptake
        REAL(KIND=8) altitude                               ![m], altitude, elevation
        REAL(KIND=8) fTexture(3)                            !!fraction of SAND,SILT,CLAY
        REAL(KIND=8) SIN_C12_C14                            !ratio of C12 to C14 in SOC input
        REAL(KIND=8) GPPref                                 ![mg C/cm3/h],reference GPP, set to mean GPP under control condition
        INTEGER iGPPscaler                             !Function to scale plant N uptake rate; iGPPscaler=0: (GPP/GPPref)**b; 1: exp[(GPP/GPPref - 1)*b]
         
        REAL(KIND=8) CN_MB_avg                              ![mg C/mg N], fixed C:N ratio for MB
        REAL(KIND=8) CN_MB_min                              ![mg C/mg N], min C:N ratio for MB
        REAL(KIND=8) CN_MB_max                              ![mg C/mg N], max C:N ratio for MB
        REAL(KIND=8) CN_ENZM                                ![mg C/mg N], fixed C:N ratio for ENZYME
        REAL(KIND=8) CN_ENZP(const_nPOM)                    ! [mg C/mg N], fixed C:N ratio for Ligninases & Cellulases
        
        REAL(KIND=8) CP_MB_avg                              !####[mg C/mg P], fixed C:P ratio for MB
        REAL(KIND=8) CP_MB_min                              !####[mg C/mg P], min C:P ratio for MB
        REAL(KIND=8) CP_MB_max                              !####[mg C/mg P], max C:P ratio for MB
        REAL(KIND=8) CP_ENZM                                !####[mg C/mg P], fixed C:P ratio for ENZYME
        REAL(KIND=8) CP_ENZP(const_nPOM)                    !#### [mg C/mg P], fixed C:P ratio for Ligninases & Cellulases
         
        REAL(KIND=8) CN_LITT_avg                             ![mg C/mg N], C:N ratio for LITTER INPUT
        REAL(KIND=8) CN_WOOD_avg                             ![mg C/mg N], C:N ratio for COARSE WOOD INPUT
        REAL(KIND=8) CN_ROOT_avg                             ![mg C/mg N], C:N ratio for ROOT INPUT
        
        REAL(KIND=8) CP_LITT_avg                             !####[mg C/mg P], C:P ratio for LITTER INPUT
        REAL(KIND=8) CP_WOOD_avg                             !####[mg C/mg P], C:P ratio for COARSE WOOD INPUT
        REAL(KIND=8) CP_ROOT_avg                             !####[mg C/mg P], C:P ratio for ROOT INPUT
         
        REAL(KIND=8) CN_LITT(3)                             ![mg C/mg N], C:N ratio for LITTER INPUT: LIGNIN, CELLULOSE, LABILE
        REAL(KIND=8) CN_WOOD(3)                             ![mg C/mg N], C:N ratio for COARSE WOOD INPUT: LIGNIN, CELLULOSE, LABILE
        REAL(KIND=8) CN_ROOT(3)                             ![mg C/mg N], C:N ratio for ROOT INPUT: LIGNIN, CELLULOSE, LABILE
        
        REAL(KIND=8) CP_LITT(3)                             !####[mg C/mg P], C:P ratio for LITTER INPUT: LIGNIN, CELLULOSE, LABILE
        REAL(KIND=8) CP_WOOD(3)                             !####[mg C/mg P], C:P ratio for COARSE WOOD INPUT: LIGNIN, CELLULOSE, LABILE
        REAL(KIND=8) CP_ROOT(3)                             !####[mg C/mg P], C:P ratio for ROOT INPUT: LIGNIN, CELLULOSE, LABILE
         
        CHARACTER(len=20) Name_Pool(const_nPOOL)     !Names for C/N Pools
        CHARACTER(len=20) Name_MNPOOL(const_nPOOL_MN)
        CHARACTER(len=20) Name_FLUX(const_nFLUX)
        CHARACTER(len=20) Name_MNFLUX(const_nFLUX_MN)
        CHARACTER(len=20) Name_PAR(const_nPAR)
        CHARACTER(len=20) Name_PAR0(const_nPAR0)   !parameter names for xx(:)
        CHARACTER(len=20) Name_RATE(const_nRATE)
         
        CHARACTER(len=200) sFile_VAR_hour                          !file unit for all response variables, user-defined output time-step
        CHARACTER(len=200) sFile_FLX_hour                          !file unit for all flux variables, user-defined output time-step
        CHARACTER(len=200) sFile_RATE_hour                        !file unit for derived rates, e.g., 1st-order decomposition rate, active fraction, phi
        CHARACTER(len=200) sFile_PAR_hour

        INTEGER iFout_all                                  !output file1: print the results of SCEUA
        INTEGER iFout_end                                  !output file2: summarize optimal parameter values from nRun
        INTEGER iFout_ini                                  !output file3: other text output
        INTEGER iFout_mcmc                                 !output file for MCMC
         
        REAL(KIND=8), ALLOCATABLE:: dINI(:)                 !Array to store initial values read from 'soil.ini'
        REAL(KIND=8), ALLOCATABLE:: rOBJ(:)                 !Nash-Sutcliffe Efficency Coefficient, rNSE(1): mean NSEC, rNSE(2:nObs_var+1): individual NSEC for nObs_var
        REAL(KIND=8), ALLOCATABLE:: rOBJw(:)                !obj weighting factor
        REAL(KIND=8), ALLOCATABLE:: dOBS_opt(:, :)          !date,obs,iVARopt
        REAL(KIND=8), ALLOCATABLE:: dSIM_opt(:, :)          !date,sim,sim_sd
        REAL(KIND=8), ALLOCATABLE:: STP(:)                  ![C], soil temperature
        REAL(KIND=8), ALLOCATABLE:: SWC(:)                  ![fraction],soil water content
        REAL(KIND=8), ALLOCATABLE:: SWP(:)                  ![MPa], soil water potential
        REAL(KIND=8), ALLOCATABLE:: SpH(:)                  !soil pH
        REAL(KIND=8), ALLOCATABLE:: SIN(:)                  ![mgC/cm3/h],SOC input, e.g., litterfall
        REAL(KIND=8), ALLOCATABLE:: SIN3(:)                 ![mgC/cm3/h],type3 SOC input, e.g., roots
        REAL(KIND=8) SIN_other(2,3)                         ![mgC/cm3/h | mgC/cm3], other constant inputs to 3 pools, e.g., coarse wood, roots, SIN_other(1,:): between sDate_beg_all & sDate_end_all; SIN_other(2,:): between sDate_beg_inp2 & sDate_end_inp2;
        REAL(KIND=8) SIN_frac(3)                            !fraction of SOC input to 3 pools (POM1, POM2, & DOC)
        REAL(KIND=8) SIN3_frac(3)                           !fraction of type3 SOC input to 3 pools (POM1, POM2, & DOC)
        REAL(KIND=8) SIN_Multiplier                         !multiplier for litter input during post-data period simulation: SIN*SIN_multiplier
        REAL(KIND=8), ALLOCATABLE:: SIN_NH4(:)              ![mgN/cm3/h],NH4 INPUT RATE
        REAL(KIND=8), ALLOCATABLE:: SIN_NO3(:)              ![mgN/cm3/h],NO3 INPUT RATE
        REAL(KIND=8), ALLOCATABLE:: SIN_P(:)                ![mgP/cm3/h],P   INPUT RATE 
         
        REAL(KIND=8), ALLOCATABLE:: VARobj_tol(:)           !VAR objective function tolerance
        REAL(KIND=8), ALLOCATABLE:: VARobjw(:)              !VAR objective function weighting factor (any number, will be normalized)
        INTEGER, ALLOCATABLE:: VARopt(:)               !if use observed VAR to calibrate model, see MEND.ini
        INTEGER, ALLOCATABLE:: VARstep(:)              !time-step, 0(hourly),1(daily),2(monthly),3(yearly)
        INTEGER, ALLOCATABLE:: VARcol(:)               !VAR data in which column of input file
        INTEGER, ALLOCATABLE:: VARopt_int(:, :)        !nrow = nVARopt, ncol = 3 (index of the VAR for opt, no. of data, tstep); tstep = 0(hourly),1(daily),2(monthly),3(yearly)
         
        INTEGER iScenario                             !Scenario design: 1(1-yr mean hourly data), 2(multiple-year hourly data)
        REAL(KIND=8) waterRetention_vg(4)                   !vg_SWCres,vg_SWCsat,vg_alpha,vg_n  !!Line-28, van-Genuchten equation
        REAL(KIND=8) STP_delta                              ![C], annual change in temperature
        REAL(KIND=8) SWC_logis(3)                           !Soil Water Content (SWC) parameter in logistic equation:  (1) p <SWC_t=infinity / SWC0>, (2) r <steepness>, (3) t0 <reference year>
        REAL(KIND=8) SIN_logis(4)                           !Litter Input (SIN) Parameters in logistic equation:       (1) beta0 <intercept>, (2)beta1 <steepness>, (3)t0 <reference year>, (4) fDOC_delta <annual change of DOC fraction in SOC input>
    END TYPE sMEND_INI
      
!-----------------------------------------------------------------------------! 
END MODULE MOD_MEND_TYPE 

