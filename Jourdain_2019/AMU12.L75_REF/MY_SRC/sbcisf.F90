MODULE sbcisf
   !!======================================================================
   !!                       ***  MODULE  sbcisf  ***
   !! Surface module :  update surface ocean boundary condition under ice
   !!                   shelf
   !!======================================================================
   !! History :  3.2   !  2011-02  (C.Harris  ) Original code isf cav
   !!            X.X   !  2006-02  (C. Wang   ) Original code bg03
   !!            3.4   !  2013-03  (P. Mathiot) Merging + parametrization
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_isf        : update sbc under ice shelf
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE eosbn2          ! equation of state
   USE sbc_oce         ! surface boundary condition: ocean fields
   USE lbclnk          !
   USE iom             ! I/O manager library
   USE in_out_manager  ! I/O manager
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing
   USE lib_fortran     ! glob_sum
   USE zdfbfr
   USE fldread         ! read input field at current time step



   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_isf, sbc_isf_div, sbc_isf_alloc  ! routine called in sbcmod and divcur

   ! public in order to be able to output then 

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   risf_tsc_b, risf_tsc   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qisf              !: net heat flux from ice shelf
   REAL(wp), PUBLIC ::   rn_hisf_tbl                 !: thickness of top boundary layer [m]
   INTEGER , PUBLIC ::   nn_isf                      !: flag to choose between explicit/param/specified  
   INTEGER , PUBLIC ::   nn_isfblk                   !: 
   INTEGER , PUBLIC ::   nn_gammablk                 !:
   REAL(wp), PUBLIC ::   rn_gammat0                  !: temperature exchange coeficient
   REAL(wp), PUBLIC ::   rn_gammas0                  !: salinity    exchange coeficient 

   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  rzisf_tbl              !:depth of calving front (shallowest point) nn_isf ==2/3
   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  rhisf_tbl, rhisf_tbl_0 !:thickness of tbl
   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  r1_hisf_tbl            !:1/thickness of tbl
   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  ralpha                 !:proportion of bottom cell influenced by tbl 
   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  risfLeff               !:effective length (Leff) BG03 nn_isf==2
   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  ttbl, stbl, utbl, vtbl !:top boundary layer variable at T point
#if defined key_agrif
   ! AGRIF can not handle these arrays as integers. The reason is a mystery but problems avoided by declaring them as reals
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  misfkt, misfkb         !:Level of ice shelf base
                                                                                          !: (first wet level and last level include in the tbl)
#else
   INTEGER,    PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  misfkt, misfkb         !:Level of ice shelf base
#endif


   REAL(wp), PUBLIC, SAVE ::   rcpi   = 2000.0_wp     ! phycst ?
   REAL(wp), PUBLIC, SAVE ::   kappa  = 1.54e-6_wp    ! phycst ?
   REAL(wp), PUBLIC, SAVE ::   rhoisf = 920.0_wp      ! phycst ?
   REAL(wp), PUBLIC, SAVE ::   tsurf  = -20.0_wp      ! phycst ?
   REAL(wp), PUBLIC, SAVE ::   lfusisf= 0.334e6_wp    ! phycst ?

!: Variable used in fldread to read the forcing file (nn_isf == 4 .OR. nn_isf == 3)
   CHARACTER(len=100), PUBLIC ::   cn_dirisf  = './'    !: Root directory for location of ssr files
   TYPE(FLD_N)       , PUBLIC ::   sn_fwfisf            !: information about the isf file to be read
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::  sf_fwfisf
   TYPE(FLD_N)       , PUBLIC ::   sn_rnfisf              !: information about the runoff file to be read
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_rnfisf           
   TYPE(FLD_N)       , PUBLIC ::   sn_depmax_isf, sn_depmin_isf, sn_Leff_isf     !: information about the runoff file to be read
   
   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.0 , LOCEAN-IPSL (2008)
   !! $Id: sbcisf.F90 5905 2015-11-20 16:59:57Z mathiot $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
 
  SUBROUTINE sbc_isf(kt)
    INTEGER, INTENT(in) :: kt                   ! ocean time step
    INTEGER             :: ji, jj, jk           ! loop index
    INTEGER             :: ikt, ikb             ! top and bottom level of the isf boundary layer
    INTEGER             :: inum, ierror
    INTEGER             :: ios                  ! Local integer output status for namelist read
    REAL(wp)            :: zhk
    REAL(wp), DIMENSION (:,:), POINTER :: zt_frz, zdep ! freezing temperature (zt_frz) at depth (zdep) 
    CHARACTER(len=256)  :: cvarzisf, cvarhisf   ! name for isf file
    CHARACTER(LEN=32 )  :: cvarLeff             ! variable name for efficient Length scale
      !
      !!---------------------------------------------------------------------
      NAMELIST/namsbc_isf/ nn_isfblk, rn_hisf_tbl, rn_gammat0, rn_gammas0, nn_gammablk, nn_isf      , &
                         & sn_fwfisf, sn_rnfisf, sn_depmax_isf, sn_depmin_isf, sn_Leff_isf
      !
      ! allocation
      CALL wrk_alloc( jpi,jpj, zt_frz, zdep  )
      !
      !                                         ! ====================== !
      IF( kt == nit000 ) THEN                   !  First call kt=nit000  !
         !                                      ! ====================== !
         REWIND( numnam_ref )              ! Namelist namsbc_rnf in reference namelist : Runoffs 
         READ  ( numnam_ref, namsbc_isf, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_isf in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist namsbc_rnf in configuration namelist : Runoffs
         READ  ( numnam_cfg, namsbc_isf, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_isf in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namsbc_isf )


         IF ( lwp ) WRITE(numout,*)
         IF ( lwp ) WRITE(numout,*) 'sbc_isf: heat flux of the ice shelf'
         IF ( lwp ) WRITE(numout,*) '~~~~~~~~~'
         IF ( lwp ) WRITE(numout,*) 'sbcisf :' 
         IF ( lwp ) WRITE(numout,*) '~~~~~~~~'
         IF ( lwp ) WRITE(numout,*) '        nn_isf      = ', nn_isf
         IF ( lwp ) WRITE(numout,*) '        nn_isfblk   = ', nn_isfblk
         IF ( lwp ) WRITE(numout,*) '        rn_hisf_tbl = ', rn_hisf_tbl
         IF ( lwp ) WRITE(numout,*) '        nn_gammablk = ', nn_gammablk 
         IF ( lwp ) WRITE(numout,*) '        rn_gammat0  = ', rn_gammat0  
         IF ( lwp ) WRITE(numout,*) '        rn_gammas0  = ', rn_gammas0  
         IF ( lwp ) WRITE(numout,*) '        rn_tfri2    = ', rn_tfri2 
         !
         ! Allocate public variable
         IF ( sbc_isf_alloc()  /= 0 )         CALL ctl_stop( 'STOP', 'sbc_isf : unable to allocate arrays' )
         !
         ! initialisation
         qisf(:,:)        = 0._wp  ; fwfisf  (:,:) = 0._wp
         risf_tsc(:,:,:)  = 0._wp  ; fwfisf_b(:,:) = 0._wp
         !
         ! define isf tbl tickness, top and bottom indice
         IF      (nn_isf == 1) THEN
            rhisf_tbl(:,:) = rn_hisf_tbl
            misfkt(:,:)    = mikt(:,:)         ! same indice for bg03 et cav => used in isfdiv
         ELSE IF ((nn_isf == 3) .OR. (nn_isf == 2)) THEN
            ALLOCATE( sf_rnfisf(1), STAT=ierror )
            ALLOCATE( sf_rnfisf(1)%fnow(jpi,jpj,1), sf_rnfisf(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_rnfisf, (/ sn_rnfisf /), cn_dirisf, 'sbc_isf_init', 'read fresh water flux isf data', 'namsbc_isf' )

            !: read effective lenght (BG03)
            IF (nn_isf == 2) THEN
               cvarLeff = 'soLeff' 
               CALL iom_open( sn_Leff_isf%clname, inum )
               CALL iom_get( inum, jpdom_data, cvarLeff, risfLeff , 1)
               CALL iom_close(inum)
               !
               risfLeff = risfLeff*1000.0_wp           !: convertion in m
            END IF

           ! read depth of the top and bottom of the isf top boundary layer (in this case, isf front depth and grounding line depth)
            CALL iom_open( sn_depmax_isf%clname, inum )
            cvarhisf = TRIM(sn_depmax_isf%clvar)
            CALL iom_get( inum, jpdom_data, cvarhisf, rhisf_tbl, 1) !: depth of deepest point of the ice shelf base
            CALL iom_close(inum)
            !
            CALL iom_open( sn_depmin_isf%clname, inum )
            cvarzisf = TRIM(sn_depmin_isf%clvar)
            CALL iom_get( inum, jpdom_data, cvarzisf, rzisf_tbl, 1) !: depth of shallowest point of the ice shelves base
            CALL iom_close(inum)
            !
            rhisf_tbl(:,:) = rhisf_tbl(:,:) - rzisf_tbl(:,:)        !: tickness isf boundary layer

           !! compute first level of the top boundary layer
           DO ji = 1, jpi
              DO jj = 1, jpj
                  jk = 2
                  DO WHILE ( jk .LE. mbkt(ji,jj) .AND. fsdepw(ji,jj,jk) < rzisf_tbl(ji,jj) ) ;  jk = jk + 1 ;  END DO
                  misfkt(ji,jj) = jk-1
               END DO
            END DO

         ELSE IF ( nn_isf == 4 ) THEN
            ! as in nn_isf == 1
            rhisf_tbl(:,:) = rn_hisf_tbl
            misfkt(:,:)    = mikt(:,:)         ! same indice for bg03 et cav => used in isfdiv
            
            ! load variable used in fldread (use for temporal interpolation of isf fwf forcing)
            ALLOCATE( sf_fwfisf(1), STAT=ierror )
            ALLOCATE( sf_fwfisf(1)%fnow(jpi,jpj,1), sf_fwfisf(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_fwfisf, (/ sn_fwfisf /), cn_dirisf, 'sbc_isf_init', 'read fresh water flux isf data', 'namsbc_isf' )
         END IF
         
         rhisf_tbl_0(:,:) = rhisf_tbl(:,:)

         ! compute bottom level of isf tbl and thickness of tbl below the ice shelf
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkt(ji,jj)
               ! thickness of boundary layer at least the top level thickness
               rhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), fse3t_n(ji,jj,ikt))

               ! determine the deepest level influenced by the boundary layer
               DO jk = ikt+1, mbkt(ji,jj)
                  IF ( (SUM(fse3t_n(ji,jj,ikt:jk-1)) .LT. rhisf_tbl(ji,jj)) .AND. (tmask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               rhisf_tbl(ji,jj) = MIN(rhisf_tbl(ji,jj), SUM(fse3t_n(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.
               misfkb(ji,jj) = ikb                                                  ! last wet level of the tbl
               r1_hisf_tbl(ji,jj) = 1._wp / rhisf_tbl(ji,jj)

               zhk           = SUM( fse3t(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj)  ! proportion of tbl cover by cell from ikt to ikb - 1
               ralpha(ji,jj) = rhisf_tbl(ji,jj) * (1._wp - zhk ) / fse3t(ji,jj,ikb)  ! proportion of bottom cell influenced by boundary layer
            END DO
         END DO
         
      END IF

      !                                            ! ---------------------------------------- !
      IF( kt /= nit000 ) THEN                      !          Swap of forcing fields          !
         !                                         ! ---------------------------------------- !
         fwfisf_b  (:,:  ) = fwfisf  (:,:  )               ! Swap the ocean forcing fields except at nit000
         risf_tsc_b(:,:,:) = risf_tsc(:,:,:)               ! where before fields are set at the end of the routine
         !
      ENDIF

      IF( MOD( kt-1, nn_fsbc) == 0 ) THEN


         ! compute salf and heat flux
         IF (nn_isf == 1) THEN
            ! realistic ice shelf formulation
            ! compute T/S/U/V for the top boundary layer
            CALL sbc_isf_tbl(tsn(:,:,:,jp_tem),ttbl(:,:),'T')
            CALL sbc_isf_tbl(tsn(:,:,:,jp_sal),stbl(:,:),'T')
            CALL sbc_isf_tbl(un(:,:,:)        ,utbl(:,:),'U')
            CALL sbc_isf_tbl(vn(:,:,:)        ,vtbl(:,:),'V')
            ! iom print
            CALL iom_put('ttbl',ttbl(:,:))
            CALL iom_put('stbl',stbl(:,:))
            CALL iom_put('utbl',utbl(:,:))
            CALL iom_put('vtbl',vtbl(:,:))
            ! compute fwf and heat flux
            CALL sbc_isf_cav (kt)

         ELSE IF (nn_isf == 2) THEN
            ! Beckmann and Goosse parametrisation 
            stbl(:,:)   = soce
            CALL sbc_isf_bg03(kt)

         ELSE IF (nn_isf == 3) THEN
            ! specified runoff in depth (Mathiot et al., XXXX in preparation)
            CALL fld_read ( kt, nn_fsbc, sf_rnfisf   )
            fwfisf(:,:) = - sf_rnfisf(1)%fnow(:,:,1)         ! fwf  flux from the isf (fwfisf <0 mean melting) 
            qisf(:,:)   = fwfisf(:,:) * lfusisf              ! heat flux
            stbl(:,:)   = soce

         ELSE IF (nn_isf == 4) THEN
            ! specified fwf and heat flux forcing beneath the ice shelf
            CALL fld_read ( kt, nn_fsbc, sf_fwfisf   )
            fwfisf(:,:) = - sf_fwfisf(1)%fnow(:,:,1)           ! fwf  flux from the isf (fwfisf <0 mean melting)
            qisf(:,:)   = fwfisf(:,:) * lfusisf                ! heat flux
            stbl(:,:)   = soce
         END IF

         ! compute tsc due to isf
         ! isf melting implemented as a volume flux and we assume that melt water is at 0 PSU.
         ! WARNING water add at temp = 0C, need to add a correction term (fwfisf * tfreez / rau0).
         ! compute freezing point beneath ice shelf (or top cell if nn_isf = 3)
         DO jj = 1,jpj
            DO ji = 1,jpi
               zdep(ji,jj)=fsdepw_n(ji,jj,misfkt(ji,jj))
            END DO
         END DO
         CALL eos_fzp( stbl(:,:), zt_frz(:,:), zdep(:,:) )
         
         risf_tsc(:,:,jp_tem) = qisf(:,:) * r1_rau0_rcp - fwfisf(:,:) * zt_frz(:,:) * r1_rau0 !
         risf_tsc(:,:,jp_sal) = 0.0_wp

         ! lbclnk
         CALL lbc_lnk(risf_tsc(:,:,jp_tem),'T',1.)
         CALL lbc_lnk(risf_tsc(:,:,jp_sal),'T',1.)
         CALL lbc_lnk(fwfisf(:,:)   ,'T',1.)
         CALL lbc_lnk(qisf(:,:)     ,'T',1.)

         IF( kt == nit000 ) THEN                         !   set the forcing field at nit000 - 1    !
            IF( ln_rstart .AND.    &                     ! Restart: read in restart file
                 & iom_varid( numror, 'fwf_isf_b', ldstop = .FALSE. ) > 0 ) THEN
               IF(lwp) WRITE(numout,*) '          nit000-1 isf tracer content forcing fields read in the restart file'
               CALL iom_get( numror, jpdom_autoglo, 'fwf_isf_b', fwfisf_b(:,:) )   ! before salt content isf_tsc trend
               CALL iom_get( numror, jpdom_autoglo, 'isf_sc_b', risf_tsc_b(:,:,jp_sal) )   ! before salt content isf_tsc trend
               CALL iom_get( numror, jpdom_autoglo, 'isf_hc_b', risf_tsc_b(:,:,jp_tem) )   ! before salt content isf_tsc trend
            ELSE
               fwfisf_b(:,:)    = fwfisf(:,:)
               risf_tsc_b(:,:,:)= risf_tsc(:,:,:)
            END IF
         ENDIF
         ! 
         ! output
         IF( iom_use('qisf'  ) )   CALL iom_put('qisf'  , qisf)
         IF( iom_use('fwfisf') )   CALL iom_put('fwfisf', fwfisf)
      END IF

      ! deallocation
      CALL wrk_dealloc( jpi,jpj, zt_frz, zdep  )
  
  END SUBROUTINE sbc_isf

  INTEGER FUNCTION sbc_isf_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION sbc_isf_rnf_alloc  ***
      !!----------------------------------------------------------------------
      sbc_isf_alloc = 0       ! set to zero if no array to be allocated
      IF( .NOT. ALLOCATED( qisf ) ) THEN
         ALLOCATE(  risf_tsc(jpi,jpj,jpts), risf_tsc_b(jpi,jpj,jpts), qisf(jpi,jpj)   , &
               &    rhisf_tbl(jpi,jpj)    , r1_hisf_tbl(jpi,jpj), rzisf_tbl(jpi,jpj)  , &
               &    ttbl(jpi,jpj)         , stbl(jpi,jpj)       , utbl(jpi,jpj)       , &
               &    vtbl(jpi, jpj)        , risfLeff(jpi,jpj)   , rhisf_tbl_0(jpi,jpj), &
               &    ralpha(jpi,jpj)       , misfkt(jpi,jpj)     , misfkb(jpi,jpj)     , &
               &    STAT= sbc_isf_alloc )
         !
         IF( lk_mpp                  )   CALL mpp_sum ( sbc_isf_alloc )
         IF( sbc_isf_alloc /= 0 )   CALL ctl_warn('sbc_isf_alloc: failed to allocate arrays.')
         !
      ENDIF
  END FUNCTION

  SUBROUTINE sbc_isf_bg03(kt)
   !!==========================================================================
   !!                 *** SUBROUTINE sbcisf_bg03  ***
   !! add net heat and fresh water flux from ice shelf melting
   !! into the adjacent ocean using the parameterisation by
   !! Beckmann and Goosse (2003), "A parameterization of ice shelf-ocean
   !!     interaction for climate models", Ocean Modelling 5(2003) 157-170.
   !!  (hereafter BG)
   !!==========================================================================
   !!----------------------------------------------------------------------
   !!   sbc_isf_bg03      : routine called from sbcmod
   !!----------------------------------------------------------------------
   !!
   !! ** Purpose   :   Add heat and fresh water fluxes due to ice shelf melting
   !! ** Reference :   Beckmann et Goosse, 2003, Ocean Modelling
   !!
   !! History :
   !!      !  06-02  (C. Wang) Original code
   !!----------------------------------------------------------------------

    INTEGER, INTENT ( in ) :: kt

    INTEGER :: ji, jj, jk !temporary integer

    REAL(wp) :: zt_sum    ! sum of the temperature between 200m and 600m
    REAL(wp) :: zt_ave    ! averaged temperature between 200m and 600m
    REAL(wp) :: zt_frz    ! freezing point temperature at depth z
    REAL(wp) :: zpress    ! pressure to compute the freezing point in depth
    
    !!----------------------------------------------------------------------
    IF ( nn_timing == 1 ) CALL timing_start('sbc_isf_bg03')
     !

    ! This test is false only in the very first time step of a run (JMM ???- Initialy build to skip 1rst year of run )
    DO ji = 1, jpi
       DO jj = 1, jpj
          jk = misfkt(ji,jj)
          !! Initialize arrays to 0 (each step)
          zt_sum = 0.e0_wp
          IF ( jk .GT. 1 ) THEN
    ! 1. -----------the average temperature between 200m and 600m ---------------------
             DO jk = misfkt(ji,jj),misfkb(ji,jj)
             ! freezing point temperature  at ice shelf base BG eq. 2 (JMM sign pb ??? +7.64e-4 !!!)
             ! after verif with UNESCO, wrong sign in BG eq. 2
             ! Calculate freezing temperature
                CALL eos_fzp(stbl(ji,jj), zt_frz, zpress) 
                zt_sum = zt_sum + (tsn(ji,jj,jk,jp_tem)-zt_frz) * fse3t(ji,jj,jk) * tmask(ji,jj,jk)  ! sum temp
             ENDDO
             zt_ave = zt_sum/rhisf_tbl(ji,jj) ! calcul mean value
    
    ! 2. ------------Net heat flux and fresh water flux due to the ice shelf
          ! For those corresponding to zonal boundary    
             qisf(ji,jj) = - rau0 * rcp * rn_gammat0 * risfLeff(ji,jj) * e1t(ji,jj) * zt_ave  &
                         & / (e1t(ji,jj) * e2t(ji,jj)) * tmask(ji,jj,jk) 
             
             fwfisf(ji,jj) = qisf(ji,jj) / lfusisf          !fresh water flux kg/(m2s)                  
             fwfisf(ji,jj) = fwfisf(ji,jj) * ( soce / stbl(ji,jj) )
             !add to salinity trend
          ELSE
             qisf(ji,jj) = 0._wp ; fwfisf(ji,jj) = 0._wp
          END IF
       ENDDO
    ENDDO
    !
    IF( nn_timing == 1 )  CALL timing_stop('sbc_isf_bg03')
  END SUBROUTINE sbc_isf_bg03

   SUBROUTINE sbc_isf_cav( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_isf_cav  ***
      !!
      !! ** Purpose :   handle surface boundary condition under ice shelf
      !!
      !! ** Method  : -
      !!
      !! ** Action  :   utau, vtau : remain unchanged
      !!                taum, wndm : remain unchanged
      !!                qns        : update heat flux below ice shelf
      !!                emp, emps  : update freshwater flux below ice shelf
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)          ::   kt         ! ocean time step
      !
      REAL(wp), DIMENSION(:,:), POINTER ::   zfrz, zsfrz
      REAL(wp), DIMENSION(:,:), POINTER ::   zgammat, zgammas 
      REAL(wp), DIMENSION(:,:), POINTER ::   zfwflx, zhtflx, zhtflx_b
      REAL(wp) ::   zlamb1, zlamb2, zlamb3
      REAL(wp) ::   zeps1,zeps2,zeps3,zeps4,zeps6,zeps7
      REAL(wp) ::   zaqe,zbqe,zcqe,zaqer,zdis,zcfac
      REAL(wp) ::   zeps = 1.e-20_wp        
      REAL(wp) ::   zerr
      INTEGER  ::   ji, jj     ! dummy loop indices
      INTEGER  ::   nit
      LOGICAL  ::   lit
      !!---------------------------------------------------------------------
      !
      ! coeficient for linearisation of potential tfreez
      ! Crude approximation for pressure (but commonly used)
      zlamb1=-0.0573_wp
      zlamb2= 0.0832_wp
      zlamb3=-7.53e-08 * grav * rau0
      IF( nn_timing == 1 )  CALL timing_start('sbc_isf_cav')
      !
      CALL wrk_alloc( jpi,jpj, zfrz  , zsfrz, zgammat, zgammas  )
      CALL wrk_alloc( jpi,jpj, zfwflx, zhtflx , zhtflx_b )

      ! initialisation
      zgammat(:,:)=rn_gammat0 ; zgammas (:,:)=rn_gammas0;
      zhtflx (:,:)=0.0_wp     ; zhtflx_b(:,:)=0.0_wp    ;
      zfwflx (:,:)=0.0_wp

      ! compute ice shelf melting
      nit = 1 ; lit = .TRUE.
      DO WHILE ( lit )    ! maybe just a constant number of iteration as in blk_core is fine
         IF (nn_isfblk == 1) THEN 
            ! ISOMIP formulation (2 equations) for volume flux (Hunter et al., 2006)
            ! Calculate freezing temperature
            CALL eos_fzp( stbl(:,:), zfrz(:,:), risfdep(:,:) )

            ! compute gammat every where (2d)
            CALL sbc_isf_gammats(zgammat, zgammas, zhtflx, zfwflx)
            
            ! compute upward heat flux zhtflx and upward water flux zwflx
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zhtflx(ji,jj) =   zgammat(ji,jj)*rcp*rau0*(ttbl(ji,jj)-zfrz(ji,jj))
                  zfwflx(ji,jj) = - zhtflx(ji,jj)/lfusisf
               END DO
            END DO

            ! Compute heat flux and upward fresh water flux
            qisf  (:,:) = - zhtflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)
            fwfisf(:,:) =   zfwflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)

         ELSE IF (nn_isfblk == 2 ) THEN
            ! ISOMIP+ formulation (3 equations) for volume flux (Asay-Davis et al., 2015) 
            ! compute gammat every where (2d)
            CALL sbc_isf_gammats(zgammat, zgammas, zhtflx, zfwflx)

            ! compute upward heat flux zhtflx and upward water flux zwflx
            ! Resolution of a 2d equation from equation 21, 22 and 23 to find Sb (Asay-Davis et al., 2015)
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! compute coeficient to solve the 2nd order equation
                  zeps1=rcp*rau0*zgammat(ji,jj)
                  zeps2=lfusisf*rau0*zgammas(ji,jj)
                  zeps3=rhoisf*rcpi*kappa/MAX(risfdep(ji,jj),zeps)
                  zeps4=zlamb2+zlamb3*risfdep(ji,jj)
                  zeps6=zeps4-ttbl(ji,jj)
                  zeps7=zeps4-tsurf
                  zaqe=zlamb1 * (zeps1 + zeps3)
                  zaqer=0.5/MIN(zaqe,-zeps)
                  zbqe=zeps1*zeps6+zeps3*zeps7-zeps2
                  zcqe=zeps2*stbl(ji,jj)
                  zdis=zbqe*zbqe-4.0*zaqe*zcqe               

                  ! Presumably zdis can never be negative because gammas is very small compared to gammat
                  ! compute s freeze
                  zsfrz(ji,jj)=(-zbqe-SQRT(zdis))*zaqer
                  IF ( zsfrz(ji,jj) .LT. 0.0_wp ) zsfrz(ji,jj)=(-zbqe+SQRT(zdis))*zaqer

                  ! compute t freeze (eq. 22)
                  zfrz(ji,jj)=zeps4+zlamb1*zsfrz(ji,jj)
  
                  ! zfwflx is upward water flux
                  ! zhtflx is upward heat flux (out of ocean)
                  ! compute the upward water and heat flux (eq. 28 and eq. 29)
                  zfwflx(ji,jj) = rau0 * zgammas(ji,jj) * (zsfrz(ji,jj)-stbl(ji,jj)) / MAX(zsfrz(ji,jj),zeps)
                  zhtflx(ji,jj) = zgammat(ji,jj) * rau0 * rcp * (ttbl(ji,jj) - zfrz(ji,jj) ) 
               END DO
            END DO

            ! compute heat and water flux
            qisf  (:,:) = - zhtflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)
            fwfisf(:,:) =   zfwflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)

         ENDIF

         ! define if we need to iterate (nn_gammablk 0/1 do not need iteration)
         IF ( nn_gammablk .LT. 2 ) THEN ; lit = .FALSE.
         ELSE                           
            ! check total number of iteration
            IF (nit .GE. 100) THEN ; CALL ctl_stop( 'STOP', 'sbc_isf_hol99 : too many iteration ...' )
            ELSE                   ; nit = nit + 1
            ENDIF

            ! compute error between 2 iterations
            ! if needed save gammat and compute zhtflx_b for next iteration
            zerr = MAXVAL(ABS(zhtflx-zhtflx_b))
            IF ( zerr .LE. 0.01 ) THEN ; lit = .FALSE.
            ELSE                       ; zhtflx_b(:,:) = zhtflx(:,:)
            ENDIF
         END IF
      END DO
      !
      ! output
      IF( iom_use('isfthermdr') ) CALL iom_put('isfthermdr', ttbl-zfrz)
      IF( iom_use('isfhalindr') ) CALL iom_put('isfhalindr', zsfrz-stbl)
      IF( iom_use('isfgammat') ) CALL iom_put('isfgammat', zgammat)
      IF( iom_use('isfgammas') ) CALL iom_put('isfgammas', zgammas)
      ! 
      CALL wrk_dealloc( jpi,jpj, zfrz  , zsfrz, zgammat, zgammas  )
      CALL wrk_dealloc( jpi,jpj, zfwflx, zhtflx , zhtflx_b )
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_isf_cav')

   END SUBROUTINE sbc_isf_cav

   SUBROUTINE sbc_isf_gammats(gt, gs, zqhisf, zqwisf )
      !!----------------------------------------------------------------------
      !! ** Purpose    : compute the coefficient echange for heat flux
      !!
      !! ** Method     : gamma assume constant or depends of u* and stability
      !!
      !! ** References : Holland and Jenkins, 1999, JPO, p1787-1800, eq 14
      !!                Jenkins et al., 2010, JPO, p2298-2312
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) :: gt, gs
      REAL(wp), DIMENSION(:,:), INTENT(in ) :: zqhisf, zqwisf

      INTEGER  :: ikt                        
      INTEGER  :: ji,jj                      ! loop index
      REAL(wp), DIMENSION(:,:), POINTER :: zustar           ! U, V at T point and friction velocity
      REAL(wp) :: zdku, zdkv                 ! U, V shear 
      REAL(wp) :: zPr, zSc, zRc              ! Prandtl, Scmidth and Richardson number 
      REAL(wp) :: zmob, zmols                ! Monin Obukov length, coriolis factor at T point
      REAL(wp) :: zbuofdep, zhnu             ! Bouyancy length scale, sublayer tickness
      REAL(wp) :: zhmax                      ! limitation of mol
      REAL(wp) :: zetastar                   ! stability parameter
      REAL(wp) :: zgmolet, zgmoles, zgturb   ! contribution of modelecular sublayer and turbulence 
      REAL(wp) :: zcoef                      ! temporary coef
      REAL(wp) :: zdep
      REAL(wp) :: zeps = 1.0e-20_wp    
      REAL(wp), PARAMETER :: zxsiN = 0.052   ! dimensionless constant
      REAL(wp), PARAMETER :: znu   = 1.95e-6 ! kinamatic viscosity of sea water (m2.s-1)
      REAL(wp), DIMENSION(2) :: zts, zab
      !!---------------------------------------------------------------------
      CALL wrk_alloc( jpi,jpj, zustar )
      !
      IF      ( nn_gammablk == 0 ) THEN
      !! gamma is constant (specified in namelist)
      !! ISOMIP formulation (Hunter et al, 2006)
         gt(:,:) = rn_gammat0
         gs(:,:) = rn_gammas0

      ELSE IF ( nn_gammablk == 1 ) THEN
      !! gamma is assume to be proportional to u* 
      !! Jenkins et al., 2010, JPO, p2298-2312
      !! Adopted by Asay-Davis et al. (2015)

      !! compute ustar (eq. 24)
         zustar(:,:) = SQRT( rn_tfri2 * (utbl(:,:) * utbl(:,:) + vtbl(:,:) * vtbl(:,:) + rn_tfeb2) )

      !! Compute gammats
         gt(:,:) = zustar(:,:) * rn_gammat0
         gs(:,:) = zustar(:,:) * rn_gammas0
      
      ELSE IF ( nn_gammablk == 2 ) THEN
      !! gamma depends of stability of boundary layer
      !! Holland and Jenkins, 1999, JPO, p1787-1800, eq 14
      !! as MOL depends of flux and flux depends of MOL, best will be iteration (TO DO)
      !! compute ustar
         zustar(:,:) = SQRT( rn_tfri2 * (utbl(:,:) * utbl(:,:) + vtbl(:,:) * vtbl(:,:) + rn_tfeb2) )

      !! compute Pr and Sc number (can be improved)
         zPr =   13.8
         zSc = 2432.0

      !! compute gamma mole
         zgmolet = 12.5 * zPr ** (2.0/3.0) - 6.0
         zgmoles = 12.5 * zSc ** (2.0/3.0) -6.0

      !! compute gamma
         DO ji=2,jpi
            DO jj=2,jpj
               ikt = mikt(ji,jj)

               IF (zustar(ji,jj) == 0._wp) THEN           ! only for kt = 1 I think
                  gt = rn_gammat0
                  gs = rn_gammas0
               ELSE
      !! compute Rc number (as done in zdfric.F90)
                  zcoef = 0.5 / fse3w(ji,jj,ikt)
                  !                                            ! shear of horizontal velocity
                  zdku = zcoef * (  un(ji-1,jj  ,ikt  ) + un(ji,jj,ikt  )  &
                     &             -un(ji-1,jj  ,ikt+1) - un(ji,jj,ikt+1)  )
                  zdkv = zcoef * (  vn(ji  ,jj-1,ikt  ) + vn(ji,jj,ikt  )  &
                     &             -vn(ji  ,jj-1,ikt+1) - vn(ji,jj,ikt+1)  )
                  !                                            ! richardson number (minimum value set to zero)
                  zRc = rn2(ji,jj,ikt+1) / MAX( zdku*zdku + zdkv*zdkv, zeps )

      !! compute bouyancy 
                  zts(jp_tem) = ttbl(ji,jj)
                  zts(jp_sal) = stbl(ji,jj)
                  zdep        = fsdepw(ji,jj,ikt)
                  !
                  CALL eos_rab( zts, zdep, zab )
                  !
      !! compute length scale 
                  zbuofdep = grav * ( zab(jp_tem) * zqhisf(ji,jj) - zab(jp_sal) * zqwisf(ji,jj) )  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! compute Monin Obukov Length
                  ! Maximum boundary layer depth
                  zhmax = fsdept(ji,jj,mbkt(ji,jj)) - fsdepw(ji,jj,mikt(ji,jj)) - 0.001_wp
                  ! Compute Monin obukhov length scale at the surface and Ekman depth:
                  zmob   = zustar(ji,jj) ** 3 / (vkarmn * (zbuofdep + zeps))
                  zmols  = SIGN(1._wp, zmob) * MIN(ABS(zmob), zhmax) * tmask(ji,jj,ikt)

      !! compute eta* (stability parameter)
                  zetastar = 1._wp / ( SQRT(1._wp + MAX(zxsiN * zustar(ji,jj) / ( ABS(ff(ji,jj)) * zmols * zRc ), 0.0_wp)))

      !! compute the sublayer thickness
                  zhnu = 5 * znu / zustar(ji,jj)

      !! compute gamma turb
                  zgturb = 1._wp / vkarmn * LOG(zustar(ji,jj) * zxsiN * zetastar * zetastar / ( ABS(ff(ji,jj)) * zhnu )) &
                  &      + 1._wp / ( 2 * zxsiN * zetastar ) - 1._wp / vkarmn

      !! compute gammats
                  gt(ji,jj) = zustar(ji,jj) / (zgturb + zgmolet)
                  gs(ji,jj) = zustar(ji,jj) / (zgturb + zgmoles)
               END IF
            END DO
         END DO
         CALL lbc_lnk(gt(:,:),'T',1.)
         CALL lbc_lnk(gs(:,:),'T',1.)
      END IF
      CALL wrk_dealloc( jpi,jpj, zustar )

   END SUBROUTINE

   SUBROUTINE sbc_isf_tbl( varin, varout, cptin )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE sbc_isf_tbl  ***
      !!
      !! ** Purpose : compute mean T/S/U/V in the boundary layer at T- point
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in) :: varin
      REAL(wp), DIMENSION(:,:)  , INTENT(out):: varout
      
      CHARACTER(len=1), INTENT(in) :: cptin ! point of variable in/out

      REAL(wp) :: ze3, zhk
      REAL(wp), DIMENSION(:,:), POINTER :: zhisf_tbl ! thickness of the tbl

      INTEGER :: ji,jj,jk                   ! loop index
      INTEGER :: ikt,ikb                    ! top and bottom index of the tbl
 
      ! allocation
      CALL wrk_alloc( jpi,jpj, zhisf_tbl)
      
      ! initialisation
      varout(:,:)=0._wp
   
      ! compute U in the top boundary layer at T- point 
      IF (cptin == 'U') THEN
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = miku(ji,jj) ; ikb = miku(ji,jj)
            ! thickness of boundary layer at least the top level thickness
               zhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), fse3u_n(ji,jj,ikt))

            ! determine the deepest level influenced by the boundary layer
               DO jk = ikt+1, mbku(ji,jj)
                  IF ( (SUM(fse3u_n(ji,jj,ikt:jk-1)) .LT. zhisf_tbl(ji,jj)) .AND. (umask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               zhisf_tbl(ji,jj) = MIN(zhisf_tbl(ji,jj), SUM(fse3u_n(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.

            ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = fse3u_n(ji,jj,jk)
                  varout(ji,jj) = varout(ji,jj) + varin(ji,jj,jk) / zhisf_tbl(ji,jj) * ze3
               END DO

            ! level partially include in ice shelf boundary layer 
               zhk = SUM( fse3u_n(ji, jj, ikt:ikb - 1)) / zhisf_tbl(ji,jj)
               varout(ji,jj) = varout(ji,jj) + varin(ji,jj,ikb) * (1._wp - zhk)
            END DO
         END DO
         DO jj = 2,jpj
            DO ji = 2,jpi
               varout(ji,jj) = 0.5_wp * (varout(ji,jj) + varout(ji-1,jj))
            END DO
         END DO
         CALL lbc_lnk(varout,'T',-1.)
      END IF

      ! compute V in the top boundary layer at T- point 
      IF (cptin == 'V') THEN
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = mikv(ji,jj) ; ikb = mikv(ji,jj)
           ! thickness of boundary layer at least the top level thickness
               zhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), fse3v_n(ji,jj,ikt))

            ! determine the deepest level influenced by the boundary layer
               DO jk = ikt+1, mbkv(ji,jj)
                  IF ( (SUM(fse3v_n(ji,jj,ikt:jk-1)) .LT. zhisf_tbl(ji,jj)) .AND. (vmask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               zhisf_tbl(ji,jj) = MIN(zhisf_tbl(ji,jj), SUM(fse3v_n(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.

            ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = fse3v_n(ji,jj,jk)
                  varout(ji,jj) = varout(ji,jj) + varin(ji,jj,jk) / zhisf_tbl(ji,jj) * ze3
               END DO

            ! level partially include in ice shelf boundary layer 
               zhk = SUM( fse3v_n(ji, jj, ikt:ikb - 1)) / zhisf_tbl(ji,jj)
               varout(ji,jj) = varout(ji,jj) + varin(ji,jj,ikb) * (1._wp - zhk)
            END DO
         END DO
         DO jj = 2,jpj
            DO ji = 2,jpi
               varout(ji,jj) = 0.5_wp * (varout(ji,jj) + varout(ji,jj-1))
            END DO
         END DO
         CALL lbc_lnk(varout,'T',-1.)
      END IF

      ! compute T in the top boundary layer at T- point 
      IF (cptin == 'T') THEN
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkb(ji,jj)

            ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = fse3t_n(ji,jj,jk)
                  varout(ji,jj) = varout(ji,jj) + varin(ji,jj,jk) * r1_hisf_tbl(ji,jj) * ze3
               END DO

            ! level partially include in ice shelf boundary layer 
               zhk = SUM( fse3t_n(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj)
               varout(ji,jj) = varout(ji,jj) + varin(ji,jj,ikb) * (1._wp - zhk)
            END DO
         END DO
      END IF

      ! mask mean tbl value
      varout(:,:) = varout(:,:) * ssmask(:,:)

      ! deallocation
      CALL wrk_dealloc( jpi,jpj, zhisf_tbl )      

   END SUBROUTINE sbc_isf_tbl
      

   SUBROUTINE sbc_isf_div( phdivn )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE sbc_isf_div  ***
      !!       
      !! ** Purpose :   update the horizontal divergence with the runoff inflow
      !!
      !! ** Method  :   
      !!                CAUTION : risf_tsc(:,:,jp_sal) is negative (outflow) increase the 
      !!                          divergence and expressed in m/s
      !!
      !! ** Action  :   phdivn   decreased by the runoff inflow
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   phdivn   ! horizontal divergence
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikt, ikb 
      REAL(wp) ::   zhk
      REAL(wp) ::   zfact     ! local scalar
      !!----------------------------------------------------------------------
      !
      zfact   = 0.5_wp
      !
      IF (lk_vvl) THEN     ! need to re compute level distribution of isf fresh water
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkt(ji,jj)
               ! thickness of boundary layer at least the top level thickness
               rhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), fse3t(ji,jj,ikt))

               ! determine the deepest level influenced by the boundary layer
               DO jk = ikt+1, mbkt(ji,jj)
                  IF ( (SUM(fse3t(ji,jj,ikt:jk-1)) .LT. rhisf_tbl(ji,jj)) .AND. (tmask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               rhisf_tbl(ji,jj) = MIN(rhisf_tbl(ji,jj), SUM(fse3t(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.
               misfkb(ji,jj) = ikb                                                  ! last wet level of the tbl
               r1_hisf_tbl(ji,jj) = 1._wp / rhisf_tbl(ji,jj)

               zhk           = SUM( fse3t(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj) ! proportion of tbl cover by cell from ikt to ikb - 1
               ralpha(ji,jj) = rhisf_tbl(ji,jj) * (1._wp - zhk ) / fse3t(ji,jj,ikb)  ! proportion of bottom cell influenced by boundary layer
            END DO
         END DO
      END IF 
      !
      !==   ice shelf melting distributed over several levels   ==!
      DO jj = 1,jpj
         DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkb(ji,jj)
               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  phdivn(ji,jj,jk) = phdivn(ji,jj,jk) + ( fwfisf(ji,jj) + fwfisf_b(ji,jj) ) &
                    &               * r1_hisf_tbl(ji,jj) * r1_rau0 * zfact
               END DO
               ! level partially include in ice shelf boundary layer 
               phdivn(ji,jj,ikb) = phdivn(ji,jj,ikb) + ( fwfisf(ji,jj) &
                  &             + fwfisf_b(ji,jj) ) * r1_hisf_tbl(ji,jj) * r1_rau0 * zfact * ralpha(ji,jj) 
         END DO
      END DO
      !
   END SUBROUTINE sbc_isf_div
   !!======================================================================
END MODULE sbcisf
