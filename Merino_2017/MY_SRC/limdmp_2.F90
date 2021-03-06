MODULE limdmp_2
   !!======================================================================
   !!                       ***  MODULE limdmp_2   ***
   !!  LIM-2 ice model : restoring Ice thickness and Fraction leads
   !!======================================================================
   !! History :   2.0  !  2004-04 (S. Theetten) Original code
   !!             3.3  !  2010-06 (J.-M. Molines) use of fldread
   !!----------------------------------------------------------------------
#if defined key_lim2
   !!----------------------------------------------------------------------
   !!   'key_lim2'                                    LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_dmp_2     : ice model damping
   !!----------------------------------------------------------------------
   USE ice_2          ! ice variables 
   USE sbc_oce, ONLY : nn_fsbc ! for fldread
   USE dom_oce         ! for mi0; mi1 etc ...
   USE fldread         ! read input fields
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE lbclnk
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_dmp_2     ! called by sbc_ice_lim2

   INTEGER  , PARAMETER :: jp_hicif = 1 , jp_frld = 2
   REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:) ::   resto_ice   ! restoring coeff. on ICE   [s-1]
   TYPE(FLD), ALLOCATABLE, DIMENSION(:)     ::   sf_icedmp   ! structure of ice damping input
   
   !! * Substitution
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/LIM 3.3 , UCL-NEMO-consortium (2010) 
   !! $Id: limdmp_2.F90 3551 2012-11-14 11:00:10Z gm $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_dmp_2( kt )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_dmp_2  ***
      !!
      !! ** purpose :   restore ice thickness and lead fraction
      !!
      !! ** method  :   restore ice thickness and lead fraction using a restoring
      !!              coefficient defined by the user in lim_dmp_init
      !!
      !! ** Action  : - update hicif and frld  
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::   zfrld, zhice  ! local scalars
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zrestoice
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'lim_dmp_2 : Ice thickness and ice concentration restoring'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
         !
         ! ice_resto_init create resto_ice (in 1/s) for restoring ice parameters near open boundaries.
         ! Double check this routine to verify if it corresponds to your config
         CALL lim_dmp_init
      ENDIF
      !
      IF( ln_limdmp ) THEN   ! ice restoring in this case
         !
         CALL fld_read( kt, nn_fsbc, sf_icedmp )
         !
!CDIR COLLAPSE
         hicif(:,:) = MAX( 0._wp,                     &        ! h >= 0         avoid spurious out of physical range
            &         hicif(:,:) - rdt_ice * resto_ice(:,:,1) * ( hicif(:,:) -    sf_icedmp(jp_hicif)%fnow(:,:,1) )     ) 
!CDIR COLLAPSE
         frld (:,:) = MAX( 0._wp, MIN( 1._wp,         &        ! 0<= frld<=1    values which blow the run up
            &         frld (:,:) - rdt_ice * resto_ice(:,:,1) * ( frld (:,:) - (1._wp - sf_icedmp(jp_frld )%fnow(:,:,1)))  )  )
         ALLOCATE(zrestoice(jpi, jpj))
         zrestoice(:,:)=rdt_ice * resto_ice(:,:,1)
         CALL iom_put('iceresto',zrestoice)
         zrestoice(:,:)=sf_icedmp(jp_hicif)%fnow(:,:,1)
         CALL iom_put('restoih',zrestoice)
         zrestoice(:,:)=1-sf_icedmp(jp_frld)%fnow(:,:,1)
         CALL iom_put('restoic',zrestoice)
         !CALL iom_put('iceresto',frld)
         DEALLOCATE(zrestoice)
         !
      ENDIF
      !
   END SUBROUTINE lim_dmp_2


   SUBROUTINE lim_dmp_init
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE lim_dmp_init  ***
      !!
      !! ** Purpose :   set the coefficient for the ice thickness and lead fraction restoring
      !!
      !! ** Method  :   restoring is used to mimic ice open boundaries.
      !!              the restoring coef. (a 2D array) has to be defined by the user.
      !!              here is given as an example a restoring along north and south boundaries
      !!      
      !! ** Action  :   define resto_ice(:,:,1)
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jj, jk       ! dummy loop indices
      INTEGER  :: irelax, ierror   ! error flag for allocation
      !
      REAL(wp) :: zdmpmax, zdmpmin, zfactor, zreltim ! temporary scalar
      !
      CHARACTER(len=100)           ::   cn_dir       ! Root directory for location of ssr files
      TYPE(FLD_N), DIMENSION (2)   ::   sl_icedmp    ! informations about the icedmp  field to be read
      TYPE(FLD_N)                  ::   sn_hicif     ! 
      TYPE(FLD_N)                  ::   sn_frld      ! 
      NAMELIST/namice_dmp/ cn_dir, sn_hicif, sn_frld
      !!----------------------------------------------------------------------
      !
      ! 1)  initialize fld read structure for input data 
      !     --------------------------------------------
      cn_dir    = './'
      ! (NB: frequency positive => hours, negative => months)
      !                !    file     ! frequency ! variable ! time intep !  clim   ! 'yearly' or ! weights  ! rotation !
      !                !    name     !  (hours)  !  name    !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs    !
      sn_hicif = FLD_N( 'ice_damping ', -1       , 'hicif'  ,  .true.    , .true.  ,  'yearly'   ,  ''      ,  ''      )
      sn_frld  = FLD_N( 'ice_damping ', -1       , 'frld'   ,  .true.    , .true.  ,  'yearly'   ,  ''      ,  ''      )

      REWIND( numnam_ice )                !* read in namelist_ice namicedmp
      READ  ( numnam_ice, namice_dmp )
      !
      IF ( lwp ) THEN                     !* control print
         WRITE (numout,*)'     lim_dmp_init : lim_dmp initialization ' 
         WRITE (numout,*)'       Namelist namicedmp read '
         WRITE (numout,*)'         Ice restoring (T) or not (F) ln_limdmp =', ln_limdmp 
         WRITE (numout,*)
         WRITE (numout,*)'     CAUTION : here hard coded ice restoring along northern and southern boundaries'
         WRITE (numout,*)'               adapt the lim_dmp_init routine to your needs'
      ENDIF

      ! 2)  initialise resto_ice    ==>  config dependant !
      !     --------------------         ++++++++++++++++
      !
      IF( ln_limdmp ) THEN                !* ice restoring is used, follow initialization
         ! 
         sl_icedmp ( jp_hicif ) = sn_hicif
         sl_icedmp ( jp_frld  ) = sn_frld
         ALLOCATE ( sf_icedmp (2) , resto_ice(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'lim_dmp_init: unable to allocate sf_icedmp structure or resto_ice array' )   ;   RETURN
         ENDIF
         ALLOCATE( sf_icedmp(jp_hicif)%fnow(jpi,jpj,1) , sf_icedmp(jp_hicif)%fdta(jpi,jpj,1,2) )
         ALLOCATE( sf_icedmp(jp_frld )%fnow(jpi,jpj,1) , sf_icedmp(jp_frld )%fdta(jpi,jpj,1,2) )
         !                         ! fill sf_icedmp with sn_icedmp and control print
         CALL fld_fill( sf_icedmp, sl_icedmp, cn_dir, 'lim_dmp_init', 'Ice  restoring input data', 'namicedmp' )
      
         resto_ice(:,:,:) = 0._wp
         !      Re-calculate the North and South boundary restoring term
         !      because those boundaries may change with the prescribed zoom area.
         !
         irelax  = 10                     ! width of buffer zone with respect to close boundary
         zdmpmax = 10._wp                  ! max restoring time scale  (days) (low restoring)
         zdmpmin = rdt_ice / 86400._wp    ! min restoring time scale  (days) (high restoring)
         !                                ! days / grid-point
         zfactor = ( zdmpmax - zdmpmin ) / REAL( irelax, wp )

         !    South boundary restoring term
         ! REM: if there is no ice in the model and in the data, 
         !      no restoring even with non zero resto_ice
         !DO jj = mj0((jpjzoom+1) - 1 + 1), mj1((jpjzoom+1) -1 + irelax)
         !   zreltim = zdmpmin + zfactor * ( mjg(jj) - (jpjzoom+1) )
         !   resto_ice(:,jj,:) = 1._wp / ( zreltim * 86400._wp )
         !END DO

         ! North boundary restoring term
         resto_ice(:,mj0(jpjglo-1):mj1(jpjglo-1),1)=1._wp / ( zdmpmin * 86400 )
         DO jj =  mj0((jpjzoom-2) -1 + jpjglo - irelax), mj1((jpjzoom-2) - 1 + jpjglo)
            zreltim = zdmpmin + zfactor * (jpjglo - ( mjg(jj) - (jpjzoom-2) +1))
            resto_ice(:,jj,:) = 1._wp / ( zreltim * 86400 )
         END DO
  
         ! East boundary restoring term
         resto_ice(mi0(jpiglo-1):mi1(jpiglo-1),:,1)=1._wp / ( zdmpmin * 86400 )
         DO jj = 1,jpj
            DO ji =  mi0((jpizoom-2) -1 + jpiglo - irelax), mi1((jpizoom-2) - 1 + jpiglo)
               zreltim = zdmpmin + zfactor * (jpiglo - ( mig(ji) - (jpizoom-2) +1 ))
               resto_ice(ji,jj,1) = MAX(1._wp / ( zreltim * 86400 ), resto_ice(ji,jj,1))
            END DO

         ! West boundary restoring term
            DO ji = mi0((jpizoom+2) -1 + 1), mi1((jpizoom+2) - 1 + irelax)
               zreltim = zdmpmin + zfactor * ( mig(ji) - (jpizoom+2) )
               resto_ice(ji,jj,1) = MAX(1._wp / ( zreltim * 86400 ), resto_ice(ji,jj,1))
            END DO
         END DO
         resto_ice(mi0(jpizoom+1):mi1(jpizoom+1),:,1)=1._wp / ( zdmpmin * 86400 )

      ENDIF
      !
   END SUBROUTINE lim_dmp_init
   
#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty Module                  No ice damping
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_dmp_2( kt )        ! Dummy routine
      WRITE(*,*) 'lim_dmp_2: You should not see this print! error? ', kt
   END SUBROUTINE lim_dmp_2
#endif

   !!======================================================================
END MODULE limdmp_2
