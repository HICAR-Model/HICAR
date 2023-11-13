!>----------------------------------------------------------
!! This module provides a wrapper to call various PBL models
!! It sets up variables specific to the physics package to be used including both
!!
!! The main entry point to the code is pbl(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!!  pbl_init->[ external initialization routines]
!!  pbl->[  external PBL routines]
!!  pbl_finalize
!!
!! High level routine descriptions / purpose
!!   pbl_init           - initializes physics package
!!   pbl                - sets up and calls main physics package
!!   pbl_finalize       - permits physics package cleanup (close files, deallocate memory)
!!
!! Inputs: domain, options, dt
!!      domain,options  = as defined in data_structures
!!      dt              = time step (seconds)
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module planetary_boundary_layer
    use data_structures
    use domain_interface,   only : domain_t
    use options_interface,  only : options_t
    use pbl_simple,    only : simple_pbl, finalize_simple_pbl, init_simple_pbl
    use pbl_diagnostic, only : diagnostic_pbl, finalize_diagnostic_pbl, init_diagnostic_pbl
    !use module_bl_ysu, only : ysuinit, ysu
    use module_bl_ysu, only : ysuinit, ysu
    use mod_wrf_constants, only : EOMEG, XLV, r_v, R_d, KARMAN, gravity, EP_1, EP_2, cp, rcp, rovg
    use icar_constants !, only : karman,stefan_boltzmann
    use mod_pbl_utilities, only : da_sfc_wtq
    use ieee_arithmetic ! for debugging
    use array_utilities, only : array_offset_x_3d, array_offset_y_3d


    implicit none
    real,allocatable, dimension(:,:)    ::  windspd, hpbl, psim, &
                                            psih, u10d, v10d, CHS, xland_real, regime
    ! integer, allocatable, dimension(:,:) :: kpbl2d
    real, allocatable, dimension(:,:,:) :: tend_u_ugrid, tend_v_vgrid, RTHRATEN

    private
    public :: pbl_var_request, pbl_init, pbl, pbl_finalize

    integer :: ids, ide, jds, jde, kds, kde,  &
               ims, ime, jms, jme, kms, kme,  &
               its, ite, jts, jte, kts, kte, j, k, i

    logical :: allowed_to_read, restart, flag_qi

contains

    subroutine pbl_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%landsurface == kPBL_SIMPLE) then
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, &
                         kVARS%cloud_water, kVARS%cloud_ice,              &
                         kVARS%rain_in_air, kVARS%snow_in_air,            &
                         kVARS%exner, kVARS%dz_interface, kVARS%density,  &
                         kVARS%u, kVARS%v, kVARS%land_mask])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, &
                         kVARS%exner, kVARS%dz_interface, kVARS%density,  &
                         kVARS%u, kVARS%v, kVARS%land_mask])
        endif
        if (options%physics%boundarylayer==kPBL_YSU) then

            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface,          &
                         kVARS%skin_temperature, kVARS%terrain, kVARS%ground_surf_temperature,              &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,                  &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%ground_heat_flux,                 &
                         kVARS%roughness_z0, kVARS%ustar, kVARS%cloud_ice,                                  &
                         kVARS%tend_th_pbl, kVARS%tend_qc_pbl, kVARS%tend_qi_pbl,  kVARS%temperature_2m,    &
                         kVARS%tend_u, kVARS%tend_v, kVARS%tend_qv_pbl, kVARS%pressure, kVARS%kpbl,         &
                         kVARS%fm, kVARS%fh, kVARS%QFX, kVARS%br,                                          &
                         kVARS%land_mask, kVARS%cloud_water, kVARS%coeff_heat_exchange_3d, kVARS%coeff_momentum_exchange_3d, kVARS%hpbl ]) !kVARS%tend_qv_adv,kVARS%tend_qv, kVARS%tend_qs, kVARS%tend_qr,, kVARS%u_mass, kVARS%v_mass,
!           kVARS%coeff_momentum_drag, ??
             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_ice, kVARS%cloud_water]) !??

             call options%restart_vars( &
                        [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                &
                        kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface,          &
                        kVARS%skin_temperature, kVARS%terrain, kVARS%ground_surf_temperature,              &
                        kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,                  &
                        kVARS%humidity_2m, kVARS%surface_pressure, kVARS%ground_heat_flux,                 &
                        kVARS%roughness_z0, kVARS%cloud_ice, kVARS%QFX,       &
                        kVARS%temperature_2m,    &
                        kVARS%pressure,         &
                        kVARS%cloud_water,kVARS%coeff_heat_exchange_3d, kVARS%coeff_momentum_exchange_3d, kVARS%hpbl  ]) !kVARS%u_mass, kVARS%v_mass,
        endif
    end subroutine pbl_var_request


    subroutine pbl_init(domain,options)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(options_t),    intent(in)      :: options

        ids = domain%ids ; ide = domain%ide ; jds = domain%jds ; jde = domain%jde ; kds = domain%kds ; kde = domain%kde
        ims = domain%ims ; ime = domain%ime ; jms = domain%jms ; jme = domain%jme ; kms = domain%kms ; kme = domain%kme
        its = domain%its ; ite = domain%ite ; jts = domain%jts ; jte = domain%jte ; kts = domain%kts ; kte = domain%kte

        allowed_to_read = .True.
        restart = .False.
        flag_qi = .true.
        if (.not.allocated(domain%tend%qv_pbl)) allocate(domain%tend%qv_pbl(ims:ime,kms:kme,jms:jme))
        domain%tend%qv_pbl=0

        if (this_image()==1) write(*,*) "Initializing PBL Scheme"

        if (options%physics%boundarylayer==kPBL_SIMPLE) then
            if (this_image()==1) write(*,*) "    Simple PBL"
            call init_simple_pbl(domain, options)
        else if (options%physics%boundarylayer==kPBL_DIAGNOSTIC) then
            if (this_image()==1) write(*,*) "    Diagnostic PBL"
            call init_diagnostic_pbl(domain, options)
        else if (options%physics%boundarylayer==kPBL_YSU) then

            if (this_image()==1) write(*,*) "    YSU PBL"

            ! allocate local vars YSU:
            allocate(windspd(ims:ime, jms:jme))
            ! allocate(hpbl(ims:ime, jms:jme))  ! this should go to domain object for convective modules!!
            allocate(psim(ims:ime, jms:jme))
            ! psim= 0.5
            allocate(psih(ims:ime, jms:jme))
            ! psih=0.5
            allocate(u10d(ims:ime, jms:jme))
            allocate(v10d(ims:ime, jms:jme))
            ! allocate(kpbl2d(ims:ime, jms:jme)) ! domain%kpbl now
            ! allocate(CHS(ims:ime,jms:jme))
            ! CHS = 0.01
            allocate(xland_real(ims:ime,jms:jme))
            xland_real=real(domain%land_mask)
            allocate(regime(ims:ime,jms:jme))
            allocate(tend_u_ugrid(ims:ime+1, kms:kme, jms:jme)) ! to add the calculated u/v tendencies to the u/v grid
            allocate(tend_v_vgrid(ims:ime, kms:kme, jms:jme+1))
            allocate(RTHRATEN(ims:ime, kms:kme, jms:jme)) !initialize radiative heating tendencies and set to 0 in case user turns on ysu radiative heating w/o radiations scheme
            RTHRATEN = 0.0
            ! initialize tendencies (this is done in ysu init but only for tiles, not mem (ie its vs ims))
            ! BK: check if this actually matters ???
            if(.not.restart)then
                do j = jms,jme
                do k = kms,kme
                do i = ims,ime
                    domain%tend%u(i,k,j) = 0.
                    domain%tend%v(i,k,j) = 0.
                    domain%tend%th_pbl(i,k,j) = 0.
                    domain%tend%qv_pbl(i,k,j) = 0.
                    domain%tend%qc_pbl(i,k,j) = 0.
                    domain%tend%qi_pbl(i,k,j) = 0.
                enddo
                enddo
                enddo
              endif


            call ysuinit(rublten=domain%tend%u                  &
                        ,rvblten=domain%tend%v                  &
                        ,rthblten=domain%tend%th_pbl            &
                        ,rqvblten=domain%tend%qv_pbl            &
                        ,rqcblten=domain%tend%qc_pbl            &
                        ,rqiblten=domain%tend%qi_pbl            &
                        ,p_qi=1                                 &
                        ,p_first_scalar=1                       &
                        ,restart=restart                        &
                        ,allowed_to_read= allowed_to_read      &
                        ,ids=ids, ide=ide, jds=jds, jde=jde     &
                        ,kds=kds, kde=kde, ims=ims, ime=ime     &
                        ,jms=jms, jme=jme, kms=kms, kme=kme     &
                        ,its=its, ite=ite, jts=jts, jte=jte     &
                        ,kts=kts, kte=kte-1)
        endif
    end subroutine pbl_init

    subroutine pbl(domain, options, dt_in)
        implicit none
        type(domain_t),  intent(inout)  :: domain
        type(options_t), intent(in)     :: options
        real,            intent(in)     :: dt_in  !  =real(dt%seconds())

        if (options%physics%boundarylayer==kPBL_SIMPLE) then
            call simple_pbl(domain% potential_temperature %data_3d,     &
                            domain% water_vapor           %data_3d,     &
                            domain% cloud_water_mass      %data_3d,     &
                            domain% cloud_ice_mass        %data_3d,     &
                            domain% rain_mass             %data_3d,     &
                            domain% snow_mass             %data_3d,     &
                            domain% u_mass                %data_3d,     &
                            domain% v_mass                %data_3d,     &
                            domain% exner                 %data_3d,     &
                            domain% density               %data_3d,     &
                            domain% z                     %data_3d,     &
                            domain% dz_mass               %data_3d,     &
                            domain% terrain               %data_2d,     &
                            domain% land_mask,                          &
                            its, ite, jts, jte, kts, kte,               &
                            dt_in)
                            ! domain% qv_pbl_tendency     %data_3d)
        endif
        if (options%physics%boundarylayer==kPBL_DIAGNOSTIC) then
            call diagnostic_pbl(domain,dt_in)
                            ! domain% qv_pbl_tendency     %data_3d)
        endif
        if (options%physics%boundarylayer==kPBL_YSU) then

            ! windspd=sqrt(  domain%u_mass%data_3d(ims:ime, 1, jms:jme)**2 +     &
            !             domain%v_mass%data_3d(ims:ime, 1, jms:jme)**2   )
            windspd = sqrt(domain%u_10m%data_2d**2 + domain%v_10m%data_2d**2) ! as it is done in lsm_driver.
            where(windspd==0) windspd=1e-5

            if (options%physics%radiation==kRA_RRTMG) then
                RTHRATEN = domain%tend%th_lwrad + domain%tend%th_swrad
            endif

            call ysu(u3d=domain%u_mass%data_3d                           & !-- u3d         3d u-velocity interpolated to theta points (m/s)
                    ,v3d=domain%v_mass%data_3d                           & !-- v3d         3d v-velocity interpolated to theta points (m/s)
                    ,th3d=domain%potential_temperature%data_3d           &
                    ,t3d=domain%temperature%data_3d                      &
                    ,qv3d=domain%water_vapor%data_3d                     &
                    ,qc3d=domain%cloud_water_mass%data_3d                & !-- qc3d        cloud water mixing ratio (kg/kg)
                    ,qi3d=domain%cloud_ice_mass%data_3d                  & !-- qi3d        cloud ice mixing ratio (kg/kg)
                    ,p3d=domain%pressure%data_3d                         & !-- p3d         3d pressure (pa)
                    ,p3di=domain%pressure_interface%data_3d              & !-- p3di        3d pressure (pa) at interface level
                    ,pi3d=domain%exner%data_3d                           & !-- pi3d        3d exner function (dimensionless)
                    ,rublten=domain%tend%u                               & ! i/o
                    ,rvblten=domain%tend%v                  & ! i/o
                    ,rthblten=domain%tend%th_pbl            & ! i/o
                    ,rqvblten=domain%tend%qv_pbl            & ! i/o
                    ,rqcblten=domain%tend%qc_pbl            & ! i/o
                    ,rqiblten=domain%tend%qi_pbl            & ! i/o
                    ,flag_qi=.True.                         & ! not used in ysu code, so can be whatever?
                    ,cp=cp                                  &
                    ,g=gravity                              &
                    ,rovcp=rcp                            & ! rovcp = Rd/cp
                    ,rd=R_d                                 &  ! J/(kg K) specific gas constant for dry air
                    ,rovg=rovg                              &
                    ,dz8w=domain%dz_interface%data_3d       & !-- dz8w        dz between full levels (m)
                    ,xlv=XLV                    & !-- xlv         latent heat of vaporization (j/kg)
                    ,rv=r_v                                  &  ! J/(kg K) specific gas constant for wet/moist air
                    ,psfc=domain%surface_pressure%data_2d   &
                    ,znt=domain%roughness_z0%data_2d       &  ! i/o -- znt		roughness length (m) (input only)
                    ,ust=domain%ustar                       & ! i/o -- ust		u* in similarity theory (m/s)
                    ,hpbl=domain%hpbl%data_2d               & ! i/o -- hpbl	pbl height (m) - intent(inout)
                    ,psim=domain%fm%data_2d               & !-- psim        similarity stability function for momentum - intent(in)
                    ,psih=domain%fh%data_2d               & !-- psih        similarity stability function for heat- intent(in)
                    ,xland=real(domain%land_mask)                               &
                    ,hfx=domain%sensible_heat%data_2d                     & !  HFX  - net upward heat flux at the surface (W/m^2)
                    ,qfx=domain%qfx%data_2d           & !  QFX  - net upward moisture flux at the surface (kg/m^2/s)
                    !,UOCE=uoce,VOCE=voce                                  & !ocean currents -- not currently used
                    !,CTOPO=ctopo,CTOPO2=ctopo2                            & !optional, only applied to momentum tendencies, not currently used
                    ,YSU_TOPDOWN_PBLMIX=options%pbl_options%ysu_topdown_pblmix                &
                    ,wspd=windspd                           & ! i/o -- wspd        wind speed at lowest model level (m/s)
                    ,br=domain%br%data_2d                   & !-- br          bulk richardson number in surface layer
                    ,dt=dt_in                               & !-- dt		time step (s)
                    ,kpbl2d=domain%kpbl                          & ! o --     ?? k layer of pbl top??
                    ,ep1=EP_1                                & !-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
                    ,ep2=EP_2                                & !-- ep2         constant for specific humidity calculation
                    ,karman=karman                          & !-- karman      von karman constant
                    ,RTHRATEN=RTHRATEN                                    &
!                    ,WSTAR=wstar,DELTA=delta                              &  !Output variables of YSU which we currently dont use
                    ,exch_h=domain%coeff_heat_exchange_3d%data_3d  & ! i/o -- exch_h ! exchange coefficient for heat, K m/s , but 3d??
                    ,exch_m=domain%coeff_momentum_exchange_3d%data_3d  & ! i/o -- exch_h ! exchange coefficient for heat, K m/s , but 3d??
                    ,u10=domain%u_10m%data_2d               &
                    ,v10=domain%v_10m%data_2d               &
                    ,ids=ids, ide=ide, jds=jds, jde=jde     &
                    ,kds=kds, kde=kde, ims=ims, ime=ime     &
                    ,jms=jms, jme=jme, kms=kms, kme=kme     &
                    ,its=its, ite=ite, jts=jts, jte=jte     &
                    ,kts=kts, kte=kte-1                     &
                !optional
                    ,regime=regime                          )!  i/o -- regime	flag indicating pbl regime (stable, unstable, etc.) - not used?

                    ! if(this_image()==1 .and. options%parameters%debug) write(*,*) "  pbl height/lev is:", maxval(domain%hpbl%data_2d ),"m/", maxval(domain%kpbl)  ! uncomment if you want to see the pbl height.

            !> ------------  add tendency terms  ------------
            !
            ! Here the tendency terms that were calculated by the ysu routine are added to the domain-wide fields.
            ! For u and v, we need to re-balance the uvw fields and re-compute dt after we change them. This is done in the
            ! step routine in time_step.f90, after the pbl call.
            !
            !> -----------------------------------------------

            ! Offset u/v tendencies to u and v grid, then add
            ! call array_offset_x_3d(domain%tend%u , tend_u_ugrid)
            ! call array_offset_y_3d(domain%tend%v , tend_v_vgrid)

            ! domain%u%data_3d   =  domain%u%data_3d  +  tend_u_ugrid  * dt_in
            ! domain%v%data_3d   =  domain%v%data_3d  +  tend_v_vgrid  * dt_in

            ! add mass grid tendencies
            domain%water_vapor%data_3d            =  domain%water_vapor%data_3d            + domain%tend%qv_pbl  * dt_in
            domain%cloud_water_mass%data_3d       =  domain%cloud_water_mass%data_3d       + domain%tend%qc_pbl  * dt_in
            domain%potential_temperature%data_3d  =  domain%potential_temperature%data_3d  + domain%tend%th_pbl  * dt_in
            domain%cloud_ice_mass%data_3d         =  domain%cloud_ice_mass%data_3d         + domain%tend%qi_pbl  * dt_in
            
            ! Reset tendencies before the next pbl call. (not sure if necessary)
            domain%tend%qv_pbl    = 0
            domain%tend%th_pbl    = 0
            domain%tend%qc_pbl    = 0
            domain%tend%qi_pbl    = 0
            domain%tend%u         = 0
            domain%tend%v         = 0



            ! -------------------- omp loop   - how to deal with offset (v) grid??   ---------------
            ! ! $omp parallel private(j) &
            ! ! $omp default(shared)
            ! ! $omp do schedule(static)
            ! do j=jts,jte ! OMP  loop

                ! domain%u%data_3d(:,:,j)  =  domain%u%data_3d(:,:,j) + tend_u_ugrid(:,:,j) * dt_in
                ! ! domain%v%data_3d(:,:,j)            =  domain%v%data_3d(:,:,j)       + domain%tend%v(:,:,j) * dt_in

                ! domain%water_vapor%data_3d(:,:,j)  =  domain%water_vapor%data_3d(:,:,j)  +  domain%tend%qv_pbl(:,:,j) * dt_in
                ! domain%cloud_water_mass%data_3d(:,:,j)      = domain%cloud_water_mass%data_3d(:,:,j)      + domain%tend%qc_pbl(:,:,j) * dt_in
                ! domain%potential_temperature%data_3d(:,:,j) = domain%potential_temperature%data_3d(:,:,j) + domain%tend%th_pbl(:,:,j) * dt_in
                ! domain%cloud_ice_mass%data_3d(:,:,j)        = domain%cloud_ice_mass%data_3d(:,:,j)        + domain%tend%qi_pbl(:,:,j) * dt_in

                ! ! Reset tendencies before the next pbl call. (necessary?)
                ! domain%tend%qv_pbl(:,:,j)   = 0
                ! domain%tend%th_pbl(:,:,j)   = 0
                ! domain%tend%qc_pbl(:,:,j)   = 0
                ! domain%tend%qi_pbl(:,:,j)   = 0

            ! enddo
            ! ! $omp end do
            ! ! $omp end parallel

        endif ! End YSU call

    end subroutine pbl

    subroutine pbl_finalize(options)
        implicit none
        type(options_t), intent(in) :: options

        if (options%physics%boundarylayer==kPBL_SIMPLE) then
            call finalize_simple_pbl()
        else if (options%physics%boundarylayer==kPBL_DIAGNOSTIC) then
            call finalize_diagnostic_pbl()
        endif

    end subroutine pbl_finalize
end module planetary_boundary_layer
