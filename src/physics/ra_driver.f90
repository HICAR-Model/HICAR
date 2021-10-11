!>----------------------------------------------------------
!! This module provides a wrapper to call various radiation models
!! It sets up variables specific to the physics package to be used
!!
!! The main entry point to the code is rad(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!!  radiation_init->[ external initialization routines]
!!  rad->[  external radiation routines]
!!
!! High level routine descriptions / purpose
!!   radiation_init     - initializes physics package
!!   rad                - sets up and calls main physics package
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
module radiation
    use module_ra_simple, only: ra_simple, ra_simple_init
    use module_ra_rrtmg_lw, only: rrtmg_lwinit, rrtmg_lwrad
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use data_structures
    use icar_constants, only : kVARS, cp, Rd, gravity

    implicit none

contains
    subroutine radiation_init(domain,options)
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in) :: options

        if (this_image()==1) write(*,*) "Initializing Radiation"

        if (options%physics%radiation==kRA_BASIC) then
            if (this_image()==1) write(*,*) "    Basic Radiation"
        endif
        if (options%physics%radiation==kRA_SIMPLE) then
            if (this_image()==1) write(*,*) "    Simple Radiation"
            call ra_simple_init(domain, options)
        endif

        if (options%physics%radiation==kRA_RRTMG) then
            if (this_image()==1) write(*,*) "    RRTMG"
            if(.not.allocated(domain%tend%th_lwrad)) &
                allocate(domain%tend%th_lwrad(domain%ims:domain%ime,domain%kms:domain%kme,domain%jms:domain%jme))
            call rrtmg_lwinit(                           &
                p_top=minval(domain%pressure_interface%data_3d(:,domain%kme,:)), allowed_to_read=.TRUE. ,                     &
                ids=domain%ids, ide=domain%ide, jds=domain%jds, jde=domain%jde, kds=domain%kds, kde=domain%kde,                &
                ims=domain%ims, ime=domain%ime, jms=domain%jms, jme=domain%jme, kms=domain%kms, kme=domain%kme,                &
                its=domain%its, ite=domain%ite, jts=domain%jts, jte=domain%jte, kts=domain%kts, kte=domain%kte                 )
        endif

    end subroutine radiation_init


    subroutine ra_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        if (options%physics%radiation == kRA_SIMPLE) then
            call ra_simple_var_request(options)
        endif

        if (options%physics%radiation == kRA_RRTMG) then
            call ra_rrtmg_var_request(options)
        endif

    end subroutine ra_var_request


    !> ----------------------------------------------
    !! Communicate to the master process requesting the variables requred to be allocated, used for restart files, and advected
    !!
    !! ----------------------------------------------
    subroutine ra_simple_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options


        ! List the variables that are required to be allocated for the simple radiation code
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%cloud_fraction,   &
                      kVARS%water_vapor, kVARS%cloud_water,             kVARS%rain_in_air,  kVARS%snow_in_air,      &
                      kVARS%shortwave,   kVARS%longwave,                kVARS%cloud_ice,    kVARS%graupel_in_air])

        ! List the variables that are required to be advected for the simple radiation code
        call options%advect_vars( &
                      [kVARS%potential_temperature] )

        ! List the variables that are required when restarting for the simple radiation code
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature, kVARS%shortwave,   kVARS%longwave, kVARS%cloud_fraction] )

    end subroutine ra_simple_var_request


!> ----------------------------------------------
    !! Communicate to the master process requesting the variables requred to be allocated, used for restart files, and advected
    !!
    !! ----------------------------------------------
    subroutine ra_rrtmg_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options


        ! List the variables that are required to be allocated for the simple radiation code
        call options%alloc_vars( &
                     [kVARS%pressure,     kVARS%pressure_interface,    kVARS%potential_temperature,   kVARS%exner,            &
                      kVARS%water_vapor,  kVARS%cloud_water,           kVARS%rain_in_air,             kVARS%snow_in_air,      &
                      kVARS%shortwave,    kVARS%longwave,              kVARS%cloud_ice,               kVARS%graupel_in_air,   &
                      kVARS%re_cloud,     kVARS%re_ice,                kVARS%re_snow,                 kVARS%out_longwave_rad, &
                      kVARS%land_mask,    kVARS%snow_water_equivalent,                                                        &
                      kVARS%dz_interface, kVARS%skin_temperature,      kVARS%temperature,             kVARS%density,          &
                      kVARS%longwave_cloud_forcing, kVARS%land_emissivity, kVARS%temperature_interface                        &
                      ])


        ! List the variables that are required when restarting for the simple radiation code
        call options%restart_vars( &
                     [kVARS%pressure,     kVARS%pressure_interface,    kVARS%potential_temperature,   kVARS%exner,            &
                      kVARS%water_vapor,  kVARS%cloud_water,           kVARS%rain_in_air,             kVARS%snow_in_air,      &
                      kVARS%shortwave,    kVARS%longwave,              kVARS%cloud_ice,               kVARS%graupel_in_air,   &
                      kVARS%re_cloud,     kVARS%re_ice,                kVARS%re_snow,                 kVARS%out_longwave_rad, &
                      kVARS%snow_water_equivalent,                                                                            &
                      kVARS%dz_interface, kVARS%skin_temperature,      kVARS%temperature,             kVARS%density,          &
                      kVARS%longwave_cloud_forcing, kVARS%land_emissivity, kVARS%temperature_interface   ] )

    end subroutine ra_rrtmg_var_request


    subroutine rad(domain, options, dt, halo, subset)
        implicit none

        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt
        integer,        intent(in), optional :: halo, subset

        integer :: ims, ime, jms, jme, kms, kme
        integer :: its, ite, jts, jte, kts, kte
        integer :: ids, ide, jds, jde, kds, kde

        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        kms = domain%grid%kms
        kme = domain%grid%kme
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte
        kts = domain%grid%kts
        kte = domain%grid%kte

        ids = domain%grid%ids
        ide = domain%grid%ide
        jds = domain%grid%jds
        jde = domain%grid%jde
        kds = domain%grid%kds
        kde = domain%grid%kde

        if (options%physics%radiation==kRA_SIMPLE) then
            call ra_simple(theta = domain%potential_temperature%data_3d,         &
                           pii= domain%exner%data_3d,                            &
                           qv = domain%water_vapor%data_3d,                      &
                           qc = domain%cloud_water_mass%data_3d,                 &
                           qs = domain%snow_mass%data_3d                         &
                                + domain%cloud_ice_mass%data_3d                  &
                                + domain%graupel_mass%data_3d,                   &
                           qr = domain%rain_mass%data_3d,                        &
                           p =  domain%pressure%data_3d,                         &
                           swdown =  domain%shortwave%data_2d,                   &
                           lwdown =  domain%longwave%data_2d,                    &
                           cloud_cover =  domain%cloud_fraction%data_2d,         &
                           lat = domain%latitude%data_2d,                        &
                           lon = domain%longitude%data_2d,                       &
                           date = domain%model_time,                             &
                           options = options,                                    &
                           dt = dt,                                              &
                           ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                           its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
        endif

        if (options%physics%radiation==kRA_RRTMG) then
            call RRTMG_LWRAD(rthratenlw=domain%tend%th_lwrad,                     &
!                lwupt, lwuptc, lwuptcln, lwdnt, lwdntc, lwdntcln, &        !if lwupt defined, all MUST be defined
!                lwupb, lwupbc, lwupbcln, lwdnb, lwdnbc, lwdnbcln, &
                            glw = domain%longwave%data_2d,                        &
                            olr = domain%out_longwave_rad%data_2d,                &
                            lwcf = domain%longwave_cloud_forcing%data_2d,         &
                            emiss = domain%land_emissivity%data_2d,               &
                            p8w = domain%pressure_interface%data_3d,              &
                            p3d = domain%pressure%data_3d,                        &
                            pi3d = domain%exner%data_3d,                          &
                            dz8w = domain%dz_interface%data_3d,                   &
                            tsk = domain%skin_temperature%data_2d,                &
                            t3d = domain%temperature%data_3d,                     &
                            t8w = domain%temperature_interface%data_3d,           &     ! temperature interface
                            rho3d = domain%density%data_3d,                       &
                            r = Rd,                                               &
                            g = gravity,                                          &
                            icloud = 0,                                           & ! set to nonzero if effective radius is available from microphysics
                            warm_rain = .False.,                                  & ! when adding WSM3scheme, add option for .True.
!                            cldfra3d  ,                                          & ! if icloud > 0, include this
                            cldovrlp=1,                                           & ! set to 1 for now. Could make this ICAR namelist option
!                            lradius,iradius,                                     & !goes with CAMMGMP (Morrison Gettelman CAM mp)
                            is_cammgmp_used = .False.,                            & !goes with CAMMGMP (Morrison Gettelman CAM mp)
!                            f_ice_phy, f_rain_phy,                               & !goes with MP option 5 (Ferrier)
                            xland=real(domain%land_mask),                         &
                            xice=real(domain%land_mask)*0,                        & ! should add a variable for sea ice fraction
                            snow=domain%snow_water_equivalent%data_2d,            &
                            qv3d=domain%water_vapor%data_3d,                      &
                            qc3d=domain%cloud_water_mass%data_3d,                 &
                            qr3d=domain%rain_mass%data_3d,                        &
                            qi3d=domain%cloud_ice_mass%data_3d,                   &
                            qs3d=domain%snow_mass%data_3d,                        &
                            qg3d=domain%graupel_mass%data_3d,                     &
!                            o3input, o33d,                             &
!                            f_qv, f_qc, f_qr, f_qi, f_qs, f_qg,                  & ! if icloud > 0, these can be set to true
                            re_cloud = domain%re_cloud%data_3d,                   &
                            re_ice   = domain%re_ice%data_3d,                     &
                            re_snow  = domain%re_snow%data_3d,                    &
                            has_reqc=1,                                           & ! use with icloud > 0
                            has_reqi=1,                                           & ! use with icloud > 0
                            has_reqs=1,                                           & ! use with icloud > 0 ! G. Thompson
!                            tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4,   & ! czhao
!                            tauaerlw5,tauaerlw6,tauaerlw7,tauaerlw8,   & ! czhao
!                            tauaerlw9,tauaerlw10,tauaerlw11,tauaerlw12,   & ! czhao
!                            tauaerlw13,tauaerlw14,tauaerlw15,tauaerlw16,   & ! czhao
!                            aer_ra_feedback,                           & !czhao
!                        !jdfcz                 progn,prescribe,                           & !czhao
!                            progn,
                            calc_clean_atm_diag=0,                               & ! used with wrf_chem !czhao
!                            qndrop3d=domain%cloud_number%data_3d,                & ! used with icould > 0
!                            f_qndrop,                         & ! if icloud > 0, use this
                        !ccc added for time varying gases.
                            yr=domain%model_time%year,                             &
                            julian=domain%model_time%day_of_year(),                &
                        !ccc
                            mp_physics=0,                                          &
                            ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,  &
                            ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,  &
!                            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte &
                            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte-1 &
!                            lwupflx, lwupflxc, lwdnflx, lwdnflxc       &
                            )
            domain%temperature%data_3d = domain%temperature%data_3d+domain%tend%th_lwrad*dt
            domain%potential_temperature%data_3d = domain%temperature%data_3d/domain%exner%data_3d
        endif

    end subroutine rad
end module radiation
