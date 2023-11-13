!>----------------------------------------------------------
!! This module provides a wrapper to call various land surface models
!!
!! It sets up variables specific to the LSM to be used including both
!! history variables not currently stored in the domain level data
!! structure, and runtime parameters
!!
!! The main entry point to the code is lsm(domain,options,dt,model_time)
!!
!! <pre>
!! Call tree graph :
!!  lsm_init->[ allocate_noah_data,
!!              external initialization routines]
!!  lsm->[  sat_mr,
!!          calc_exchange_coefficient,
!!          external LSM routines]
!!
!! High level routine descriptions / purpose
!!   lsm_init           - allocates module data and initializes physics package
!!   lsm                - sets up and calls main physics package
!!  calc_exchange_coefficient - calculates surface exchange coefficient (for Noah)
!!  allocate_noah_data  - allocate module level data for Noah LSM
!!  apply_fluxes        - apply LSM fluxes (e.g. sensible and latent heat fluxes) to atmosphere
!!
!! Inputs: domain, options, dt, model_time
!!      domain,options  = as defined in data_structures
!!      dt              = time step (seconds)
!!      model_time      = time since beginning date (seconds)
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module land_surface
    use module_sf_noahdrv,   only : lsm_noah, lsm_noah_init
    use module_sf_noahmpdrv, only : noahmplsm, noahmp_init
    ! use module_lsm_basic,    only : lsm_basic
    ! use module_lsm_simple,   only : lsm_simple, lsm_simple_init
    use module_water_simple, only : water_simple
    use module_water_lake,   only : lake, lakeini, nlevsoil, nlevsnow, nlevlake
    use mod_atm_utilities,   only : sat_mr, calc_Richardson_nr
    use time_object,         only : Time_type
    use data_structures
    use icar_constants,      only : kVARS, kLSM_SIMPLE, kLSM_NOAH, kLSM_NOAHMP, kPBL_DIAGNOSTIC, kSM_FSM
    use options_interface,   only : options_t
    use domain_interface,    only : domain_t
    use module_ra_simple,    only : calc_solar_elevation_corr
    use ieee_arithmetic
    use mod_wrf_constants,   only : gravity, KARMAN, cp, R_d, XLV, rcp, STBOLT, epsilon
    use module_sf_FSMdrv,   only : sm_FSM_init, sm_FSM

    implicit none

    private
    public :: lsm_init, lsm, lsm_var_request

    ! Noah LSM required variables.  Some of these should be stored in domain, but tested here for now
    integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions

    ! LOTS of variables required by Noah, placed here "temporarily", this may be where some stay.
    ! Keeping them at the module level prevents having to allocate/deallocate every call
    ! also avoids adding LOTS of LSM variables to the main domain datastrucurt
    real,allocatable, dimension(:,:)    :: SMSTAV,SFCRUNOFF,UDRUNOFF, SW,                           &
                                           SNOW,SNOWC,SNOWH, ACSNOW, ACSNOM, SNOALB,           &
                                           QGH, GSW, ALBEDO, ALBBCK, Z0, XICE,               &
                                           EMBCK, QSFC, CPM, SR, CHS, CHS2, CQS2,                   &
                                           CHKLOWQ, LAI, QZ0, VEGFRAC, SHDMIN,SHDMAX,SNOTIME,SNOPCX,&
                                           POTEVP,RIB, NOAHRES,FLX4_2D,FVB_2D,FBUR_2D,              &
                                           FGSN_2D, z_atm,Ri,       &
                                           current_precipitation, nmp_snow, nmp_snowh
    double precision,allocatable, dimension(:,:)    :: RAINBL

    integer,allocatable, dimension(:,:) :: rain_bucket ! used to start the previous time step rain bucket

    logical :: MYJ, FRPCPN,ua_phys,RDLAI2D,USEMONALB
    real,allocatable, dimension(:,:,:)  :: SH2O,SMCREL
    real,allocatable, dimension(:,:)    :: dTemp,lhdQV, windspd
    real,allocatable, dimension(:)      :: Zs,DZs
    real :: XICE_THRESHOLD
    integer,allocatable, dimension(:,:) :: IVGTYP,ISLTYP ! IVGTYP not used?
    integer :: ITIMESTEP, update_interval, cur_vegmonth

!     real, parameter :: kappa=0.4 ! this should be karman from data_structure
    real, parameter :: freezing_threshold=273.15
    real, parameter :: SMALL_PRESSURE=0.1 ! note: 0.1Pa is very small 1e-10 wouldn't affect a single-precision float
    real, parameter :: SMALL_QV=1e-10
    real, parameter :: MAX_EXCHANGE_C = 0.5
    real, parameter :: MIN_EXCHANGE_C = 0.004

    character(len=MAXVARLENGTH) :: MMINLU
    logical :: FNDSOILW,FNDSNOWH,RDMAXALB
    integer :: num_soil_layers,ISURBAN,ISICE,ISWATER, ISLAKE
    integer :: exchange_term
    real*8  :: last_model_time
    real*8  :: Delta_t !! MJ added to detect the time for outputting 
    real    :: lh_feedback_fraction, sh_feedback_fraction
    real    :: sfc_layer_thickness


    !Noah-MP specific
    real    :: NMP_SOILTSTEP
    integer :: IDVEG,IOPT_CRS,IOPT_BTR,IOPT_RUN,IOPT_SFC,IOPT_FRZ,IOPT_INF,IOPT_RAD,IOPT_ALB,IOPT_SNF,IOPT_TBOT, IOPT_TDRN, IOPT_NMPOUTPUT
    integer :: IOPT_STC, IOPT_GLA, IOPT_RSF, IOPT_SOIL, IOPT_PEDO, IOPT_CROP, IOPT_IRR, IOPT_IRRM, IZ0TLND, SF_URBAN_PHYSICS
    real,allocatable,dimension(:,:) :: chstarxy
    character(len=MAXVARLENGTH) :: landuse_name
    real, allocatable :: day_frac(:), solar_elevation(:)

    ! MJ added for FSM
    double precision,allocatable, dimension(:,:)    :: SNOWBL    
    real,allocatable, dimension(:,:) :: current_snow, current_rain
    integer,allocatable, dimension(:,:) :: snow_bucket      

    ! Lake model: (allocated on lake init)
    real, allocatable, dimension (:,:)      ::      TH2 !, savedtke12d  lakedepth2d,
    integer :: lakeflag, lake_depth_flag, use_lakedepth, lake_count
    LOGICAL, allocatable, DIMENSION( :,: ) :: lake_or_not

contains


    subroutine lsm_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%landsurface == kLSM_NOAH .or. options%physics%landsurface == kLSM_BASIC) then
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature,       &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%shortwave,     &
                         kVARS%longwave, kVARS%vegetation_fraction, kVARS%canopy_water, kVARS%snow_water_equivalent,    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0, kVARS%ustar,        &
                         kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2,                                                  &
                         kVARS%snow_height, kVARS%lai, kVARS%temperature_2m_veg, kVARS%albedo,                          &
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask, kVARS%land_emissivity])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature,       &
                         kVARS%density, kVARS%pressure_interface, kVARS%shortwave,                                      &
                         kVARS%longwave, kVARS%canopy_water, kVARS%snow_water_equivalent,                               &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%snow_height, kVARS%QFX,  kVARS%land_emissivity,                                          &  ! BK 2020/10/26
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0])!, kVARS%veg_type])    ! BK uncommented 2021/03/20
                         ! kVARS%soil_type, kVARS%land_mask, kVARS%vegetation_fraction]
        endif

        if (options%physics%landsurface == kLSM_NOAHMP) then
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature,       &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%shortwave,     &
                         kVARS%shortwave_direct, kVARS%shortwave_diffuse, kVARS%albedo,                                 &
                         kVARS%longwave, kVARS%vegetation_fraction, kVARS%canopy_water, kVARS%snow_water_equivalent,    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0, kVARS%ustar,        &
                         kVARS%snow_height, kVARS%canopy_vapor_pressure, kVARS%canopy_temperature,                      &
                         kVARS%veg_leaf_temperature, kVARS%coeff_momentum_drag, kVARS%hpbl,                             &
                         kVARS%canopy_fwet, kVARS%snow_water_eq_prev, kVARS%water_table_depth, kVARS%water_aquifer,     &
                         kVARS%mass_leaf, kVARS%mass_root, kVARS%mass_stem, kVARS%mass_wood, kVARS%soil_carbon_fast,    &
                         kVARS%soil_carbon_stable, kVARS%eq_soil_moisture, kVARS%smc_watertable_deep, kVARS%recharge,   &
                         kVARS%recharge_deep, kVARS%storage_lake, kVARS%storage_gw, kVARS%mass_ag_grain,                &
                         kVARS%growing_degree_days, kVARS%plant_growth_stage, kVARS%temperature_2m_veg,                 &
                         kVARS%temperature_2m_bare, kVARS%mixing_ratio_2m_veg, kVARS%mixing_ratio_2m_bare,              &
                         kVARS%surface_rad_temperature, kVARS%net_ecosystem_exchange, kVARS%gross_primary_prod,         &
                         kVARS%net_primary_prod, kVARS%runoff_surface, kVARS%runoff_subsurface,                         &
                         kVARS%evap_canopy, kVARS%evap_soil_surface, kVARS%rad_absorbed_total, kVARS%rad_net_longwave,  &
                         kVARS%apar, kVARS%photosynthesis_total, kVARS%rad_absorbed_veg, kVARS%rad_absorbed_bare,       &
                         kVARS%stomatal_resist_total, kVARS%stomatal_resist_sun, kVARS%stomatal_resist_shade,           &
                         kVARS%lai, kVARS%sai, kVARS%snow_albedo_prev, kVARS%snow_age_factor, kVARS%canopy_water_ice,   &
                         kVARS%canopy_water_liquid, kVARS%vegetation_fraction_max, kVARS%crop_category,                 &
                         kVARS%date_planting, kVARS%date_harvest, kVARS%growing_season_gdd, kVARS%transpiration_rate,   &
                         kVARS%frac_within_gap, kVARS%frac_between_gap, kVARS%ground_temperature_canopy,                &
                         kVARS%ground_temperature_bare, kVARS%ch_veg, kVARS%ch_veg_2m, kVARS%ch_bare, kVARS%ch_bare_2m, &
                         kVARS%ch_under_canopy, kVARS%ch_leaf, kVARS%sensible_heat_veg, kVARS%sensible_heat_bare,       &
                         kVARS%sensible_heat_canopy, kVARS%evap_heat_veg, kVARS%evap_heat_bare, kVARS%evap_heat_canopy, &
                         kVARS%transpiration_heat, kVARS%ground_heat_veg, kVARS%ground_heat_bare, kVARS%snow_nlayers,   &
                         kVARS%net_longwave_veg, kVARS%net_longwave_bare, kVARS%net_longwave_canopy,                    &
                         kVARS%irr_frac_total, kVARS%irr_frac_sprinkler, kVARS%irr_frac_micro, kVARS%irr_frac_flood,    &
                         kVARS%irr_eventno_sprinkler, kVARS%irr_eventno_micro, kVARS%irr_eventno_flood,                 &
                         kVARS%irr_alloc_sprinkler, kVARS%irr_alloc_micro, kVARS%irr_alloc_flood, kVARS%irr_amt_flood,  &
                         kVARS%irr_evap_loss_sprinkler, kVARS%irr_amt_sprinkler, kVARS%irr_amt_micro,                   &
                         kVARS%evap_heat_sprinkler, kVARS%snowfall_ground, kVARS%rainfall_ground, kVARS%crop_type,      &
                         kVARS%ground_surf_temperature, kVARS%snow_temperature, kVARS%snow_layer_depth,                 &
                         kVARS%snow_layer_ice, kVARS%snow_layer_liquid_water, kVARS%soil_texture_1, kVARS%gecros_state, &
                         kVARS%soil_texture_2, kVARS%soil_texture_3, kVARS%soil_texture_4, kVARS%soil_sand_and_clay,    &
                         kVARS%vegetation_fraction_out, kVARS%latitude, kVARS%longitude, kVARS%cosine_zenith_angle,     &
                         kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2,                                                  &
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask, kVARS%land_emissivity])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature,       &
                         kVARS%density, kVARS%pressure_interface, kVARS%shortwave,  kVARS%hpbl, kVARS%land_emissivity,  &
                         kVARS%longwave, kVARS%canopy_water, kVARS%snow_water_equivalent, kVARS%QFX,                    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%snow_height, kVARS%canopy_water_ice, kVARS%canopy_vapor_pressure, kVARS%canopy_temperature,    &  ! BK 2020/10/26
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0])!, kVARS%veg_type])    ! BK uncommented 2021/03/20
                         ! kVARS%soil_type, kVARS%land_mask, kVARS%vegetation_fraction]
        endif



        if (options%physics%snowmodel == kSM_FSM) then
            call options%alloc_vars( &
                         [kVARS%sst, kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature, &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%shortwave,     &
                         kVARS%longwave, kVARS%vegetation_fraction, kVARS%canopy_water, kVARS%snow_water_equivalent,    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0, kVARS%ustar,        &
                         kVARS%snow_height, kVARS%lai, kVARS%temperature_2m_veg,                                        &
                         kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2, kVARS%land_emissivity,                           &
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask, kVARS%snowfall, kVARS%albedo,                &
                         kVARS%runoff_tstep, kVARS%Tsnow, kVARS%Sice, kVARS%Sliq, kVARS%Ds, kVARS%fsnow, kVARS%Nsnow,   &
                         kVARS%rainfall_tstep, kVARS%shd, kVARS%snowfall_tstep, kVARS%meltflux_out_tstep, kVARS%Sliq_out, &
                         kVARS%windspd_10m, kVARS%factor_p, kVARS%dm_salt, kVARS%dm_susp, kVARS%dm_subl, kVARS%dm_slide])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%exch_vars([kVARS%Ds, kVARS%Nsnow, kVARS%fsnow, kVARS%Sice, kVARS%Sliq, kVARS%Tsnow])
             
             call options%restart_vars( &
                         [kVARS%sst, kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature, &
                         kVARS%density, kVARS%pressure_interface, kVARS%shortwave,                                      &
                         kVARS%longwave, kVARS%canopy_water, kVARS%snow_water_equivalent,                               &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%snow_height,  kVARS%snowfall, kVARS%albedo, kVARS%QFX, kVARS%land_emissivity,            &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0,                     &
                         kVARS%runoff_tstep, kVARS%Tsnow, kVARS%Sice, kVARS%Sliq, kVARS%Ds, kVARS%fsnow, kVARS%Nsnow  ])
        endif

       if (options%physics%watersurface > 1) then
            call options%alloc_vars( &
                         [kVARS%sst, kVARS%ustar, kVARS%surface_pressure, kVARS%water_vapor,            &
                         kVARS%temperature, kVARS%sensible_heat, kVARS%latent_heat, kVARS%land_mask,    &
                         kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2,                                            &
                         kVARS%humidity_2m, kVARS%temperature_2m, kVARS%skin_temperature, kVARS%u_10m, kVARS%v_10m])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%sst, kVARS%potential_temperature, kVARS%water_vapor, kVARS%skin_temperature,        &
                         kVARS%surface_pressure, kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,  &
                         kVARS%QFX,                                                                                 &
                         kVARS%humidity_2m, kVARS%temperature_2m])
        endif

        if (options%physics%watersurface == kWATER_LAKE ) then
            call options%alloc_vars( &
            [kVARS%lake_depth,kVARS%veg_type,kVARS%soil_type, kVARS%land_mask,kVARS%terrain,                        &
            kVARS%temperature,kVARS%pressure_interface, kVARS%dz_interface, kVARS%shortwave,  kVARS%longwave,       &
            kVARS%water_vapor, kVARS%latitude, kVARS%longitude, kVARS%sensible_heat, kVARS%latent_heat,             &
            kVARS%ground_heat_flux, kVARS%snow_water_equivalent, kVARS%t_lake3d, kVARS%dz_lake3d,                   &
            kVARS%t_soisno3d, kVARS%h2osoi_ice3d, kVARS%h2osoi_liq3d, kVARS%h2osoi_vol3d, kVARS%z3d,                &
            kVARS%dz3d, kVARS%watsat3d, kVARS%csol3d, kVARS%tkmg3d, kVARS%lakemask, kVARS%zi3d,                     &
            kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2, kVARS%land_emissivity,                                    &
            kVARS%tksatu3d, kVARS%tkdry3d, kVARS%snl2d, kVARS%t_grnd2d,  kVARS%savedtke12d, kVARS%lakedepth2d,      & !  kVARS%snowdp2d, kVARS%h2osno2d,
            kVARS%lake_icefrac3d, kVARS%z_lake3d,kVARS%water_vapor, kVARS%potential_temperature     ])

            ! call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

            call options%restart_vars( &
            [kVARS%lake_depth,kVARS%veg_type,kVARS%soil_type, kVARS%land_mask,kVARS%terrain,                        &
            kVARS%temperature,kVARS%pressure_interface, kVARS%dz_interface, kVARS%shortwave,  kVARS%longwave,       &
            kVARS%water_vapor, kVARS%latitude, kVARS%longitude, kVARS%sensible_heat, kVARS%latent_heat,             &
            kVARS%ground_heat_flux, kVARS%snow_water_equivalent, kVARS%t_lake3d, kVARS%dz_lake3d,                   &
            kVARS%t_soisno3d, kVARS%h2osoi_ice3d, kVARS%h2osoi_liq3d, kVARS%h2osoi_vol3d, kVARS%z3d,                &
            kVARS%dz3d, kVARS%watsat3d, kVARS%csol3d, kVARS%tkmg3d, kVARS%lakemask, kVARS%zi3d,                     &
            kVARS%QFX, kVARS%land_emissivity,                                                                       &
            kVARS%tksatu3d, kVARS%tkdry3d, kVARS%snl2d, kVARS%t_grnd2d, kVARS%savedtke12d, kVARS%lakedepth2d,       & !kVARS%snowdp2d, kVARS%h2osno2d,
            kVARS%lake_icefrac3d, kVARS%z_lake3d,kVARS%water_vapor, kVARS%potential_temperature ])
        endif



    end subroutine lsm_var_request

    subroutine calc_exchange_coefficient(wind,tskin,z0,airt,exchange_C)
        implicit none
        real, dimension(ims:ime,jms:jme),intent(in) :: wind, tskin, z0
        real, dimension(ims:ime,kms:kme,jms:jme),intent(in) :: airt
        real,dimension(ims:ime,jms:jme),intent(inout) :: exchange_C
        
        real, dimension(ims:ime,jms:jme)  :: lnz_atm_term, base_exchange_term
        
        lnz_atm_term = log((z_atm+z0)/z0)
        base_exchange_term=(75*karman**2 * sqrt((z_atm+z0)/z0)) / (lnz_atm_term**2)
        lnz_atm_term=(karman/lnz_atm_term)**2

        ! Richardson number
        exchange_C = 0

        Ri = gravity/airt(:,1,:) * (airt(:,1,:)-tskin)*z_atm/(wind**2+epsilon)
        ! Ri now is a function in atm_utlilities:
        ! calc_Richardson_nr(Ri, airt, tskin, z_atm, wind)

        ! "--------------------------------------------------"
        !  "Surface Richardson number"
        where(Ri<0)  exchange_C = lnz_atm_term * (1.0-(15.0*Ri)/(1.0+(base_exchange_term * sqrt((-1.0)*Ri))))
        where(Ri>=0) exchange_C = lnz_atm_term * 1.0/((1.0+15.0*Ri)*sqrt(1.0+5.0*Ri))

        where(exchange_C > MAX_EXCHANGE_C) exchange_C=MAX_EXCHANGE_C
        where(exchange_C < MIN_EXCHANGE_C) exchange_C=MIN_EXCHANGE_C
    end subroutine calc_exchange_coefficient


! eqn A11 in Appendix A.2 of Chen et al 1997 (see below for reference)
    subroutine F2_formula(F2, z0, Ri)
        real, dimension(ims:ime,jms:jme), intent(inout) :: F2
        real, dimension(ims:ime,jms:jme), intent(in)    :: z0, Ri
        
        real, dimension(ims:ime,jms:jme)  :: lnz_atm_term
        
        lnz_atm_term = log((z_atm+z0)/z0)

        ! for the stable case from Mahrt (1987)
        where(Ri>=0) F2=exp(-Ri)
        ! for the unstable case from Holtslag and Beljaars (1989)
        where(Ri<0)  F2=(1-(15*Ri)/(1+((70.5*karman**2 * sqrt(-Ri * z_atm/z0))/(lnz_atm_term**2))) )

    end subroutine F2_formula
    
    subroutine psi_m_stable_formula(rat,psi_m)
        real, dimension(:,:),intent(in) :: rat
        real, dimension(:,:),intent(inout) :: psi_m

        real :: a, b

        a = 6.1
        b = 2.5
        
        psi_m = -a * log(rat + (1 + rat**(b))**(1/b))
    
    end subroutine psi_m_stable_formula
    
    subroutine psi_h_stable_formula(rat,psi_h)
        real, dimension(:,:),intent(in) :: rat
        real, dimension(:,:),intent(inout) :: psi_h

        real :: c, d

        c = 5.3
        d = 1.1
        
        psi_h = -c * log(rat + (1 + rat**(d))**(1/d))
    
    end subroutine psi_h_stable_formula
    
    
!From Appendix A.2 in Chen et al 1997
! Impact of Atmospheric Surface-layer Parameterizations in the new Land-surface Scheme of the Ncep Mesoscale ETA Model
! Boundary-Layer Meteorology 85:391-421
    subroutine calc_mahrt_holtslag_exchange_coefficient(wind,tskin,airt,znt,exchange_C)
        implicit none
        real, dimension(:,:),intent(in) :: wind,tskin,znt
        real, dimension(:,:,:),intent(in) :: airt
        real,dimension(:,:),intent(inout) :: exchange_C
        
        real, dimension(ims:ime,jms:jme)  :: lnz_atm_term, base_exchange_term
        
        lnz_atm_term = log((z_atm+znt)/znt)
        
        ! Richardson number
        Ri = gravity/airt(:,1,:) * (airt(:,1,:)-tskin)*z_atm/(wind**2+epsilon)

        call F2_formula(base_exchange_term,znt,Ri)

        exchange_C = karman**2 * base_exchange_term / lnz_atm_term**2

        where(exchange_C > MAX_EXCHANGE_C) exchange_C=MAX_EXCHANGE_C
        where(exchange_C < MIN_EXCHANGE_C) exchange_C=MIN_EXCHANGE_C
    end subroutine calc_mahrt_holtslag_exchange_coefficient


    subroutine surface_diagnostics(HFX, QFX, TSK, QSFC, CHS2, CQS2,T2, Q2, PSFC, &
                                    VEGFRAC, veg_type, land_mask, T2veg, T2bare, Q2veg, Q2bare)
        ! taken almost directly / exactly from WRF's module_sf_sfcdiags.F
        !-- HFX           net upward heat flux at the surface (W/m^2)
        !-- QFX           net upward moisture flux at the surface (kg/m^2/s)
        !-- TSK           surface temperature (K)
        !-- qsfc          specific humidity at lower boundary (kg/kg)
        implicit none
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  HFX, QFX, TSK, QSFC
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(INOUT) ::  Q2, T2
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  PSFC, CHS2, CQS2
        REAL, DIMENSION( : , : ),  POINTER, INTENT(IN)    ::  T2veg, T2bare, Q2veg, Q2bare
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  VEGFRAC
        INTEGER, DIMENSION(ims:ime, jms:jme ), INTENT(IN) ::  land_mask, veg_type
        integer :: i,j
        real :: rho

        ! !$omp parallel default(shared), private(i,j,rho)
        ! !$omp do
        do j=jts,jte
            do i=its,ite

                ! if ((domain%veg_type(i,j)/=13).and.(domain%veg_type(i,j)/=15).and.(domain%veg_type(i,j)/=16).and.(domain%veg_type(i,j)/=21)) then
                ! over glacier, urban and barren, noahmp veg 2m T is 0 or -9999e35
                if ((T2veg(i,j) > 200).and.(land_mask(i,j)==kLC_LAND).and.(associated(T2bare))) then
                    T2(i,j) = VEGFRAC(i,j) * T2veg(i,j) &
                        + (1-VEGFRAC(i,j)) * T2bare(i,j)
                    Q2(i,j) = VEGFRAC(i,j) * Q2veg(i,j) &
                        + (1-VEGFRAC(i,j)) * Q2bare(i,j)
                else
                    ! over glacier we don't want to use the bare ground temperature though
                    if ((veg_type(i,j)/=ISICE)           &  ! was /=15  (15=snow/ice in MODIFIED_IGBP_MODIS_NOAH)
                        .and.(veg_type(i,j)/=ISLAKE)     &  ! was /=21  (ISLAKE  = options%lsm_options%lake_category)   # 17 is water, 21 is lakes (in MODIFIED_IGBP_MODIS_NOAH ) MODIFY FOR GENERIC LU TYPES!!
                        .and.(land_mask(i,j)==kLC_LAND)  &
                        .and.(associated(T2bare))) then
                        T2(i,j) = T2bare(i,j)
                        Q2(i,j) = Q2bare(i,j)
                    else
                        RHO = PSFC(I,J)/(R_d * TSK(I,J))

                        if(CQS2(I,J).lt.1.E-3) then
                           Q2(I,J) = QSFC(I,J)
                        else
                           Q2(I,J) = QSFC(I,J) - QFX(I,J)/(RHO*CQS2(I,J))
                        endif

                        if(CHS2(I,J).lt.1.E-3) then
                           T2(I,J) = TSK(I,J)
                        else
                           T2(I,J) = TSK(I,J) - HFX(I,J)/(RHO*CP*CHS2(I,J))
                        endif
                    endif
                endif
                if (Q2(i,j) < SMALL_QV) Q2(i,j) = SMALL_QV

                ! TH2(I,J) = T2(I,J)*(1.E5/PSFC(I,J))**ROVCP
            enddo
        enddo
        ! !$omp end do
        ! !$omp end parallel
    end subroutine surface_diagnostics
    

    subroutine surface_diagnostics_FSM(HFX, QFX, TSK, QSFC, CHS2, CQS2,T2, Q2, PSFC)
        ! taken almost directly / exactly from WRF's module_sf_sfcdiags.F
        !-- HFX           net upward heat flux at the surface (W/m^2)
        !-- QFX           net upward moisture flux at the surface (kg/m^2/s)
        !-- TSK           surface temperature (K)
        !-- qsfc          specific humidity at lower boundary (kg/kg)

        implicit none
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  HFX, QFX, TSK, QSFC
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(INOUT) ::  Q2, T2
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  PSFC, CHS2, CQS2
        integer :: i,j
        real :: rho
        real :: CP_=1005.

        !$omp parallel default(shared), private(i,j,rho)
        !$omp do
        do j=jts,jte
            do i=its,ite
                RHO = PSFC(I,J)/(R_d * TSK(I,J))
                !if(CQS2(I,J).lt.1.E-3) then
                !   Q2(I,J) = QSFC(I,J)
                !else
                   Q2(I,J) = QSFC(I,J) - QFX(I,J)/(RHO*CQS2(I,J))
                !endif
                !if(CHS2(I,J).lt.1.E-3) then
                !   T2(I,J) = TSK(I,J)
                !else
                   T2(I,J) = TSK(I,J) - HFX(I,J)/(RHO*CP_*CHS2(I,J))
                !endif
                if (Q2(i,j) < SMALL_QV) Q2(i,j) = SMALL_QV
                ! TH2(I,J) = T2(I,J)*(1.E5/PSFC(I,J))**ROVCP
                 if ( isnan(HFX(I,J)) .or. abs(TSK(I,J)-T2(I,J))>30 ) write(*,*) "img-H222",i,j,this_image(), HFX(I,J), TSK(I,J), T2(I,J), CHS2(I,J) 
                 !if (this_image()==1) write(*,*),"img-H222",i,j,this_image(), RHO, TSK(I,J), T2(I,J),PSFC(I,J)
            enddo
        enddo
        !$omp end do
        !$omp end parallel
         !write(*,*),"RHO*CP",this_image(), RHO,CP,rd
    end subroutine surface_diagnostics_FSM
        

    subroutine apply_fluxes(domain,dt)
        ! add sensible and latent heat fluxes to the first atm level
        implicit none
        type(domain_t), intent(inout) :: domain
        real, intent(in) :: dt
        integer :: i,j,k
        integer, SAVE :: nz = 0
        real :: layer_fraction

        if (nz==0) then
            layer_fraction = 0
            do k=kts, kte
                layer_fraction = maxval(domain%dz_interface%data_3d(:,k,:)) + layer_fraction
                if (layer_fraction < sfc_layer_thickness) nz=k
            end do
        end if

        associate(density       => domain%density%data_3d,             &
                  sensible_heat => domain%sensible_heat%data_2d,           &
                  latent_heat   => domain%latent_heat%data_2d,             &
                  dz            => domain%dz_interface%data_3d,        &
                  pii           => domain%exner%data_3d,               &
                  th            => domain%potential_temperature%data_3d, &
                  qv            => domain%water_vapor%data_3d          &
            )

        do j = jts, jte
        do k = kts, kts + nz
        do i = its, ite
            ! compute the fraction of the current gridcell that is within the surface layer
            if (k==kts) Then
                layer_fraction = min(1.0, sfc_layer_thickness / dz(i,k,j))
            else
                layer_fraction = max(0.0, min(1.0, (sfc_layer_thickness - sum(dz(i,kts:k-1,j))) / dz(i,k,j) ) )
            endif

            ! convert sensible heat flux to a temperature delta term
            ! (J/(s*m^2) * s / (J/(kg*K)) => kg*K/m^2) ... /((kg/m^3) * m) => K
            dTemp(i,j) = (sh_feedback_fraction * sensible_heat(i,j) * dt/cp)  &
                     / (density(i,k,j) * sfc_layer_thickness)
            ! add temperature delta converted back to potential temperature
            th(i,k,j) = th(i,k,j) + (dTemp(i,j) / pii(i,k,j)) * layer_fraction

            ! convert latent heat flux to a mixing ratio tendancy term
            ! (J/(s*m^2) * s / (J/kg) => kg/m^2) ... / (kg/m^3 * m) => kg/kg
            lhdQV(i,j) = (lh_feedback_fraction * latent_heat(i,j) / XLV * dt) &
                    / (density(i,k,j) * sfc_layer_thickness)
            ! add water vapor in kg/kg
            qv(i,k,j) = qv(i,k,j) + lhdQV(i,j) * layer_fraction

        end do ! i
        end do ! k
        end do ! j

        ! write(*,*) MINVAL(lhdQV), MAXVAL(lhdQV), 'kg/kg (min/max) added to QV at', domain%model_time%hour

        ! enforce some minimum water vapor content... just in case
        where(qv < SMALL_QV) qv = SMALL_QV
        end associate

    end subroutine apply_fluxes

    subroutine allocate_noah_data(num_soil_layers)
        implicit none
        integer, intent(in) :: num_soil_layers
        integer :: i

        ITIMESTEP=1
        
        allocate(SMSTAV(ims:ime,jms:jme))
        SMSTAV = 0.5 !average soil moisture available for transp (between SMCWLT and SMCMAX)
        allocate(SFCRUNOFF(ims:ime,jms:jme))
        SFCRUNOFF = 0
        allocate(UDRUNOFF(ims:ime,jms:jme))
        UDRUNOFF = 0
        allocate(SNOW(ims:ime,jms:jme))
        SNOW = 0
        allocate(SNOWC(ims:ime,jms:jme))
        SNOWC = 0
        allocate(SNOWH(ims:ime,jms:jme))
        SNOWH = 0
        allocate(ACSNOW(ims:ime,jms:jme))
        ACSNOW = 0
        allocate(ACSNOM(ims:ime,jms:jme))
        ACSNOM = 0
        allocate(SNOALB(ims:ime,jms:jme))
        SNOALB = 0.8

        allocate(QGH(ims:ime,jms:jme))
        QGH = 0.02 ! saturated mixing ratio at ~20C
        allocate(GSW(ims:ime,jms:jme))
        GSW = 0

        allocate(ALBEDO(ims:ime,jms:jme))
        ALBEDO = 0.17
        allocate(ALBBCK(ims:ime,jms:jme))
        ALBBCK = 0.17 !?
        allocate(XICE(ims:ime,jms:jme))
        XICE = 0
        allocate(EMBCK(ims:ime,jms:jme))
        EMBCK = 0.99
        allocate(CPM(ims:ime,jms:jme))
        CPM = 0
        allocate(SR(ims:ime,jms:jme))
        SR = 0
        allocate(CHKLOWQ(ims:ime,jms:jme))
        CHKLOWQ = 0
        allocate(QZ0(ims:ime,jms:jme))
        QZ0 = 0 ! used to check for saturation? but only relevant if myj == True

        allocate(FLX4_2D(ims:ime,jms:jme))
        allocate(FVB_2D(ims:ime,jms:jme))
        allocate(FBUR_2D(ims:ime,jms:jme))
        allocate(FGSN_2D(ims:ime,jms:jme))

        allocate(SHDMIN(ims:ime,jms:jme))
        SHDMIN = 0
        allocate(SHDMAX(ims:ime,jms:jme))
        SHDMAX = 100
        allocate(SNOTIME(ims:ime,jms:jme))
        SNOTIME = 0
        allocate(SNOPCX(ims:ime,jms:jme))
        SNOPCX = 0
        allocate(POTEVP(ims:ime,jms:jme))
        POTEVP = 0
        allocate(SMCREL(ims:ime,num_soil_layers,jms:jme))
        SMCREL = 0
        allocate(RIB(ims:ime,jms:jme))
        RIB = 0
        allocate(NOAHRES(ims:ime,jms:jme))
        NOAHRES = 0
        allocate(VEGFRAC(ims:ime,jms:jme))
        VEGFRAC = 50

        allocate(nmp_snow(ims:ime,jms:jme))
        allocate(nmp_snowh(ims:ime,jms:jme))

        allocate(day_frac(ims:ime))
        allocate(solar_elevation(ims:ime))

        XICE_THRESHOLD = 1
        RDLAI2D = .false. !TLE check back on this one
        USEMONALB = .false.
        MYJ = .false.
        FRPCPN = .false. ! set this to true and calculate snow ratio to use microphysics based snow/rain partitioning
        ua_phys = .false.

        allocate(SH2O(ims:ime,num_soil_layers,jms:jme))
        SH2O = 0.25

        allocate(Zs(num_soil_layers))
        allocate(DZs(num_soil_layers))
        !DZs = [0.1,0.3,0.6,1.0]
        !if (options%physics%landsurface==kSM_FSM) then !! MJ added to adapt with FSM for the soil layer thickness
            DZs = [0.1,0.2,0.4,0.8]
        !endif
        Zs(1) = DZs(1)/2
        do i = 2,num_soil_layers
            Zs(i) = Zs(i-1) + DZs(i)/2 + DZs(i-1)/2
        end do

    end subroutine allocate_noah_data

    subroutine lsm_init(domain,options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        integer :: i,j

        if (options%physics%landsurface == 0) return   !! So we cannot (currently) run without lsm but with water.

        if (this_image()==1) write(*,*) "Initializing LSM"

        if (this_image()==1) write(*,*) "    max soil_deep_temperature on init: ", maxval(domain%soil_deep_temperature%data_2d)
        if (this_image()==1) write(*,*) "    max skin_temperature on init: ", maxval(domain%skin_temperature%data_2d)

        exchange_term = 1

        ! module level variables for easy access... need to think about tiling to permit halo processing separately.
        ids = domain%grid%ids
        ide = domain%grid%ide
        jds = domain%grid%jds
        jde = domain%grid%jde
        kds = domain%grid%kds
        kde = domain%grid%kde
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

        lh_feedback_fraction = options%lsm_options%lh_feedback_fraction
        sh_feedback_fraction = options%lsm_options%sh_feedback_fraction
        sfc_layer_thickness = options%lsm_options%sfc_layer_thickness

        allocate(SW(ims:ime,jms:jme))

        allocate(dTemp(its:ite,jts:jte))
        dTemp = 0
        allocate(lhdQV(its:ite,jts:jte))
        lhdQV = 0
        allocate(Z0(ims:ime,jms:jme))
        Z0 = domain%roughness_z0%data_2d ! this should get updated by the LSM(?)
        allocate(QSFC(ims:ime,jms:jme))
        QSFC = domain%water_vapor%data_3d(:,kms,:) ! this should get updated by the lsm    (BK: BUT does not fed back to domain%water_vapor%data_3d ?)
        allocate(Ri(ims:ime,jms:jme))
        Ri = 0
        allocate(z_atm(ims:ime,jms:jme))
        z_atm = 50

        allocate(current_precipitation(ims:ime,jms:jme))
        current_precipitation = 0

        allocate(windspd(ims:ime,jms:jme))
        windspd = 3

        ! NOTE, these fields have probably not been initialized yet...
        ! windspd = sqrt(domain%u10**2+domain%v10**2)

        if (options%physics%landsurface > kLSM_BASIC) then
            if (options%physics%microphysics == 0) then
                write(*,*) "Land Surface models need microphysics to run"
                stop "Land Surface models need microphysics to run"
            endif
            allocate(RAINBL(ims:ime,jms:jme))
            RAINBL = domain%accumulated_precipitation%data_2dd  ! used to store last time step accumulated precip so that it can be subtracted from the current step
                                ! set to domain%rain incase this is a restart run and rain is non-zero to start
            allocate(rain_bucket(ims:ime,jms:jme))
            rain_bucket = domain%precipitation_bucket  ! used to store last time step accumulated precip so that it can be subtracted from the current step

            ! MJ added:
            allocate(current_snow(ims:ime,jms:jme)) ! MJ added 
            current_snow = 0

            allocate(current_rain(ims:ime,jms:jme)) ! MJ added 
            current_rain = 0

            allocate(SNOWBL(ims:ime,jms:jme))! for snowfall:
            SNOWBL = domain%accumulated_snowfall%data_2dd  ! used to store last time step accumulated precip so that it can be subtracted from the current step

            allocate(snow_bucket(ims:ime,jms:jme))
            snow_bucket = domain%snowfall_bucket
        endif
        
        if (options%physics%landsurface==kLSM_NOAH .or. options%physics%landsurface==kLSM_NOAHMP .or. options%physics%snowmodel==kSM_FSM) then
            num_soil_layers=4 ! to .nml?
            call allocate_noah_data(num_soil_layers)
        endif
        ! initial guesses (not needed?)
        domain%temperature_2m%data_2d = domain%temperature%data_3d(:,kms,:)
        domain%humidity_2m%data_2d = domain%water_vapor%data_3d(:,kms,:)

        if (options%physics%landsurface==kLSM_SIMPLE) then
            write(*,*) "    Simple LSM (may not work?)"
            stop "Simple LSM not settup, choose a different LSM options"
            ! call lsm_simple_init(domain,options)
        endif
        ! Noah Land Surface Model
        if (options%physics%landsurface==kLSM_NOAH) then
            if (this_image()==1) write(*,*) "    Noah LSM"

            ! if (this_image()==1) then
            !     write(*,*) "    options%parameters%external_files: ", trim(options%parameters%external_files)
            !     write(*,*) "    options%parameters%restart: ", options%parameters%restart
            !     write(*,*) "    options%parameters%rho_snow_ext ", trim(options%parameters%rho_snow_ext)
            !     write(*,*) "    options%parameters%swe_ext ", trim(options%parameters%swe_ext )
            ! endif

            if (options%parameters%rho_snow_ext /="" .AND. options%parameters%swe_ext /="") then ! calculate snowheight from external swe and density, but only if both are provided. (Swe alone will give FNDSNW = F)
                FNDSNOWH = .True.
                if (this_image()==1) write(*,*) "    Find snow height in file i.s.o. calculating them from SWE: FNDSNOWH=", FNDSNOWH
            elseif(options%parameters%hsnow_ext /="" ) then  ! read in external snowheight if supplied
                FNDSNOWH= .True.
                if (this_image()==1) write(*,*) "    Find snow height in file i.s.o. calculating them from SWE: FNDSNOWH=", FNDSNOWH
            elseif(options%parameters%restart) then   ! If restarting read in snow height, but only if this is in restart file?
                FNDSNOWH= .True.
                if (this_image()==1) write(*,*) "    Find snow height in file i.s.o. calculating them from SWE: FNDSNOWH=", FNDSNOWH
            else
                FNDSNOWH=.False. ! calculate SNOWH from SNOW
            endif

            !These are still needed for NoahLSM, since it assumes that the exchange coefficients are multiplied by windspeed
            allocate(CHS(ims:ime,jms:jme))
            allocate(CHS2(ims:ime,jms:jme))
            allocate(CQS2(ims:ime,jms:jme))

            FNDSOILW=.False. ! calculate SOILW (this parameter is ignored in LSM_NOAH_INIT)
            RDMAXALB=.False.

            ISURBAN = options%lsm_options%urban_category
            ISICE   = options%lsm_options%ice_category
            ISWATER = options%lsm_options%water_category
            MMINLU  = options%lsm_options%LU_Categories !"MODIFIED_IGBP_MODIS_NOAH"
            ISLAKE  = options%lsm_options%lake_category


            if (options%lsm_options%monthly_albedo) then
                if (.not.options%lsm_options%monthly_vegfrac) Then
                    print*, "ERROR, monthly albedo requires monthly vegfrac"
                    error stop
                endif
                ALBEDO = domain%albedo%data_3d(:, domain%model_time%month, :)
            else
                ALBEDO = domain%albedo%data_3d(:, 1, :)
            endif
            if (options%lsm_options%monthly_vegfrac) then
                VEGFRAC = domain%vegetation_fraction%data_3d(:, domain%model_time%month, :)
            else
                VEGFRAC = domain%vegetation_fraction%data_3d(:, 1, :)
            endif
            cur_vegmonth = domain%model_time%month

            ! save the canopy water in a temporary variable in case this is a restart run because lsm_init resets it to 0
            domain%cqs2%data_2d = domain%canopy_water%data_2d
            ! prevents init from failing when processing water points that may have "soil_t"=0
            where(domain%soil_temperature%data_3d<200) domain%soil_temperature%data_3d=200
            where(domain%soil_water_content%data_3d<0.0001) domain%soil_water_content%data_3d=0.0001

            call LSM_NOAH_INIT( VEGFRAC,                             &
                                domain%snow_water_equivalent%data_2d,& !SNOW, &  BK 18/03/2021
                                SNOWC,                               &
                                domain%snow_height%data_2d,          & !SNOWH, &   BK 18/03/2021
                                domain%canopy_water%data_2d,         &
                                domain%soil_temperature%data_3d,     & !-- SMSTAV      Soil moisture availability for evapotranspiration ( fraction between SMCWLT and SMCMXA)
                                domain%soil_water_content%data_3d,   &
                                SFCRUNOFF,                           &
                                UDRUNOFF,                            &
                                ACSNOW,                              &
                                ACSNOM,                              &
                                domain%veg_type,                     &
                                domain%soil_type,                    &
                                domain%soil_temperature%data_3d,     &
                                domain%soil_water_content%data_3d,   &
                                SH2O,                                &
                                ZS,                                  &
                                DZS,                                 &
                                MMINLU,                              &
                                SNOALB,                              &
                                FNDSOILW,                            &
                                FNDSNOWH,                            &
                                RDMAXALB,                            &
                                num_soil_layers,                     &
                                .False.,                             & ! nlayers, is_restart (can't yet)
                                .True. ,                             & ! allowed_to_read (e.g. soilparm.tbl)
                                ids,ide, jds,jde, kds,kde,           &
                                ims,ime, jms,jme, kms,kme,           &
                                its,ite, jts,jte, kts,kte  )

            domain%canopy_water%data_2d = domain%cqs2%data_2d
            domain%cqs2%data_2d=0.01
            ! where(domain%veg_type==ISWATER) domain%land_mask=kLC_WATER ! ensure VEGTYPE (land cover) and land-sea mask are consistent
            where((domain%veg_type==ISWATER) .OR. (domain%veg_type==ISLAKE)) domain%land_mask=kLC_WATER  ! include lakes.
        endif

        ! Noah-MP Land Surface Model
        if (options%physics%landsurface==kLSM_NOAHMP) then
            if (this_image()==1) write(*,*) "    Noah-MP LSM"

            ! if (this_image()==1) then
            !     write(*,*) "    options%parameters%external_files: ", trim(options%parameters%external_files)
            !     write(*,*) "    options%parameters%restart: ", options%parameters%restart
            !     write(*,*) "    options%parameters%rho_snow_ext ", trim(options%parameters%rho_snow_ext)
            !     write(*,*) "    options%parameters%swe_ext ", trim(options%parameters%swe_ext )
            ! endif

            if (options%parameters%rho_snow_ext /="" .AND. options%parameters%swe_ext /="") then ! calculate snowheight from external swe and density, but only if both are provided. (Swe alone will give FNDSNW = F)
                FNDSNOWH = .True.
                if (this_image()==1) write(*,*) "    Find snow height in file i.s.o. calculating them from SWE: FNDSNOWH=", FNDSNOWH
            elseif(options%parameters%hsnow_ext /="" ) then  ! read in external snowheight if supplied
                FNDSNOWH= .True.
                if (this_image()==1) write(*,*) "    Find snow height in file i.s.o. calculating them from SWE: FNDSNOWH=", FNDSNOWH
            elseif(options%parameters%restart) then   ! If restarting read in snow height, but only if this is in restart file?
                FNDSNOWH= .True.
                if (this_image()==1) write(*,*) "    Find snow height in file i.s.o. calculating them from SWE: FNDSNOWH=", FNDSNOWH
            else
                FNDSNOWH=.False. ! calculate SNOWH from SNOW
            endif

            FNDSOILW=.False. ! calculate SOILW (this parameter is ignored in LSM_NOAH_INIT)
            RDMAXALB=.False.

            ISURBAN = options%lsm_options%urban_category
            ISICE   = options%lsm_options%ice_category
            ISWATER = options%lsm_options%water_category
            MMINLU  = options%lsm_options%LU_Categories !"MODIFIED_IGBP_MODIS_NOAH"
            ISLAKE  = options%lsm_options%lake_category

            if (options%lsm_options%monthly_albedo) then
                if (.not.options%lsm_options%monthly_vegfrac) Then
                    print*, "ERROR, monthly albedo requires monthly vegfrac"
                    error stop
                endif
                ALBEDO = domain%albedo%data_3d(:, domain%model_time%month, :)
            else
                ALBEDO = domain%albedo%data_3d(:, 1, :)
            endif
            if (options%lsm_options%monthly_vegfrac) then
                VEGFRAC = domain%vegetation_fraction%data_3d(:, domain%model_time%month, :)
            else
                VEGFRAC = domain%vegetation_fraction%data_3d(:, 1, :)
            endif
            cur_vegmonth = domain%model_time%month

            ! save the canopy water in a temporary variable in case this is a restart run because lsm_init resets it to 0
            domain%cqs2%data_2d = domain%canopy_water%data_2d
            ! prevents init from failing when processing water points that may have "soil_t"=0
            where(domain%soil_temperature%data_3d<200) domain%soil_temperature%data_3d=200
            where(domain%soil_water_content%data_3d<0.0001) domain%soil_water_content%data_3d=0.0001

            IDVEG = options%lsm_options%nmp_dveg
            IOPT_CRS = options%lsm_options%nmp_opt_crs
            IOPT_SFC = options%lsm_options%nmp_opt_sfc
            IOPT_BTR = options%lsm_options%nmp_opt_btr
            IOPT_RUN = options%lsm_options%nmp_opt_run
            IOPT_INF = options%lsm_options%nmp_opt_inf
            IOPT_FRZ = options%lsm_options%nmp_opt_frz
            IOPT_INF = options%lsm_options%nmp_opt_inf
            IOPT_RAD = options%lsm_options%nmp_opt_rad
            IOPT_ALB = options%lsm_options%nmp_opt_alb
            IOPT_SNF = options%lsm_options%nmp_opt_snf
            IOPT_TBOT = options%lsm_options%nmp_opt_tbot
            IOPT_STC = options%lsm_options%nmp_opt_stc
            IOPT_GLA = options%lsm_options%nmp_opt_gla
            IOPT_RSF = options%lsm_options%nmp_opt_rsf
            IOPT_SOIL = options%lsm_options%nmp_opt_soil
            IOPT_PEDO = options%lsm_options%nmp_opt_pedo
            IOPT_CROP = options%lsm_options%nmp_opt_crop
            IOPT_IRR = options%lsm_options%nmp_opt_irr
            IOPT_IRRM = options%lsm_options%nmp_opt_irrm
            IOPT_TDRN = options%lsm_options%nmp_opt_tdrn
            IOPT_NMPOUTPUT = options%lsm_options%noahmp_output
            IZ0TLND = options%sfc_options%iz0tlnd
            NMP_SOILTSTEP = options%lsm_options%nmp_soiltstep
            SF_URBAN_PHYSICS = options%lsm_options%sf_urban_phys


            if (options%physics%snowmodel==kSM_FSM) IOPT_SNF = 4 !This will turn off precipitation partitioning in NoahMP, letting us remove snow from NoahMP


            !allocate dummy variable that doesn't do anything
            allocate(chstarxy(ims:ime,jms:jme))
            chstarxy = 0

            call NOAHMP_INIT ( MMINLU,                                  &
                                domain%snow_water_equivalent%data_2d,   &
                                domain%snow_height%data_2d,             &
                                domain%canopy_water%data_2d,            &
                                domain%soil_type,                       &
                                domain%veg_type,                        &
                                domain%latitude%data_2d,                &
                                domain%soil_temperature%data_3d,        &
                                domain%soil_water_content%data_3d,      &
                                SH2O , DZS ,                            &
                                FNDSOILW , FNDSNOWH ,                   &
                                domain%skin_temperature%data_2d,        &
                                domain%snow_nlayers,                    &
                                domain%veg_leaf_temperature%data_2d,    &
                                domain%ground_surf_temperature%data_2d, &
                                domain%canopy_water_ice%data_2d,        &
                                domain%soil_deep_temperature%data_2d,   &
                                XICE,                                   &
                                domain%canopy_water_liquid%data_2d,     &
                                domain%canopy_vapor_pressure%data_2d,   &
                                domain%canopy_temperature%data_2d,      &
                                domain%coeff_momentum_drag%data_2d,     &
                                domain%chs%data_2d,                     &
                                domain%canopy_fwet%data_2d,             &
                                domain%snow_water_eq_prev%data_2d,      &
                                domain%snow_albedo_prev%data_2d,        &
                                domain%snowfall_ground%data_2d,         &
                                domain%rainfall_ground%data_2d,         &
                                domain%storage_lake%data_2d,            &
                                domain%water_table_depth%data_2d,       &
                                domain%water_aquifer%data_2d,           &
                                domain%storage_gw%data_2d,              &
                                domain%snow_temperature%data_3d,        &
                                domain%snow_layer_depth%data_3d,        &
                                domain%snow_layer_ice%data_3d,          &
                                domain%snow_layer_liquid_water%data_3d, &
                                domain%mass_leaf%data_2d,               &
                                domain%mass_root%data_2d,               &
                                domain%mass_stem%data_2d,               &
                                domain%mass_wood%data_2d,               &
                                domain%soil_carbon_stable%data_2d,      &
                                domain%soil_carbon_fast%data_2d,        &
                                domain%lai%data_2d,                     &
                                domain%sai%data_2d,                     &
                                domain%mass_ag_grain%data_2d,           &
                                domain%growing_degree_days%data_2d,     &
                                domain%crop_type%data_3d,               &
                                domain%crop_category,                   &
                                domain%irr_eventno_sprinkler,           &
                                domain%irr_eventno_micro,               &
                                domain%irr_eventno_flood,               &
                                domain%irr_alloc_sprinkler%data_2d,     &
                                domain%irr_alloc_micro%data_2d,         &
                                domain%irr_alloc_flood%data_2d,         &
                                domain%irr_evap_loss_sprinkler%data_2d, &
                                domain%irr_amt_sprinkler%data_2d,       &
                                domain%irr_amt_micro%data_2d,           &
                                domain%irr_amt_flood%data_2d,           &
                                domain%evap_heat_sprinkler%data_2d,     &
                                domain%temperature_2m_veg%data_2d,      &
                                domain%temperature_2m_bare%data_2d,     &
                                chstarxy,                               &   !doesn't do anything -_-
                                num_soil_layers,                        &
                                .False.,                                &    !restart
                                .True.,                                 &    !allowed_to_read
                                IOPT_RUN,  IOPT_CROP, IOPT_IRR, IOPT_IRRM, &
                                SF_URBAN_PHYSICS,                         &  ! urban scheme
                                ids,ide, jds,jde, kds,kde,                &
                                ims,ime, jms,jme, kms,kme,                &
                                its,ite, jts,jte, kts,kte)

  !                           TLE: GROUNDWATER OFF FOR NOW
  !                                   smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy, dx, dy, msftx, msfty,&     ! Optional groundwater
  !                                   wtddt    ,stepwtd  ,dt       ,qrfsxy     ,qspringsxy  , qslatxy    ,  &      ! Optional groundwater
  !                                   fdepthxy ,ht     ,riverbedxy ,eqzwt     ,rivercondxy ,pexpxy       ,  &      ! Optional groundwater
  !                                   rechclim,                                                             &      ! Optional groundwater
  !                                   gecros_state)                                                                ! Optional gecros crop

            domain%canopy_water%data_2d = domain%cqs2%data_2d
            domain%cqs2%data_2d=0.01
            ! where(domain%veg_type==ISWATER) domain%land_mask=kLC_WATER ! ensure VEGTYPE (land cover) and land-sea mask are consistent (BK 202208: this does not include lakes!!)
            where((domain%veg_type==ISWATER) .OR. (domain%veg_type==ISLAKE)) domain%land_mask=kLC_WATER

        endif

        if(options%physics%watersurface==kWATER_LAKE) then
        ! ____________ Lake model ______________________
        ! From WRF's /run/README.namelist:  These could at some point become namelist options in ICAR?
        ! lakedepth_default(max_dom)          = 50,      ! default lake depth (If there is no lake_depth information in the input data, then lake depth is assumed to be 50m)
        ! lake_min_elev(max_dom)              = 5,       ! minimum elevation of lakes. May be used to determine whether a water point is a lake in the absence of lake
        !                                                  category. If the landuse type includes 'lake' (i.e. Modis_lake and USGS_LAKE), this variable is of no effects.
        ! use_lakedepth (max_dom)             = 1,       ! option to use lake depth data. Lake depth data is available from 3.6 geogrid program. If one didn't process
        !                                                    the lake depth data, but this switch is set to 1, the program will stop and tell one to go back to geogrid
        !                                                     program.
        !                                     = 0, do not use lake depth data.

            if (this_image()==1) write(*,*) "Initializing Lake model"

            ! allocate arrays:
            allocate( lake_or_not(ims:ime, jms:jme))
            allocate( TH2( ims:ime, jms:jme ))
            if( .not.(allocated(XICE))) then
                allocate(XICE(ims:ime,jms:jme))   ! already allocated for NoahMP, so check?
                XICE = 0
            endif

            ! ISURBAN = options%lsm_options%urban_category
            ISICE   = options%lsm_options%ice_category
            ISWATER = options%lsm_options%water_category
            ! MMINLU  = options%lsm_options%LU_Categories !"MODIFIED_IGBP_MODIS_NOAH"
            ISLAKE  = options%lsm_options%lake_category

            ! allocate_noah_data already sets xice_threshold, so if we are using noah (mp/lsm) leave as is.
            if(.not.(options%physics%landsurface==kLSM_NOAHMP .OR. options%physics%landsurface==kLSM_NOAH)) then
                xice_threshold = 1.0  ! allocate_noah_data sets it to 1., BUT WRF's module_physics_init.F sets xice_threshold to 0.5 .... so?
            endif

            lake_count=0
            if(ISLAKE==-1) then
                if(this_image()==1) write(*,*)  "   WARNING: no lake category in LU data: The model will try to guess lake-gridpoints. This option has not been properly tested!"
                lakeflag=0  ! If no lake cat is provided, the lake model will determine lakes based
                            ! on the criterion (ivgtyp(i,j)==iswater .and. ht(i,j)>=lake_min_elev))
            else
                lakeflag=1
                ! from WRF's module_initialize_real.F:
                DO j = jts, MIN(jde-1,jte)
                    DO i = its, MIN(ide-1,ite)
                    !    IF ( grid%lu_index(i,j) .NE. grid%islake ) THEN
                        if(domain%veg_type(i,j) .NE. ISLAKE ) then
                            domain%lakemask%data_2d(i,j) = 0       ! grid%lakemask(i,j) = 0
                        ELSE
                            domain%lakemask%data_2d(i,j) = 1       ! grid%lakemask(i,j) = 1
                            lake_count= lake_count + 1
                        end if
                    END DO
                END DO
            endif
            ! if(options%parameters%debug) write(*,*)"   ",lake_count, " lake cells in image ", this_image()

            ! setlake_depth_flag and use_lakedepth flag. (They seem to be redundant, but whatever):
            if( associated(domain%lake_depth%data_2d) ) then
                if(this_image()==1) write(*,*) "   Using Lake depth data "
                use_lakedepth = 1
                lake_depth_flag = 1
            else
                use_lakedepth = 0
                lake_depth_flag = 0
            endif

            call lakeini( &
                IVGTYP = domain%veg_type                        &
                ,ISLTYP = domain%soil_type                      &
                ,HT=domain%terrain%data_2d                      & ! terrain height [m] if ht(i,j)>=lake_min_elev -> lake  (grid%ht in WRF)
                ,SNOW=domain%snow_water_equivalent%data_2d      & !i  ! SNOW in kg/m^2  (NoahLSM: SNOW liquid water-equivalent snow depth (m)
                ,lake_min_elev=5.                               & ! minimum elevation of lakes. May be used to determine whether a water point is a lake in the absence of lake category. If the landuse type includes 'lake' (i.e. Modis_lake and USGS_LAKE), this variable is of no effects.
                ,restart=options%parameters%restart             & ! if restart, this (lakeini) subroutine is simply skipped.
                ,lakedepth_default=50.                          & ! default lake depth (If there is no lake_depth information in the input data, then lake depth is assumed to be 50m)
                ,lake_depth=domain%lake_depth%data_2d           & !INTENT(IN)
                ,lakedepth2d=domain%lakedepth2d%data_2d         & !INTENT(OUT) (will be equal to lake_depth if lake_depth data is provided in hi-res input, otherwise lakedepth_default)
                ,savedtke12d=domain%savedtke12d%data_2d         & !INTENT(OUT)
                ,snowdp2d=domain%snow_height%data_2d            & ! domain%snowdp2d%data_2d
                ,h2osno2d=domain%snow_water_equivalent%data_2d  & !domain%h2osno2d%data_2d
                ,snl2d=domain%snl2d%data_2d                     & ! snowlevel 2d?
                ,t_grnd2d=domain%t_grnd2d%data_2d               & ! ground temperature?
                ,t_lake3d=domain%t_lake3d%data_3d               & ! lake temperature 3d
                ,lake_icefrac3d=domain%lake_icefrac3d%data_3d   & ! lake ice fraction ?
                ,z_lake3d=domain%z_lake3d%data_3d               & !
                ,dz_lake3d=domain%dz_lake3d%data_3d             &
                ,t_soisno3d=domain%t_soisno3d%data_3d           & ! temperature of both soil and snow
                ,h2osoi_ice3d=domain%h2osoi_ice3d%data_3d       & !  ice lens (kg/m2)
                ,h2osoi_liq3d=domain%h2osoi_liq3d%data_3d       & ! liquid water (kg/m2)
                ,h2osoi_vol3d=domain%h2osoi_vol3d%data_3d       & ! volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
                ,z3d=domain%z3d%data_3d                         & ! layer depth for snow & soil (m)
                ,dz3d=domain%dz3d%data_3d                       & ! layer thickness for soil or snow (m)
                ,zi3d=domain%zi3d%data_3d                                              &
                ,watsat3d=domain%watsat3d%data_3d               &
                ,csol3d=domain%csol3d%data_3d                   &
                ,tkmg3d=domain%tkmg3d%data_3d                   &
                ,iswater=iswater,       xice=xice,           xice_threshold=xice_threshold                                              &
                ,xland=domain%land_mask                         & !-- XLAND         land mask (1 for land, 2 for water)  i/o
                ,tsk=domain%skin_temperature%data_2d            &
                ,lakemask=domain%lakemask%data_2d               & ! 2d var that says lake(1) or not lake(0)
                ,lakeflag=lakeflag                              & ! flag to read in lakemask (lakeflag=1), or to determine lakemask from ivgtyp(i,j)==iswater.and.ht(i,j)>=lake_min_elev (lakeflag=0)
                ,lake_depth_flag=lake_depth_flag,   use_lakedepth=use_lakedepth               & ! flags to use the provided lake depth data (in hi-res input domain file) or not.
                ,tkdry3d=domain%tkdry3d%data_3d                  &
                ,tksatu3d=domain%tksatu3d%data_3d                  &
                ,lake=lake_or_not                               & ! Logical (:,:) if gridpoint is lake or not (INTENT(OUT)) not used further?
                ,its=its, ite=ite, jts=jts, jte=jte             &
                ,ims=ims, ime=ime, jms=jms, jme=jme             &
                )
        endif

        if (options%physics%snowmodel == kSM_FSM) then
            if (this_image()==1) write(*,*) "    SnowModel: FSM2"
            call sm_FSM_init(domain,options)
            
            !DZs already allocated in alloc_noah above
            !Hard code DZs here just for FSM, since this is also hard-coded in FSM
            DZs = [0.1,0.2,0.4,0.8]
            Zs(1) = DZs(1)/2
            do i = 2,num_soil_layers
                Zs(i) = Zs(i-1) + DZs(i)/2 + DZs(i-1)/2
            end do
            !!
            !allocate(Zs(num_soil_layers))
            !allocate(DZs(num_soil_layers))
            !!
        endif
        
        ! defines the height of the middle of the first model level
        z_atm = domain%z%data_3d(:,kts,:) - domain%terrain%data_2d

        update_interval=options%lsm_options%update_interval
        last_model_time=-999

    end subroutine lsm_init


    subroutine lsm(domain,options,dt)
        implicit none

        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real, intent(in) :: dt
        real :: lsm_dt
        integer :: i,j
        
        if (options%physics%landsurface == 0) return

        if (last_model_time==-999) then
            last_model_time = domain%model_time%seconds()-update_interval
        endif

        if ((domain%model_time%seconds() - last_model_time) >= update_interval) then
            lsm_dt = domain%model_time%seconds() - last_model_time
            last_model_time = domain%model_time%seconds()

            if (options%physics%radiation_downScaling==1  .and. options%physics%radiation>1) then
                SW = domain%shortwave_total%data_2d
            else
                SW = domain%shortwave%data_2d
            endif

            ! if (this_image()==1) write(*,*) "    lsm start: snow_water_equivalent max:", MAXVAL(domain%snow_water_equivalent%data_2d)

            ! exchange coefficients
            windspd = sqrt(domain%u_10m%data_2d**2 + domain%v_10m%data_2d**2)
            
            !If there is no surface layer scheme, then we should compute exchange coefficients here. Assign CHS to CHS2 and CQS2 as well...
            if (options%physics%surfacelayer ==  0) then
                if (exchange_term==1) then
                    call calc_exchange_coefficient(windspd,domain%skin_temperature%data_2d,domain%roughness_z0%data_2d,domain%temperature%data_3d,domain%chs%data_2d)
                elseif (exchange_term==2) then
                    call calc_mahrt_holtslag_exchange_coefficient(windspd,domain%skin_temperature%data_2d,domain%temperature%data_3d,domain%roughness_z0%data_2d,domain%chs%data_2d)
                endif
                
                domain%chs2%data_2d = domain%chs%data_2d
                domain%cqs2%data_2d = domain%chs%data_2d
            endif


            ! --------------------------------------------------
            ! First handle the open water surface options
            ! --------------------------------------------------
            ! if (options%physics%watersurface==kWATER_BASIC) then
                ! Note, do nothing because QFX and QSFC are only used for to calculate diagnostic
                !    T2m and Q2m.  However, the fluxes and stability terms are not coordinated, so
                !    This leads to problems in the current formulation and this has been removed.
                ! do j=1,ny
                !     do i=1,nx
                !         if (domain%landmask(i,j)==kLC_WATER) then
                !             QFX(i,j) = domain%latent_heat(i,j) / LH_vaporization
                !             QSFC(i,j)=sat_mr(domain%T2m(i,j),domain%psfc(i,j))
                !         endif
                !     enddo
                ! enddo
            ! else

            if((options%physics%watersurface==kWATER_SIMPLE) .or.      &
                (options%physics%watersurface==kWATER_LAKE) ) then
                    call water_simple(options,                              &
                                      domain%sst%data_2d(its:ite,jts:jte),                   &
                                      domain%surface_pressure%data_2d(its:ite,jts:jte),      &
                                      windspd(its:ite,jts:jte),                              &
                                      domain%ustar(its:ite,jts:jte),                         &
                                      domain%water_vapor%data_3d(its:ite,kms,jts:jte),       &
                                      domain%temperature%data_3d(its:ite,kms,jts:jte),       &
                                      domain%sensible_heat%data_2d(its:ite,jts:jte),         &
                                      domain%latent_heat%data_2d(its:ite,jts:jte),           &
                                      z_atm(its:ite,jts:jte),                                &
                                      domain%roughness_z0%data_2d(its:ite,jts:jte),          &
                                      domain%land_mask(its:ite,jts:jte),                     &
                                      QSFC(its:ite,jts:jte),                                 &
                                      domain%qfx%data_2d(its:ite,jts:jte),                   &
                                      domain%skin_temperature%data_2d(its:ite,jts:jte),      &
                                      domain%chs%data_2d(its:ite,jts:jte),   &
                                      domain%veg_type(its:ite,jts:jte),                      &
                                      its, ite, kts, kte, jts, jte)
                                !   ,domain%terrain%data_2d               & ! terrain height [m] if ht(i,j)>=lake_min_elev -> lake (in case no lake category is provided, but lake model is selected, we need to not run the simple water as well - left comment in for future reference)
            endif

            !___________________ Lake model _____________________
            ! This lake model (ported from WRF V4.4) is run for the grid cells that are defined as lake in the hi-res input file.
            ! It also is advised to supply a lake_depth parameter in the hi-res input, otherwise the default depth of 50m is used (see lakeini above)
            ! It requires the VEGPARM.TBL landuse category to be one which has a separate lake category (i.e. MODIFIED_IGBP_MODIS_NOAH, USGS-RUC or MODI-RUC).
            ! For the grid cells that are defined as water, but not as lake (i.e. oceans), the simple water model above will be run.
            if (options%physics%watersurface==kWATER_LAKE) then    ! WRF's lake model

                ! current_precipitation = (domain%accumulated_precipitation%data_2dd-RAINBL)+(domain%precipitation_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE  ! analogous to noah calls

                call lake( &
                    t_phy=domain%temperature%data_3d                            & !-- t_phy         temperature (K)     !Temprature at the mid points (K)
                    ,p8w=domain%pressure_interface%data_3d(:,kms:kme,:)         & !-- p8w           pressure at full levels (Pa) ! Naming convention: 8~at => p8w reads as "p-at-w" (w=full levels)
                    ,dz8w=domain%dz_interface%data_3d                           & !-- dz8w          dz between full levels (m)
                    ,qvcurr=domain%water_vapor%data_3d                          &  !i
                    ,u_phy=domain%u_mass%data_3d                                & !-- u_phy         u-velocity interpolated to theta points (m/s)
                    ,v_phy=domain%v_mass%data_3d                                & !-- v_phy         v-velocity interpolated to theta points (m/s)
                    ,glw=domain%longwave%data_2d                                & !-- GLW           downward long wave flux at ground surface (W/m^2)
                    ,emiss=domain%land_emissivity%data_2d                       & !-- EMISS         surface emissivity (between 0 and 1)
                    ,rainbl=current_precipitation                               & ! RAINBL in mm (Accumulation between PBL calls)
                    ,dtbl=lsm_dt                                                & !-- dtbl          timestep (s) or ITIMESTEP?
                    ,swdown=SW                                                  & !-- SWDOWN        downward short wave flux at ground surface (W/m^2)
                    ,albedo=ALBEDO                                              & ! albedo? fixed at 0.17?
                    ,xlat_urb2d=domain%latitude%data_2d                         & ! optional ?
                    ,z_lake3d=domain%z_lake3d%data_3d                           &
                    ,dz_lake3d=domain%dz_lake3d%data_3d                         &
                    ,lakedepth2d=domain%lakedepth2d%data_2d                     &
                    ,watsat3d=domain%watsat3d%data_3d                           &
                    ,csol3d=domain%csol3d%data_3d                               &
                    ,tkmg3d=domain%tkmg3d%data_3d                               &
                    ,tkdry3d=domain%tkdry3d%data_3d        &
                    ,tksatu3d=domain%tksatu3d%data_3d                  &
                    ,ivgtyp=domain%veg_type                                     &
                    ,HT=domain%terrain%data_2d                                  &
                    ,xland=real(domain%land_mask)                               & !-- XLAND         land mask (1 for land, 2 OR 0 for water)  i/o
                    ,iswater=iswater,  xice=xice,   xice_threshold=xice_threshold   &
                    ,lake_min_elev=5.                                           & ! if this value is changed, also change it in lake_ini
                    ,ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde       &
                    ,ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme       &
                    ,its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte       &
                    ,h2osno2d=domain%snow_water_equivalent%data_2d             & !domain%h2osno2d%data_2d
                    ,snowdp2d=domain%snow_height%data_2d                        & ! domain%snowdp2d%data_2d
                    ,snl2d=domain%snl2d%data_2d                                 &
                    ,z3d=domain%z3d%data_3d                                     &
                    ,dz3d=domain%dz3d%data_3d                                   &
                    ,zi3d=domain%zi3d%data_3d                                   &
                    ,h2osoi_vol3d=domain%h2osoi_vol3d%data_3d           &
                    ,h2osoi_liq3d=domain%h2osoi_liq3d%data_3d       &
                    ,h2osoi_ice3d=domain%h2osoi_ice3d%data_3d           &
                    ,t_grnd2d=domain%t_grnd2d%data_2d               &
                    ,t_soisno3d=domain%t_soisno3d%data_3d                                      &
                    ,t_lake3d=domain%t_lake3d%data_3d                           & ! 3d lake temperature (K)
                    ,savedtke12d=domain%savedtke12d%data_2d                &
                    ,lake_icefrac3d=domain%lake_icefrac3d%data_3d    &
                    ,lakemask=domain%lakemask%data_2d                                        &
                    ,lakeflag=lakeflag                                          &
                    ,hfx= domain%sensible_heat%data_2d                          & !(OUT)-- HFX         upward heat flux at the surface (W/m^2)   (INTENT:OUT)
                    ,lh=domain%latent_heat%data_2d                              & !(OUT)-- LH          net upward latent heat flux at surface (W/m^2)
                    ,grdflx=domain%ground_heat_flux%data_2d                     & !(OUT)-- GRDFLX(I,J) ground heat flux (W m-2)
                    ,tsk=domain%skin_temperature%data_2d                        & !(OUT)-- TSK          skin temperature [K]
                    ,qfx=domain%qfx%data_2d                                     & !(OUT)-- QFX        upward moisture flux at the surface (kg/m^2/s) in
                    ,t2= domain%temperature_2m%data_2d                          & !(OUT)-- t2         diagnostic 2-m temperature from surface layer and lsm
                    ,th2=TH2                                                    & !(OUT)-- th2        diagnostic 2-m theta from surface layer and lsm
                    ,q2=domain%humidity_2m%data_2d                              & !(OUT)-- q2         diagnostic 2-m mixing ratio from surface layer and lsm
                )

            endif


            where(windspd<1) windspd=1 ! minimum wind speed to prevent the exchange coefficient from blowing up
            
            
            ! --------------------------------------------------
            ! Now handle the land surface options
            ! --------------------------------------------------
            ! if (options%physics%landsurface==kLSM_BASIC) then
                ! call lsm_basic(domain,options,lsm_dt)
                ! Note, do nothing because QFX and QSFC are only used for to calculate diagnostic
                !    T2m and Q2m.  However, the fluxes and stability terms are not coordinated, so
                !    This leads to problems in the current formulation and this has been removed.
                ! do j=1,ny
                !     do i=1,nx
                !         if (domain%landmask(i,j)==kLC_LAND) then
                !             QFX(i,j) = domain%latent_heat(i,j) / LH_vaporization
                !             QSFC(i,j)=max(domain%water_vapor%data_3d(i,1,j),0.5*sat_mr(domain%T2m(i,j),domain%psfc(i,j)))
                !         endif
                !     enddo
                ! enddo


            ! else
            if (options%physics%landsurface == kLSM_SIMPLE) then
                write(*,*) "--------------------------"
                stop "Simple LSM not implemented yet"
                ! call lsm_simple(domain%th,domain%pii,domain%qv,domain%current_rain, domain%current_snow,domain%p_inter, &
                !                 domain%swdown,domain%lwdown, sqrt(domain%u10**2+domain%v10**2), &
                !                 domain%sensible_heat%data_2d, domain%latent_heat, domain%ground_heat_flux, &
                !                 domain%skin_t, domain%soil_t, domain%soil_vwc, domain%snow_swe, &
                !                 options,lsm_dt)

            else if (options%physics%landsurface == kLSM_NOAH) then
                ! Call the Noah Land Surface Model

                ! 2m saturated mixing ratio
                do j=jms,jme
                    do i=ims,ime
                        if (domain%land_mask(i,j) == kLC_LAND) then
                            QGH(i,j) = sat_mr(domain%temperature_2m%data_2d(i,j),domain%surface_pressure%data_2d(i,j))
                        endif
                    enddo
                enddo
                if (options%lsm_options%monthly_albedo) then
                    if (cur_vegmonth /= domain%model_time%month) then
                        ALBEDO = domain%albedo%data_3d(:, domain%model_time%month, :)
                    endif
                endif
                if (options%lsm_options%monthly_vegfrac) then
                    if (cur_vegmonth /= domain%model_time%month) then
                        VEGFRAC = domain%vegetation_fraction%data_3d(:, domain%model_time%month, :)
                        cur_vegmonth = domain%model_time%month
                    endif
                endif

                ! if (this_image()==1) write(*,*) "    lsm start: accumulated_precipitation max:", MAXVAL(domain%accumulated_precipitation%data_2dd)
                ! if (this_image()==1) write(*,*) "    lsm start: RAINBL max:", MAXVAL(RAINBL)
                ! if (this_image()==1) write(*,*) "    lsm start: domain%precipitation_bucket max:", MAXVAL(domain%precipitation_bucket)
                ! if (this_image()==1) write(*,*) "    lsm start: rain_bucket max:", MAXVAL(rain_bucket)


                ! RAINBL(i,j) = [kg m-2]   RAINBL = domain%accumulated_precipitation%data_2dd  ! used to store last time step accumulated precip so that it can be subtracted from the current step
                current_precipitation = (domain%accumulated_precipitation%data_2dd - RAINBL) !+(domain%precipitation_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE
                if (allocated(domain%rain_fraction)) current_precipitation = current_precipitation * domain%rain_fraction(:,:,domain%model_time%get_month())

                CHS = domain%chs%data_2d*windspd
                CHS2 = domain%chs2%data_2d*windspd
                CQS2 = domain%cqs2%data_2d*windspd


                call lsm_noah(domain%dz_interface%data_3d,                &
                            domain%water_vapor%data_3d,                   &
                            domain%pressure_interface%data_3d(:,kms:kme,:),&
                            domain%temperature%data_3d,                   &
                            domain%skin_temperature%data_2d,              &
                            domain%sensible_heat%data_2d,                 &
                            domain%qfx%data_2d,                           &
                            domain%latent_heat%data_2d,                   &
                            domain%ground_heat_flux%data_2d,              &
                            QGH,                                          &
                            GSW,                                          &
                            SW,                                           &
                            domain%longwave%data_2d,                      &
                            SMSTAV,                                       &
                            domain%soil_totalmoisture%data_2d,            &  ! this is not defined on lsm_init (BK 2021/03/20)
                            SFCRUNOFF,                                    &
                            UDRUNOFF,                                     &
                            domain%veg_type,                              &
                            domain%soil_type,                             &
                            ISURBAN,                                      &
                            ISICE,                                        &
                            VEGFRAC,                                      &
                            ALBEDO,                                       &
                            ALBBCK,                                       &
                            domain%roughness_z0%data_2d,                  &
                            Z0,                                           &
                            domain%soil_deep_temperature%data_2d,         &
                            real(domain%land_mask),                       &
                            XICE,                                         &
                            domain%land_emissivity%data_2d,               &
                            EMBCK,                                        &
                            SNOWC,                                        &
                            QSFC,                                         &
                            current_precipitation,                        &  ! RAINBL
                            MMINLU,                                       &
                            num_soil_layers,                              &
                            lsm_dt,                                       &
                            DZS,                                          &
                            ITIMESTEP,                                    &
                            domain%soil_water_content%data_3d,            &
                            domain%soil_temperature%data_3d,              &
                            domain%snow_water_equivalent%data_2d,         &
                            domain%canopy_water%data_2d,                  &
                            CHS,                                          &
                            CHS2,                                         &
                            CQS2,                                         &
                            CPM,                                          &
                            rcp,                                          &
                            SR,                                           &
                            chklowq,                                      &
                            domain%lai%data_2d,                           &
                            qz0,                                          & !H
                            myj,frpcpn,                                   &
                            SH2O,                                         &
                            domain%snow_height%data_2d,                   &     !SNOWH,                                   & !H
                            SNOALB,SHDMIN,SHDMAX,                         & !I
                            SNOTIME,                                      & !?
                            ACSNOM,ACSNOW,                                & !O
                            SNOPCX,                                       & !O
                            POTEVP,                                       & !O
                            SMCREL,                                       & !O
                            XICE_THRESHOLD,                               &
                            RDLAI2D,USEMONALB,                            &
                            Ri,                                           & !I
                            NOAHRES,                                      &
                            ua_phys,flx4_2d,fvb_2d,fbur_2d,fgsn_2d,       & ! Noah UA changes
                            ids,ide, jds,jde, kds,kde,                    &
                            ims,ime, jms,jme, kms,kme,                    &
                            its,ite, jts,jte, kts,kte)

                where(domain%snow_water_equivalent%data_2d > options%lsm_options%max_swe) domain%snow_water_equivalent%data_2d = options%lsm_options%max_swe
                ! now that znt (roughness_z0) has been updated, we need to recalculate terms
            else if (options%physics%landsurface == kLSM_NOAHMP) then
            ! Call the Noah-MP Land Surface Model

                ! 2m saturated mixing ratio
                do j=jms,jme
                    do i=ims,ime
                        if (domain%land_mask(i,j) == kLC_LAND) then
                            QGH(i,j) = sat_mr(domain%temperature_2m%data_2d(i,j),domain%surface_pressure%data_2d(i,j))
                        endif
                    enddo
                enddo
                if (options%lsm_options%monthly_albedo) then
                    ALBEDO = domain%albedo%data_3d(:, domain%model_time%month, :)
                else
                    ALBEDO = domain%albedo%data_3d(:, 1, :)
                endif
                if (options%lsm_options%monthly_vegfrac) then
                    if (cur_vegmonth /= domain%model_time%month) then
                        VEGFRAC = domain%vegetation_fraction%data_3d(:, domain%model_time%month, :)
                        cur_vegmonth = domain%model_time%month
                    endif
                endif

                !more parameters
                landuse_name = options%lsm_options%LU_Categories            !test whether this works or if we need something separate

                ! if (this_image()==1) write(*,*) "    lsm start: accumulated_precipitation max:", MAXVAL(domain%accumulated_precipitation%data_2d)
                ! if (this_image()==1) write(*,*) "    lsm start: RAINBL max:", MAXVAL(RAINBL)
                ! if (this_image()==1) write(*,*) "    lsm start: domain%precipitation_bucket max:", MAXVAL(domain%precipitation_bucket)
                ! if (this_image()==1) write(*,*) "    lsm start: rain_bucket max:", MAXVAL(rain_bucket)
                
                current_snow = (domain%accumulated_snowfall%data_2dd-SNOWBL)!+(domain%snowfall_bucket-snow_bucket)*kPRECIP_BUCKET_SIZE !! MJ: snowfall in kg m-2
                current_precipitation = (domain%accumulated_precipitation%data_2dd - RAINBL) !+(domain%precipitation_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE
                if (allocated(domain%rain_fraction)) current_precipitation = current_precipitation * domain%rain_fraction(:,:,domain%model_time%get_month())

!                do I = ims,ime
!                  do J = jms,jme
!                    call calc_declin(domain%model_time%day_of_year(),real(domain%model_time%hour),real(domain%model_time%minute),real(domain%model_time%second),domain%latitude%data_2d(I,J),domain%longitude%data_2d(I,J),domain%cos_zenith%data_2d(I,J))
!                  enddo
!                enddo


                do j = jms,jme
                    !! MJ commented as it does not work in Erupe
                    !solar_elevation  = calc_solar_elevation(date=domain%model_time, lon=domain%longitude%data_2d, &
                    !                j=j, ims=ims,ime=ime,jms=jms,jme=jme,its=its,ite=ite,day_frac=day_frac)

                    solar_elevation  = calc_solar_elevation_corr(date=domain%model_time, lon=domain%longitude%data_2d, &
                                     j=j, ims=ims,ime=ime,jms=jms,jme=jme,its=its,ite=ite,day_frac=day_frac)
                    domain%cosine_zenith_angle%data_2d(ims:ime,j)=sin(solar_elevation(ims:ime))
                enddo

                if (options%physics%snowmodel==kSM_FSM) then
                    SR = 0.0 ! This, in combination with setting OPT_SNF to 4 in the LSM_init, will turn off snowfall partitioning in NoahMP
                    current_precipitation = current_precipitation-current_snow ! Now remove snowfall from precipitation, so we only have liquid precip going into NMP
                    
                    nmp_snowh = 0.0
                    nmp_snow = 0.0
                else
                    SR = current_snow/(current_precipitation+epsilon)
                    nmp_snowh = domain%snow_height%data_2d
                    nmp_snow = domain%snow_water_equivalent%data_2d
                endif

                call noahmplsm(ITIMESTEP,                              &
                             domain%model_time%year,                   &
                             domain%model_time%day_of_year(),          &
                             domain%cosine_zenith_angle%data_2d,       &
                             domain%latitude%data_2d,                  &
                             domain%longitude%data_2d,                 &
                             domain%dz_interface%data_3d * options%lsm_options%dz_lsm_modification, & ! domain%dz_interface%data_3d,              & !
                             lsm_dt,                                   &
                             DZS,                                      &
                             num_soil_layers,                          &
                             domain%dx,                                &
                             domain%veg_type,                          &
                             domain%soil_type,                         &
                             VEGFRAC,                                  &
                             domain%vegetation_fraction_max%data_2d,   &
                             domain%soil_deep_temperature%data_2d,     &
                             real(domain%land_mask),                   &
                             XICE,                                     &
                             XICE_THRESHOLD,                           &
                             domain%crop_category,                     &  !only used if iopt_crop>0; not currently set up
                             domain%date_planting%data_2d,             &  !only used if iopt_crop>0; not currently set up
                             domain%date_harvest%data_2d,              &  !only used if iopt_crop>0; not currently set up
                             domain%growing_season_gdd%data_2d,        &  !only used if iopt_crop>0; not currently set up
                             IDVEG, IOPT_CRS,  IOPT_BTR, IOPT_RUN,     &
                             IOPT_SFC, IOPT_FRZ, IOPT_INF, IOPT_RAD,   &
                             IOPT_ALB, IOPT_SNF, IOPT_TBOT, IOPT_STC,  &
                             IOPT_GLA, IOPT_RSF, IOPT_SOIL,IOPT_PEDO,  &
                             IOPT_CROP, IOPT_IRR, IOPT_IRRM, IZ0TLND,  &
                             SF_URBAN_PHYSICS,                         &
                             domain%soil_sand_and_clay%data_3d,        &  ! only used if iopt_soil = 3
                             domain%soil_texture_1%data_2d,            &  ! only used if iopt_soil = 2
                             domain%soil_texture_2%data_2d,            &  ! only used if iopt_soil = 2
                             domain%soil_texture_3%data_2d,            &  ! only used if iopt_soil = 2
                             domain%soil_texture_4%data_2d,            &  ! only used if iopt_soil = 2
                             domain%temperature%data_3d,               &
                             domain%water_vapor%data_3d,               &
                             domain%u_mass%data_3d * options%lsm_options%wind_enhancement, &
                             domain%v_mass%data_3d * options%lsm_options%wind_enhancement, &
                             SW,                                       &
                             domain%shortwave_direct%data_2d,          &  ! only used in urban modules, which are currently disabled
                             domain%shortwave_diffuse%data_2d,         &  ! only used in urban modules, which are currently disabled
                             domain%longwave%data_2d,                  &
                             domain%pressure_interface%data_3d(:,kms:kme,:),&
                             current_precipitation,                    &
                             SR,                                       &
                             domain%irr_frac_total%data_2d,            &  ! only used if iopt_irr > 0
                             domain%irr_frac_sprinkler%data_2d,        &  ! only used if iopt_irr > 0
                             domain%irr_frac_micro%data_2d,            &  ! only used if iopt_irr > 0
                             domain%irr_frac_flood%data_2d,            &  ! only used if iopt_irr > 0
                             domain%skin_temperature%data_2d,          &  ! TSK
                             domain%sensible_heat%data_2d,             &  !  HFX
                             domain%qfx%data_2d,                       &
                             domain%latent_heat%data_2d,               &  ! LH
                             domain%ground_heat_flux%data_2d,          &  ! GRDFLX
                             SMSTAV,                                   &
                             domain%soil_totalmoisture%data_2d,        &
                             SFCRUNOFF, UDRUNOFF,                      &
                             ALBEDO, SNOWC,                            &
                             domain%soil_water_content%data_3d,        &
                             SH2O,                                     &
                             domain%soil_temperature%data_3d,          &
                             nmp_snow,                                 &
                             nmp_snowh,                                &
                             domain%canopy_water%data_2d,              &
                             ACSNOM, ACSNOW,                           &
                             domain%land_emissivity%data_2d, QSFC, Z0, &
                             domain%roughness_z0%data_2d,              &
                             domain%hpbl%data_2d,                      &
                             domain%irr_eventno_sprinkler,             &  ! only used if iopt_irr > 0
                             domain%irr_eventno_micro,                 &  ! only used if iopt_irr > 0
                             domain%irr_eventno_flood,                 &  ! only used if iopt_irr > 0
                             domain%irr_alloc_sprinkler%data_2d,       &  ! only used if iopt_irr > 0
                             domain%irr_alloc_micro%data_2d,           &  ! only used if iopt_irr > 0
                             domain%irr_alloc_flood%data_2d,           &  ! only used if iopt_irr > 0
                             domain%irr_evap_loss_sprinkler%data_2d,   &  ! only used if iopt_irr > 0
                             domain%irr_amt_sprinkler%data_2d,         &  ! only used if iopt_irr > 0
                             domain%irr_amt_micro%data_2d,             &  ! only used if iopt_irr > 0
                             domain%irr_amt_flood%data_2d,             &  ! only used if iopt_irr > 0
                             domain%evap_heat_sprinkler%data_2d,       &  ! only used if iopt_irr > 0
                             landuse_name,                             &
                             domain%snow_nlayers,                      &
                             domain%veg_leaf_temperature%data_2d,      &
                             domain%ground_surf_temperature%data_2d,   &
                             domain%canopy_water_ice%data_2d,          &
                             domain%canopy_water_liquid%data_2d,       &
                             domain%canopy_vapor_pressure%data_2d,     &
                             domain%canopy_temperature%data_2d,        &
                             domain%coeff_momentum_drag%data_2d,       &
                             domain%chs%data_2d,                       &
                             domain%canopy_fwet%data_2d,               &
                             domain%snow_water_eq_prev%data_2d,        &
                             domain%snow_albedo_prev%data_2d,          &
                             domain%snowfall_ground%data_2d,           &
                             domain%rainfall_ground%data_2d,           &
                             domain%storage_lake%data_2d,              &
                             domain%water_table_depth%data_2d,         &
                             domain%water_aquifer%data_2d,             &
                             domain%storage_gw%data_2d,                &
                             domain%snow_temperature%data_3d,          &
                             domain%snow_layer_depth%data_3d,          &
                             domain%snow_layer_ice%data_3d,            &
                             domain%snow_layer_liquid_water%data_3d,   &
                             domain%mass_leaf%data_2d,                 &
                             domain%mass_root%data_2d,                 &
                             domain%mass_stem%data_2d,                 &
                             domain%mass_wood%data_2d,                 &
                             domain%soil_carbon_stable%data_2d,        &
                             domain%soil_carbon_fast%data_2d,          &
                             domain%lai%data_2d,                       &
                             domain%sai%data_2d,                       &
                             domain%snow_age_factor%data_2d,           &
                             domain%eq_soil_moisture%data_3d,          &
                             domain%smc_watertable_deep%data_2d,       &
                             domain%recharge_deep%data_2d,             &
                             domain%recharge%data_2d,                  &
                             domain%mass_ag_grain%data_2d,             &  ! currently left as zeroes; not used if iopt_crop = 0?
                             domain%growing_degree_days%data_2d,       &  ! currently left as zeroes; not used if iopt_crop = 0?
                             domain%plant_growth_stage,                &  ! currently left as zeroes; not used if iopt_crop = 0?
                             domain%gecros_state%data_3d,              &  ! not set up; only used if iopt_crop = 2
                             domain%temperature_2m_veg%data_2d,        &
                             domain%temperature_2m_bare%data_2d,       &
                             domain%mixing_ratio_2m_veg%data_2d,       &
                             domain%mixing_ratio_2m_bare%data_2d,      &
                             domain%surface_rad_temperature%data_2d,   &
              	             domain%net_ecosystem_exchange%data_2d,    &
                             domain%gross_primary_prod%data_2d,        &
                             domain%net_primary_prod%data_2d,          &
                             domain%vegetation_fraction_out%data_2d,   &
                             domain%runoff_surface%data_2d,            &
                             domain%runoff_subsurface%data_2d,         &
                             domain%evap_canopy%data_2d,               &
                             domain%evap_soil_surface%data_2d,         &
                             domain%transpiration_rate%data_2d,        &
                             domain%rad_absorbed_total%data_2d,        &
                             domain%rad_net_longwave%data_2d,          &
                             domain%apar%data_2d,                      &
                             domain%photosynthesis_total%data_2d,      &
                             domain%rad_absorbed_veg%data_2d,          &
                             domain%rad_absorbed_bare%data_2d,         &
                             domain%stomatal_resist_sun%data_2d,       &
                             domain%stomatal_resist_shade%data_2d,     &
                             domain%frac_between_gap%data_2d,          &
                             domain%frac_within_gap%data_2d,           &
                             domain%ground_temperature_canopy%data_2d, &
                             domain%ground_temperature_bare%data_2d,   &
                             domain%ch_veg%data_2d,                    &
                             domain%ch_bare%data_2d,                   &
                             domain%sensible_heat_veg%data_2d,         &
                             domain%sensible_heat_canopy%data_2d,      &
                             domain%sensible_heat_bare%data_2d,        &
                             domain%evap_heat_veg%data_2d,             &
                             domain%evap_heat_bare%data_2d,            &
                             domain%ground_heat_veg%data_2d,           &
                             domain%ground_heat_bare%data_2d,          &
                             domain%net_longwave_veg%data_2d,          &
                             domain%net_longwave_canopy%data_2d,       &
                             domain%net_longwave_bare%data_2d,         &
                             domain%transpiration_heat%data_2d,        &
                             domain%evap_heat_canopy%data_2d,          &
                             domain%ch_leaf%data_2d,                   &
                             domain%ch_under_canopy%data_2d,           &
                             domain%ch_veg_2m%data_2d,                 &
                             domain%ch_bare_2m%data_2d,                &
                             domain%stomatal_resist_total%data_2d,     &
                             ids,ide,  jds,jde,  kds,kde,              &
                             ims,ime,  jms,jme,  kms,kme,              &
                             its,ite,  jts,jte,  kts,kte)
                             
                             
                if (options%lsm_options%monthly_albedo) then
                    domain%albedo%data_3d(:, domain%model_time%month, :) = ALBEDO
                else
                    domain%albedo%data_3d(:, 1, :) = ALBEDO
                endif

                VEGFRAC = domain%vegetation_fraction_out%data_2d(:, :)*100.0

                if ( .not.(options%physics%snowmodel==kSM_FSM)) then
                    domain%snow_height%data_2d = nmp_snowh
                    domain%snow_water_equivalent%data_2d = nmp_snow
                endif

    !         TLE: OMITTING OPTIONAL PRECIP INPUTS FOR NOW
    !                         MP_RAINC, MP_RAINNC, MP_SHCV, MP_SNOW, MP_GRAUP, MP_HAIL     )
                where(domain%snow_water_equivalent%data_2d > options%lsm_options%max_swe) domain%snow_water_equivalent%data_2d = options%lsm_options%max_swe
            endif

            !! MJ added: this block is for FSM as sm.
            if (options%physics%snowmodel == kSM_FSM) then
                current_precipitation = (domain%accumulated_precipitation%data_2dd-RAINBL)!+(domain%precipitation_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE !! MJ: this is the total prep=rainfall+snowfall in kg m-2
                current_snow = (domain%accumulated_snowfall%data_2dd-SNOWBL)!+(domain%snowfall_bucket-snow_bucket)*kPRECIP_BUCKET_SIZE !! MJ: snowfall in kg m-2
                current_rain = max(current_precipitation-current_snow,0.) !! MJ: rainfall in kg m-2
                !!
                domain%windspd_10m%data_2d(its:ite,jts:jte)=windspd(its:ite,jts:jte)
                !!
                call sm_FSM(domain,options,lsm_dt,current_rain,current_snow,windspd)
                !!
                
                do i = 1, num_soil_layers              ! soil
                   SH2O(its:ite,i,jts:jte) =  domain%soil_water_content%data_3d(its:ite,i,jts:jte)    / (1000. * DZS(i))
                end do
                
                !if (.not. options%lsm_options%surface_diagnostics) then
                !    domain%temperature_2m%data_2d = domain%temperature%data_3d(:,kms,:)
                !    domain%humidity_2m%data_2d = domain%water_vapor%data_3d(:,kms,:)
                !endif
            endif
            !!
            if (options%physics%landsurface > kLSM_BASIC) then
            
                RAINBL = domain%accumulated_precipitation%data_2dd
                rain_bucket = domain%precipitation_bucket
                SNOWBL = domain%accumulated_snowfall%data_2dd
                snow_bucket = domain%snowfall_bucket
                
                domain%longwave_up%data_2d = STBOLT * domain%land_emissivity%data_2d * domain%skin_temperature%data_2d**4
                ! accumulate soil moisture over the entire column
                domain%soil_totalmoisture%data_2d = domain%soil_water_content%data_3d(:,1,:) * DZS(1) * 1000
                do i = 2,num_soil_layers
                    domain%soil_totalmoisture%data_2d = domain%soil_totalmoisture%data_2d + domain%soil_water_content%data_3d(:,i,:) * DZS(i) * 1000
                enddo

                ! 2m Air T and Q are not well defined if Tskin is not coupled with the surface fluxes
                call surface_diagnostics(domain%sensible_heat%data_2d,          &
                                         domain%qfx%data_2d,                    &
                                         domain%skin_temperature%data_2d,       &
                                         QSFC,                                  &
                                         domain%chs2%data_2d,                   &
                                         domain%cqs2%data_2d,                   &
                                         domain%temperature_2m%data_2d,         &
                                         domain%humidity_2m%data_2d,            &
                                         domain%surface_pressure%data_2d,       &
                                         (VEGFRAC/100.0),                       &
                                         domain%veg_type,                       &
                                         domain%land_mask,                      &
                                         domain%temperature_2m_veg%data_2d,     &
                                         domain%temperature_2m_bare%data_2d,    &
                                         domain%mixing_ratio_2m_veg%data_2d,    &
                                         domain%mixing_ratio_2m_bare%data_2d)

            endif
            !!
        endif
        
        ! PBL scheme should handle the distribution of sensible and latent heat fluxes. If we are
        ! running the LSM without a PBL scheme, as may be done for High-resolution runs, then 
        ! run apply fluxes to still apply heat fluxes calculated by LSM
        if ( (options%physics%landsurface>0 .or. options%physics%watersurface>0 ) .and. (options%physics%boundarylayer==0)) then
            call apply_fluxes(domain, dt)
        endif

    end subroutine lsm
end module land_surface
