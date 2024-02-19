module domain_interface
  use iso_c_binding
  use options_interface,        only : options_t
  use boundary_interface,       only : boundary_t
  use exchangeable_interface,   only : exchangeable_t
  use grid_interface,           only : grid_t
  use variable_interface,       only : variable_t
  use variable_dict_interface,  only : var_dict_t
  use meta_data_interface,      only : meta_data_t
  use time_object,              only : Time_type
  use time_delta_object,        only : time_delta_t
  use data_structures,          only : interpolable_type, tendencies_type
  use timer_interface,          only : timer_t

  implicit none

  private
  public :: domain_t

  type domain_t
    type(meta_data_t)    :: info
    type(grid_t)         :: grid,   grid8w,  u_grid,   v_grid
    type(grid_t)         :: grid2d, u_grid2d, v_grid2d
    type(grid_t)         :: grid_monthly, grid_soil
    type(grid_t)         :: grid_snow, grid_snowsoil
    type(grid_t)         :: grid_soilcomp, grid_gecros, grid_croptype
    type(grid_t)         :: grid_hlm !! MJ added
    type(grid_t)         :: grid_lake , grid_lake_soisno, grid_lake_soi, grid_lake_soisno_1

    type(Time_type) :: model_time

    ! note that not all variables are allocated at runtime, physics packages must request a variable be created
    ! though variables considered "required" are requested by the domain object itself (e.g. terrain)
    ! core model species to be advected

    ! wind field to control advection
    type(exchangeable_t) :: u
    type(exchangeable_t) :: v
    type(exchangeable_t) :: w

    type(exchangeable_t) :: water_vapor
    type(exchangeable_t) :: potential_temperature
    type(exchangeable_t) :: cloud_water_mass
    type(exchangeable_t) :: cloud_number
    type(exchangeable_t) :: cloud_ice_mass
    type(exchangeable_t) :: cloud_ice_number
    type(exchangeable_t) :: rain_mass
    type(exchangeable_t) :: rain_number
    type(exchangeable_t) :: snow_mass
    type(exchangeable_t) :: snow_number
    type(exchangeable_t) :: graupel_mass
    type(exchangeable_t) :: graupel_number
    type(exchangeable_t) :: ice1_a
    type(exchangeable_t) :: ice1_c
    type(exchangeable_t) :: ice2_mass
    type(exchangeable_t) :: ice2_number
    type(exchangeable_t) :: ice2_a
    type(exchangeable_t) :: ice2_c
    type(exchangeable_t) :: ice3_mass
    type(exchangeable_t) :: ice3_number
    type(exchangeable_t) :: ice3_a
    type(exchangeable_t) :: ice3_c

    ! other model variables (not advected)
    type(variable_t) :: exner
    type(variable_t) :: density
    type(variable_t) :: pressure
    type(variable_t) :: pressure_d
    type(variable_t) :: pressure_interface
    type(variable_t) :: temperature
    type(variable_t) :: z
    type(variable_t) :: dz_interface
    type(variable_t) :: z_interface
    type(variable_t) :: dz_mass
    type(variable_t) :: nsquared
    type(variable_t) :: graupel
    type(variable_t) :: accumulated_precipitation
    integer,allocatable :: precipitation_bucket(:,:)
    type(variable_t) :: accumulated_convective_pcp
    integer,allocatable :: cu_precipitation_bucket(:,:)
    type(variable_t) :: accumulated_snowfall
    type(variable_t) :: precip_in_total
    type(variable_t) :: snowfall_ground
    type(variable_t) :: rainfall_ground
    integer,allocatable :: snowfall_bucket(:,:)
    type(variable_t) :: external_precipitation
    type(variable_t) :: cloud_fraction
    type(variable_t) :: longwave
    type(variable_t) :: shortwave
    type(variable_t) :: shortwave_direct
    type(variable_t) :: shortwave_diffuse
    type(variable_t) :: shortwave_direct_above !! MJ aded
    type(variable_t) :: shortwave_total !! MJ added
    type(variable_t) :: terrain
    type(variable_t) :: forcing_terrain  ! BK 05/2020: The forcing terrain interpolated 2d to the hi-res grid. In order to calculate difference in slope
    type(variable_t) :: forcing_terrain2 ! test 9-6-2020
    type(variable_t) :: u_10m
    type(variable_t) :: v_10m
    type(variable_t) :: windspd_10m
    type(variable_t) :: coeff_momentum_drag
    type(variable_t) :: chs
    type(variable_t) :: chs2
    type(variable_t) :: cqs2
    type(variable_t) :: br
    type(variable_t) :: qfx
    type(variable_t) :: psim
    type(variable_t) :: psih
    type(variable_t) :: fm
    type(variable_t) :: fh
    type(variable_t) :: coeff_momentum_exchange_3d ! used in YSU pbl
    type(variable_t) :: coeff_heat_exchange_3d ! used in YSU pbl
    integer,allocatable :: kpbl(:,:)  ! used in YSU pbl / BMJ cu
    type(variable_t) :: hpbl          ! used in YSU pbl /NSAS cu
    type(variable_t) :: surface_rad_temperature
    type(variable_t) :: temperature_2m
    type(variable_t) :: humidity_2m
    type(variable_t) :: temperature_2m_veg
    type(variable_t) :: temperature_2m_bare
    type(variable_t) :: mixing_ratio_2m_veg
    type(variable_t) :: mixing_ratio_2m_bare
    type(variable_t) :: surface_pressure
    type(variable_t) :: rad_absorbed_total
    type(variable_t) :: rad_absorbed_veg
    type(variable_t) :: rad_absorbed_bare
    type(variable_t) :: rad_net_longwave
    type(variable_t) :: longwave_up
    type(variable_t) :: ground_heat_flux
    type(variable_t) :: sensible_heat
    type(variable_t) :: latent_heat
    integer,allocatable :: veg_type(:,:)
    type(variable_t) :: mass_leaf
    type(variable_t) :: mass_root
    type(variable_t) :: mass_stem
    type(variable_t) :: mass_wood
    integer,allocatable :: soil_type(:,:)
    type(variable_t) :: soil_texture_1
    type(variable_t) :: soil_texture_2
    type(variable_t) :: soil_texture_3
    type(variable_t) :: soil_texture_4
    type(variable_t) :: soil_sand_and_clay
    type(variable_t) :: soil_carbon_stable
    type(variable_t) :: soil_carbon_fast
    type(variable_t) :: roughness_z0
    type(variable_t) :: albedo
    type(variable_t) :: vegetation_fraction
    type(variable_t) :: vegetation_fraction_max
    type(variable_t) :: vegetation_fraction_out
    type(variable_t) :: lai
    type(variable_t) :: sai
    integer,allocatable :: crop_category(:,:)
    type(variable_t) :: crop_type
    type(variable_t) :: date_planting
    type(variable_t) :: date_harvest
    type(variable_t) :: growing_season_gdd
    type(variable_t) :: irr_frac_total
    type(variable_t) :: irr_frac_sprinkler
    type(variable_t) :: irr_frac_micro
    type(variable_t) :: irr_frac_flood
    integer,allocatable :: irr_eventno_sprinkler(:,:)
    integer,allocatable :: irr_eventno_micro(:,:)
    integer,allocatable :: irr_eventno_flood(:,:)
    type(variable_t) :: irr_alloc_sprinkler
    type(variable_t) :: irr_alloc_micro
    type(variable_t) :: irr_alloc_flood
    type(variable_t) :: irr_evap_loss_sprinkler
    type(variable_t) :: irr_amt_sprinkler
    type(variable_t) :: irr_amt_micro
    type(variable_t) :: irr_amt_flood
    type(variable_t) :: evap_heat_sprinkler
    type(variable_t) :: mass_ag_grain
    type(variable_t) :: growing_degree_days
    integer,allocatable :: plant_growth_stage(:,:)
    type(variable_t) :: net_ecosystem_exchange
    type(variable_t) :: gross_primary_prod
    type(variable_t) :: net_primary_prod
    type(variable_t) :: apar
    type(variable_t) :: photosynthesis_total
    type(variable_t) :: stomatal_resist_total
    type(variable_t) :: stomatal_resist_sun
    type(variable_t) :: stomatal_resist_shade
    type(variable_t) :: gecros_state
    type(variable_t) :: canopy_water
    type(variable_t) :: canopy_water_ice
    type(variable_t) :: canopy_water_liquid
    type(variable_t) :: canopy_vapor_pressure
    type(variable_t) :: canopy_temperature
    type(variable_t) :: canopy_fwet
    type(variable_t) :: veg_leaf_temperature
    type(variable_t) :: ground_surf_temperature
    type(variable_t) :: frac_between_gap
    type(variable_t) :: frac_within_gap
    type(variable_t) :: ground_temperature_bare
    type(variable_t) :: ground_temperature_canopy
    type(variable_t) :: snow_water_equivalent
    type(variable_t) :: snow_water_eq_prev
    type(variable_t) :: snow_albedo_prev
    type(variable_t) :: snow_temperature
    type(variable_t) :: snow_layer_depth
    type(variable_t) :: snow_layer_ice
    type(variable_t) :: snow_layer_liquid_water
    type(variable_t) :: snow_age_factor
    type(variable_t) :: snow_height
    type(variable_t) :: snow_nlayers
    type(variable_t) :: skin_temperature
    type(variable_t) :: sst
    type(variable_t) :: soil_water_content
    type(variable_t) :: eq_soil_moisture
    type(variable_t) :: smc_watertable_deep
    type(variable_t) :: recharge
    type(variable_t) :: recharge_deep
    type(variable_t) :: soil_temperature
    type(variable_t) :: runoff_subsurface
    type(variable_t) :: runoff_surface
    type(variable_t) :: evap_canopy
    type(variable_t) :: evap_soil_surface
    type(variable_t) :: transpiration_rate
    type(variable_t) :: ch_veg
    type(variable_t) :: ch_veg_2m
    type(variable_t) :: ch_bare
    type(variable_t) :: ch_bare_2m
    type(variable_t) :: ch_under_canopy
    type(variable_t) :: ch_leaf
    type(variable_t) :: sensible_heat_veg
    type(variable_t) :: sensible_heat_bare
    type(variable_t) :: sensible_heat_canopy
    type(variable_t) :: evap_heat_veg
    type(variable_t) :: evap_heat_bare
    type(variable_t) :: evap_heat_canopy
    type(variable_t) :: transpiration_heat
    type(variable_t) :: ground_heat_veg
    type(variable_t) :: ground_heat_bare
    type(variable_t) :: net_longwave_veg
    type(variable_t) :: net_longwave_bare
    type(variable_t) :: net_longwave_canopy
    type(variable_t) :: soil_totalmoisture
    type(variable_t) :: soil_water_content_liq
    type(variable_t) :: soil_deep_temperature
    type(variable_t) :: water_table_depth
    type(variable_t) :: water_aquifer
    type(variable_t) :: storage_gw
    type(variable_t) :: storage_lake
    ! lake model vars:
    type(variable_t) :: lake_depth
    type(variable_t) :: t_lake3d
    type(variable_t) :: snl2d
    type(variable_t) :: t_grnd2d
    type(variable_t) :: lake_icefrac3d
    type(variable_t) :: z_lake3d
    type(variable_t) :: dz_lake3d
    type(variable_t) :: t_soisno3d
    type(variable_t) :: h2osoi_ice3d
    type(variable_t) :: h2osoi_liq3d! liquid water (kg/m2)
    type(variable_t) :: h2osoi_vol3d! volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
    type(variable_t) :: z3d ! layer depth for snow & soil (m)
    type(variable_t) :: dz3d
    type(variable_t) :: watsat3d
    type(variable_t) :: csol3d
    type(variable_t) :: tkmg3d
    type(variable_t) :: lakemask
    type(variable_t) :: xice
    type(variable_t) :: tksatu3d
    type(variable_t) :: tkdry3d
    type(variable_t) :: zi3d
    type(variable_t) :: savedtke12d
    type(variable_t) :: lakedepth2d
    ! diagnostics
    type(variable_t) :: ivt
    type(variable_t) :: iwv
    type(variable_t) :: iwl
    type(variable_t) :: iwi

    ! link effective radius from microphysics to radiation scheme
    type(variable_t) :: re_cloud
    type(variable_t) :: re_ice
    type(variable_t) :: re_snow
    
    ! ice hydrometeor properties diagnosted by ISHMAEL
    type(variable_t) :: ice1_rho
    type(variable_t) :: ice1_phi
    type(variable_t) :: ice1_vmi
    type(variable_t) :: ice2_rho
    type(variable_t) :: ice2_phi
    type(variable_t) :: ice2_vmi
    type(variable_t) :: ice3_rho
    type(variable_t) :: ice3_phi
    type(variable_t) :: ice3_vmi

    type(variable_t) :: out_longwave_rad
    type(variable_t) :: longwave_cloud_forcing
    type(variable_t) :: shortwave_cloud_forcing
    type(variable_t) :: land_emissivity
    type(variable_t) :: temperature_interface
    type(variable_t) :: cosine_zenith_angle
    type(variable_t) :: tend_swrad

    integer,allocatable :: land_mask(:,:)
    type(variable_t) :: latitude
    type(variable_t) :: longitude
    real, allocatable :: latitude_global(:,:)
    real, allocatable :: longitude_global(:,:)
    type(variable_t) :: u_latitude
    type(variable_t) :: u_longitude
    type(variable_t) :: v_latitude
    type(variable_t) :: v_longitude

    type(variable_t) :: w_real
    type(variable_t) :: u_mass
    type(variable_t) :: v_mass

    type(variable_t) :: alpha  !wind-alpha
    type(variable_t) :: froude !Froude number
    type(variable_t) :: Ri     !Bulk richardson number


    type(tendencies_type) :: tend

    type(var_dict_t) :: variables_to_force
    type(var_dict_t) :: vars_to_out
    
    ! Array listing variables to advect with pointers to local data
    type(var_dict_t) :: adv_vars
    type(var_dict_t) :: exch_vars

    type(interpolable_type) :: geo
    type(interpolable_type) :: geo_agl
    type(interpolable_type) :: geo_u
    type(interpolable_type) :: geo_v

    real :: smooth_height, dx
    integer :: nsmooth

    complex(C_DOUBLE_COMPLEX),  allocatable :: terrain_frequency(:,:) ! FFT(terrain)
    double precision,           allocatable :: costheta(:,:)
    double precision,           allocatable :: sintheta(:,:)
    real,                       allocatable :: relax_filter_2d(:,:)
    real,                       allocatable :: relax_filter_3d(:,:,:)
    real,                       allocatable :: advection_dz(:,:,:)
    real,                       allocatable :: rain_fraction(:,:,:) ! monthly varying fraction to multiple precipitation  [-]
    ! store the ratio between the average dz and each grid cells topographically modified dz (for space varying dz only)
    real,                       allocatable :: jacobian(:,:,:)
    real,                       allocatable :: jacobian_u(:,:,:)
    real,                       allocatable :: jacobian_v(:,:,:)
    real,                       allocatable :: jacobian_w(:,:,:)
    real,                       allocatable :: dzdx(:,:,:) ! change in height with change in x/y position (used to calculate w_real vertical motions)
    real,                       allocatable :: dzdy(:,:,:) ! change in height with change in x/y position (used to calculate w_real vertical motions)
    real,                       allocatable :: dzdx_u(:,:,:) ! change in height with change in x/y position on u-grid
    real,                       allocatable :: dzdy_v(:,:,:) ! change in height with change in x/y position on v-grid
    ! BK 2020/05
    real,                       allocatable :: froude_terrain(:,:,:,:) ! Terrain length-scale to use at each point for Froude Number calculation
    real,                       allocatable :: h1(:,:)     ! the large-scale terrain (h1) for the SLEVE coordinate (achieved by smoothin the org terrain)
    real,                       allocatable :: h2(:,:)     ! the small-scale terrain (h2) for the SLEVE coordinate (difference org and h1 terrain)
    real,                       allocatable :: h1_u(:,:)     ! the large-scale terrain (h1) on the u grid
    real,                       allocatable :: h1_v(:,:)     ! the large-scale terrain (h1) on the v grid
    real,                       allocatable :: h2_u(:,:)     ! the small-scale terrain (h2) on the u grid
    real,                       allocatable :: h2_v(:,:)     ! the small-scale terrain (h2) on the v grid
    real,                       allocatable :: dz_scl(:)  ! the scaled dz levels, required for delta terrain calculation    
    real,                       allocatable :: Sx(:,:,:,:)
    real,                       allocatable :: TPI(:,:)
    real,                       allocatable :: neighbor_TPI(:,:)
    real,                       allocatable :: ustar(:,:)
    real,                       allocatable :: znu(:)
    real,                       allocatable :: znw(:)
    
    ! these data are stored on the domain wide grid even if this process is only looking at a subgrid
    ! these variables are necessary with linear winds, especially with spatially variable dz, to compute the LUT
    real,                       allocatable :: global_terrain(:,:)
    real,                       allocatable :: neighbor_terrain(:,:)
    real,                       allocatable :: global_z_interface(:,:,:)
    real,                       allocatable :: global_dz_interface(:,:,:)


    ! these coarrays are used to send all data to/from a master image for IO... ?
    ! For now this will be taken care of in the boundary conditions object
    ! real, allocatable :: transfer_3d(:,:,:)[:]
    ! real, allocatable :: transfer_2d(:,:)[:]
    
    ! Neighboring images of this image
    integer, allocatable :: neighbors(:)
    integer, allocatable :: corner_neighbors(:)

    real, allocatable :: south_in_3d(:,:,:,:)[:]
    real, allocatable :: north_in_3d(:,:,:,:)[:]
    real, allocatable :: west_in_3d(:,:,:,:)[:]
    real, allocatable :: east_in_3d(:,:,:,:)[:]

    real, allocatable :: north_buffer_3d(:,:,:,:)
    real, allocatable :: south_buffer_3d(:,:,:,:)
    real, allocatable :: east_buffer_3d(:,:,:,:)
    real, allocatable :: west_buffer_3d(:,:,:,:)

    real, allocatable :: south_in_2d(:,:,:)[:]
    real, allocatable :: north_in_2d(:,:,:)[:]
    real, allocatable :: west_in_2d(:,:,:)[:]
    real, allocatable :: east_in_2d(:,:,:)[:]

    real, allocatable :: north_buffer_2d(:,:,:)
    real, allocatable :: south_buffer_2d(:,:,:)
    real, allocatable :: east_buffer_2d(:,:,:)
    real, allocatable :: west_buffer_2d(:,:,:)

    ! MPI communicator object for doing parallel communications among domain objects
    integer, public :: IO_comms

    ! contains the size of the domain (or the local tile?)
    integer :: nx, ny, nz, nx_global, ny_global
    integer :: ximg, ximages, yimg, yimages
    integer :: north_neighbor, south_neighbor, east_neighbor, west_neighbor
    integer :: northwest_neighbor, southwest_neighbor, northeast_neighbor, southeast_neighbor


    logical :: north_boundary = .True.
    logical :: south_boundary = .True.
    logical :: east_boundary = .True.
    logical :: west_boundary = .True.

    ! store the start (s) and end (e) for the i,j,k dimensions
    integer ::  ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
                ims,ime, jms,jme, kms,kme, & ! for the memory in these arrays (m)
                its,ite, jts,jte, kts,kte, & ! for the data tile to process   (t)
                ihs,ihe, jhs,jhe, khs,khe    ! for the neighborhood arrays for non-local calculations (h)

    integer :: neighborhood_max ! The maximum neighborhood radius in indices
    
    !! MJ added new vars needed for FSM
    !real,allocatable :: FSM_slopemu(:,:)
    type(variable_t) :: runoff_tstep 
    type(variable_t) :: Tsnow
    type(variable_t) :: Sice
    type(variable_t) :: Sliq
    type(variable_t) :: Ds
    type(variable_t) :: fsnow
    type(variable_t) :: dm_salt
    type(variable_t) :: dm_susp
    type(variable_t) :: dm_subl
    type(variable_t) :: dm_slide
    type(variable_t) :: Nsnow
    type(variable_t) :: rainfall_tstep
    type(variable_t) :: snowfall_tstep
    type(variable_t) :: meltflux_out_tstep
    type(variable_t) :: slope
    type(variable_t) :: slope_angle
    type(variable_t) :: aspect_angle
    type(variable_t) :: svf
    type(variable_t) :: factor_p
    type(variable_t) :: Sliq_out
    type(variable_t) :: hlm
    type(variable_t) :: ridge_dist
    type(variable_t) :: valley_dist
    type(variable_t) :: ridge_drop
    type(variable_t) :: shd


  contains
    procedure :: init
    procedure :: var_request
    
    procedure :: halo_send
    procedure :: halo_retrieve
    procedure :: halo_exchange
    procedure :: halo_3d_send_batch
    procedure :: halo_3d_retrieve_batch
    procedure :: halo_3d_exchange_batch
    procedure :: halo_2d_send_batch
    procedure :: halo_2d_retrieve_batch
    procedure :: halo_2d_exchange_batch
    procedure :: enforce_limits

    procedure :: get_initial_conditions
    procedure :: diagnostic_update
    procedure :: interpolate_forcing
    procedure :: interpolate_external
    procedure :: update_delta_fields
    procedure :: apply_forcing

  end type

  integer, parameter :: space_dimension=3

  interface

    ! Set default component values
    module subroutine init(this, options)
        implicit none
        class(domain_t), intent(inout) :: this
        type(options_t), intent(inout) :: options
    end subroutine
    
    ! read initial atmospheric conditions from forcing data
    module subroutine get_initial_conditions(this, forcing, options, external_conditions)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(inout) :: forcing
        type(boundary_t), intent(inout), optional :: external_conditions  ! external data such as SWE
        type(options_t),  intent(in)    :: options
    end subroutine

    module subroutine diagnostic_update(this,options)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(options_t),  intent(in)    :: options
    end subroutine

    module subroutine interpolate_external(this, external_conditions, options)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(in)    :: external_conditions
        type(options_t),  intent(in)    :: options
    end subroutine

    module subroutine interpolate_forcing(this, forcing, update)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(in)    :: forcing
        logical,          intent(in),   optional :: update
    end subroutine

    module subroutine var_request(this, options)
        implicit none
        class(domain_t), intent(inout) :: this
        type(options_t), intent(inout)  :: options
    end subroutine

    module subroutine halo_send(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine halo_retrieve(this, wait_timer)
        implicit none
        class(domain_t), intent(inout) :: this
        type(timer_t),   intent(inout) :: wait_timer
    end subroutine

    ! Exchange subdomain boundary information
    module subroutine halo_exchange(this, send_timer, ret_timer, wait_timer)
        implicit none
        class(domain_t), intent(inout) :: this
        type(timer_t),   intent(inout) :: send_timer, ret_timer, wait_timer
    end subroutine

    module subroutine halo_3d_send_batch(this, exch_var_only)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, intent(in) :: exch_var_only
    end subroutine

    module subroutine halo_3d_retrieve_batch(this, exch_var_only, wait_timer)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, intent(in) :: exch_var_only
        type(timer_t), optional,   intent(inout) :: wait_timer
    end subroutine

    ! Exchange subdomain boundary information as a batched exchange
    module subroutine halo_3d_exchange_batch(this, send_timer, ret_timer, wait_timer, exch_var_only)
        implicit none
        class(domain_t), intent(inout) :: this
        type(timer_t),   optional, intent(inout) :: send_timer, ret_timer, wait_timer
        logical, optional, intent(in) :: exch_var_only
    end subroutine

    module subroutine halo_2d_send_batch(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine halo_2d_retrieve_batch(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    ! Exchange subdomain boundary information as a batched exchange
    module subroutine halo_2d_exchange_batch(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    ! Make sure no hydrometeors are getting below 0
    module subroutine enforce_limits(this,update_in)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, optional, intent(in)  :: update_in
    end subroutine

    module subroutine update_delta_fields(this, dt)
        implicit none
        class(domain_t),    intent(inout) :: this
        type(time_delta_t), intent(in)    :: dt
    end subroutine

    module subroutine apply_forcing(this, forcing, options, dt)
        implicit none
        class(domain_t),    intent(inout) :: this
        class(boundary_t),  intent(inout) :: forcing
        type(options_t), intent(in)       :: options
        real, intent(in)                  :: dt
    end subroutine

  end interface

end module
