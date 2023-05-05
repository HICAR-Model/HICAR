!-----------------------------------------------------------------------
! Physical constants
!-----------------------------------------------------------------------
module CONSTANTS
real, parameter :: &
  cp = 1005,         &! Specific heat capacity of air (J/K/kg)
  eps = 0.622,       &! Ratio of molecular weights of water and dry air
  e0 = 610.78,       &! Saturation vapour pressure at Tm (Pa)
  grav = 9.81,       &! Acceleration due to gravity (m/s^2)
  hcap_ice = 2100,   &! Specific heat capacity of ice (J/K/kg)
  hcap_wat = 4180,   &! Specific heat capacity of water (J/K/kg)
  hcon_air = 0.025,  &! Thermal conductivity of air (W/m/K)
  hcon_clay = 1.16,  &! Thermal conductivity of clay (W/m/K)
  hcon_ice = 2.24,   &! Thermal conducivity of ice (W/m/K)
  hcon_sand = 1.57,  &! Thermal conductivity of sand (W/m/K)
  hcon_wat = 0.56,   &! Thermal conductivity of water (W/m/K)
  I0 = 1367,         &! Solar constant (W/m^2)
  Lf = 0.334e6,      &! Latent heat of fusion (J/kg)
  Lv = 2.501e6,      &! Latent heat of vapourisation (J/kg)
  Ls = Lf + Lv,      &! Latent heat of sublimation (J/kg)
  pi = 3.14159,      &! pi
  Rair = 287,        &! Gas constant for air (J/K/kg)
  Rwat = 462,        &! Gas constant for water vapour (J/K/kg)
  Runi = 8313,       &! Universal gas constant (J/kmol/K)
  rho_ice = 917,     &! Density of ice (kg/m^3)
  rho_wat = 1000,    &! Density of water (kg/m^3)
  rho_air = 1.275,   &! Density of air (kg/m^3)
  sb = 5.67e-8,      &! Stefan-Boltzmann constant (W/m^2/K^4)
  em_snow = 0.99,    &! Emissivity snow for Stefan-Boltzmann
  em_soil = 0.90,    &! Emissivity soil for Stefan-Boltzmann
  Tm = 273.15,       &! Melting point (K)
  vkman = 0.4,       &! Von Karman constant
  xM = 18.01,        &! Molecular weight of water (kg/kmol)
  visc_air = 13.e-6, &! Kinematic viscosity of air (m^2/s)
  undef = 1.e+6       ! Initialization value for allocatables
integer, parameter :: &
  iundef = 1.e+6      ! Initialization value for integer allocatables
end module CONSTANTS

!-----------------------------------------------------------------------
! Model configuration
!-----------------------------------------------------------------------
module MODCONF
! Process options                                         : Possible values
integer :: &
  ALBEDO,     &! snow albedo                              : 0, 1, 2
  CANMOD,     &! forest canopy                            : 0, 1
  CONDCT,     &! snow thermal conductivity                : 0, 1
  DENSTY,     &! snow density                             : 0, 1, 2, 3
  EXCHNG,     &! turbulent exchange                       : 0, 1
  HYDROL,     &! snow hydraulics                          : 0, 1, 2
  SNFRAC,     &! snow cover fraction                      : 0, 1, 2, 3, 4
  RADSBG,     &! subgrid radiation param                  : 0, 1
! Driving data options
  ZOFFST,     &! Measurement height offset                : 0, 1
  SNTRAN,     &! snow transport                           : 0, 1
  SNSLID,     &! snow slides                              : 0, 1
  SNOLAY,     &! density-dependent layering               : 0, 1
  HISWET,     &! use historical variable of layer wetting : 0, 1
  CHECKS       ! check state variables at every time step : 0, 1, 2
logical :: &
  HN_ON,      &!activate the new snow model
  FOR_HN       !write 18h states for the hn model.
end module MODCONF

!-----------------------------------------------------------------------
! Model output configuration
!-----------------------------------------------------------------------
module MODOUTPUT
! list of diagnostic variables that can be written to output bin files.
! at the moment, only 2d real variables are handled.
character(len=4), dimension(40) :: &
  WRITE_DIAG_VARS = (/'rotc', &  ! Roff      Total runoff, snow and bare soil (kg/m^2)
                      'hsnt', &  ! snowdepth Total snowdepth (m)
                      'swet', &  ! SWE       Total SWE (kg/m^2)
                      'slqt', &  ! Sliq_out  Total LWC (kg/m^2)
                      'swtb', &  ! Sdirt     Incoming direct beam radiation corrected for subgrid topography (W/m^2)
                      'swtd', &  ! Sdift     Incoming diffuse beam radiation corrected for subgrid topography (W/m^2)
                      'lwtr', &  ! Lwt       Incoming longwave radiation corrected for subgrid topography (W/m^2)
                      'romc', &  ! meltflux_out Runoff from snowmelt at the base of snow (kg/m^2)
                      'sbsc', &  ! Sbsrf     Snow sublimation rate (kg/m^2/s)
                      'asrf', &  ! asrf_out  Surface albedo
                      'emlt', &  ! Melt      Surface melt rate (kg/m^2/s)
                      'esrf', &  ! Esrf      Moisture flux from the surface (kg/m^2/s)
                      'eveg', &  ! Eveg      Moisture flux from vegetation (kg/m^2/s)
                      'ghsl', &  ! Gsoil     Heat flux into soil (W/m^2)
                      'hesr', &  ! Hsrf      Sensible heat flux from the surface (W/m^2)
                      'intc', &  ! intcpt    Canopy interception (kg/m^2)
                      'khag', &  ! KH        Eddy diffusivity for heat to the atmosphere (m/s)
                      'khac', &  ! KHa       Eddy diffusivity for heat from the canopy air space (m/s)
                      'khgr', &  ! KHg       Eddy diffusivity for heat from the ground (m/s)
                      'khve', &  ! KHv       Eddy diffusivity for heat from vegetation (m/s)
                      'kwgr', &  ! KWg       Eddy diffusivity for water from the ground (m/s)
                      'kwve', &  ! KWv       Eddy diffusivity for water from vegetation (m/s)
                      'lahe', &  ! LE        Latent heat flux to the atmosphere (W/m^2)
                      'lesr', &  ! LEsrf     Latent heat flux from the surface (W/m^2)
                      'lwsc', &  ! LWsci     Subcanopy incoming longwave radiation (W/m^2)
                      'lwve', &  ! LWveg     Net longwave radiation absorbed by vegetation (W/m^2)
                      'rnet', &  ! Rnet      Net radiation (W/m^2)
                      'rnsr', &  ! Rsrf      Net radiation at surface (W/m^2)
                      'sbve', &  ! Sbveg     Sublimation from vegetation (kg/m^2)
                      'sehe', &  ! H         Sensible heat flux to the atmosphere (W/m^2)
                      'swsc', &  ! SWsci     Subcanopy incoming shortwave radiation (W/m^2)
                      'swsr', &  ! SWsrf     Net SW radiation absorbed by the surface (W/m^2)
                      'swve', &  ! SWveg     Net SW radiation absorbed by vegetation (W/m^2)
                      'uasc', &  ! Usc       Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)
                      'unld', &  ! unload    Snow mass unloaded from canopy (kg/m^2)
                      'sltc', &  ! dm_salt   SWE change due to saltation (kg/m^2)
                      'sspc', &  ! dm_susp   SWE change due to suspension (kg/m^2)
                      'sblc', &  ! dm_subl   SWE change due to sublimation (kg/m^2)
                      'sbgc', &  ! dm_subgrid SWE change due to subgrid redistribution (kg/m^2)
                      'sldc'  &  ! dm_slide  SWE change due to snow slides (kg/m^2)
                      /)
   ! GM Note: Energy-balance relevant diagnostics omitted in this list can be derived from existing variables:
                    !Hveg = H-Hsrf
                    !LEveg = LE-LEsrf
                    !LWsrf = Tss^4*sb
                    !Rveg  = SWveg + LWveg      
! list of state variables that can be written to output bin files.
! at the moment, only 2d real variables are handled.
character(len=4), dimension(14):: &
  WRITE_STATE_VARS = (/&
                  'alse', &  ! albs
                ! 'hsnl', &  ! Ds(3d var)
                  'scfe', &  ! fsnow
                ! 'nsne', &  ! Nsnow
                ! 'sicl', &  ! Sice (3d var)
                ! 'slql', &  ! Sliq (3d var)
                  'tsfe', &  ! Tsrf
                ! 'tsnl', &  ! Tsnow (3d var)
                ! 'tsll', &  ! Tsoil (3d var)
                  'hsmn', &  ! snowdepthmin
                  'hsmx', &  ! snowdepthmax
                ! 'hshs', &  ! snowdepthhist (3d var, 3rd dimension: time.)
                  'swmn', &  ! swemin
                  'swmx', &  ! swemax
                ! 'swhs', &  ! swehist (3d var, 3rd dimension: time.)
                  'qcan', &  ! Qcan
                  'sveg', &  ! Sveg
                  'tcan', &  ! Tcan
                  'tveg', &  ! Tveg
                  'sblt', &  ! dm_tot_subl
                  'trat', &  ! dm_tot_trans
                  'sldt'  &  ! dm_tot_slide
                  /)
character(len=4), dimension(40) :: &
  LIST_DIAG_RESULTS  ! List of diagnostic variables that the user wants to write into output bin files (subset of WRITE_DIAG_VARS)
character(len=4), dimension(14) :: &
  LIST_STATE_RESULTS ! List of result variables that the user wants to write into output bin files (subset of WRITE_STATE_VARS)
logical, dimension(40)  :: WRITE_DIAG_TABLE  ! table specifying whether the corresponding variables in WRITE_DIAG_VARS should be written or not.
logical, dimension(14) :: WRITE_STATE_TABLE ! table specifying whether the corresponding variables in WRITE_STATE_VARS should be written or not.
end module MODOUTPUT
!-----------------------------------------------------------------------
! Output diagnostics
!-----------------------------------------------------------------------
module DIAGNOSTICS
integer :: &
  Nave                ! Number of timesteps in average outputs
end module DIAGNOSTICS

!-----------------------------------------------------------------------
! Meteorological driving variables
!-----------------------------------------------------------------------
module DRIVING
integer :: &
  year,              &! Year
  month,             &! Month of year
  day                 ! Day of month
real :: &
  hour                ! Hour of day
real :: &
  dt,                &! Timestep (s)
  zT,                &! Temperature measurement height (m)
  zU,                &! Wind speed measurement height (m)
  zRH                 ! Relative humidity measurement height (m)
real, allocatable :: &
  LW(:,:),           &! Incoming longwave radiation (W/m^2)
  Ps(:,:),           &! Surface pressure (Pa)
  Qa(:,:),           &! Specific humidity (kg/kg)
  RH(:,:),           &! Relative humidity (%)
  Rf(:,:),           &! Rainfall rate (kg/m^2/s)
  Sf(:,:),           &! Snowfall rate (kg/m^2/s)
  Sf24h(:,:),        &! Snowfall 24hr (kg/m^2)
  Sdif(:,:),         &! Diffuse shortwave radiation (W/m^2)
  Sdir(:,:),         &! Direct-beam shortwave radiation (W/m^2)
  Ta(:,:),           &! Air temperature (K)
  Tv(:,:),           &! Time-varying transmissivity for direct SWR (-)
  Ua(:,:),           &! Wind speed (m/s)
  Udir(:,:)           ! Wind direction (degrees, clockwise from N)
end module DRIVING

!-----------------------------------------------------------------------
! Grid parameters
!-----------------------------------------------------------------------
module GRID
integer :: &
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions
real :: &
  Ds_min,            &! Minimum possible snow layer thickness (m)
  Ds_surflay          ! Maximum thickness of surface fine snow layering (m)
real, allocatable :: &
  Dzsnow(:),         &! Minimum snow layer thicknesses (m)
  Dzsoil(:)           ! Soil layer thicknesses (m)
end module GRID

!-----------------------------------------------------------------------
! Input / output unit numbers
!-----------------------------------------------------------------------
! module IOUNITS
! integer, parameter :: &
!   umta = 51           ! Metadata output file unit number
! end module IOUNITS

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
module PARAMETERS
! Numerical solution parameter
integer :: &
  Nt,                &! Number of iterations of physics within an input timestep
                      ! (Nt = 4 implies a 15 min timestep for an hourly input timestep)
  Nitr                ! Number of iterations in energy balance calulation
! Vegetation parameters
real :: &
  avg0,              &! Snow-free vegetation albedo
  avgs,              &! Snow-covered vegetation albedo
  cden,              &! Dense canopy turbulent transfer coefficient
  cvai,              &! Canopy snow capacity per unit VAI (kg/m^2)
  Gcn1,              &! Leaf angle distribution parameter
  Gcn2,              &! Leaf angle distribution parameter
  gsnf,              &! Snow-free vegetation moisture conductance (m/s)
  kdif,              &! Diffuse radiation extinction coefficient
  kveg,              &! Canopy cover coefficient
  cveg,              &! Vegetation turbulent transfer coefficient
  rchd,              &! Ratio of displacement height to canopy height
  rchz,              &! Ratio of roughness length to canopy height
  tcnc,              &! Canopy unloading time scale for cold snow (s)
  tcnm                ! Canopy unloading time scale for melting snow (s)
! Snow parameters
real :: &
  a_eta,             &! Temperature factor for Crocus B92 compaction (K^-1)
  asmx,              &! Maximum albedo for fresh snow
  asmn,              &! Minimum albedo for melting snow
  b_eta,             &! First density factor for Crocus B92 compaction (m^3/kg)
  bstb,              &! Atmospheric stability parameter
  bthr,              &! Snow thermal conductivity exponent
  c_eta,             &! Second density factor for Crocus B92 compaction (kg/m^3)
  eta0,              &! Reference snow viscosity (Pa s)
  eta1,              &! Reference snow viscosity for Crocus B92 compaction (Pa s)
  hfsn,              &! Snowcover fraction depth scale (m)
  kfix,              &! Fixed thermal conductivity of snow (W/m/K)
  rgr0,              &! Fresh snow grain radius (m)
  rho0,              &! Fixed snow density (kg/m^3)
  rhob,              &! Temperature factor in fresh snow density (kg/m^3/K)
  rhoc,              &! Wind factor in fresh snow density (kg s^0.5/m^3.5)
  rhof,              &! Fresh snow density (kg/m^3)
  rhos_min,          &! Minimum snow density (kg/m^3)
  rhos_max,          &! Maximum snow density (kg/m^3)
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  Salb,              &! Snowfall to refresh albedo (kg/m^2)
  snda,              &! Thermal metamorphism parameter (1/s)
  Talb,              &! Albedo decay temperature threshold (C)
  tcld,              &! Cold snow albedo decay time scale (s)
  tmlt,              &! Melting snow albedo decay time scale (s)
  trho,              &! Snow compaction time scale (s)
  Wirr,              &! Irreducible liquid water content of snow
  z0sn                ! Snow roughness length (m)
! Surface parameters
real :: &
  gsat                ! Surface conductance for saturated soil (m/s)
! Additional parameters used in forest snow process parametrizations
real :: &      
  adfs,              &! Snow albedo adjustment dependent on SWR
  adfl,              &! Snow albedo adjustment dependent on LWR
  fsar,              &! Snow albedo adjustment range dependent on vegetation fraction
  psf,               &! Scaling factor for solid precipitation (within forest stand, at min CC)
  psr,               &! Range of solid precipitation (within forest stand, spread min-max CC)
  wcan,              &! Parameter of exponential wind profile 
  zsub,              &! Sub-canopy reference height (m)
  zgf,               &! Roughness length adjustment factor depending on vegetation fraction  
  zgr,               &! Roughness length adjustment range depending on vegetation fraction  
  khcf                ! Diffusivity adjustment for canopy effects (Finnigan 2000)
! Parameters characterizing grid cell tiles
real :: & 
  fthresh             ! Forest fraction required for forest tile to be considered 
end module PARAMETERS

!-----------------------------------------------------------------------
! Spatial surface characteristics
!-----------------------------------------------------------------------
module PARAMMAPS
real, allocatable :: &
  alb0(:,:),         &! Snow-free ground albedo
  canh(:,:),         &! Canopy heat capacity (J/K/m^2)
  fcly(:,:),         &! Soil clay fraction
  fsnd(:,:),         &! Soil sand fraction
  fsky(:,:),         &! Sky view fraction
  fves(:,:),         &! Stand-scale canopy cover fraction
  scap(:,:),         &! Canopy snow capacity (kg/m^2)
  trcn(:,:),         &! Canopy transmissivity
  VAI(:,:),          &! Vegetation area index
  z0sf(:,:),         &! Snow-free roughness length (m)
  vegsnowd_xy(:,:)    ! Vegetation snow holding capacity (m)
end module PARAMMAPS

!-----------------------------------------------------------------------
! Soil properties
!-----------------------------------------------------------------------
module SOILPARAMS
real, allocatable :: &
  b(:,:),            &! Clapp-Hornberger exponent
  hcap_soil(:,:),    &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil(:,:),    &! Thermal conductivity of dry soil (W/m/K)
  sathh(:,:),        &! Saturated soil water pressure (m)
  Vcrit(:,:),        &! Volumetric soil moisture at critical point
  Vsat(:,:)           ! Volumetric soil moisture at saturation
end module SOILPARAMS

!-----------------------------------------------------------------------
! Model state variables  
!-----------------------------------------------------------------------
module STATE_VARIABLES
! First iteration identifier
integer, allocatable :: &
  firstit             ! First iteration identifier

! Canopy properties
real, allocatable :: &
  Qcan(:,:),         &! Canopy air space humidity
  Tcan(:,:),         &! Canopy air space temperature (K)
  Sveg(:,:),         &! Snow mass on vegetation (kg/m^2)
  Tveg(:,:)           ! Vegetation temperature (K)
  
! Surface state variables
real, allocatable :: &
  Tsrf(:,:),         &! Surface skin temperature (K)
  fsnow(:,:)          ! Snow cover fraction terrain

! Snow state variables
integer, allocatable :: &
  Nsnow(:,:)          ! Number of snow layers
real, allocatable ::     &
  albs(:,:),             &! Snow albedo
  Ds(:,:,:),             &! Snow layer thicknesses (m)
  histowet(:,:,:),       &! Historical variable for past wetting of a layer (0-1)
  rgrn(:,:,:),           &! Snow layer grain radius (m)
  Sice(:,:,:),           &! Ice content of snow layers (kg/m^2)
  Sliq(:,:,:),           &! Liquid content of snow layers (kg/m^2)
  Tsnow(:,:,:),          &! Snow layer temperatures (K)
  swehist(:,:,:),        &! History of SWE during last 14 days (kg/m^2). Most recent entries first.
  swemin(:,:),           &! Minimum SWE during the season (m)
  swemax(:,:),           &! Maximum SWE during the season (m)
  snowdepthhist(:,:,:),  &! History of snow depth during last 14 days (m). Most recent entries first.
  snowdepthmin(:,:),     &! Minimum snow depth at time step of swemin (m)
  snowdepthmax(:,:)       ! Maximum snow depth at time stemp of swemax(m)

! Soil state variables
real, allocatable :: &
  theta(:,:,:),          &! Volumetric moisture content of soil layers
  Tsoil(:,:,:)            ! Soil layer temperatures (K)
  
! SnowTran3D state variables
real, allocatable ::     &
  dm_tot_subl(:,:),      &! Cumulated SWE change due to sublimation (kg/m^2)
  dm_tot_trans(:,:)       ! Cumulated transported SWE (kg/m^2)
  
! SnowSlide state variables
real, allocatable ::     &
  dm_tot_slide(:,:)       ! Cumulated SWE change due to snow slides (kg/m^2)
  
integer, allocatable ::     &
  index_grid_dem_sorted(:,:)   ! Location (i,j) of sorted grid points
  
end module STATE_VARIABLES

!-----------------------------------------------------------------------
! Landuse information 
!-----------------------------------------------------------------------
module LANDUSE
! Canopy properties
real, allocatable :: &
  fveg(:,:),         &! Canopy cover fraction
  lai(:,:),          &! Leaf area index 
  vfhp(:,:),         &! Hemispherical sky-view fraction including canopy
  hcan(:,:)           ! Canopy height (m)

! Terrain properties
real, allocatable :: &
  fsky_terr(:,:),    &! Sky view fraction terrain
  slopemu(:,:),      &! slope parameter 
  xi(:,:),           &! terrain correlation length
  Ld(:,:),           &! grid cell size or domain size (m)
  lat(:,:),          &! latitude of each grid cell (center?)
  lon(:,:),          &! longitude of each grid cell (center?) 
  dem(:,:),          &! grid elevation (m)
  slope(:,:),        &! slope (deg)
  Shd(:,:),          &! snow holding depth, for SnowSlide (m)
  forest(:,:)         ! grid cell forest fraction
  
! Terrain properties
real, allocatable :: &
  dem_sorted(:)      ! sorted grid elevation (m)
  
end module LANDUSE

!-----------------------------------------------------------------------
! Parameters for SnowTran3D
!-----------------------------------------------------------------------
module PARAM_SNOWTRAN3D
real, parameter :: &

  ! Define whether the threshold surface shear velocity will be
  ! constant during the simulation (Utau_t_flag = 0.0), or whether
  ! it will evolve as a funtion of air temperature and wind speed
  ! (Utau_t_flag = 1.0).
  !!! (recommended default value: Utau_t_flag = 0.0)
  Utau_t_flag = 1.0,    &! Flag for variable threshold friction velocitiy (1.0) or constant (0.0)

  ! For the case of Utau_t_flag = 0.0, define what the threshold
  ! surface shear velocity (m/s) will be during the simulation.
  !!! (recommended default value: Utau_t_const = 0.25)
  Utau_t_const = 0.25,  &! Constant threshold friction velocity (m/s) used if Utau_t_flag = 0.0

  ! Define whether you want to apply the Tabler subgrid algorithm.
  ! Flag = 0.0 is the original model (no subgrid redistribution, as
  ! in the original model), flag = 1.0 uses Tabler (1975) to define
  ! equilibrium profiles for topographic drift catchments.  Note
  ! that the model has not been tested for grid increments less than
  ! 1.0 m deltax and deltay.  Also note that this subgrid algorithm
  ! is not really appropriate for grid increments greater than
  ! 30.0 m.  There is also a computational penalty for running with
  ! this turned on: the simulations can take as much as 35% longer
  ! when subgrid_flag = 1.0 (but it does produce more realistic
  ! drift profiles at the higher resolutions).  If subgrid_flag = 1.0,
  ! you need to set topoflag = 0.0.
  !!! (recommended default value: subgrid_flag = 0.0)
  subgrid_flag = 0.0,   &! Flag for subgrid parameterization

  ! The implementation of Tabler (1975) allows the user to define
  ! the area immediately upwind of any Tabler drift trap to be an
  ! erosion area where no snow is allowed to accumulate.  This area
  ! is defined in meters (the distance upwind of the drift trap lip
  ! that will not be allowed to accumulate snow).  A value of 1.0 m
  ! or greater is required to simulate realistic drift profiles.
  !!! (recommended default value: erosion_dist = 1.0)
  erosion_dist = 1.0,   &! Erosion area (m)

  ! Define the Tabler snow accumulation surfaces to include the drift
  ! traps and some fraction of the initial snow depth and
  ! precipitation inputs.  The "tabler precipitation scale",
  ! tp_scale, defines the fraction of the initial snow depth and
  ! precipitation that is added to the Tabler Surfaces.  This is
  ! required so that precipitation accumulating on flat ground, for
  ! example, is not removed because it is higher than the Tabler
  ! Surface.  Valid numbers are between 0.0 and 1.0, but values
  ! greater than 0.9 typically lead to unrealistic drift profiles.
  ! Values closer to 0.0 lead to drift profiles that equal the
  ! Tabler equilibrium drift profiles.  In general, the higher the
  ! tp_scale value, the rougher the drift profile and the more
  ! sensitive the simulation is to the model wind parameters like
  ! curve_len_scale and slopewt and curvewt.
  !!! (recommended default value: tp_scale = 0.9)
  tp_scale = 0.9,       &! Tabler precipitation scale

  ! For subgrid_flag = 1.0, you have to define the dominant wind
  ! direction (one of the 8 primary wind directions) you want
  ! Tabler surfaces and drifts to be generated for.
  !!! (recommended default value: user input required)
  !tabler_dir = 0.0,    &! Dominant wind direction N
  !tabler_dir = 45.0,   &! Dominant wind direction NE
  !tabler_dir = 90.0,   &! Dominant wind direction E
  !tabler_dir = 135.0,  &! Dominant wind direction SE
  !tabler_dir = 180.0,  &! Dominant wind direction S
  tabler_dir = 225.0,  &! Dominant wind direction SW
  !tabler_dir = 270.0,  &! Dominant wind direction W
  !tabler_dir = 315.0,   &! Dominant wind direction NW

  ! You can make adjustments to the slopes of the Tabler surfaces
  ! using the 'slope_adjust' parameter.  This factor is multiplied by
  ! the calculated Tabler equilibrium drift surface slope to get the
  ! simulated slope.  So, for example, if you want the simulated
  ! slope to be 50% of what Tabler calculated, use a value of 0.5; a
  ! value of 1.5 would give half again steeper than Tabler's original
  ! model.
  !!! (recommended default value: slope_adjust = 1.0)
  slope_adjust = 1.0,   &! Slope adjustment factor

  ! Define whether you want two-layer accounting for soft (movable)
  ! and hard (unmovable) layers (flag = 0.0 = off, flag = 1.0 = on).
  !!! (recommended default value: twolayer_flag = 1.0)
  twolayer_flag = 1.0,  &! Flag for soft/hard layers distinction

  ! Define whether you want the upwind boundaries to have zero
  ! incoming transport flux, or the equilibrium transport flux.
  ! bc_flag = 0.0 corresponds to a zero flux, and bc_flag = 1.0
  ! corresponds to the equilibrium flux.
  !!! (recommended default value: bc_flag = 0.0)
  bc_flag = 0.0,        &! Boundary condition flag

  ! Density of snow deposited by snowdrift
  !!! (recommended default value: rho_snow = 300.0)
  rho_snow = 300.0       ! Constant snow density (kg/m^3)
  
end module PARAM_SNOWTRAN3D

!-----------------------------------------------------------------------
! Constants for SnowTran3D
!-----------------------------------------------------------------------
module CONSTANTS_SNOWTRAN3D
real, parameter :: &

  ! Constants related to surface shear stress and saltation transport.
  fetch = 500.0,        &! Equilibrium fetch distance (m)
  xmu = 3.0,            &! Scaling constant for non-equilibrium saltation transport
  C_z = 0.12,           &! Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500
  h_const = 1.6,        &! Coefficient 1.6 in Liston and Sturm (1998) eq. 14 p. 501
  wind_min = 4.0,       &! Minimum wind speed to compute snow transport (m/s)
  
  ! Constants related to suspended snow profile.
  Up_const = 2.8,       &! Constant coefficient for calculation of U_p
  dz_susp = 0.20,       &! dz in suspension layer (m)
  ztop_susp = 2.0,      &! Height of the top of suspension layer (m)
  fall_vel = 0.3,       &! Particle-settling velocity, s in L&S 1998 eq. 12 (m/s)
  Ur_const = 0.5         ! Constant coefficient beta in L&S 1998, eq. 13
  
end module CONSTANTS_SNOWTRAN3D

!-----------------------------------------------------------------------
! Parameters for SnowSlide
!-----------------------------------------------------------------------
module PARAM_SNOWSLIDE
real, parameter :: &

  ! Dynamic ratio for the snow holding depth during an avalanche
  !!! (recommended default value: dyn_ratio = 1.0)
  dyn_ratio = 1.0,  &! Dynamic snow holding depth ratio

  ! Density of snow deposits by avalanches
  !!! (recommended default value: rho_snow = 300.0)
  rho_deposit = 300.0,  &! Constant snow avalanche deposit density (kg/m^3)  

  ! Snow slide can only be triggered above a certain slope
  !!! (recommended default value: slope_min = 25 deg)
  slope_min = 25.0,     &! Minimum slope for snow slide occurrence (deg)
  
  ! Snow holding depth has a min threshold 
  !!! (recommended default value: Shd_min = 0.05 m)
  Shd_min = 0.05         ! Minimum snow holding depth (m)
  
end module PARAM_SNOWSLIDE
