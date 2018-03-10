##ICAR settings file documentation

###Main Namelists
REQUIRED:
 * model_version
 * physics
 * parameters
 * files_list
 * z_info (technically optional)
 
 OPTIONAL:
 * restart_info
 * lt_parameters
 * lsm_parameters
 * mp_parameters
 * adv_parameters


###model_version
This namelist specifies the model version associated with the settings file, and provides a place to enter a comment into the output files. 

###physics
This namelist specifies which physics options will be used when the model runs. 

###parameters
This namelist specifies various model settings. 

###files_list
This namelist specifies the name of various files that will be used by ICAR.  The most important files are the forcing data files, and the ICAR hi-res grid file (akin to the WRF geogrid).  Additional files can be specified for spatially variable calibration input, and a file to list the input forcing file names. 

###z_info
This namelist specifies the thickness of each model layer.  Technically this is optional as it will default to something reasonable, but it is recommended that it be set as the default uses too fine a discretization (and thus runs slower) than should be used. 

###restart_info
This optional namelist specifies the necessary input for a restart run. 

###lt_parameters
This optional namelist specifies various parameters that are used in the linear theory calculations. 

###lsm_parameters
This optional namelist specifies various parameters that are used in the land surface model. 

###mp_parameters
This optional namelist specifies various parameters that are used in the microphysics scheme. 

###adv_parameters
This optional namelist specifies various parameters that are used in the advection calculations. 


<!-- !---------------------------------------------------------
!   Model and run meta-data
!---------------------------------------------------------
&model_version
    version = "0.9.3",                    ! This must match the version of the compiled code
    comment = "Add your comment here"     ! This will be stored in output files
/

!---------------------------------------------------------
!   Model levels specification (may be optional, but should be specified)
!---------------------------------------------------------
&z_info
    !   Sample model level thickness [m]  Bottom levels could be thicker.
    dz_levels = 50.,   75.,  125.,  200.,  300.,  400.,  500.,  500.,  500.,  500.,    !  1-10
               500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,    ! 10-20
               500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,    ! 20-30
               500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.     ! 30-40

    !   If you want to line up model level with common forcing data
    !   ERAi levels
    !dz_levels= 24.8,  36.5,  51.8,  70.1,  90.8, 113.5, 137.9, 163.7, 190.5, 218.1,   !  1-10
    !          246.4, 275.1, 304.3, 333.6, 363.0, 392.4, 421.7, 450.8, 479.6, 508.0,   ! 10-20
    !          535.9, 563.2, 589.8, 615.7, 640.9, 665.5, 689.8, 714.1, 739.4, 767.2,   ! 20-30
    !          796.8, 826.6, 856.2, 885.1, 912.5, 937.9, 961.4, 979.4, 990.1, 976.6    ! 30-40
    !   WRF levels from Headwaters 36km runs
    !dz_levels= 36.,   51.,   58.,   73.,   74.,  111.,  113.,  152.,  155.,  157.,    !  1-10
    !          160.,  245.,  251.,  258.,  265.,  365.,  379.,  395.,  413.,  432.,    ! 10-20
    !          453.,  476.,  503.,  533.,  422.,  443.,  467.,  326.,  339.,  353.,    ! 20-30
    !          369.,  386.,  405.,  426.,  450.,  477.,  455.,  429.,  396.,  357.     ! 30-40
/

!---------------------------------------------------------
!   Specify physics options to use for the model run
!---------------------------------------------------------
&physics
    ! Common precipitation downscaling run use pbl=0 lsm=0 mp=1 rad=0 conv=0 adv=1 wind=1
    ! For a FASTER run (simpler physics), set mp=2
    ! If surface air temperature is important use pbl=2 lsm=3 rad=2 water=2 this requires Noah LSM data
    ! N/A = Not Available or Not fully implemented
    ! wishlist = No Code Present yet

    pbl = 0,  ! 1=legacy (deprecated)      2=Simple (Local HP96)        3=YSU             (N/A)
    lsm = 0,  ! 1=use prescribed fluxes    2=Simple LSM (N/A)           3=Noah LSM
    water=2,  ! 1=use prescribed (w/lsm=1) 2=Simple sea surface fluxes
    mp  = 1,  ! 1=Thompson                 2=Simple (SB04)              3=Morrison        (wishlist)
    rad = 0,  ! 1=use prescribed fluxes    2=Simple (empirical)         3=RRTMG           (wishlist)
    conv= 0,  ! 1=Tiedke Scheme            2=Simple Scheme (wishlist)   3=Kain-Fritsch
    adv = 1,  ! 1=Upwind                   2=MPDATA                     3=Adams-Bashforth (wishlist)
    wind= 1   ! 1=Linear Theory            2=INFORM style (wishlist)    3=Dynamical?      (wishlist)
/

!---------------------------------------------------------
!   Files to be used by the run
!---------------------------------------------------------
&files_list
    !   This is the high-resolution input filename
    !   primary inputs from this file are lat, lon, and terrain, optionally soil and veg types
    init_conditions_file="baseline/geo_4km_conus.nc",

    !   This is the prefix for all output files (any directories must be created prior to running)
    output_file="output/icar_out_",

    !   This is a list of the boundary conditions files number of files must match nfiles variable above
    boundary_files= "forcing/wrfout_d01_2001-04-01_03:00:00", "forcing/wrfout_d01_2001-06-30_03:00:00"

    !   Alternatively a separate file containing one forcing file name per line may be specified
    !   This file may be generated by : 
    !       ls -1 forcing/* | sed 's/$/"/g;s/^/"/g'>file_list.txt
    !   sed is used to add " around the filename. 
    !   The quotes are probably only necessary if there are special characters or spaces on the line
    ! forcing_file_list = "file_list.txt"

    !   Files to read "calibration" data from
    ! nsq_calibration_file = "nsq_calibration.nc",
    ! linear_mask_file = "linear_weights.nc"
/

!---------------------------------------------------------
!   Main List of Parameters
!---------------------------------------------------------
&parameters
    !   Set this to the starting date of the first low-resolution forcing file
    forcing_start_date = '2001-04-01 03:00:00',
    !   Set this to the date to start running the model (defaults to the forcing_start_date)
    start_date = "2001-04-02 00:00:00",
    !   Set this to the date to stop running the model 
    end_date = "2001-04-10 00:00:00",
    !   Calendar used by the forcing data "gregorian", "standard", "noleap", "365-day", "360-day"
    calendar = "standard",

    !   The length of an input forcing time step
    inputinterval = 3600,   ! [s]
    !   The output interval
    outputinterval = 3600,  ! [s]

    !   Limit output data to near surface variables
    !   WARNING if true it is impossible to restart the run (for now)
    ! surface_io_only = False,

    !   The grid spacing of the high-resolution data
    dx = 4000.0,        ! [m]
    !   The approximate grid spacing of the forcing data
    !   only used in rm_linear_winds?
    ! dxlow = 20000.0,    ! [m]

    !   Read dz from the namelist file (below)
    readdz = True,

    !   The number of vertical levels to run (suggest ~10-30 levels with a model top around 4-8km)
    !   this is now optional, if not supplied, ICAR will determine it from the number of levels specified
    !   if it is supplied it must be less than or equal to the number of levels specified below
    !   but it can be used to subset the number of levels used.
    nz = 15, ! []
    !   Set this to true of the zvar in the input data is actually in units of geopotential height (m/s^2)
    z_is_geopotential = False,
    !   Specify that the height of the forcing data will change through the simulation (common for atmospheric model-level output)
    time_varying_z = True,
    !   Use height above ground layer to interpolate the wind field instead of height above sea level.
    use_agl_height = False,
    
    !   Multiplier on CFL number to increase stability if unstable (make it <1 to increase stability, >1 to increase model speed)
    ! cfl_reduction_factor = 1.0
    !   CFL method 1 = max(1D winds) * sqrt(3), 2=max(1D,ave.3D wind)*sqrt(3), 3=max(sum.3D wind), 4=max(sum.3D wind)*sqrt(3), 5 = sum(max.3d)
    !   Note that 4 is probably the safest, but 3 has always been stable and is left as the default. 
    !   5 is the value that used to be used. 
    !   Simulations with 4 will run 1.7x slower. 
    ! cfl_strictness = 3

    !   If the forcing data come from WRF, the temperature data probably have an offset applied
    !   t_offset will be added to the forcing temperature data.  Defaults to 0
    ! t_offset = 300, ! [K]

    !   Distance to smooth winds over [m] ~100000 is reasonable
    !   larger values result in less large scale convergence/divergence in the flow field
    !   smaller value result in more and can destroy orographic precip and result in odd spatial coherence
    !   depending on the forcing data resolution. At a minimum, this should be ~dxlow
    smooth_wind_distance = 72000, ! [m]

    !   To run an ideal simulation in which the boundary conditions are held constant
    ! ideal = false,
    !   To use an externally supplied high-resolution wind field (ignore)
    ! external_winds = false,
    !   Number of external wind files (ignore)
    ! n_ext_winds = 1,  ! [n-files]
    !   Run with a horizontally averaged wind field
    ! mean_winds = false,
    !   Run with a horizontally averaged boundary conditions
    ! mean_fields = false,

    !   Use this to restart the model restart_info must be supplied below
    restart = false,

    !   Use density in the advection step (violates linear theory assumptions)
    advect_density = false,

    !   The number of grid cells to remove from all sides of the high-resolution grid
    !   used primarily for faster test runs over a smaller domain
    ! buffer = 0,   ! [n-gridcells]

    !   Doesn't do much at the moment, increases output print at runtime
    debug = true,
    warning_level = 4, ! 0-10 increases the level of errors it warns about and quits over (slightly)

    !   If the following are true, their respective namelists (below) will also be read in. 
    !   Read parameters for advection
    use_adv_options = true,
    !   Read parameters for linear theory
    use_lt_options = true,
    !   Read parameters for microphysics (thompson only at this point)
    use_mp_options = true
    !   Read parameters for land surface model
    use_lsm_options = true,
/



!---------------------------------------------------------
!   Specification of variable names in input files
!---------------------------------------------------------
&var_list
    ! These are the names of the variables in the forcing data files
    ! variables on the mass / center grid
    pvar    = "P",          ! pressure                  [Pa]
    pbvar   = "PB",         ! base pressure state       [Pa]        OPTIONAL
    tvar    = "T",          ! temperature               [K]   (with optional offset)
    qvvar   = "QVAPOR",     ! water vapor mixing ratio  [kg/kg]
    qcvar   = "QCLOUD",     ! cloud water mixing ratio  [kg/kg]     OPTIONAL
    qivar   = "QICE",       ! cloud ice mixing ratio    [kg/kg]     OPTIONAL
    hgtvar  = "HGT",        ! surface elevation         [m]
    zvar    = "PH",         ! model level elevations    [m or m/s^2 if z_is_geopotential]
    zbvar   = "PHB",        ! base height state         [m or m/s^2] OPTIONAL
    latvar  = "XLAT",       ! latitude                  [degrees]
    lonvar  = "XLONG",      ! longitude                 [degrees]
    sst_var = "TSK"         ! Water surface temperature [K]          OPTIONAL (used with water=2)

    ! variables on the ew staggered (U) grid
    uvar    = "U",          ! East-West wind speed      [m/s]
    ulat    = "XLAT_U",     ! latitude                  [degrees]
    ulon    = "XLONG_U",    ! longitude                 [degrees]

    ! variables on the NS staggered (V) grid
    vvar    = "V",          ! North-South wind speed    [m/s]
    vlat    = "XLAT_V",     ! latitude                  [degrees]
    vlon    = "XLONG_V",    ! longitude                 [degrees]

    ! these are only used with lsm=1 (pbl should also be >0)
    ! shvar = "HFX",        ! sensible heat flux        [W/m^2]
    ! lhvar = "LH",         ! latent heat flux          [W/m^2]

    ! for lsm=1,pbl=1
    ! pblhvar = "PBLH",     ! Planetary boundary layer height [m]

    ! Radiative fluxes at the surface required with physics:rad=1
    swdown_var = "SWDOWN",  ! Shortwave down            [W/m^2]
    lwdown_var = "GLW",     ! Longwave down             [W/m^2]

    ! only required for some physics code (Noah LSM, water, Tiedke, KF(?))
    landvar = "LANDMASK",   ! land-water mask (as in WRF) 1=land, 0 or 2=water

    ! NOTE, these variables should be in the high-resolution initial conditions netcdf file
    lat_hi  = "XLAT_M",     ! latitude  (mass grid)         [degrees]
    lon_hi  = "XLONG_M",    ! longitude (mass grid)         [degrees]
    ulat_hi = "XLAT_U",     ! latitude  (ew-staggered grid) [degrees]
    ulon_hi = "XLONG_U",    ! longitude (ew-staggered grid) [degrees]
    vlat_hi = "XLAT_V",     ! latitude  (ns-staggered grid) [degrees]
    vlon_hi = "XLONG_V",    ! longitude (ns-staggered grid) [degrees]
    hgt_hi  = "HGT_M"       ! surface elevation             [m]
    
    ! to use the Noah LSM the following fields should also be specified on the high-res grid
    ! vegtype_var    = "IVGTYP",    ! vegetation type index (classification to match VEGPARM.TBL file)
    ! vegfrac_var    = "VEGFRA",    ! vegetation cover fraction
    ! soiltype_var   = "ISLTYP",    ! soil type index (classification to match SOILPARM.TBL file)
    ! soil_deept_var = "SOILTEMP",  ! deep soil temperature         [K]
                                    ! if soil_t_var is not specified this is used
                                    ! throughout the soil column, not just at the bottom.
    ! soil_t_var   = "TSLB",        ! soil temperature (4 levels)   [K]
    ! soil_vwc_var = "SMOIS",       ! soil water content (4 levels) [m^3/m^3]

    ! variables to read from calibration files, both default to "data"
    ! nsq_calibration_var = "data",
    ! linear_mask_var = "data"
/


!---------------------------------------------------------
!   Optionally specified Microphysics parameters (mostly for Thompson)
!---------------------------------------------------------
&mp_parameters
    update_interval = 60 ! maximum update interval allowed
                         ! MP only updated when this interval will be exceeded in the next step

    Nt_c  = 100.e6      !  50, 100,500,1000
    TNO   = 5.0         !  0.5, 5, 50
    am_s  = 0.069       ! 0.052 (Heymsfield), 0.02 (Mitchell), 0.01.
                        ! Note that these values are converted to mks units. Was given as cgs units in Morrison p3 code
    rho_g = 500.0       ! 800, 500, 200
    av_s  = 40.0        ! 11.72 (Locatelli and Hobbs)
    bv_s  = 0.55        ! 0.41
    fv_s  = 100.0       ! 0
    av_g  = 442.0       ! 19.3   from "Cloud-Resolving Modelling of Convective Processes, by Gao and Li,
    bv_g  = 0.89        ! 0.37
    av_i  = 1847.5      ! 700 (Ikawa and Saito)
    Ef_si = 0.05
    Ef_rs = 0.95        ! 1
    Ef_rg = 0.75        ! 1
    Ef_ri = 0.95        ! 1
    C_cubes = 0.5       ! 0.25 Based on Thesis paper "Validation and Improvements of Simulated
                        !      Cloud Microphysics and Orographic Precipitation over the Pacific Northwest"
    C_sqrd  = 0.3
    mu_r    = 0.        ! 1, 2, 5
    t_adjust= 0.0       ! -5, 10, 15
    Ef_rw_l = .False.   ! True sets ef_rw = 1, insted of max 0.95
    Ef_sw_l = .False.   ! True sets ef_rw = 1, insted of max 0.95

    top_mp_level = 0    ! if <=0 just use the actual model top
    local_precip_fraction = 1.0 ! Fraction of micrphysics derived precipitation to deposit in the local grid cell
                                ! the remaining precip is distributed to the surrounding grid cells.
/

!---------------------------------------------------------
!   Optionally specified advection parameters (only used by MPDATA right now)
!---------------------------------------------------------
&adv_parameters
    flux_corrected_transport = true ! Use a flux correction in the transport calculations to prevent ringing and overshoots
                                    ! this should keep MPDATA stable enough for use with the linear winds

    mpdata_order = 2                ! Int: Closure order to use (IORD in MPDATA papers)
                                    ! order=1 equivalent to simple upwind
                                    ! order=2 is standard MPDATA
                                    ! order>2 is a higher order correction that will be very expensive with relatively little gain

    boundary_buffer = False         ! smooth a one grid cell buffer around the boundary
                                    ! to avoid ringing artifacts in non-flux-corrected advection
                                    ! better just to use flux correction as it may crash without it.
/

!---------------------------------------------------------
!   Optionally specified land surface model parameters (mostly for Noah)
!---------------------------------------------------------
&lsm_parameters
    update_interval = 600             ! Int : Seconds to wait before updating land surface fluxes again (default=300)

    LU_Categories = "MODIFIED_IGBP_MODIS_NOAH"   ! Land Use Category definitions
                                    ! Note, this must match a category in VEGPARM.TBL and correspond to
                                    ! the values stored in vegtype_var in the hi-res input var (default="MODIFIED_IGBP_MODIS_NOAH")
                                    ! common values are USGS, USGS-RUC, MODI-RUC, and NLCD40

    monthly_vegfrac = true            ! read / use a 12 month phenology of vegetation fraction

    ! These all default to values defined in common LU_Categories
    ! urban_category = -1             ! Int: index that defines the urban category in LU_Categories
    ! ice_category   = -1             ! Int: index that defines the ice category in LU_Categories
    ! water_category = -1             ! Int: index that defines the water category in LU_Categories
/


!---------------------------------------------------------
!   Optionally specified Linear Theory parameters
!---------------------------------------------------------
&lt_parameters
    buffer = 50                     ! The number of grid cells of buffer to use around the topography for the fft calculations
    stability_window_size = 2       ! The number of grid cells in all directions to average Nsq over for variable_N
    vert_smooth = 2,                ! The number of vertical levels to look up and down when calculating brunt vaisalla frequency
    max_stability = 6e-4            ! The maximum Brunt-Vaisalla frequency to allow
    min_stability = 1e-7            ! The minimum Brunt-Vaisalla frequency to allow

    ! If you want to run with a constant BV instead of a time varying one, it can be set here (and set variable_N to false)
    ! NOTE this will be used for the dry BV,  moist will be dry/10
    ! N_squared = 3.0e-5            ! set this to use a fixed brunt-vaisalla frequency in linear wind calculations
    variable_N = true,              ! use a time varying Nsq (e.g. calculate it from the data don't use the above fixed value)
    linear_update_fraction = 0.5    ! set this to the fraction of the current linear calculation to add to a time-varying perturbation
                                    ! setting to 1 means that waves form instantly, setting it to 0 means they will never form
                                    ! anything in between is the contribution from the current input forcing time step thus it should 
                                    ! change if inputinterval changes. 

    ! linear_contribution = 1.0,    ! set this to the fraction of the linear perturbation you wish to use (1.0 = full/standard linear field)
    spatial_linear_fields = true,   ! use a spatially variable wind field when calculating the linear wind field
    smooth_nsq = .False.,           ! set to true to provide additional vertical smoothing of Nsq within spatial_linear_winds

    ! NOTE THIS DOES NOT WORK RIGHT NOW
    ! rm_N_squared = 9e-5,          ! set this to use a fixed brunt-vaisalla frequency in linear wind calculations
    ! remove_lowres_linear = false, ! attempt to "remove" the low resolution linear winds from the forcing data
    ! rm_linear_contribution = 0.4, ! fraction of linear perturbation to remove from the low-res wind field (if rm_lowres_linear==true)

    ! Used to test possible model calibration... not sure what these will do longer term.
    ! To use these, you must also specify a filename and variable name to be read in for these fields. 
    ! nsq_calibration = false,
    ! linear_mask = false,

    ! Linear theory Look Up Table generation parameters
    ! NOTE: if you have memory problems running the model, decrease the n_X_values below
    ! direction ranges and number of bins
    dirmax = 6.283185307 ! 2*pi
    dirmin = 0
    n_dir_values = 36

    ! wind speed ranges and number of bins
    spdmax = 30
    spdmin = 0
    n_spd_values = 10

    ! BV frequency ranges (in log space) and number of bins
    nsqmax = -7.42  ! ln(6e-4) defaults to ln(max_stability)
    nsqmin = -16.12 ! ln(1e-7) defaults to ln(min_stability)
    n_nsq_values = 10

    ! NOTE: this look up table requires a LOT of RAM.  (and you still need some for the rest of the model)
    !   Calculate bytes of RAM required as nx * ny * nz * n_dir * n_spd * n_nsq * 2 * 4
    !   e.g. 320 * 250 * 14 * 36 * 10 * 10 * 2 * 4 = 30GB!
    !   (* 2 is for having to store U and V)
    !   (* 4 is for the number of bytes per float)

    ! To speed up model initialization, the look up table can be saved to disk (write_LUT=True)
    ! On the next run, the LUT can be read back in (read_LUT=True).
    ! Error checking will be performed and if the attributes above, or the spatial domain
    ! used to generate the LUT does not match the values for the current run it will regenerate the LUT.
    read_LUT  = True    ! read the Look up table from a specified file
    write_LUT = True    ! write the look up table to the specified file
    LUT_filename = "Linear_Theory_LUT.nc"
/


!---------------------------------------------------------
!   Optionally specified Restart information
!---------------------------------------------------------
&restart_info
    ! file to read for initial conditions (an ICAR output file will work)
    restart_file = "restart/icar_1990_09_30_01-01.nc",
    
    ! date to start from, used to calculate position in both restart file and forcing file
    restart_date =  1990, 9, 30, 23, 0, 0
/ -->
