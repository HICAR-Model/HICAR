!>------------------------------------------------------------
!!  Implementation of domain object
!!
!!  implements all domain type bound procedures
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
submodule(domain_interface) domain_implementation
    use assertions_mod,       only : assert, assertions
    use mod_atm_utilities,    only : exner_function, update_pressure, compute_ivt, compute_iq
    use icar_constants
    use string,               only : str
    use co_util,              only : broadcast
    use io_routines,          only : io_write, io_read
    use geo,                  only : geo_lut, geo_interp, geo_interp2d, standardize_coordinates
    use array_utilities,      only : array_offset_x, array_offset_y, smooth_array, smooth_array_2d, make_2d_x, make_2d_y
    use vertical_interpolation,only : vinterp, vLUT
    use wind_surf,            only : calc_Sx, calc_TPI
    use output_metadata,            only : get_varname
    use mod_wrf_constants,    only : gravity, R_d, KARMAN
    implicit none

    interface setup
        module procedure setup_var, setup_exch
    end interface
    
    real, allocatable :: mod_temp_3d(:,:,:), surf_temp_1(:,:), surf_temp_2(:,:)
    real, parameter::deg2rad=0.017453293 !2*pi/360
    
    ! primary public routines : init, get_initial_conditions, halo_send, halo_retrieve, or halo_exchange
contains


    !> -------------------------------
    !! Initialize the size of the domain
    !!
    !! -------------------------------
    module subroutine init(this, options)
        class(domain_t), intent(inout) :: this
        type(options_t), intent(inout) :: options

        this%dx = options%parameters%dx

        call this%var_request(options)
        
        call read_domain_shape(this, options)
        
        call create_variables(this, options)
        
        call initialize_core_variables(this, options)  ! split into several subroutines?

        call read_land_variables(this, options)
        
        call setup_meta_data(this, options)

        call set_var_lists(this, options)

        if (options%parameters%batched_exch) call setup_batch_exch(this)

        call init_relax_filters(this)

    end subroutine
    
    subroutine set_var_lists(this, options)
        class(domain_t), intent(inout) :: this
        type(options_t), intent(in)    :: options

        integer :: var_list(kMAX_STORAGE_VARS)

        !Advection variables -- these are exchanged AND advected
        if (options%vars_to_advect(kVARS%water_vapor)>0) call this%adv_vars%add_var('qv', this%water_vapor%meta_data)
        if (options%vars_to_advect(kVARS%potential_temperature)>0) call this%adv_vars%add_var('theta', this%potential_temperature%meta_data) 
        if (options%vars_to_advect(kVARS%cloud_water)>0) call this%adv_vars%add_var('qc', this%cloud_water_mass%meta_data)                  
        if (options%vars_to_advect(kVARS%rain_in_air)>0) call this%adv_vars%add_var('qr', this%rain_mass%meta_data)                    
        if (options%vars_to_advect(kVARS%snow_in_air)>0) call this%adv_vars%add_var('qs', this%snow_mass%meta_data)                    
        if (options%vars_to_advect(kVARS%cloud_ice)>0) call this%adv_vars%add_var('qi', this%cloud_ice_mass%meta_data)                      
        if (options%vars_to_advect(kVARS%graupel_in_air)>0) call this%adv_vars%add_var('qg', this%graupel_mass%meta_data)                 
        if (options%vars_to_advect(kVARS%ice_number_concentration)>0)  call this%adv_vars%add_var('ni', this%cloud_ice_number%meta_data)       
        if (options%vars_to_advect(kVARS%rain_number_concentration)>0) call this%adv_vars%add_var('nr', this%rain_number%meta_data)      
        if (options%vars_to_advect(kVARS%snow_number_concentration)>0) call this%adv_vars%add_var('ns', this%snow_number%meta_data)      
        if (options%vars_to_advect(kVARS%graupel_number_concentration)>0) call this%adv_vars%add_var('ng', this%graupel_number%meta_data)   
        if (options%vars_to_advect(kVARS%ice1_a)>0) call this%adv_vars%add_var('ice1_a', this%ice1_a%meta_data)   
        if (options%vars_to_advect(kVARS%ice1_c)>0) call this%adv_vars%add_var('ice1_c', this%ice1_c%meta_data)   
        if (options%vars_to_advect(kVARS%ice2_mass)>0) call this%adv_vars%add_var('ice2_mass', this%ice2_mass%meta_data)   
        if (options%vars_to_advect(kVARS%ice2_number)>0) call this%adv_vars%add_var('ice2_number', this%ice2_number%meta_data)   
        if (options%vars_to_advect(kVARS%ice2_a)>0) call this%adv_vars%add_var('ice2_a', this%ice2_a%meta_data)   
        if (options%vars_to_advect(kVARS%ice2_c)>0) call this%adv_vars%add_var('ice2_c', this%ice2_c%meta_data)   
        if (options%vars_to_advect(kVARS%ice3_mass)>0) call this%adv_vars%add_var('ice3_mass', this%ice3_mass%meta_data)   
        if (options%vars_to_advect(kVARS%ice3_number)>0) call this%adv_vars%add_var('ice3_number', this%ice3_number%meta_data)   
        if (options%vars_to_advect(kVARS%ice3_a)>0) call this%adv_vars%add_var('ice3_a', this%ice3_a%meta_data)   
        if (options%vars_to_advect(kVARS%ice3_c)>0) call this%adv_vars%add_var('ice3_c', this%ice3_c%meta_data)   

        !Exchange-only variables
        if (options%vars_to_exch(kVARS%sensible_heat)>0) call this%exch_vars%add_var('hfss', this%sensible_heat) 
        if (options%vars_to_exch(kVARS%skin_temperature)>0) call this%exch_vars%add_var('tsfe', this%skin_temperature)   
        if (options%vars_to_exch(kVARS%Ds)>0) call this%exch_vars%add_var('Ds', this%Ds)   
        if (options%vars_to_exch(kVARS%fsnow)>0) call this%exch_vars%add_var('fsnow', this%fsnow)   
        if (options%vars_to_exch(kVARS%Sice)>0) call this%exch_vars%add_var('Sice', this%Sice)   
        if (options%vars_to_exch(kVARS%Sliq)>0) call this%exch_vars%add_var('Sliq', this%Sliq)   
        if (options%vars_to_exch(kVARS%Tsnow)>0) call this%exch_vars%add_var('Tsnow', this%Tsnow)   
        if (options%vars_to_exch(kVARS%Nsnow)>0) call this%exch_vars%add_var('Nsnow', this%Nsnow)   

        var_list = options%io_options%vars_for_output + options%vars_for_restart
        if (0<var_list( kVARS%u) )                          call this%vars_to_out%add_var( trim( get_varname( kVARS%u                            )), this%u%meta_data)
        if (0<var_list( kVARS%v) )                          call this%vars_to_out%add_var( trim( get_varname( kVARS%v                            )), this%v%meta_data)
        if (0<var_list( kVARS%w) )                          call this%vars_to_out%add_var( trim( get_varname( kVARS%w                            )), this%w%meta_data)
        if (0<var_list( kVARS%w_real) )                     call this%vars_to_out%add_var( trim( get_varname( kVARS%w_real                       )), this%w_real)
        if (0<var_list( kVARS%nsquared) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%nsquared                     )), this%nsquared)
        if (0<var_list( kVARS%water_vapor) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%water_vapor                  )), this%water_vapor%meta_data)
        if (0<var_list( kVARS%potential_temperature) )      call this%vars_to_out%add_var( trim( get_varname( kVARS%potential_temperature        )), this%potential_temperature%meta_data)
        if (0<var_list( kVARS%cloud_water) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%cloud_water                  )), this%cloud_water_mass%meta_data)
        if (0<var_list( kVARS%cloud_number_concentration))  call this%vars_to_out%add_var( trim( get_varname( kVARS%cloud_number_concentration   )), this%cloud_number%meta_data)
        if (0<var_list( kVARS%cloud_ice) )                  call this%vars_to_out%add_var( trim( get_varname( kVARS%cloud_ice                    )), this%cloud_ice_mass%meta_data)
        if (0<var_list( kVARS%ice_number_concentration))    call this%vars_to_out%add_var( trim( get_varname( kVARS%ice_number_concentration     )), this%cloud_ice_number%meta_data)
        if (0<var_list( kVARS%rain_in_air) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%rain_in_air                  )), this%rain_mass%meta_data)
        if (0<var_list( kVARS%rain_number_concentration))   call this%vars_to_out%add_var( trim( get_varname( kVARS%rain_number_concentration    )), this%rain_number%meta_data)
        if (0<var_list( kVARS%snow_in_air) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_in_air                  )), this%snow_mass%meta_data)
        if (0<var_list( kVARS%snow_number_concentration) )  call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_number_concentration    )), this%snow_number%meta_data)
        if (0<var_list( kVARS%graupel_in_air) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%graupel_in_air               )), this%graupel_mass%meta_data)
        if (0<var_list( kVARS%graupel_number_concentration))call this%vars_to_out%add_var( trim( get_varname( kVARS%graupel_number_concentration )), this%graupel_number%meta_data)
        if (0<var_list( kVARS%ice1_a))                      call this%vars_to_out%add_var( trim( get_varname( kVARS%ice1_a                       )), this%ice1_a%meta_data)
        if (0<var_list( kVARS%ice1_c))                      call this%vars_to_out%add_var( trim( get_varname( kVARS%ice1_c                       )), this%ice1_c%meta_data)
        if (0<var_list( kVARS%ice2_mass))                   call this%vars_to_out%add_var( trim( get_varname( kVARS%ice2_mass                    )), this%ice2_mass%meta_data)
        if (0<var_list( kVARS%ice2_number))                 call this%vars_to_out%add_var( trim( get_varname( kVARS%ice2_number                  )), this%ice2_number%meta_data)
        if (0<var_list( kVARS%ice2_a))                      call this%vars_to_out%add_var( trim( get_varname( kVARS%ice2_a                       )), this%ice2_a%meta_data)
        if (0<var_list( kVARS%ice2_c))                      call this%vars_to_out%add_var( trim( get_varname( kVARS%ice2_c                       )), this%ice2_c%meta_data)
        if (0<var_list( kVARS%ice3_mass))                   call this%vars_to_out%add_var( trim( get_varname( kVARS%ice3_mass                    )), this%ice3_mass%meta_data)
        if (0<var_list( kVARS%ice3_number))                 call this%vars_to_out%add_var( trim( get_varname( kVARS%ice3_number                  )), this%ice3_number%meta_data)
        if (0<var_list( kVARS%ice3_a))                      call this%vars_to_out%add_var( trim( get_varname( kVARS%ice3_a                       )), this%ice3_a%meta_data)
        if (0<var_list( kVARS%ice3_c))                      call this%vars_to_out%add_var( trim( get_varname( kVARS%ice3_c                       )), this%ice3_c%meta_data)

        if (0<var_list( kVARS%precipitation) )              call this%vars_to_out%add_var( trim( get_varname( kVARS%precipitation                )), this%accumulated_precipitation)
        if (0<var_list( kVARS%convective_precipitation) )   call this%vars_to_out%add_var( trim( get_varname( kVARS%convective_precipitation     )), this%accumulated_convective_pcp)
        if (0<var_list( kVARS%snowfall) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%snowfall                     )), this%accumulated_snowfall)
        if (0<var_list( kVARS%graupel) )                    call this%vars_to_out%add_var( trim( get_varname( kVARS%graupel                      )), this%graupel)
        if (0<var_list( kVARS%pressure) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%pressure                     )), this%pressure)
        if (0<var_list( kVARS%temperature) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%temperature                  )), this%temperature)
        if (0<var_list( kVARS%exner) )                      call this%vars_to_out%add_var( trim( get_varname( kVARS%exner                        )), this%exner)
        if (0<var_list( kVARS%z) )                          call this%vars_to_out%add_var( trim( get_varname( kVARS%z                            )), this%z)
        if (0<var_list( kVARS%dz_interface) )               call this%vars_to_out%add_var( trim( get_varname( kVARS%dz_interface                 )), this%dz_interface)
        if (0<var_list( kVARS%z_interface) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%z_interface                  )), this%z_interface)
        if (0<var_list( kVARS%dz) )                         call this%vars_to_out%add_var( trim( get_varname( kVARS%dz                           )), this%dz_mass)
        if (0<var_list( kVARS%density) )                    call this%vars_to_out%add_var( trim( get_varname( kVARS%density                      )), this%density)
        if (0<var_list( kVARS%pressure_interface) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%pressure_interface           )), this%pressure_interface)
        if (0<var_list( kVARS%cloud_fraction) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%cloud_fraction               )), this%cloud_fraction)
        if (0<var_list( kVARS%shortwave) )                  call this%vars_to_out%add_var( trim( get_varname( kVARS%shortwave                    )), this%shortwave)
        if (0<var_list( kVARS%shortwave_direct) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%shortwave_direct             )), this%shortwave_direct)
        if (0<var_list( kVARS%shortwave_diffuse) )          call this%vars_to_out%add_var( trim( get_varname( kVARS%shortwave_diffuse            )), this%shortwave_diffuse)
        if (0<var_list( kVARS%shortwave_direct_above) )     call this%vars_to_out%add_var( trim( get_varname( kVARS%shortwave_direct_above       )), this%shortwave_direct_above) !! MJ aded
        if (0<var_list( kVARS%shortwave_total) )            call this%vars_to_out%add_var( trim( get_varname( kVARS%shortwave_total              )), this%shortwave_total) !! MJ aded
        if (0<var_list( kVARS%longwave) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%longwave                     )), this%longwave)
        if (0<var_list( kVARS%vegetation_fraction) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%vegetation_fraction          )), this%vegetation_fraction)
        if (0<var_list( kVARS%vegetation_fraction_max) )    call this%vars_to_out%add_var( trim( get_varname( kVARS%vegetation_fraction_max      )), this%vegetation_fraction_max)
        if (0<var_list( kVARS%vegetation_fraction_out) )    call this%vars_to_out%add_var( trim( get_varname( kVARS%vegetation_fraction_out      )), this%vegetation_fraction_out)
        if (0<var_list( kVARS%lai) )                        call this%vars_to_out%add_var( trim( get_varname( kVARS%lai                          )), this%lai)
        if (0<var_list( kVARS%sai) )                        call this%vars_to_out%add_var( trim( get_varname( kVARS%sai                          )), this%sai)
        if (0<var_list( kVARS%crop_type) )                  call this%vars_to_out%add_var( trim( get_varname( kVARS%crop_type                    )), this%crop_type)
        if (0<var_list( kVARS%date_planting) )              call this%vars_to_out%add_var( trim( get_varname( kVARS%date_planting                )), this%date_planting)
        if (0<var_list( kVARS%date_harvest) )               call this%vars_to_out%add_var( trim( get_varname( kVARS%date_harvest                 )), this%date_harvest)
        if (0<var_list( kVARS%growing_season_gdd) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%growing_season_gdd           )), this%growing_season_gdd)
        if (0<var_list( kVARS%irr_frac_total) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_frac_total               )), this%irr_frac_total)
        if (0<var_list( kVARS%irr_frac_sprinkler) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_frac_sprinkler           )), this%irr_frac_sprinkler)
        if (0<var_list( kVARS%irr_frac_micro) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_frac_micro               )), this%irr_frac_micro)
        if (0<var_list( kVARS%irr_frac_flood) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_frac_flood               )), this%irr_frac_flood)
        if (0<var_list( kVARS%irr_alloc_sprinkler) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_alloc_sprinkler          )), this%irr_alloc_sprinkler)
        if (0<var_list( kVARS%irr_alloc_micro) )            call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_alloc_micro              )), this%irr_alloc_micro)
        if (0<var_list( kVARS%irr_alloc_flood) )            call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_alloc_flood              )), this%irr_alloc_flood)
        if (0<var_list( kVARS%irr_evap_loss_sprinkler) )    call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_evap_loss_sprinkler      )), this%irr_evap_loss_sprinkler)
        if (0<var_list( kVARS%irr_amt_sprinkler) )          call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_amt_sprinkler            )), this%irr_amt_sprinkler)
        if (0<var_list( kVARS%irr_amt_micro) )              call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_amt_micro                )), this%irr_amt_micro)
        if (0<var_list( kVARS%irr_amt_flood) )              call this%vars_to_out%add_var( trim( get_varname( kVARS%irr_amt_flood                )), this%irr_amt_flood)
        if (0<var_list( kVARS%evap_heat_sprinkler) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%evap_heat_sprinkler          )), this%evap_heat_sprinkler)
        if (0<var_list( kVARS%mass_ag_grain) )              call this%vars_to_out%add_var( trim( get_varname( kVARS%mass_ag_grain                )), this%mass_ag_grain)
        if (0<var_list( kVARS%growing_degree_days) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%growing_degree_days          )), this%growing_degree_days)
        if (0<var_list( kVARS%net_ecosystem_exchange) )     call this%vars_to_out%add_var( trim( get_varname( kVARS%net_ecosystem_exchange       )), this%net_ecosystem_exchange)
        if (0<var_list( kVARS%gross_primary_prod) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%gross_primary_prod           )), this%gross_primary_prod)
        if (0<var_list( kVARS%net_primary_prod) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%net_primary_prod             )), this%net_primary_prod)
        if (0<var_list( kVARS%apar) )                       call this%vars_to_out%add_var( trim( get_varname( kVARS%apar                         )), this%apar)
        if (0<var_list( kVARS%photosynthesis_total) )       call this%vars_to_out%add_var( trim( get_varname( kVARS%photosynthesis_total         )), this%photosynthesis_total)
        if (0<var_list( kVARS%stomatal_resist_total) )      call this%vars_to_out%add_var( trim( get_varname( kVARS%stomatal_resist_total        )), this%stomatal_resist_total)
        if (0<var_list( kVARS%stomatal_resist_sun) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%stomatal_resist_sun          )), this%stomatal_resist_sun)
        if (0<var_list( kVARS%stomatal_resist_shade) )      call this%vars_to_out%add_var( trim( get_varname( kVARS%stomatal_resist_shade        )), this%stomatal_resist_shade)
        if (0<var_list( kVARS%gecros_state) )               call this%vars_to_out%add_var( trim( get_varname( kVARS%gecros_state                 )), this%gecros_state)
        if (0<var_list( kVARS%canopy_water) )               call this%vars_to_out%add_var( trim( get_varname( kVARS%canopy_water                 )), this%canopy_water)
        if (0<var_list( kVARS%canopy_water_ice) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%canopy_water_ice             )), this%canopy_water_ice)
        if (0<var_list( kVARS%canopy_water_liquid) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%canopy_water_liquid          )), this%canopy_water_liquid)
        if (0<var_list( kVARS%canopy_vapor_pressure) )      call this%vars_to_out%add_var( trim( get_varname( kVARS%canopy_vapor_pressure        )), this%canopy_vapor_pressure)
        if (0<var_list( kVARS%canopy_temperature) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%canopy_temperature           )), this%canopy_temperature)
        if (0<var_list( kVARS%canopy_fwet) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%canopy_fwet                  )), this%canopy_fwet)
        if (0<var_list( kVARS%veg_leaf_temperature) )       call this%vars_to_out%add_var( trim( get_varname( kVARS%veg_leaf_temperature         )), this%veg_leaf_temperature)
        if (0<var_list( kVARS%ground_surf_temperature) )    call this%vars_to_out%add_var( trim( get_varname( kVARS%ground_surf_temperature      )), this%ground_surf_temperature)
        if (0<var_list( kVARS%frac_within_gap) )            call this%vars_to_out%add_var( trim( get_varname( kVARS%frac_within_gap              )), this%frac_within_gap)
        if (0<var_list( kVARS%frac_between_gap) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%frac_between_gap             )), this%frac_between_gap)
        if (0<var_list( kVARS%ground_temperature_bare) )    call this%vars_to_out%add_var( trim( get_varname( kVARS%ground_temperature_bare      )), this%ground_temperature_bare)
        if (0<var_list( kVARS%ground_temperature_canopy) )  call this%vars_to_out%add_var( trim( get_varname( kVARS%ground_temperature_canopy    )), this%ground_temperature_canopy)
        if (0<var_list( kVARS%snowfall_ground) )            call this%vars_to_out%add_var( trim( get_varname( kVARS%snowfall_ground              )), this%snowfall_ground)
        if (0<var_list( kVARS%rainfall_ground) )            call this%vars_to_out%add_var( trim( get_varname( kVARS%rainfall_ground              )), this%rainfall_ground)
        if (0<var_list( kVARS%snow_water_equivalent) )      call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_water_equivalent        )), this%snow_water_equivalent)
        if (0<var_list( kVARS%snow_water_eq_prev) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_water_eq_prev           )), this%snow_water_eq_prev)
        if (0<var_list( kVARS%albedo) )                     call this%vars_to_out%add_var( trim( get_varname( kVARS%albedo                       )), this%albedo)
        if (0<var_list( kVARS%snow_albedo_prev) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_albedo_prev             )), this%snow_albedo_prev)
        if (0<var_list( kVARS%snow_temperature) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_temperature             )), this%snow_temperature)
        if (0<var_list( kVARS%snow_layer_depth) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_layer_depth             )), this%snow_layer_depth)
        if (0<var_list( kVARS%snow_layer_ice) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_layer_ice               )), this%snow_layer_ice)
        if (0<var_list( kVARS%snow_layer_liquid_water) )    call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_layer_liquid_water      )), this%snow_layer_liquid_water)
        if (0<var_list( kVARS%snow_age_factor) )            call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_age_factor              )), this%snow_age_factor)
        if (0<var_list( kVARS%snow_height) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%snow_height                  )), this%snow_height)
        if (0<var_list( kVARS%skin_temperature) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%skin_temperature             )), this%skin_temperature)
        if (0<var_list( kVARS%sst) )                        call this%vars_to_out%add_var( trim( get_varname( kVARS%sst                          )), this%sst)
        if (0<var_list( kVARS%soil_water_content) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_water_content           )), this%soil_water_content)
        if (0<var_list( kVARS%eq_soil_moisture) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%eq_soil_moisture             )), this%eq_soil_moisture)
        if (0<var_list( kVARS%smc_watertable_deep) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%smc_watertable_deep          )), this%smc_watertable_deep)
        if (0<var_list( kVARS%recharge) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%recharge                     )), this%recharge)
        if (0<var_list( kVARS%recharge_deep) )              call this%vars_to_out%add_var( trim( get_varname( kVARS%recharge_deep                )), this%recharge_deep)
        if (0<var_list( kVARS%soil_temperature) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_temperature             )), this%soil_temperature)
        if (0<var_list( kVARS%latitude) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%latitude                     )), this%latitude)
        if (0<var_list( kVARS%longitude) )                  call this%vars_to_out%add_var( trim( get_varname( kVARS%longitude                    )), this%longitude)
        if (0<var_list( kVARS%u_latitude) )                 call this%vars_to_out%add_var( trim( get_varname( kVARS%u_latitude                   )), this%u_latitude)
        if (0<var_list( kVARS%u_longitude) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%u_longitude                  )), this%u_longitude)
        if (0<var_list( kVARS%v_latitude) )                 call this%vars_to_out%add_var( trim( get_varname( kVARS%v_latitude                   )), this%v_latitude)
        if (0<var_list( kVARS%v_longitude) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%v_longitude                  )), this%v_longitude)
        if (0<var_list( kVARS%terrain) )                    call this%vars_to_out%add_var( trim( get_varname( kVARS%terrain                      )), this%terrain)
        if (0<var_list( kVARS%sensible_heat) )              call this%vars_to_out%add_var( trim( get_varname( kVARS%sensible_heat                )), this%sensible_heat)
        if (0<var_list( kVARS%latent_heat) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%latent_heat                  )), this%latent_heat)
        if (0<var_list( kVARS%u_10m) )                      call this%vars_to_out%add_var( trim( get_varname( kVARS%u_10m                        )), this%u_10m)
        if (0<var_list( kVARS%v_10m) )                      call this%vars_to_out%add_var( trim( get_varname( kVARS%v_10m                        )), this%v_10m)
        if (0<var_list( kVARS%windspd_10m) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%windspd_10m                  )), this%windspd_10m) !! MJ added
        if (0<var_list( kVARS%coeff_momentum_drag) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%coeff_momentum_drag          )), this%coeff_momentum_drag)
        if (0<var_list( kVARS%chs) )                        call this%vars_to_out%add_var( trim( get_varname( kVARS%chs                          )), this%chs)
        if (0<var_list( kVARS%chs2) )                       call this%vars_to_out%add_var( trim( get_varname( kVARS%chs2                         )), this%chs2)
        if (0<var_list( kVARS%cqs2) )                       call this%vars_to_out%add_var( trim( get_varname( kVARS%cqs2                         )), this%cqs2)
        if (0<var_list( kVARS%br) )                         call this%vars_to_out%add_var( trim( get_varname( kVARS%br                           )), this%br)
        if (0<var_list( kVARS%QFX) )                        call this%vars_to_out%add_var( trim( get_varname( kVARS%QFX                          )), this%qfx)
        if (0<var_list( kVARS%surface_rad_temperature) )    call this%vars_to_out%add_var( trim( get_varname( kVARS%surface_rad_temperature      )), this%surface_rad_temperature)
        if (0<var_list( kVARS%temperature_2m) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%temperature_2m               )), this%temperature_2m)
        if (0<var_list( kVARS%humidity_2m) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%humidity_2m                  )), this%humidity_2m)
        if (0<var_list( kVARS%temperature_2m_veg) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%temperature_2m_veg           )), this%temperature_2m_veg)
        if (0<var_list( kVARS%temperature_2m_bare) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%temperature_2m_bare          )), this%temperature_2m_bare)
        if (0<var_list( kVARS%mixing_ratio_2m_veg) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%mixing_ratio_2m_veg          )), this%mixing_ratio_2m_veg)
        if (0<var_list( kVARS%mixing_ratio_2m_bare) )       call this%vars_to_out%add_var( trim( get_varname( kVARS%mixing_ratio_2m_bare         )), this%mixing_ratio_2m_bare)
        if (0<var_list( kVARS%surface_pressure) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%surface_pressure             )), this%surface_pressure)
        if (0<var_list( kVARS%rad_absorbed_total) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%rad_absorbed_total           )), this%rad_absorbed_total)
        if (0<var_list( kVARS%rad_absorbed_veg) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%rad_absorbed_veg             )), this%rad_absorbed_veg)
        if (0<var_list( kVARS%rad_absorbed_bare) )          call this%vars_to_out%add_var( trim( get_varname( kVARS%rad_absorbed_bare            )), this%rad_absorbed_bare)
        if (0<var_list( kVARS%rad_net_longwave) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%rad_net_longwave             )), this%rad_net_longwave)
        if (0<var_list( kVARS%longwave_up) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%longwave_up                  )), this%longwave_up)
        if (0<var_list( kVARS%ground_heat_flux) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%ground_heat_flux             )), this%ground_heat_flux)
        if (0<var_list( kVARS%soil_deep_temperature) )      call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_deep_temperature        )), this%soil_deep_temperature)
        if (0<var_list( kVARS%evap_canopy) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%evap_canopy                  )), this%evap_canopy)
        if (0<var_list( kVARS%evap_soil_surface) )          call this%vars_to_out%add_var( trim( get_varname( kVARS%evap_soil_surface            )), this%evap_soil_surface)
        if (0<var_list( kVARS%transpiration_rate) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%transpiration_rate           )), this%transpiration_rate)
        if (0<var_list( kVARS%ch_veg) )                     call this%vars_to_out%add_var( trim( get_varname( kVARS%ch_veg                       )), this%ch_veg)
        if (0<var_list( kVARS%ch_veg_2m) )                  call this%vars_to_out%add_var( trim( get_varname( kVARS%ch_veg_2m                    )), this%ch_veg_2m)
        if (0<var_list( kVARS%ch_bare) )                    call this%vars_to_out%add_var( trim( get_varname( kVARS%ch_bare                      )), this%ch_bare)
        if (0<var_list( kVARS%ch_bare_2m) )                 call this%vars_to_out%add_var( trim( get_varname( kVARS%ch_bare_2m                   )), this%ch_bare_2m)
        if (0<var_list( kVARS%ch_under_canopy) )            call this%vars_to_out%add_var( trim( get_varname( kVARS%ch_under_canopy              )), this%ch_under_canopy)
        if (0<var_list( kVARS%ch_leaf) )                    call this%vars_to_out%add_var( trim( get_varname( kVARS%ch_leaf                      )), this%ch_leaf)
        if (0<var_list( kVARS%sensible_heat_veg) )          call this%vars_to_out%add_var( trim( get_varname( kVARS%sensible_heat_veg            )), this%sensible_heat_veg)
        if (0<var_list( kVARS%sensible_heat_bare) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%sensible_heat_bare           )), this%sensible_heat_bare)
        if (0<var_list( kVARS%sensible_heat_canopy) )       call this%vars_to_out%add_var( trim( get_varname( kVARS%sensible_heat_canopy         )), this%sensible_heat_canopy)
        if (0<var_list( kVARS%evap_heat_veg) )              call this%vars_to_out%add_var( trim( get_varname( kVARS%evap_heat_veg                )), this%evap_heat_veg)
        if (0<var_list( kVARS%evap_heat_bare) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%evap_heat_bare               )), this%evap_heat_bare)
        if (0<var_list( kVARS%evap_heat_canopy) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%evap_heat_canopy             )), this%evap_heat_canopy)
        if (0<var_list( kVARS%transpiration_heat) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%transpiration_heat           )), this%transpiration_heat)
        if (0<var_list( kVARS%ground_heat_veg) )            call this%vars_to_out%add_var( trim( get_varname( kVARS%ground_heat_veg              )), this%ground_heat_veg)
        if (0<var_list( kVARS%ground_heat_bare) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%ground_heat_bare             )), this%ground_heat_bare)
        if (0<var_list( kVARS%net_longwave_veg) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%net_longwave_veg             )), this%net_longwave_veg)
        if (0<var_list( kVARS%net_longwave_bare) )          call this%vars_to_out%add_var( trim( get_varname( kVARS%net_longwave_bare            )), this%net_longwave_bare)
        if (0<var_list( kVARS%net_longwave_canopy) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%net_longwave_canopy          )), this%net_longwave_canopy)
        if (0<var_list( kVARS%runoff_surface) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%runoff_surface               )), this%runoff_surface)
        if (0<var_list( kVARS%runoff_subsurface) )          call this%vars_to_out%add_var( trim( get_varname( kVARS%runoff_subsurface            )), this%runoff_subsurface)
        if (0<var_list( kVARS%soil_totalmoisture) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_totalmoisture           )), this%soil_totalmoisture)
        if (0<var_list( kVARS%water_table_depth) )          call this%vars_to_out%add_var( trim( get_varname( kVARS%water_table_depth            )), this%water_table_depth)
        if (0<var_list( kVARS%water_aquifer) )              call this%vars_to_out%add_var( trim( get_varname( kVARS%water_aquifer                )), this%water_aquifer)
        if (0<var_list( kVARS%storage_gw) )                 call this%vars_to_out%add_var( trim( get_varname( kVARS%storage_gw                   )), this%storage_gw)
        if (0<var_list( kVARS%storage_lake) )               call this%vars_to_out%add_var( trim( get_varname( kVARS%storage_lake                 )), this%storage_lake)
        if (0<var_list( kVARS%roughness_z0) )               call this%vars_to_out%add_var( trim( get_varname( kVARS%roughness_z0                 )), this%roughness_z0)
        !if (0<var_list( kVARS%veg_type) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%veg_type                     )), this%veg_type)
        if (0<var_list( kVARS%mass_leaf) )                  call this%vars_to_out%add_var( trim( get_varname( kVARS%mass_leaf                    )), this%mass_leaf)
        if (0<var_list( kVARS%mass_root) )                  call this%vars_to_out%add_var( trim( get_varname( kVARS%mass_root                    )), this%mass_root)
        if (0<var_list( kVARS%mass_stem) )                  call this%vars_to_out%add_var( trim( get_varname( kVARS%mass_stem                    )), this%mass_stem)
        if (0<var_list( kVARS%mass_wood) )                  call this%vars_to_out%add_var( trim( get_varname( kVARS%mass_wood                    )), this%mass_wood)
        if (0<var_list( kVARS%soil_carbon_fast) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_carbon_fast             )), this%soil_carbon_fast)
        if (0<var_list( kVARS%soil_carbon_stable) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_carbon_stable           )), this%soil_carbon_stable)
        if (0<var_list( kVARS%soil_texture_1) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_texture_1               )), this%soil_texture_1)
        if (0<var_list( kVARS%soil_texture_2) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_texture_2               )), this%soil_texture_2)
        if (0<var_list( kVARS%soil_texture_3) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_texture_3               )), this%soil_texture_3)
        if (0<var_list( kVARS%soil_texture_4) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_texture_4               )), this%soil_texture_4)
        if (0<var_list( kVARS%soil_sand_and_clay) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%soil_sand_and_clay           )), this%soil_sand_and_clay)
        if (0<var_list( kVARS%re_cloud) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%re_cloud                     )), this%re_cloud)
        if (0<var_list( kVARS%re_ice) )                     call this%vars_to_out%add_var( trim( get_varname( kVARS%re_ice                       )), this%re_ice)
        if (0<var_list( kVARS%re_snow) )                    call this%vars_to_out%add_var( trim( get_varname( kVARS%re_snow                      )), this%re_snow)
        if (0<var_list( kVARS%ice1_rho) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%ice1_rho                      )), this%ice1_rho)
        if (0<var_list( kVARS%ice1_phi) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%ice1_phi                      )), this%ice1_phi)
        if (0<var_list( kVARS%ice2_rho) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%ice2_rho                      )), this%ice2_rho)
        if (0<var_list( kVARS%ice2_phi) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%ice2_phi                      )), this%ice2_phi)
        if (0<var_list( kVARS%ice3_rho) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%ice3_rho                      )), this%ice3_rho)
        if (0<var_list( kVARS%ice3_phi) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%ice3_phi                      )), this%ice3_phi)
        if (0<var_list( kVARS%out_longwave_rad) )           call this%vars_to_out%add_var( trim( get_varname( kVARS%out_longwave_rad             )), this%out_longwave_rad)
        if (0<var_list( kVARS%longwave_cloud_forcing) )     call this%vars_to_out%add_var( trim( get_varname( kVARS%longwave_cloud_forcing       )), this%longwave_cloud_forcing)
        if (0<var_list( kVARS%shortwave_cloud_forcing) )    call this%vars_to_out%add_var( trim( get_varname( kVARS%shortwave_cloud_forcing      )), this%shortwave_cloud_forcing)
        if (0<var_list( kVARS%cosine_zenith_angle) )        call this%vars_to_out%add_var( trim( get_varname( kVARS%cosine_zenith_angle          )), this%cosine_zenith_angle)
        if (0<var_list( kVARS%land_emissivity) )            call this%vars_to_out%add_var( trim( get_varname( kVARS%land_emissivity              )), this%land_emissivity)
        if (0<var_list( kVARS%temperature_interface) )      call this%vars_to_out%add_var( trim( get_varname( kVARS%temperature_interface        )), this%temperature_interface)
        if (0<var_list( kVARS%tend_swrad) )                 call this%vars_to_out%add_var( trim( get_varname( kVARS%tend_swrad                   )), this%tend_swrad)
        !! MJ added for FSM:
        if (0<var_list( kVARS%runoff_tstep) )               call this%vars_to_out%add_var( trim( get_varname( kVARS%runoff_tstep                 )), this%runoff_tstep)
        if (0<var_list( kVARS%Tsnow) )                      call this%vars_to_out%add_var( trim( get_varname( kVARS%Tsnow                        )), this%Tsnow)
        if (0<var_list( kVARS%Sice) )                       call this%vars_to_out%add_var( trim( get_varname( kVARS%Sice                         )), this%Sice)
        if (0<var_list( kVARS%Sliq) )                       call this%vars_to_out%add_var( trim( get_varname( kVARS%Sliq                         )), this%Sliq)
        if (0<var_list( kVARS%Ds) )                         call this%vars_to_out%add_var( trim( get_varname( kVARS%Ds                           )), this%Ds)
        if (0<var_list( kVARS%fsnow) )                      call this%vars_to_out%add_var( trim( get_varname( kVARS%fsnow                        )), this%fsnow)
        if (0<var_list( kVARS%Nsnow) )                      call this%vars_to_out%add_var( trim( get_varname( kVARS%Nsnow                        )), this%Nsnow)   
        if (0<var_list( kVARS%dm_salt))                     call this%vars_to_out%add_var( trim( get_varname( kVARS%dm_salt                      )), this%dm_salt) 
        if (0<var_list( kVARS%dm_susp))                     call this%vars_to_out%add_var( trim( get_varname( kVARS%dm_susp                      )), this%dm_susp) 
        if (0<var_list( kVARS%dm_subl))                     call this%vars_to_out%add_var( trim( get_varname( kVARS%dm_subl                      )), this%dm_subl) 
        if (0<var_list( kVARS%dm_slide))                    call this%vars_to_out%add_var( trim( get_varname( kVARS%dm_slide                     )), this%dm_slide) 

        !!
        if (0<var_list( kVARS%rainfall_tstep) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%rainfall_tstep               )), this%rainfall_tstep)        
        if (0<var_list( kVARS%snowfall_tstep) )             call this%vars_to_out%add_var( trim( get_varname( kVARS%snowfall_tstep               )), this%snowfall_tstep)        
        if (0<var_list( kVARS%meltflux_out_tstep) )         call this%vars_to_out%add_var( trim( get_varname( kVARS%meltflux_out_tstep           )), this%meltflux_out_tstep)        
        if (0<var_list( kVARS%Sliq_out) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%Sliq_out                     )), this%Sliq_out)        
        if (0<var_list( kVARS%slope) )                      call this%vars_to_out%add_var( trim( get_varname( kVARS%slope                        )), this%slope)        
        if (0<var_list( kVARS%slope_angle) )                call this%vars_to_out%add_var( trim( get_varname( kVARS%slope_angle                  )), this%slope_angle)        
        if (0<var_list( kVARS%aspect_angle) )               call this%vars_to_out%add_var( trim( get_varname( kVARS%aspect_angle                 )), this%aspect_angle)        
        if (0<var_list( kVARS%svf) )                        call this%vars_to_out%add_var( trim( get_varname( kVARS%svf                          )), this%svf)        
        if (0<var_list( kVARS%factor_p) )                   call this%vars_to_out%add_var( trim( get_varname( kVARS%factor_p                     )), this%factor_p)        
        if (0<var_list( kVARS%hlm) )                        call this%vars_to_out%add_var( trim( get_varname( kVARS%hlm                          )), this%hlm)        
        if (0<var_list( kVARS%hpbl) )                       call this%vars_to_out%add_var( trim( get_varname( kVARS%hpbl                         )), this%hpbl)
        if (0<var_list( kVARS%coeff_heat_exchange_3d) )     call this%vars_to_out%add_var( trim( get_varname( kVARS%coeff_heat_exchange_3d       )), this%coeff_heat_exchange_3d)
        if (0<var_list( kVARS%coeff_momentum_exchange_3d) ) call this%vars_to_out%add_var( trim( get_varname( kVARS%coeff_momentum_exchange_3d   )), this%coeff_momentum_exchange_3d)
        if (0<var_list( kVARS%wind_alpha) )                 call this%vars_to_out%add_var( trim( get_varname( kVARS%wind_alpha                   )), this%alpha)
        if (0<var_list( kVARS%froude) )                     call this%vars_to_out%add_var( trim( get_varname( kVARS%froude                       )), this%froude)
        if (0<var_list( kVARS%blk_ri) )                     call this%vars_to_out%add_var( trim( get_varname( kVARS%blk_ri                       )), this%Ri)

    end subroutine set_var_lists

    !> -------------------------------
    !! Initialize the arrays and co-arrays needed to perform a batch exchange
    !!
    !! -------------------------------
    module subroutine setup_batch_exch(this)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(variable_t) :: var

        integer :: n_2d, n_3d = 0
      
        ! Loop over all adv_vars and count how many are 3D
        call this%adv_vars%reset_iterator()
        
        do while (this%adv_vars%has_more_elements())
            var = this%adv_vars%next()
            if (var%three_d) n_3d = n_3d + 1
        end do
        
        ! Loop over all exch vars and count how many are 3D
        call this%exch_vars%reset_iterator()
        
        do while (this%exch_vars%has_more_elements())
            var = this%exch_vars%next()
            if (var%three_d) n_3d = n_3d + 1
        end do

        ! Determine number of 2D and 3D vars present
        n_2d = (this%adv_vars%n_vars+this%exch_vars%n_vars)-n_3d

        allocate(this%north_in_3d(n_3d,1:(this%grid%ns_halo_nx+this%grid%halo_size*2),&
                        this%kms:this%kme,1:this%grid%halo_size)[*])
        allocate(this%south_in_3d(n_3d,1:(this%grid%ns_halo_nx+this%grid%halo_size*2),&
                        this%kms:this%kme,1:this%grid%halo_size)[*])
        allocate(this%east_in_3d(n_3d,1:this%grid%halo_size,&
                        this%kms:this%kme,1:(this%grid%ew_halo_ny+this%grid%halo_size*2))[*])
        allocate(this%west_in_3d(n_3d,1:this%grid%halo_size,&
                        this%kms:this%kme,1:(this%grid%ew_halo_ny+this%grid%halo_size*2))[*])

        if (.not.(this%north_boundary)) allocate(this%north_buffer_3d(n_3d,1:(this%grid%ns_halo_nx+this%grid%halo_size*2),&
                        this%kms:this%kme,1:this%grid%halo_size))
        if (.not.(this%south_boundary)) allocate(this%south_buffer_3d(n_3d,1:(this%grid%ns_halo_nx+this%grid%halo_size*2),&
                        this%kms:this%kme,1:this%grid%halo_size))
        if (.not.(this%east_boundary)) allocate(this%east_buffer_3d(n_3d,1:this%grid%halo_size,&
                        this%kms:this%kme,1:(this%grid%ew_halo_ny+this%grid%halo_size*2)))
        if (.not.(this%west_boundary)) allocate(this%west_buffer_3d(n_3d,1:this%grid%halo_size,&
                        this%kms:this%kme,1:(this%grid%ew_halo_ny+this%grid%halo_size*2)))
      

        ! If no 2D vars present, don't allocate arrays (nothing should be calling exch 2D then)
        if (n_2d > 0) then
            allocate(this%north_in_2d(n_2d,1:(this%grid%ns_halo_nx+this%grid%halo_size*2),1:this%grid%halo_size)[*])
            allocate(this%south_in_2d(n_2d,1:(this%grid%ns_halo_nx+this%grid%halo_size*2),1:this%grid%halo_size)[*])
            allocate(this%east_in_2d(n_2d,1:this%grid%halo_size,1:(this%grid%ew_halo_ny+this%grid%halo_size*2))[*])
            allocate(this%west_in_2d(n_2d,1:this%grid%halo_size,1:(this%grid%ew_halo_ny+this%grid%halo_size*2))[*])

            if (.not.(this%north_boundary)) allocate(this%north_buffer_2d(n_2d,1:(this%grid%ns_halo_nx+this%grid%halo_size*2),1:this%grid%halo_size))
            if (.not.(this%south_boundary)) allocate(this%south_buffer_2d(n_2d,1:(this%grid%ns_halo_nx+this%grid%halo_size*2),1:this%grid%halo_size))
            if (.not.(this%east_boundary)) allocate(this%east_buffer_2d(n_2d,1:this%grid%halo_size,1:(this%grid%ew_halo_ny+this%grid%halo_size*2)))
            if (.not.(this%west_boundary)) allocate(this%west_buffer_2d(n_2d,1:this%grid%halo_size,1:(this%grid%ew_halo_ny+this%grid%halo_size*2)))
        endif
    
    end subroutine setup_batch_exch

    !> -------------------------------
    !! Set up the initial conditions for the domain
    !!
    !! This includes setting up all of the geographic interpolation now that we have the forcing grid
    !! and interpolating the first time step of forcing data on to the high res domain grids
    !!
    !! -------------------------------
    module subroutine get_initial_conditions(this, forcing, options, external_conditions)
      implicit none
      class(domain_t),  intent(inout) :: this
      type(boundary_t), intent(inout) :: forcing
      type(boundary_t), intent(inout), optional :: external_conditions
      type(options_t),  intent(in)    :: options

      integer :: i

      ! create geographic lookup table for domain
      call setup_geo_interpolation(this, forcing, options)

      ! for all variables with a forcing_var /= "", get forcing, interpolate to local domain
      call this%interpolate_forcing(forcing)
      
      if (allocated(this%znw).or.allocated(this%znu)) call init_znu(this)


      ! - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! read in any external variables, such as SWE or snow height.
      ! if (present(external_conditions).AND. (options%parameters%restart .eqv. .False.))    then
      if (present(external_conditions).AND. (options%parameters%restart .neqv. .True.))    then
        ! if (this_image()==1) write(*,*) "   (Dom) - Setting up ext files.  "
        ! create geographic lookup table for domain
        call setup_geo_interpolation(this, external_conditions, options)

        ! interpolate external variables to the hi-res grid
        call this%interpolate_external( external_conditions, options)
        ! if (this_image()==1) write(*,*) " interpolating exteral conditions"
      endif

      ! - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call diagnostic_update(this,options)
      call this%enforce_limits()
      this%model_time = options%parameters%start_time

    end subroutine


    !>------------------------------------------------------------
    !! Update model diagnostic fields
    !!
    !! Calculates most model diagnostic fields such as Psfc, 10m height winds and ustar
    !!
    !! @param domain    Model domain data structure to be updated
    !! @param options   Model options (not used at present)
    !!
    !!------------------------------------------------------------
    module subroutine diagnostic_update(this, options)
        implicit none
        class(domain_t),  intent(inout)   :: this
        type(options_t), intent(in)      :: options
        integer :: i, j, k
        real :: qsum
        
        associate(ims => this%ims, ime => this%ime,                             &
                  jms => this%jms, jme => this%jme,                             &
                  kms => this%kms, kme => this%kme,                             &
                  its => this%its, ite => this%ite,                             &
                  jts => this%jts, jte => this%jte,                             &
                  exner                 => this%exner%data_3d,                  &
                  pressure              => this%pressure%data_3d,               &
                  pressure_i            => this%pressure_interface%data_3d,     &
                  dz_i                  => this%dz_interface%data_3d,           &
                  dz_mass               => this%dz_mass%data_3d,                &
                  psfc                  => this%surface_pressure%data_2d,       &
                  density               => this%density%data_3d,                &
                  temperature           => this%temperature%data_3d,            &
                  qv                    => this%water_vapor%data_3d,            &  
                  cloud_water           => this%cloud_water_mass%data_3d,       &
                  rain_water            => this%rain_mass%data_3d,              &
                  cloud_ice             => this%cloud_ice_mass%data_3d,         &
                  snow_ice              => this%snow_mass%data_3d,              &
                  graupel_ice           => this%graupel_mass%data_3d,           &
                  temperature_i         => this%temperature_interface%data_3d,  &
                  u                     => this%u%data_3d,                      &
                  v                     => this%v%data_3d,                      &
                  u_mass                => this%u_mass%data_3d,                 &
                  v_mass                => this%v_mass%data_3d,                 &
                  potential_temperature => this%potential_temperature%data_3d )

        exner = exner_function(pressure)

        !Calculation of density
        if (associated(this%density%data_3d)) then
            do j = jms,jme
                do k = kms,kme
                    do i = ims,ime
                        qsum = qv(i,k,j)
                        if(associated(this%cloud_water_mass%data_3d)) qsum = qsum + this%cloud_water_mass%data_3d(i,k,j)
                        if(associated(this%cloud_ice_mass%data_3d)) qsum = qsum + this%cloud_ice_mass%data_3d(i,k,j)
                        if(associated(this%ice2_mass%data_3d)) qsum = qsum + this%ice2_mass%data_3d(i,k,j)
                        if(associated(this%ice3_mass%data_3d)) qsum = qsum + this%ice3_mass%data_3d(i,k,j)
                        if(associated(this%rain_mass%data_3d)) qsum = qsum + this%rain_mass%data_3d(i,k,j)
                        if(associated(this%snow_mass%data_3d)) qsum = qsum + this%snow_mass%data_3d(i,k,j)
                        if(associated(this%graupel_mass%data_3d)) qsum = qsum + this%graupel_mass%data_3d(i,k,j)
                        
                        temperature(i,k,j) = potential_temperature(i,k,j) * exner(i,k,j)
                        density(i,k,j) =  pressure(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(i,k,j))) ! kg/m^3
                    enddo
                enddo
            enddo
        endif
        
        temperature_i(ims:ime,kms,jms:jme) = temperature(ims:ime,kms,jms:jme) + (temperature(ims:ime,kms,jms:jme) - temperature(ims:ime,kms+1,jms:jme)) * 0.5
        pressure_i(ims:ime,kms,jms:jme) = pressure(ims:ime,kms,jms:jme) + (pressure(ims:ime,kms,jms:jme) - pressure(ims:ime,kms+1,jms:jme)) * 0.5

        do j = jms,jme
            do k = kms+1,kme
                do i = ims,ime
                    pressure_i(i,k,j) = (dz_i(i,k,j)*pressure(i,k-1,j)+dz_i(i,k-1,j)*pressure(i,k,j))/((dz_i(i,k-1,j)+dz_i(i,k,j)))
                    temperature_i(i,k,j) = (dz_i(i,k,j)*temperature(i,k-1,j)+dz_i(i,k-1,j)*temperature(i,k,j))/((dz_i(i,k-1,j)+dz_i(i,k,j)))
                enddo
            enddo
        enddo
        temperature_i(ims:ime,kme+1,jms:jme) = temperature(ims:ime,kme,jms:jme) + (temperature(ims:ime,kme,jms:jme) - temperature(ims:ime,kme-1,jms:jme)) * 0.5
        pressure_i(ims:ime,kme+1,jms:jme) = pressure(ims:ime,kme,jms:jme) + (pressure(ims:ime,kme,jms:jme) - pressure(ims:ime,kme-1,jms:jme)) * 0.5

        if (associated(this%u_mass%data_3d)) then
            do j = jms,jme
                do k = kms,kme
                    do i = ims,ime
                        u_mass(i,k,j) = (u(i+1,k,j) + u(i,k,j)) * 0.5
                        v_mass(i,k,j) = (v(i,k,j+1) + v(i,k,j)) * 0.5
                    enddo
                enddo
            enddo
        endif
                
        
        if (associated(this%surface_pressure%data_2d)) then
            psfc = pressure_i(ims:ime, kms, jms:jme)
        endif
        if (associated(this%ivt%data_2d)) then
            call compute_ivt(this%ivt%data_2d, qv, u_mass, v_mass, pressure_i(:,kms:kme,:))
        endif
        if (associated(this%iwv%data_2d)) then
            call compute_iq(this%iwv%data_2d, qv, pressure_i(:,kms:kme,:))
        endif
        if (associated(this%iwl%data_2d)) then
            mod_temp_3d = 0
            if (associated(this%cloud_water_mass%data_3d)) mod_temp_3d = mod_temp_3d + cloud_water
            if (associated(this%rain_mass%data_3d)) mod_temp_3d = mod_temp_3d + rain_water
            call compute_iq(this%iwl%data_2d, mod_temp_3d, pressure_i(:,kms:kme,:))
        endif
        if (associated(this%iwi%data_2d)) then
            mod_temp_3d = 0
            if (associated(this%cloud_ice_mass%data_3d)) mod_temp_3d = mod_temp_3d + cloud_ice
            if (associated(this%snow_mass%data_3d)) mod_temp_3d = mod_temp_3d + snow_ice
            if (associated(this%graupel_mass%data_3d)) mod_temp_3d = mod_temp_3d + graupel_ice
            call compute_iq(this%iwi%data_2d, mod_temp_3d, pressure_i(:,kms:kme,:))
        endif
        
        if (options%physics%surfacelayer == 0) then
            ! temporary constant
            if (associated(this%roughness_z0%data_2d)) then
                ! use log-law of the wall to convert from first model level to surface
                surf_temp_1 = karman / log((this%z%data_3d(ims:ime,kms,jms:jme) - this%terrain%data_2d(ims:ime,jms:jme)) / this%roughness_z0%data_2d(ims:ime,jms:jme))
                ! use log-law of the wall to convert from surface to 10m height
                surf_temp_2 = log(10.0 / this%roughness_z0%data_2d(ims:ime,jms:jme)) / karman
            endif

            if (associated(this%u_10m%data_2d)) then
                this%ustar        (ims:ime,jms:jme) = u_mass      (ims:ime,kms,jms:jme) * surf_temp_1
                this%u_10m%data_2d(ims:ime,jms:jme) = this%ustar(ims:ime,jms:jme)     * surf_temp_2
                this%ustar        (ims:ime,jms:jme) = v_mass      (ims:ime,kms,jms:jme) * surf_temp_1
                this%v_10m%data_2d(ims:ime,jms:jme) = this%ustar(ims:ime,jms:jme)     * surf_temp_2
            endif

            if (allocated(this%ustar)) then
                ! now calculate master ustar based on U and V combined in quadrature
                this%ustar(its:ite,jts:jte) = sqrt(u_mass(its:ite,kms,jts:jte)**2 + v_mass(its:ite,kms,jts:jte)**2) * surf_temp_1(its:ite,jts:jte)
            endif
        endif
        
        end associate

    end subroutine diagnostic_update



    !> -------------------------------
    !! Send the halos from all exchangable objects to their neighbors
    !!
    !! -------------------------------
    module subroutine halo_send(this)
      class(domain_t), intent(inout) :: this
      if (associated(this%water_vapor%data_3d))           call this%water_vapor%send()
      if (associated(this%potential_temperature%data_3d)) call this%potential_temperature%send()
      if (associated(this%cloud_water_mass%data_3d))      call this%cloud_water_mass%send()
      if (associated(this%cloud_number%data_3d))          call this%cloud_number%send()
      if (associated(this%cloud_ice_mass%data_3d))        call this%cloud_ice_mass%send()
      if (associated(this%cloud_ice_number%data_3d))      call this%cloud_ice_number%send()
      if (associated(this%rain_mass%data_3d))             call this%rain_mass%send()
      if (associated(this%rain_number%data_3d))           call this%rain_number%send()
      if (associated(this%snow_mass%data_3d))             call this%snow_mass%send()
      if (associated(this%snow_number%data_3d))           call this%snow_number%send()
      if (associated(this%graupel_mass%data_3d))          call this%graupel_mass%send()
      if (associated(this%graupel_number%data_3d))        call this%graupel_number%send()
      if (associated(this%ice1_a%data_3d))                call this%ice1_a%send()
      if (associated(this%ice1_c%data_3d))                call this%ice1_c%send()
      if (associated(this%ice2_mass%data_3d))             call this%ice2_mass%send()
      if (associated(this%ice2_number%data_3d))           call this%ice2_number%send()
      if (associated(this%ice2_a%data_3d))                call this%ice2_a%send()
      if (associated(this%ice2_c%data_3d))                call this%ice2_c%send()
      if (associated(this%ice3_mass%data_3d))             call this%ice3_mass%send()
      if (associated(this%ice3_number%data_3d))           call this%ice3_number%send()
      if (associated(this%ice3_a%data_3d))                call this%ice3_a%send()
      if (associated(this%ice3_c%data_3d))                call this%ice3_c%send()

    end subroutine

    !> -------------------------------
    !! Get the halos from all exchangable objects from their neighbors
    !!
    !! -------------------------------
    module subroutine halo_retrieve(this, wait_timer)
      class(domain_t), intent(inout) :: this
      type(timer_t),   intent(inout) :: wait_timer
      
      call wait_timer%start()
      if (associated(this%potential_temperature%data_3d)) call this%potential_temperature%retrieve()! the first retrieve call will sync all
      call wait_timer%stop()
      
      if (associated(this%water_vapor%data_3d))           call this%water_vapor%retrieve(no_sync=.True.)
      if (associated(this%cloud_water_mass%data_3d))      call this%cloud_water_mass%retrieve(no_sync=.True.)
      if (associated(this%cloud_number%data_3d))          call this%cloud_number%retrieve(no_sync=.True.)
      if (associated(this%cloud_ice_mass%data_3d))        call this%cloud_ice_mass%retrieve(no_sync=.True.)
      if (associated(this%cloud_ice_number%data_3d))      call this%cloud_ice_number%retrieve(no_sync=.True.)
      if (associated(this%rain_mass%data_3d))             call this%rain_mass%retrieve(no_sync=.True.)
      if (associated(this%rain_number%data_3d))           call this%rain_number%retrieve(no_sync=.True.)
      if (associated(this%snow_mass%data_3d))             call this%snow_mass%retrieve(no_sync=.True.)
      if (associated(this%snow_number%data_3d))           call this%snow_number%retrieve(no_sync=.True.)
      if (associated(this%graupel_mass%data_3d))          call this%graupel_mass%retrieve(no_sync=.True.)
      if (associated(this%graupel_number%data_3d))        call this%graupel_number%retrieve(no_sync=.True.)
      if (associated(this%ice1_a%data_3d))                call this%ice1_a%retrieve(no_sync=.True.)
      if (associated(this%ice1_c%data_3d))                call this%ice1_c%retrieve(no_sync=.True.)
      if (associated(this%ice2_mass%data_3d))             call this%ice2_mass%retrieve(no_sync=.True.)
      if (associated(this%ice2_number%data_3d))           call this%ice2_number%retrieve(no_sync=.True.)
      if (associated(this%ice2_a%data_3d))                call this%ice2_a%retrieve(no_sync=.True.)
      if (associated(this%ice2_c%data_3d))                call this%ice2_c%retrieve(no_sync=.True.)
      if (associated(this%ice3_mass%data_3d))             call this%ice3_mass%retrieve(no_sync=.True.)
      if (associated(this%ice3_number%data_3d))           call this%ice3_number%retrieve(no_sync=.True.)
      if (associated(this%ice3_a%data_3d))                call this%ice3_a%retrieve(no_sync=.True.)
      if (associated(this%ice3_c%data_3d))                call this%ice3_c%retrieve(no_sync=.True.)

    end subroutine

    !> -------------------------------
    !! Send and get the halos from all exchangable objects to/from their neighbors
    !!
    !! -------------------------------
    module subroutine halo_exchange(this, send_timer, ret_timer, wait_timer)
      class(domain_t), intent(inout) :: this
      type(timer_t),   intent(inout) :: send_timer, ret_timer, wait_timer
      
      call send_timer%start()
      call this%halo_send()
      call send_timer%stop()

      call ret_timer%start()
      call this%halo_retrieve(wait_timer)
      call ret_timer%stop()

    end subroutine
    
    module subroutine halo_3d_send_batch(this, exch_var_only)
        class(domain_t), intent(inout) :: this
        logical, intent(in) :: exch_var_only
        
        type(variable_t) :: var
        integer :: n, k_max

        call this%adv_vars%reset_iterator()
        call this%exch_vars%reset_iterator()
        n = 1
        ! Now iterate through the dictionary as long as there are more elements present
        if (.not.(exch_var_only)) then
            do while (this%adv_vars%has_more_elements())
                ! get the next variable
                var = this%adv_vars%next()
                if (var%three_d) then
                    if (.not.(this%north_boundary)) this%north_buffer_3d(n,1:(this%ite-this%its+3),:,:) = &
                            var%data_3d(this%its-1:this%ite+1,:,(this%jte-this%grid%halo_size+1):this%jte)
                    if (.not.(this%south_boundary)) this%south_buffer_3d(n,1:(this%ite-this%its+3),:,:) = &
                            var%data_3d(this%its-1:this%ite+1,:,this%jts:(this%jts+this%grid%halo_size-1))
                    if (.not.(this%east_boundary)) this%east_buffer_3d(n,:,:,1:(this%jte-this%jts+3)) = &
                            var%data_3d((this%ite-this%grid%halo_size+1):this%ite,:,this%jts-1:this%jte+1)
                    if (.not.(this%west_boundary)) this%west_buffer_3d(n,:,:,1:(this%jte-this%jts+3)) = &
                            var%data_3d(this%its:(this%its+this%grid%halo_size)-1,:,this%jts-1:this%jte+1)

                    n = n+1
                endif
            enddo
        endif

        ! Now iterate through the exchange-only objects as long as there are more elements present
        do while (this%exch_vars%has_more_elements())
            ! get the next variable
            var = this%exch_vars%next()
            if (var%three_d) then
                k_max = ubound(var%data_3d,2)
                if (.not.(this%north_boundary)) this%north_buffer_3d(n,1:(this%ite-this%its+3),1:k_max,:) = &
                        var%data_3d(this%its-1:this%ite+1,1:k_max,(this%jte-this%grid%halo_size+1):this%jte)
                if (.not.(this%south_boundary)) this%south_buffer_3d(n,1:(this%ite-this%its+3),1:k_max,:) = &
                        var%data_3d(this%its-1:this%ite+1,1:k_max,this%jts:(this%jts+this%grid%halo_size-1))
                if (.not.(this%east_boundary)) this%east_buffer_3d(n,:,1:k_max,1:(this%jte-this%jts+3)) = &
                        var%data_3d((this%ite-this%grid%halo_size+1):this%ite,1:k_max,this%jts-1:this%jte+1)
                if (.not.(this%west_boundary)) this%west_buffer_3d(n,:,1:k_max,1:(this%jte-this%jts+3)) = &
                        var%data_3d(this%its:(this%its+this%grid%halo_size)-1,1:k_max,this%jts-1:this%jte+1)

                n = n+1
            endif
        enddo

        if (.not.(this%north_boundary)) then
            !DIR$ PGAS DEFER_SYNC
            this%south_in_3d(:,:,:,:)[this%north_neighbor] = this%north_buffer_3d(:,:,:,:)
        endif
        if (.not.(this%south_boundary)) then
            !DIR$ PGAS DEFER_SYNC
            this%north_in_3d(:,:,:,:)[this%south_neighbor] = this%south_buffer_3d(:,:,:,:)
        endif
        if (.not.(this%east_boundary)) then
            !DIR$ PGAS DEFER_SYNC
            this%west_in_3d(:,:,:,:)[this%east_neighbor] = this%east_buffer_3d(:,:,:,:)
        endif
        if (.not.(this%west_boundary)) then
            !DIR$ PGAS DEFER_SYNC
            this%east_in_3d(:,:,:,:)[this%west_neighbor] = this%west_buffer_3d(:,:,:,:)
        endif

    end subroutine halo_3d_send_batch

    module subroutine halo_3d_retrieve_batch(this, exch_var_only, wait_timer)
        class(domain_t), intent(inout) :: this
        logical, intent(in) :: exch_var_only
        type(timer_t), optional,  intent(inout) :: wait_timer
        
        type(variable_t) :: var
        integer :: n, k_max

        if (present(wait_timer)) call wait_timer%start()
        sync images( this%neighbors )
        if (present(wait_timer)) call wait_timer%stop()

        call this%adv_vars%reset_iterator()
        call this%exch_vars%reset_iterator()
        n = 1
        ! Now iterate through the dictionary as long as there are more elements present
        if (.not.(exch_var_only)) then
            do while (this%adv_vars%has_more_elements())
                ! get the next variable
                var = this%adv_vars%next()
                if (var%three_d) then
                    if (.not.(this%north_boundary)) var%data_3d(this%its-1:this%ite+1,:,(this%jte+1):this%jme) = &
                            this%north_in_3d(n,1:(this%ite-this%its+3),:,:)
                    if (.not.(this%south_boundary)) var%data_3d(this%its-1:this%ite+1,:,this%jms:(this%jts-1)) = &
                            this%south_in_3d(n,1:(this%ite-this%its+3),:,:)
                    if (.not.(this%east_boundary)) var%data_3d((this%ite+1):this%ime,:,this%jts-1:this%jte+1) = &
                            this%east_in_3d(n,:,:,1:(this%jte-this%jts+3))
                    if (.not.(this%west_boundary)) var%data_3d(this%ims:(this%its-1),:,this%jts-1:this%jte+1) = &
                            this%west_in_3d(n,:,:,1:(this%jte-this%jts+3))
                    n = n+1
                endif
            enddo
        endif
    
        ! Now iterate through the exchange-only objects as long as there are more elements present
        do while (this%exch_vars%has_more_elements())
            ! get the next variable
            var = this%exch_vars%next()
            if (var%three_d) then
                k_max = ubound(var%data_3d,2)
                if (.not.(this%north_boundary)) var%data_3d(this%its-1:this%ite+1,1:k_max,(this%jte+1):this%jme) = &
                        this%north_in_3d(n,1:(this%ite-this%its+3),1:k_max,:)
                if (.not.(this%south_boundary)) var%data_3d(this%its-1:this%ite+1,1:k_max,this%jms:(this%jts-1)) = &
                        this%south_in_3d(n,1:(this%ite-this%its+3),1:k_max,:)
                if (.not.(this%east_boundary)) var%data_3d((this%ite+1):this%ime,1:k_max,this%jts-1:this%jte+1) = &
                        this%east_in_3d(n,:,1:k_max,1:(this%jte-this%jts+3))
                if (.not.(this%west_boundary)) var%data_3d(this%ims:(this%its-1),1:k_max,this%jts-1:this%jte+1) = &
                        this%west_in_3d(n,:,1:k_max,1:(this%jte-this%jts+3))
                n = n+1
            endif
        enddo

    
    end subroutine halo_3d_retrieve_batch

    module subroutine halo_2d_send_batch(this)
        class(domain_t), intent(inout) :: this
        type(variable_t) :: var
        integer :: n

        call this%exch_vars%reset_iterator()
        n = 1
        ! Now iterate through the exchange-only objects as long as there are more elements present
        do while (this%exch_vars%has_more_elements())
            ! get the next variable
            var = this%exch_vars%next()
            if (var%two_d) then
                if (.not.(this%north_boundary)) this%north_buffer_2d(n,1:(this%ite-this%its+3),:) = &
                        var%data_2d(this%its-1:this%ite+1,(this%jte-this%grid%halo_size+1):this%jte)
                if (.not.(this%south_boundary)) this%south_buffer_2d(n,1:(this%ite-this%its+3),:) = &
                        var%data_2d(this%its-1:this%ite+1,this%jts:(this%jts+this%grid%halo_size-1))
                if (.not.(this%east_boundary)) this%east_buffer_2d(n,:,1:(this%jte-this%jts+3)) = &
                        var%data_2d((this%ite-this%grid%halo_size+1):this%ite,this%jts-1:this%jte+1)
                if (.not.(this%west_boundary)) this%west_buffer_2d(n,:,1:(this%jte-this%jts+3)) = &
                        var%data_2d(this%its:(this%its+this%grid%halo_size)-1,this%jts-1:this%jte+1)

                n = n+1
            endif
        enddo

        if (.not.(this%north_boundary)) then
            !DIR$ PGAS DEFER_SYNC
            this%south_in_2d(:,:,:)[this%north_neighbor] = this%north_buffer_2d(:,:,:)
        endif
        if (.not.(this%south_boundary)) then
            !DIR$ PGAS DEFER_SYNC
            this%north_in_2d(:,:,:)[this%south_neighbor] = this%south_buffer_2d(:,:,:)
        endif
        if (.not.(this%east_boundary)) then
            !DIR$ PGAS DEFER_SYNC
            this%west_in_2d(:,:,:)[this%east_neighbor] = this%east_buffer_2d(:,:,:)
        endif
        if (.not.(this%west_boundary)) then
            !DIR$ PGAS DEFER_SYNC
            this%east_in_2d(:,:,:)[this%west_neighbor] = this%west_buffer_2d(:,:,:)
        endif

    end subroutine halo_2d_send_batch

    module subroutine halo_2d_retrieve_batch(this)
        class(domain_t), intent(inout) :: this
        type(variable_t) :: var
        integer :: n

        sync images( this%neighbors )

        call this%exch_vars%reset_iterator()
        n = 1    
        ! Now iterate through the exchange-only objects as long as there are more elements present
        do while (this%exch_vars%has_more_elements())
            ! get the next variable
            var = this%exch_vars%next()
            if (var%two_d) then
                if (.not.(this%north_boundary)) var%data_2d(this%its-1:this%ite+1,(this%jte+1):this%jme) = this%north_in_2d(n,1:(this%ite-this%its+3),:)
                if (.not.(this%south_boundary)) var%data_2d(this%its-1:this%ite+1,this%jms:(this%jts-1)) = this%south_in_2d(n,1:(this%ite-this%its+3),:)
                if (.not.(this%east_boundary)) var%data_2d((this%ite+1):this%ime,this%jts-1:this%jte+1) = this%east_in_2d(n,:,1:(this%jte-this%jts+3))
                if (.not.(this%west_boundary)) var%data_2d(this%ims:(this%its-1),this%jts-1:this%jte+1) = this%west_in_2d(n,:,1:(this%jte-this%jts+3))
                n = n+1
            endif
        enddo

    
    end subroutine halo_2d_retrieve_batch

    !> -------------------------------
    !! Send and get the data from all exch+adv objects to/from their neighbors (3D)
    !!
    !! -------------------------------
    module subroutine halo_3d_exchange_batch(this, send_timer, ret_timer, wait_timer, exch_var_only)
        class(domain_t), intent(inout) :: this
        type(timer_t),   optional, intent(inout) :: send_timer, ret_timer, wait_timer
        logical, optional, intent(in) :: exch_var_only
        
        logical :: exch_only
        
        exch_only = .False.
        if (present(exch_var_only)) exch_only = exch_var_only

        if (present(send_timer)) call send_timer%start()
        call this%halo_3d_send_batch(exch_var_only=exch_only)
        if (present(send_timer)) call send_timer%stop()

        if (present(ret_timer)) call ret_timer%start()
        if (present(wait_timer)) then
            call this%halo_3d_retrieve_batch(exch_only,wait_timer)
        else
            call this%halo_3d_retrieve_batch(exch_only)
        endif
        if (present(ret_timer)) call ret_timer%stop()
    end subroutine


    !> -------------------------------
    !! Send and get the data from all exch+adv objects to/from their neighbors (2D)
    !!
    !! -------------------------------
    module subroutine halo_2d_exchange_batch(this)
        class(domain_t), intent(inout) :: this
        
        call this%halo_2d_send_batch()

        call this%halo_2d_retrieve_batch()
    end subroutine


    !> -------------------------------
    !! Allocate and or initialize all domain variables if they have been requested
    !!
    !! -------------------------------
    subroutine create_variables(this, opt)
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: opt
        integer :: i,j

        integer :: ims, ime, jms, jme, kms, kme

        ims = this%grid%ims
        ime = this%grid%ime
        kms = this%grid%kms
        kme = this%grid%kme
        jms = this%grid%jms
        jme = this%grid%jme

        allocate( mod_temp_3d( ims:ime, kms:kme, jms:jme))
        allocate( surf_temp_1( ims:ime, jms:jme))
        allocate( surf_temp_2( ims:ime, jms:jme))

        if (this_image()==1) print *,"  Initializing variables"

        if (0<opt%vars_to_allocate( kVARS%u) )                          call setup(this%u,                        this%u_grid,   forcing_var=opt%parameters%uvar,       list=this%variables_to_force, force_boundaries=.False.)
        if (0<opt%vars_to_allocate( kVARS%u) )                          call setup(this%u_mass,                   this%grid)
        if (0<opt%vars_to_allocate( kVARS%v) )                          call setup(this%v,                        this%v_grid,   forcing_var=opt%parameters%vvar,       list=this%variables_to_force, force_boundaries=.False.)
        if (0<opt%vars_to_allocate( kVARS%v) )                          call setup(this%v_mass,                   this%grid)
        if (0<opt%vars_to_allocate( kVARS%w) )                          call setup(this%w,                        this%grid)
        if (0<opt%vars_to_allocate( kVARS%w_real) )                     call setup(this%w_real,                   this%grid,   forcing_var=opt%parameters%wvar,       list=this%variables_to_force, force_boundaries=.False. )
        if (0<opt%vars_to_allocate( kVARS%nsquared) )                   call setup(this%nsquared,                 this%grid )
        if (0<opt%vars_to_allocate( kVARS%water_vapor) )                call setup(this%water_vapor,              this%grid,     forcing_var=opt%parameters%qvvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%potential_temperature) )      call setup(this%potential_temperature,    this%grid,     forcing_var=opt%parameters%tvar,       list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%cloud_water) )                call setup(this%cloud_water_mass,         this%grid,     forcing_var=opt%parameters%qcvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%cloud_number_concentration))  call setup(this%cloud_number,             this%grid,     forcing_var=opt%parameters%qncvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%cloud_ice) )                  call setup(this%cloud_ice_mass,           this%grid,     forcing_var=opt%parameters%qivar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice_number_concentration))    call setup(this%cloud_ice_number,         this%grid,     forcing_var=opt%parameters%qnivar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%rain_in_air) )                call setup(this%rain_mass,                this%grid,     forcing_var=opt%parameters%qrvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%rain_number_concentration))   call setup(this%rain_number,              this%grid,     forcing_var=opt%parameters%qnrvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%snow_in_air) )                call setup(this%snow_mass,                this%grid,     forcing_var=opt%parameters%qsvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%snow_number_concentration) )  call setup(this%snow_number,              this%grid,     forcing_var=opt%parameters%qnsvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%graupel_in_air) )             call setup(this%graupel_mass,             this%grid,     forcing_var=opt%parameters%qgvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%graupel_number_concentration))call setup(this%graupel_number,           this%grid,     forcing_var=opt%parameters%qngvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice1_a))                      call setup(this%ice1_a,           this%grid,     forcing_var=opt%parameters%i1avar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice1_c))                      call setup(this%ice1_c,           this%grid,     forcing_var=opt%parameters%i1cvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice2_mass))                   call setup(this%ice2_mass,        this%grid,     forcing_var=opt%parameters%i2mvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice2_number))                 call setup(this%ice2_number,      this%grid,     forcing_var=opt%parameters%i2nvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice2_a))                      call setup(this%ice2_a,           this%grid,     forcing_var=opt%parameters%i2avar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice2_c))                      call setup(this%ice2_c,           this%grid,     forcing_var=opt%parameters%i2cvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice3_mass))                   call setup(this%ice3_mass,        this%grid,     forcing_var=opt%parameters%i3mvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice3_number))                 call setup(this%ice3_number,      this%grid,     forcing_var=opt%parameters%i3nvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice3_a))                      call setup(this%ice3_a,           this%grid,     forcing_var=opt%parameters%i3avar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%ice3_c))                      call setup(this%ice3_c,           this%grid,     forcing_var=opt%parameters%i3cvar,      list=this%variables_to_force, force_boundaries=.True.)
        if (0<opt%vars_to_allocate( kVARS%precipitation) )              call setup(this%accumulated_precipitation,this%grid2d, dtype=kDOUBLE )
        if (0<opt%vars_to_allocate( kVARS%convective_precipitation) )   call setup(this%accumulated_convective_pcp,this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%external_precipitation) )     call setup(this%external_precipitation,   this%grid2d,   forcing_var=opt%parameters%rain_var,  list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%snowfall) )                   call setup(this%accumulated_snowfall,     this%grid2d, dtype=kDOUBLE )
        if (0<opt%vars_to_allocate( kVARS%pressure) )                   call setup(this%pressure,                 this%grid,     forcing_var=opt%parameters%pvar,       list=this%variables_to_force, force_boundaries=.False.)
        if (0<opt%vars_to_allocate( kVARS%temperature) )                call setup(this%temperature,              this%grid )
        if (0<opt%vars_to_allocate( kVARS%exner) )                      call setup(this%exner,                    this%grid )
        if (0<opt%vars_to_allocate( kVARS%z) )                          call setup(this%z,                        this%grid )
        if (0<opt%vars_to_allocate( kVARS%dz_interface) )               call setup(this%dz_interface,             this%grid )
        if (0<opt%vars_to_allocate( kVARS%z_interface) )                call setup(this%z_interface,              this%grid8w )
        if (0<opt%vars_to_allocate( kVARS%dz) )                         call setup(this%dz_mass,                  this%grid )
        if (0<opt%vars_to_allocate( kVARS%density) )                    call setup(this%density,                  this%grid )
        if (0<opt%vars_to_allocate( kVARS%pressure_interface) )         call setup(this%pressure_interface,       this%grid8w )
        if (0<opt%vars_to_allocate( kVARS%graupel) )                    call setup(this%graupel,                  this%grid2d, dtype=kDOUBLE )
        if (0<opt%vars_to_allocate( kVARS%cloud_fraction) )             call setup(this%cloud_fraction,           this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%shortwave) )                  call setup(this%shortwave,                this%grid2d,   forcing_var=opt%parameters%swdown_var,  list=this%variables_to_force, force_boundaries=.False.)
        if (0<opt%vars_to_allocate( kVARS%shortwave_direct) )           call setup(this%shortwave_direct,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%shortwave_diffuse) )          call setup(this%shortwave_diffuse,        this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%shortwave_direct_above) )     call setup(this%shortwave_direct_above,   this%grid2d) !! MJ added
        if (0<opt%vars_to_allocate( kVARS%shortwave_total) )            call setup(this%shortwave_total,          this%grid2d) !! MJ added
        if (0<opt%vars_to_allocate( kVARS%longwave) )                   call setup(this%longwave,                 this%grid2d,   forcing_var=opt%parameters%lwdown_var,  list=this%variables_to_force, force_boundaries=.False.)
        if (0<opt%vars_to_allocate( kVARS%albedo) )                     call setup(this%albedo,                   this%grid_monthly )
        if (0<opt%vars_to_allocate( kVARS%vegetation_fraction) )        call setup(this%vegetation_fraction,      this%grid_monthly )
        if (0<opt%vars_to_allocate( kVARS%vegetation_fraction_max) )    call setup(this%vegetation_fraction_max,  this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%vegetation_fraction_out) )    call setup(this%vegetation_fraction_out,  this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%lai) )                        call setup(this%lai,                      this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%sai) )                        call setup(this%sai,                      this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%crop_type) )                  call setup(this%crop_type,                this%grid_croptype )
        if (0<opt%vars_to_allocate( kVARS%date_planting) )              call setup(this%date_planting,            this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%date_harvest) )               call setup(this%date_harvest,             this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%growing_season_gdd) )         call setup(this%growing_season_gdd,       this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_frac_total) )             call setup(this%irr_frac_total,           this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_frac_sprinkler) )         call setup(this%irr_frac_sprinkler,       this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_frac_micro) )             call setup(this%irr_frac_micro,           this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_frac_flood) )             call setup(this%irr_frac_flood,           this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_alloc_sprinkler) )        call setup(this%irr_alloc_sprinkler,      this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_alloc_micro) )            call setup(this%irr_alloc_micro,          this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_alloc_flood) )            call setup(this%irr_alloc_flood,          this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_evap_loss_sprinkler) )    call setup(this%irr_evap_loss_sprinkler,  this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_amt_sprinkler) )          call setup(this%irr_amt_sprinkler,        this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_amt_micro) )              call setup(this%irr_amt_micro,            this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%irr_amt_flood) )              call setup(this%irr_amt_flood,            this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%evap_heat_sprinkler) )        call setup(this%evap_heat_sprinkler,      this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%mass_ag_grain) )              call setup(this%mass_ag_grain,            this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%growing_degree_days) )        call setup(this%growing_degree_days,      this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%net_ecosystem_exchange) )     call setup(this%net_ecosystem_exchange,   this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%gross_primary_prod) )         call setup(this%gross_primary_prod,       this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%net_primary_prod) )           call setup(this%net_primary_prod,         this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%apar) )                       call setup(this%apar,                     this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%photosynthesis_total) )       call setup(this%photosynthesis_total,     this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%stomatal_resist_total) )      call setup(this%stomatal_resist_total,    this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%stomatal_resist_sun) )        call setup(this%stomatal_resist_sun,      this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%stomatal_resist_shade) )      call setup(this%stomatal_resist_shade,    this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%gecros_state) )               call setup(this%gecros_state,             this%grid_gecros )
        if (0<opt%vars_to_allocate( kVARS%canopy_water) )               call setup(this%canopy_water,             this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%canopy_water_ice) )           call setup(this%canopy_water_ice,         this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%canopy_water_liquid) )        call setup(this%canopy_water_liquid,      this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%canopy_vapor_pressure) )      call setup(this%canopy_vapor_pressure,    this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%canopy_temperature) )         call setup(this%canopy_temperature,       this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%canopy_fwet) )                call setup(this%canopy_fwet,              this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%veg_leaf_temperature) )       call setup(this%veg_leaf_temperature,     this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%ground_surf_temperature) )    call setup(this%ground_surf_temperature,  this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%frac_within_gap) )            call setup(this%frac_within_gap,          this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%frac_between_gap) )           call setup(this%frac_between_gap,         this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%ground_temperature_bare) )    call setup(this%ground_temperature_bare,  this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%ground_temperature_canopy) )  call setup(this%ground_temperature_canopy,this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snowfall_ground) )            call setup(this%snowfall_ground,          this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%rainfall_ground) )            call setup(this%rainfall_ground,          this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snow_water_equivalent) )      call setup(this%snow_water_equivalent,    this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snow_water_eq_prev) )         call setup(this%snow_water_eq_prev,       this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snow_albedo_prev) )           call setup(this%snow_albedo_prev,         this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snow_temperature) )           call setup(this%snow_temperature,         this%grid_snow )
        if (0<opt%vars_to_allocate( kVARS%snow_layer_depth) )           call setup(this%snow_layer_depth,         this%grid_snowsoil )
        if (0<opt%vars_to_allocate( kVARS%snow_layer_ice) )             call setup(this%snow_layer_ice,           this%grid_snow )
        if (0<opt%vars_to_allocate( kVARS%snow_layer_liquid_water) )    call setup(this%snow_layer_liquid_water,  this%grid_snow )
        if (0<opt%vars_to_allocate( kVARS%snow_age_factor) )            call setup(this%snow_age_factor,          this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snow_height) )                call setup(this%snow_height,              this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%sst) )                        call setup(this%sst,                      this%grid2d,   forcing_var=opt%parameters%sst_var,     list=this%variables_to_force, force_boundaries=.False.)
        if (0<opt%vars_to_allocate( kVARS%skin_temperature) )           call setup(this%skin_temperature,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_water_content) )         call setup(this%soil_water_content,       this%grid_soil)
        if (0<opt%vars_to_allocate( kVARS%eq_soil_moisture) )           call setup(this%eq_soil_moisture,         this%grid_soil)
        if (0<opt%vars_to_allocate( kVARS%smc_watertable_deep) )        call setup(this%smc_watertable_deep,      this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%recharge) )                   call setup(this%recharge,                 this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%recharge_deep) )              call setup(this%recharge_deep,            this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_temperature) )           call setup(this%soil_temperature,         this%grid_soil)
        if (0<opt%vars_to_allocate( kVARS%latitude) )                   call setup(this%latitude,                 this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%longitude) )                  call setup(this%longitude,                this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%u_latitude) )                 call setup(this%u_latitude,               this%u_grid2d)
        if (0<opt%vars_to_allocate( kVARS%u_longitude) )                call setup(this%u_longitude,              this%u_grid2d)
        if (0<opt%vars_to_allocate( kVARS%v_latitude) )                 call setup(this%v_latitude,               this%v_grid2d)
        if (0<opt%vars_to_allocate( kVARS%v_longitude) )                call setup(this%v_longitude,              this%v_grid2d)
        if (0<opt%vars_to_allocate( kVARS%terrain) )                    call setup(this%terrain,                  this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%terrain) )                    call setup(this%forcing_terrain,          this%grid2d) !,    forcing_var=opt%parameters%hgtvar, list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%sensible_heat) )              call setup(this%sensible_heat,            this%grid2d,   forcing_var=opt%parameters%shvar,     list=this%variables_to_force, force_boundaries=.False.)
        if (0<opt%vars_to_allocate( kVARS%latent_heat) )                call setup(this%latent_heat,              this%grid2d,   forcing_var=opt%parameters%lhvar,     list=this%variables_to_force, force_boundaries=.False.)
        if (0<opt%vars_to_allocate( kVARS%u_10m) )                      call setup(this%u_10m,                    this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%v_10m) )                      call setup(this%v_10m,                    this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%blk_ri) )                     call setup(this%Ri,                       this%grid)
        if (0<opt%vars_to_allocate( kVARS%froude) )                     call setup(this%froude,                   this%grid)
        if (0<opt%vars_to_allocate( kVARS%wind_alpha) )                 call setup(this%alpha,                    this%grid)

        if (0<opt%vars_to_allocate( kVARS%windspd_10m) )                call setup(this%windspd_10m,              this%grid2d) !! MJ added
        if (0<opt%vars_to_allocate( kVARS%coeff_momentum_drag) )        call setup(this%coeff_momentum_drag,      this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%chs) )                        call setup(this%chs,                      this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%chs2) )                       call setup(this%chs2,                     this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%cqs2) )                       call setup(this%cqs2,                     this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%br) )                         call setup(this%br,                       this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%QFX) )                        call setup(this%qfx,                      this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%psim) )                       call setup(this%psim,                     this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%psih) )                       call setup(this%psih,                     this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%fm) )                         call setup(this%fm,                       this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%fh) )                        call setup(this%fh,                      this%grid2d)
        
        if (0<opt%vars_to_allocate( kVARS%coeff_heat_exchange_3d) )     call setup(this%coeff_heat_exchange_3d,   this%grid)    ! for pbl ysu
        if (0<opt%vars_to_allocate( kVARS%coeff_momentum_exchange_3d) ) call setup(this%coeff_momentum_exchange_3d,this%grid)    ! for pbl ysu
        if (0<opt%vars_to_allocate( kVARS%hpbl) )                       call setup(this%hpbl,                     this%grid2d)    ! for pbl ysu
        if (0<opt%vars_to_allocate( kVARS%surface_rad_temperature) )    call setup(this%surface_rad_temperature,  this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%temperature_2m) )             call setup(this%temperature_2m,           this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%humidity_2m) )                call setup(this%humidity_2m,              this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%temperature_2m_veg) )         call setup(this%temperature_2m_veg,       this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%temperature_2m_bare) )        call setup(this%temperature_2m_bare,      this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%mixing_ratio_2m_veg) )        call setup(this%mixing_ratio_2m_veg,      this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%mixing_ratio_2m_bare) )       call setup(this%mixing_ratio_2m_bare,     this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%surface_pressure) )           call setup(this%surface_pressure,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%rad_absorbed_total) )         call setup(this%rad_absorbed_total,       this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%rad_absorbed_veg) )           call setup(this%rad_absorbed_veg,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%rad_absorbed_bare) )          call setup(this%rad_absorbed_bare,        this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%rad_net_longwave) )           call setup(this%rad_net_longwave,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%longwave_up) )                call setup(this%longwave_up,              this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%ground_heat_flux) )           call setup(this%ground_heat_flux,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%evap_canopy) )                call setup(this%evap_canopy,              this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%evap_soil_surface) )          call setup(this%evap_soil_surface,        this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%transpiration_rate) )         call setup(this%transpiration_rate,       this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%ch_veg) )                     call setup(this%ch_veg,                   this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%ch_veg_2m) )                  call setup(this%ch_veg_2m,                this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%ch_bare) )                    call setup(this%ch_bare,                  this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%ch_bare_2m) )                 call setup(this%ch_bare_2m,               this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%ch_under_canopy) )            call setup(this%ch_under_canopy,          this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%ch_leaf) )                    call setup(this%ch_leaf,                  this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%sensible_heat_veg) )          call setup(this%sensible_heat_veg,        this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%sensible_heat_bare) )         call setup(this%sensible_heat_bare,       this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%sensible_heat_canopy) )       call setup(this%sensible_heat_canopy,     this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%evap_heat_veg) )              call setup(this%evap_heat_veg,            this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%evap_heat_bare) )             call setup(this%evap_heat_bare,           this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%evap_heat_canopy) )           call setup(this%evap_heat_canopy,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%transpiration_heat) )         call setup(this%transpiration_heat,       this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%ground_heat_veg) )            call setup(this%ground_heat_veg,          this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%ground_heat_bare) )           call setup(this%ground_heat_bare,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%net_longwave_veg) )           call setup(this%net_longwave_veg,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%net_longwave_bare) )          call setup(this%net_longwave_bare,        this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%net_longwave_canopy) )        call setup(this%net_longwave_canopy,      this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%runoff_surface) )             call setup(this%runoff_surface,           this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%runoff_subsurface) )          call setup(this%runoff_subsurface,        this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_totalmoisture) )         call setup(this%soil_totalmoisture,       this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_deep_temperature) )      call setup(this%soil_deep_temperature,    this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%water_table_depth) )          call setup(this%water_table_depth,        this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%water_aquifer) )              call setup(this%water_aquifer,            this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%storage_gw) )                 call setup(this%storage_gw,               this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%storage_lake) )               call setup(this%storage_lake,             this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%roughness_z0) )               call setup(this%roughness_z0,             this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%mass_leaf) )                  call setup(this%mass_leaf,                this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%mass_root) )                  call setup(this%mass_root,                this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%mass_stem) )                  call setup(this%mass_stem,                this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%mass_wood) )                  call setup(this%mass_wood,                this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_carbon_fast) )           call setup(this%soil_carbon_fast,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_carbon_stable) )         call setup(this%soil_carbon_stable,       this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_texture_1) )             call setup(this%soil_texture_1,           this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_texture_2) )             call setup(this%soil_texture_2,           this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_texture_3) )             call setup(this%soil_texture_3,           this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_texture_4) )             call setup(this%soil_texture_4,           this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_sand_and_clay) )         call setup(this%soil_sand_and_clay,       this%grid_soilcomp)
        if (0<opt%vars_to_allocate( kVARS%re_cloud) )                   call setup(this%re_cloud,                 this%grid)
        if (0<opt%vars_to_allocate( kVARS%re_ice) )                     call setup(this%re_ice,                   this%grid)
        if (0<opt%vars_to_allocate( kVARS%re_snow) )                    call setup(this%re_snow,                  this%grid)
        if (0<opt%vars_to_allocate( kVARS%ice1_rho) )                   call setup(this%ice1_rho,                 this%grid)
        if (0<opt%vars_to_allocate( kVARS%ice1_phi) )                   call setup(this%ice1_phi,                 this%grid)
        if (0<opt%vars_to_allocate( kVARS%ice2_rho) )                   call setup(this%ice2_rho,                 this%grid)
        if (0<opt%vars_to_allocate( kVARS%ice2_phi) )                   call setup(this%ice2_phi,                 this%grid)
        if (0<opt%vars_to_allocate( kVARS%ice3_rho) )                   call setup(this%ice3_rho,                 this%grid)
        if (0<opt%vars_to_allocate( kVARS%ice3_phi) )                   call setup(this%ice3_phi,                 this%grid)
        if (0<opt%vars_to_allocate( kVARS%out_longwave_rad) )           call setup(this%out_longwave_rad,         this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%longwave_cloud_forcing) )     call setup(this%longwave_cloud_forcing,   this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%shortwave_cloud_forcing) )    call setup(this%shortwave_cloud_forcing,  this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%cosine_zenith_angle) )        call setup(this%cosine_zenith_angle,      this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%temperature_interface) )      call setup(this%temperature_interface,    this%grid8w)
        if (0<opt%vars_to_allocate( kVARS%land_emissivity) )            call setup(this%land_emissivity,          this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%tend_swrad) )                 call setup(this%tend_swrad,               this%grid)
        ! lake vars:
        if (0<opt%vars_to_allocate( kVARS%lake_depth) )                 call setup(this%lake_depth,              this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%t_lake3d) )                   call setup(this%t_lake3d,                this%grid_lake )
        if (0<opt%vars_to_allocate( kVARS%lake_icefrac3d) )             call setup(this%lake_icefrac3d,          this%grid_lake )
        if (0<opt%vars_to_allocate( kVARS%z_lake3d) )                   call setup(this%z_lake3d,                this%grid_lake )
        if (0<opt%vars_to_allocate( kVARS%dz_lake3d) )                  call setup(this%dz_lake3d,               this%grid_lake )
        if (0<opt%vars_to_allocate( kVARS%snl2d) )                      call setup(this%snl2d,                   this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%t_grnd2d) )                   call setup(this%t_grnd2d,                this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%t_soisno3d) )                 call setup(this%t_soisno3d,              this%grid_lake_soisno )
        if (0<opt%vars_to_allocate( kVARS%h2osoi_ice3d) )               call setup(this%h2osoi_ice3d,            this%grid_lake_soisno )
        if (0<opt%vars_to_allocate( kVARS%h2osoi_liq3d) )               call setup(this%h2osoi_liq3d,            this%grid_lake_soisno )
        if (0<opt%vars_to_allocate( kVARS%h2osoi_vol3d) )               call setup(this%h2osoi_vol3d,            this%grid_lake_soisno )
        if (0<opt%vars_to_allocate( kVARS%z3d) )                        call setup(this%z3d,                     this%grid_lake_soisno )
        if (0<opt%vars_to_allocate( kVARS%dz3d) )                       call setup(this%dz3d,                    this%grid_lake_soisno )
        if (0<opt%vars_to_allocate( kVARS%zi3d) )                       call setup(this%zi3d,                    this%grid_lake_soisno_1 )
        if (0<opt%vars_to_allocate( kVARS%watsat3d) )                   call setup(this%watsat3d,                this%grid_lake_soi )
        if (0<opt%vars_to_allocate( kVARS%csol3d) )                     call setup(this%csol3d,                  this%grid_lake_soi )
        if (0<opt%vars_to_allocate( kVARS%tkmg3d) )                     call setup(this%tkmg3d,                  this%grid_lake_soi )
        if (0<opt%vars_to_allocate( kVARS%tksatu3d) )                   call setup(this%tksatu3d,                this%grid_lake_soi )
        if (0<opt%vars_to_allocate( kVARS%tkdry3d) )                    call setup(this%tkdry3d,                 this%grid_lake_soi )
        if (0<opt%vars_to_allocate( kVARS%lakemask) )                   call setup(this%lakemask,                this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%savedtke12d) )                call setup(this%savedtke12d,             this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%lakedepth2d) )                call setup(this%lakedepth2d,             this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%ivt) )                        call setup(this%ivt,                     this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%iwv) )                        call setup(this%iwv,                     this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%iwl) )                        call setup(this%iwl,                     this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%iwi) )                        call setup(this%iwi,                     this%grid2d )

        ! integer variable_t types aren't available (yet...)
        if (0<opt%vars_to_allocate( kVARS%convective_precipitation) )   allocate(this%cu_precipitation_bucket  (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%precipitation) )              allocate(this%precipitation_bucket     (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%snowfall) )                   allocate(this%snowfall_bucket          (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%veg_type) )                   allocate(this%veg_type                 (ims:ime, jms:jme),          source=7)
        if (0<opt%vars_to_allocate( kVARS%soil_type) )                  allocate(this%soil_type                (ims:ime, jms:jme),          source=3)
        if (0<opt%vars_to_allocate( kVARS%land_mask) )                  allocate(this%land_mask                (ims:ime, jms:jme),          source=kLC_LAND)
        if (0<opt%vars_to_allocate( kVARS%snow_nlayers) )               allocate(this%snow_nlayers             (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%crop_category) )              allocate(this%crop_category            (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%irr_eventno_sprinkler) )      allocate(this%irr_eventno_sprinkler    (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%irr_eventno_micro) )          allocate(this%irr_eventno_micro        (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%irr_eventno_flood) )          allocate(this%irr_eventno_flood        (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%plant_growth_stage) )         allocate(this%plant_growth_stage       (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%kpbl) )                       allocate(this%kpbl                     (ims:ime, jms:jme),          source=0) ! for pbl ysu

        ! tendency variables that don't need to be output... maybe these should be set up the same way
        if (0<opt%vars_to_allocate( kVARS%tend_qv_adv) )                allocate(this%tend%qv_adv(ims:ime, kms:kme, jms:jme),   source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qv_pbl) )                allocate(this%tend%qv_pbl(ims:ime, kms:kme, jms:jme),   source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qv) )                    allocate(this%tend%qv(ims:ime, kms:kme, jms:jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_th) )                    allocate(this%tend%th(ims:ime, kms:kme, jms:jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_th_pbl) )                allocate(this%tend%th_pbl(ims:ime, kms:kme, jms:jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qc_pbl) )                allocate(this%tend%qc_pbl(ims:ime, kms:kme, jms:jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qi_pbl) )                allocate(this%tend%qi_pbl(ims:ime, kms:kme, jms:jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qc) )                    allocate(this%tend%qc(ims:ime, kms:kme, jms:jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qi) )                    allocate(this%tend%qi(ims:ime, kms:kme, jms:jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qs) )                    allocate(this%tend%qs(ims:ime, kms:kme, jms:jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qr) )                    allocate(this%tend%qr(ims:ime, kms:kme, jms:jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_u) )                     allocate(this%tend%u(ims:ime, kms:kme, jms:jme),        source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_v) )                     allocate(this%tend%v(ims:ime, kms:kme, jms:jme),        source=0.0)

        if (0<opt%vars_to_allocate( kVARS%ustar) )                      allocate(this%ustar(ims:ime, jms:jme),   source=0.1)

        if (0<opt%vars_to_allocate( kVARS%znu) )                        allocate(this%znu(kms:kme),   source=0.0)
        if (0<opt%vars_to_allocate( kVARS%znw) )                        allocate(this%znw(kms:kme),   source=0.0)

        !! MJ added for needed new vars for FSM
        !! note that, in lsm_driver, it is alreaddy decieded if we need these vars or not..so it should be here..
        !if (0<opt%vars_to_allocate( kVARS%FSM_slopemu) )               allocate(this%FSM_slopemu                (ims:ime, jms:jme))!,          source=0.0)        
        if (0<opt%vars_to_allocate( kVARS%runoff_tstep) )               call setup(this%runoff_tstep,     this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%Tsnow) )                      call setup(this%Tsnow,      this%grid_snow)        
        if (0<opt%vars_to_allocate( kVARS%Sice) )                       call setup(this%Sice,       this%grid_snow)        
        if (0<opt%vars_to_allocate( kVARS%Sliq) )                       call setup(this%Sliq,       this%grid_snow)        
        if (0<opt%vars_to_allocate( kVARS%Ds) )                         call setup(this%Ds,         this%grid_snow)        
        if (0<opt%vars_to_allocate( kVARS%fsnow) )                      call setup(this%fsnow,      this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%Nsnow) )                      call setup(this%Nsnow,      this%grid2d)   
        if (0<opt%vars_to_allocate( kVARS%dm_salt) )                    call setup(this%dm_salt,    this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%dm_susp) )                    call setup(this%dm_susp,    this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%dm_subl) )                    call setup(this%dm_subl,    this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%dm_slide) )                   call setup(this%dm_slide,   this%grid2d)        

        !!
        if (0<opt%vars_to_allocate( kVARS%rainfall_tstep) )             call setup(this%rainfall_tstep,     this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%snowfall_tstep) )             call setup(this%snowfall_tstep,     this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%meltflux_out_tstep) )         call setup(this%meltflux_out_tstep, this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%Sliq_out) )                   call setup(this%Sliq_out,           this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%slope) )                      call setup(this%slope,              this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%slope_angle) )                call setup(this%slope_angle,        this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%aspect_angle) )               call setup(this%aspect_angle,       this%grid2d)        
        if (0<opt%vars_to_allocate( kVARS%svf) )                        call setup(this%svf,                this%grid2d) 
        if (0<opt%vars_to_allocate( kVARS%factor_p) )                   call setup(this%factor_p,           this%grid2d) 
        if (0<opt%vars_to_allocate( kVARS%ridge_dist) )                 call setup(this%ridge_dist,         this%grid2d) 
        if (0<opt%vars_to_allocate( kVARS%valley_dist) )                call setup(this%valley_dist,        this%grid2d) 
        if (0<opt%vars_to_allocate( kVARS%ridge_drop) )                 call setup(this%ridge_drop,         this%grid2d) 
        if (0<opt%vars_to_allocate( kVARS%hlm) )                        call setup(this%hlm,                this%grid_hlm) 
        if (0<opt%vars_to_allocate( kVARS%shd) )                        call setup(this%shd,                this%grid2d) 

    end subroutine

    !> -------------------------------
    !! Setup a regular variable.
    !!
    !! Initializes the variable
    !! including the forcing_variable if it was set
    !! and adds that variable to the list of variables that has forcing data if the list is supplied
    !! and the forcing_var is both present and not blank ("")
    !!
    !! -------------------------------
    subroutine setup_var(var, grid, forcing_var, list, force_boundaries, dtype)
        implicit none
        type(variable_t),   intent(inout) :: var
        type(grid_t),       intent(in)    :: grid
        character(len=*),   intent(in),   optional :: forcing_var
        type(var_dict_t),   intent(inout),optional :: list
        logical,            intent(in),   optional :: force_boundaries
        integer,            intent(in),   optional :: dtype

        if (present(forcing_var)) then
            if (present(dtype)) then
                call var%initialize(grid, forcing_var=forcing_var, dtype=dtype)
            else
                call var%initialize(grid, forcing_var=forcing_var)
            endif

            if (present(list)) then
                if (Len(Trim(forcing_var)) /= 0) then
                    if (present(force_boundaries)) var%force_boundaries = force_boundaries
                    call list%add_var(forcing_var, var)
                endif
            endif
        else

            if (present(dtype)) then
                call var%initialize(grid, dtype=dtype)
            else
                call var%initialize(grid)
            endif
        endif

    end subroutine


    !> -------------------------------
    !! Setup an exchangeable variable.
    !!
    !! Initializes the variable
    !! including the forcing_variable if it was set
    !! and adds that variable to the list of variables that has forcing data if the list is supplied
    !! and the forcing_var is both present and not blank ("")
    !!
    !! -------------------------------
    subroutine setup_exch(var, grid, forcing_var, list, force_boundaries)
        implicit none
        type(exchangeable_t),   intent(inout) :: var
        type(grid_t),           intent(in)    :: grid
        character(len=*),       intent(in),   optional :: forcing_var
        type(var_dict_t),       intent(inout),optional :: list
        logical,                intent(in),   optional :: force_boundaries

        if (present(forcing_var)) then
            call var%initialize(grid, forcing_var=forcing_var)

            if (present(list)) then
                if (Len(Trim(forcing_var)) /= 0) then
                    if (present(force_boundaries)) var%meta_data%force_boundaries = force_boundaries
                    call list%add_var(forcing_var, var%meta_data)
                endif
            endif
        else

            call var%initialize(grid)
        endif

    end subroutine


    !> ---------------------------------
    !! Read the core model variables from disk
    !!
    !! Reads Terrain, lat, lon and u/v lat/lon on the high-res domain grid
    !! Passing data between images and disk is handled by io_read
    !!
    !! ---------------------------------
    subroutine read_core_variables(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options
        real, allocatable :: temporary_data(:,:), temp_offset(:,:)

        ! Read the terrain data
        call io_read(options%parameters%init_conditions_file,   &
                       options%parameters%hgt_hi,                 &
                       temporary_data)
        this%terrain%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
        
        !while we have global terrain loaded, pass to split_topography
        if (options%parameters%sleve) call split_topography(this, temporary_data, options)  ! here h1 and h2 are calculated
        
        allocate(this%neighbor_terrain(this%ihs:this%ihe, this%jhs:this%jhe), &
                    source=temporary_data(this%ihs:this%ihe, this%jhs:this%jhe))
        allocate(temp_offset(1:this%grid%ide+1,1:this%grid%jde+1))

        if ( (options%physics%windtype == kWIND_LINEAR) .or. (options%physics%windtype == kLINEAR_OBRIEN_WINDS) .or. &
             (options%physics%windtype == kLINEAR_ITERATIVE_WINDS) ) then
            this%global_terrain = temporary_data ! save the global terrain map for the linear wind solution
        end if


        ! Read the latitude data
        call io_read(options%parameters%init_conditions_file,   &
                       options%parameters%lat_hi,                 &
                       temporary_data)

        call make_2d_y(temporary_data, this%grid%ims, this%grid%ime)
        this%latitude%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
        ! allocate(this%latitude_global, source=temporary_data)

        ! Read the longitude data
        call io_read(options%parameters%init_conditions_file,   &
                       options%parameters%lon_hi,                 &
                       temporary_data)
        call make_2d_x(temporary_data, this%grid%jms, this%grid%jme)
        this%longitude%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
        ! allocate(this%longitude_global, source=temporary_data)

        !-----------------------------------------
        !
        ! Handle staggered lat/lon grids, straightfoward if ulat/ulon are supplied
        ! If not, then read in mass grid lat/lon and stagger them
        !
        !-----------------------------------------
        ! Read the u-grid longitude data if specified, other wise interpolate from mass grid
        if (options%parameters%ulon_hi /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%ulon_hi,                &
                           temporary_data)

            call make_2d_y(temporary_data, 1, this%jde)
            this%u_longitude%data_2d = temporary_data(this%u_grid%ims:this%u_grid%ime,this%u_grid%jms:this%u_grid%jme)
        else
            ! load the mass grid data again to get the full grid
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%lon_hi,                 &
                           temporary_data)

            call make_2d_y(temporary_data, 1, this%jde)
            call array_offset_x(temporary_data, temp_offset)
            this%u_longitude%data_2d = temp_offset(this%u_grid%ims:this%u_grid%ime,this%u_grid%jms:this%u_grid%jme)
        endif

        ! Read the u-grid latitude data if specified, other wise interpolate from mass grid
        if (options%parameters%ulat_hi /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%ulat_hi,                &
                           temporary_data)

            call make_2d_x(temporary_data, 1, this%ide+1)
            this%u_latitude%data_2d = temporary_data(this%u_grid%ims:this%u_grid%ime,this%u_grid%jms:this%u_grid%jme)
        else
            ! load the mass grid data again to get the full grid
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%lat_hi,                 &
                           temporary_data)

            call make_2d_x(temporary_data, 1, this%ide+1)
            call array_offset_x(temporary_data, temp_offset)
            this%u_latitude%data_2d = temp_offset(this%u_grid%ims:this%u_grid%ime,this%u_grid%jms:this%u_grid%jme)
        endif

        ! Read the v-grid longitude data if specified, other wise interpolate from mass grid
        if (options%parameters%vlon_hi /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%vlon_hi,                &
                           temporary_data)

            call make_2d_y(temporary_data, 1, this%jde+1)
            this%v_longitude%data_2d = temporary_data(this%v_grid%ims:this%v_grid%ime,this%v_grid%jms:this%v_grid%jme)
        else
            ! load the mass grid data again to get the full grid
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%lon_hi,                 &
                           temporary_data)

            call make_2d_y(temporary_data, 1, this%jde+1)
            call array_offset_y(temporary_data, temp_offset)
            this%v_longitude%data_2d = temp_offset(this%v_grid%ims:this%v_grid%ime,this%v_grid%jms:this%v_grid%jme)
        endif

        ! Read the v-grid latitude data if specified, other wise interpolate from mass grid
        if (options%parameters%vlat_hi /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%vlat_hi,                &
                           temporary_data)

            call make_2d_x(temporary_data, 1, this%ide)
            this%v_latitude%data_2d = temporary_data(this%v_grid%ims:this%v_grid%ime,this%v_grid%jms:this%v_grid%jme)
        else
            ! load the mass grid data again to get the full grid
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%lat_hi,                 &
                           temporary_data)

            call make_2d_x(temporary_data, 1, this%ide)
            call array_offset_y(temporary_data, temp_offset)
            this%v_latitude%data_2d = temp_offset(this%v_grid%ims:this%v_grid%ime,this%v_grid%jms:this%v_grid%jme)
        endif

        if (this_image()==1) write(*,*) "  Finished reading core domain variables"

    end subroutine


    !> -------------------------------
    !! Setup a single Geographic structure given a latitude, longitude, and z array
    !!
    !! -------------------------------
    subroutine setup_geo(geo, latitude, longitude, longitude_system, z)
        implicit none
        type(interpolable_type),  intent(inout) :: geo
        real,                     intent(in)    :: latitude(:,:)
        real,                     intent(in)    :: longitude(:,:)
        integer,                  intent(in)    :: longitude_system
        real, optional,           intent(in)    :: z(:,:,:)
        if (allocated(geo%lat)) deallocate(geo%lat)
        allocate( geo%lat, source=latitude)

        if (allocated(geo%lon)) deallocate(geo%lon)
        allocate( geo%lon, source=longitude)

        if (present(z)) then
            if (allocated(geo%z)) deallocate(geo%z)
            allocate( geo%z, source=z)
        endif
        ! This makes 2D variables out of lat/lon if they come in as 1D variables
        ! This also puts the longitudes onto a 0-360 if they are -180-180 (important for Alaska)
        ! Though if working in Europe the -180-180 grid is better ideally the optimal value should be checked.
        ! and good luck if you want to work over the poles...
        call standardize_coordinates(geo, longitude_system)

    end subroutine


    function find_flat_model_level(options, nz, dz) result(max_level)
        implicit none
        type(options_t), intent(in) :: options
        integer,         intent(in) :: nz
        real,            intent(in) :: dz(:)
        integer :: max_level

        integer :: j
        real :: height

        if (options%parameters%flat_z_height > nz) then
            if (this_image()==1) write(*,*) "    Treating flat_z_height as specified in meters above mean terrain height: ", options%parameters%flat_z_height," meters"
            height = 0
            do j = 1, nz
                if (height <= options%parameters%flat_z_height) then
                    height = height + dz(j)
                    max_level = j
                endif
            enddo

        elseif (options%parameters%flat_z_height <= 0) then
            if (this_image()==1) write(*,*) "    Treating flat_z_height as counting levels down from the model top: ", options%parameters%flat_z_height," levels"
            max_level = nz + options%parameters%flat_z_height

        else
            if (this_image()==1) write(*,*) "    Treating flat_z_height as counting levels up from the ground: ", options%parameters%flat_z_height," levels"
            max_level = options%parameters%flat_z_height
        endif

    end function find_flat_model_level





    subroutine allocate_z_arrays(this, options)
        implicit none
        class(domain_t),  intent(inout)  :: this
        class(options_t), intent(in)     :: options

        allocate(this%jacobian(this% ims : this% ime, &
                                    this% kms : this% kme, &
                                    this% jms : this% jme) )

        allocate(this%jacobian_u(this% ims : this% ime+1, &
                                    this% kms : this% kme, &
                                    this% jms : this% jme) )

        allocate(this%jacobian_v(this% ims : this% ime, &
                                    this% kms : this% kme, &
                                    this% jms : this% jme+1) )

        allocate(this%jacobian_w(this% ims : this% ime, &
                                    this% kms : this% kme, &
                                    this% jms : this% jme) )
                                                                
        allocate(this%dzdx(this% ims : this% ime, &
                           this% kms : this% kme, &
                           this% jms : this% jme) )

        allocate(this%dzdy(this% ims : this% ime, &
                           this% kms : this% kme, &
                           this% jms : this% jme) )
                           
        allocate(this%dzdx_u(this% ims : this% ime+1, &
                           this% kms : this% kme, &
                           this% jms : this% jme) )

        allocate(this%dzdy_v(this% ims : this% ime, &
                           this% kms : this% kme, &
                           this% jms : this% jme+1) )

        allocate(this%dz_scl( this%kms : this%kme))

        if (options%physics%windtype == kLINEAR_ITERATIVE_WINDS .or. options%physics%windtype == kITERATIVE_WINDS) then
            allocate(this%dzdxz(this% ims : this% ime, &
                                this% kms : this% kme, &
                                this% jms : this% jme) )
            allocate(this%dzdyz(this% ims : this% ime, &
                                this% kms : this% kme, &
                                this% jms : this% jme) )
        endif

        if ( (options%physics%windtype == kWIND_LINEAR) .or. (options%physics%windtype == kLINEAR_OBRIEN_WINDS) .or. &
             (options%physics%windtype == kLINEAR_ITERATIVE_WINDS) ) then
             
            allocate(this%global_z_interface(this% ids : this% ide,   &
                                             this% kds : this% kde+1, &
                                             this% jds : this% jde)   )

            allocate(this%global_dz_interface(this% ids : this% ide,   &
                                              this% kds : this% kde,   &
                                              this% jds : this% jde)   )
        else
            allocate(this%global_z_interface(this% ihs : this% ihe,   &
                                             this% khs : this% khe+1, &
                                             this% jhs : this% jhe)   )

            allocate(this%global_dz_interface(this% ihs : this% ihe,   &
                                              this% khs : this% khe,   &
                                              this% jhs : this% jhe)   )
        endif


        allocate(this%sintheta( this% ims : this% ime, &
                                this% jms : this% jme) )

        allocate(this%costheta( this% ims : this% ime, &
                                this% jms : this% jme) )

        allocate(this%advection_dz( this%ims:this%ime,  &
                                    this%kms:this%kme,  &
                                    this%jms:this%jme) )


    end subroutine allocate_z_arrays


    !> -------------------------------
    !! Setup the SLEVE vertical grid structure.
    !!   This basically entails 2 transformations: First a linear one so that sum(dz) ranges from 0 to smooth_height H.
    !!   (boundary cnd (3) in Schär et al 2002)  Next, the nonlinear SLEVE transformation
    !!    eqn (2) from Leuenberger et al 2009 z_sleve = Z + terrain * sinh((H/s)**n - (Z/s)**n) / SINH((H/s)**n) (for both smallscale and largescale terrain)
    !!   Here H is the model top or (flat_z_height in m), s controls how fast the terrain decays
    !!   and n controls the compression throughout the column (this last factor was added by Leuenberger et al 2009)
    !!   References: Leuenberger et al 2009 "A Generalization of the SLEVE Vertical Coordinate"
    !!               Schär et al 2002 "A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models"
    !!
    !! N.B. flat dz height != 0 makes little sense here? But works (?)
    !! -------------------------------
    subroutine setup_sleve(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options

        real, allocatable :: temp(:,:,:), gamma_n(:), neighbor_jacobian(:,:,:), neighbor_z(:,:,:)
        integer :: i, max_level
        real :: s, n, s1, s2, gamma, gamma_min
        real :: b1_i, b1_mass, db1_i, db1_mass, b2_i, b2_mass, db2_i, db2_mass
        
        associate(ims => this%ims,      ime => this%ime,                        &
            jms => this%jms,      jme => this%jme,                        &
            kms => this%kms,      kme => this%kme,                        &
            z                     => this%z%data_3d,                      &
            z_u                   => this%geo_u%z,                        &
            z_v                   => this%geo_v%z,                        &
            z_interface           => this%z_interface%data_3d,            &
            nz                    => options%parameters%nz,               &
            dz                    => options%parameters%dz_levels,        &
            dz_mass               => this%dz_mass%data_3d,                &
            dz_interface          => this%dz_interface%data_3d,           &
            terrain               => this%terrain%data_2d,                &
            h1                    => this%h1,                &
            h2                    => this%h2,                &
            h1_u                  => this%h1_u,                &
            h2_u                  => this%h2_u,                &
            h1_v                  => this%h1_v,                &
            h2_v                  => this%h2_v,                &
            global_z_interface    => this%global_z_interface,             &
            global_dz_interface   => this%global_dz_interface,            &
            neighbor_terrain      => this%neighbor_terrain,               &
            jacobian_u            => this%jacobian_u,                     &
            jacobian_v            => this%jacobian_v,                     &
            jacobian_w            => this%jacobian_w,                     &
            dzdx                  => this%dzdx,                           &
            dzdy                  => this%dzdy,                           &
            dzdx_u                => this%dzdx_u,                         &
            dzdy_v                => this%dzdy_v,                         &
            jacobian              => this%jacobian,                       &
            smooth_height         => this%smooth_height,                  &
            dz_scl                => this%dz_scl)

            ! Still not 100% convinced this works well in cases other than flat_z_height = 0 (w sleve). So for now best to keep at 0 when using sleve?
            max_level = find_flat_model_level(options, nz, dz)

            smooth_height = sum(dz(1:max_level))!+dz(max_level)*0.5

            ! Terminology from Schär et al 2002, Leuenberger 2009: (can be simpliied later on, but for clarity)
            s1 = smooth_height / options%parameters%decay_rate_L_topo
            s2 = smooth_height / options%parameters%decay_rate_S_topo
            n  =  options%parameters%sleve_n 

            ! Scale dz with smooth_height/sum(dz(1:max_level)) before calculating sleve levels.
            dz_scl(:)   =   dz(1:nz) !*  smooth_height / sum(dz(1:max_level))  ! this leads to a jump in dz thickness at max_level+1. Not sure if this is a problem.

            ! - - -   calculate invertibility parameter gamma (Schär et al 2002 eqn 20):  - - - - - -
            gamma  =  1  -  MAXVAL(h1)/s1 * COSH(smooth_height/s1)/SINH(smooth_height/s1) &
                          - MAXVAL(h2)/s2 * COSH(smooth_height/s2)/SINH(smooth_height/s2)

            ! with the new (leuenberger et al 2010) Sleve formulation, the inveribiltiy criterion is as follows:
            ! ( Although an argument could be made to calculate this on the offset (u/v) grid b/c that is most
            !   relevant for advection? In reality this is probably a sufficient approximation, as long as we
            !   aren't pushing the gamma factor too close to zero )
            allocate(gamma_n(this%kds : this%kde+1))
            i=kms
            gamma_n(i) =  1                                                     &
                - MAXVAL(h1) * n/(s1**n)                                        &
                * COSH((smooth_height/s1)**n) / SINH((smooth_height/s1)**n)     &
                - MAXVAL(h2) * n/(s2**n)                                        &
                * COSH((smooth_height/s2)**n) / SINH((smooth_height/s2)**n)

            do i = this%grid%kds, this%grid%kde
                gamma_n(i+1)  =  1                                    &    ! # for i != kds !!
                - MAXVAL(h1) * n/(s1**n) * sum(dz_scl(1:i))**(n-1)                                             &
                * COSH((smooth_height/s1)**n -(sum(dz_scl(1:i))/s1)**n ) / SINH((smooth_height/s1)**n)    &
                - MAXVAL(h2) * n/(s2**n) *  sum(dz_scl(1:i))**(n-1)                                            &
                * COSH((smooth_height/s2)**n -(sum(dz_scl(1:i))/s2)**n ) / SINH((smooth_height/s2)**n)
            enddo

            if (n==1) then
                gamma_min = gamma
            else
                gamma_min = MINVAL(gamma_n)
            endif


            ! For reference: COSMO1 operational setting (but model top is at ~22000 masl):
            !    Decay Rate for Large-Scale Topography: svc1 = 10000.0000
            !    Decay Rate for Small-Scale Topography: svc2 =  3300.0000
            if ((this_image()==1)) then
                write(*,*) "    Using a SLEVE coordinate with a Decay height for Large-Scale Topography: (s1) of ", s1, " m."
                write(*,*) "    Using a SLEVE coordinate with a Decay height for Small-Scale Topography: (s2) of ", s2, " m."
                write(*,*) "    Using a sleve_n of ", options%parameters%sleve_n
                write(*,*) "    Smooth height is ", smooth_height, "m.a.s.l     (model top ", sum(dz(1:nz)), "m.a.s.l.)"
                write(*,*) "    invertibility parameter gamma is: ", gamma_min
                if(gamma_min <= 0) write(*,*) " CAUTION: coordinate transformation is not invertible (gamma <= 0 ) !!! reduce decay rate(s), and/or increase flat_z_height!"
                ! if(options%parameters%debug)  write(*,*) "   (for (debugging) reference: 'gamma(n=1)'= ", gamma,")"
                write(*,*) ""
            endif

            ! use temp to store global z-interface so that global-jacobian can be calculated

            if ( (options%physics%windtype == kWIND_LINEAR) .or. (options%physics%windtype == kLINEAR_OBRIEN_WINDS) .or. &
                 (options%physics%windtype == kLINEAR_ITERATIVE_WINDS) ) then
                 
                allocate(temp(this%ids:this%ide, this%kds:this%kde, this%jds:this%jde))
                temp(:,kms,:)   = this%global_terrain
            else
                allocate(temp(this%ihs:this%ihe, this%khs:this%khe, this%jhs:this%jhe))
                temp(:,kms,:)   = neighbor_terrain
            endif
            
            allocate(neighbor_jacobian(this%ihs:this%ihe, this%khs:this%khe, this%jhs:this%jhe))
            allocate(neighbor_z(this%ihs:this%ihe, this%khs:this%khe, this%jhs:this%jhe))
                                    
            ! - - - - -  k levels  - - - - -
            do i = this%grid%kms, this%grid%kme

                if (i==kms) then
                    b1_i = SINH( (smooth_height/s1)**n - (dz_scl(i)/s1)**n ) / SINH((smooth_height/s1)**n)
                    b2_i = SINH( (smooth_height/s2)**n - (dz_scl(i)/s2)**n ) / SINH((smooth_height/s2)**n)
                    b1_mass = SINH( (smooth_height/s1)**n -  ( (dz_scl(i)/2) /s1)**n ) / SINH((smooth_height/s1)**n)
                    b2_mass = SINH( (smooth_height/s2)**n -  ( (dz_scl(i)/2) /s2)**n ) / SINH((smooth_height/s2)**n)

                    db1_i = -n/(s1**n) * dz_scl(i)**(n-1) * COSH((smooth_height/s1)**n - & 
                            (dz_scl(i)/s1)**n ) / SINH((smooth_height/s1)**n)
                    db2_i = -n/(s2**n) * dz_scl(i)**(n-1) * COSH((smooth_height/s2)**n - & 
                            (dz_scl(i)/s2)**n ) / SINH((smooth_height/s2)**n)

                    db1_mass = -n/(s1**n) * (dz_scl(i)/2)**(n-1) * COSH((smooth_height/s1)**n - &
                            ((dz_scl(i)/2)/s1)**n ) / SINH((smooth_height/s1)**n)
                    db2_mass = -n/(s2**n) * (dz_scl(i)/2)**(n-1) * COSH((smooth_height/s2)**n - &
                            ((dz_scl(i)/2)/s2)**n ) / SINH((smooth_height/s2)**n)

                    temp(:,i+1,:)  = dz_scl(i) + h1*b1_i + h2*b2_i

                    global_dz_interface(:,i,:)  =  temp(:,i+1,:) - temp(:,i,:)  ! same for higher k
                    global_z_interface(:,i,:)  = temp(:,i,:)

                    dz_mass(:,i,:)       = global_dz_interface(ims:ime,i,jms:jme) / 2           ! Diff for k=1            

                    ! ! - - - - -   u/v grid calculations for lowest level (i=kms)  - - - - -
                    ! ! for the u and v grids, z(1) was already initialized with terrain.
                    ! ! but the first level needs to be offset, and the rest of the levels need to be created
                    ! ! BK: So if z_u is already offset in the u dir, but not in the z dir, we can say that
                    ! !     z_u(:,1,:) is the terrain on the u grid, and it needs to be offset in the z-dir
                    ! !     to reach mass levels (so by dz[i]/2)

                    neighbor_z(:,i,:)  = (dz_scl(i)/2)  + h1*b1_mass + h2*b2_mass
                    z_u(:,i,:)   = (dz_scl(i)/2) + h1_u*b1_mass + h2_u*b2_mass
                    z_v(:,i,:)   = (dz_scl(i)/2) + h1_v*b1_mass + h2_v*b2_mass

                else if(i>kms) then
                    if(i<=max_level) then

                        b1_i = SINH( (smooth_height/s1)**n - (sum(dz_scl(1:i))/s1)**n ) / SINH((smooth_height/s1)**n)
                        b2_i = SINH( (smooth_height/s2)**n - (sum(dz_scl(1:i))/s2)**n ) / SINH((smooth_height/s2)**n)
                        b1_mass = SINH( (smooth_height/s1)**n -  ( (sum(dz_scl(1:(i-1)))+dz_scl(i)/2) /s1)**n ) / SINH((smooth_height/s1)**n)
                        b2_mass = SINH( (smooth_height/s2)**n -  ( (sum(dz_scl(1:(i-1)))+dz_scl(i)/2) /s2)**n ) / SINH((smooth_height/s2)**n)


                        db1_i = -n/(s1**n) * sum(dz_scl(1:i))**(n-1) * COSH((smooth_height/s1)**n - & 
                                (sum(dz_scl(1:i))/s1)**n ) / SINH((smooth_height/s1)**n)
                        db2_i = -n/(s2**n) * sum(dz_scl(1:i))**(n-1) * COSH((smooth_height/s2)**n - & 
                                (sum(dz_scl(1:i))/s2)**n ) / SINH((smooth_height/s2)**n)
                                
                        db1_mass = -n/(s1**n) * (sum(dz_scl(1:(i-1)))+dz_scl(i)/2)**(n-1) * COSH((smooth_height/s1)**n - &
                                ((sum(dz_scl(1:(i-1)))+dz_scl(i)/2)/s1)**n ) / SINH((smooth_height/s1)**n)
                        db2_mass = -n/(s2**n) * (sum(dz_scl(1:(i-1)))+dz_scl(i)/2)**(n-1) * COSH((smooth_height/s2)**n - &
                                ((sum(dz_scl(1:(i-1)))+dz_scl(i)/2)/s2)**n ) / SINH((smooth_height/s2)**n)

                        if (i==this%grid%kme) then  ! if we are at the model top i+1 is not defined
                            global_dz_interface(:,i,:)  =  smooth_height - temp(:,i,:)
                        else
                            temp(:,i+1,:)  = sum(dz_scl(1:i)) + h1*b1_i + h2*b2_i 
                            global_dz_interface(:,i,:)  =  temp(:,i+1,:) - temp(:,i,:)
                        endif

                        global_z_interface(:,i,:)  = global_z_interface(:,i-1,:) + global_dz_interface(:,i-1,:)

                        neighbor_z(:,i,:)  = (sum(dz_scl(1:(i-1))) + dz_scl(i)/2)  + h1*b1_mass + h2*b2_mass  
                        z_u(:,i,:)   = (sum(dz_scl(1:(i-1))) + dz_scl(i)/2) + h1_u*b1_mass + h2_u*b2_mass  
                        z_v(:,i,:)   = (sum(dz_scl(1:(i-1))) + dz_scl(i)/2) + h1_v*b1_mass + h2_v*b2_mass  

                        if ( ANY(global_z_interface(:,i,:)<0) ) then   ! Eror catching. Probably good to engage.
                        if (this_image()==1) then
                            write(*,*) "Error: dz_interface below zero (for level  ",i,")"
                            write(*,*)  "min max dz_interface: ",MINVAL(global_z_interface(:,i,:)),MAXVAL(global_z_interface(:,i,:))
                            error stop
                        endif
                        else if ( ANY(global_z_interface(:,i,:)<=0.01) ) then
                            write(*,*) "WARNING: dz_interface very low (at level ",i,")"
                        endif

                    else ! above the flat_z_height
                        b1_i = 0
                        b2_i = 0
                        b1_mass = 0
                        b2_mass = 0

                        db1_i = 0
                        db2_i = 0
                        db1_mass = 0
                        db2_mass = 0

                        global_dz_interface(:,i,:) =  dz_scl(i)
                        global_z_interface(:,i,:)  = global_z_interface(:,i-1,:) + global_dz_interface(:,i-1,:)
                        !  if (i/=this%grid%kme)   z_interface(:,i+1,:) = z_interface(:,i,:) + dz(i) ! (dz(i) + dz_scl( i) )/2 !test in icar_s5T
                        z_u(:,i,:)  = z_u(:,i-1,:)  + (dz_scl(i) + dz_scl(i-1))*0.5 ! zr_u only relevant for first i above max level, aferwards both zr_u(i) AND zr_u(i-1) ar
                        z_v(:,i,:)  = z_v(:,i-1,:)  + (dz_scl(i) + dz_scl(i-1))*0.5
                        neighbor_z(:,i,:) =  neighbor_z(:,i-1,:) + (dz_scl(i-1) + dz_scl(i))*0.5

                    endif
                    dz_mass(:,i,:)   =  global_dz_interface(ims:ime,i-1,jms:jme) / 2  +  global_dz_interface(ims:ime,i,jms:jme) / 2
                endif ! if (i>kms)
                
                neighbor_jacobian(:,i,:) = 1 + h1*db1_mass + h2*db2_mass
                jacobian_u(:,i,:) = 1 + h1_u*db1_mass + h2_u*db2_mass
                jacobian_v(:,i,:) = 1 + h1_v*db1_mass + h2_v*db2_mass

                jacobian_w(:,i,:) = 1 + h1(ims:ime,jms:jme)*db1_i + h2(ims:ime,jms:jme)*db2_i
                
                dzdx(ims+1:ime-1,i,:) = (b1_mass*(h1(ims+2:ime,jms:jme)-h1(ims:ime-2,jms:jme)) + b2_mass*(h2(ims+2:ime,jms:jme)-h2(ims:ime-2,jms:jme)))/(2*this%dx)
                !dzdx(ims,i,:)   = (-neighbor_z(ims+2,i,:) + 4*neighbor_z(ims+1,i,:) - 3*neighbor_z(ims,i,:) )/(2*this%dx)
                !dzdx(ime,i,:)   = (neighbor_z(ime-2,i,:) - 4*neighbor_z(ime-1,i,:) + 3*neighbor_z(ime,i,:) )/(2*this%dx)
                dzdx(ims,i,:)   = dzdx(ims+1,i,:)
                dzdx(ime,i,:)   = dzdx(ime-1,i,:)

                dzdy(:,i,jms+1:jme-1) = (b1_mass*(h1(ims:ime,jms+2:jme)-h1(ims:ime,jms:jme-2)) + b2_mass*(h2(ims:ime,jms+2:jme)-h2(ims:ime,jms:jme-2)))/(2*this%dx)
                !dzdy(:,i,jms)   = (-neighbor_z(:,i,jms+2) + 4*neighbor_z(:,i,jms+1) - 3*neighbor_z(:,i,jms) )/(2*this%dx)
                !dzdy(:,i,jme)   = (neighbor_z(:,i,jme-2) - 4*neighbor_z(:,i,jme-1) + 3*neighbor_z(:,i,jme) )/(2*this%dx)
                dzdy(:,i,jms)   = dzdy(:,i,jms+1)
                dzdy(:,i,jme)   = dzdy(:,i,jme-1)
                !dzdx(:,i,:) = (b1_mass*(h1_u(ims+1:ime+1,jms:jme)-h1_u(ims:ime,jms:jme)) + b2_mass*(h2_u(ims+1:ime+1,jms:jme)-h2_u(ims:ime,jms:jme)))/(this%dx)
                !dzdy(:,i,:) = (b1_mass*(h1_v(ims:ime,jms+1:jme+1)-h1_v(ims:ime,jms:jme)) + b2_mass*(h2_v(ims:ime,jms+1:jme+1)-h2_v(ims:ime,jms:jme)))/(this%dx)
                
                dzdx_u(ims+1:ime,i,:) = (b1_mass*(h1(ims+1:ime,jms:jme)-h1(ims:ime-1,jms:jme))   + b2_mass*(h2(ims+1:ime,jms:jme)-h2(ims:ime-1,jms:jme)))/(this%dx)
                dzdy_v(:,i,jms+1:jme) = (b1_mass*(h1(ims:ime,jms+1:jme)-h1(ims:ime,jms:jme-1))   + b2_mass*(h2(ims:ime,jms+1:jme)-h2(ims:ime,jms:jme-1)))/(this%dx)

                dzdx_u(ims+1:ime-1,i,:) = (b1_mass*(h1_u(ims+2:ime,jms:jme)-h1_u(ims:ime-2,jms:jme))   + b2_mass*(h2_u(ims+2:ime,jms:jme)-h2_u(ims:ime-2,jms:jme)))/(2*this%dx) !(dzdx(ims+1:ime,i,:)+dzdx(ims:ime-1,i,:))*0.5
                dzdx_u(ims+1:ime,i,:) = (dzdx(ims:ime-1,i,:)+dzdx(ims+1:ime,i,:))*0.5
                dzdx_u(ims,i,:)   = dzdx(ims,i,:)*1.5 - dzdx(ims+1,i,:)*0.5
                !dzdx_u(ime,i,:)   = dzdx(ime,i,:)
                dzdx_u(ime+1,i,:)   = dzdx(ime,i,:)*1.5 - dzdx(ime-1,i,:)*0.5

                dzdy_v(:,i,jms+1:jme-1) = (b1_mass*(h1_v(ims:ime,jms+2:jme)-h1_v(ims:ime,jms:jme-2))   + b2_mass*(h2_v(ims:ime,jms+2:jme)-h2_v(ims:ime,jms:jme-2)))/(2*this%dx)!(dzdy(:,i,jms+1:jme)+dzdy(:,i,jms:jme-1))*0.5
                dzdy_v(:,i,jms+1:jme) = (dzdy(:,i,jms:jme-1)+dzdy(:,i,jms+1:jme))*0.5
                dzdy_v(:,i,jms)   = dzdy(:,i,jms)*1.5 - dzdy(:,i,jms+1)*0.5
                !dzdy_v(:,i,jme)   = dzdy(:,i,jme)
                dzdy_v(:,i,jme+1)   = dzdy(:,i,jme)*1.5 - dzdy(:,i,jme-1)*0.5
                if (options%physics%windtype == kLINEAR_ITERATIVE_WINDS .or. options%physics%windtype == kITERATIVE_WINDS) then
                    this%dzdxz(:,i,:) = (db1_mass*(h1_u(ims+1:ime+1,jms:jme)-h1_u(ims:ime,jms:jme)) + db2_mass*(h2_u(ims+1:ime+1,jms:jme)-h2_u(ims:ime,jms:jme)))/this%dx
                    this%dzdyz(:,i,:) = (db1_mass*(h1_v(ims:ime,jms+1:jme+1)-h1_v(ims:ime,jms:jme)) + db2_mass*(h2_v(ims:ime,jms+1:jme+1)-h2_v(ims:ime,jms:jme)))/this%dx
                    
                    this%dzdxz(ims+1:ime-1,i,:) = (db1_mass*(h1(ims+2:ime,jms:jme)-h1(ims:ime-2,jms:jme)) + db2_mass*(h2(ims+2:ime,jms:jme)-h2(ims:ime-2,jms:jme)))/(2*this%dx)
                    this%dzdxz(ims,i,:)   = this%dzdxz(ims+1,i,:)
                    this%dzdxz(ime,i,:)   = this%dzdxz(ime-1,i,:)

                    this%dzdyz(:,i,jms+1:jme-1) = (db1_mass*(h1(ims:ime,jms+2:jme)-h1(ims:ime,jms:jme-2)) + db2_mass*(h2(ims:ime,jms+2:jme)-h2(ims:ime,jms:jme-2)))/(2*this%dx)
                    this%dzdyz(:,i,jms)   = this%dzdyz(:,i,jms+1)
                    this%dzdyz(:,i,jme)   = this%dzdyz(:,i,jme-1)

                endif
            enddo
            
            !Finishing touch
            i=kme+1
            global_z_interface(:,i,:)  = global_z_interface(:,i-1,:) + global_dz_interface(:,i-1,:)
            
            ! this is on the subset grid:
            dz_interface = global_dz_interface(ims:ime,:,jms:jme)
            z_interface  = global_z_interface(ims:ime,:,jms:jme)
            z            = neighbor_z(ims:ime,:,jms:jme)
            jacobian     = neighbor_jacobian(ims:ime,:,jms:jme)
                        
            !call setup_dzdxy(this, options, neighbor_jacobian)


        end associate

    end subroutine setup_sleve



    !> -------------------------------
    !! Setup the vertical grid structure, in case SLEVE coordinates are not used.
    !!    This means either constant vertical height, or a simple terrain following coordinate (Gal-Chen)
    !!
    !! --------------------------------
    subroutine setup_simple_z(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options

        real, allocatable :: temp(:,:,:), temp_offset(:,:), global_jacobian(:,:,:)
        integer :: i, max_level

        associate(  ims => this%ims,      ime => this%ime,                        &
                    jms => this%jms,      jme => this%jme,                        &
                    kms => this%kms,      kme => this%kme,                        &
                    z                     => this%z%data_3d,                      &
                    z_u                   => this%geo_u%z,                        &
                    z_v                   => this%geo_v%z,                        &
                    z_interface           => this%z_interface%data_3d,            &
                    nz                    => options%parameters%nz,               &
                    dz                    => options%parameters%dz_levels,        &
                    dz_mass               => this%dz_mass%data_3d,                &
                    dz_interface          => this%dz_interface%data_3d,           &
                    terrain               => this%terrain%data_2d,                &
                    global_z_interface    => this%global_z_interface,             &
                    global_dz_interface   => this%global_dz_interface,            &
                    neighbor_terrain      => this%neighbor_terrain,               &
                    jacobian_u            => this%jacobian_u,                     &
                    jacobian_v            => this%jacobian_v,                     &
                    jacobian_w            => this%jacobian_w,                     &
                    dzdy                  => this%dzdy,                           &
                    jacobian              => this%jacobian,                       &
                    smooth_height         => this%smooth_height,                  &
                    dz_scl                => this%dz_scl)

            ! Start with a separate calculation for the lowest model level z=1
            i = this%grid%kms

            max_level = nz
            
            
            if ( (options%physics%windtype == kWIND_LINEAR) .or. (options%physics%windtype == kLINEAR_OBRIEN_WINDS) .or. &
                 (options%physics%windtype == kLINEAR_ITERATIVE_WINDS) ) then
                global_z_interface(:,i,:)   = this%global_terrain
                allocate(global_jacobian(this%ids:this%ide, this%kds:this%kde, this%jds:this%jde))
            else
                global_z_interface(:,i,:)   = neighbor_terrain
                allocate(global_jacobian(this%ihs:this%ihe, this%khs:this%khe, this%jhs:this%jhe))
            endif


            if (options%parameters%space_varying_dz) then
                max_level = find_flat_model_level(options, nz, dz)

                smooth_height = sum(dz(1:max_level))

                jacobian(:,i,:) = (smooth_height - terrain) / smooth_height ! sum(dz(1:max_level))
                global_jacobian(:,i,:) = (smooth_height - global_z_interface(:,i,:) ) /smooth_height !sum(dz(1:max_level))

            else
                jacobian = 1
                global_jacobian = 1
            endif

            dz_mass(:,i,:)      = dz(i) / 2 * jacobian(:,i,:)
            dz_interface(:,i,:) = dz(i) * jacobian(:,i,:)
            z(:,i,:)            = terrain + dz_mass(:,i,:)
            z_interface(:,i,:)  = terrain

            global_dz_interface(:,i,:) = dz(i) * global_jacobian(:,i,:)
            
            ! Now the higher (k!=1) levels can be calculated:
            do i = this%grid%kms+1, this%grid%kme
                if (i<=max_level) then
                    jacobian(:,i,:) = jacobian(:,i-1,:)
                    global_jacobian(:,i,:) = global_jacobian(:,i-1,:)

                else
                    jacobian(:,i,:) = 1
                    global_jacobian(:,i,:) = 1
                endif

                dz_mass(:,i,:)     = (dz(i)/2 * jacobian(:,i,:) + dz(i-1)/2 * jacobian(:,i-1,:))
                dz_interface(:,i,:)= dz(i) * jacobian(:,i,:)
                z(:,i,:)           = z(:,i-1,:)           + dz_mass(:,i,:)
                z_interface(:,i,:) = z_interface(:,i-1,:) + dz_interface(:,i-1,:)

                global_dz_interface(:,i,:) = dz(i) * global_jacobian(:,i,:)
                global_z_interface(:,i,:)  = global_z_interface(:,i-1,:) + global_dz_interface(:,i-1,:)

                jacobian(:,i,:) = dz_interface(:,i,:)/dz(i)
                global_jacobian(:,i,:) = global_dz_interface(:,i,:)/dz(i)

            enddo

            i = this%grid%kme + 1
            global_z_interface(:,i,:) = global_z_interface(:,i-1,:) + global_dz_interface(:,i-1,:)
            
            if (allocated(temp)) deallocate(temp)
            allocate(temp(this%ihs:this%ihe+1, this%khs:this%khe, this%jhs:this%jhe+1))

            temp(this%ihs,:,this%jhs:this%jhe) = global_jacobian(this%ihs,:,this%jhs:this%jhe)
            temp(this%ihe+1,:,this%jhs:this%jhe) = global_jacobian(this%ihe,:,this%jhs:this%jhe)
            temp(this%ihs+1:this%ihe,:,this%jhs:this%jhe) = (global_jacobian(this%ihs+1:this%ihe,:,this%jhs:this%jhe) + &
                                                                global_jacobian(this%ihs:this%ihe-1,:,this%jhs:this%jhe))/2
            jacobian_u = temp(ims:ime+1,:,jms:jme)

            temp(this%ihs:this%ihe,:,this%jhs) = global_jacobian(this%ihs:this%ihe,:,this%jhs)
            temp(this%ihs:this%ihe,:,this%jhe+1) = global_jacobian(this%ihs:this%ihe,:,this%jhe)
            temp(this%ihs:this%ihe,:,this%jhs+1:this%jhe) = (global_jacobian(this%ihs:this%ihe,:,this%jhs+1:this%jhe) + &
                                                global_jacobian(this%ihs:this%ihe,:,this%jhs:this%jhe-1))/2
            jacobian_v = temp(ims:ime,:,jms:jme+1)

            jacobian_w(:,this%kme,:) = jacobian(:,this%kme,:)
            jacobian_w(:,this%kms:this%kme-1,:) = (dz_interface(:,this%kms+1:this%kme,:)* jacobian(:,this%kms:this%kme-1,:) + &
                                                   dz_interface(:,this%kms:this%kme-1,:)* jacobian(:,this%kms+1:this%kme,:))/ &
                                                                                (dz_interface(:,this%kms:this%kme-1,:)+dz_interface(:,this%kms+1:this%kme,:))
                                                                                

            call array_offset_x(terrain, temp_offset)
            z_u(:,1,:) = temp_offset
            call array_offset_y(terrain, temp_offset)
            z_v(:,1,:) = temp_offset

            z_u(:,1,:)          = z_u(:,1,:) + dz(1) / 2 * jacobian_u(:,1,:)
            z_v(:,1,:)          = z_v(:,1,:) + dz(1) / 2 * jacobian_v(:,1,:)

            do i = this%grid%kms+1, this%grid%kme
                z_u(:,i,:) = z_u(:,i-1,:)  + ((dz(i)/2 * jacobian_u(:,i,:) + dz(i-1)/2 * jacobian_u(:,i-1,:)))
                z_v(:,i,:) = z_v(:,i-1,:)  + ((dz(i)/2 * jacobian_v(:,i,:) + dz(i-1)/2 * jacobian_v(:,i-1,:)))  
            enddo
                                                                                
            call setup_dzdxy(this, options, global_jacobian)
            
            jacobian = jacobian*smooth_height
            jacobian_u = jacobian_u*smooth_height
            jacobian_v = jacobian_v*smooth_height
            jacobian_w = jacobian_w*smooth_height
            where(jacobian==smooth_height) jacobian=1
            where(jacobian_u==smooth_height) jacobian_u=1
            where(jacobian_v==smooth_height) jacobian_v=1
            where(jacobian_w==smooth_height) jacobian_w=1

        end associate

    end subroutine setup_simple_z



    !> -------------------------------
    !! Initialize various domain variables, mostly z, dz, etc.
    !!
    !! -------------------------------
    subroutine initialize_core_variables(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options

        real, allocatable :: temp(:,:,:)
        integer :: i, j

        call read_core_variables(this, options)

        !setup geo_u/v here, because their z arrays will be calculated in the setup methods below
        call setup_geo(this%geo_u,   this%u_latitude%data_2d,   this%u_longitude%data_2d, options%parameters%longitude_system)
        call setup_geo(this%geo_v,   this%v_latitude%data_2d,   this%v_longitude%data_2d, options%parameters%longitude_system)
        allocate( this%geo_u%z(this%u_grid%ims:this%u_grid%ime, this%u_grid%nz, this%u_grid%jms:this%u_grid%jme))
        allocate( this%geo_v%z(this%v_grid%ims:this%v_grid%ime, this%v_grid%nz, this%v_grid%jms:this%v_grid%jme))

        call allocate_z_arrays(this, options)

        do i=this%grid%kms, this%grid%kme
            this%advection_dz(:,i,:) = options%parameters%dz_levels(i)
        enddo

        ! Setup the vertical grid structure, either as a SLEVE coordinate, or a more 'simple' vertical structure:
        if (options%parameters%sleve) then
            call setup_sleve(this, options)
        else
            ! This will set up either a Gal-Chen terrainfollowing coordinate, or no terrain following.
            call setup_simple_z(this, options)
        endif
        
        call setup_geo(this%geo,   this%latitude%data_2d,   this%longitude%data_2d, options%parameters%longitude_system,  this%z%data_3d)
        call setup_geo(this%geo_agl,   this%latitude%data_2d,   this%longitude%data_2d, options%parameters%longitude_system,  this%z%data_3d)

        call setup_grid_rotations(this, options)
        
        
        ! Setup variables applicable to near-surface wind modifications
        if (options%wind%Sx) then
            call setup_Sx(this, options)
        endif


    end subroutine initialize_core_variables
        
    subroutine setup_grid_rotations(this,options)
        type(domain_t),  intent(inout) :: this
        type(options_t), intent(in)    :: options

        integer :: i, j, ims, ime, jms, jme, smooth_loops
        integer :: starti, endi
        double precision :: dist, dlat, dlon

        real, allocatable :: lat(:,:), lon(:,:), costheta(:,:), sintheta(:,:)

        if (options%parameters%sinalpha_var /= "") then
            ims = this%ims
            ime = this%ime
            jms = this%jms
            jme = this%jme

            if (this_image()==1) print*, "Reading Sinalpha/cosalpha"

            call io_read(options%parameters%init_conditions_file, options%parameters%sinalpha_var, lon)
            this%sintheta = lon(ims:ime, jms:jme)

            call io_read(options%parameters%init_conditions_file, options%parameters%cosalpha_var, lon)
            this%costheta = lon(ims:ime, jms:jme)

        else

            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%lat_hi,                 &
                           lat)
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%lon_hi,                 &
                           lon)

            ims = this%ims
            ime = this%ime
            jms = this%jms
            jme = this%jme

            allocate(sintheta(ims:ime,jms:jme))
            allocate(costheta(ims:ime,jms:jme))

            do j = jms, jme
                do i = ims, ime
                    ! in case we are in the first or last grid, reset boundaries
                    starti = max(this%ids, i-2)
                    endi   = min(this%ide, i+2)

                    ! change in latitude
                    dlat = DBLE(lat(endi,j) - lat(starti,j))
                    ! change in longitude
                    dlon = DBLE(lon(endi,j) - lon(starti,j)) * cos(deg2rad*DBLE(lat(i,j)))
                    !if (abs(dlat) > 1) write(*,*) 'dlat:  ', dlat, '  ', ims, '  ', ime, '  ', jms, '  ', jme
                    !if (abs(dlon) > 1) write(*,*) 'dlon:  ', dlon, '  ', ims, '  ', ime, '  ', jms, '  ', jme
                    
                    ! distance between two points
                    dist = sqrt(DBLE(dlat)**2 + DBLE(dlon)**2) 

                    ! sin/cos of angles for use in rotating fields later
                    costheta(i, j) = abs(dlon / dist)
                    sintheta(i, j) =  (-1) * dlat / dist

                enddo
            enddo
            
            !Smooth cos/sin in case there are jumps from the lat/lon grids (more likely at low resolutions)
            smooth_loops = int(1000/this%dx)
            
            do i=1,smooth_loops
             call smooth_array_2d( costheta , windowsize  =  int((ime-ims)/5))
             call smooth_array_2d( sintheta , windowsize  =  int((ime-ims)/5))
            enddo
            this%costheta(ims:ime,jms:jme) = costheta(ims:ime,jms:jme)
            this%sintheta(ims:ime,jms:jme) = sintheta(ims:ime,jms:jme)
             
        endif
        if (options%parameters%debug .and.(this_image()==1)) then
            print*, ""
            print*, "Domain Geometry"
            print*, "MAX / MIN SIN(theta) (ideally 0)"
            print*, "   ", maxval(this%sintheta), minval(this%sintheta)
            print*, "MAX / MIN COS(theta) (ideally 1)"
            print*, "   ", maxval(this%costheta), minval(this%costheta)
            print*, ""
        endif


    end subroutine setup_grid_rotations

    subroutine setup_dzdxy(this, options, neighbor_jacobian)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options
        real, allocatable, intent(in)   :: neighbor_jacobian(:,:,:)
        
        real, allocatable :: neighbor_z(:,:,:)
        real, allocatable :: neighbor_dzdx(:,:,:)
        real, allocatable :: neighbor_dzdy(:,:,:)
        integer :: i

        allocate(neighbor_z( this% ihs : this% ihe, this% khs : this% khe, this% jhs : this% jhe) )
        allocate(neighbor_dzdx( this% ihs : this% ihe+1, this% khs : this% khe, this% jhs : this% jhe) )
        allocate(neighbor_dzdy( this% ihs : this% ihe, this% khs : this% khe, this% jhs : this% jhe+1) )

        neighbor_z(:,1,:) = this%neighbor_terrain + (options%parameters%dz_levels(1)/2)*neighbor_jacobian(:,1,:)

        do i=2,this%khe
            neighbor_z(:,i,:) = neighbor_z(:,i-1,:) + (((options%parameters%dz_levels(i)) / 2)*neighbor_jacobian(:,i,:)) + &
                                                  (((options%parameters%dz_levels(i-1)) / 2)*neighbor_jacobian(:,i-1,:))
        enddo

        neighbor_dzdx = 0
        neighbor_dzdy = 0

        !For dzdx
        ! neighbor_dzdx(this%ihs+1:this%ihe-1,:,:) = (neighbor_z(this%ihs+2:this%ihe,:,:) - &
        !                                                    neighbor_z(this%ihs:this%ihe-2,:,:))/(2*this%dx)
        !                                                                                                   
        ! neighbor_dzdx(this%ihs,:,:) = (-3*neighbor_z(this%ihs,:,:) + &
        !                                   4*neighbor_z(this%ihs+1,:,:) - neighbor_z(this%ihs+2,:,:)) / (2*this%dx)
        !                                   
        ! neighbor_dzdx(this%ihe,:,:) = (3*neighbor_z(this%ihe,:,:) - &
        !                                  4*neighbor_z(this%ihe-1,:,:) + neighbor_z(this%ihe-2,:,:)) / (2*this%dx)
        ! this%dzdx(:,:,:) = neighbor_dzdx(this%ims:this%ime,:,this%jms:this%jme)
        ! 
        
        
        ! neighbor_dzdx(this%ims+1:this%ime-1,:,this%jms:this%jme) = (this%z%data_3d(this%ims+2:this%ime,:,:) - &
        !                                                   this%z%data_3d(this%ims:this%ime-2,:,:))/(2*this%dx)
        !                                                                                                   
        ! neighbor_dzdx(this%ims,:,this%jms:this%jme) = (-3*this%z%data_3d(this%ims,:,:) + &
        !                                   4*this%z%data_3d(this%ims+1,:,:) - this%z%data_3d(this%ims+2,:,:)) / (2*this%dx)
        !                                   
        ! neighbor_dzdx(this%ime,:,this%jms:this%jme) = (3*this%z%data_3d(this%ime,:,:) - &
        !                                  4*this%z%data_3d(this%ime-1,:,:) + this%z%data_3d(this%ime-2,:,:)) / (2*this%dx)
        ! this%dzdx(:,:,:) = neighbor_dzdx(this%ims:this%ime,:,this%jms:this%jme)
        ! 
! 
        ! neighbor_dzdx(this%ims+1:this%ime,:,this%jms:this%jme) = (this%geo_u%z(this%ims+2:this%ime+1,:,:) - &
        !                                                    this%geo_u%z(this%ims:this%ime-1,:,:))/(2*this%dx)
        !                                                                                                   
        ! neighbor_dzdx(this%ims,:,this%jms:this%jme) = (-3*this%geo_u%z(this%ims,:,:) + &
        !                                   4*this%geo_u%z(this%ims+1,:,:) - this%geo_u%z(this%ims+2,:,:)) / (2*this%dx)
        !                                   
        ! neighbor_dzdx(this%ime+1,:,this%jms:this%jme) = (3*this%geo_u%z(this%ime+1,:,:) - &
        !                                  4*this%geo_u%z(this%ime,:,:) + this%geo_u%z(this%ime-1,:,:)) / (2*this%dx)
        ! this%dzdx_u(:,:,:) = neighbor_dzdx(this%ims:this%ime+1,:,this%jms:this%jme)
        
        
        this%dzdx(:,:,:) = (this%geo_u%z(this%ims+1:this%ime+1,:,:) - this%geo_u%z(this%ims:this%ime,:,:))/(this%dx)

        neighbor_dzdx(this%ihs+1:this%ihe,:,:) = (neighbor_z(this%ihs+1:this%ihe,:,:) - neighbor_z(this%ihs:this%ihe-1,:,:))/this%dx
        neighbor_dzdx(this%ihs,:,:) = neighbor_dzdx(this%ihs+1,:,:) 
        neighbor_dzdx(this%ihe+1,:,:) = neighbor_dzdx(this%ihe,:,:)
        this%dzdx_u(this%ims+1:this%ime,:,:) = neighbor_dzdx(this%ims+1:this%ime,:,this%jms:this%jme)
        
        this%dzdx_u(this%ims,:,:)   = this%dzdx(this%ims,:,:)*1.5 - this%dzdx(this%ims+1,:,:)*0.5
        this%dzdx_u(this%ime+1,:,:) = this%dzdx(this%ime,:,:)*1.5 - this%dzdx(this%ime-1,:,:)*0.5
        
        
        !For dzdy
        ! neighbor_dzdy(:,:,this%jhs+1:this%jhe-1) = (neighbor_z(:,:,this%jhs+2:this%jhe) - &
        !                                                    neighbor_z(:,:,this%jhs:this%jhe-2))/(2*this%dx)
        ! neighbor_dzdy(:,:,this%jhs) = (-3*neighbor_z(:,:,this%jms) + &
        !                                   4*neighbor_z(:,:,this%jms+1) - neighbor_z(:,:,this%jms+2)) / (2*this%dx)
        !                                   
        ! neighbor_dzdy(:,:,this%jhe) = (3*neighbor_z(:,:,this%jhe) - &
        !                                  4*neighbor_z(:,:,this%jhe-1) + neighbor_z(:,:,this%jhe-2)) / (2*this%dx)
        ! this%dzdy(:,:,:) = neighbor_dzdy(this%ims:this%ime,:,this%jms:this%jme)

        !  neighbor_dzdy(this%ims:this%ime,:,this%jms+1:this%jme-1) = (this%z%data_3d(:,:,this%jms+2:this%jme) - &
        !                                                     this%z%data_3d(:,:,this%jms:this%jme-2))/(2*this%dx)
        !  neighbor_dzdy(this%ims:this%ime,:,this%jms) = (-3*this%z%data_3d(:,:,this%jms) + &
        !                                    4*this%z%data_3d(:,:,this%jms+1) - this%z%data_3d(:,:,this%jms+2)) / (2*this%dx)
        !                                    
        !  neighbor_dzdy(this%ims:this%ime,:,this%jme) = (3*this%z%data_3d(:,:,this%jme) - &
        !                                   4*this%z%data_3d(:,:,this%jme-1) + this%z%data_3d(:,:,this%jme-2)) / (2*this%dx)
        !  this%dzdy(:,:,:) = neighbor_dzdy(this%ims:this%ime,:,this%jms:this%jme)
!  
!  
        !  neighbor_dzdy(this%ims:this%ime,:,this%jms+1:this%jme) = (this%geo_v%z(:,:,this%jms+2:this%jme+1) - &
        !                                                     this%geo_v%z(:,:,this%jms:this%jme-1))/(2*this%dx)
        !  neighbor_dzdy(this%ims:this%ime,:,this%jms) = (-3*this%geo_v%z(:,:,this%jms) + &
        !                                    4*this%geo_v%z(:,:,this%jms+1) - this%geo_v%z(:,:,this%jms+2)) / (2*this%dx)
        !                                    
        !  neighbor_dzdy(this%ims:this%ime,:,this%jme+1) = (3*this%geo_v%z(:,:,this%jme+1) - &
        !                                   4*this%geo_v%z(:,:,this%jme) + this%geo_v%z(:,:,this%jme-1)) / (2*this%dx)
        !  this%dzdy_v(:,:,:) = neighbor_dzdy(this%ims:this%ime,:,this%jms:this%jme+1)

        if (this%east_boundary .and. this%grid%yimg==8) then
            call io_write("dzdx_u_in_domain.nc", "dzdx_u", this%dzdx_u)
            call io_write("jacobian_u_in_domain.nc", "jacobian_u", this%jacobian_u)
            call io_write("geo_u_z_in_domain.nc", "geo_u_z", this%geo_u%z)
        endif

        this%dzdy(:,:,:) = (this%geo_v%z(:,:,this%jms+1:this%jme+1) - this%geo_v%z(:,:,this%jms:this%jme))/(this%dx)

        neighbor_dzdy(:,:,this%jhs+1:this%jhe) = (neighbor_z(:,:,this%jhs+1:this%jhe) - neighbor_z(:,:,this%jhs:this%jhe-1))/this%dx
        neighbor_dzdy(:,:,this%jhs) = neighbor_dzdy(:,:,this%jhs+1) 
        neighbor_dzdy(:,:,this%jhe+1) = neighbor_dzdy(:,:,this%jhe)
                
        this%dzdy_v(:,:,this%jms+1:this%jme) = neighbor_dzdy(this%ims:this%ime,:,this%jms+1:this%jme)
        this%dzdy_v(:,:,this%jms)   = this%dzdy(:,:,this%jms)*1.5 - this%dzdy(:,:,this%jms+1)*0.5
        this%dzdy_v(:,:,this%jme+1) = this%dzdy(:,:,this%jme)*1.5 - this%dzdy(:,:,this%jme-1)*0.5

    end subroutine setup_dzdxy


    !> -------------------------------
    !!  Separate the terrain into large scale and small scale terrain for SLEVE coordinate calculation
    !!  h(x,y) = h_1(x,y) + h_2(x,y) ;
    !!  where the subscripts 1 and 2 refer to large-scale and small-scale contributions, respectively.
    !!  The large-scale contribution h1 can be obtained from the full topography by an appropriate smoothing operation.
    !!
    !!  The smoothing is done over the entire (non-parallelized terrain, i.e. ids-ide). Afterwards the relevant variables
    !!  are subset to the respective paralellized grids. This is not the most efficient, but it makes the smoothing easier.
    !!
    !> -------------------------------

    subroutine split_topography(this, global_terr, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        real, dimension(this%ids:this%ide,this%jds:this%jde), intent(in)   :: global_terr
        type(options_t), intent(in)     :: options

        real, allocatable :: temp(:,:)
        integer :: i


        allocate(this%h1( this% ihs : this% ihe, &
                          this% jhs : this% jhe) )

        allocate(this%h2( this% ihs : this% ihe, &
                          this% jhs : this% jhe) )

        allocate(this%h1_u( this%u_grid2d% ims : this%u_grid2d% ime,   &
                            this%u_grid2d% jms : this%u_grid2d% jme) )

        allocate(this%h1_v( this%v_grid2d% ims : this%v_grid2d% ime,   &
                            this%v_grid2d% jms : this%v_grid2d% jme) )

        allocate(this%h2_u( this%u_grid2d% ims : this%u_grid2d% ime,   &
                            this%u_grid2d% jms : this%u_grid2d% jme) )

        allocate(this%h2_v( this%v_grid2d% ims : this%v_grid2d% ime,   &
                            this%v_grid2d% jms : this%v_grid2d% jme) )
                            
        associate(ims => this%ims,      ime => this%ime,                        &
                  jms => this%jms,      jme => this%jme,                        &
                  kms => this%kms,      kme => this%kme,                        &
                  h1                    => this%h1,                             &
                  h2                    => this%h2,                             &
                  h1_u                  => this%h1_u,                           &
                  h2_u                  => this%h2_u,                           &
                  h1_v                  => this%h1_v,                           &
                  h2_v                  => this%h2_v,                           &
                  terrain               => this%terrain%data_2d)


        if ((this_image()==1)) then
          write(*,*) "  Setting up the SLEVE vertical coordinate:"
          write(*,*) "    Smoothing large-scale terrain (h1) with a windowsize of ", &
                  options%parameters%terrain_smooth_windowsize, " for ",        &
                  options%parameters%terrain_smooth_cycles, " smoothing cylces."
        endif


        ! create a separate variable that will be smoothed later on:
        allocate(temp(this%ids:this%ide,this%jds:this%jde))
        temp =  global_terr


        h2 = temp(this%ihs:this%ihe,this%jhs:this%jhe)
        ! Smooth the terrain to attain the large-scale contribution h1 (_u/v):
        
        do i =1,options%parameters%terrain_smooth_cycles
          call smooth_array( temp, windowsize  =  options%parameters%terrain_smooth_windowsize)
        enddo
        
        h1   =  temp(this%ihs:this%ihe,this%jhs:this%jhe)
        h2   =  h2 - h1

        ! offset the global terrain for the h_(u/v) calculations:
        deallocate(temp)
        allocate(temp(this%ids:this%ide+1,this%jds:this%jde))
        call array_offset_x(global_terr, temp)
        temp(this%ids,this%jds:this%jde) = temp(this%ids+1,this%jds:this%jde)
        temp(this%ide+1,this%jds:this%jde) = temp(this%ide,this%jds:this%jde)
        
        h2_u = temp(this%u_grid2d%ims:this%u_grid2d%ime, this%u_grid2d%jms:this%u_grid2d%jme)
        do i =1,options%parameters%terrain_smooth_cycles
          call smooth_array( temp, windowsize = options%parameters%terrain_smooth_windowsize)
        enddo
        
        h1_u = temp(this%u_grid2d%ims:this%u_grid2d%ime, this%u_grid2d%jms:this%u_grid2d%jme)
        h2_u =  h2_u  - h1_u

        
        ! offset the global terrain for the h_(u/v) calculations:
        deallocate(temp)
        allocate(temp(this%ids:this%ide,this%jds:this%jde+1))
        call array_offset_y(global_terr, temp)
        temp(this%ids:this%ide,this%jds) = temp(this%ids:this%ide,this%jds+1)
        temp(this%ids:this%ide,this%jde+1) = temp(this%ids:this%ide,this%jde)
        h2_v = temp(this%v_grid2d%ims:this%v_grid2d%ime, this%v_grid2d%jms:this%v_grid2d%jme)
        
        do i =1,options%parameters%terrain_smooth_cycles
          call smooth_array( temp, windowsize = options%parameters%terrain_smooth_windowsize)
        enddo
        
        h1_v = temp(this%v_grid2d%ims:this%v_grid2d%ime, this%v_grid2d%jms:this%v_grid2d%jme)
        h2_v =  h2_v  - h1_v
        
        !if (this_image()==1) then
        !   write(*,*) "       Max of full topography", MAXVAL(neighbor_terrain )
        !   write(*,*) "       Max of large-scale topography (h1)  ", MAXVAL(h1)
        !   write(*,*) "       Max of small-scale topography (h2)  ", MAXVAL(h2)
        !end if

        end associate

    end subroutine split_topography


    subroutine setup_Sx(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options
        
        character(len=(len(trim(options%parameters%init_conditions_file))+3)) :: filename
        logical          :: fexist
        real, allocatable :: temporary_data(:,:,:,:)
        
        filename = trim(options%parameters%init_conditions_file)
        filename = filename(1:(len(filename)-6))//'_Sx.nc'
        !Check if we already have an Sx file
        inquire(file=filename,exist=fexist)
        
        !If we don't have an Sx file saved, build one
        !if (.not.(fexist)) then
        if (this_image()==1) write(*,*) "    Calculating Sx and TPI for wind modification"
        call calc_TPI(this, options)
        call calc_Sx(this, options, filename)
        !endif
    
        !Load Sx from file into domain_array
        !Assume 8 wind-directions for Sx
        !allocate(this%Sx( 1:72, this%grid2d% ims : this%grid2d% ime, 1:30,&
        !                       this%grid2d% jms : this%grid2d% jme) )        

        !call io_read(filename, 'Sx', temporary_data)

        !this%Sx = temporary_data(:,this%grid%ims:this%grid%ime, 1:30, this%grid%jms:this%grid%jme)
        
        !deallocate(temporary_data)
        
    end subroutine setup_Sx

    
    !>------------------------------------------------------------
    !! Calculate the ZNU and ZNW variables
    !!
    !! @param domain    Model domain structure
    !!
    !!------------------------------------------------------------
    subroutine init_znu(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        integer :: i, xpt, ypt
        real    :: ptop
        integer :: kms, kme

        kms = domain%kms
        kme = domain%kme

        ! one grid point into the domain gets a non-boundary point
        xpt = domain%ims + domain%grid%halo_size
        ypt = domain%jms + domain%grid%halo_size

        associate(p     => domain%pressure%data_3d,                         &
                  nz    => domain%nz,                                       &
                  psfc  => domain%surface_pressure%data_2d(xpt, ypt))

        ptop = p(xpt,kme,ypt) - (p(xpt,kme-1,ypt) - p(xpt,kme,ypt))/2.0 !NOT CORRECT
        ptop = max(ptop,1.0)

        if (allocated(domain%znu)) then
            do i=kms, kme
                domain%znu(i) = (p(xpt,i,ypt) - ptop) / (psfc - ptop)
            enddo
        endif

        if (allocated(domain%znw)) then
            do i = kms, kme
                if (i > kms) then
                    domain%znw(i) = ((p(xpt,i,ypt) + p(xpt,i-1,ypt)) / 2 - ptop) / (psfc-ptop)
                else
                    domain%znw(i) = 1
                endif
            enddo
        endif

        end associate

    end subroutine init_znu



    subroutine read_land_variables(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options

        integer :: i, nsoil
        real, allocatable :: temporary_data(:,:), temporary_data_3d(:,:,:)
        real :: soil_thickness(20)
        real :: init_surf_temp

        soil_thickness = 1.0
        soil_thickness(1:4) = [0.1, 0.2, 0.5, 1.0]
        init_surf_temp = 280

        if (this_image()==1) write (*,*) "Reading Land Variables"
        if (associated(this%soil_water_content%data_3d)) then
            nsoil = size(this%soil_water_content%data_3d, 2)
        elseif (associated(this%soil_temperature%data_3d)) then
            nsoil = size(this%soil_temperature%data_3d, 2)
        endif

        if (options%parameters%landvar /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%landvar,         &
                           temporary_data)
            if (allocated(this%land_mask)) then
                this%land_mask = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                where(this%land_mask==0) this%land_mask = kLC_WATER  ! To ensure conisitency. land_mask can be 0 or 2 for water, enforce a single value.
            endif
        endif

        if ((options%physics%watersurface==kWATER_LAKE) .AND.(options%parameters%lakedepthvar /= "")) then
            if (this_image()==1) write(*,*) "   reading lake depth data from hi-res file"

            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%lakedepthvar,         &
                           temporary_data)
            if (associated(this%lake_depth%data_2d)) then
                this%lake_depth%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif

        endif

        if (options%parameters%soiltype_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%soiltype_var,         &
                           temporary_data)
            if (allocated(this%soil_type)) then
                this%soil_type = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        endif

        if (options%parameters%soil_deept_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%soil_deept_var,       &
                           temporary_data)
            if (associated(this%soil_deep_temperature%data_2d)) then
                this%soil_deep_temperature%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)

                if (minval(temporary_data)< 200) then
                    if (this_image()==1) write(*,*) "WARNING, VERY COLD SOIL TEMPERATURES SPECIFIED:", minval(temporary_data)
                    if (this_image()==1) write(*,*) trim(options%parameters%init_conditions_file),"  ",trim(options%parameters%soil_deept_var)
                endif
                if (minval(this%soil_deep_temperature%data_2d)< 200) then
                    where(this%soil_deep_temperature%data_2d<200) this%soil_deep_temperature%data_2d=init_surf_temp ! <200 is just broken, set to mean annual air temperature at mid-latidudes
                endif
            endif
        else
            if (associated(this%soil_deep_temperature%data_2d)) then
                this%soil_deep_temperature%data_2d = init_surf_temp
            endif
        endif

        if (options%parameters%soil_t_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%soil_t_var,           &
                           temporary_data_3d)
            if (associated(this%soil_temperature%data_3d)) then
                do i=1,nsoil
                    this%soil_temperature%data_3d(:,i,:) = temporary_data_3d(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme, i)
                enddo
                if (options%parameters%soil_deept_var == "") then
                    if (associated(this%soil_deep_temperature%data_2d)) then
                        this%soil_deep_temperature%data_2d = this%soil_temperature%data_3d(:,nsoil,:)
                    endif
                endif
            endif

        else
            if (associated(this%soil_temperature%data_3d)) then
                if (associated(this%soil_deep_temperature%data_2d)) then
                    do i=1,nsoil
                        this%soil_temperature%data_3d(:,i,:) = this%soil_deep_temperature%data_2d
                    enddo
                endif
            endif
        endif


        if (options%parameters%swe_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%swe_var,         &
                           temporary_data)
            if (associated(this%snow_water_equivalent%data_2d)) then
                this%snow_water_equivalent%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif

        else
            if (associated(this%snow_water_equivalent%data_2d)) then
                this%snow_water_equivalent%data_2d = 0
            endif
        endif

        if (options%parameters%snowh_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%snowh_var,         &
                           temporary_data)
            if (associated(this%snow_height%data_2d)) then
                this%snow_height%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif

        else
            if (associated(this%snow_height%data_2d)) then
                this%snow_height%data_2d = 0
            endif
        endif

        if (options%parameters%soil_vwc_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%soil_vwc_var,         &
                           temporary_data_3d)
            if (associated(this%soil_water_content%data_3d)) then
                do i=1,nsoil
                    this%soil_water_content%data_3d(:,i,:) = temporary_data_3d(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme, i)
                enddo
            endif

        else
            if (associated(this%soil_water_content%data_3d)) then
                this%soil_water_content%data_3d = 0.4
            endif
        endif

        if (options%parameters%vegtype_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%vegtype_var,          &
                           temporary_data)
            if (allocated(this%veg_type)) then
                this%veg_type = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        endif

        if (options%parameters%albedo_var /= "") then
            if (options%lsm_options%monthly_albedo) then
                call io_read(options%parameters%init_conditions_file,   &
                            options%parameters%albedo_var,          &
                            temporary_data_3d)

                if (associated(this%albedo%data_3d)) then
                    do i=1,size(this%albedo%data_3d, 2)
                        this%albedo%data_3d(:,i,:) = temporary_data_3d(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme,i)
                    enddo

                    if (maxval(temporary_data_3d) > 1) then
                        if (this_image()==1) write(*,*) "Changing input ALBEDO % to fraction"
                        this%albedo%data_3d = this%albedo%data_3d / 100
                    endif
                endif
            else
                call io_read(options%parameters%init_conditions_file,   &
                               options%parameters%albedo_var,          &
                               temporary_data)
                if (associated(this%albedo%data_3d)) then
                    do i=1,size(this%albedo%data_3d, 2)
                        this%albedo%data_3d(:,i,:) = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                    enddo

                    if (maxval(temporary_data) > 1) then
                        if (this_image()==1) write(*,*) "Changing input ALBEDO % to fraction"
                        this%albedo%data_3d = this%albedo%data_3d / 100
                    endif
                endif
            endif

        else
            if (associated(this%albedo%data_3d)) then
                this%albedo%data_3d = 0.17
            endif
        endif


        if (options%parameters%vegfrac_var /= "") then
            if (options%lsm_options%monthly_albedo) then
                call io_read(options%parameters%init_conditions_file,   &
                            options%parameters%vegfrac_var,          &
                            temporary_data_3d)

                if (associated(this%vegetation_fraction%data_3d)) then
                    do i=1,size(this%vegetation_fraction%data_3d, 2)
                        this%vegetation_fraction%data_3d(:,i,:) = temporary_data_3d(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme,i)
                    enddo
                endif
            else
                call io_read(options%parameters%init_conditions_file,   &
                               options%parameters%vegfrac_var,          &
                               temporary_data)
                if (associated(this%vegetation_fraction%data_3d)) then
                    do i=1,size(this%vegetation_fraction%data_3d, 2)
                        this%vegetation_fraction%data_3d(:,i,:) = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                    enddo
                endif
            endif

        else
            if (associated(this%vegetation_fraction%data_3d)) then
                this%vegetation_fraction%data_3d = 60.
            endif
        endif

        if (associated(this%soil_totalmoisture%data_2d)) then
            this%soil_totalmoisture%data_2d = 0
            if (associated(this%soil_water_content%data_3d)) then
                do i=1, nsoil
                    this%soil_totalmoisture%data_2d = this%soil_totalmoisture%data_2d + this%soil_water_content%data_3d(:,i,:) * soil_thickness(i) * 1000 !! MJ added
                enddo
            endif
        endif

        if (options%parameters%vegfracmax_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%vegfracmax_var,       &
                           temporary_data)
            if (associated(this%vegetation_fraction_max%data_2d)) then
                this%vegetation_fraction_max%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (associated(this%vegetation_fraction_max%data_2d)) then
                if (this_image()==1) write(*,*) "    VEGMAX not specified; using default value of 0.8"
                this%vegetation_fraction_max%data_2d = 80.
            endif
        endif

        if (options%parameters%lai_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%lai_var,              &
                           temporary_data)
            if (associated(this%lai%data_2d)) then
                this%lai%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (associated(this%lai%data_2d)) then
                if (this_image()==1) write(*,*) "    LAI not specified; using default value of 1"
                this%lai%data_2d = 1
            endif
        endif

        if (options%parameters%canwat_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%canwat_var,              &
                           temporary_data)
            if (associated(this%canopy_water%data_2d)) then
                this%canopy_water%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (associated(this%canopy_water%data_2d)) then
                if (this_image()==1) write(*,*) "    CANWAT not specified; using default value of 0"
                this%canopy_water%data_2d = 0
            endif
        endif

        if (options%parameters%ridge_dist_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%ridge_dist_var,       &
                           temporary_data)
            if (associated(this%ridge_dist%data_2d)) then
                this%ridge_dist%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        endif
        
        if (options%parameters%valley_dist_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%valley_dist_var,       &
                           temporary_data)
            if (associated(this%valley_dist%data_2d)) then
                this%valley_dist%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        endif
        
        if (options%parameters%ridge_drop_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%ridge_drop_var,       &
                           temporary_data)
            if (associated(this%ridge_drop%data_2d)) then
                this%ridge_drop%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        endif
        
        if (options%parameters%shd_var /= "") then
            call io_read(options%parameters%init_conditions_file,   &
                           options%parameters%shd_var,       &
                           temporary_data)
            if (associated(this%shd%data_2d)) then
                this%shd%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        endif

        if (options%physics%radiation_downScaling==1) then
            !!
            if (options%parameters%hlm_var /= "") then
                call io_read(options%parameters%init_conditions_file,   &
                               options%parameters%hlm_var,           &
                               temporary_data_3d)
                if (associated(this%hlm%data_3d)) then
                    do i=1,90
                        this%hlm%data_3d(:,i,:) = temporary_data_3d(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme, i)
                        !if (this_image()==1) write(*,*),"hlm ", i, this%hlm%data_3d(this%grid%its,i,this%grid%jts)
                    enddo
                endif
            endif
            !!
            if (options%parameters%svf_var /= "") then
                call io_read(options%parameters%init_conditions_file,   &
                               options%parameters%svf_var,         &
                               temporary_data)
                if (associated(this%svf%data_2d)) then
                    this%svf%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                endif
            endif            
            !!
            if (options%parameters%slope_var /= "") then
                call io_read(options%parameters%init_conditions_file,   &
                               options%parameters%slope_var,         &
                               temporary_data)
                if (associated(this%slope%data_2d)) then
                    this%slope%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                endif
            endif            
            !!
            if (options%parameters%slope_angle_var /= "") then
                call io_read(options%parameters%init_conditions_file,   &
                               options%parameters%slope_angle_var,         &
                               temporary_data)
                if (associated(this%slope_angle%data_2d)) then
                    this%slope_angle%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                endif
            endif
            !!
            if (options%parameters%aspect_angle_var /= "") then
                call io_read(options%parameters%init_conditions_file,   &
                               options%parameters%aspect_angle_var,         &
                               temporary_data)
                if (associated(this%aspect_angle%data_2d)) then
                    this%aspect_angle%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                endif
            endif
            if (options%parameters%factor_p_var /= "") then
                 !if (this_image()==1) write(*,*) "facto_p is read...domain"
                call io_read(options%parameters%init_conditions_file,   &
                               options%parameters%factor_p_var,         &
                               temporary_data)
                if (associated(this%factor_p%data_2d)) then
                    this%factor_p%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                endif
            endif
        endif

        ! these will all be udpated by either forcing data or the land model, but initialize to sensible values to avoid breaking other initialization routines
        if (associated(this%skin_temperature%data_2d)) this%skin_temperature%data_2d = init_surf_temp
        if (associated(this%sst%data_2d)) this%sst%data_2d = init_surf_temp
        if (associated(this%roughness_z0%data_2d)) this%roughness_z0%data_2d = 0.001
        if (associated(this%sensible_heat%data_2d)) this%sensible_heat%data_2d=0
        if (associated(this%latent_heat%data_2d)) this%latent_heat%data_2d=0
        if (associated(this%u_10m%data_2d)) this%u_10m%data_2d=0
        if (associated(this%v_10m%data_2d)) this%v_10m%data_2d=0
        if (associated(this%windspd_10m%data_2d)) this%windspd_10m%data_2d=0
        if (associated(this%temperature_2m%data_2d)) this%temperature_2m%data_2d=init_surf_temp
        if (associated(this%humidity_2m%data_2d)) this%humidity_2m%data_2d=0.001
        if (associated(this%surface_pressure%data_2d)) this%surface_pressure%data_2d=102000
        if (associated(this%longwave_up%data_2d)) this%longwave_up%data_2d=0
        if (associated(this%ground_heat_flux%data_2d)) this%ground_heat_flux%data_2d=0
        if (associated(this%veg_leaf_temperature%data_2d)) this%veg_leaf_temperature%data_2d=init_surf_temp
        if (associated(this%ground_surf_temperature%data_2d)) this%ground_surf_temperature%data_2d=init_surf_temp
        if (associated(this%canopy_vapor_pressure%data_2d)) this%canopy_vapor_pressure%data_2d=2000
        if (associated(this%canopy_temperature%data_2d)) this%canopy_temperature%data_2d=init_surf_temp
        if (associated(this%coeff_momentum_drag%data_2d)) this%coeff_momentum_drag%data_2d=0.01
        if (associated(this%chs%data_2d)) this%chs%data_2d=0.01
        if (associated(this%chs2%data_2d)) this%chs2%data_2d=0.01
        if (associated(this%cqs2%data_2d)) this%cqs2%data_2d=0.01
        if (associated(this%qfx%data_2d)) this%qfx%data_2d=0.0
        if (associated(this%br%data_2d)) this%br%data_2d=0.0
        if (associated(this%psim%data_2d)) this%psim%data_2d=0.0
        if (associated(this%psih%data_2d)) this%psih%data_2d=0.0
        if (associated(this%fm%data_2d)) this%fm%data_2d=0.0
        if (associated(this%fh%data_2d)) this%fh%data_2d=0.0
        if (associated(this%hpbl%data_2d)) this%hpbl%data_2d=100.0
        if (associated(this%coeff_heat_exchange_3d%data_3d)) this%coeff_heat_exchange_3d%data_3d=0.01
        if (associated(this%coeff_momentum_exchange_3d%data_3d)) this%coeff_momentum_exchange_3d%data_3d=0.01
        if (associated(this%canopy_fwet%data_2d)) this%canopy_fwet%data_2d=0
        if (associated(this%snow_water_eq_prev%data_2d)) this%snow_water_eq_prev%data_2d=0
        if (associated(this%snow_albedo_prev%data_2d)) this%snow_albedo_prev%data_2d=0.65
        if (associated(this%storage_lake%data_2d)) this%storage_lake%data_2d=0
        
        if (associated(this%runoff_tstep%data_2d))        this%runoff_tstep%data_2d=0.
        if (associated(this%Tsnow%data_3d))               this%Tsnow%data_3d=273.15
        if (associated(this%Sice%data_3d))                this%Sice%data_3d=0.
        if (associated(this%Sliq%data_3d))                this%Sliq%data_3d=0.
        if (associated(this%Ds%data_3d))                  this%Ds%data_3d=0.
        if (associated(this%fsnow%data_2d))               this%fsnow%data_2d=0.
        if (associated(this%Nsnow%data_2d))               this%Nsnow%data_2d=0.
        if (associated(this%dm_salt%data_2d))             this%dm_salt%data_2d=0.
        if (associated(this%dm_susp%data_2d))             this%dm_susp%data_2d=0.
        if (associated(this%dm_subl%data_2d))             this%dm_subl%data_2d=0.
        if (associated(this%dm_slide%data_2d))            this%dm_slide%data_2d=0.

        !!
        if (associated(this%rainfall_tstep%data_2d))      this%rainfall_tstep%data_2d=0.
        if (associated(this%snowfall_tstep%data_2d))      this%snowfall_tstep%data_2d=0.
        if (associated(this%meltflux_out_tstep%data_2d))  this%meltflux_out_tstep%data_2d=0.
        if (associated(this%svf%data_2d))                 this%svf%data_2d=1.
        if (associated(this%hlm%data_2d))                 this%hlm%data_2d=90.
        if (associated(this%shortwave_direct%data_2d))  this%shortwave_direct%data_2d=0.
        if (associated(this%shortwave_diffuse%data_2d))  this%shortwave_diffuse%data_2d=0.
        if (associated(this%shortwave_direct_above%data_2d))  this%shortwave_direct_above%data_2d=0.
        if (associated(this%shortwave_total%data_2d))  this%shortwave_total%data_2d=0.
        if (associated(this%Sliq_out%data_2d))  this%Sliq_out%data_2d=0.
    end subroutine read_land_variables



    !> -------------------------------
    !! Populare the metadata structure in the domain for later output
    !!
    !! -------------------------------
    subroutine setup_meta_data(this, options)
        implicit none
        class(domain_t), intent(inout) :: this
        type(options_t), intent(in)    :: options
        character*60 :: a_string

        call this%info%add_attribute("comment",options%parameters%comment)
        call this%info%add_attribute("source","ICAR version:"//trim(options%parameters%version))

        ! Add info on grid setting:
        write(a_string,*) options%parameters%space_varying_dz
        call this%info%add_attribute("space_varying_dz",a_string)
        write(a_string,*) options%parameters%sleve
        call this%info%add_attribute("sleve",a_string)
        if (options%parameters%sleve) then
          write(a_string,*) options%parameters%terrain_smooth_windowsize
          call this%info%add_attribute("terrain_smooth_windowsize",a_string )
          write(a_string,*) options%parameters%terrain_smooth_cycles
          call this%info%add_attribute("terrain_smooth_cycles",a_string )
          write(a_string,*) options%parameters%decay_rate_L_topo
          call this%info%add_attribute("decay_rate_L_topo",a_string )
          write(a_string,*) options%parameters%decay_rate_s_topo
          call this%info%add_attribute("decay_rate_S_topo",a_string )
          write(a_string,*) options%parameters%sleve_n
          call this%info%add_attribute("sleve_n",a_string )
        endif
        ! Add some more info on physics settings:
        write(a_string,*) options%physics%boundarylayer
        call this%info%add_attribute("pbl", a_string )
        write(a_string,*) options%physics%landsurface
        call this%info%add_attribute("lsm", a_string )
        write(a_string,*) options%physics%watersurface
        call this%info%add_attribute("water", a_string )
        write(a_string,*) options%physics%microphysics
        call this%info%add_attribute("mp", a_string )
        write(a_string,*) options%physics%radiation
        call this%info%add_attribute("rad", a_string )
        write(a_string,*) options%physics%convection
        call this%info%add_attribute("conv", a_string )
        write(a_string,*) options%physics%advection
        call this%info%add_attribute("adv", a_string )
        write(a_string,*) options%physics%windtype
        call this%info%add_attribute("wind", a_string )


        call this%info%add_attribute("ids",str(this%grid%ids))
        call this%info%add_attribute("ide",str(this%grid%ide))
        call this%info%add_attribute("jds",str(this%grid%jds))
        call this%info%add_attribute("jde",str(this%grid%jde))
        call this%info%add_attribute("kds",str(this%grid%kds))
        call this%info%add_attribute("kde",str(this%grid%kde))

        !call this%info%add_attribute("ims",str(this%grid%ims))
        !call this%info%add_attribute("ime",str(this%grid%ime))
        !call this%info%add_attribute("jms",str(this%grid%jms))
        !call this%info%add_attribute("jme",str(this%grid%jme))
        !call this%info%add_attribute("kms",str(this%grid%kms))
        !call this%info%add_attribute("kme",str(this%grid%kme))

        !call this%info%add_attribute("its",str(this%grid%its))
        !call this%info%add_attribute("ite",str(this%grid%ite))
        !call this%info%add_attribute("jts",str(this%grid%jts))
        !call this%info%add_attribute("jte",str(this%grid%jte))
        !call this%info%add_attribute("kts",str(this%grid%kts))
        !call this%info%add_attribute("kte",str(this%grid%kte))

    end subroutine setup_meta_data


    !> -------------------------------
    !! Add variables needed by all domains to the list of requested variables
    !!
    !! -------------------------------
    module subroutine var_request(this, options)
        class(domain_t), intent(inout) :: this
        type(options_t), intent(inout) :: options

        integer :: i
        
        ! List the variables that are required to be allocated for any domain
        call options%alloc_vars(                                                    &
                     [kVARS%z,                      kVARS%z_interface,              &
                      kVARS%dz,                     kVARS%dz_interface,             &
                      kVARS%u,                      kVARS%v,                        &
                      kVARS%w,                      kVARS%w_real,                   &
                      kVARS%surface_pressure,       kVARS%roughness_z0,             &
                      kVARS%terrain,                kVARS%pressure,                 &
                      kVARS%temperature,            kVARS%pressure_interface,       &
                      kVARS%exner,                  kVARS%potential_temperature,    &
                      kVARS%latitude,               kVARS%longitude,                &
                      kVARS%u_latitude,             kVARS%u_longitude,              &
                      kVARS%v_latitude,             kVARS%v_longitude,              &
                      kVARS%temperature_interface,  kVars%density])

        if (trim(options%parameters%rain_var) /= "") call options%alloc_vars([kVARS%external_precipitation])

        ! List the variables that are required for any restart
        call options%restart_vars(                                                  &
                     [kVARS%z,                                                      &
                      kVARS%terrain,                kVARS%potential_temperature,    &
                      kVARS%latitude,               kVARS%longitude,                &
                      kVARS%u_latitude,             kVARS%u_longitude,              &
                      kVARS%v_latitude,             kVARS%v_longitude               ])

        !clean output var list
        do i=1, size(options%io_options%vars_for_output)
            if ((options%io_options%vars_for_output(i)+options%vars_for_restart(i) > 0) .and. (options%vars_to_allocate(i) <= 0)) then
                if (this_image()==1) write(*,*) 'var for kVARS index: ',options%vars_to_allocate(i),' requested for output/restart, but was not allocated by one of the modules'
                options%io_options%vars_for_output(i) = 0
            endif
        enddo

    end subroutine var_request

    !> -------------------------------
    !! Read in the shape of the domain required and setup the grid objects
    !!
    !! -------------------------------
    subroutine read_domain_shape(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options

        real, allocatable :: temporary_data(:,:)
        integer :: nx_global, ny_global, nz_global, nsmooth, n_neighbors, current, adv_order, my_index

        nsmooth = max(1, int(options%parameters%smooth_wind_distance / options%parameters%dx))
        if (options%parameters%smooth_wind_distance == 0.0) nsmooth = 0
        this%nsmooth = nsmooth
        if ((this_image()==1).and.(options%parameters%debug)) write(*,*) "number of gridcells to smooth = ",nsmooth
        ! This doesn't need to read in this variable, it could just request the dimensions
        ! but this is not a performance sensitive part of the code (for now)
        call io_read(options%parameters%init_conditions_file,   &
                     options%parameters%hgt_hi,                 &
                     temporary_data)

        nx_global = size(temporary_data,1)
        ny_global = size(temporary_data,2)
        nz_global = options%parameters%nz
        
        adv_order = max(options%adv_options%h_order,options%adv_options%v_order)
        
        !If we are using the monotonic flux limiter, it is necesarry to calculate the fluxes one location deep into the
        !halo. Thus, we need one extra cell in each halo direction to support the finite difference stencil
        !This is achieved here by artificially inflating the adv_order which is passed to the grid setup
        if (options%adv_options%flux_corr==kFLUXCOR_MONO) adv_order = adv_order+2
        
        !If using MPDATA, we need a halo of size 2 to support the difference stencil
        if (options%physics%advection==kADV_MPDATA) adv_order = 4

        call this%grid%set_grid_dimensions(     nx_global, ny_global, nz_global,adv_order=adv_order)
        call this%grid8w%set_grid_dimensions(   nx_global, ny_global, nz_global+1,adv_order=adv_order)

        call this%u_grid%set_grid_dimensions( nx_global, ny_global, nz_global,adv_order=adv_order, nx_extra = 1)
        call this%v_grid%set_grid_dimensions( nx_global, ny_global, nz_global,adv_order=adv_order, ny_extra = 1)

        ! for 2D mass variables
        call this%grid2d%set_grid_dimensions( nx_global, ny_global, 0,adv_order=adv_order)

        ! setup a 2D lat/lon grid extended by nsmooth grid cells so that smoothing can take place "across" images
        ! This just sets up the fields to interpolate u and v to so that the input data are handled on an extended
        ! grid.  They are then subset to the u_grid and v_grids above before actual use.
        call this%u_grid2d%set_grid_dimensions(     nx_global, ny_global, 0,adv_order=adv_order, nx_extra = 1)

        ! handle the v-grid too
        call this%v_grid2d%set_grid_dimensions(     nx_global, ny_global, 0,adv_order=adv_order, ny_extra = 1)
        
        call this%grid_soil%set_grid_dimensions(         nx_global, ny_global, kSOIL_GRID_Z,adv_order=adv_order)
        call this%grid_snow%set_grid_dimensions(         nx_global, ny_global, kSNOW_GRID_Z,adv_order=adv_order)
        call this%grid_snowsoil%set_grid_dimensions(     nx_global, ny_global, kSNOWSOIL_GRID_Z,adv_order=adv_order)
        call this%grid_soilcomp%set_grid_dimensions(     nx_global, ny_global, kSOILCOMP_GRID_Z,adv_order=adv_order)
        call this%grid_gecros%set_grid_dimensions(       nx_global, ny_global, kGECROS_GRID_Z,adv_order=adv_order)
        call this%grid_croptype%set_grid_dimensions(     nx_global, ny_global, kCROP_GRID_Z,adv_order=adv_order)
        call this%grid_monthly%set_grid_dimensions(      nx_global, ny_global, kMONTH_GRID_Z,adv_order=adv_order)
        call this%grid_lake%set_grid_dimensions(         nx_global, ny_global, kLAKE_Z,adv_order=adv_order)
        call this%grid_lake_soisno%set_grid_dimensions(  nx_global, ny_global, kLAKE_SOISNO_Z,adv_order=adv_order)
        call this%grid_lake_soi%set_grid_dimensions(     nx_global, ny_global, kLAKE_SOI_Z,adv_order=adv_order)
        call this%grid_lake_soisno_1%set_grid_dimensions(nx_global, ny_global, kLAKE_SOISNO_1_Z,adv_order=adv_order)
        call this%grid_hlm%set_grid_dimensions(     nx_global, ny_global, 90,adv_order=adv_order) !! MJ added


        this%north_boundary = (this%grid%yimg == this%grid%yimages)
        this%south_boundary = (this%grid%yimg == 1)
        this%east_boundary  = (this%grid%ximg == this%grid%ximages)
        this%west_boundary  = (this%grid%ximg == 1)

        this%ims = this%grid%ims; this%its = this%grid%its; this%ids = this%grid%ids
        this%ime = this%grid%ime; this%ite = this%grid%ite; this%ide = this%grid%ide
        this%kms = this%grid%kms; this%kts = this%grid%kts; this%kds = this%grid%kds
        this%kme = this%grid%kme; this%kte = this%grid%kte; this%kde = this%grid%kde
        this%jms = this%grid%jms; this%jts = this%grid%jts; this%jds = this%grid%jds
        this%jme = this%grid%jme; this%jte = this%grid%jte; this%jde = this%grid%jde
        
        !Calculate neighborhood indexes. These are used to store terrain in the local neighborhood for non-local wind calculations
        this%neighborhood_max = max(nsmooth,8)
        
        !Considering blocking terrain...
        if (options%physics%windtype == kLINEAR_ITERATIVE_WINDS .or. options%physics%windtype ==  kITERATIVE_WINDS) then
            this%neighborhood_max = int(max(4000.0/this%dx,1.0))
        endif
        
        !Considering TPI...
        if (options%wind%Sx) then
            this%neighborhood_max = max(this%neighborhood_max,floor(max(1.0,(options%wind%TPI_dmax+options%wind%Sx_dmax)/this%dx)))
        endif
        
        this%ihs=max(this%grid%ims-this%neighborhood_max,this%grid%ids); this%ihe=min(this%grid%ime+this%neighborhood_max,this%grid%ide)
        this%jhs=max(this%grid%jms-this%neighborhood_max,this%grid%jds); this%jhe=min(this%grid%jme+this%neighborhood_max,this%grid%jde)
        this%khs=this%grid%kms;                                          this%khe=this%grid%kme


        my_index = FINDLOC(DOM_IMG_INDX,this_image(),dim=1)
        !If we were found/are a compute process
        if (my_index > 0) then
            !Compute cardinal direction neighbors
            if (.not.(this%south_boundary)) this%south_neighbor = DOM_IMG_INDX(my_index - this%grid%ximages)
            if (.not.(this%north_boundary)) this%north_neighbor = DOM_IMG_INDX(my_index + this%grid%ximages)
            if (.not.(this%east_boundary)) this%east_neighbor  = DOM_IMG_INDX(my_index + 1)
            if (.not.(this%west_boundary)) this%west_neighbor  = DOM_IMG_INDX(my_index - 1)

            n_neighbors = merge(0,1,this%south_boundary)  &
                     +merge(0,1,this%north_boundary)  &
                     +merge(0,1,this%east_boundary)   &
                     +merge(0,1,this%west_boundary)
            n_neighbors = max(1, n_neighbors)

            allocate(this%neighbors(n_neighbors))

            current = 1
            if (.not. this%south_boundary) then
                this%neighbors(current) = this%south_neighbor
                current = current+1
            endif
            if (.not. this%north_boundary) then
                this%neighbors(current) = this%north_neighbor
                current = current+1
            endif
            if (.not. this%east_boundary) then
                this%neighbors(current) = this%east_neighbor
                current = current+1
            endif
            if (.not. this%west_boundary) then
                this%neighbors(current) = this%west_neighbor
                current = current+1
            endif
            
        ! if current = 1 then all of the boundaries were set, just store ourself as our "neighbor"
            if (current == 1) then
                this%neighbors(current) = this_image()
            endif

            !Compute diagonal direction neighbors
            if (.not.(this%south_boundary) .and. .not.(this%west_boundary)) this%southwest_neighbor = DOM_IMG_INDX(my_index - this%grid%ximages - 1)
            if (.not.(this%north_boundary) .and. .not.(this%west_boundary)) this%northwest_neighbor = DOM_IMG_INDX(my_index + this%grid%ximages - 1)
            if (.not.(this%south_boundary) .and. .not.(this%east_boundary)) this%southeast_neighbor  = DOM_IMG_INDX(my_index - this%grid%ximages + 1)
            if (.not.(this%north_boundary) .and. .not.(this%east_boundary)) this%northeast_neighbor  = DOM_IMG_INDX(my_index + this%grid%ximages + 1)

            n_neighbors = merge(0,1,(this%south_boundary .or. this%west_boundary))  &
                     +merge(0,1,(this%north_boundary .or. this%west_boundary))  &
                     +merge(0,1,(this%south_boundary .or. this%east_boundary))   &
                     +merge(0,1,(this%north_boundary .or. this%east_boundary))
            n_neighbors = max(1, n_neighbors)

            allocate(this%corner_neighbors(n_neighbors))

            current = 1
            if (.not.(this%south_boundary) .and. .not.(this%west_boundary)) then
                this%corner_neighbors(current) = this%southwest_neighbor
                current = current+1
            endif
            if (.not.(this%north_boundary) .and. .not.(this%west_boundary)) then
                this%corner_neighbors(current) = this%northwest_neighbor
                current = current+1
            endif
            if (.not.(this%south_boundary) .and. .not.(this%east_boundary)) then
                this%corner_neighbors(current) = this%southeast_neighbor
                current = current+1
            endif
            if (.not.(this%north_boundary) .and. .not.(this%east_boundary)) then
                this%corner_neighbors(current) = this%northeast_neighbor
                current = current+1
            endif

        ! if current = 1 then all of the boundaries were set, just store ourself as our "neighbor"
            if (current == 1) then
                this%corner_neighbors(current) = this_image()
            endif
        endif
    end subroutine

    !> -------------------------------
    !! Check that a set of variables is within realistic bounds (i.e. >0)
    !!
    !! Need to add more variables to the list
    !!
    !! -------------------------------
    module subroutine enforce_limits(this,update_in)
      class(domain_t), intent(inout) :: this
      logical, optional,  intent(in) :: update_in

      logical update
      update = .False.
      if (present(update_in)) update = update_in

      if (update) then
        if (associated(this%water_vapor%meta_data%dqdt_3d)      ) where(this%water_vapor%meta_data%dqdt_3d < 0)           this%water_vapor%meta_data%dqdt_3d = 0
        if (associated(this%potential_temperature%meta_data%dqdt_3d) ) where(this%potential_temperature%meta_data%dqdt_3d < 0) this%potential_temperature%meta_data%dqdt_3d = 0
        if (associated(this%cloud_water_mass%meta_data%dqdt_3d) ) where(this%cloud_water_mass%meta_data%dqdt_3d < 0)      this%cloud_water_mass%meta_data%dqdt_3d = 0
        if (associated(this%cloud_number%meta_data%dqdt_3d)     ) where(this%cloud_number%meta_data%dqdt_3d < 0)          this%cloud_number%meta_data%dqdt_3d = 0
        if (associated(this%cloud_ice_mass%meta_data%dqdt_3d)   ) where(this%cloud_ice_mass%meta_data%dqdt_3d < 0)        this%cloud_ice_mass%meta_data%dqdt_3d = 0
        if (associated(this%cloud_ice_number%meta_data%dqdt_3d) ) where(this%cloud_ice_number%meta_data%dqdt_3d < 0)      this%cloud_ice_number%meta_data%dqdt_3d = 0
        if (associated(this%rain_mass%meta_data%dqdt_3d)        ) where(this%rain_mass%meta_data%dqdt_3d < 0)             this%rain_mass%meta_data%dqdt_3d = 0
        if (associated(this%rain_number%meta_data%dqdt_3d)      ) where(this%rain_number%meta_data%dqdt_3d < 0)           this%rain_number%meta_data%dqdt_3d = 0
        if (associated(this%snow_mass%meta_data%dqdt_3d)        ) where(this%snow_mass%meta_data%dqdt_3d < 0)             this%snow_mass%meta_data%dqdt_3d = 0
        if (associated(this%snow_number%meta_data%dqdt_3d)      ) where(this%snow_number%meta_data%dqdt_3d < 0)           this%snow_number%meta_data%dqdt_3d = 0
        if (associated(this%graupel_mass%meta_data%dqdt_3d)     ) where(this%graupel_mass%meta_data%dqdt_3d < 0)          this%graupel_mass%meta_data%dqdt_3d = 0
        if (associated(this%graupel_number%meta_data%dqdt_3d)   ) where(this%graupel_number%meta_data%dqdt_3d < 0)        this%graupel_number%meta_data%dqdt_3d = 0
        if (associated(this%ice1_a%meta_data%dqdt_3d)           ) where(this%ice1_a%meta_data%dqdt_3d < 0)                this%ice1_a%meta_data%dqdt_3d = 0
        if (associated(this%ice1_c%meta_data%dqdt_3d)           ) where(this%ice1_c%meta_data%dqdt_3d < 0)                this%ice1_c%meta_data%dqdt_3d = 0
        if (associated(this%ice2_mass%meta_data%dqdt_3d)        ) where(this%ice2_mass%meta_data%dqdt_3d < 0)             this%ice2_mass%meta_data%dqdt_3d = 0
        if (associated(this%ice2_number%meta_data%dqdt_3d)      ) where(this%ice2_number%meta_data%dqdt_3d < 0)           this%ice2_number%meta_data%dqdt_3d = 0
        if (associated(this%ice2_a%meta_data%dqdt_3d)           ) where(this%ice2_a%meta_data%dqdt_3d < 0)                this%ice2_a%meta_data%dqdt_3d = 0
        if (associated(this%ice2_c%meta_data%dqdt_3d)           ) where(this%ice2_c%meta_data%dqdt_3d < 0)                this%ice2_c%meta_data%dqdt_3d = 0
        if (associated(this%ice3_mass%meta_data%dqdt_3d)        ) where(this%ice3_mass%meta_data%dqdt_3d < 0)             this%ice3_mass%meta_data%dqdt_3d = 0
        if (associated(this%ice3_number%meta_data%dqdt_3d)      ) where(this%ice3_number%meta_data%dqdt_3d < 0)           this%ice3_number%meta_data%dqdt_3d = 0
        if (associated(this%ice3_a%meta_data%dqdt_3d)           ) where(this%ice3_a%meta_data%dqdt_3d < 0)                this%ice3_a%meta_data%dqdt_3d = 0
        if (associated(this%ice3_c%meta_data%dqdt_3d)           ) where(this%ice3_c%meta_data%dqdt_3d < 0)                this%ice3_c%meta_data%dqdt_3d = 0

      else
        if (associated(this%water_vapor%data_3d)           ) where(this%water_vapor%data_3d < 0)             this%water_vapor%data_3d = 0
        if (associated(this%potential_temperature%data_3d) ) where(this%potential_temperature%data_3d < 0)   this%potential_temperature%data_3d = 0
        if (associated(this%cloud_water_mass%data_3d)      ) where(this%cloud_water_mass%data_3d < 0)        this%cloud_water_mass%data_3d = 0
        if (associated(this%cloud_number%data_3d)          ) where(this%cloud_number%data_3d < 0)            this%cloud_number%data_3d = 0
        if (associated(this%cloud_ice_mass%data_3d)        ) where(this%cloud_ice_mass%data_3d < 0)          this%cloud_ice_mass%data_3d = 0
        if (associated(this%cloud_ice_number%data_3d)      ) where(this%cloud_ice_number%data_3d < 0)        this%cloud_ice_number%data_3d = 0
        if (associated(this%rain_mass%data_3d)             ) where(this%rain_mass%data_3d < 0)               this%rain_mass%data_3d = 0
        if (associated(this%rain_number%data_3d)           ) where(this%rain_number%data_3d < 0)             this%rain_number%data_3d = 0
        if (associated(this%snow_mass%data_3d)             ) where(this%snow_mass%data_3d < 0)               this%snow_mass%data_3d = 0
        if (associated(this%snow_number%data_3d)           ) where(this%snow_number%data_3d < 0)             this%snow_number%data_3d = 0
        if (associated(this%graupel_mass%data_3d)          ) where(this%graupel_mass%data_3d < 0)            this%graupel_mass%data_3d = 0
        if (associated(this%graupel_number%data_3d)        ) where(this%graupel_number%data_3d < 0)          this%graupel_number%data_3d = 0
        if (associated(this%ice1_a%data_3d)                ) where(this%ice1_a%data_3d < 0)                  this%ice1_a%data_3d = 0
        if (associated(this%ice1_c%data_3d)                ) where(this%ice1_c%data_3d < 0)                  this%ice1_c%data_3d = 0
        if (associated(this%ice2_mass%data_3d)             ) where(this%ice2_mass%data_3d < 0)               this%ice2_mass%data_3d = 0
        if (associated(this%ice2_number%data_3d)           ) where(this%ice2_number%data_3d < 0)             this%ice2_number%data_3d = 0
        if (associated(this%ice2_a%data_3d)                ) where(this%ice2_a%data_3d < 0)                  this%ice2_a%data_3d = 0
        if (associated(this%ice2_c%data_3d)                ) where(this%ice2_c%data_3d < 0)                  this%ice2_c%data_3d = 0
        if (associated(this%ice3_mass%data_3d)             ) where(this%ice3_mass%data_3d < 0)               this%ice3_mass%data_3d = 0
        if (associated(this%ice3_number%data_3d)           ) where(this%ice3_number%data_3d < 0)             this%ice3_number%data_3d = 0
        if (associated(this%ice3_a%data_3d)                ) where(this%ice3_a%data_3d < 0)                  this%ice3_a%data_3d = 0
        if (associated(this%ice3_c%data_3d)                ) where(this%ice3_c%data_3d < 0)                  this%ice3_c%data_3d = 0

      endif
    end subroutine


    !> -------------------------------
    !! Setup the Geographic look up tables for interpolating a given forcing data set to each of the grids
    !!
    !! -------------------------------
    subroutine setup_geo_interpolation(this, forcing, options)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(inout) :: forcing
        type(options_t), intent(in)     :: options

        type(interpolable_type) :: forc_u_from_mass, forc_v_from_mass

        integer :: nx, ny, nz, i, j, k, ims, ime, jms, jme
        real, allocatable, dimension(:,:) :: AGL_cap, AGL_u_cap, AGL_v_cap, AGL_n, AGL_u_n, AGL_v_n

        ! this%geo and forcing%geo have to be of class interpolable
        ! which means they must contain lat, lon, z, geolut, and vLUT components

        call geo_LUT(this%geo_agl,forcing%geo_agl)
        call geo_LUT(this%geo,forcing%geo)
        call geo_LUT(this%geo_u,  forcing%geo_u)
        call geo_LUT(this%geo_v,  forcing%geo_v)

        if (allocated(forcing%z)) then  ! In case of external 2D forcing data, skip the VLUTs.

            forc_u_from_mass%lat = forcing%geo%lat
            forc_u_from_mass%lon = forcing%geo%lon
            forc_v_from_mass%lat = forcing%geo%lat
            forc_v_from_mass%lon = forcing%geo%lon

            call geo_LUT(this%geo_u, forc_u_from_mass)
            
            call geo_LUT(this%geo_v, forc_v_from_mass)
            
            nz = ubound(forcing%z,  2)
            ims = lbound(this%geo_u%z,1)
            ime = ubound(this%geo_u%z,1)
            jms = lbound(this%geo_u%z,3)
            jme = ubound(this%geo_u%z,3)
            allocate(forcing%geo_u%z(ims:ime, forcing%kts:forcing%kte, jms:jme))            

            ims = lbound(this%geo_v%z,1)
            ime = ubound(this%geo_v%z,1)
            jms = lbound(this%geo_v%z,3)
            jme = ubound(this%geo_v%z,3)            
            allocate(forcing%geo_v%z(ims:ime, forcing%kts:forcing%kte, jms:jme))
            

            ims = lbound(this%geo%z,1)
            ime = ubound(this%geo%z,1)
            jms = lbound(this%geo%z,3)
            jme = ubound(this%geo%z,3)
            allocate(forcing%geo%z(ims:ime, forcing%kts:forcing%kte, jms:jme))            
            allocate(forcing%geo_agl%z(ims:ime, forcing%kts:forcing%kte, jms:jme))            

            call geo_interp(forcing%geo%z, forcing%z, forcing%geo%geolut)
            call vLUT(this%geo,   forcing%geo)

            call geo_interp(forcing%geo_agl%z, forcing%z, forcing%geo%geolut)
            call geo_interp(forcing%geo_u%z, forcing%z, forc_u_from_mass%geolut)
            call geo_interp(forcing%geo_v%z, forcing%z, forc_v_from_mass%geolut)
            

            if (options%parameters%use_agl_height) then
                
                nx = size(this%geo_agl%z, 1)
                ny = size(this%geo_agl%z, 3)
                allocate(AGL_n(nx,ny))
                allocate(AGL_cap(nx,ny))

                nx = size(this%geo_u%z, 1)
                ny = size(this%geo_u%z, 3)
                allocate(AGL_u_n(nx,ny))
                allocate(AGL_u_cap(nx,ny))
                
                nx = size(this%geo_v%z, 1)
                ny = size(this%geo_v%z, 3)
                allocate(AGL_v_n(nx,ny))
                allocate(AGL_v_cap(nx,ny))


                AGL_cap = forcing%geo_agl%z(:,1,:)+real(options%parameters%agl_cap)
                where (AGL_cap <= (this%geo_agl%z(:,1,:)+200)) AGL_cap = this%geo_agl%z(:,1,:)+200
                
                AGL_u_cap = forcing%geo_u%z(:,1,:)+real(options%parameters%agl_cap)
                where (AGL_u_cap <= (this%geo_u%z(:,1,:)+200)) AGL_u_cap = this%geo_u%z(:,1,:)+200
                
                AGL_v_cap = forcing%geo_v%z(:,1,:)+real(options%parameters%agl_cap)
                where (AGL_v_cap <= (this%geo_v%z(:,1,:)+200)) AGL_v_cap = this%geo_v%z(:,1,:)+200

                !Do AGL interpolation for forcing geo z's
                do k=size(forcing%geo_agl%z, 2),1,-1
                    AGL_u_n = (AGL_u_cap-forcing%geo_u%z(:,k,:))/max(abs(AGL_u_cap-forcing%geo_u%z(:,1,:)),0.00001)
                    AGL_v_n = (AGL_v_cap-forcing%geo_v%z(:,k,:))/max(abs(AGL_v_cap-forcing%geo_v%z(:,1,:)),0.00001)
                    AGL_n = (AGL_cap-forcing%geo_agl%z(:,k,:))/max(abs(AGL_cap-forcing%geo_agl%z(:,1,:)),0.00001)

                    where (AGL_n < 0.0) AGL_n = 0.0
                    where (AGL_u_n < 0.0) AGL_u_n = 0.0
                    where (AGL_v_n < 0.0) AGL_v_n = 0.0
                    
                    forcing%geo_u%z(:,k,:) = forcing%geo_u%z(:,k,:)-forcing%geo_u%z(:,1,:)*AGL_u_n
                    forcing%geo_v%z(:,k,:) = forcing%geo_v%z(:,k,:)-forcing%geo_v%z(:,1,:)*AGL_v_n
                    forcing%geo_agl%z(:,k,:) = forcing%geo_agl%z(:,k,:)-forcing%geo_agl%z(:,1,:)*AGL_n
                enddo
                ! Step in reverse so that the bottom level is preserved until it is no longer needed
                ! Do AGL interpolation for domain grid
                
                do k=size(this%geo_agl%z,   2),1,-1
                    ! Multiply subtraction of base-topography by a factor that scales from 1 at surface to 0 at AGL_cap height
                    AGL_u_n = (AGL_u_cap-this%geo_u%z(:,k,:))/max(abs(AGL_u_cap-this%geo_u%z(:,1,:)),0.00001)
                    AGL_v_n = (AGL_v_cap-this%geo_v%z(:,k,:))/max(abs(AGL_v_cap-this%geo_v%z(:,1,:)),0.00001)
                    AGL_n = (AGL_cap-this%geo_agl%z(:,k,:))/max(abs(AGL_cap-this%geo_agl%z(:,1,:)),0.00001)

                    where (AGL_n < 0.0) AGL_n = 0.0
                    where (AGL_u_n < 0.0) AGL_u_n = 0.0
                    where (AGL_v_n < 0.0) AGL_v_n = 0.0
                    
                    this%geo_u%z(:,k,:) = this%geo_u%z(:,k,:)-this%geo_u%z(:,1,:)*AGL_u_n
                    this%geo_v%z(:,k,:) = this%geo_v%z(:,k,:)-this%geo_v%z(:,1,:)*AGL_v_n
                    this%geo_agl%z(:,k,:) = this%geo_agl%z(:,k,:)-this%geo_agl%z(:,1,:)*AGL_n
                enddo
            endif

            call vLUT(this%geo_agl,   forcing%geo_agl)
            call vLUT(this%geo_u, forcing%geo_u)
            call vLUT(this%geo_v, forcing%geo_v)
                        
            
            !if (options%parameters%use_agl_height) then
            !    do k=size(forcing%z,  2),1,-1
            !         forcing%z(:,k,:) = forcing%z(:,k,:)+forcing%original_geo%z(:,1,:)*AGL_forcing_n(:,k,:)
            !    enddo
            !endif

        end if
        

    end subroutine

    subroutine init_relax_filters(this)
        implicit none
        class(domain_t),    intent(inout) :: this
        integer :: hs, nr, k, i
        real, dimension(7) :: rs, rs_r
        logical :: corner
        !Setup relaxation filters, start with 2D then expand for 3D version
        
        allocate(this%relax_filter_2d(this%ims:this%ime,this%jms:this%jme))
        allocate(this%relax_filter_3d(this%ims:this%ime,this%kms:this%kme,this%jms:this%jme))
        
        associate( relax_filter => this%relax_filter_2d, relax_filter_3d => this%relax_filter_3d)

        corner = ((this%west_boundary .or. this%east_boundary) .and. (this%north_boundary .or. this%south_boundary))

        hs = this%grid%halo_size

        !relaxation boundary -- set to be 7 for default
        nr = min(7,(this%ime-this%ims-hs),(this%jme-this%jms-hs))
        
        rs = (/0.9, 0.75, 0.6, 0.5, 0.4, 0.25, 0.1 /)
        rs_r = (/0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9/)
        relax_filter = 0.0
        
        if (.not.(corner)) then
            if (this%west_boundary) then
                relax_filter(this%ims:this%ims+hs-1,this%jms:this%jme) = 1.0
                do k=this%jms,this%jme
                    relax_filter(this%ims+hs:this%ims+hs+nr-1,k) = rs(1:nr)
                enddo
            else if (this%east_boundary) then
                relax_filter(this%ime-hs+1:this%ime,this%jms:this%jme) = 1.0
                do k=this%jms,this%jme
                    relax_filter(this%ime-hs-nr+1:this%ime-hs,k) = rs_r(1:nr)        
                enddo
            else if (this%north_boundary) then
                relax_filter(this%ims:this%ime,this%jme-hs+1:this%jme) = 1.0
                do k=this%ims,this%ime
                    relax_filter(k,this%jme-hs-nr+1:this%jme-hs) = rs_r(1:nr)
                enddo
            else if (this%south_boundary) then
                relax_filter(this%ims:this%ime,this%jms:this%jms+hs-1) = 1.0
                do k=this%ims,this%ime
                    relax_filter(k,this%jms+hs:this%jms+hs+nr-1) = rs(1:nr)
                enddo
            endif
        else
            if (this%north_boundary .and. this%west_boundary) then
                relax_filter(this%ims:this%ims+hs-1,this%jms:this%jme) = 1.0
                relax_filter(this%ims:this%ime,this%jme-hs+1:this%jme) = 1.0

                do k=this%jms,this%jme-hs
                    relax_filter(this%ims+hs:this%ims+hs+nr-1,k) = rs(1:nr)
                enddo
                do k=this%ims+hs,this%ime
                    relax_filter(k,this%jme-hs-nr+1:this%jme-hs) = rs_r(1:nr)
                enddo
                do i = 1, nr
                    do k = 1, nr
                        relax_filter(this%ims+hs+i-1,this%jme-hs-k+1) = rs(min(i,k))
                    enddo
                enddo
            else if (this%north_boundary .and. this%east_boundary) then
                relax_filter(this%ime-hs+1:this%ime,this%jms:this%jme) = 1.0   
                relax_filter(this%ims:this%ime,this%jme-hs+1:this%jme) = 1.0

                do k=this%jms,this%jme-hs
                    relax_filter(this%ime-hs-nr+1:this%ime-hs,k) = rs_r(1:nr)        
                enddo        
                do k=this%ims,this%ime-hs
                    relax_filter(k,this%jme-hs-nr+1:this%jme-hs) = rs_r(1:nr)
                enddo
                do i = 1, nr
                    do k = 1, nr
                        relax_filter(this%ime-hs-i+1,this%jme-hs-k+1) = rs(min(i,k))
                    enddo
                enddo
            else if (this%south_boundary .and. this%west_boundary) then
                relax_filter(this%ims:this%ims+hs-1,this%jms:this%jme) = 1.0
                relax_filter(this%ims:this%ime,this%jms:this%jms+hs-1) = 1.0

                do k=this%jms+hs,this%jme
                    relax_filter(this%ims+hs:this%ims+hs+nr-1,k) = rs(1:nr)
                enddo
                do k=this%ims+hs,this%ime
                    relax_filter(k,this%jms+hs:this%jms+hs+nr-1) = rs(1:nr)
                enddo
                do i = 1, nr
                    do k = 1, nr
                        relax_filter(this%ims+hs+i-1,this%jms+hs+k-1) = rs(min(i,k))
                    enddo
                enddo
            else if (this%south_boundary .and. this%east_boundary) then
                relax_filter(this%ime-hs+1:this%ime,this%jms:this%jme) = 1.0   
                relax_filter(this%ims:this%ime,this%jms:this%jms+hs-1) = 1.0
            
                do k=this%jms+hs,this%jme
                    relax_filter(this%ime-hs-nr+1:this%ime-hs,k) = rs_r(1:nr)        
                enddo        
                do k=this%ims,this%ime-hs
                    relax_filter(k,this%jms+hs:this%jms+hs+nr-1) = rs(1:nr)
                enddo
                do i = 1, nr
                    do k = 1, nr
                        relax_filter(this%ime-hs-i+1,this%jms+hs+k-1) = rs(min(i,k))
                    enddo
                enddo
            endif
        endif

        do k=this%kms,this%kme
            relax_filter_3d(this%ims:this%ime,k,this%jms:this%jme) = relax_filter
        enddo
        
        end associate
    end subroutine init_relax_filters
    
    
    !> -------------------------------
    !! Update the dQdt fields for all forced variables which force the whole domain
    !! Forced variables which force just the boundary are handeled by a similar function called on the boundary object
    !! 
    !! For domain-forced variables, this routine is the partner of apply_forcing below.
    !! update_delta_fields normalizes the difference by the time step of that difference field
    !! apply_forcing multiplies that /second value and multiplies it by the current time step before adding it
    !!
    !! -------------------------------
    module subroutine update_delta_fields(this, dt)
        implicit none
        class(domain_t),    intent(inout) :: this
        type(time_delta_t), intent(in)    :: dt

        ! temporary to hold the variable to be interpolated to
        type(variable_t) :: var_to_update
                
        ! make sure the dictionary is reset to point to the first variable
        call this%variables_to_force%reset_iterator()

        ! Now iterate through the dictionary as long as there are more elements present
        do while (this%variables_to_force%has_more_elements())
            ! get the next variable
            var_to_update = this%variables_to_force%next()
            
            if (var_to_update%two_d) then
                if (.not.(var_to_update%force_boundaries)) var_to_update%dqdt_2d = (var_to_update%dqdt_2d - var_to_update%data_2d) / dt%seconds()
            else if (var_to_update%three_d) then
                if (.not.(var_to_update%force_boundaries)) var_to_update%dqdt_3d = (var_to_update%dqdt_3d - var_to_update%data_3d) / dt%seconds()
            endif

        enddo

        ! w has to be handled separately because it is the only variable that can be updated using the delta fields but is not
        ! actually read from disk. Note that if we move to balancing winds every timestep, then it doesn't matter.
        var_to_update = this%w%meta_data
        var_to_update%dqdt_3d = (var_to_update%dqdt_3d - var_to_update%data_3d) / dt%seconds()

    end subroutine


    !> -------------------------------
    !! Add the forcing update to boundaries and internal diagnosed fields
    !!
    !! This routine is the partner of update_delta_fields above.
    !! update_delta_fields normalizes the difference by the time step of that difference field
    !! apply forcing multiplies that /second value and multiplies it by the current time step before adding it
    !!
    !! -------------------------------
    module subroutine apply_forcing(this, forcing, options, dt)
        implicit none
        class(domain_t),    intent(inout) :: this
        class(boundary_t),  intent(inout) :: forcing
        type(options_t), intent(in)       :: options
        real, intent(in)                  :: dt
        integer :: ims, ime, jms, jme
        ! temporary to hold the variable to be interpolated to
        type(variable_t) :: var_to_update
        type(variable_t) :: forcing_hi
        integer :: i, k, j
        real    :: dt_h
        
        !calculate dt in units of hours
        dt_h = dt/3600.0
        
        ! make sure the dictionary is reset to point to the first variable
        call this%variables_to_force%reset_iterator()

        ! Now iterate through the dictionary as long as there are more elements present
        do while (this%variables_to_force%has_more_elements())
            ! get the next variable
            var_to_update = this%variables_to_force%next()
            
            forcing_hi = forcing%variables_hi%get_var(var_to_update%forcing_var)

            ims = var_to_update%grid%ims
            ime = var_to_update%grid%ime
            jms = var_to_update%grid%jms
            jme = var_to_update%grid%jme

            if (var_to_update%two_d) then
                ! apply forcing throughout the domain for 2D diagnosed variables (e.g. SST, SW)
                if (.not.(var_to_update%force_boundaries)) then
                    do j = jms,jme
                        do i = ims,ime
                            var_to_update%data_2d(i,j) = var_to_update%data_2d(i,j) + (var_to_update%dqdt_2d(i,j) * dt)
                        enddo
                    enddo
                else if (any(this%relax_filter_2d > 0.0)) then
                    !Update forcing data to current time step
                    do j = jms,jme
                        do i = ims,ime
                            if (this%relax_filter_2d(i,j) > 0.0) then
                                forcing_hi%data_2d(i,j) = forcing_hi%data_2d(i,j) + (forcing_hi%dqdt_2d(i,j) * dt)
                                if (this%relax_filter_2d(i,j) == 1.0) then
                                    var_to_update%data_2d(i,j) = forcing_hi%data_2d(i,j)
                                else
                                    var_to_update%data_2d(i,j) = var_to_update%data_2d(i,j) + &
                                                    (this%relax_filter_2d(i,j) * dt_h) * &
                                                    (forcing_hi%data_2d(i,j) - var_to_update%data_2d(i,j))
                                endif
                            endif
                        enddo
                    enddo
                endif 

            else if (var_to_update%three_d) then
                ! only apply forcing data on the boundaries for advected scalars (e.g. temperature, humidity)
                ! applying forcing to the edges has already been handeled when updating dqdt using the relaxation filter
                if (.not.(var_to_update%force_boundaries)) then
                    do j = jms,jme
                        do k = this%kms,this%kme
                            do i = ims,ime
                                forcing_hi%data_3d(i,k,j) = forcing_hi%data_3d(i,k,j) + (forcing_hi%dqdt_3d(i,k,j) * dt)
                                var_to_update%data_3d(i,k,j) = var_to_update%data_3d(i,k,j) + &
                                                              (var_to_update%dqdt_3d(i,k,j) * dt)
                            enddo
                        enddo
                    enddo
                else if (any(this%relax_filter_3d > 0.0)) then
                    !Update forcing data to current time step
                    do j = jms,jme
                        do k = this%kms,this%kme
                            do i = ims,ime
                                if (this%relax_filter_3d(i,k,j) > 0.0) then
                                    forcing_hi%data_3d(i,k,j) = forcing_hi%data_3d(i,k,j) + (forcing_hi%dqdt_3d(i,k,j) * dt)
                                    if (this%relax_filter_3d(i,k,j) == 1.0) then
                                        var_to_update%data_3d(i,k,j) = forcing_hi%data_3d(i,k,j)
                                    else
                                        var_to_update%data_3d(i,k,j) = var_to_update%data_3d(i,k,j) + &
                                                        (this%relax_filter_3d(i,k,j) * dt_h) * &
                                                        (forcing_hi%data_3d(i,k,j) - var_to_update%data_3d(i,k,j))
                                    endif
                                endif
                            enddo
                        enddo
                    enddo
                endif
            endif

        enddo

        ! w has to be handled separately because it is the only variable that can be updated using the delta fields but is not
        ! actually read from disk. Note that if we move to balancing winds every timestep, then it doesn't matter.
        if (.not.(options%parameters%advect_density)) then
            do j = jms,jme
                do k = this%kms,this%kme
                    do i = ims,ime
                        this%w%meta_data%data_3d(i,k,j) = this%w%meta_data%data_3d(i,k,j) + (this%w%meta_data%dqdt_3d(i,k,j) * dt)
                    enddo
                enddo
            enddo
        endif
        if (associated(this%external_precipitation%data_2d)) then
            if (associated(this%accumulated_precipitation%data_2d)) then
                this%accumulated_precipitation%data_2d = this%accumulated_precipitation%data_2d + (this%external_precipitation%data_2d * dt)
            endif
            if (associated(this%accumulated_precipitation%data_2dd)) then
                this%accumulated_precipitation%data_2dd = this%accumulated_precipitation%data_2dd + (this%external_precipitation%data_2d * dt)
            endif
        endif


    end subroutine

    !> -----------------------------------------------------------------------------------------------------------------
    !! Loop through external variables if supplied and interpolate the external data to the domain
    !!
    !! -----------------------------------------------------------------------------------------------------------------
    module subroutine interpolate_external(this, external_conditions, options)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(in)    :: external_conditions
        type(options_t), intent(in)     :: options

        integer :: i, nsoil=4
        character(len=99) :: varname
        type(variable_t) :: external_var, external_var2
        ! real, allocatable :: ext_snowheight_int(:,:)

        if(options%parameters%external_files/="MISSING") then
          ! -------  repeat this code block for other external variables?   -----------------
          if(options%parameters%swe_ext/="") then

            varname = options%parameters%swe_ext   !   options%ext_var_list(j)

            if (this_image()==1) write(*,*) "    interpolating external var ", trim(varname) , " for initial conditions"
            external_var =external_conditions%variables%get_var(trim(varname))  ! the external variable

            if (associated(this%snow_water_equivalent%data_2d)) then
                call geo_interp2d(  this%snow_water_equivalent%data_2d, & ! ( this%grid2d% ids : this%grid2d% ide, this%grid2d% jds : this%grid2d% jde)   ,            &
                                    external_var%data_2d,               &
                                    external_conditions%geo%geolut )
            endif

          endif

          ! -------  external snow height   -----------------
          if (options%parameters%hsnow_ext/="") then

            varname = options%parameters%hsnow_ext   !   options%ext_var_list(j)

            if (this_image()==1) write(*,*) "    interpolating external var ", trim(varname) , " for initial conditions"
            external_var =external_conditions%variables%get_var(trim(varname))  ! the external variable

            if (associated(this%snow_height%data_2d)) then
                call geo_interp2d(  this%snow_height%data_2d, & ! ( this%grid2d% ids : this%grid2d% ide, this%grid2d% jds : this%grid2d% jde)   ,            &
                                    external_var%data_2d,               &
                                    external_conditions%geo%geolut )
            endif
          ! -------  external snow height from external swe and density  -----------------
          elseif (options%parameters%swe_ext/="" .AND. options%parameters%rho_snow_ext/="") then

            varname = options%parameters%rho_snow_ext   !   options%ext_var_list(j)

            if (this_image()==1) write(*,*) "    interpolating external var ", trim(varname) , " to calculate initial snow height"
            external_var =external_conditions%variables%get_var(trim(varname))  ! the external variable
            external_var2 =external_conditions%variables%get_var(trim(options%parameters%swe_ext))  ! the external swe
            if (associated(this%snow_height%data_2d)) then
                call geo_interp2d(  this%snow_height%data_2d, &
                                    external_var2%data_2d / external_var%data_2d,               &  ! ext_swe / rho_snow_swe = hsnow_ext
                                    external_conditions%geo%geolut )
            endif
          endif

          ! ------ soil temperature  (2D or 3D)_______________________
          if (options%parameters%tsoil2D_ext/="") then

            varname = options%parameters%tsoil2D_ext   !   options%ext_var_list(j)

            if (this_image()==1) write(*,*) "    interpolating external var ", trim(varname) , " for initial conditions"
            external_var =external_conditions%variables%get_var(trim(varname))  ! the external variable

            if (associated(this%soil_deep_temperature%data_2d)) then

                call geo_interp2d(  this%soil_deep_temperature%data_2d, &
                                    external_var%data_2d,               &
                                    external_conditions%geo%geolut )
                if (associated(this%soil_temperature%data_3d)) then
                    do i=1,nsoil
                        this%soil_temperature%data_3d(:,i,:) = this%soil_deep_temperature%data_2d
                    enddo
                endif
            endif

          elseif (options%parameters%tsoil3D_ext/="") then  ! if 3D soil is provided we take the lowest level only. (can/should be expanded later)

            varname = options%parameters%tsoil3D_ext

            if (this_image()==1) write(*,*) "    interpolating external var ", trim(varname) , " for initial conditions"
            external_var =external_conditions%variables%get_var(trim(varname))  ! the external variable

            if (associated(this%soil_deep_temperature%data_2d)) then
                call geo_interp2d(  this%soil_deep_temperature%data_2d, &
                                    external_var%data_3d(:,size(external_var%data_3d,2),:)  ,               &
                                    external_conditions%geo%geolut )
                if (associated(this%soil_temperature%data_3d)) then
                    do i=1,nsoil
                        this%soil_temperature%data_3d(:,i,:) = this%soil_deep_temperature%data_2d
                    enddo
                endif
            endif
          endif
        endif
    end subroutine


    !> -------------------------------
    !! Loop through all variables for which forcing data have been supplied and interpolate the forcing data to the domain
    !!
    !! -------------------------------
    module subroutine interpolate_forcing(this, forcing, update)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(in)    :: forcing
        logical,          intent(in),   optional :: update

        ! internal field always present for value of optional "update"
        logical :: update_only
        ! temporary to hold the variable to be interpolated to
        type(variable_t) :: var_to_interpolate
        ! temporary to hold pressure and temperature for later below-grid adjustments
        type(variable_t) :: pressure, potential_temp
        ! temporary to hold the forcing variable to be interpolated from
        type(variable_t) :: input_data
        type(variable_t) :: forcing_hi

        ! number of layers has to be used when subsetting for update_pressure (for now)
        integer :: nz
        logical :: var_is_u, var_is_v, var_is_pressure, var_is_potential_temp, agl_interp

        update_only = .False.
        if (present(update)) update_only = update

        ! make sure the dictionary is reset to point to the first variable
        call this%variables_to_force%reset_iterator()

        ! Now iterate through the dictionary as long as there are more elements present
        do while (this%variables_to_force%has_more_elements())
            ! get the next variable
            var_to_interpolate = this%variables_to_force%next()
            
            forcing_hi = forcing%variables_hi%get_var(var_to_interpolate%forcing_var)

            ! get the associated forcing data
            input_data = forcing%variables%get_var(var_to_interpolate%forcing_var)

            ! interpolate
            if (var_to_interpolate%two_d) then
                if (update_only) then
                    call geo_interp2d(forcing_hi%dqdt_2d, input_data%data_2d, forcing%geo%geolut)
                    !If this variable is forcing the whole domain, we can copy the next forcing step directly over to domain
                    if (.not.(var_to_interpolate%force_boundaries)) var_to_interpolate%dqdt_2d = forcing_hi%dqdt_2d
                else
                    call geo_interp2d(forcing_hi%data_2d, input_data%data_2d, forcing%geo%geolut)
                    !If this is an initialization step, copy high res directly over to domain
                    var_to_interpolate%data_2d = forcing_hi%data_2d
                endif

            else
                var_is_pressure = (trim(var_to_interpolate%forcing_var) == trim(this%pressure%forcing_var))
                var_is_potential_temp = (trim(var_to_interpolate%forcing_var) == trim(this%potential_temperature%meta_data%forcing_var))
                var_is_u = (trim(var_to_interpolate%forcing_var) == trim(this%u%meta_data%forcing_var))
                var_is_v = (trim(var_to_interpolate%forcing_var) == trim(this%v%meta_data%forcing_var))
                
                !If we are dealing with anything but pressure and temperature (basically mass/number species), consider height above ground
                !for interpolation. If the user has not selected AGL interpolation in the namelist, this will result in standard z-interpolation
                agl_interp = .not.(var_is_pressure .or. var_is_potential_temp)

                ! if just updating, use the dqdt variable otherwise use the 3D variable
                if (update_only) then
                    call interpolate_variable(forcing_hi%dqdt_3d, input_data, forcing, this, &
                                    interpolate_agl_in=agl_interp, var_is_u=var_is_u, var_is_v=var_is_v, nsmooth=this%nsmooth)
                    !If this variable is forcing the whole domain, we can copy the next forcing step directly over to domain
                    if (.not.(var_to_interpolate%force_boundaries).and..not.var_is_u.and..not.var_is_v) var_to_interpolate%dqdt_3d = forcing_hi%dqdt_3d
                else
                    call interpolate_variable(forcing_hi%data_3d, input_data, forcing, this, &
                                    interpolate_agl_in=agl_interp, var_is_u=var_is_u, var_is_v=var_is_v, nsmooth=this%nsmooth)
                    !If this is an initialization step, copy high res directly over to domain
                    var_to_interpolate%data_3d = forcing_hi%data_3d
                endif
                
                if (var_is_pressure) pressure = forcing_hi
                if (var_is_potential_temp) potential_temp = forcing_hi
            endif
            
        enddo

        !Adjust potential temperature (first) and pressure (second) to account for points below forcing grid
        !Only domain-wide-forced variables are updated with the domain dqdt_3d
        if (update_only) then
            call adjust_pressure_temp(pressure%dqdt_3d,potential_temp%dqdt_3d, forcing%geo%z, this%geo%z)
            this%pressure%dqdt_3d = pressure%dqdt_3d
            this%potential_temperature%meta_data%dqdt_3d = potential_temp%dqdt_3d
            !this%w_real%dqdt_3d = 0
        else
            call adjust_pressure_temp(pressure%data_3d,potential_temp%data_3d, forcing%geo%z, this%geo%z)
            this%pressure%data_3d = pressure%data_3d
            this%potential_temperature%data_3d = potential_temp%data_3d
            !this%w_real%data_3d = 0
        endif
        
        !Ensure that input data for hydrometeors after interpolation have been forced to 0-minimum
        call this%enforce_limits(update_in=update_only)

    end subroutine

    subroutine adjust_pressure_temp(pressure, potential_temp, input_z, output_z)
        implicit none
        real, intent(inout), dimension(:,:,:) :: pressure, potential_temp
        real, intent(in), dimension(:,:,:) :: input_z, output_z !> z on the forcing and ICAR model levels [m]
        integer :: i,j,k, nz, nx, ny
        real    :: t, p_guess, dz
        
        !For all output_z less than input_z, extrapolate downwards based on lapse rate of -6.5C/km
        
        nx = size(potential_temp, 1)
        nz = size(potential_temp, 2)
        ny = size(potential_temp, 3)

        do j = 1, ny
            do i = 1, nx
                do k = 1, nz
                    if (input_z(i,1,j) > output_z(i,k,j)) then
                        
                        !From vertical interpolation, potential_temperature and pressure will be kept constant when below the grid
                        !So the current values at these below-indices reflect the temp/pressure of the closest forcing grid cell
                    
                        dz = input_z(i,1,j)-output_z(i,k,j)
                        
                        !Assume lapse rate of -6.5ºC/1km
                        !potential_temp(i,k,j) = potential_temp(i,k,j) + 6.5*dz/1000.0
                        
                        !estimate pressure difference 1100 Pa for each 100m difference for exner function
                        p_guess = pressure(i,k,j) + 1100*dz/100.0
                        t = exner_function(p_guess) * potential_temp(i,k,j)
                        pressure(i,k,j) = pressure(i,k,j) * exp( ((gravity/R_d) * dz) / t )
                    else
                        exit
                    endif
                end do
            enddo
        enddo


    end subroutine adjust_pressure_temp

    !> -------------------------------
    !! Adjust a 3d pressure field from the forcing data to the ICAR model grid
    !!
    !! Because the GCM grid can be very different from the ICAR grid, we first roughly match up
    !! the GCM level that is closest to the ICAR level. This has to be done grid cell by gridcell.
    !! This still is not ideal, in that it has already subset the GCM levels to the same number as are in ICAR
    !! If the GCM has a LOT of fine layers ICAR will not be getting layers higher up in the atmosphere.
    !! It would be nice to first use vinterp to get as close as we can, then update pressure only for grid cells below.
    !! Uses update_pressure to make a final adjustment (including below the lowest model level).
    !!
    !! -------------------------------
    subroutine adjust_pressure(pressure, input_z, output_z, potential_temperature)
        implicit none
        real, intent(inout), dimension(:,:,:) :: pressure !> Pressure on the forcing model levels [Pa]
        real, intent(in), dimension(:,:,:) :: input_z, output_z !> z on the forcing and ICAR model levels [m]
        real, intent(in), dimension(:,:,:) :: potential_temperature !> potential temperature of the forcing data [K]

        ! store a temporary copy of P and Z from the forcing data after selecting the closest GCM level to the ICAR data
        real, allocatable, dimension(:,:,:) :: temp_z, temp_p, temp_t
        ! loop counter variables
        integer :: k, nz, in_z_idx
        integer :: i,j, nx, ny

        allocate(temp_z, temp_p, temp_t, mold=pressure)

        nx = size(pressure, 1)
        nz = size(pressure, 2)
        ny = size(pressure, 3)

        do j = 1, ny
            do i = 1, nx
                ! keep track of the nearest z level from the forcing data
                in_z_idx = 1
                do k = 1, nz
                    ! if the ICAR z level is more than half way to the next forcing z level, then increment the GCM z
                    findz: do while (output_z(i,k,j) > ((input_z(i,in_z_idx,j) + input_z(i,min(nz,in_z_idx+1),j)) / 2))
                        in_z_idx = min(nz, in_z_idx + 1)

                        if (in_z_idx == nz) then
                            exit findz
                        endif
                    end do findz
                    ! make a new copy of the pressure and z data from the closest GCM model level
                    temp_z(i,k,j) = input_z(i,in_z_idx,j)
                    temp_p(i,k,j) = pressure(i,in_z_idx,j)
                    temp_t(i,k,j) = exner_function(pressure(i,in_z_idx,j)) * potential_temperature(i,in_z_idx,j)
                end do
            enddo
        enddo

        ! put the updated pressure data into the pressure variable prior to adjustments
        pressure = temp_p


        ! update pressure for the change in height between the closest GCM model level and each ICAR level.
        call update_pressure(pressure, temp_z, output_z, temp_t)

    end subroutine

    !> -------------------------------
    !! Interpolate one variable by requesting the forcing data from the boundary data structure then
    !! calling the appropriate interpolation routine (2D vs 3D) with the appropriate grid (mass, u, v)
    !!
    !! -------------------------------
    subroutine interpolate_variable(var_data, input_data, forcing, dom, interpolate_agl_in, var_is_u, var_is_v, nsmooth)
        implicit none
        real,            intent(inout) :: var_data(:,:,:)
        type(variable_t),   intent(in) :: input_data
        type(boundary_t),   intent(in)    :: forcing
        type(domain_t),     intent(in)    :: dom
        logical,            intent(in),   optional :: interpolate_agl_in
        logical,            intent(in),   optional :: var_is_u, var_is_v
        integer,            intent(in),   optional :: nsmooth

        ! note that 3D variables have a different number of vertical levels, so they have to first be interpolated
        ! to the high res horizontal grid, then vertically interpolated to the actual icar domain
        real, allocatable :: temp_3d(:,:,:)
        logical :: interpolate_agl, uvar, vvar
        integer :: nx, ny, nz, ims, ime, jms, jme
        integer :: windowsize, z

        interpolate_agl=.False.
        if (present(interpolate_agl_in)) interpolate_agl = interpolate_agl_in
        uvar = .False.
        if (present(var_is_u)) uvar = var_is_u
        vvar = .False.
        if (present(var_is_v)) vvar = var_is_v
        windowsize = 0
        if (present(nsmooth)) windowsize = nsmooth

        ims = lbound(var_data,1)
        ime = ubound(var_data,1)
        jms = lbound(var_data,3)
        jme = ubound(var_data,3)

        ! allocate a temporary variable to hold the horizontally interpolated data before vertical interpolation
        allocate(temp_3d(ims:ime, size(input_data%data_3d,2), jms:jme ))

        ! Sequence of if statements to test if this variable needs to be interpolated onto the staggared grids
        ! This could all be combined by passing in the geo data to use, along with a smoothing flag.

        ! Interpolate to the Mass grid
        if ((size(var_data,1) == size(forcing%geo%geolut%x,2)).and.(size(var_data,3) == size(forcing%geo%geolut%x,3))) then

            call geo_interp(temp_3d, input_data%data_3d, forcing%geo%geolut)

            if (interpolate_agl) then
                call vinterp(var_data, temp_3d, forcing%geo_agl%vert_lut)
            else
                call vinterp(var_data, temp_3d, forcing%geo%vert_lut)
            endif
            
        ! Interpolate to the u staggered grid
        else if (uvar) then

            ! One grid cell smoothing of original input data
            if (windowsize > 0) call smooth_array(input_data%data_3d, windowsize=1, ydim=3)
            call geo_interp(temp_3d, input_data%data_3d, forcing%geo_u%geolut)

            call vinterp(var_data, temp_3d, forcing%geo_u%vert_lut)
            ! temp_3d = pre_smooth(:,:nz,:) ! no vertical interpolation option
            if (windowsize > 0) call smooth_array(var_data, windowsize=windowsize, ydim=3)
                        
        ! Interpolate to the v staggered grid
        else if (vvar) then

            ! One grid cell smoothing of original input data
            if (windowsize > 0) call smooth_array(input_data%data_3d, windowsize=1, ydim=3)
            call geo_interp(temp_3d, input_data%data_3d, forcing%geo_v%geolut)
            
            call vinterp(var_data, temp_3d, forcing%geo_v%vert_lut)
            ! temp_3d = pre_smooth(:,:nz,:) ! no vertical interpolation option
            if (windowsize > 0) call smooth_array(var_data, windowsize=windowsize, ydim=3)
        endif
        
    end subroutine



end submodule
