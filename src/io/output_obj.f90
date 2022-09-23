submodule(output_interface) output_implementation
  use output_metadata,          only : get_metadata
  use debug_module,             only : check_ncdf
  implicit none

contains

    module subroutine set_domain(this, domain)
        class(output_t),  intent(inout)  :: this
        type(domain_t),   intent(in)     :: domain
        integer :: i

        if (.not.this%is_initialized) call this%init()

        do i=1,domain%info%n_attrs
            call this%add_attribute(domain%info%attributes(i)%name, domain%info%attributes(i)%value)
        enddo

    end subroutine


    module subroutine add_to_output(this, variable)
        class(output_t),   intent(inout) :: this
        type(variable_t),  intent(in)     :: variable

        if (.not.this%is_initialized) call this%init()

        if (associated(variable%data_2d).or.associated(variable%data_3d)) then

            if (this%n_variables == size(this%variables)) call this%increase_var_capacity()

            this%n_variables = this%n_variables + 1

            this%variables(this%n_variables) = variable
        endif

    end subroutine


    module subroutine save_file(this, filename, current_step, time)
        class(output_t),  intent(inout) :: this
        character(len=*), intent(in)    :: filename
        integer,          intent(in)    :: current_step
        type(Time_type),  intent(in)    :: time
        integer :: err, oldmode

        if (.not.this%is_initialized) call this%init()
        
        ! open file
        this%filename = filename
        err = nf90_open(filename, IOR(NF90_WRITE,NF90_NETCDF4), this%ncfile_id, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
        if (err /= NF90_NOERR) then
            call check_ncdf( nf90_create(filename, IOR(NF90_CLOBBER,NF90_NETCDF4), this%ncfile_id, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "Opening:"//trim(filename))
            this%creating=.True.
        else
            ! in case we need to add a new variable when setting up variables
            call check_ncdf(nf90_redef(this%ncfile_id), "Setting redefine mode for: "//trim(filename))
        endif
        
        call check_ncdf( nf90_set_fill(this%ncfile_id, nf90_nofill, oldmode), "Setting fill mode to none")

        ! define variables or find variable IDs (and dimensions)
        call setup_variables(this, time)
        if (this%creating) then
            ! add global attributes such as the image number, domain dimension, creation time
            call add_global_attributes(this)
        endif
        ! End define mode. This tells netCDF we are done defining metadata.
        call check_ncdf( nf90_enddef(this%ncfile_id), "end define mode" )

        ! store output
        call save_data(this, current_step, time)

        this%creating = .false.
        ! close file
        call check_ncdf(nf90_close(this%ncfile_id), "Closing file "//trim(filename))
    end subroutine

    module subroutine add_variables(this, var_list, domain)
        class(output_t), intent(inout)  :: this
        integer,         intent(in)     :: var_list(:)
        type(domain_t),  intent(in)     :: domain

        if (0<var_list( kVARS%u) )                          call this%add_to_output( get_metadata( kVARS%u                            , domain%u))
        if (0<var_list( kVARS%v) )                          call this%add_to_output( get_metadata( kVARS%v                            , domain%v))
        if (0<var_list( kVARS%w) )                          call this%add_to_output( get_metadata( kVARS%w                            , domain%w))
        if (0<var_list( kVARS%w) )                          call this%add_to_output( get_metadata( kVARS%w_real                       , domain%w_real))
        if (0<var_list( kVARS%nsquared) )                   call this%add_to_output( get_metadata( kVARS%nsquared                     , domain%nsquared))
        if (0<var_list( kVARS%water_vapor) )                call this%add_to_output( get_metadata( kVARS%water_vapor                  , domain%water_vapor))
        if (0<var_list( kVARS%potential_temperature) )      call this%add_to_output( get_metadata( kVARS%potential_temperature        , domain%potential_temperature))
        if (0<var_list( kVARS%cloud_water) )                call this%add_to_output( get_metadata( kVARS%cloud_water                  , domain%cloud_water_mass))
        if (0<var_list( kVARS%cloud_number_concentration))  call this%add_to_output( get_metadata( kVARS%cloud_number_concentration   , domain%cloud_number))
        if (0<var_list( kVARS%cloud_ice) )                  call this%add_to_output( get_metadata( kVARS%cloud_ice                    , domain%cloud_ice_mass))
        if (0<var_list( kVARS%ice_number_concentration))    call this%add_to_output( get_metadata( kVARS%ice_number_concentration     , domain%cloud_ice_number))
        if (0<var_list( kVARS%rain_in_air) )                call this%add_to_output( get_metadata( kVARS%rain_in_air                  , domain%rain_mass))
        if (0<var_list( kVARS%rain_number_concentration))   call this%add_to_output( get_metadata( kVARS%rain_number_concentration    , domain%rain_number))
        if (0<var_list( kVARS%snow_in_air) )                call this%add_to_output( get_metadata( kVARS%snow_in_air                  , domain%snow_mass))
        if (0<var_list( kVARS%snow_number_concentration) )  call this%add_to_output( get_metadata( kVARS%snow_number_concentration    , domain%snow_number))
        if (0<var_list( kVARS%graupel_in_air) )             call this%add_to_output( get_metadata( kVARS%graupel_in_air               , domain%graupel_mass))
        if (0<var_list( kVARS%graupel_number_concentration))call this%add_to_output( get_metadata( kVARS%graupel_number_concentration , domain%graupel_number))
        if (0<var_list( kVARS%ice1_a))                      call this%add_to_output( get_metadata( kVARS%ice1_a                       , domain%ice1_a))
        if (0<var_list( kVARS%ice1_c))                      call this%add_to_output( get_metadata( kVARS%ice1_c                       , domain%ice1_c))
        if (0<var_list( kVARS%ice2_mass))                   call this%add_to_output( get_metadata( kVARS%ice2_mass                    , domain%ice2_mass))
        if (0<var_list( kVARS%ice2_number))                 call this%add_to_output( get_metadata( kVARS%ice2_number                  , domain%ice2_number))
        if (0<var_list( kVARS%ice2_a))                      call this%add_to_output( get_metadata( kVARS%ice2_a                       , domain%ice2_a))
        if (0<var_list( kVARS%ice2_c))                      call this%add_to_output( get_metadata( kVARS%ice2_c                       , domain%ice2_c))
        if (0<var_list( kVARS%ice3_mass))                   call this%add_to_output( get_metadata( kVARS%ice3_mass                    , domain%ice3_mass))
        if (0<var_list( kVARS%ice3_number))                 call this%add_to_output( get_metadata( kVARS%ice3_number                  , domain%ice3_number))
        if (0<var_list( kVARS%ice3_a))                      call this%add_to_output( get_metadata( kVARS%ice3_a                       , domain%ice3_a))
        if (0<var_list( kVARS%ice3_c))                      call this%add_to_output( get_metadata( kVARS%ice3_c                       , domain%ice3_c))

        if (0<var_list( kVARS%precipitation) )              call this%add_to_output( get_metadata( kVARS%precipitation                , domain%accumulated_precipitation))
        if (0<var_list( kVARS%convective_precipitation) )   call this%add_to_output( get_metadata( kVARS%convective_precipitation     , domain%accumulated_convective_pcp))
        if (0<var_list( kVARS%snowfall) )                   call this%add_to_output( get_metadata( kVARS%snowfall                     , domain%accumulated_snowfall))
        if (0<var_list( kVARS%graupel) )                    call this%add_to_output( get_metadata( kVARS%graupel                      , domain%graupel))
        if (0<var_list( kVARS%pressure) )                   call this%add_to_output( get_metadata( kVARS%pressure                     , domain%pressure))
        if (0<var_list( kVARS%temperature) )                call this%add_to_output( get_metadata( kVARS%temperature                  , domain%temperature))
        if (0<var_list( kVARS%exner) )                      call this%add_to_output( get_metadata( kVARS%exner                        , domain%exner))
        if (0<var_list( kVARS%z) )                          call this%add_to_output( get_metadata( kVARS%z                            , domain%z))
        if (0<var_list( kVARS%dz_interface) )               call this%add_to_output( get_metadata( kVARS%dz_interface                 , domain%dz_interface))
        if (0<var_list( kVARS%z_interface) )                call this%add_to_output( get_metadata( kVARS%z_interface                  , domain%z_interface))
        if (0<var_list( kVARS%dz) )                         call this%add_to_output( get_metadata( kVARS%dz                           , domain%dz_mass))
        if (0<var_list( kVARS%density) )                    call this%add_to_output( get_metadata( kVARS%density                      , domain%density))
        if (0<var_list( kVARS%pressure_interface) )         call this%add_to_output( get_metadata( kVARS%pressure_interface           , domain%pressure_interface))
        if (0<var_list( kVARS%cloud_fraction) )             call this%add_to_output( get_metadata( kVARS%cloud_fraction               , domain%cloud_fraction))
        if (0<var_list( kVARS%shortwave) )                  call this%add_to_output( get_metadata( kVARS%shortwave                    , domain%shortwave))
        if (0<var_list( kVARS%shortwave_direct) )           call this%add_to_output( get_metadata( kVARS%shortwave_direct             , domain%shortwave_direct))
        if (0<var_list( kVARS%shortwave_diffuse) )          call this%add_to_output( get_metadata( kVARS%shortwave_diffuse            , domain%shortwave_diffuse))
        if (0<var_list( kVARS%longwave) )                   call this%add_to_output( get_metadata( kVARS%longwave                     , domain%longwave))
        if (0<var_list( kVARS%vegetation_fraction) )        call this%add_to_output( get_metadata( kVARS%vegetation_fraction          , domain%vegetation_fraction))
        if (0<var_list( kVARS%vegetation_fraction_max) )    call this%add_to_output( get_metadata( kVARS%vegetation_fraction_max      , domain%vegetation_fraction_max))
        if (0<var_list( kVARS%vegetation_fraction_out) )    call this%add_to_output( get_metadata( kVARS%vegetation_fraction_out      , domain%vegetation_fraction_out))
        if (0<var_list( kVARS%lai) )                        call this%add_to_output( get_metadata( kVARS%lai                          , domain%lai))
        if (0<var_list( kVARS%sai) )                        call this%add_to_output( get_metadata( kVARS%sai                          , domain%sai))
        if (0<var_list( kVARS%crop_type) )                  call this%add_to_output( get_metadata( kVARS%crop_type                    , domain%crop_type))
        if (0<var_list( kVARS%date_planting) )              call this%add_to_output( get_metadata( kVARS%date_planting                , domain%date_planting))
        if (0<var_list( kVARS%date_harvest) )               call this%add_to_output( get_metadata( kVARS%date_harvest                 , domain%date_harvest))
        if (0<var_list( kVARS%growing_season_gdd) )         call this%add_to_output( get_metadata( kVARS%growing_season_gdd           , domain%growing_season_gdd))
        if (0<var_list( kVARS%irr_frac_total) )             call this%add_to_output( get_metadata( kVARS%irr_frac_total               , domain%irr_frac_total))
        if (0<var_list( kVARS%irr_frac_sprinkler) )         call this%add_to_output( get_metadata( kVARS%irr_frac_sprinkler           , domain%irr_frac_sprinkler))
        if (0<var_list( kVARS%irr_frac_micro) )             call this%add_to_output( get_metadata( kVARS%irr_frac_micro               , domain%irr_frac_micro))
        if (0<var_list( kVARS%irr_frac_flood) )             call this%add_to_output( get_metadata( kVARS%irr_frac_flood               , domain%irr_frac_flood))
        if (0<var_list( kVARS%irr_alloc_sprinkler) )        call this%add_to_output( get_metadata( kVARS%irr_alloc_sprinkler          , domain%irr_alloc_sprinkler))
        if (0<var_list( kVARS%irr_alloc_micro) )            call this%add_to_output( get_metadata( kVARS%irr_alloc_micro              , domain%irr_alloc_micro))
        if (0<var_list( kVARS%irr_alloc_flood) )            call this%add_to_output( get_metadata( kVARS%irr_alloc_flood              , domain%irr_alloc_flood))
        if (0<var_list( kVARS%irr_evap_loss_sprinkler) )    call this%add_to_output( get_metadata( kVARS%irr_evap_loss_sprinkler      , domain%irr_evap_loss_sprinkler))
        if (0<var_list( kVARS%irr_amt_sprinkler) )          call this%add_to_output( get_metadata( kVARS%irr_amt_sprinkler            , domain%irr_amt_sprinkler))
        if (0<var_list( kVARS%irr_amt_micro) )              call this%add_to_output( get_metadata( kVARS%irr_amt_micro                , domain%irr_amt_micro))
        if (0<var_list( kVARS%irr_amt_flood) )              call this%add_to_output( get_metadata( kVARS%irr_amt_flood                , domain%irr_amt_flood))
        if (0<var_list( kVARS%evap_heat_sprinkler) )        call this%add_to_output( get_metadata( kVARS%evap_heat_sprinkler          , domain%evap_heat_sprinkler))
        if (0<var_list( kVARS%mass_ag_grain) )              call this%add_to_output( get_metadata( kVARS%mass_ag_grain                , domain%mass_ag_grain))
        if (0<var_list( kVARS%growing_degree_days) )        call this%add_to_output( get_metadata( kVARS%growing_degree_days          , domain%growing_degree_days))
        if (0<var_list( kVARS%net_ecosystem_exchange) )     call this%add_to_output( get_metadata( kVARS%net_ecosystem_exchange       , domain%net_ecosystem_exchange))
        if (0<var_list( kVARS%gross_primary_prod) )         call this%add_to_output( get_metadata( kVARS%gross_primary_prod           , domain%gross_primary_prod))
        if (0<var_list( kVARS%net_primary_prod) )           call this%add_to_output( get_metadata( kVARS%net_primary_prod             , domain%net_primary_prod))
        if (0<var_list( kVARS%apar) )                       call this%add_to_output( get_metadata( kVARS%apar                         , domain%apar))
        if (0<var_list( kVARS%photosynthesis_total) )       call this%add_to_output( get_metadata( kVARS%photosynthesis_total         , domain%photosynthesis_total))
        if (0<var_list( kVARS%stomatal_resist_total) )      call this%add_to_output( get_metadata( kVARS%stomatal_resist_total        , domain%stomatal_resist_total))
        if (0<var_list( kVARS%stomatal_resist_sun) )        call this%add_to_output( get_metadata( kVARS%stomatal_resist_sun          , domain%stomatal_resist_sun))
        if (0<var_list( kVARS%stomatal_resist_shade) )      call this%add_to_output( get_metadata( kVARS%stomatal_resist_shade        , domain%stomatal_resist_shade))
        if (0<var_list( kVARS%gecros_state) )               call this%add_to_output( get_metadata( kVARS%gecros_state                 , domain%gecros_state))
        if (0<var_list( kVARS%canopy_water) )               call this%add_to_output( get_metadata( kVARS%canopy_water                 , domain%canopy_water))
        if (0<var_list( kVARS%canopy_water_ice) )           call this%add_to_output( get_metadata( kVARS%canopy_water_ice             , domain%canopy_water_ice))
        if (0<var_list( kVARS%canopy_water_liquid) )        call this%add_to_output( get_metadata( kVARS%canopy_water_liquid          , domain%canopy_water_liquid))
        if (0<var_list( kVARS%canopy_vapor_pressure) )      call this%add_to_output( get_metadata( kVARS%canopy_vapor_pressure        , domain%canopy_vapor_pressure))
        if (0<var_list( kVARS%canopy_temperature) )         call this%add_to_output( get_metadata( kVARS%canopy_temperature           , domain%canopy_temperature))
        if (0<var_list( kVARS%canopy_fwet) )                call this%add_to_output( get_metadata( kVARS%canopy_fwet                  , domain%canopy_fwet))
        if (0<var_list( kVARS%veg_leaf_temperature) )       call this%add_to_output( get_metadata( kVARS%veg_leaf_temperature         , domain%veg_leaf_temperature))
        if (0<var_list( kVARS%ground_surf_temperature) )    call this%add_to_output( get_metadata( kVARS%ground_surf_temperature      , domain%ground_surf_temperature))
        if (0<var_list( kVARS%frac_within_gap) )            call this%add_to_output( get_metadata( kVARS%frac_within_gap              , domain%frac_within_gap))
        if (0<var_list( kVARS%frac_between_gap) )           call this%add_to_output( get_metadata( kVARS%frac_between_gap             , domain%frac_between_gap))
        if (0<var_list( kVARS%ground_temperature_bare) )    call this%add_to_output( get_metadata( kVARS%ground_temperature_bare      , domain%ground_temperature_bare))
        if (0<var_list( kVARS%ground_temperature_canopy) )  call this%add_to_output( get_metadata( kVARS%ground_temperature_canopy    , domain%ground_temperature_canopy))
        if (0<var_list( kVARS%snowfall_ground) )            call this%add_to_output( get_metadata( kVARS%snowfall_ground              , domain%snowfall_ground))
        if (0<var_list( kVARS%rainfall_ground) )            call this%add_to_output( get_metadata( kVARS%rainfall_ground              , domain%rainfall_ground))
        if (0<var_list( kVARS%snow_water_equivalent) )      call this%add_to_output( get_metadata( kVARS%snow_water_equivalent        , domain%snow_water_equivalent))
        if (0<var_list( kVARS%snow_water_eq_prev) )         call this%add_to_output( get_metadata( kVARS%snow_water_eq_prev           , domain%snow_water_eq_prev))
        if (0<var_list( kVARS%snow_albedo_prev) )           call this%add_to_output( get_metadata( kVARS%snow_albedo_prev             , domain%snow_albedo_prev))
        if (0<var_list( kVARS%snow_temperature) )           call this%add_to_output( get_metadata( kVARS%snow_temperature             , domain%snow_temperature))
        if (0<var_list( kVARS%snow_layer_depth) )           call this%add_to_output( get_metadata( kVARS%snow_layer_depth             , domain%snow_layer_depth))
        if (0<var_list( kVARS%snow_layer_ice) )             call this%add_to_output( get_metadata( kVARS%snow_layer_ice               , domain%snow_layer_ice))
        if (0<var_list( kVARS%snow_layer_liquid_water) )    call this%add_to_output( get_metadata( kVARS%snow_layer_liquid_water      , domain%snow_layer_liquid_water))
        if (0<var_list( kVARS%snow_age_factor) )            call this%add_to_output( get_metadata( kVARS%snow_age_factor              , domain%snow_age_factor))
        if (0<var_list( kVARS%snow_height) )                call this%add_to_output( get_metadata( kVARS%snow_height                  , domain%snow_height))
        if (0<var_list( kVARS%skin_temperature) )           call this%add_to_output( get_metadata( kVARS%skin_temperature             , domain%skin_temperature))
        if (0<var_list( kVARS%soil_water_content) )         call this%add_to_output( get_metadata( kVARS%soil_water_content           , domain%soil_water_content))
        if (0<var_list( kVARS%eq_soil_moisture) )           call this%add_to_output( get_metadata( kVARS%eq_soil_moisture             , domain%eq_soil_moisture))
        if (0<var_list( kVARS%smc_watertable_deep) )        call this%add_to_output( get_metadata( kVARS%smc_watertable_deep          , domain%smc_watertable_deep))
        if (0<var_list( kVARS%recharge) )                   call this%add_to_output( get_metadata( kVARS%recharge                     , domain%recharge))
        if (0<var_list( kVARS%recharge_deep) )              call this%add_to_output( get_metadata( kVARS%recharge_deep                , domain%recharge_deep))
        if (0<var_list( kVARS%soil_temperature) )           call this%add_to_output( get_metadata( kVARS%soil_temperature             , domain%soil_temperature))
        if (0<var_list( kVARS%latitude) )                   call this%add_to_output( get_metadata( kVARS%latitude                     , domain%latitude))
        if (0<var_list( kVARS%longitude) )                  call this%add_to_output( get_metadata( kVARS%longitude                    , domain%longitude))
        if (0<var_list( kVARS%u_latitude) )                 call this%add_to_output( get_metadata( kVARS%u_latitude                   , domain%u_latitude))
        if (0<var_list( kVARS%u_longitude) )                call this%add_to_output( get_metadata( kVARS%u_longitude                  , domain%u_longitude))
        if (0<var_list( kVARS%v_latitude) )                 call this%add_to_output( get_metadata( kVARS%v_latitude                   , domain%v_latitude))
        if (0<var_list( kVARS%v_longitude) )                call this%add_to_output( get_metadata( kVARS%v_longitude                  , domain%v_longitude))
        if (0<var_list( kVARS%terrain) )                    call this%add_to_output( get_metadata( kVARS%terrain                      , domain%terrain))
        if (0<var_list( kVARS%sensible_heat) )              call this%add_to_output( get_metadata( kVARS%sensible_heat                , domain%sensible_heat))
        if (0<var_list( kVARS%latent_heat) )                call this%add_to_output( get_metadata( kVARS%latent_heat                  , domain%latent_heat))
        if (0<var_list( kVARS%u_10m) )                      call this%add_to_output( get_metadata( kVARS%u_10m                        , domain%u_10m))
        if (0<var_list( kVARS%v_10m) )                      call this%add_to_output( get_metadata( kVARS%v_10m                        , domain%v_10m))
        if (0<var_list( kVARS%coeff_momentum_drag) )        call this%add_to_output( get_metadata( kVARS%coeff_momentum_drag          , domain%coeff_momentum_drag))
        if (0<var_list( kVARS%coeff_heat_exchange) )        call this%add_to_output( get_metadata( kVARS%coeff_heat_exchange          , domain%coeff_heat_exchange))
        if (0<var_list( kVARS%surface_rad_temperature) )    call this%add_to_output( get_metadata( kVARS%surface_rad_temperature      , domain%surface_rad_temperature))
        if (0<var_list( kVARS%temperature_2m) )             call this%add_to_output( get_metadata( kVARS%temperature_2m               , domain%temperature_2m))
        if (0<var_list( kVARS%humidity_2m) )                call this%add_to_output( get_metadata( kVARS%humidity_2m                  , domain%humidity_2m))
        if (0<var_list( kVARS%temperature_2m_veg) )         call this%add_to_output( get_metadata( kVARS%temperature_2m_veg           , domain%temperature_2m_veg))
        if (0<var_list( kVARS%temperature_2m_bare) )        call this%add_to_output( get_metadata( kVARS%temperature_2m_bare          , domain%temperature_2m_bare))
        if (0<var_list( kVARS%mixing_ratio_2m_veg) )        call this%add_to_output( get_metadata( kVARS%mixing_ratio_2m_veg          , domain%mixing_ratio_2m_veg))
        if (0<var_list( kVARS%mixing_ratio_2m_bare) )       call this%add_to_output( get_metadata( kVARS%mixing_ratio_2m_bare         , domain%mixing_ratio_2m_bare))
        if (0<var_list( kVARS%surface_pressure) )           call this%add_to_output( get_metadata( kVARS%surface_pressure             , domain%surface_pressure))
        if (0<var_list( kVARS%rad_absorbed_total) )         call this%add_to_output( get_metadata( kVARS%rad_absorbed_total           , domain%rad_absorbed_total))
        if (0<var_list( kVARS%rad_absorbed_veg) )           call this%add_to_output( get_metadata( kVARS%rad_absorbed_veg             , domain%rad_absorbed_veg))
        if (0<var_list( kVARS%rad_absorbed_bare) )          call this%add_to_output( get_metadata( kVARS%rad_absorbed_bare            , domain%rad_absorbed_bare))
        if (0<var_list( kVARS%rad_net_longwave) )           call this%add_to_output( get_metadata( kVARS%rad_net_longwave             , domain%rad_net_longwave))
        if (0<var_list( kVARS%longwave_up) )                call this%add_to_output( get_metadata( kVARS%longwave_up                  , domain%longwave_up))
        if (0<var_list( kVARS%ground_heat_flux) )           call this%add_to_output( get_metadata( kVARS%ground_heat_flux             , domain%ground_heat_flux))
        if (0<var_list( kVARS%soil_deep_temperature) )      call this%add_to_output( get_metadata( kVARS%soil_deep_temperature        , domain%soil_deep_temperature))
        if (0<var_list( kVARS%evap_canopy) )                call this%add_to_output( get_metadata( kVARS%evap_canopy                  , domain%evap_canopy))
        if (0<var_list( kVARS%evap_soil_surface) )          call this%add_to_output( get_metadata( kVARS%evap_soil_surface            , domain%evap_soil_surface))
        if (0<var_list( kVARS%transpiration_rate) )         call this%add_to_output( get_metadata( kVARS%transpiration_rate           , domain%transpiration_rate))
        if (0<var_list( kVARS%ch_veg) )                     call this%add_to_output( get_metadata( kVARS%ch_veg                       , domain%ch_veg))
        if (0<var_list( kVARS%ch_veg_2m) )                  call this%add_to_output( get_metadata( kVARS%ch_veg_2m                    , domain%ch_veg_2m))
        if (0<var_list( kVARS%ch_bare) )                    call this%add_to_output( get_metadata( kVARS%ch_bare                      , domain%ch_bare))
        if (0<var_list( kVARS%ch_bare_2m) )                 call this%add_to_output( get_metadata( kVARS%ch_bare_2m                   , domain%ch_bare_2m))
        if (0<var_list( kVARS%ch_under_canopy) )            call this%add_to_output( get_metadata( kVARS%ch_under_canopy              , domain%ch_under_canopy))
        if (0<var_list( kVARS%ch_leaf) )                    call this%add_to_output( get_metadata( kVARS%ch_leaf                      , domain%ch_leaf))
        if (0<var_list( kVARS%sensible_heat_veg) )          call this%add_to_output( get_metadata( kVARS%sensible_heat_veg            , domain%sensible_heat_veg))
        if (0<var_list( kVARS%sensible_heat_bare) )         call this%add_to_output( get_metadata( kVARS%sensible_heat_bare           , domain%sensible_heat_bare))
        if (0<var_list( kVARS%sensible_heat_canopy) )       call this%add_to_output( get_metadata( kVARS%sensible_heat_canopy         , domain%sensible_heat_canopy))
        if (0<var_list( kVARS%evap_heat_veg) )              call this%add_to_output( get_metadata( kVARS%evap_heat_veg                , domain%evap_heat_veg))
        if (0<var_list( kVARS%evap_heat_bare) )             call this%add_to_output( get_metadata( kVARS%evap_heat_bare               , domain%evap_heat_bare))
        if (0<var_list( kVARS%evap_heat_canopy) )           call this%add_to_output( get_metadata( kVARS%evap_heat_canopy             , domain%evap_heat_canopy))
        if (0<var_list( kVARS%transpiration_heat) )         call this%add_to_output( get_metadata( kVARS%transpiration_heat           , domain%transpiration_heat))
        if (0<var_list( kVARS%ground_heat_veg) )            call this%add_to_output( get_metadata( kVARS%ground_heat_veg              , domain%ground_heat_veg))
        if (0<var_list( kVARS%ground_heat_bare) )           call this%add_to_output( get_metadata( kVARS%ground_heat_bare             , domain%ground_heat_bare))
        if (0<var_list( kVARS%net_longwave_veg) )           call this%add_to_output( get_metadata( kVARS%net_longwave_veg             , domain%net_longwave_veg))
        if (0<var_list( kVARS%net_longwave_bare) )          call this%add_to_output( get_metadata( kVARS%net_longwave_bare            , domain%net_longwave_bare))
        if (0<var_list( kVARS%net_longwave_canopy) )        call this%add_to_output( get_metadata( kVARS%net_longwave_canopy          , domain%net_longwave_canopy))
        if (0<var_list( kVARS%runoff_surface) )             call this%add_to_output( get_metadata( kVARS%runoff_surface               , domain%runoff_surface))
        if (0<var_list( kVARS%runoff_subsurface) )          call this%add_to_output( get_metadata( kVARS%runoff_subsurface            , domain%runoff_subsurface))
        if (0<var_list( kVARS%soil_totalmoisture) )         call this%add_to_output( get_metadata( kVARS%soil_totalmoisture           , domain%soil_totalmoisture))
        if (0<var_list( kVARS%water_table_depth) )          call this%add_to_output( get_metadata( kVARS%water_table_depth            , domain%water_table_depth))
        if (0<var_list( kVARS%water_aquifer) )              call this%add_to_output( get_metadata( kVARS%water_aquifer                , domain%water_aquifer))
        if (0<var_list( kVARS%storage_gw) )                 call this%add_to_output( get_metadata( kVARS%storage_gw                   , domain%storage_gw))
        if (0<var_list( kVARS%storage_lake) )               call this%add_to_output( get_metadata( kVARS%storage_lake                 , domain%storage_lake))
        if (0<var_list( kVARS%roughness_z0) )               call this%add_to_output( get_metadata( kVARS%roughness_z0                 , domain%roughness_z0))
        if (0<var_list( kVARS%mass_leaf) )                  call this%add_to_output( get_metadata( kVARS%mass_leaf                    , domain%mass_leaf))
        if (0<var_list( kVARS%mass_root) )                  call this%add_to_output( get_metadata( kVARS%mass_root                    , domain%mass_root))
        if (0<var_list( kVARS%mass_stem) )                  call this%add_to_output( get_metadata( kVARS%mass_stem                    , domain%mass_stem))
        if (0<var_list( kVARS%mass_wood) )                  call this%add_to_output( get_metadata( kVARS%mass_wood                    , domain%mass_wood))
        if (0<var_list( kVARS%soil_carbon_fast) )           call this%add_to_output( get_metadata( kVARS%soil_carbon_fast             , domain%soil_carbon_fast))
        if (0<var_list( kVARS%soil_carbon_stable) )         call this%add_to_output( get_metadata( kVARS%soil_carbon_stable           , domain%soil_carbon_stable))
        if (0<var_list( kVARS%soil_texture_1) )             call this%add_to_output( get_metadata( kVARS%soil_texture_1               , domain%soil_texture_1))
        if (0<var_list( kVARS%soil_texture_2) )             call this%add_to_output( get_metadata( kVARS%soil_texture_2               , domain%soil_texture_2))
        if (0<var_list( kVARS%soil_texture_3) )             call this%add_to_output( get_metadata( kVARS%soil_texture_3               , domain%soil_texture_3))
        if (0<var_list( kVARS%soil_texture_4) )             call this%add_to_output( get_metadata( kVARS%soil_texture_4               , domain%soil_texture_4))
        if (0<var_list( kVARS%soil_sand_and_clay) )         call this%add_to_output( get_metadata( kVARS%soil_sand_and_clay           , domain%soil_sand_and_clay))
        if (0<var_list( kVARS%re_cloud) )                   call this%add_to_output( get_metadata( kVARS%re_cloud                     , domain%re_cloud))
        if (0<var_list( kVARS%re_ice) )                     call this%add_to_output( get_metadata( kVARS%re_ice                       , domain%re_ice))
        if (0<var_list( kVARS%re_snow) )                    call this%add_to_output( get_metadata( kVARS%re_snow                      , domain%re_snow))
        if (0<var_list( kVARS%out_longwave_rad) )           call this%add_to_output( get_metadata( kVARS%out_longwave_rad             , domain%out_longwave_rad))
        if (0<var_list( kVARS%longwave_cloud_forcing) )     call this%add_to_output( get_metadata( kVARS%longwave_cloud_forcing       , domain%longwave_cloud_forcing))
        if (0<var_list( kVARS%shortwave_cloud_forcing) )    call this%add_to_output( get_metadata( kVARS%shortwave_cloud_forcing      , domain%shortwave_cloud_forcing))
        if (0<var_list( kVARS%cosine_zenith_angle) )        call this%add_to_output( get_metadata( kVARS%cosine_zenith_angle          , domain%cosine_zenith_angle))
        if (0<var_list( kVARS%land_emissivity) )            call this%add_to_output( get_metadata( kVARS%land_emissivity              , domain%land_emissivity))
        if (0<var_list( kVARS%temperature_interface) )      call this%add_to_output( get_metadata( kVARS%temperature_interface        , domain%temperature_interface))
        if (0<var_list( kVARS%tend_swrad) )                 call this%add_to_output( get_metadata( kVARS%tend_swrad                  , domain%tend_swrad))

    end subroutine

    subroutine add_global_attributes(this)
        implicit none
        class(output_t), intent(inout)  :: this
        integer :: i

        character(len=19)       :: todays_date_time
        integer,dimension(8)    :: date_time
        character(len=49)       :: date_format
        character(len=5)        :: UTCoffset
        character(len=64)       :: err
        integer                 :: ncid

        ncid = this%ncfile_id

        err="Creating global attributes"
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"Conventions","CF-1.6"), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"title","Intermediate Complexity Atmospheric Research (ICAR) model output"), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"institution","National Center for Atmospheric Research"), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"history","Created:"//todays_date_time//UTCoffset), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"references", &
                    "Gutmann et al. 2016: The Intermediate Complexity Atmospheric Model (ICAR). J.Hydrometeor. doi:10.1175/JHM-D-15-0155.1, 2016."), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"contact","Ethan Gutmann : gutmann@ucar.edu"), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"git",VERSION), trim(err))

        if (this%n_attrs > 0) then
            do i=1,this%n_attrs
                call check_ncdf( nf90_put_att(   this%ncfile_id,             &
                                            NF90_GLOBAL,                &
                                            trim(this%attributes(i)%name),    &
                                            trim(this%attributes(i)%value)),  &
                                            "global attr:"//trim(this%attributes(i)%name))
            enddo
        endif

        call date_and_time(values=date_time,zone=UTCoffset)
        date_format='(I4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)'
        write(todays_date_time,date_format) date_time(1:3),date_time(5:7)

        call check_ncdf(nf90_put_att(this%ncfile_id, NF90_GLOBAL,"history","Created:"//todays_date_time//UTCoffset), "global attr")
        !call check_ncdf(nf90_put_att(this%ncfile_id, NF90_GLOBAL, "image", this_image()))

    end subroutine add_global_attributes

    subroutine setup_variables(this, time)
        implicit none
        class(output_t), intent(inout) :: this
        type(Time_type), intent(in)    :: time
        integer :: i

        ! iterate through variables creating or setting up variable IDs if they exist, also dimensions
        do i=1,this%n_variables
            ! create all dimensions or find dimension IDs if they exist already

            call setup_dims_for_var(this, this%variables(i))

            call setup_variable(this, this%variables(i))
            
        end do

        call setup_time_variable(this, time)

    end subroutine setup_variables

    subroutine setup_time_variable(this, time)
        implicit none
        class(output_t), intent(inout) :: this
        type(Time_type), intent(in)    :: time
        integer :: err
        character(len=kMAX_NAME_LENGTH) :: calendar

        associate(var => this%time)
        var%name = "time"
        var%dimensions = [ "time" ]
        var%n_dimensions = 1

        select case (time%calendar)
            case(GREGORIAN)
                calendar = "proleptic_gregorian"
            case(NOLEAP)
                calendar = "noleap"
            case(THREESIXTY)
                calendar = "360-day"
            case default
                calendar = "standard"
        end select


        err = nf90_inq_varid(this%ncfile_id, var%name, var%var_id)

        ! if the variable was not found in the netcdf file then we will define it.
        if (err /= NF90_NOERR) then

            if (allocated(var%dim_ids)) deallocate(var%dim_ids)
            allocate(var%dim_ids(1))

            ! Try to find the dimension ID if it exists already.
            err = nf90_inq_dimid(this%ncfile_id, trim(var%dimensions(1)), var%dim_ids(1))

            ! if the dimension doesn't exist in the file, create it.
            if (err/=NF90_NOERR) then
                call check_ncdf( nf90_def_dim(this%ncfile_id, trim(var%dimensions(1)), NF90_UNLIMITED, &
                            var%dim_ids(1) ), "def_dim"//var%dimensions(1) )
            endif

            call check_ncdf( nf90_def_var(this%ncfile_id, var%name, NF90_DOUBLE, var%dim_ids(1), var%var_id), "Defining time" )
            call check_ncdf( nf90_put_att(this%ncfile_id, var%var_id,"standard_name","time"))
            call check_ncdf( nf90_put_att(this%ncfile_id, var%var_id,"calendar",trim(calendar)))
            call check_ncdf( nf90_put_att(this%ncfile_id, var%var_id,"units",time%units()))
            call check_ncdf( nf90_put_att(this%ncfile_id, var%var_id,"UTCoffset","0"))

        endif
        call check_ncdf( nf90_var_par_access(this%ncfile_id, var%var_id, nf90_collective))
        end associate

    end subroutine setup_time_variable


    subroutine save_data(this, current_step, time)
        implicit none
        class(output_t), intent(in) :: this
        integer,         intent(in) :: current_step
        type(Time_type), intent(in) :: time
        integer :: i, i_s, i_e, k_s, k_e, j_s, j_e

        integer :: start_three_D_t(4)
        integer :: start_two_D_t(3)
        integer :: cnt_3d(3)
        integer :: cnt_2d(2)


        do i=1,this%n_variables
            associate(var => this%variables(i))
                i_s = var%grid%its
                i_e = var%grid%ite
                j_s = var%grid%jts
                j_e = var%grid%jte
        
                if (var%grid%ims == var%grid%ids) i_s = var%grid%ids
                if (var%grid%ime == var%grid%ide) i_e = var%grid%ide
                if (var%grid%jms == var%grid%jds) j_s = var%grid%jds
                if (var%grid%jme == var%grid%jde) j_e = var%grid%jde


                if (var%three_d) then
                
                    k_s = var%grid%kts
                    k_e = var%grid%kte
                    if (var%grid%kms == var%grid%kds) k_s = var%grid%kds
                    if (var%grid%kme == var%grid%kde) k_e = var%grid%kde
                    start_three_D_t = (/ i_s, j_s, k_s, current_step /)
                    cnt_3d = (/ (i_e-i_s+1), (j_e-j_s+1), (k_e-k_s+1)  /)


                    if (var%unlimited_dim) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id, &
                            reshape(var%data_3d(i_s:i_e,k_s:k_e,j_s:j_e), shape=cnt_3d, order=[1,3,2]), &
                                        start_three_D_t, count=(/cnt_3d(1), cnt_3d(2), cnt_3d(3), 1/)), "saving:"//trim(var%name) )
                    elseif (this%creating) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id, &
                            reshape(var%data_3d(i_s:i_e,k_s:k_e,j_s:j_e),  &
                            shape=cnt_3d, order=[1,3,2]), start=(/ start_three_D_t(1), start_three_D_t(2), start_three_D_t(3) /),&
                            count=cnt_3d ), "saving:"//trim(var%name) )
                    endif

                elseif (var%two_d) then
                    start_two_D_t = (/ i_s, j_s, current_step /)
                    cnt_2d = (/ (i_e-i_s+1), (j_e-j_s+1) /)

                    if (var%unlimited_dim) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d(i_s:i_e,j_s:j_e), &
                                start_two_D_t,count=(/ cnt_2d(1), cnt_2d(2), 1/)), "saving:"//trim(var%name) )
                    elseif (this%creating) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d(i_s:i_e,j_s:j_e), &
                                    start=(/ start_two_D_t(1), start_two_D_t(2) /), &
                                    count=cnt_2d), "saving:"//trim(var%name) )
                    endif
                endif
            end associate
        end do

        call check_ncdf( nf90_put_var(this%ncfile_id, this%time%var_id, dble(time%mjd()), [current_step]),   &
                   "saving:"//trim(this%time%name) )


    end subroutine save_data

    subroutine setup_dims_for_var(this, var)
        implicit none
        class(output_t),    intent(inout) :: this
        type(variable_t),   intent(inout) :: var
        integer :: i, err
        

        if (allocated(var%dim_ids)) deallocate(var%dim_ids)

        allocate(var%dim_ids(var%n_dimensions))

        do i = 1, size(var%dim_ids)

            ! Try to find the dimension ID if it exists already.
            err = nf90_inq_dimid(this%ncfile_id, trim(var%dimensions(i)), var%dim_ids(i))

            ! probably the dimension doesn't exist in the file, so we will create it.
            if (err/=NF90_NOERR) then
                ! assume that the last dimension should be the unlimited dimension (generally a good idea...)
                if (var%unlimited_dim .and. (i==size(var%dim_ids))) then
                    call check_ncdf( nf90_def_dim(this%ncfile_id, trim(var%dimensions(i)), NF90_UNLIMITED, &
                                var%dim_ids(i) ), "def_dim"//var%dimensions(i) )
                else
                    call check_ncdf( nf90_def_dim(this%ncfile_id, var%dimensions(i), var%global_dim_len(i),       &
                                var%dim_ids(i) ), "def_dim"//var%dimensions(i) )
                endif
            endif
        end do

    end subroutine setup_dims_for_var

    subroutine setup_variable(this, var)
        implicit none
        class(output_t),   intent(inout) :: this
        type(variable_t),  intent(inout) :: var
        integer :: i, n, err
        integer, allocatable :: chunks(:)

        err = nf90_inq_varid(this%ncfile_id, var%name, var%var_id)

        ! if the variable was not found in the netcdf file then we will define it.
        if (err /= NF90_NOERR) then
                    
            allocate(chunks(var%n_dimensions))
            chunks(1) = var%global_dim_len(1); chunks(2) = var%global_dim_len(2);
            if (var%three_d) chunks(3) = var%global_dim_len(3)
            if (var%n_dimensions > size(var%dim_len)) chunks(var%n_dimensions) = 1

            
            call check_ncdf( nf90_def_var(this%ncfile_id, var%name, NF90_REAL, var%dim_ids, var%var_id, chunksizes=chunks), &
                        "Defining variable:"//trim(var%name) )

            ! setup attributes
            do i=1,size(var%attributes)
                call check_ncdf( nf90_put_att(this%ncfile_id,                &
                                         var%var_id,                    &
                                         trim(var%attributes(i)%name),        &
                                         trim(var%attributes(i)%value)),      &
                            "saving attribute"//trim(var%attributes(i)%name))
            enddo
        endif
        
        call check_ncdf( nf90_var_par_access(this%ncfile_id, var%var_id, nf90_collective))


    end subroutine setup_variable

    module subroutine init(this)
        implicit none
        class(output_t),   intent(inout)  :: this
        
        allocate(this%variables(kINITIAL_VAR_SIZE))
        this%n_variables = 0
        this%n_dims      = 0
        this%is_initialized = .True.


    end subroutine

    module subroutine increase_var_capacity(this)
        implicit none
        class(output_t),   intent(inout)  :: this
        type(variable_t),  allocatable :: new_variables(:)

        ! assert allocated(this%variables)
        allocate(new_variables, source=this%variables)
        ! new_variables = this%variables

        deallocate(this%variables)

        allocate(this%variables(size(new_variables)*2))
        this%variables(:size(new_variables)) = new_variables

        deallocate(new_variables)

    end subroutine


end submodule
