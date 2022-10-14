submodule(output_interface) output_implementation
  use output_metadata,          only : get_metadata
  use debug_module,             only : check_ncdf
  use iso_fortran_env,          only : output_unit

  implicit none

contains


    module subroutine init(this,domain,options, its, ite, kts, kte, jts, jte)
        implicit none
        class(output_t),  intent(inout)  :: this
        type(domain_t),   intent(in)     :: domain
        type(options_t),  intent(in)     :: options
        integer,          intent(in)     :: its, ite, kts, kte, jts, jte
        
        integer :: i
        
        allocate(this%variables(kINITIAL_VAR_SIZE))
        
        this%n_vars = 0
        this%n_dims      = 0
        this%is_initialized = .True.
        this%output_counter = 1
        this%output_count = options%io_options%frames_per_outfile
        this%its = its; this%ite = ite; this%kts = kts; this%kte = kte; this%jts = jts; this%jte = jte
        this%global_dim_len = (/domain%ide, domain%jde, domain%kde /)
        
        call set_attrs(this, domain)
        call add_variables(this, options%io_options%vars_for_output)
        
        
        
        !do i=1, size(options%io_options%vars_for_output)
        !    if (options%io_options%vars_for_output(i) > 0) then
        !        call this%add_to_output( get_metadata( i ))
        !    endif
        !enddo

        !call open_file(this,options%parameters%start_time)

    end subroutine

    module subroutine set_attrs(this, domain)
        class(output_t),  intent(inout)  :: this
        type(domain_t),   intent(in)     :: domain
        integer :: i

        do i=1,domain%info%n_attrs
            call this%add_attribute(domain%info%attributes(i)%name, domain%info%attributes(i)%value)
        enddo

    end subroutine


    module subroutine add_to_output(this, in_variable)
        class(output_t),   intent(inout) :: this
        type(variable_t),  intent(in) :: in_variable
        
        type(variable_t) :: variable
        
        variable = in_variable
        
        if (variable%dim_len(3)<=0) variable%dim_len(3) = this%kte
            
        if (this%n_vars == size(this%variables)) call this%increase_var_capacity()

        this%n_vars = this%n_vars + 1

        this%variables(this%n_vars) = variable

        !if (associated(variable%data_2d).or.associated(variable%data_3d)) then

        !endif

    end subroutine


    module subroutine save_file(this, time, par_comms)
        class(output_t),  intent(inout) :: this
        type(Time_type),  intent(in)  :: time
        integer,          intent(in)  :: par_comms

        if (this%ncfile_id < 0) call open_file(this,time, par_comms)
        call flush(output_unit)

        if (.not.this%block_checked) call block_hunter(this)

        ! store output
        call save_data(this, this%output_counter, time)
        
        !In case we had creating set to true, set to false
        this%creating = .false.

        !Check if we should close the file
        if (this%output_counter >= this%output_count) then
            this%output_counter = 1
            call this%close_file()
        else
            this%output_counter = this%output_counter + 1
        endif

    end subroutine
    
    !See if this outputer has 'blocks', or a piece of it which it should not write
    !This can occur due to using only images on a given node, which results in stepped
    !patterns if the output is not appropriately treated
    !
    !The goal of the routine is to set all start and cnt object fields so that either
    !(1) if there is no block, the standard fields contain the whole bounds of the object
    !(2) if there is a block, the standard fields contain the continuous block of the object
    ! and the 'b' fields contain the discontinuous block
    subroutine block_hunter(this)
        class(output_t),  intent(inout) :: this

        integer :: i_s, i_e, k_s, k_e, j_s, j_e, nx, ny
        integer :: i_s_b, i_e_b, j_s_b, j_e_b
        integer :: i_s_b2, i_e_b2, j_s_b2, j_e_b2

        real, allocatable, dimension(:,:) :: datas
        logical :: blocked_LL, blocked_UR
        
        this%block_checked = .True.
        this%is_blocked = .False.
        blocked_LL = .False.
        blocked_UR = .False.

        i_s_b = 2; i_e_b = 1; j_s_b = 2; j_e_b = 1;
        i_s_b2 = 2; i_e_b2 = 1; j_s_b2 = 2; j_e_b2 = 1;

        i_s = this%its
        i_e = this%ite
        j_s = this%jts
        j_e = this%jte
        k_s = this%kts
        k_e = this%kte

        nx = i_e - i_s + 1
        ny = j_e - j_s + 1

        if (this%variables(1)%three_d) then
            datas = this%variables(1)%data_3d(1:nx,1,1:ny)
        else if (this%variables(1)%two_d) then
            datas = this%variables(1)%data_2d(1:nx,1:ny)
        endif
        !Check each corner to see where this starts
        !LL
        if (datas(1,1) == kEMPT_BUFF) then
            i_s_b = findloc(datas(:,1),kEMPT_BUFF,dim=1,back=.True.) + 1
            i_e_b = i_e
            j_e_b = findloc(datas(1,:),kEMPT_BUFF,dim=1,back=.True.) + j_s - 1
            j_s_b = j_s
            blocked_LL = .True.
        endif
        !UR
        if (datas(nx,ny) == kEMPT_BUFF) then
            i_s_b2 = i_s 
            i_e_b2 = findloc(datas(:,ny),kEMPT_BUFF,dim=1,back=.False.) - 1
            j_e_b2 = j_e 
            j_s_b2 = findloc(datas(nx,:),kEMPT_BUFF,dim=1,back=.False.) + j_s + 1
            blocked_UR = .True.
        endif
        
        if (blocked_LL) j_s = j_e_b+1
        if (blocked_UR) j_e = j_s_b2-1
        
        this%start_3d = (/ i_s, j_s, k_s /)
        this%cnt_3d = (/ (i_e-i_s+1), (j_e-j_s+1), (k_e-k_s+1)  /)
        this%cnt_2d = (/ (i_e-i_s+1), (j_e-j_s+1) /)

        !Compute block start and cnts accordingly
        if (blocked_LL) then
            this%start_3d_b = (/ i_s_b, j_s_b, k_s /)
            this%cnt_3d_b = (/ (i_e_b-i_s_b+1), (j_e_b-j_s_b+1), (k_e-k_s+1)  /)
            this%cnt_2d_b = (/ (i_e_b-i_s_b+1), (j_e_b-j_s_b+1) /)
        else
            this%start_3d_b = (/ 1, 1, 1 /)
            this%cnt_3d_b = (/ 0, 0, 0  /)
            this%cnt_2d_b = (/ 0, 0 /)
        endif 
        if (blocked_UR) then
            this%start_3d_b2 = (/ i_s_b2, j_s_b2, k_s /)
            this%cnt_3d_b2 = (/ (i_e_b2-i_s_b2+1), (j_e_b2-j_s_b2+1), (k_e-k_s+1)  /)
            this%cnt_2d_b2 = (/ (i_e_b2-i_s_b2+1), (j_e_b2-j_s_b2+1) /)
        else
            this%start_3d_b2 = (/ 1, 1, 1 /)
            this%cnt_3d_b2 = (/ 0, 0, 0  /)
            this%cnt_2d_b2 = (/ 0, 0 /)
        endif
    end subroutine block_hunter
    
    subroutine open_file(this, time, par_comms)
        class(output_t),  intent(inout) :: this
        type(Time_type),  intent(in)    :: time
        integer,          intent(in)    :: par_comms

        integer :: err
        
        write(this%filename, '(A,A,".nc")')    &
            trim(this%base_file_name),   &
            trim(time%as_string(this%file_date_format))

        ! open file
        err = nf90_open(this%filename, IOR(NF90_WRITE,NF90_NETCDF4), this%ncfile_id, comm = par_comms, info = MPI_INFO_NULL)
        if (err /= NF90_NOERR) then
            call check_ncdf( nf90_create(this%filename, IOR(NF90_CLOBBER,NF90_NETCDF4), this%ncfile_id, comm = par_comms, info = MPI_INFO_NULL), "Opening:"//trim(this%filename))
            this%creating=.True.
        else
            ! in case we need to add a new variable when setting up variables
            call check_ncdf(nf90_redef(this%ncfile_id), "Setting redefine mode for: "//trim(this%filename))
        endif

        !This command cuts down on write time and must be done once the file is opened
        call check_ncdf( nf90_set_fill(this%ncfile_id, nf90_nofill, err), "Setting fill mode to none")

        ! define variables or find variable IDs (and dimensions)
        call setup_variables(this, time)

        if (this%creating) then
            ! add global attributes such as the image number, domain dimension, creation time
            call add_global_attributes(this)
        endif
        ! End define mode. This tells netCDF we are done defining metadata.
        call check_ncdf( nf90_enddef(this%ncfile_id), "end define mode" )
    

    end subroutine

    module subroutine add_variables(this, var_list)
        class(output_t), intent(inout)  :: this
        integer,         intent(in)     :: var_list(:)



        if (0<var_list( kVARS%u) )                          call this%add_to_output( get_metadata( kVARS%u ))                        
        if (0<var_list( kVARS%v) )                          call this%add_to_output( get_metadata( kVARS%v ))                        
        if (0<var_list( kVARS%w) )                          call this%add_to_output( get_metadata( kVARS%w ))                        
        if (0<var_list( kVARS%w) )                          call this%add_to_output( get_metadata( kVARS%w_real))                    
        if (0<var_list( kVARS%nsquared) )                   call this%add_to_output( get_metadata( kVARS%nsquared  ))     
        if (0<var_list( kVARS%water_vapor) )                call this%add_to_output( get_metadata( kVARS%water_vapor  )) 
        if (0<var_list( kVARS%potential_temperature) )      call this%add_to_output( get_metadata( kVARS%potential_temperature ))
        if (0<var_list( kVARS%cloud_water) )                call this%add_to_output( get_metadata( kVARS%cloud_water ))
        if (0<var_list( kVARS%cloud_number_concentration))  call this%add_to_output( get_metadata( kVARS%cloud_number_concentration))
        if (0<var_list( kVARS%cloud_ice) )                  call this%add_to_output( get_metadata( kVARS%cloud_ice ))
        if (0<var_list( kVARS%ice_number_concentration))    call this%add_to_output( get_metadata( kVARS%ice_number_concentration ))
        if (0<var_list( kVARS%rain_in_air) )                call this%add_to_output( get_metadata( kVARS%rain_in_air               ))    
        if (0<var_list( kVARS%rain_number_concentration))   call this%add_to_output( get_metadata( kVARS%rain_number_concentration ))    
        if (0<var_list( kVARS%snow_in_air) )                call this%add_to_output( get_metadata( kVARS%snow_in_air               ))    
        if (0<var_list( kVARS%snow_number_concentration) )  call this%add_to_output( get_metadata( kVARS%snow_number_concentration ))    
        if (0<var_list( kVARS%graupel_in_air) )             call this%add_to_output( get_metadata( kVARS%graupel_in_air            ))
        if (0<var_list( kVARS%graupel_number_concentration))call this%add_to_output( get_metadata( kVARS%graupel_number_concentration ))
        if (0<var_list( kVARS%ice1_a))                      call this%add_to_output( get_metadata( kVARS%ice1_a ))                   
        if (0<var_list( kVARS%ice1_c))                      call this%add_to_output( get_metadata( kVARS%ice1_c ))                   
        if (0<var_list( kVARS%ice2_mass))                   call this%add_to_output( get_metadata( kVARS%ice2_mass  ))   
        if (0<var_list( kVARS%ice2_number))                 call this%add_to_output( get_metadata( kVARS%ice2_number ))       
        if (0<var_list( kVARS%ice2_a))                      call this%add_to_output( get_metadata( kVARS%ice2_a ))                   
        if (0<var_list( kVARS%ice2_c))                      call this%add_to_output( get_metadata( kVARS%ice2_c ))                   
        if (0<var_list( kVARS%ice3_mass))                   call this%add_to_output( get_metadata( kVARS%ice3_mass  ))   
        if (0<var_list( kVARS%ice3_number))                 call this%add_to_output( get_metadata( kVARS%ice3_number ))        
        if (0<var_list( kVARS%ice3_a))                      call this%add_to_output( get_metadata( kVARS%ice3_a ))                   
        if (0<var_list( kVARS%ice3_c))                      call this%add_to_output( get_metadata( kVARS%ice3_c ))                   

        if (0<var_list( kVARS%precipitation) )              call this%add_to_output( get_metadata( kVARS%precipitation  ))
        if (0<var_list( kVARS%convective_precipitation) )   call this%add_to_output( get_metadata( kVARS%convective_precipitation  ))
        if (0<var_list( kVARS%snowfall) )                   call this%add_to_output( get_metadata( kVARS%snowfall   ))
        if (0<var_list( kVARS%graupel) )                    call this%add_to_output( get_metadata( kVARS%graupel         ))
        if (0<var_list( kVARS%pressure) )                   call this%add_to_output( get_metadata( kVARS%pressure        ))
        if (0<var_list( kVARS%temperature) )                call this%add_to_output( get_metadata( kVARS%temperature        ))
        if (0<var_list( kVARS%exner) )                      call this%add_to_output( get_metadata( kVARS%exner        ))
        if (0<var_list( kVARS%z) )                          call this%add_to_output( get_metadata( kVARS%z        ))
        if (0<var_list( kVARS%dz_interface) )               call this%add_to_output( get_metadata( kVARS%dz_interface        ))
        if (0<var_list( kVARS%z_interface) )                call this%add_to_output( get_metadata( kVARS%z_interface        ))
        if (0<var_list( kVARS%dz) )                         call this%add_to_output( get_metadata( kVARS%dz             ))
        if (0<var_list( kVARS%density) )                    call this%add_to_output( get_metadata( kVARS%density        ))
        if (0<var_list( kVARS%pressure_interface) )         call this%add_to_output( get_metadata( kVARS%pressure_interface ))
        if (0<var_list( kVARS%cloud_fraction) )             call this%add_to_output( get_metadata( kVARS%cloud_fraction   ))
        if (0<var_list( kVARS%shortwave) )                  call this%add_to_output( get_metadata( kVARS%shortwave        ))
        if (0<var_list( kVARS%shortwave_direct) )           call this%add_to_output( get_metadata( kVARS%shortwave_direct   ))
        if (0<var_list( kVARS%shortwave_diffuse) )          call this%add_to_output( get_metadata( kVARS%shortwave_diffuse   ))
        if (0<var_list( kVARS%longwave) )                   call this%add_to_output( get_metadata( kVARS%longwave        ))
        if (0<var_list( kVARS%vegetation_fraction) )        call this%add_to_output( get_metadata( kVARS%vegetation_fraction ))
        if (0<var_list( kVARS%vegetation_fraction_max) )    call this%add_to_output( get_metadata( kVARS%vegetation_fraction_max    ))
        if (0<var_list( kVARS%vegetation_fraction_out) )    call this%add_to_output( get_metadata( kVARS%vegetation_fraction_out    ))
        if (0<var_list( kVARS%lai) )                        call this%add_to_output( get_metadata( kVARS%lai        ))
        if (0<var_list( kVARS%sai) )                        call this%add_to_output( get_metadata( kVARS%sai        ))
        if (0<var_list( kVARS%crop_type) )                  call this%add_to_output( get_metadata( kVARS%crop_type   ))
        if (0<var_list( kVARS%date_planting) )              call this%add_to_output( get_metadata( kVARS%date_planting   ))
        if (0<var_list( kVARS%date_harvest) )               call this%add_to_output( get_metadata( kVARS%date_harvest     ))
        if (0<var_list( kVARS%growing_season_gdd) )         call this%add_to_output( get_metadata( kVARS%growing_season_gdd ))
        if (0<var_list( kVARS%irr_frac_total) )             call this%add_to_output( get_metadata( kVARS%irr_frac_total   ))
        if (0<var_list( kVARS%irr_frac_sprinkler) )         call this%add_to_output( get_metadata( kVARS%irr_frac_sprinkler ))
        if (0<var_list( kVARS%irr_frac_micro) )             call this%add_to_output( get_metadata( kVARS%irr_frac_micro   ))
        if (0<var_list( kVARS%irr_frac_flood) )             call this%add_to_output( get_metadata( kVARS%irr_frac_flood   ))
        if (0<var_list( kVARS%irr_alloc_sprinkler) )        call this%add_to_output( get_metadata( kVARS%irr_alloc_sprinkler    ))
        if (0<var_list( kVARS%irr_alloc_micro) )            call this%add_to_output( get_metadata( kVARS%irr_alloc_micro   ))
        if (0<var_list( kVARS%irr_alloc_flood) )            call this%add_to_output( get_metadata( kVARS%irr_alloc_flood   ))
        if (0<var_list( kVARS%irr_evap_loss_sprinkler) )    call this%add_to_output( get_metadata( kVARS%irr_evap_loss_sprinkler    ))
        if (0<var_list( kVARS%irr_amt_sprinkler) )          call this%add_to_output( get_metadata( kVARS%irr_amt_sprinkler   ))
        if (0<var_list( kVARS%irr_amt_micro) )              call this%add_to_output( get_metadata( kVARS%irr_amt_micro     ))
        if (0<var_list( kVARS%irr_amt_flood) )              call this%add_to_output( get_metadata( kVARS%irr_amt_flood     ))
        if (0<var_list( kVARS%evap_heat_sprinkler) )        call this%add_to_output( get_metadata( kVARS%evap_heat_sprinkler ))
        if (0<var_list( kVARS%mass_ag_grain) )              call this%add_to_output( get_metadata( kVARS%mass_ag_grain     ))
        if (0<var_list( kVARS%growing_degree_days) )        call this%add_to_output( get_metadata( kVARS%growing_degree_days ))
        if (0<var_list( kVARS%net_ecosystem_exchange) )     call this%add_to_output( get_metadata( kVARS%net_ecosystem_exchange    ))
        if (0<var_list( kVARS%gross_primary_prod) )         call this%add_to_output( get_metadata( kVARS%gross_primary_prod ))
        if (0<var_list( kVARS%net_primary_prod) )           call this%add_to_output( get_metadata( kVARS%net_primary_prod   ))
        if (0<var_list( kVARS%apar) )                       call this%add_to_output( get_metadata( kVARS%apar        ))
        if (0<var_list( kVARS%photosynthesis_total) )       call this%add_to_output( get_metadata( kVARS%photosynthesis_total    ))
        if (0<var_list( kVARS%stomatal_resist_total) )      call this%add_to_output( get_metadata( kVARS%stomatal_resist_total    ))
        if (0<var_list( kVARS%stomatal_resist_sun) )        call this%add_to_output( get_metadata( kVARS%stomatal_resist_sun    ))
        if (0<var_list( kVARS%stomatal_resist_shade) )      call this%add_to_output( get_metadata( kVARS%stomatal_resist_shade    ))
        if (0<var_list( kVARS%gecros_state) )               call this%add_to_output( get_metadata( kVARS%gecros_state     ))
        if (0<var_list( kVARS%canopy_water) )               call this%add_to_output( get_metadata( kVARS%canopy_water     ))
        if (0<var_list( kVARS%canopy_water_ice) )           call this%add_to_output( get_metadata( kVARS%canopy_water_ice   ))
        if (0<var_list( kVARS%canopy_water_liquid) )        call this%add_to_output( get_metadata( kVARS%canopy_water_liquid ))
        if (0<var_list( kVARS%canopy_vapor_pressure) )      call this%add_to_output( get_metadata( kVARS%canopy_vapor_pressure    ))
        if (0<var_list( kVARS%canopy_temperature) )         call this%add_to_output( get_metadata( kVARS%canopy_temperature ))
        if (0<var_list( kVARS%canopy_fwet) )                call this%add_to_output( get_metadata( kVARS%canopy_fwet     ))
        if (0<var_list( kVARS%veg_leaf_temperature) )       call this%add_to_output( get_metadata( kVARS%veg_leaf_temperature    ))
        if (0<var_list( kVARS%ground_surf_temperature) )    call this%add_to_output( get_metadata( kVARS%ground_surf_temperature    ))
        if (0<var_list( kVARS%frac_within_gap) )            call this%add_to_output( get_metadata( kVARS%frac_within_gap     ))
        if (0<var_list( kVARS%frac_between_gap) )           call this%add_to_output( get_metadata( kVARS%frac_between_gap   ))
        if (0<var_list( kVARS%ground_temperature_bare) )    call this%add_to_output( get_metadata( kVARS%ground_temperature_bare    ))
        if (0<var_list( kVARS%ground_temperature_canopy) )  call this%add_to_output( get_metadata( kVARS%ground_temperature_canopy    ))
        if (0<var_list( kVARS%snowfall_ground) )            call this%add_to_output( get_metadata( kVARS%snowfall_ground   ))
        if (0<var_list( kVARS%rainfall_ground) )            call this%add_to_output( get_metadata( kVARS%rainfall_ground   ))
        if (0<var_list( kVARS%snow_water_equivalent) )      call this%add_to_output( get_metadata( kVARS%snow_water_equivalent    ))
        if (0<var_list( kVARS%snow_water_eq_prev) )         call this%add_to_output( get_metadata( kVARS%snow_water_eq_prev ))
        if (0<var_list( kVARS%snow_albedo_prev) )           call this%add_to_output( get_metadata( kVARS%snow_albedo_prev   ))
        if (0<var_list( kVARS%snow_temperature) )           call this%add_to_output( get_metadata( kVARS%snow_temperature   ))
        if (0<var_list( kVARS%snow_layer_depth) )           call this%add_to_output( get_metadata( kVARS%snow_layer_depth   ))
        if (0<var_list( kVARS%snow_layer_ice) )             call this%add_to_output( get_metadata( kVARS%snow_layer_ice     ))
        if (0<var_list( kVARS%snow_layer_liquid_water) )    call this%add_to_output( get_metadata( kVARS%snow_layer_liquid_water    ))
        if (0<var_list( kVARS%snow_age_factor) )            call this%add_to_output( get_metadata( kVARS%snow_age_factor     ))
        if (0<var_list( kVARS%snow_height) )                call this%add_to_output( get_metadata( kVARS%snow_height   ))
        if (0<var_list( kVARS%skin_temperature) )           call this%add_to_output( get_metadata( kVARS%skin_temperature   ))
        if (0<var_list( kVARS%soil_water_content) )         call this%add_to_output( get_metadata( kVARS%soil_water_content ))
        if (0<var_list( kVARS%eq_soil_moisture) )           call this%add_to_output( get_metadata( kVARS%eq_soil_moisture   ))
        if (0<var_list( kVARS%smc_watertable_deep) )        call this%add_to_output( get_metadata( kVARS%smc_watertable_deep ))
        if (0<var_list( kVARS%recharge) )                   call this%add_to_output( get_metadata( kVARS%recharge      ))
        if (0<var_list( kVARS%recharge_deep) )              call this%add_to_output( get_metadata( kVARS%recharge_deep     ))
        if (0<var_list( kVARS%soil_temperature) )           call this%add_to_output( get_metadata( kVARS%soil_temperature   ))
        if (0<var_list( kVARS%latitude) )                   call this%add_to_output( get_metadata( kVARS%latitude      ))
        if (0<var_list( kVARS%longitude) )                  call this%add_to_output( get_metadata( kVARS%longitude      ))
        if (0<var_list( kVARS%u_latitude) )                 call this%add_to_output( get_metadata( kVARS%u_latitude      ))
        if (0<var_list( kVARS%u_longitude) )                call this%add_to_output( get_metadata( kVARS%u_longitude      ))
        if (0<var_list( kVARS%v_latitude) )                 call this%add_to_output( get_metadata( kVARS%v_latitude      ))
        if (0<var_list( kVARS%v_longitude) )                call this%add_to_output( get_metadata( kVARS%v_longitude      ))
        if (0<var_list( kVARS%terrain) )                    call this%add_to_output( get_metadata( kVARS%terrain      ))
        if (0<var_list( kVARS%sensible_heat) )              call this%add_to_output( get_metadata( kVARS%sensible_heat     ))
        if (0<var_list( kVARS%latent_heat) )                call this%add_to_output( get_metadata( kVARS%latent_heat      ))
        if (0<var_list( kVARS%u_10m) )                      call this%add_to_output( get_metadata( kVARS%u_10m        ))
        if (0<var_list( kVARS%v_10m) )                      call this%add_to_output( get_metadata( kVARS%v_10m        ))
        if (0<var_list( kVARS%coeff_momentum_drag) )        call this%add_to_output( get_metadata( kVARS%coeff_momentum_drag ))
        if (0<var_list( kVARS%coeff_heat_exchange) )        call this%add_to_output( get_metadata( kVARS%coeff_heat_exchange ))
        if (0<var_list( kVARS%surface_rad_temperature) )    call this%add_to_output( get_metadata( kVARS%surface_rad_temperature    ))
        if (0<var_list( kVARS%temperature_2m) )             call this%add_to_output( get_metadata( kVARS%temperature_2m     ))
        if (0<var_list( kVARS%humidity_2m) )                call this%add_to_output( get_metadata( kVARS%humidity_2m      ))
        if (0<var_list( kVARS%temperature_2m_veg) )         call this%add_to_output( get_metadata( kVARS%temperature_2m_veg ))
        if (0<var_list( kVARS%temperature_2m_bare) )        call this%add_to_output( get_metadata( kVARS%temperature_2m_bare ))
        if (0<var_list( kVARS%mixing_ratio_2m_veg) )        call this%add_to_output( get_metadata( kVARS%mixing_ratio_2m_veg ))
        if (0<var_list( kVARS%mixing_ratio_2m_bare) )       call this%add_to_output( get_metadata( kVARS%mixing_ratio_2m_bare ))
        if (0<var_list( kVARS%surface_pressure) )           call this%add_to_output( get_metadata( kVARS%surface_pressure   ))
        if (0<var_list( kVARS%rad_absorbed_total) )         call this%add_to_output( get_metadata( kVARS%rad_absorbed_total ))
        if (0<var_list( kVARS%rad_absorbed_veg) )           call this%add_to_output( get_metadata( kVARS%rad_absorbed_veg   ))
        if (0<var_list( kVARS%rad_absorbed_bare) )          call this%add_to_output( get_metadata( kVARS%rad_absorbed_bare   ))
        if (0<var_list( kVARS%rad_net_longwave) )           call this%add_to_output( get_metadata( kVARS%rad_net_longwave   ))
        if (0<var_list( kVARS%longwave_up) )                call this%add_to_output( get_metadata( kVARS%longwave_up      ))
        if (0<var_list( kVARS%ground_heat_flux) )           call this%add_to_output( get_metadata( kVARS%ground_heat_flux   ))
        if (0<var_list( kVARS%soil_deep_temperature) )      call this%add_to_output( get_metadata( kVARS%soil_deep_temperature    ))
        if (0<var_list( kVARS%evap_canopy) )                call this%add_to_output( get_metadata( kVARS%evap_canopy      ))
        if (0<var_list( kVARS%evap_soil_surface) )          call this%add_to_output( get_metadata( kVARS%evap_soil_surface   ))
        if (0<var_list( kVARS%transpiration_rate) )         call this%add_to_output( get_metadata( kVARS%transpiration_rate ))
        if (0<var_list( kVARS%ch_veg) )                     call this%add_to_output( get_metadata( kVARS%ch_veg      ))
        if (0<var_list( kVARS%ch_veg_2m) )                  call this%add_to_output( get_metadata( kVARS%ch_veg_2m      ))
        if (0<var_list( kVARS%ch_bare) )                    call this%add_to_output( get_metadata( kVARS%ch_bare      ))
        if (0<var_list( kVARS%ch_bare_2m) )                 call this%add_to_output( get_metadata( kVARS%ch_bare_2m      ))
        if (0<var_list( kVARS%ch_under_canopy) )            call this%add_to_output( get_metadata( kVARS%ch_under_canopy     ))
        if (0<var_list( kVARS%ch_leaf) )                    call this%add_to_output( get_metadata( kVARS%ch_leaf      ))
        if (0<var_list( kVARS%sensible_heat_veg) )          call this%add_to_output( get_metadata( kVARS%sensible_heat_veg   ))
        if (0<var_list( kVARS%sensible_heat_bare) )         call this%add_to_output( get_metadata( kVARS%sensible_heat_bare ))
        if (0<var_list( kVARS%sensible_heat_canopy) )       call this%add_to_output( get_metadata( kVARS%sensible_heat_canopy    ))
        if (0<var_list( kVARS%evap_heat_veg) )              call this%add_to_output( get_metadata( kVARS%evap_heat_veg     ))
        if (0<var_list( kVARS%evap_heat_bare) )             call this%add_to_output( get_metadata( kVARS%evap_heat_bare     ))
        if (0<var_list( kVARS%evap_heat_canopy) )           call this%add_to_output( get_metadata( kVARS%evap_heat_canopy   ))
        if (0<var_list( kVARS%transpiration_heat) )         call this%add_to_output( get_metadata( kVARS%transpiration_heat ))
        if (0<var_list( kVARS%ground_heat_veg) )            call this%add_to_output( get_metadata( kVARS%ground_heat_veg     ))
        if (0<var_list( kVARS%ground_heat_bare) )           call this%add_to_output( get_metadata( kVARS%ground_heat_bare     ))
        if (0<var_list( kVARS%net_longwave_veg) )           call this%add_to_output( get_metadata( kVARS%net_longwave_veg     ))
        if (0<var_list( kVARS%net_longwave_bare) )          call this%add_to_output( get_metadata( kVARS%net_longwave_bare   ))
        if (0<var_list( kVARS%net_longwave_canopy) )        call this%add_to_output( get_metadata( kVARS%net_longwave_canopy ))
        if (0<var_list( kVARS%runoff_surface) )             call this%add_to_output( get_metadata( kVARS%runoff_surface     ))
        if (0<var_list( kVARS%runoff_subsurface) )          call this%add_to_output( get_metadata( kVARS%runoff_subsurface   ))
        if (0<var_list( kVARS%soil_totalmoisture) )         call this%add_to_output( get_metadata( kVARS%soil_totalmoisture ))
        if (0<var_list( kVARS%water_table_depth) )          call this%add_to_output( get_metadata( kVARS%water_table_depth   ))
        if (0<var_list( kVARS%water_aquifer) )              call this%add_to_output( get_metadata( kVARS%water_aquifer     ))
        if (0<var_list( kVARS%storage_gw) )                 call this%add_to_output( get_metadata( kVARS%storage_gw      ))
        if (0<var_list( kVARS%storage_lake) )               call this%add_to_output( get_metadata( kVARS%storage_lake      ))
        if (0<var_list( kVARS%roughness_z0) )               call this%add_to_output( get_metadata( kVARS%roughness_z0      ))
        if (0<var_list( kVARS%mass_leaf) )                  call this%add_to_output( get_metadata( kVARS%mass_leaf      ))
        if (0<var_list( kVARS%mass_root) )                  call this%add_to_output( get_metadata( kVARS%mass_root      ))
        if (0<var_list( kVARS%mass_stem) )                  call this%add_to_output( get_metadata( kVARS%mass_stem      ))
        if (0<var_list( kVARS%mass_wood) )                  call this%add_to_output( get_metadata( kVARS%mass_wood      ))
        if (0<var_list( kVARS%soil_carbon_fast) )           call this%add_to_output( get_metadata( kVARS%soil_carbon_fast   ))
        if (0<var_list( kVARS%soil_carbon_stable) )         call this%add_to_output( get_metadata( kVARS%soil_carbon_stable ))
        if (0<var_list( kVARS%soil_texture_1) )             call this%add_to_output( get_metadata( kVARS%soil_texture_1     ))
        if (0<var_list( kVARS%soil_texture_2) )             call this%add_to_output( get_metadata( kVARS%soil_texture_2     ))
        if (0<var_list( kVARS%soil_texture_3) )             call this%add_to_output( get_metadata( kVARS%soil_texture_3     ))
        if (0<var_list( kVARS%soil_texture_4) )             call this%add_to_output( get_metadata( kVARS%soil_texture_4     ))
        if (0<var_list( kVARS%soil_sand_and_clay) )         call this%add_to_output( get_metadata( kVARS%soil_sand_and_clay ))
        if (0<var_list( kVARS%re_cloud) )                   call this%add_to_output( get_metadata( kVARS%re_cloud        ))
        if (0<var_list( kVARS%re_ice) )                     call this%add_to_output( get_metadata( kVARS%re_ice        ))
        if (0<var_list( kVARS%re_snow) )                    call this%add_to_output( get_metadata( kVARS%re_snow        ))
        if (0<var_list( kVARS%out_longwave_rad) )           call this%add_to_output( get_metadata( kVARS%out_longwave_rad     ))
        if (0<var_list( kVARS%longwave_cloud_forcing) )     call this%add_to_output( get_metadata( kVARS%longwave_cloud_forcing    ))
        if (0<var_list( kVARS%shortwave_cloud_forcing) )    call this%add_to_output( get_metadata( kVARS%shortwave_cloud_forcing    ))
        if (0<var_list( kVARS%cosine_zenith_angle) )        call this%add_to_output( get_metadata( kVARS%cosine_zenith_angle ))
        if (0<var_list( kVARS%land_emissivity) )            call this%add_to_output( get_metadata( kVARS%land_emissivity     ))
        if (0<var_list( kVARS%temperature_interface) )      call this%add_to_output( get_metadata( kVARS%temperature_interface ))
        if (0<var_list( kVARS%tend_swrad) )                 call this%add_to_output( get_metadata( kVARS%tend_swrad      ))

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
        do i=1,this%n_vars
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
        end associate

    end subroutine setup_time_variable


    subroutine save_data(this, current_step, time)
        implicit none
        class(output_t), intent(in) :: this
        integer,         intent(in) :: current_step
        type(Time_type), intent(in) :: time
        
        integer :: i, k_s, k_e

        integer :: start_three_D_t(4), start_three_D_t_b(4), start_three_D_t_b2(4)
        integer :: start_two_D_t(3), start_two_D_t_b(3), start_two_D_t_b2(3)
        integer :: cnt_3d(3), cnt_3d_b(3), cnt_3d_b2(3)
        integer :: cnt_2d(2)
        integer :: v_i_s, v_i_e, v_j_s, v_j_e
        integer :: v_i_s_b, v_i_e_b, v_j_s_b, v_j_e_b
        integer :: v_i_s_b2, v_i_e_b2, v_j_s_b2, v_j_e_b2


        do i=1,this%n_vars
            associate(var => this%variables(i))
                k_s = this%kts
                k_e = var%dim_len(3)
                
                start_three_D_t = (/ this%start_3d(1), this%start_3d(2), this%start_3d(3), current_step /)
                cnt_3d = (/ this%cnt_3d(1), this%cnt_3d(2), (k_e-k_s+1) /)
                    
                start_three_D_t_b = (/ this%start_3d_b(1), this%start_3d_b(2), this%start_3d_b(3), current_step /)
                cnt_3d_b = (/ this%cnt_3d_b(1), this%cnt_3d_b(2), (k_e-k_s+1) /)
                        
                start_three_D_t_b2 = (/ this%start_3d_b2(1), this%start_3d_b2(2), this%start_3d_b2(3), current_step /)
                cnt_3d_b2 = (/ this%cnt_3d_b2(1), this%cnt_3d_b2(2), (k_e-k_s+1) /)
                    
                start_two_D_t = (/ this%start_3d(1), this%start_3d(2), current_step /)
                start_two_D_t_b = (/ this%start_3d_b(1), this%start_3d_b(2), current_step /)
                start_two_D_t_b2 = (/ this%start_3d_b2(1), this%start_3d_b2(2), current_step /)

                if (this%ite == this%global_dim_len(1)) cnt_3d(1) = cnt_3d(1) + var%xstag
                if ((this%start_3d_b(1) - this%its + cnt_3d_b(1)) == this%global_dim_len(1)) cnt_3d_b(1) = cnt_3d_b(1) + var%xstag
                if ((this%start_3d_b2(1) - this%its + cnt_3d_b2(1)) == this%global_dim_len(1)) cnt_3d_b2(1) = cnt_3d_b2(1) + var%xstag
                if (this%jte == this%global_dim_len(2)) cnt_3d(2) = cnt_3d(2) + var%ystag
                if ((this%start_3d_b2(2) - this%jts + cnt_3d_b2(2)) == this%global_dim_len(2)) cnt_3d_b2(2) = cnt_3d_b2(2) + var%ystag

                v_i_s = this%start_3d(1) - this%its + 1
                v_i_e = v_i_s + cnt_3d(1) - 1
        
                v_j_s = this%start_3d(2) - this%jts + 1
                v_j_e = v_j_s + cnt_3d(2) - 1

                v_i_s_b = this%start_3d_b(1) - this%its + 1
                v_i_e_b = v_i_s_b + cnt_3d_b(1) - 1
        
                v_j_s_b = this%start_3d_b(2) - this%jts + 1
                v_j_e_b = v_j_s_b + cnt_3d_b(2) - 1
        
                v_i_s_b2 = this%start_3d_b2(1) - this%its + 1
                v_i_e_b2 = v_i_s_b2 + cnt_3d_b2(1) - 1
        
                v_j_s_b2 = this%start_3d_b2(2) - this%jts + 1
                v_j_e_b2 = v_j_s_b2 + cnt_3d_b2(2) - 1
                
                call check_ncdf( nf90_var_par_access(this%ncfile_id, var%var_id, nf90_collective))

                if (var%three_d) then
                    if (var%unlimited_dim) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id, &
                            reshape(var%data_3d(v_i_s:v_i_e,:,v_j_s:v_j_e), shape=cnt_3d, order=[1,3,2]), &
                                        start_three_D_t, count=(/cnt_3d(1), cnt_3d(2), cnt_3d(3), 1/)), "saving:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id, &
                            reshape(var%data_3d(v_i_s_b:v_i_e_b,:,v_j_s_b:v_j_e_b), shape=cnt_3d_b, order=[1,3,2]), &
                                        start_three_D_t_b, count=(/cnt_3d_b(1), cnt_3d_b(2), (k_e-k_s+1), 1/)), "saving:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id, &
                            reshape(var%data_3d(v_i_s_b2:v_i_e_b2,:,v_j_s_b2:v_j_e_b2), shape=cnt_3d_b2, order=[1,3,2]), &
                                        start_three_D_t_b2, count=(/cnt_3d_b2(1), cnt_3d_b2(2), (k_e-k_s+1), 1/)), "saving:"//trim(var%name) )
                    elseif (this%creating) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id, reshape(var%data_3d(v_i_s:v_i_e,:,v_j_s:v_j_e),  &
                            shape=cnt_3d, order=[1,3,2]), start=this%start_3d,&
                            count=cnt_3d ), "saving:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id, reshape(var%data_3d(v_i_s_b:v_i_e_b,:,v_j_s_b:v_j_e_b),  &
                            shape=cnt_3d_b, order=[1,3,2]), start=this%start_3d_b,&
                            count=cnt_3d_b ), "saving:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id, reshape(var%data_3d(v_i_s_b2:v_i_e_b2,:,v_j_s_b2:v_j_e_b2),  &
                            shape=cnt_3d_b2, order=[1,3,2]), start=this%start_3d_b2,&
                            count=cnt_3d_b2 ), "saving:"//trim(var%name) )
                    endif
                elseif (var%two_d) then
                    if (var%unlimited_dim) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d(v_i_s:v_i_e,v_j_s:v_j_e), &
                                start_two_D_t,count=(/ cnt_3d(1), cnt_3d(2), 1/)), "saving:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d(v_i_s_b:v_i_e_b,v_j_s_b:v_j_e_b), &
                                start_two_D_t_b,count=(/ cnt_3d_b(1), cnt_3d_b(2), 1/)), "saving:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d(v_i_s_b2:v_i_e_b2,v_j_s_b2:v_j_e_b2), &
                                start_two_D_t_b2,count=(/ cnt_3d_b2(1), cnt_3d_b2(2), 1/)), "saving:"//trim(var%name) )
                    elseif (this%creating) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d(v_i_s:v_i_e,v_j_s:v_j_e), &
                                    start=(/ this%start_3d(1), this%start_3d(2) /), &
                                    count=(/ cnt_3d(1), cnt_3d(2) /)), "saving:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d(v_i_s_b:v_i_e_b,v_j_s_b:v_j_e_b), &
                                    start=(/ this%start_3d_b(1), this%start_3d_b(2) /), &
                                    count=(/ cnt_3d_b(1), cnt_3d_b(2) /)), "saving:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d(v_i_s_b2:v_i_e_b2,v_j_s_b2:v_j_e_b2), &
                                    start=(/ this%start_3d_b2(1), this%start_3d_b2(2) /), &
                                    count=(/ cnt_3d_b2(1), cnt_3d_b2(2) /)), "saving:"//trim(var%name) )
                    endif
                endif
            end associate
        end do
        
        call check_ncdf( nf90_var_par_access(this%ncfile_id, this%time%var_id, nf90_collective))

        call check_ncdf( nf90_put_var(this%ncfile_id, this%time%var_id, dble(time%mjd()), [current_step]),   &
                   "saving:"//trim(this%time%name) )



    end subroutine save_data

    subroutine setup_dims_for_var(this, var)
        implicit none
        class(output_t),    intent(inout) :: this
        type(variable_t),   intent(inout) :: var
        integer :: i, err, dim_len
        

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
                    dim_len = this%global_dim_len(i)
                    if (i == 1) dim_len = dim_len+var%xstag
                    if (i == 2) dim_len = dim_len+var%ystag
                    if (i == 3) dim_len = var%dim_len(3)

                    call check_ncdf( nf90_def_dim(this%ncfile_id, var%dimensions(i), dim_len,       &
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

        err = nf90_inq_varid(this%ncfile_id, var%name, var%var_id)

        ! if the variable was not found in the netcdf file then we will define it.
        if (err /= NF90_NOERR) then
            
            call check_ncdf( nf90_def_var(this%ncfile_id, var%name, NF90_REAL, var%dim_ids, var%var_id), &
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
        

    end subroutine setup_variable

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

    module subroutine close_file(this)
        implicit none
        class(output_t),   intent(inout)  :: this

        if (this%ncfile_id > 0) then
            call check_ncdf(nf90_close(this%ncfile_id), "Closing file ")
            this%ncfile_id = -1
        endif

    end subroutine

end submodule
