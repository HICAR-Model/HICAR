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
        do i=1, size(options%io_options%vars_for_output)
            if (options%io_options%vars_for_output(i) > 0) then
                call this%add_to_output( get_metadata( i ))
            endif
        enddo

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
        integer :: i, i_s, i_e, k_s, k_e, j_s, j_e


        integer :: start_three_D_t(4)
        integer :: start_two_D_t(3)
        integer :: cnt_3d(3)
        integer :: cnt_2d(2)


        do i=1,this%n_vars
            associate(var => this%variables(i))
                i_s = this%its
                i_e = this%ite
                j_s = this%jts
                j_e = this%jte
        
                !if (i_e == this%global_dim_len(1)) i_e = i_e + var%xstag
                !if (j_e == this%global_dim_len(2)) j_e = j_e + var%ystag

                call check_ncdf( nf90_var_par_access(this%ncfile_id, var%var_id, nf90_collective))

                if (var%three_d) then
                    k_s = this%kts
                    k_e = var%dim_len(3)
                    
                    start_three_D_t = (/ i_s, j_s, k_s, current_step /)
                    cnt_3d = (/ (i_e-i_s+1), (j_e-j_s+1), (k_e-k_s+1)  /)
                    !write(*,*) start_three_D_t
                    !write(*,*) cnt_3d
                    !write(*,*) 'lbound data_3d, ax1:  ',lbound(var%data_3d,1)
                    !write(*,*) 'ubound data_3d, ax1:  ',ubound(var%data_3d,1)
                    !write(*,*) 'lbound data_3d, ax2:  ',lbound(var%data_3d,2)
                    !write(*,*) 'ubound data_3d, ax2:  ',ubound(var%data_3d,2)
                    !write(*,*) 'lbound data_3d, ax3:  ',lbound(var%data_3d,3)
                    !write(*,*) 'ubound data_3d, ax3:  ',ubound(var%data_3d,3)

                    if (var%unlimited_dim) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id, &
                            reshape(var%data_3d, shape=cnt_3d, order=[1,3,2]), &
                                        start_three_D_t, count=(/cnt_3d(1), cnt_3d(2), cnt_3d(3), 1/)), "saving:"//trim(var%name) )
                    elseif (this%creating) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id, &
                            reshape(var%data_3d,  &
                            shape=cnt_3d, order=[1,3,2]), start=(/ start_three_D_t(1), start_three_D_t(2), start_three_D_t(3) /),&
                            count=cnt_3d ), "saving:"//trim(var%name) )
                    endif

                elseif (var%two_d) then
                    start_two_D_t = (/ i_s, j_s, current_step /)
                    cnt_2d = (/ (i_e-i_s+1), (j_e-j_s+1) /)

                    if (var%unlimited_dim) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d, &
                                start_two_D_t,count=(/ cnt_2d(1), cnt_2d(2), 1/)), "saving:"//trim(var%name) )
                    elseif (this%creating) then
                        call check_ncdf( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d, &
                                    start=(/ start_two_D_t(1), start_two_D_t(2) /), &
                                    count=cnt_2d), "saving:"//trim(var%name) )
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
