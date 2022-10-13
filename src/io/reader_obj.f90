
!>----------------------------------------------------------
!!  Define the interface for the output object
!!
!!  Output objects store all of the data and references to data necessary to write
!!  an output file.  This includes primarily internal netcdf related IDs.
!!  Output objects also store an array of variables to output.
!!  These variables maintain pointers to the data to be output as well as
!!  Metadata (e.g. dimension names, units, other attributes)
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!  Ethan Gutmann (gutmann@ucar.edu)
!! 
!!----------------------------------------------------------
submodule(reader_interface) reader_implementation
  use debug_module,           only : check_ncdf
  use variable_interface,     only : variable_t
  use timer_interface,    only : timer_t

  use time_io,                only : read_times, find_timestep_in_file
  implicit none

contains


    module subroutine init(this, its, ite, kts, kte, jts, jte, options)
        class(reader_t), intent(inout) :: this
        integer, intent(in) :: its, ite, kts, kte, jts, jte
        type(options_t), intent(in) :: options

        
        integer, allocatable                         :: var_dimensions(:)
        character(len=kMAX_NAME_LENGTH), allocatable :: vars_to_read(:)
        type(variable_t) :: new_variable
        integer :: i, dims(3)
        

        this%file_list = options%parameters%boundary_files
        this%time_var  = options%parameters%time_var
        this%ncfile_id = -1
        ! figure out while file and timestep contains the requested start_time
        call set_curfile_curstep(this, options%parameters%start_time, this%file_list, this%time_var)
        !Now setup dimensions of reader_object so we know what part of input file that we should read
        this%its = its; this%ite = ite
        this%kts = kts; this%kte = kte
        this%jts = jts; this%jte = jte
        dims = (/ (this%ite-this%its+1),(this%kte-this%kts+1),(this%jte-this%jts+1) /)
        
        ! the parameters option type can't contain allocatable arrays because it is a coarray
        ! so we need to allocate the vars_to_read and var_dimensions outside of the options type
        call setup_variable_lists(options%parameters%vars_to_read, options%parameters%dim_list, vars_to_read, var_dimensions)

        this%n_vars = size(vars_to_read)

        do i=1, this%n_vars
            if (var_dimensions(i)==2) then
                call new_variable%initialize( [dims(1),dims(3)] )
                call this%variables%add_var(vars_to_read(i), new_variable)
            elseif (var_dimensions(i)==3) then
                call new_variable%initialize( dims )
                call this%variables%add_var(vars_to_read(i), new_variable)
            endif
        end do
        
    end subroutine init

    !>------------------------------------------------------------
    !! Reads a new set of forcing data for the next time step
    !!
    !!------------------------------------------------------------
    module subroutine read_next_step(this, buffer, par_comms)
        class(reader_t), intent(inout) :: this
        real, allocatable, intent(inout) :: buffer(:,:,:,:)
        integer, intent(in)              :: par_comms

        real, allocatable :: data3d(:,:,:,:)
        type(variable_t)  :: var
        character(len=kMAX_NAME_LENGTH) :: name
        integer :: nx, ny, nz, err, varid, n, start_3d(4), cnt_3d(4), start_2d(3), cnt_2d(3)
        type(timer_t)   :: timer

        if (allocated(buffer)) deallocate(buffer)
        allocate(buffer(this%n_vars,this%its:this%ite,this%kts:this%kte,this%jts:this%jte))

        !See if we must open the file
        if (this%ncfile_id < 0) then
            err = nf90_open(this%file_list(this%curfile), IOR(nf90_nowrite,NF90_NETCDF4), this%ncfile_id, comm = par_comms, info = MPI_INFO_NULL)
        endif
        
        ! setup start/count arrays accordingly
        start_3d = (/ this%its,this%jts,this%kts,this%curstep /)
        start_2d = (/ this%its,this%jts,this%curstep /)
        cnt_3d = (/ (this%ite-this%its+1),(this%jte-this%jts+1),(this%kte-this%kts+1),1 /)
        cnt_2d = (/ (this%ite-this%its+1),(this%jte-this%jts+1),1 /)

        associate(list => this%variables)
            
        ! loop through the list of variables that need to be read in
        call list%reset_iterator()
        n = 1
        do while (list%has_more_elements())
            ! get the next variable in the structure
            var = list%next(name)
            if (var%var_id < 0) call check_ncdf( nf90_inq_varid(this%ncfile_id, name, var%var_id), " Getting var ID for "//trim(name))
            call check_ncdf( nf90_var_par_access(this%ncfile_id, var%var_id, nf90_collective))
            if (var%three_d) then
                nx = size(var%data_3d, 1)
                ny = size(var%data_3d, 3)
                nz = size(var%data_3d, 2)
                if (allocated(data3d)) deallocate(data3d)
                allocate(data3d(nx,ny,nz,1))
                call check_ncdf( nf90_get_var(this%ncfile_id, var%var_id, data3d, start=start_3d, count=cnt_3d), " Getting 3D var "//trim(name))

                buffer(n,this%its:this%ite,this%kts:this%kte,this%jts:this%jte) = reshape(data3d(:,:,:,1), shape=[nx,nz,ny], order=[1,3,2])
            else if (var%two_d) then
                call check_ncdf( nf90_get_var(this%ncfile_id, var%var_id, buffer(n,:,1,:), start=start_2d, count=cnt_2d), " Getting 2D "//trim(name))
            endif

            n = n+1
        end do
        
        end associate

        call update_forcing_step(this)

    end subroutine
    
    
    !>------------------------------------------------------------
    !! Update the curstep and curfile (increments curstep and curfile if necessary)
    !!
    !!------------------------------------------------------------
    subroutine update_forcing_step(this)
        implicit none
        type(reader_t),   intent(inout) :: this


        integer :: steps_in_file

        this%curstep = this%curstep + 1 ! this may be all we have to do most of the time

        ! check that we haven't stepped passed the end of the current file
        steps_in_file = get_n_timesteps(this, this%time_var, 0)

        if (steps_in_file < this%curstep) then
            ! if we have, use the next file
            this%curfile = this%curfile + 1
            ! and the first timestep in the next file
            this%curstep = 1
            
            ! close current file
            call check_ncdf(nf90_close(this%ncfile_id), "Closing file "//trim(this%file_list(this%curfile)))
            this%ncfile_id = -1


            ! if we have run out of input files, stop with an error message
            if (this%curfile > size(this%file_list)) then
                stop "Ran out of files to process while searching for matching time variable!"
            endif

        endif
    end subroutine

    !>------------------------------------------------------------
    !! Figure out how many time steps are in a file based on a specified variable
    !!
    !! By default assumes that the variable has three dimensions.  If not, var_space_dims must be set
    !!------------------------------------------------------------
    function get_n_timesteps(this, varname, var_space_dims) result(steps_in_file)
        implicit none
        type(reader_t),   intent(in) :: this
        character(len=*), intent(in) :: varname
        integer,          intent(in), optional :: var_space_dims
        integer :: steps_in_file, i

        integer :: space_dims, varid, ndims, numDims, dimlen, dims(100), dimIds(100)

        space_dims=3
        if (present(var_space_dims)) space_dims = var_space_dims

        ! Get the varid of the variable, based on its name.
        call check_ncdf(nf90_inq_varid(this%ncfile_id, varname, varid),varname)
        ! find the number of dimensions
        call check_ncdf(nf90_inquire_variable(this%ncfile_id, varid, ndims = numDims),varname)
        ! find the dimension IDs
        call check_ncdf(nf90_inquire_variable(this%ncfile_id, varid, dimids = dimIds(:numDims)),varname)
        dims(1)=numDims
        ! finally, find the length of each dimension
        do i=1,numDims
            call check_ncdf(nf90_inquire_dimension(this%ncfile_id, dimIds(i), len = dimlen))
            dims(i+1)=dimlen
        end do

        if (dims(1) == space_dims) then
            steps_in_file = 1
        else
            steps_in_file = dims(dims(1)+1)
        endif

    end function


    !>------------------------------------------------------------
    !! Set the boundary data structure to the correct time step / file in the list of files
    !!
    !! Reads the time_var from each file successively until it finds a timestep that matches time
    !!------------------------------------------------------------
    subroutine set_curfile_curstep(this, time, file_list, time_var)
        implicit none
        class(reader_t), intent(inout) :: this
        type(Time_type),    intent(in) :: time
        character(len=*),   intent(in) :: file_list(:)
        character(len=*),   intent(in) :: time_var


        integer :: error

        ! these are module variables that should be correctly set when the subroutine returns
        this%curstep = 1
        this%curfile = 0
        error = 1
        do while ( (error/=0) .and. (this%curfile < size(file_list)) )
            this%curfile = this%curfile + 1
            this%curstep = find_timestep_in_file(file_list(this%curfile), time_var, time, error=error)
        enddo

        if (error==1) then
            stop "Ran out of files to process while searching for matching time variable!"
        endif

    end subroutine

    !>------------------------------------------------------------
    !! Setup the vars_to_read and var_dimensions arrays given a master set of variables
    !!
    !! Count the number of variables specified, then allocate and store those variables in a list just their size.
    !! The master list will have all variables, but not all will be set
    !!------------------------------------------------------------
    subroutine setup_variable_lists(master_var_list, master_dim_list, vars_to_read, var_dimensions)
        implicit none
        character(len=kMAX_NAME_LENGTH), intent(in)                 :: master_var_list(:)
        integer,                         intent(in)                 :: master_dim_list(:)
        character(len=kMAX_NAME_LENGTH), intent(inout), allocatable :: vars_to_read(:)
        integer,                         intent(inout), allocatable :: var_dimensions(:)

        integer :: n_valid_vars
        integer :: i, curvar, err

        n_valid_vars = 0
        do i=1, size(master_var_list)
            if (trim(master_var_list(i)) /= '') then
                n_valid_vars = n_valid_vars + 1
            endif
        enddo

        allocate(vars_to_read(  n_valid_vars), stat=err)
        if (err /= 0) stop "vars_to_read: Allocation request denied"

        allocate(var_dimensions(  n_valid_vars), stat=err)
        if (err /= 0) stop "var_dimensions: Allocation request denied"

        curvar = 1
        do i=1, size(master_var_list)
            if (trim(master_var_list(i)) /= '') then
                vars_to_read(curvar) = master_var_list(i)
                var_dimensions(curvar) = master_dim_list(i)
                ! if (this_image()==1) print *, "in variable list: ", vars_to_read(curvar)
                curvar = curvar + 1
            endif
        enddo
    end subroutine

    module subroutine close_file(this)
        implicit none
        class(reader_t),   intent(inout)  :: this

        if (this%ncfile_id > 0) then
            call check_ncdf(nf90_close(this%ncfile_id), "Closing file ")
            this%ncfile_id = -1
        endif

    end subroutine

    
end submodule
