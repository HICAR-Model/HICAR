
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
!!
!!----------------------------------------------------------
submodule(ioclient_interface) ioclient_implementation
  use debug_module,             only : check_ncdf
  use iso_fortran_env
  use output_metadata,          only : get_varindx


  implicit none

contains        
    
    module subroutine init(this, domain, forcing, options)
        implicit none
        class(ioclient_t),  intent(inout)  :: this
        type(domain_t),     intent(inout)  :: domain
        type(boundary_t),   intent(in)     :: forcing
        type(options_t),    intent(in)     :: options

        type(variable_t) :: var
        integer :: n, out_i, rst_i, var_indx

        if (this_image()==1) write(*,*) 'Initializing I/O Clients'
    
        this%i_s_r = forcing%its; this%i_e_r = forcing%ite
        this%k_s_r = forcing%kts; this%k_e_r = forcing%kte
        this%j_s_r = forcing%jts; this%j_e_r = forcing%jte
        
        
        this%i_s_w = domain%its; this%i_e_w = domain%ite
        this%k_s_w = domain%kts; this%k_e_w = domain%kte
        this%j_s_w = domain%jts; this%j_e_w = domain%jte
        if (domain%ims == domain%ids) this%i_s_w = domain%ids
        if (domain%ime == domain%ide) this%i_e_w = domain%ide !Add extra to accomodate staggered vars
        if (domain%jms == domain%jds) this%j_s_w = domain%jds
        if (domain%jme == domain%jde) this%j_e_w = domain%jde !Add extra to accomodate staggered vars

        !Setup arrays for information about accessing variables from write buffer
        allocate(this%out_var_indices(count(options%io_options%vars_for_output > 0)))
        allocate(this%rst_var_indices(count(options%vars_for_restart > 0)))
        allocate(this%rst_var_names(count(options%vars_for_restart > 0)))

        out_i = 1
        rst_i = 1

        associate(list => domain%vars_to_out)
        
        n = 1
        ! loop through the list of variables that need to be written out
        call list%reset_iterator()
        
        do while (list%has_more_elements())
            var = list%next()
            
            var_indx = get_varindx(var%name)
            if (options%io_options%vars_for_output(var_indx) > 0) then
                this%out_var_indices(out_i) = n
                out_i = out_i + 1
            endif
            if (options%vars_for_restart(var_indx) > 0) then
                this%rst_var_indices(rst_i) = n
                this%rst_var_names(rst_i) = var%name
                rst_i = rst_i + 1
            endif
            n = n+1
        enddo
        end associate
    end subroutine init

    ! This subroutine pushes the output fields from the domain object
    ! to the write buffer for the IO processes to use
    module subroutine push(this, domain, write_buffer)
        implicit none
        class(ioclient_t),   intent(inout) :: this
        type(domain_t),   intent(inout)    :: domain
        real, intent(inout), allocatable   :: write_buffer(:,:,:,:)[:]
        
        type(variable_t) :: var
        integer :: i, n, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w
                
        associate(list => domain%vars_to_out)
        
        n = 1

        ! loop through the list of variables that need to be written out
        call list%reset_iterator()
        
        do while (list%has_more_elements())
            ! get the next variable in the structure
            var = list%next()
            
            i_s_w = this%i_s_w; i_e_w = this%i_e_w
            j_s_w = this%j_s_w; j_e_w = this%j_e_w
            if (domain%ime == domain%ide) i_e_w = i_e_w+var%xstag !Add extra to accomodate staggered vars
            if (domain%jme == domain%jde) j_e_w = j_e_w+var%ystag !Add extra to accomodate staggered vars
            nx = i_e_w - i_s_w + 1
            ny = j_e_w - j_s_w + 1
            if (var%two_d) then
                write_buffer(n,1:nx,1,1:ny) = &
                        var%data_2d(i_s_w:i_e_w,j_s_w:j_e_w)
            else
                write_buffer(n,1:nx,1:var%dim_len(2),1:ny) = &
                        var%data_3d(i_s_w:i_e_w,1:var%dim_len(2),j_s_w:j_e_w)
            endif
            n = n+1
        enddo
        end associate  
                
    end subroutine 
    
    ! This subroutine receives the input fields from the IO buffer
    ! for assignment to the forcing object
    module subroutine receive(this, forcing, read_buffer)
        implicit none
        class(ioclient_t), intent(inout) :: this
        type(boundary_t), intent(inout)  :: forcing
        real, intent(in), allocatable    :: read_buffer(:,:,:,:)[:]

        type(variable_t)     :: var
        integer :: i, n, nx, ny
                
        associate(list => forcing%variables)
        n = 1
        nx = this%i_e_r - this%i_s_r + 1
        ny = this%j_e_r - this%j_s_r + 1

        ! loop through the list of variables that need to be read in
        call list%reset_iterator()
        
        !If the parent I/O server has not yet written all of our input vars, wait

        do while (list%has_more_elements())
            ! get the next variable in the structure
            var = list%next()
            if (var%computed) then
                cycle
            else
                    if (var%two_d) then
                        var%data_2d = read_buffer(n,1:nx,1,1:ny)
                    else
                        var%data_3d = read_buffer(n,1:nx,1:var%dim_len(2),1:ny)
                    endif
                    n = n+1
            endif
        enddo
        end associate
    
    end subroutine 
    
    ! This subroutine assigns the data from the write_buffer to the appropriate fields
    ! of the domain object. It is a "reverse write"
    module subroutine receive_rst(this, domain, write_buffer)
        implicit none
        class(ioclient_t), intent(inout) :: this
        type(domain_t),   intent(inout)  :: domain
        real, intent(in), allocatable    :: write_buffer(:,:,:,:)[:]

        type(variable_t)     :: var
        integer :: i, n, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w

        do i = 1,size(this%rst_var_indices)
            n = this%rst_var_indices(i)
            var = domain%vars_to_out%get_var(this%rst_var_names(i))

            i_s_w = this%i_s_w; i_e_w = this%i_e_w
            j_s_w = this%j_s_w; j_e_w = this%j_e_w
            i_e_w = i_e_w+var%xstag !Add extra to accomodate staggered vars
            j_e_w = j_e_w+var%ystag !Add extra to accomodate staggered vars
            nx = i_e_w - i_s_w + 1
            ny = j_e_w - j_s_w + 1
            if (var%two_d) then
                var%data_2d(i_s_w:i_e_w,j_s_w:j_e_w) = &
                     write_buffer(n,1:nx,1,1:ny)
            else
                var%data_3d(i_s_w:i_e_w,1:var%dim_len(2),j_s_w:j_e_w) = &
                    write_buffer(n,1:nx,1:var%dim_len(2),1:ny)
            endif
        enddo

    end subroutine 

    
end submodule
