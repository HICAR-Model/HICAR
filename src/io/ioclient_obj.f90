
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


  implicit none

contains        
    
    module subroutine init(this, domain, forcing)
        implicit none
        class(ioclient_t),  intent(inout)  :: this
        type(domain_t),     intent(in)     :: domain
        type(boundary_t),   intent(in)     :: forcing


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

    
    end subroutine init

    ! This subroutine has two functions, depending on the process which calls it
    ! For child processes, they send their domain output vars to the output coarray
    ! buffer of the parent process, and set a flag on themselves that they are ready for output
    ! For parent processes, they wait until all children report that they are ready for output
    ! Then, the parent process hands this coarray to the outputer object, which performs
    ! parallel I/O with the other parent processes.
    module subroutine push(this, domain, write_buffer)
        implicit none
        class(ioclient_t),   intent(inout) :: this
        type(domain_t),   intent(inout)    :: domain
        real, intent(inout), allocatable   :: write_buffer(:,:,:,:)[:]
        
        type(variable_t) :: var
        integer :: i, n, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w
        
        !Send domain fields to parent ca_write_buffer
        associate(list => domain%vars_to_out)
        
        n = 1

        ! loop through the list of variables that need to be read in
        call list%reset_iterator()
        
        !If the parent I/O server has not yet written all of our input vars, wait

        do while (list%has_more_elements())
            ! get the next variable in the structure
            var = list%next()
            if (var%computed) then
                cycle
            else
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
                    write_buffer(n,1:nx,:,1:ny) = &
                            var%data_3d(i_s_w:i_e_w,:,j_s_w:j_e_w)
                endif
                n = n+1
            endif
        enddo
        end associate        

    end subroutine 
    
    ! This subroutine calls the read file function from the input object
    ! and then passes the read-in data to the coarray fields on the 
    ! child forcing objects
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
                        var%data_3d = read_buffer(n,1:nx,:,1:ny)
                    endif
                    n = n+1
            endif
        enddo
        end associate
    
    end subroutine 
    
end submodule
