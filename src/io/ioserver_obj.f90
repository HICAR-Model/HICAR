
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
submodule(ioserver_interface) ioserver_implementation
  use debug_module,             only : check_ncdf
  use iso_fortran_env
  use timer_interface,    only : timer_t

  implicit none

contains


    module subroutine init(this, domain, options, isrc, ierc, ksrc, kerc, jsrc, jerc, iswc, iewc, kswc, kewc, jswc, jewc)
        class(ioserver_t),   intent(inout)  :: this
        type(domain_t),   intent(in)        :: domain
        type(options_t), intent(in) :: options
        integer, allocatable, dimension(:), intent(in) :: isrc, ierc, ksrc, kerc, jsrc, jerc, iswc, iewc, kswc, kewc, jswc, jewc

        integer ::  n, globalComm, splitComm, color, ierr
        type(options_t) :: fake_options
        
        this%server_id = (this_image()/(num_images()/kNUM_SERVERS))
        if (this%server_id==1) write(*,*) 'Initializing I/O Server'

        this%outputer%base_file_name = options%io_options%output_file
        this%restart_count = options%io_options%restart_count
        this%rster%base_file_name = options%io_options%restart_out_file
        
        call decompose(this, isrc, ierc, ksrc, kerc, jsrc, jerc, iswc, iewc, kswc, kewc, jswc, jewc)

        call this%reader%init(this%i_s_r,this%i_e_r,this%k_s_r,this%k_e_r,this%j_s_r,this%j_e_r,options)
        this%n_children = size(this%children)
        this%n_r = this%reader%n_vars

        !call this%rster%set_domain(domain)
        !call this%rster%add_variables(options%vars_for_restart, domain)

        
        call this%outputer%init(domain,options,this%i_s_w,this%i_e_w,this%k_s_w,this%k_e_w,this%j_s_w,this%j_e_w)
        !write(*,*) 'outputter inited:   ', size(options%io_options%vars_for_output)

        this%n_w = this%outputer%n_vars

        !determine if we need to increase our k index due to some very large soil field
        do n = 1,this%n_w
            if(this%outputer%variables(n)%dim_len(3) > this%k_e_w) this%k_e_w = this%outputer%variables(n)%dim_len(3)
        enddo

        
        allocate(this%parent_write_buffer(this%n_w,this%i_s_w:this%i_e_w+1,this%k_s_w:this%k_e_w,this%j_s_w:this%j_e_w+1))
        
        do n = 1,this%n_w
            !if(this%outputer%variables(n)%name == options%io_options%vars_for_output) then
            if(this%outputer%variables(n)%three_d) then
                this%outputer%variables(n)%data_3d => this%parent_write_buffer(n,:,:,:)
            else
                this%outputer%variables(n)%data_2d => this%parent_write_buffer(n,:,1,:)
            endif
        enddo

        if (options%parameters%restart) then
            if (this%server_id==1) write(*,*) "Reading restart data"
                ! call restart_model(domain, this%rster, options)
        endif

    end subroutine
    


    subroutine decompose(this,i_s_r,i_e_r,k_s_r,k_e_r,j_s_r,j_e_r,i_s_w,i_e_w,k_s_w,k_e_w,j_s_w,j_e_w)
        class(ioserver_t),   intent(inout)  :: this
        integer, allocatable, dimension(:), intent(in) :: i_s_r,i_e_r,k_s_r,k_e_r,j_s_r,j_e_r,i_s_w,i_e_w,k_s_w,k_e_w,j_s_w,j_e_w
        logical, allocatable :: passed_children(:)
        integer :: n, k, i, col_slice, row_slice
        integer :: i_s, i_e, j_s, j_e
        integer :: xsplit, ysplit, xs, ys, nx, ny
        real    :: x, y, best, current
        integer, allocatable :: mins(:)
                        
        nx = maxval(i_e_w)
        ny = maxval(j_e_w)
                    
        xsplit = 1
        ysplit = kNUM_SERVERS
        xs = xsplit
        ys = ysplit

        x = nx/real(xsplit)
        y = ny/real(ysplit)

        if (y > x) then
            best = abs(1 - ( y / x ))
        else
            best = abs(1 - ( x / y ))
        endif
  
        !do i=kNUM_SERVERS,1,-1
        !    if (mod(kNUM_SERVERS,i)==0) then
        !        ysplit = i
        !        xsplit = kNUM_SERVERS / i

        !        x = (nx/float(xsplit))
        !        y = (ny/float(ysplit))

        !        if (y > x) then
        !            current = abs(1 - ( y / x ))
        !        else
        !            current = abs(1 - ( x / y ))
        !        endif

        !        if (current < best) then
        !            best = current
        !            xs = xsplit
        !            ys = ysplit
        !        endif
        !    endif
        !enddo
        !On exit, xs and ys have the number of servers in each direction
        
        row_slice = nint(ny*1.0/ys)
        j_s = 1 + row_slice*(ceiling(this%server_id*1.0/xs)-1)
        j_e = j_s + row_slice
        if (ceiling(this%server_id*1.0/xs)==ys) j_e = ny
    
        col_slice = nint(nx*1.0/xs)
        i_s = 1 + col_slice*mod((this%server_id-1),xs)
        i_e = i_s + col_slice
        if (mod(this%server_id,xs)==0) i_e = nx     
    
        !Clip extent of parent process to the closest child whose write bounds fit within
        allocate(mins(size(i_s_w)))
        allocate(passed_children(size(i_s_w)))

        mins = abs(i_s-i_s_w)
        i_s = i_s_w(minloc(mins,dim=1))
        mins = abs(i_e-i_e_w)
        i_e = i_e_w(minloc(mins,dim=1))

        mins = abs(j_s-j_s_w)
        j_s = j_s_w(minloc(mins,dim=1))
        mins = abs(j_e-j_e_w)
        j_e = j_e_w(minloc(mins,dim=1))
    
        !do n = 1,size(i_s_w)
        !    passed_children(n) = .False.
            
        !    if ( (i_s_w(n) >= i_s) .and. (i_e_w(n) <= i_e) .and. &
        !         (j_s_w(n) >= j_s) .and. (j_e_w(n) <= j_e)) then
        !        passed_children(n) = .True.
        !    endif
        !enddo
        passed_children = .False.
        passed_children(max((this_image()-(kNUM_PROC_PER_NODE-1)),1):(this_image()-1)) = .True.
        
        allocate(this%children(count(passed_children)))
        allocate(this%iswc(count(passed_children)))
        allocate(this%iewc(count(passed_children)))
        allocate(this%kswc(count(passed_children)))
        allocate(this%kewc(count(passed_children)))
        allocate(this%jswc(count(passed_children)))
        allocate(this%jewc(count(passed_children)))

        allocate(this%isrc(count(passed_children)))
        allocate(this%ierc(count(passed_children)))
        allocate(this%ksrc(count(passed_children)))
        allocate(this%kerc(count(passed_children)))
        allocate(this%jsrc(count(passed_children)))
        allocate(this%jerc(count(passed_children)))

        k = 1

        do n = 1,size(i_s_w)
            if (passed_children(n)) then
                this%children(k) = n
                this%iswc(k) = i_s_w(n)
                this%iewc(k) = i_e_w(n)
                this%kswc(k) = k_s_w(n)
                this%kewc(k) = k_e_w(n)
                this%jswc(k) = j_s_w(n)
                this%jewc(k) = j_e_w(n)
                
                this%isrc(k) = i_s_r(n)
                this%ierc(k) = i_e_r(n)
                this%ksrc(k) = k_s_r(n)
                this%kerc(k) = k_e_r(n)
                this%jsrc(k) = j_s_r(n)
                this%jerc(k) = j_e_r(n)

                k = k + 1
            endif
        enddo
        
        this%i_s_r = minval(this%isrc)
        this%i_e_r = maxval(this%ierc)
        this%i_s_w = minval(this%iswc)
        this%i_e_w = maxval(this%iewc)

        this%j_s_r = minval(this%jsrc)
        this%j_e_r = maxval(this%jerc)
        this%j_s_w = minval(this%jswc)
        this%j_e_w = maxval(this%jewc)
        
        this%k_s_r = minval(this%ksrc)
        this%k_e_r = maxval(this%kerc)
        this%k_s_w = minval(this%kswc)
        this%k_e_w = maxval(this%kewc)

    end subroutine decompose
    
    ! This subroutine has two functions, depending on the process which calls it
    ! For child processes, they send their domain output vars to the output coarray
    ! buffer of the parent process, and set a flag on themselves that they are ready for output
    ! For parent processes, they wait until all children report that they are ready for output
    ! Then, the parent process hands this coarray to the outputer object, which performs
    ! parallel I/O with the other parent processes.
    module subroutine write_file(this, domain, time, write_buffer)
        implicit none
        class(ioserver_t), intent(inout)  :: this
        type(domain_t),   intent(in)      :: domain
        type(Time_type),  intent(in)      :: time
        real, allocatable, intent(in)     :: write_buffer(:,:,:,:)[:]

        integer :: i, n, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w
        ! Loop through child images and send chunks of buffer array to each one
        
        this%parent_write_buffer = kEMPT_BUFF
        do i=1,size(this%children)
            n = this%children(i)
            
            i_s_w = this%iswc(i); i_e_w = this%iewc(i)
            j_s_w = this%jswc(i); j_e_w = this%jewc(i)
            if (domain%ide == i_e_w) i_e_w = i_e_w+1 !Add extra to accomodate staggered vars
            if (domain%jde == j_e_w) j_e_w = j_e_w+1 !Add extra to accomodate staggered vars
            nx = i_e_w - i_s_w + 1
            ny = j_e_w - j_s_w + 1

            this%parent_write_buffer(:,i_s_w:i_e_w,:,j_s_w:j_e_w) = &
                write_buffer(:,1:nx,:,1:ny)[n]
        enddo

        call this%outputer%save_file(time,this%IO_comms)        

    end subroutine 
    
    ! This subroutine calls the read file function from the input object
    ! and then passes the read-in data to the coarray fields on the 
    ! child forcing objects
    module subroutine read_file(this, forcing, read_buffer)
        class(ioserver_t), intent(inout) :: this
        type(boundary_t), intent(inout)  :: forcing
        real, intent(inout), allocatable :: read_buffer(:,:,:,:)[:]

        real, allocatable, dimension(:,:,:,:) :: parent_read_buffer
        integer :: i, n, nx, ny

        ! read file into buffer array
        call this%reader%read_next_step(parent_read_buffer,this%IO_comms)
        
        ! Loop through child images and send chunks of buffer array to each one
        do i=1,size(this%children)
            n = this%children(i)
            nx = this%ierc(i) - this%isrc(i) + 1
            ny = this%jerc(i) - this%jsrc(i) + 1
            read_buffer(:,1:nx,:,1:ny)[n] = &
                    parent_read_buffer(:,this%isrc(i):this%ierc(i),:,this%jsrc(i):this%jerc(i))
        enddo

    end subroutine 

    ! This function closes all open file handles. Files are left open by default
    ! to minimize I/O calls. When the program exits, this must be called
    module subroutine close_files(this)
        class(ioserver_t), intent(inout) :: this
        
        this%creating = .false.
        
        ! close files
        call this%reader%close_file()
        call this%outputer%close_file()
        call this%rster%close_file()

    end subroutine 

    
    
end submodule
