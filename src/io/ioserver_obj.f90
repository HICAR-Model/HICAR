
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
  use output_metadata,          only : get_metadata, get_varindx

  implicit none

contains


    module subroutine init(this, domain, options, isrc, ierc, ksrc, kerc, jsrc, jerc, iswc, iewc, kswc, kewc, jswc, jewc)
        class(ioserver_t),  intent(inout)  :: this
        type(domain_t),     intent(inout)  :: domain
        type(options_t),    intent(in)     :: options
        integer, allocatable, dimension(:), intent(in) :: isrc, ierc, ksrc, kerc, jsrc, jerc, iswc, iewc, kswc, kewc, jswc, jewc
        integer ::  n, var_indx, out_i, rst_i
        
        
        this%server_id = (this_image()/(num_images()/kNUM_SERVERS))
        this%io_time = options%parameters%start_time
        if (this%server_id==1) write(*,*) 'Initializing I/O Server'
        this%ide = domain%ide
        this%kde = domain%kde
        this%jde = domain%jde

        call decompose(this, isrc, ierc, ksrc, kerc, jsrc, jerc, iswc, iewc, kswc, kewc, jswc, jewc)


        !Setup reading capability
        call this%reader%init(this%i_s_r,this%i_e_r,this%k_s_r,this%k_e_r,this%j_s_r,this%j_e_r,options)
        this%n_children = size(this%children)
        this%n_r = this%reader%n_vars


        !Setup writing capability
        call this%outputer%init(domain,options,this%i_s_w,this%i_e_w,this%k_s_w,this%k_e_w,this%j_s_w,this%j_e_w)
        
        this%n_w = this%outputer%n_vars

        !determine if we need to increase our k index due to some very large soil field
        do n = 1,this%n_w
            if(this%outputer%variables(n)%dim_len(3) > this%k_e_w) this%k_e_w = this%outputer%variables(n)%dim_len(3)
        enddo
        
        !Link local buffer to the outputer variables
        allocate(this%parent_write_buffer(this%n_w,this%i_s_w:this%i_e_w+1,this%k_s_w:this%k_e_w,this%j_s_w:this%j_e_w+1))
        
        do n = 1,this%n_w
            if(this%outputer%variables(n)%three_d) then
                this%outputer%variables(n)%data_3d => this%parent_write_buffer(n,:,:,:)
            else
                this%outputer%variables(n)%data_2d => this%parent_write_buffer(n,:,1,:)
            endif
        enddo

        !Setup arrays for information about accessing variables from write buffer
        allocate(this%out_var_indices(count(options%io_options%vars_for_output > 0)))
        allocate(this%rst_var_indices(count(options%vars_for_restart > 0)))
        allocate(this%rst_var_names(count(options%vars_for_restart > 0)))

        out_i = 1
        rst_i = 1
        
        do n=1,this%n_w
            var_indx = get_varindx(this%outputer%variables(n)%name)
            if (options%io_options%vars_for_output(var_indx) > 0) then
                this%out_var_indices(out_i) = n
                out_i = out_i + 1
            endif
            if (options%vars_for_restart(var_indx) > 0) then
                this%rst_var_indices(rst_i) = n
                this%rst_var_names(rst_i) = this%outputer%variables(n)%name
                rst_i = rst_i + 1
            endif
        enddo


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
    
    ! This subroutine gathers the write buffers of its children 
    ! compute processes and then writes them to the output file
    module subroutine write_file(this, time, write_buffer)
        implicit none
        class(ioserver_t), intent(inout)  :: this
        type(Time_type),  intent(in)      :: time
        real, allocatable, intent(in)     :: write_buffer(:,:,:,:)[:]

        integer :: i, n, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w
        ! Loop through child images and send chunks of buffer array to each one
        
        this%parent_write_buffer = kEMPT_BUFF
        do i=1,this%n_children
            n = this%children(i)
            
            i_s_w = this%iswc(i); i_e_w = this%iewc(i)
            j_s_w = this%jswc(i); j_e_w = this%jewc(i)
            !If a particular child process is at the boundaries
            if (this%ide == i_e_w) i_e_w = i_e_w+1 !Add extra to accomodate staggered vars
            if (this%jde == j_e_w) j_e_w = j_e_w+1 !Add extra to accomodate staggered vars
            nx = i_e_w - i_s_w + 1
            ny = j_e_w - j_s_w + 1

            this%parent_write_buffer(:,i_s_w:i_e_w,:,j_s_w:j_e_w) = &
                write_buffer(:,1:nx,:,1:ny)[n]
        enddo

        call this%outputer%save_out_file(time,this%IO_comms,this%out_var_indices,this%rst_var_indices)        

    end subroutine 

    
    ! This subroutine calls the read file function from the input object
    ! and then passes the read-in data to the read buffer
    module subroutine read_file(this, read_buffer)
        class(ioserver_t), intent(inout) :: this
        real, intent(inout), allocatable :: read_buffer(:,:,:,:)[:]

        real, allocatable, dimension(:,:,:,:) :: parent_read_buffer
        integer :: i, n, nx, ny

        ! read file into buffer array
        call this%reader%read_next_step(parent_read_buffer,this%IO_comms)
        
        ! Loop through child images and send chunks of buffer array to each one
        do i=1,this%n_children
            n = this%children(i)
            nx = this%ierc(i) - this%isrc(i) + 1
            ny = this%jerc(i) - this%jsrc(i) + 1
            read_buffer(:,1:nx,:,1:ny)[n] = &
                    parent_read_buffer(:,this%isrc(i):this%ierc(i),:,this%jsrc(i):this%jerc(i))
        enddo

    end subroutine 

    ! Same as above, but for restart file
    module subroutine read_restart_file(this, options, write_buffer)
        class(ioserver_t),   intent(inout) :: this
        type(options_t),     intent(in)    :: options
        real, intent(inout), allocatable   :: write_buffer(:,:,:,:)[:]

        integer :: i, n, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w
        integer :: ncid, var_id, dimid_3d(4), nz, err, varid, start_3d(4), cnt_3d(4), start_2d(3), cnt_2d(3)
        real, allocatable :: data3d(:,:,:,:)
        type(variable_t)  :: var
        character(len=kMAX_NAME_LENGTH) :: name
        
        err = nf90_open(options%io_options%restart_in_file, IOR(nf90_nowrite,NF90_NETCDF4), ncid, comm = this%IO_comms, info = MPI_INFO_NULL)
        
        ! setup start/count arrays accordingly
        start_3d = (/ this%i_s_w,this%j_s_w,this%k_s_w,options%io_options%restart_step_in_file /)
        start_2d = (/ this%i_s_w,this%j_s_w,options%io_options%restart_step_in_file /)
        cnt_3d = (/ (this%i_e_w-this%i_s_w+1),(this%j_e_w-this%j_s_w+1),(this%k_e_w-this%k_s_w+1),1 /)
        cnt_2d = (/ (this%i_e_w-this%i_s_w+1),(this%j_e_w-this%j_s_w+1),1 /)

        this%parent_write_buffer = kEMPT_BUFF

        do i = 1,size(this%rst_var_indices)
            n = this%rst_var_indices(i)
            name = this%rst_var_names(i)
            var = get_metadata(get_varindx(name))
            
            call check_ncdf( nf90_inq_varid(ncid, name, var_id), " Getting var ID for "//trim(name))
            call check_ncdf( nf90_var_par_access(ncid, var_id, nf90_collective))
            
            
            nx = cnt_3d(1) + var%xstag
            ny = cnt_3d(2) + var%ystag

            if (var%three_d) then
                ! Get length of z dim
                call check_ncdf( nf90_inquire_variable(ncid, var_id, dimids = dimid_3d), " Getting dim IDs for "//trim(name))
                call check_ncdf( nf90_inquire_dimension(ncid, dimid_3d(3), len = nz), " Getting z dim len for "//trim(name))
                
                if (allocated(data3d)) deallocate(data3d)
                allocate(data3d(nx,ny,nz,1))
                call check_ncdf( nf90_get_var(ncid, var_id, data3d, start=start_3d, count=(/ nx, ny, nz /)), " Getting 3D var "//trim(name))

                this%parent_write_buffer(n,this%i_s_w:this%i_e_w+var%xstag,1:nz,this%j_s_w:this%j_e_w+var%ystag) = &
                        reshape(data3d(:,:,:,1), shape=[nx,nz,ny], order=[1,3,2])
            else if (var%two_d) then
                call check_ncdf( nf90_get_var(ncid, var_id, this%parent_write_buffer(n,this%i_s_w:this%i_e_w+var%xstag,1,this%j_s_w:this%j_e_w+var%ystag), &
                        start=start_2d, count=(/ nx, ny /)), " Getting 2D "//trim(name))
            endif
        end do
        
        call check_ncdf(nf90_close(ncid), "Closing file "//trim(options%io_options%restart_in_file))
        
        
        ! Loop through child images and send chunks of buffer array to each one
        do i=1,this%n_children
            n = this%children(i)

            i_s_w = this%iswc(i); i_e_w = this%iewc(i)
            j_s_w = this%jswc(i); j_e_w = this%jewc(i)
            i_e_w = i_e_w+1 !Add extra to accomodate staggered vars
            j_e_w = j_e_w+1 !Add extra to accomodate staggered vars
            nx = i_e_w - i_s_w + 1
            ny = j_e_w - j_s_w + 1

            write_buffer(:,1:nx,:,1:ny)[n] = &
                    this%parent_write_buffer(:,i_s_w:i_e_w,:,j_s_w:j_e_w)
        enddo

    end subroutine 


    ! This function closes all open file handles. Files are left open by default
    ! to minimize I/O calls. When the program exits, this must be called
    module subroutine close_files(this)
        class(ioserver_t), intent(inout) :: this
        
        this%creating = .false.
        
        ! close files
        call this%reader%close_file()
        call this%outputer%close_files()

    end subroutine 
    
end submodule
