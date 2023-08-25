!>-----------------------------------------
!! Main Program
!!
!! Initialize options and memory in init_model
!! Read initial conditions in bc_init (from a restart file if requested)
!! initialize physics packages in init_physics (e.g. tiedke and thompson if used)
!! If this run is a restart run, then set start to the restart timestep
!!      in otherwords, ntimesteps is the number of BC updates from the beginning of the entire model
!!      run, not just from the begining of this restart run
!! calculate model time in seconds based on the time between BC updates (in_dt)
!! Calculate the next model output time from current model time + output time delta (out_dt)
!!
!! Finally, loop until ntimesteps are reached updating boundary conditions and stepping the model forward
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!-----------------------------------------
program icar
    use mpi
    use iso_fortran_env
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use output_interface,   only : output_t
    use time_step,          only : step, dt_reduce                ! Advance the model forward in time
    use initialization,     only : init_model, init_physics, init_model_state
    use timer_interface,    only : timer_t
    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use icar_constants
    use wind_iterative,     only : finalize_iter_winds
    use ioserver_interface, only : ioserver_t
    use ioclient_interface, only : ioclient_t

    use land_surface,               only : lsm_init

    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    type(boundary_t):: boundary, add_cond
    type(event_type)  :: written_ev[*], write_ev[*], read_ev[*], child_read_ev[*], end_ev[*]
    type(output_t)  :: restart_dataset
    type(output_t)  :: output_dataset
    type(ioserver_t)  :: ioserver
    type(ioclient_t)  :: ioclient
    real, allocatable :: write_buffer(:,:,:,:)[:], read_buffer(:,:,:,:)[:]

    type(timer_t)   :: initialization_timer, total_timer, input_timer, output_timer, physics_timer, wind_timer, mp_timer, adv_timer, rad_timer, lsm_timer, pbl_timer, exch_timer, send_timer, ret_timer, wait_timer, forcing_timer, diagnostic_timer, wind_bal_timer
    type(Time_type) :: next_output, next_input
    type(time_delta_t) :: small_time_delta
    
    integer :: i, ierr, exec_team
    integer :: sleep_cnt, ev_cnt, end_ev_cnt
    real :: t_val
    logical :: init_flag

    !Initialize MPI if needed
    init_flag = .False.
    call MPI_initialized(init_flag, ierr)
    if (.not.(init_flag)) then
        call MPI_INIT(ierr)
        init_flag = .True.
    endif

    !Determine split of processes which will become I/O servers and which will be compute tasks
    !Also sets constants for the program to keep track of this splitting
    call split_processes(exec_team)

    call small_time_delta%set(1)
    
    call total_timer%start()
    call initialization_timer%start()
    
    !-----------------------------------------
    !  Model Initialization
    !
    ! Reads config options and initializes domain and boundary conditions
    ! First the compute domain has to be established
    
    ! Model structure must be initialized first. 
    ! The domain dimensions will be used when establishing I/O client-server relations
    call init_model(options, domain, boundary)

    call init_IO(exec_team, domain, boundary, options, ioclient, ioserver, write_buffer, read_buffer)        

    select case(exec_team)
    case(kCOMPUTE_TEAM)
    
        call EVENT_QUERY(read_ev,ev_cnt)
        sleep_cnt = 0
        do while(ev_cnt == 0)
            call sleep(1)
            sleep_cnt = sleep_cnt + 1
            call EVENT_QUERY(read_ev,ev_cnt)
            if (sleep_cnt >= kTIMEOUT) stop 'Timeout on waiting for reading of initial conditions'
        enddo
        !Finally, do an event wait to decrement the event counter
        EVENT WAIT (read_ev)

        call ioclient%receive(boundary, read_buffer)
        EVENT POST (child_read_ev[ioclient%server])
        
        call boundary%update_computed_vars(options, update=options%parameters%time_varying_z)

        call init_model_state(options, domain, boundary, add_cond) ! added boundary structure for external files (additional conditions)

        if (options%parameters%restart) then
            sync all !Matching IO sync call
            call EVENT_QUERY(read_ev,ev_cnt)
            sleep_cnt = 0
            do while(ev_cnt == 0)
                call sleep(1)
                sleep_cnt = sleep_cnt + 1
                call EVENT_QUERY(read_ev,ev_cnt)
                if (sleep_cnt >= kTIMEOUT) stop 'Timeout on waiting for reading of restart data'
            enddo
            !Finally, do an event wait to decrement the event counter
            EVENT WAIT (read_ev)

            if (this_image()==1) write(*,*) "Reading restart data"
            call ioclient%receive_rst(domain, write_buffer)
        endif
        
        ! physics drivers need to be initialized after restart data are potentially read in.
        call init_physics(options, domain, boundary)

        call output_timer%start()
        call ioclient%push(domain, write_buffer)
        EVENT POST (write_ev[ioclient%server])
        call output_timer%stop()
        next_output = options%parameters%start_time + options%io_options%output_dt
    case(kIO_TEAM)
    
        call ioserver%read_file(read_buffer)
        do i = 1,ioserver%n_children
            EVENT POST (read_ev[ioserver%children(i)])
        enddo
        
        if (options%parameters%restart) then
            sync all !Make sure we don't overwrite the file previously read in before the compute task has used it
            call ioserver%read_restart_file(options, write_buffer)
            do i = 1,ioserver%n_children
                EVENT POST (read_ev[ioserver%children(i)])
            enddo
        endif
                
        !We do not do output here so that we can continue straight to reading the next input file and not block the compute processes
        next_output = options%parameters%start_time
    end select
    
    call initialization_timer%stop()
    next_input = options%parameters%start_time + options%io_options%input_dt


    select case(exec_team)
    case(kCOMPUTE_TEAM)

        if (this_image()==1) write(*,*) "Initialization complete, beginning physics integration."
        do while (domain%model_time + small_time_delta < options%parameters%end_time)

            ! -----------------------------------------------------
            !
            !  Read input data if necessary
            !
            ! -----------------------------------------------------
            if ((domain%model_time + small_time_delta + options%io_options%input_dt) >= next_input) then
                if (this_image()==1) write(*,*) ""
                if (this_image()==1) write(*,*) " ----------------------------------------------------------------------"
                if (this_image()==1) write(*,*) "Updating Boundary conditions"
                
                call input_timer%start()

                call EVENT_QUERY(read_ev,ev_cnt)
                sleep_cnt = 0
                do while(ev_cnt == 0)
                    call sleep(1)
                    sleep_cnt = sleep_cnt + 1
                    call EVENT_QUERY(read_ev,ev_cnt)
                    if (sleep_cnt >= kTIMEOUT) stop 'Timeout on waiting for reading of input'
                enddo
                !Finally, do an event wait to decrement the event counter
                EVENT WAIT (read_ev)

                call ioclient%receive(boundary, read_buffer)
                EVENT POST (child_read_ev[ioclient%server])
                
                ! after reading all variables that can be read, not compute any remaining variables (e.g. z from p+ps)
                call boundary%update_computed_vars(options, update=.True.)

                ! if the vertical levels of the forcing data change over time, they need to be interpolated to the original levels here.
                if (options%parameters%time_varying_z) then
                    call boundary%interpolate_original_levels(options)
                endif
                call domain%interpolate_forcing(boundary, update=.True.)
                
                ! Make the boundary condition dXdt values into units of [X]/s
                call boundary%update_delta_fields(next_input - domain%model_time)
                call domain%update_delta_fields(next_input - domain%model_time)
                
                next_input = next_input + options%io_options%input_dt
                call input_timer%stop()
            endif


            ! -----------------------------------------------------
            !
            !  Integrate physics forward in time
            !
            ! -----------------------------------------------------
            if (this_image()==1) write(*,*) "Running Physics"
            if (this_image()==1) write(*,*) "  Model time = ", trim(domain%model_time%as_string())
            if (this_image()==1) write(*,*) "   End  time = ", trim(options%parameters%end_time%as_string())
            if (this_image()==1) write(*,*) "  Next Input = ", trim(next_input%as_string())
            if (this_image()==1) write(*,*) "  Next Output= ", trim(next_output%as_string())

            ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
            if (.not.(options%wind%wind_only)) then
                call physics_timer%start()
                call step(domain, boundary, step_end(next_input, next_output), options,           &
                                mp_timer, adv_timer, rad_timer, lsm_timer, pbl_timer, exch_timer, &
                                send_timer, ret_timer, wait_timer, forcing_timer, diagnostic_timer, wind_bal_timer, wind_timer)
                call physics_timer%stop()
            elseif (options%wind%wind_only) then
                call domain%apply_forcing(boundary, options, real(options%io_options%output_dt%seconds()))
                domain%model_time = next_output
            endif

            ! -----------------------------------------------------
            !  Write output data if it is time
            ! -----------------------------------------------------
            if ((domain%model_time + small_time_delta) >= next_output) then
                if (this_image()==1) write(*,*) "Writing output file"
                call output_timer%start()
                call EVENT_QUERY(written_ev,ev_cnt)
                sleep_cnt = 0
                do while(ev_cnt == 0)
                    call sleep(1)
                    sleep_cnt = sleep_cnt + 1
                    call EVENT_QUERY(written_ev,ev_cnt)
                    if (sleep_cnt >= kTIMEOUT) stop 'Timeout on waiting for previous output file to be written'
                enddo
                !Finally, do an event wait to decrement the event counter
                EVENT WAIT (written_ev)

                call ioclient%push(domain, write_buffer)
                EVENT POST (write_ev[ioclient%server])
                next_output = next_output + options%io_options%output_dt
                call output_timer%stop()
            endif
        end do
        
        !Wait for final write before ending program
        call EVENT_QUERY(written_ev,ev_cnt)
        sleep_cnt = 0
        do while(ev_cnt == 0)
            call sleep(1)
            sleep_cnt = sleep_cnt + 1
            call EVENT_QUERY(written_ev,ev_cnt)
            if (sleep_cnt >= kTIMEOUT) stop 'Timeout on waiting for last output file to be written'
        enddo
        !Finally, do an event wait to decrement the event counter
        EVENT WAIT (written_ev)
        
        EVENT POST (end_ev[ioclient%server])

        if (options%physics%windtype==kITERATIVE_WINDS) call finalize_iter_winds() 
    case (kIO_TEAM)
    
        call EVENT_QUERY(end_ev,end_ev_cnt)

        do while (end_ev_cnt < ioserver%n_children)
            !See of it is time to read
            call EVENT_QUERY(child_read_ev,ev_cnt)
            if (ev_cnt == ioserver%n_children .and. next_input<=options%parameters%end_time) then
                !Do an event wait to decrement the event counter
                EVENT WAIT (child_read_ev, UNTIL_COUNT=ioserver%n_children)
            
                call ioserver%read_file(read_buffer)
                do i = 1,ioserver%n_children
                    EVENT POST (read_ev[ioserver%children(i)])
                enddo
                next_input = next_input + options%io_options%input_dt
            endif

            
            !See if all of our children are ready for a write
            call EVENT_QUERY(write_ev,ev_cnt)
            if (ev_cnt == ioserver%n_children .and. next_output<=options%parameters%end_time) then
                !Do an event wait to decrement the event counter
                EVENT WAIT (write_ev, UNTIL_COUNT=ioserver%n_children)

                call ioserver%write_file(next_output, write_buffer)
                do i = 1,ioserver%n_children
                    EVENT POST (written_ev[ioserver%children(i)])
                enddo 
                next_output = next_output + options%io_options%output_dt
            endif
            
            !See if it is time to end
            call EVENT_QUERY(end_ev,end_ev_cnt)

        enddo
        !If we are done with the program
        call ioserver%close_files()
    end select
    !
    !-----------------------------------------
    call total_timer%stop()
    if (this_image()==1) then
        write(*,*) ""
        write(*,*) "Model run from : ",trim(options%parameters%start_time%as_string())
        write(*,*) "           to  : ",trim(options%parameters%end_time%as_string())
        write(*,*) "Domain : ",trim(options%parameters%init_conditions_file)
        write(*,*) "Number of images:",num_images()
        write(*,*) ""
        write(*,*) "Average timing across compute images:"
    endif
    t_val = timer_mean(total_timer)
    if (this_image()==1) write(*,*) "total          : ", t_val 
    t_val = timer_mean(initialization_timer)
    if (this_image()==1) write(*,*) "init           : ", t_val
    t_val = timer_mean(input_timer)
    if (this_image()==1) write(*,*) "input          : ", t_val
    t_val = timer_mean(output_timer)
    if (this_image()==1) write(*,*) "output         : ", t_val
    t_val = timer_mean(physics_timer)
    if (this_image()==1) write(*,*) "physics        : ", t_val
    t_val = timer_mean(mp_timer)
    if (this_image()==1) write(*,*) "microphysics   : ", t_val
    t_val = timer_mean(adv_timer)
    if (this_image()==1) write(*,*) "advection      : ", t_val
    t_val = timer_mean(rad_timer)
    if (this_image()==1) write(*,*) "radiation      : ", t_val
    t_val = timer_mean(lsm_timer)
    if (this_image()==1) write(*,*) "LSM            : ", t_val
    t_val = timer_mean(pbl_timer)
    if (this_image()==1) write(*,*) "PBL            : ", t_val
    t_val = timer_mean(forcing_timer)
    if (this_image()==1) write(*,*) "forcing        : ", t_val 
    t_val = timer_mean(wind_bal_timer)
    if (this_image()==1) write(*,*) "wind bal       : ", t_val
    t_val = timer_mean(diagnostic_timer)
    if (this_image()==1) write(*,*) "diagnostic     : ", t_val
    t_val = timer_mean(exch_timer)
    if (this_image()==1) write(*,*) "halo-exchange  : ", t_val 
    t_val = timer_mean(wind_timer)
    if (this_image()==1) write(*,*) "winds          : ", t_val 

contains

    !subroutine safe_wait(event, err_msg, until_count)
    !    implicit none
    !    type(event_type), intent(inout)  :: event[*]
    !    character(len=*), intent(in)     :: err_msg
    !    integer, intent(in), optional    :: until_count
    !    integer :: sleep_cnt, ev_cnt, until_cnt
    !        timer_mean
    !    call EVENT_QUERY(event,ev_cnt)
    !    sleep_cnt = 0
    !    do while(.not.(ev_cnt == 1))
    !        call sleep(1)
    !        sleep_cnt = sleep_cnt + 1
    !        call EVENT_QUERY(event,ev_cnt)
    !        if (sleep_cnt >= kTIMEOUT) stop trim(err_msg)
    !    enddo
    !    !Finally, do an event wait to decrement the event counter
    !    EVENT WAIT (event)!, UNTIL_COUNT=ioserver%n_children)
    !end subroutine safe_wait

    function timer_mean(timer) result(mean_t)
        implicit none
        type(timer_t), intent(inout) :: timer
        
        real :: mean_t, t_sum
            
        t_sum = timer%get_time()
        call CO_SUM(t_sum)
        mean_t = t_sum/kNUM_COMPUTE
    
    end function

    function step_end(time1, time2) result(min_time)
        implicit none
        type(Time_type), intent(in) :: time1
        type(Time_type), intent(in) :: time2
        type(Time_type) :: min_time

        if (time1 <= time2 ) then
            min_time = time1
        else
            min_time = time2
        endif

    end function


    subroutine split_processes(exec_team)
        implicit none
        integer, intent(inout) :: exec_team
        integer :: n, k, name_len, ierr, node_name_i, node_names(num_images())
        character(len=MPI_MAX_PROCESSOR_NAME) :: node_name
        
        node_names = 0
        node_name_i = 0
        
        call MPI_Get_processor_name(node_name, name_len, ierr)
        do n = 1,name_len
            node_name_i = node_name_i + ichar(node_name(n:n))*n*10
        enddo
        node_names(this_image()) = node_name_i
        call co_max(node_names)
        
        kNUM_PROC_PER_NODE = count(node_names==node_names(1))

        !Assign one io process per node, this results in best co-array transfer times
        kNUM_SERVERS = ceiling(num_images()*1.0/kNUM_PROC_PER_NODE)
        kNUM_COMPUTE = num_images()-kNUM_SERVERS
        
        if ((mod(kNUM_COMPUTE,2) /= 0) .and. this_image()==1) then
            write(*,*) 'WARNING: number of compute processes is odd-numbered.' 
            write(*,*) 'One process per node is used for I/O.'
            write(*,*) 'If the total number of compute processes is odd-numbered,'
            write(*,*) 'this may lead to errors with domain decomposition'
        endif
        
        if (mod(this_image(),(num_images()/kNUM_SERVERS)) == 0) then
            exec_team = kIO_TEAM
        else
            exec_team = kCOMPUTE_TEAM
        endif
        
        k = 1
        allocate(DOM_IMG_INDX(kNUM_COMPUTE))
        do n = 1,kNUM_COMPUTE
            if (mod(k,(num_images()/kNUM_SERVERS)) == 0) k = k+1
            DOM_IMG_INDX(n) = k
            k = k+1
        enddo
        !form team(team_num,exec_team)
        
    end subroutine split_processes
    
    subroutine init_IO(exec_team, domain, boundary, options, ioclient, ioserver, write_buffer, read_buffer)
        implicit none
        integer, intent(in) :: exec_team
        type(domain_t),  intent(inout)  :: domain
        type(boundary_t),intent(in)  :: boundary
        type(options_t), intent(in)  :: options
        type(ioclient_t), intent(inout) :: ioclient
        type(ioserver_t), intent(inout) :: ioserver
        real, intent(inout), allocatable :: write_buffer(:,:,:,:)[:], read_buffer(:,:,:,:)[:]

        integer, allocatable, dimension(:) :: i_s_w, i_e_w, k_s_w, k_e_w, j_s_w, j_e_w, i_s_r, i_e_r, k_s_r, k_e_r, j_s_r, j_e_r
        integer, allocatable, dimension(:,:) :: childrens
        integer :: nx_w, ny_w, nz_w, n_w, nx_r, ny_r, nz_r, n_r, n_restart, i, globalComm, splitComm, color, ierr
        integer :: out_i, rst_i, var_indx
        
        allocate(i_s_w(num_images())); allocate(i_e_w(num_images()))
        allocate(i_s_r(num_images())); allocate(i_e_r(num_images()))
        
        allocate(j_s_w(num_images())); allocate(j_e_w(num_images()))
        allocate(j_s_r(num_images())); allocate(j_e_r(num_images()))
        
        allocate(k_s_w(num_images())); allocate(k_e_w(num_images()))
        allocate(k_s_r(num_images())); allocate(k_e_r(num_images()))

        allocate(childrens(kNUM_SERVERS,num_images()))

        i_s_w = 0; i_e_w = 0; i_s_r = 0; i_e_r = 0
        j_s_w = 0; j_e_w = 0; j_s_r = 0; j_e_r = 0
        k_s_w = 0; k_e_w = 0; k_s_r = 0; k_e_r = 0
        
        nx_w = 0; ny_w = 0; nz_w = 0; nx_r = 0; ny_r = 0; nz_r = 0;
        n_w = 0; n_r = 0; n_restart = 0;
        
        childrens = 0
        
        select case (exec_team)
        case (kCOMPUTE_TEAM)
            color = 0
        case(kIO_TEAM)
            color = 1
        end select
        
        CALL MPI_COMM_DUP( MPI_COMM_WORLD, globalComm, ierr )
        ! Assign all images in the IO team to the IO_comms MPI communicator. Use image indexing within initial team to get indexing of global MPI ranks
        CALL MPI_COMM_SPLIT( globalComm, color, (this_image()-1), splitComm, ierr )
        
        select case (exec_team)
        case (kCOMPUTE_TEAM)
            CALL MPI_COMM_DUP( splitComm, domain%IO_comms, ierr )
            call ioclient%init(domain, boundary, options)
            
            i_s_w(this_image()) = ioclient%i_s_w; i_e_w(this_image()) = ioclient%i_e_w
            i_s_r(this_image()) = ioclient%i_s_r; i_e_r(this_image()) = ioclient%i_e_r
            
            j_s_w(this_image()) = ioclient%j_s_w; j_e_w(this_image()) = ioclient%j_e_w
            j_s_r(this_image()) = ioclient%j_s_r; j_e_r(this_image()) = ioclient%j_e_r
            
            k_s_w(this_image()) = ioclient%k_s_w; k_e_w(this_image()) = ioclient%k_e_w
            k_s_r(this_image()) = ioclient%k_s_r; k_e_r(this_image()) = ioclient%k_e_r
        end select
                
        call co_max(i_s_w); call co_max(i_e_w); call co_max(i_s_r); call co_max(i_e_r)
        call co_max(j_s_w); call co_max(j_e_w); call co_max(j_s_r); call co_max(j_e_r)
        call co_max(k_s_w); call co_max(k_e_w); call co_max(k_s_r); call co_max(k_e_r)            
        
    
        select case(exec_team)
        case(kIO_TEAM)
            CALL MPI_COMM_DUP( splitComm, ioserver%IO_comms, ierr )
            call ioserver%init(domain, options,i_s_r,i_e_r,k_s_r,k_e_r,j_s_r,j_e_r,i_s_w,i_e_w,k_s_w,k_e_w,j_s_w,j_e_w)
            
            if (ioserver%server_id==1) then
                write(*,*) "Setting up output files"
                if (options%io_options%frames_per_outfile<2) then
                    print*,"  frames per output file should be 2 or more. Currently: ", options%io_options%frames_per_outfile
                else
                    print*,"  frames per output file= ", options%io_options%frames_per_outfile
                end if
            endif

            !Contribute k-extents for wrie variables, in case we had to expand beyond the domain k extent due to a large soil variable
            if ( ioserver%k_e_w > maxval(k_e_w)) k_e_w(1) = ioserver%k_e_w

            n_r = ioserver%n_r
            n_w = ioserver%n_w
            n_restart = ioserver%n_restart
            childrens(ioserver%server_id,1:ioserver%n_children) = ioserver%children
        end select

        call co_max(childrens)
        
        select case(exec_team)
        case(kCOMPUTE_TEAM)
            do i = 1,size(childrens,1)
                if (findloc(childrens(i,:),this_image(),dim=1) > 0) then
                    ioclient%server = i*(num_images()/kNUM_SERVERS)
                    exit
                endif
            enddo
            
            !Buffer arrays must be initialized in the global image team so both teams have access
            nx_w = i_e_w(this_image()) - i_s_w(this_image()) + 1
            ny_w = j_e_w(this_image()) - j_s_w(this_image()) + 1
            nz_w = k_e_w(this_image()) - k_s_w(this_image()) + 1

            nx_r = i_e_r(this_image()) - i_s_r(this_image()) + 1
            ny_r = j_e_r(this_image()) - j_s_r(this_image()) + 1
            nz_r = k_e_r(this_image()) - k_s_r(this_image()) + 1
        end select

        call co_max(nx_w)
        call co_max(ny_w)
        call co_max(nz_w)

        call co_max(nx_r)
        call co_max(ny_r)
        call co_max(nz_r)

        call co_max(n_w)
        call co_max(n_r)
        call co_max(n_restart)
        
        
        !Initialize read/write buffers. They are coarrays, defined by the comax of the read and write bounds for each compute image
        allocate(write_buffer(n_w,1:nx_w+1,1:nz_w,1:ny_w+1)[*])
        allocate(read_buffer(n_r,1:nx_r+1,1:nz_r,1:ny_r+1)[*])

    end subroutine init_IO

    
end program

! This is the Doxygen mainpage documentation.  This should be moved to another file at some point.

!>------------------------------------------
!!  @mainpage
!!
!!  @section Introduction
!!  ICAR is a simplified atmospheric model designed primarily for climate downscaling, atmospheric sensitivity tests,
!!  and hopefully educational uses. At this early stage, the model is still undergoing rapid development, and users
!!  are encouraged to get updates frequently.
!!
!!  @section Running_ICAR
!!  To run the model 3D time-varying atmospheric data are required, though an ideal test case can be generated for
!!  simple simulations as well. There are some sample python scripts to help make input forcing files, but the WRF
!!  pre-processing system can also be used. Low-resolution WRF output files can be used directly, various reanalysis
!!  and GCM output files can be used with minimal pre-processing (just get all the variables in the same netcdf file.)
!!  In addition, a high-resolution netCDF topography file is required. This will define the grid that ICAR will run on.
!!  Finally an ICAR options file is used to specify various parameters for the model. A sample options file is provided
!!  in the run/ directory.
!!
!!  @section Developing
!!  This document provides the primary API and code structure documentation. The code is based on github.com/NCAR/icar
!!  Developers are encouraged to fork the main git repository and maintain their own git repository from which to
!!  issue pull requests.
!!
!!  @section Reference
!!  Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016),
!!  The Intermediate Complexity Atmospheric Research Model, J. Hydrometeor, doi:<a href="http://dx.doi.org/10.1175/JHM-D-15-0155.1">10.1175/JHM-D-15-0155.1</a>.
!!
!!------------------------------------------
