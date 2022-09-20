!> ----------------------------------------------------------------------------
!!  Main time stepping module.
!!  Calculates a stable time step (dt) and loops over physics calls
!!  Also updates boundaries every time step.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module time_step
    use data_structures             ! *_type  types and kCONSTANTS
    use icar_constants,             only : Rd, cp
    use microphysics,               only : mp
    use advection,                  only : advect
    use mod_atm_utilities,          only : exner_function
    use convection,                 only : convect
    use land_surface,               only : lsm
    use planetary_boundary_layer,   only : pbl
    use radiation,                  only : rad
    use wind,                       only : balance_uvw
    use domain_interface,           only : domain_t
    use boundary_interface,         only : boundary_t
    use options_interface,          only : options_t
    use debug_module,               only : domain_check
    use string,                     only : str
    use timer_interface,            only : timer_t

    implicit none
    private
    public :: step, update_dt

contains


    !!  Calculate the maximum stable time step given some CFL criteria
    !!
    !!  For each grid cell, find the mean of the wind speeds from each
    !!  direction * sqrt(3) for the 3D advection CFL limited time step
    !!  Also find the maximum wind speed anywhere in the domain to check
    !!  against a 1D advection limit.
    !!
    !! @param dx  [ scalar ]        horizontal grid cell width  [m]
    !! @param u   [nx+1 x nz x ny]  east west wind speeds       [m/s]
    !! @param v   [nx x nz x ny+1]  North South wind speed      [m/s]
    !! @param w   [nx x nz x ny]    vertical wind speed         [m/s]
    !! @param CFL [ scalar ]        CFL limit to use (e.g. 1.0)
    !! @return dt [ scalar ]        Maximum stable time step    [s]
    !!
    !!------------------------------------------------------------
    function compute_dt(dx, u, v, w, rho, dz, ims, ime, kms, kme, jms, jme, its, ite, jts, jte, CFL, cfl_strictness, use_density) result(dt)
        real,       intent(in)                   :: dx
        real,       intent(in), dimension(ims:ime+1,kms:kme,jms:jme) :: u 
        real,       intent(in), dimension(ims:ime,kms:kme,jms:jme+1) :: v
        real,       intent(in), dimension(ims:ime,kms:kme,jms:jme)   :: w, rho
        real,       intent(in), dimension(kms:kme)     :: dz
        integer,    intent(in)                   :: ims, ime, kms, kme, jms, jme, its, ite, jts, jte
        real,       intent(in)                   :: CFL
        integer,    intent(in)                   :: cfl_strictness
        logical,    intent(in)                   :: use_density
        
        ! output value
        real :: dt
        ! locals
        real :: three_d_cfl = 0.577350269 ! = sqrt(3)/3
        integer :: i, j, k, zoffset
        real :: maxwind3d, maxwind1d, current_wind, sqrt3

        sqrt3 = sqrt(3.0) * 1.001 ! with a safety factor

        maxwind1d = 0
        maxwind3d = 0

        if (cfl_strictness==1) then
            ! to ensure we are stable for 1D advection:
            if (use_density) then
                !maxwind1d = max( maxval(abs(u(2:,:,:) / (rho*dz*dx) )), maxval(abs(v(:,:,2:) / (rho*dz*dx))) )
                !maxwind1d = max( maxwind1d, maxval(abs(w/(rho*dz*dx))) )
            else
                maxwind1d = max( maxval(abs(u)), maxval(abs(v)))
                maxwind1d = max( maxwind1d, maxval(abs(w)))
            endif

            maxwind3d = maxwind1d * sqrt3
        else if (cfl_strictness==5) then

            if (use_density) then
                !maxwind1d = maxval(abs(u(2:,:,:) / (rho*dz*dx) )) &
                !          + maxval(abs(v(:,:,2:) / (rho*dz*dx) )) &
                !          + maxval(abs(w(:,:, :) / (rho*dz*dx) ))
            else
                maxwind3d = maxval(abs(u)) + maxval(abs(v)) + maxval(abs(w))
            endif

        else
            ! to ensure we are stable for 3D advection we'll use the average "max" wind speed
            ! but that average has to be divided by sqrt(3) for stability in 3 dimensional advection
            do j=jts,jte
                do k=kms,kme
                    if (k==kms) then
                        zoffset = 0
                    else
                        zoffset = -1
                    endif

                    do i=its,ite
                        ! just compute the sum of the wind speeds, but take the max over the two
                        ! faces of the grid cell (e.g. east and west sides)
                        ! this will be divided by 3 later by three_d_cfl
                        if (use_density) then
                            !current_wind = (max(abs(u(i,k,j)), abs(u(i+1,k,j))) &
                            !              + max(abs(v(i,k,j)), abs(v(i,k,j+1))) &
                            !              + max(abs(w(i,k,j)), abs(w(i,k+zoffset,j))) ) &
                            !              / (rho(i,k,j) * dz(i,k,j) * dx)
                        else
                            !current_wind = max(abs(u(i,k,j)), abs(u(i+1,k,j))) / dx &
                            !              +max(abs(v(i,k,j)), abs(v(i,k,j+1))) / dx &
                            !              +max(abs(w(i,k,j)), abs(w(i,k+zoffset,j))) / dz(k)
                                          
                            current_wind = max(( max( abs(u(i,k,j)), abs(u(i+1,k,j)) ) / dx), &
                                               ( max( abs(v(i,k,j)), abs(v(i,k,j+1)) ) / dx), &
                                               ( max( abs(w(i,k,j)), abs(w(i,k+zoffset,j)) ) / dz(k) ))
                        endif
                        maxwind3d = max(maxwind3d, current_wind)
                    ENDDO
                ENDDO
            ENDDO

            if (cfl_strictness==2) then
                ! effectively divides by 3 to take the mean and multiplies by the sqrt(3) for the 3D advection limit
                maxwind3d = maxwind3d * three_d_cfl

                ! to ensure we are stable for 1D advection:
                if (use_density) then
                    !maxwind1d = max( maxval(abs(u(2:,:,:) / (rho*dz*dx) )), maxval(abs(v(:,:,2:) / (rho*dz*dx))) )
                    !maxwind1d = max( maxwind1d, maxval(abs(w/(rho*dz*dx))) )
                else
                    maxwind1d = max( maxval(abs(u)), maxval(abs(v)))
                    maxwind1d = max( maxwind1d, maxval(abs(w)))
                endif
                ! also insure stability for 1D advection
                maxwind3d = max(maxwind1d,maxwind3d)

            ! else if (cfl_strictness==3) then
            !   leave maxwind3d as the sum of the max winds
            ! This should be the default, does it need to be multiplied by sqrt(3)?
            elseif (cfl_strictness==4) then
                maxwind3d = maxwind3d * sqrt3
            endif

        endif

        !TESTING: Do we need to multiply maxwind3d by sqrt3 as the comment above suggests?
        ! maxwind3d = maxwind3d * sqrt3

        dt = CFL / maxwind3d

        ! If we have too small a time step throw an error
        ! something is probably wrong in the physics or input data
        if (dt<1e-1) then
            write(*,*) "dt   = ", dt
            write(*,*) "Umax = ", maxval(abs(u))
            write(*,*) "Vmax = ", maxval(abs(v))
            write(*,*) "Wmax = ", maxval(abs(w))
            stop "ERROR time step too small"
        endif

    end function compute_dt


    !>------------------------------------------------------------
    !!  Prints progress to the terminal if requested
    !!
    !! @param current_time  the current state of the model time
    !! @param end_time      the end of the current full time step (when step will return)
    !! @param time_step     length of full time step to be integrated by step
    !! @param dt            numerical timestep to print for information
    !!
    !!------------------------------------------------------------
    subroutine print_progress(current_time, end_time, time_step, dt, last_time)
        implicit none
        type(Time_type),    intent(in)    :: current_time,    end_time
        type(time_delta_t), intent(in)    :: time_step,       dt
        real,               intent(inout) :: last_time

        type(time_delta_t) :: progress_dt
        real :: time_percent

        ! first compute the current time until reaching the end
        progress_dt  = (end_time - current_time)

        ! convert that to a percentage of the total time required
        time_percent = 100 - progress_dt%seconds() / time_step%seconds()  * 100

        ! finally if it has been at least 5% of the time since the last time we printed output, print output
        if (time_percent > (last_time + 5.0)) then
            last_time = last_time + 5.0
            ! this used to just use the nice $ (or advance="NO") trick, but at least with some mpi implementations, it buffers this output until it crashes
            write(*,"(A,f5.1,A,A)") char(13), max(0.0, time_percent)," %  dt=",trim(dt%as_string())
        endif

    end subroutine print_progress

    !>------------------------------------------------------------
    !! Update the numerical timestep to use
    !!
    !! @param dt            numerical timestep to use
    !! @param options       set options for how to update the time step
    !! @param domain        the full domain structure (need winds to compute dt)
    !! @param end_time      the end of the current full time step (when step will return)
    !!
    !!------------------------------------------------------------
    subroutine update_dt(dt, future_seconds, options, domain, update)
        implicit none
        type(time_delta_t), intent(inout) :: dt
        double precision,   intent(inout) :: future_seconds
        type(options_t),    intent(in)    :: options
        type(domain_t),     intent(in)    :: domain
        logical, optional,  intent(in)    :: update
        
        double precision                  :: seconds
        ! compute internal timestep dt to maintain stability
        ! courant condition for 3D advection. 
                
        !First, set dt to whatever the dt for the current time step was calculated to be
        call dt%set(seconds=min(future_seconds,120.0D0))

        if ((present(update)) .and. (update)) then
            seconds = compute_dt(domain%dx, domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, &
                            domain%w%meta_data%dqdt_3d, domain%density%data_3d, options%parameters%dz_levels, &
                            domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, &
                            domain%its, domain%ite, domain%jts, domain%jte, &
                            options%time_options%cfl_reduction_factor, &
                            cfl_strictness=options%time_options%cfl_strictness, use_density=.false.)
        else
            seconds = compute_dt(domain%dx, domain%u%data_3d, domain%v%data_3d, &
                            domain%w%data_3d, domain%density%data_3d, options%parameters%dz_levels, &
                            domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, &
                            domain%its, domain%ite, domain%jts, domain%jte, &
                            options%time_options%cfl_reduction_factor, &
                            cfl_strictness=options%time_options%cfl_strictness, use_density=.false.)
        endif
        ! perform a reduction across all images to find the minimum time step required
!#ifndef __INTEL_COMPILER
!        call CO_MIN(seconds)
!#endif
!#ifdef __INTEL_COMPILER
!         seconds = dx / 100
!#endif
        call CO_MIN(seconds)

        future_seconds = seconds
        ! set an upper bound on dt to keep microphysics and convection stable (?)
        ! store this back in the dt time_delta data structure
        
        ! If we are faster than the next time step, then adjust to the slower of the two
        if (dt%seconds() > future_seconds) call dt%set(seconds=min(future_seconds,120.0D0))
        if (this_image()==1) write(*,*) 'time_step: ',dt%seconds()

    end subroutine update_dt
    
!    subroutine update_dt(dt, options, domain, end_time)
!        implicit none
!        type(time_delta_t), intent(inout) :: dt
!        type(options_t),    intent(in)    :: options
!        type(domain_t),     intent(in)    :: domain
!        type(Time_type),    intent(in)    :: end_time

!        double precision :: seconds

        ! compute internal timestep dt to maintain stability
        ! courant condition for 3D advection. Note that w is normalized by dx/dz

        ! Note this needs to be performed when advect_density is enabled
        ! if (options%parameters%advect_density) then
            ! call dt%set(seconds=compute_dt(domain%dx, domain%ur, domain%vr, domain%wr, domain%rho, domain%dz_inter,&
            !                 options%time_options%cfl_reduction_factor, cfl_strictness=options%time_options%cfl_strictness,                   &
            !                 use_density=.True.))

        ! else
        ! compute the dt to meet the CFL criteria specified given dx, u, v, w, dz
!        associate(dx         => domain%dx,                              &
!                  u          => domain%u%data_3d(domain%its:domain%ite+1,:,domain%jts:domain%jte),                     &
!                  v          => domain%v%data_3d(domain%its:domain%ite,:,domain%jts:domain%jte+1),                     &
!                  w          => domain%w%data_3d(domain%its:domain%ite,:,domain%jts:domain%jte),                       &
!                  density    => domain%density%data_3d(domain%its:domain%ite,:,domain%jts:domain%jte),                 &
!                  dz         => options%parameters%dz_levels,           &
!                  reduction  => options%time_options%cfl_reduction_factor,&
!                  strictness => options%time_options%cfl_strictness       &
!            )

!            seconds = compute_dt(dx, u, v, w, density, dz, reduction, &
!                                 cfl_strictness=strictness, use_density=.false.)

!        end associate
        ! endif

!        seconds = compute_dt(dx, u, v, w, density, options%parameters%dz_levels, options%time_options%cfl_reduction_factor, &
!                                 cfl_strictness=options%time_options%cfl_strictness, use_density=.false.)
        ! perform a reduction across all images to find the minimum time step required
!#ifndef __INTEL_COMPILER
!        call CO_MIN(seconds)
!#endif
!#ifdef __INTEL_COMPILER
!         seconds = dx / 100
!#endif
!        call co_min(seconds)

        ! set an upper bound on dt to keep microphysics and convection stable (?)
        ! store this back in the dt time_delta data structure
!        call dt%set(seconds=min(seconds,120.0D0))

!    end subroutine update_dt


    !>------------------------------------------------------------
    !!  Step forward one IO time step.
    !!
    !!  Calculated the internal model time step to satisfy the CFL criteria,
    !!  then updates all forcing update increments for that dt and loops through
    !!  time calling physics modules.
    !!  Also checks to see if it is time to write a model output file.
    !!
    !! @param domain    domain data structure containing model state
    !! @param options   model options structure
    !! @param bc        model boundary conditions data structure
    !! @param next_output   Next time to write an output file (in "model_time")
    !!
    !!------------------------------------------------------------
    subroutine step(domain, forcing, end_time, dt_in, options, mp_timer, adv_timer, exch_timer)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(boundary_t),   intent(in)      :: forcing
        type(Time_type),    intent(in)      :: end_time
        type(time_delta_t), intent(in)      :: dt_in
        type(options_t),    intent(in)      :: options
        type(timer_t),      intent(inout)   :: mp_timer, adv_timer, exch_timer

        real :: last_print_time
        type(time_delta_t) :: time_step_size, dt

        last_print_time = 0.0
        time_step_size = end_time - domain%model_time

        ! now just loop over internal timesteps computing all physics in order (operator splitting...)
        dt = dt_in
        do while (domain%model_time < end_time)

            !call update_dt(dt, options, domain, end_time)
            ! Make sure we don't over step the forcing or output period
            if ((domain%model_time + dt) > end_time) then
                dt = end_time - domain%model_time
            endif

            ! ! apply/update boundary conditions including internal wind and pressure changes.
            call domain%apply_forcing(forcing,dt)

            ! if using advect_density winds need to be balanced at each update
            if (options%parameters%advect_density) call balance_uvw(domain,options)

            ! if an interactive run was requested than print status updates everytime at least 5% of the progress has been made
            if (options%parameters%interactive .and. (this_image()==1)) then
                call print_progress(domain%model_time, end_time, time_step_size, dt, last_print_time)
            endif
            ! this if is to avoid round off errors causing an additional physics call that won't really do anything

            if (real(dt%seconds()) > 1e-3) then

                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" init", fix=.True.)

                ! first process the halo section of the domain (currently hard coded at 1 should come from domain?)
                call rad(domain, options, real(dt%seconds()))
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" rad(domain", fix=.True.)

                call lsm(domain, options, real(dt%seconds()))!, halo=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" lsm")

                call pbl(domain, options, real(dt%seconds()))!, halo=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" pbl")

                call convect(domain, options, real(dt%seconds()))!, halo=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" convect")

                
                !Call, 'unify tendencies', where all physics tendencies are applied and forcing tendency is applied. Advection tendencies do not need to be calculated as long as
                ! the input arrays to advect are the same as the start of the time step.
                !After advection, apply tendencies to cells and call diagnostic update. Now MP can be called
                
                call exch_timer%start()
                call domain%halo_exchange_big()
                call exch_timer%stop()

                call adv_timer%start()

                call advect(domain, options, real(dt%seconds()))
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" advect(domain", fix=.True.)
                call adv_timer%stop()

                call domain%diagnostic_update(options)
                
                call mp_timer%start()

                call mp(domain, options, real(dt%seconds()))!, halo=domain%grid%halo_size)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" mp_halo", fix=.True.)
                call mp_timer%stop()
                !call domain%halo_send()
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%halo_send", fix=.True.)

                ! call rad(domain, options, real(dt%seconds()), subset=1)
                ! call lsm(domain, options, real(dt%seconds()))!, subset=1)
                ! call pbl(domain, options, real(dt%seconds()))!, subset=1)
                ! call convect(domain, options, real(dt%seconds()), subset=1)

                !call mp(domain, options, real(dt%seconds()), subset=domain%grid%halo_size)

                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" mp(domain", fix=.True.)
                !call domain%halo_retrieve()
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%halo_retrieve", fix=.True.)


                !If we are in the last ~10 updates of a time step and a variable drops below 0, we have probably over-shot a value of 0. Force back to 0
                if ((end_time%seconds() - domain%model_time%seconds()) < (dt%seconds()*10)) then

                    call domain%enforce_limits()
                endif


                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%apply_forcing", fix=.True.)

            endif
            ! step model_time forward
            domain%model_time = domain%model_time + dt
            
        enddo

    end subroutine step
end module time_step
