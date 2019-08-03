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
    use icar_constants,             only : Rd
    use microphysics,               only : mp
    use wind,                       only : balance_uvw, update_winds
    use advection,                  only : advect
    use mod_atm_utilities,          only : exner_function
    ! use convection,                 only : convect
    use land_surface,               only : lsm
    ! use planetary_boundary_layer,   only : pbl
    use radiation,                  only : rad

    use domain_interface,           only : domain_t
    use options_interface,          only : options_t
    use debug_module,               only : domain_check
    use string,                     only : str

    implicit none
    private
    public :: step

contains


    !>------------------------------------------------------------
    !! Update model diagnostic fields
    !!
    !! Calculates most model diagnostic fields such as Psfc, 10m height winds and ustar
    !!
    !! @param domain    Model domain data structure to be updated
    !! @param options   Model options (not used at present)
    !!
    !!------------------------------------------------------------
    subroutine diagnostic_update(domain, options)
        implicit none
        type(domain_t),  intent(inout)   :: domain
        type(options_t), intent(in)      :: options

        associate(exner                 => domain%exner%data_3d,    &
                  pressure              => domain%pressure%data_3d, &
                  density               => domain%density%data_3d,  &
                  potential_temperature => domain%potential_temperature%data_3d)

        exner = exner_function(pressure)

        density =  pressure / &
                    (Rd * potential_temperature * exner) ! kg/m^3

        end associate

    ! NOTE: all code below is not implemented in ICAR 2.0 yet
    ! it is left as a reminder of what needs to be done, and example when the time comes
    !
    !     integer :: nx,ny,nz, y, z
    !
    !     nx=size(domain%p,1)
    !     nz=size(domain%p,2)
    !     ny=size(domain%p,3)
    !
    !     ! update p_inter, psfc, ptop, Um, Vm, mut
    !     domain%Um = 0.5*(domain%u(1:nx-1,:,:)+domain%u(2:nx,:,:))
    !     domain%Vm = 0.5*(domain%v(:,:,1:ny-1)+domain%v(:,:,2:ny))
    !     if (options%physics%convection>0) then
    !         domain%Um = domain%Um + 0.5*(domain%u_cu(1:nx-1,:,:)+domain%u_cu(2:nx,:,:))
    !         domain%Vm = domain%Vm + 0.5*(domain%v_cu(:,:,1:ny-1)+domain%v_cu(:,:,2:ny))
    !     endif
    !     domain%t  = domain%th * domain%pii
    !
    !     domain%p_inter=domain%p
    !     call update_pressure(domain%p_inter, domain%z, domain%z_inter, domain%t)
    !     domain%psfc = domain%p_inter(:,1,:)
    !     ! technically this isn't correct, we should be using update_pressure or similar to solve this
    !     domain%ptop = 2*domain%p(:,nz,:) - domain%p(:,nz-1,:)
    !
    !     ! dry mass in the gridcell is equivalent to the difference in pressure from top to bottom
    !     domain%mut(:,1:nz-1,:) = domain%p_inter(:,1:nz-1,:) - domain%p_inter(:,2:nz,:)
    !     domain%mut(:,nz,:) = domain%p_inter(:,nz,:) - domain%ptop
    !
    !     if (.not.allocated(lastw)) then
    !         allocate(lastw(nx-2,ny-2))
    !         allocate(currw(nx-2,ny-2))
    !         allocate(uw(nx-1,ny-2))
    !         allocate(vw(nx-2,ny-1))
    !     endif
    !
    !     ! temporary constant
    !     ! use log-law of the wall to convert from first model level to surface
    !     currw = karman / log((domain%z(2:nx-1,1,2:ny-1)-domain%terrain(2:nx-1,2:ny-1)) / domain%znt(2:nx-1,2:ny-1))
    !     ! use log-law of the wall to convert from surface to 10m height
    !     lastw = log(10.0 / domain%znt(2:nx-1,2:ny-1)) / karman
    !     domain%ustar(2:nx-1,2:ny-1) = domain%Um   (2:nx-1,1,2:ny-1) * currw
    !     domain%u10  (2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1)   * lastw
    !     domain%ustar(2:nx-1,2:ny-1) = domain%Vm   (2:nx-1,1,2:ny-1) * currw
    !     domain%v10  (2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1)   * lastw
    !
    !     ! now calculate master ustar based on U and V combined in quadrature
    !     domain%ustar(2:nx-1,2:ny-1) = sqrt(domain%Um(2:nx-1,1,2:ny-1)**2 + domain%Vm(2:nx-1,1,2:ny-1)**2) * currw
    !
    !     ! finally, calculate the real vertical motions (including U*dzdx + V*dzdy)
    !     lastw=0
    !     do z=1,nz
    !         ! compute the U * dz/dx component of vertical motion
    !         uw    = domain%u(2:nx,  z,2:ny-1) * domain%dzdx(:,2:ny-1)
    !         ! compute the V * dz/dy component of vertical motion
    !         vw    = domain%v(2:nx-1,z,2:ny  ) * domain%dzdy(2:nx-1,:)
    !         ! convert the W grid relative motion to m/s
    !         currw = domain%w(2:nx-1,z,2:ny-1) * domain%dz_inter(2:nx-1,z,2:ny-1) / domain%dx
    !
    !         if (options%physics%convection>0) then
    !             currw = currw + domain%w_cu(2:nx-1,z,2:ny-1) * domain%dz_inter(2:nx-1,z,2:ny-1) / domain%dx
    !         endif
    !
    !         ! compute the real vertical velocity of air by combining the different components onto the mass grid
    !         ! includes vertical interpolation between w_z-1/2 and w_z+1/2
    !         domain%w_real(2:nx-1,z,2:ny-1) = (uw(1:nx-2,:) + uw(2:nx-1,:))*0.5 &
    !                                         +(vw(:,1:ny-2) + vw(:,2:ny-1))*0.5 &
    !                                         +(lastw + currw) * 0.5
    !
    !         lastw=currw ! could avoid this memcopy cost using pointers or a single manual loop unroll
    !     end do
    !
    end subroutine diagnostic_update


    !>------------------------------------------------------------
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
    function compute_dt(dx, u, v, w, rho, dz, CFL, cfl_strictness, use_density) result(dt)
        real,       intent(in)                   :: dx
        real,       intent(in), dimension(:,:,:) :: u, v, w, rho, dz
        real,       intent(in)                   :: CFL
        integer,    intent(in)                   :: cfl_strictness
        logical,    intent(in)                   :: use_density
        ! output value
        real :: dt
        ! locals
        real :: three_d_cfl = 0.577350269 ! = sqrt(3)/3
        integer :: i, j, k, nx, nz, ny, zoffset
        real :: maxwind3d, maxwind1d, current_wind, sqrt3

        sqrt3 = sqrt(3.0) * 1.001 ! with a safety factor

        maxwind1d = 0
        maxwind3d = 0

        nx = size(rho,1)
        nz = size(rho,2)
        ny = size(rho,3)

        if (cfl_strictness==1) then
            ! to ensure we are stable for 1D advection:
            if (use_density) then
                maxwind1d = max( maxval(abs(u(2:,:,:) / (rho*dz*dx) )), maxval(abs(v(:,:,2:) / (rho*dz*dx))) )
                maxwind1d = max( maxwind1d, maxval(abs(w/(rho*dz*dx))) )
            else
                maxwind1d = max( maxval(abs(u)), maxval(abs(v)))
                maxwind1d = max( maxwind1d, maxval(abs(w)))
            endif

            maxwind3d = maxwind1d * sqrt3
        else if (cfl_strictness==5) then

            if (use_density) then
                maxwind1d = maxval(abs(u(2:,:,:) / (rho*dz*dx) )) &
                          + maxval(abs(v(:,:,2:) / (rho*dz*dx) )) &
                          + maxval(abs(w(:,:, :) / (rho*dz*dx) ))
            else
                maxwind3d = maxval(abs(u)) + maxval(abs(v)) + maxval(abs(w))
            endif

        else
            ! to ensure we are stable for 3D advection we'll use the average "max" wind speed
            ! but that average has to be divided by sqrt(3) for stability in 3 dimensional advection
            do j=1,ny
                do k=1,nz
                    if (k==1) then
                        zoffset = 0
                    else
                        zoffset = -1
                    endif

                    do i=1,nx
                        ! just compute the sum of the wind speeds, but take the max over the two
                        ! faces of the grid cell (e.g. east and west sides)
                        ! this will be divided by 3 later by three_d_cfl
                        if (use_density) then
                            current_wind = (max(abs(u(i,k,j)), abs(u(i+1,k,j))) &
                                          + max(abs(v(i,k,j)), abs(v(i,k,j+1))) &
                                          + max(abs(w(i,k,j)), abs(w(i,k+zoffset,j))) ) &
                                          / (rho(i,k,j) * dz(i,k,j) * dx)
                        else
                            current_wind = max(abs(u(i,k,j)), abs(u(i+1,k,j))) &
                                          +max(abs(v(i,k,j)), abs(v(i,k,j+1))) &
                                          +max(abs(w(i,k,j)), abs(w(i,k+zoffset,j)))
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
                    maxwind1d = max( maxval(abs(u(2:,:,:) / (rho*dz*dx) )), maxval(abs(v(:,:,2:) / (rho*dz*dx))) )
                    maxwind1d = max( maxwind1d, maxval(abs(w/(rho*dz*dx))) )
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

        dt = CFL * dx / maxwind3d

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
    subroutine update_dt(dt, options, domain, end_time)
        implicit none
        type(time_delta_t), intent(inout) :: dt
        type(options_t),    intent(in)    :: options
        type(domain_t),     intent(in)    :: domain
        type(Time_type),    intent(in)    :: end_time

        double precision :: seconds

        ! compute internal timestep dt to maintain stability
        ! courant condition for 3D advection. Note that w is normalized by dx/dz

        ! Note this needs to be performed when advect_density is enabled
        ! if (options%parameters%advect_density) then
            ! call dt%set(seconds=compute_dt(domain%dx, domain%ur, domain%vr, domain%wr, domain%rho, domain%dz_inter,&
            !                 options%parameters%cfl_reduction_factor, cfl_strictness=options%parameters%cfl_strictness,                   &
            !                 use_density=.True.))

        ! else
        ! compute the dt to meet the CFL criteria specified given dx, u, v, w, dz
        associate(dx         => domain%dx,                              &
                  u          => domain%u%data_3d,                       &
                  v          => domain%v%data_3d,                       &
                  w          => domain%w%data_3d,                       &
                  density    => domain%density%data_3d,                 &
                  dz         => domain%dz_interface%data_3d,            &
                  reduction  => options%parameters%cfl_reduction_factor,&
                  strictness => options%parameters%cfl_strictness       &
            )

            seconds = compute_dt(dx, u, v, w, density, dz, reduction, &
                                 cfl_strictness=strictness, use_density=.false.)

        end associate
        ! endif

        ! perform a reduction across all images to find the minimum time step required
#ifndef __INTEL_COMPILER
        call co_min(seconds)
#endif
        seconds = seconds/1.25
        ! seconds=60
        ! print*, seconds
        ! set an upper bound on dt to keep microphysics and convection stable (?)
        ! store this back in the dt time_delta data structure
        call dt%set(seconds=min(seconds,120.0D0))

    end subroutine update_dt


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
    subroutine step(domain, end_time, options)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(Time_type),    intent(in)      :: end_time
        type(options_t),    intent(in)      :: options

        real :: last_print_time
        type(time_delta_t) :: dt, time_step_size
        integer :: ims, jms, kms
        integer :: test_x, test_y, test_z, test_image

        test_x = 87
        test_z = 3
        test_y = 198
        test_image = 33

        ims = domain%grid%ims
        jms = domain%grid%jms
        kms = domain%grid%kms

        last_print_time = 0.0
        time_step_size = end_time - domain%model_time


        call update_dt(dt, options, domain, end_time)

        ! now just loop over internal timesteps computing all physics in order (operator splitting...)
        do while (domain%model_time < end_time)
            ! Make sure we don't over step the forcing or output period
            if ((domain%model_time + dt) > end_time) then
                dt = end_time - domain%model_time
            endif

            ! ensure internal model consistency
            call diagnostic_update(domain, options)


            ! if an interactive run was requested than print status updates everytime at least 5% of the progress has been made
            if (options%parameters%interactive .and. (this_image()==1)) then
                call print_progress(domain%model_time, end_time, time_step_size, dt, last_print_time)
            endif

            ! this if is to avoid round off errors causing an additional physics call that won't really do anything
            if (dt%seconds() > 1e-3) then

                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" init")

                ! first process the halo section of the domain (currently hard coded at 1 should come from domain?)
                call rad(domain, options, real(dt%seconds()))
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" rad(domain")

                ! call lsm(domain, options, real(dt%seconds()), halo=1)
                ! call pbl(domain, options, real(dt%seconds()), halo=1)
                ! call convect(domain, options, real(dt%seconds()), halo=1)

                call mp(domain, options, real(dt%seconds()), halo=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" mp_halo")

                call domain%halo_send()
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%halo_send")

                ! call rad(domain, options, real(dt%seconds()), subset=1)
                ! call lsm(domain, options, real(dt%seconds()))!, subset=1)
                ! call pbl(domain, options, real(dt%seconds()), subset=1)
                ! call convect(domain, options, real(dt%seconds()), subset=1)

                call mp(domain, options, real(dt%seconds()), subset=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" mp(domain")

                call domain%halo_retrieve()
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%halo_retrieve")

                call advect(domain, options, real(dt%seconds()))
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" advect(domain")


                ! ! apply/update boundary conditions including internal wind and pressure changes.
                call domain%apply_forcing(dt)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%apply_forcing")

            endif

            ! step model_time forward
            domain%model_time = domain%model_time + dt

        enddo

    end subroutine step
end module time_step
