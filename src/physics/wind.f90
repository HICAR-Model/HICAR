!>------------------------------------------------------------
!! Module to manage the ICAR wind field, including calls to linear winds
!! importantly it also rotates the wind field into the ICAR grid and
!! balances the U, V, and W fields for "mass" conservation
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module wind
    use linear_theory_winds, only : linear_perturb
    !use mod_blocking,        only : update_froude_number, initialize_blocking
    use data_structures
    use exchangeable_interface,   only : exchangeable_t
    use domain_interface,  only : domain_t
    use options_interface, only : options_t
    use grid_interface,    only : grid_t
    use wind_surf, only         : apply_Sx
    use io_routines, only : io_read, io_write
    use mod_atm_utilities,   only : calc_froude, calc_Ri, calc_dry_stability
    
    implicit none
    private
    public::update_winds, init_winds
    real, parameter::deg2rad=0.017453293 !2*pi/360
contains

    !>------------------------------------------------------------
    !! Forces u,v, and w fields to balance
    !!       du/dx + dv/dy = dw/dz
    !!
    !! Starts by setting w out of the ground=0 then works through layers
    !!
    !!------------------------------------------------------------
    subroutine balance_uvw(u,v,w, jaco_u,jaco_v,jaco_w,dz,dx,jaco,rho,smooth_height, options,vert_weight)
        implicit none
        real,           intent(inout) :: u(:,:,:), v(:,:,:), w(:,:,:)
        real,           intent(in)    :: jaco_u(:,:,:), jaco_v(:,:,:), jaco_w(:,:,:), dz(:,:,:), jaco(:,:,:), rho(:,:,:)
        real,           intent(in)    :: dx, smooth_height
        type(options_t),intent(in)    :: options
        real, optional, intent(in)    :: vert_weight(:,:,:)

        real, allocatable, dimension(:,:) :: rhou, rhov, rhow
        real, allocatable, dimension(:,:,:) :: divergence
        real, allocatable, dimension(:,:,:) :: vert_div_weight
        integer :: ims, ime, jms, jme, kms, kme, k

        !if (this_image()==1) write(*,*) "vert_div_weight: ",vert_div_weight

        ! associate(u => domain%u%data_3d,  &
        !           v => domain%v%data_3d,  &
        !           w => domain%w%data_3d  )

        ims = lbound(w,1)
        ime = ubound(w,1)
        kms = lbound(w,2)
        kme = ubound(w,2)
        jms = lbound(w,3)
        jme = ubound(w,3)

        w = 0
        
        allocate(vert_div_weight(ims:ime,kms:kme,jms:jme))
        
        !Base assumption is to use same vertical splitting as Sherman 1978 (perfectly neuteral atmosphere)
        vert_div_weight = 1.0000000000
        if (present(vert_weight)) vert_div_weight=vert_weight

        !------------------------------------------------------------
        ! These could be module level to prevent lots of allocation/deallocation/reallocations
        ! but these are relatively small, and should be allocated efficiently on the stack
        !------------------------------------------------------------
        ! if (options%advect_density) then
        !     allocate(rhou(nx-1,ny-2))
        !     allocate(rhov(nx-2,ny-1))
        !     allocate(rhow(nx-2,ny-2))
        ! endif

        allocate(divergence(ims:ime,kms:kme,jms:jme))

        call calc_divergence(divergence,u,v,w,jaco_u,jaco_v,jaco_w,dz,dx,jaco,rho,smooth_height,options,horz_only=.True.)
        divergence = divergence * vert_div_weight

        ! If this becomes a bottle neck in the code it could be parallelized over y
        ! loop over domain levels
        do k = kms,kme
            !------------------------------------------------------------
            ! If we are incorporating density into the advection equation
            ! then it needs to be incorporated when balancing the wind field
            !
            ! Note that the "else" case below does the same thing without density
            ! and is much easier to understand
            !------------------------------------------------------------

            ! this is the else, advect density is not supported at the moment
            !------------------------------------------------------------
            ! If we are not incorporating density this is simpler
            !------------------------------------------------------------
            ! calculate horizontal divergence
            !   in the North-South direction
            !   in the East-West direction
            !   in net

            ! Then calculate w to balance
                ! if this is the first model level start from 0 at the ground
                ! note the are out for w is dx^2, but there is a dx in the divergence term that is dropped to balance
                if (options%parameters%advect_density) then
                    if (k==kms) then
                        w(:,k,:) = 0 - divergence(:,k,:) * dz(:,k,:) / (jaco_w(:,k,:) * (rho(:,k,:)+rho(:,k+1,:))/2 )
                    elseif (k==kme) then
                        w(:,k,:) = (w(:,k-1,:) * jaco_w(:,k-1,:) - divergence(:,k,:) * dz(:,k,:))/ (jaco_w(:,k,:) * rho(:,k,:))
                    else
                        w(:,k,:) = (w(:,k-1,:) * jaco_w(:,k-1,:) - divergence(:,k,:) * dz(:,k,:))/ (jaco_w(:,k,:) *  (rho(:,k,:)+rho(:,k+1,:))/2 )
                    endif
                else
                    if (k==kms) then
                        w(:,k,:) = 0 - divergence(:,k,:) * dz(:,k,:) / (jaco_w(:,k,:) )
                    else 
                        w(:,k,:) = (w(:,k-1,:) * jaco_w(:,k-1,:) - divergence(:,k,:) * dz(:,k,:))/ (jaco_w(:,k,:) )
                    end if
                    !u(:,k,:) = (u(:,k,:) * jaco_u(:,k,:) - divergence(:,k,:) * 0.25 * dx) / jaco_u(:,k,:)
                    !v(:,k,:) = (v(:,k,:) * jaco_v(:,k,:) - divergence(:,k,:) * 0.25 * dx) / jaco_v(:,k,:)


                
                end if
            end do
            if (present(vert_weight) .and. this_image()==1) then
                !call io_write("ideal_vert_div.nc","vert_div",divergence)
                !call io_write("ideal_w_grid.nc","w_grid",w)
            endif
            !------------------------------------------------------------
            ! Now do the same for the convective wind field if needed
            !------------------------------------------------------------
            ! if (options%physics%convection > 0) then
            !     ! calculate horizontal divergence
            !     dv = domain%v_cu(2:nx-1,i,3:ny) - domain%v_cu(2:nx-1,i,2:ny-1)
            !     du = domain%u_cu(3:nx,i,2:ny-1) - domain%u_cu(2:nx-1,i,2:ny-1)
            !     divergence = du + dv
            !     ! Then calculate w to balance
            !     if (i==1) then
            !         ! if this is the first model level start from 0 at the ground
            !         domain%w_cu(2:nx-1,i,2:ny-1) = 0 - divergence
            !     else
            !         ! else calculate w as a change from w at the level below
            !         domain%w_cu(2:nx-1,i,2:ny-1) = domain%w_cu(2:nx-1,i-1,2:ny-1) - divergence
            !     endif
            ! endif

        ! end associate

    end subroutine balance_uvw


    subroutine calc_divergence(div, u, v, w, jaco_u, jaco_v, jaco_w, dz, dx, jaco, rho, smooth_height,options,horz_only)
        implicit none
        real,           intent(inout) :: div(:,:,:)
        real,           intent(in)    :: u(:,:,:), v(:,:,:), w(:,:,:), dz(:,:,:), jaco_u(:,:,:), jaco_v(:,:,:), jaco_w(:,:,:), jaco(:,:,:), rho(:,:,:)
        real,           intent(in)    :: dx, smooth_height
        logical, optional, intent(in)  :: horz_only
        type(options_t),intent(in)    :: options

        real, allocatable, dimension(:,:,:) :: diff_U, diff_V, u_met, v_met, w_met
        integer :: ims, ime, jms, jme, kms, kme, k
        logical :: horz

        horz = .False.
        if (present(horz_only)) horz=horz_only

        ims = lbound(w,1)
        ime = ubound(w,1)
        kms = lbound(w,2)
        kme = ubound(w,2)
        jms = lbound(w,3)
        jme = ubound(w,3)

        allocate(diff_U(ims:ime,kms:kme,jms:jme))
        allocate(diff_V(ims:ime,kms:kme,jms:jme))
        allocate(u_met(ims:ime+1,kms:kme,jms:jme))
        allocate(v_met(ims:ime,kms:kme,jms:jme+1))
        allocate(w_met(ims:ime,kms:kme,jms:jme))

        !Multiplication of U/V by metric terms, converting jacobian to staggered-grid where possible, otherwise making assumption of 
        !Constant jacobian at edges
        
        if (options%parameters%advect_density) then
            u_met = u * jaco_u * (rho(ims:ime-1,:,jms:jme) + rho(ims+1:ime,:,jms:jme))/2
            v_met = v * jaco_v * (rho(ims:ime,:,jms:jme-1) + rho(ims:ime,:,jms+1:jme))/2
        else
            u_met = u * jaco_u
            v_met = v * jaco_v
        end if
        
        diff_U = u_met(ims+1:ime+1, :, jms:jme) - u_met(ims:ime, :, jms:jme)
        diff_V = v_met(ims:ime, :, jms+1:jme+1) - v_met(ims:ime, :, jms:jme)

        div(ims:ime,kms:kme,jms:jme) = (diff_U+diff_V)/(dx)

        if (.NOT.(horz)) then
            if (options%parameters%advect_density) then
                w_met(:,kme,:) = w(:,kme,:) * jaco_w(:,kme,:) * rho(:,kme,:)
                w_met(:,kms:kme-1,:) = w(:,kms:kme-1,:) * jaco_w(:,kms:kme-1,:) * (rho(:,kms+1:kme,:) + rho(:,kms:kme-1,:))/2 
            else
                w_met = w*jaco_w
            end if
            
            do k = kms,kme
                if (k == kms) then
                    div(ims:ime, k, jms:jme) = div(ims:ime, k, jms:jme) + w_met(ims:ime, k, jms:jme)/(dz(ims:ime, k, jms:jme))
                else
                    div(ims:ime, k, jms:jme) = div(ims:ime, k, jms:jme) + &
                                   (w_met(ims:ime,k,jms:jme)-w_met(ims:ime,k-1,jms:jme))/(dz(ims:ime,k,jms:jme))
                endif
            enddo
            !if (options%parameters%advect_density) then
            !    div = div/(jaco*rho)
            !else
            !    div = div/jaco
            !end if
        endif

    end subroutine calc_divergence
    


    !>------------------------------------------------------------
    !! Correct for a grid that is locally rotated with respect to EW,NS
    !!
    !! Assumes forcing winds are EW, NS relative, not grid relative.
    !!
    !!------------------------------------------------------------
    subroutine make_winds_grid_relative(u, v, w, grid, sintheta, costheta)
        real, intent(inout) :: u(:,:,:), v(:,:,:), w(:,:,:)
        type(grid_t),   intent(in)    :: grid
        double precision, intent(in)    :: sintheta(:,:), costheta(:,:)

        real, dimension(:,:), allocatable :: cos_shifted_u,sin_shifted_u,cos_shifted_v,sin_shifted_v,u_shifted_v,v_shifted_u
        real, dimension(:), allocatable :: u_local,v_local
        integer :: k, j, ims, ime, jms, jme, kms, kme, ids, ide, jds, jde

        ims = lbound(w,1)
        ime = ubound(w,1)
        kms = lbound(w,2)
        kme = ubound(w,2)
        jms = lbound(w,3)
        jme = ubound(w,3)

        ids = lbound(costheta,1)
        ide = ubound(costheta,1)
        jds = lbound(costheta,2)
        jde = ubound(costheta,2)

        allocate(u_local(ims:ime))
        allocate(v_local(ims:ime))

        !assumes u and v come in on a staggered Arakawa C-grid with one additional grid point in x/y for u/v respectively
        ! destagger to a centered grid (the mass grid)
        u(:ime,:,:) = (u(:ime,:,:) + u(ims+1:,:,:))/2
        v(:,:,:jme) = (v(:,:,:jme) + v(:,:,jms+1:))/2

        do j = jms, jme
            do k = kms, kme
                ! rotate wind field to the real grid
                u_local = u(ims:ime,k,j) * costheta(grid%its-1:grid%ite+1,j+grid%jts-2) &
                        + v(ims:ime,k,j) * sintheta(grid%its-1:grid%ite+1,j+grid%jts-2)
                v_local = v(ims:ime,k,j) * costheta(grid%its-1:grid%ite+1,j+grid%jts-2) &
                        + u(ims:ime,k,j) * sintheta(grid%its-1:grid%ite+1,j+grid%jts-2)
                u(:ime,k,j) = u_local
                v(:ime,k,j) = v_local
           enddo
        enddo
        deallocate(u_local,v_local)

        ! put the fields back onto a staggered grid, having effectively lost two grid cells in the staggered directions
        ! estimate the "lost" grid cells by extrapolating beyond the remaining
        u(ims+1:ime,:,:) =  (u(ims:ime-1,:,:) + u(ims+1:ime,:,:))/2
        u(ims,:,:)       = 2*u(ims,:,:)       - u(ims+1,:,:)
        u(ime+1,:,:)     = 2*u(ime,:,:)       - u(ime-1,:,:)

        v(:,:,jms+1:jme) =  (v(:,:,jms:jme-1) + v(:,:,jms+1:jme))/2
        v(:,:,jms)       = 2*v(:,:,jms)       - v(:,:,jms+1)
        v(:,:,jme+1)     = 2*v(:,:,jme)       - v(:,:,jme-1)

    end subroutine


    !>------------------------------------------------------------
    !! Apply wind field physics and adjustments
    !!
    !! This will call the linear wind module if necessary, otherwise it just updates for
    !! This should ONLY be called once for each forcing step, otherwise effects will be additive.
    !!
    !!------------------------------------------------------------
    subroutine update_winds(domain, options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options

        real, allocatable, dimension(:,:,:) :: temparray
        integer :: nx, ny, nz, i, j
        

        if (.not.allocated(domain%advection_dz)) then

            call init_winds(domain, options)
            !call initialize_blocking(domain, options)
            call update_stability(domain)

            ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
            call make_winds_grid_relative(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, domain%grid, domain%sintheta, domain%costheta)

            ! flow blocking parameterization
            ! if (options%block_options%block_flow) then
            !     call add_blocked_flow(domain, options)
            ! endif

            if (options%wind%Sx) call apply_Sx(domain)!%Sx,domain%TPI,domain%u%data_3d, domain%v%data_3d, domain%w%data_3d)

            ! linear winds
            if (options%physics%windtype==kWIND_LINEAR) then
                call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%parameters%advect_density)
            ! simple acceleration over topography
            elseif (options%physics%windtype==kCONSERVE_MASS) then
                if (options%parameters%use_terrain_difference) then  ! 
                !! use the ratio between hi-res and lo-res grid deformation (i.e. due to 'additional' terrain) for speedup
                    call mass_conservative_acceleration(domain%u%data_3d, domain%v%data_3d, domain%zfr_u, domain%zfr_v)
                else    
                    call mass_conservative_acceleration(domain%u%data_3d, domain%v%data_3d, domain%zr_u, domain%zr_v)
                endif    
            elseif (options%physics%windtype==kITERATIVE_WINDS) then
                call iterative_winds(domain, options)

            endif
            ! else assumes even flow over the mountains

            ! use horizontal divergence (convergence) to calculate vertical convergence (divergence)
            call balance_uvw(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, domain%jacobian_u, domain%jacobian_v, domain%jacobian_w, domain%advection_dz, domain%dx, domain%jacobian, domain%density%data_3d, domain%smooth_height, options)
            write(*,*) "Just in initilaization of winds"
        else

            call update_stability(domain)

            ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
            !call make_winds_grid_relative(domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, domain%w%meta_data%dqdt_3d, domain%grid, domain%sintheta, domain%costheta)
            
            if (options%wind%Sx) call apply_Sx(domain)!%Sx,domain%TPI,domain%u%meta_data%dqdt_3d,domain%v%meta_data%dqdt_3d, domain%w%meta_data%dqdt_3d)

            ! linear winds
            if (options%physics%windtype==kWIND_LINEAR) then
                call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%parameters%advect_density, update=.True.)
            ! simple acceleration over topography
            elseif (options%physics%windtype==kCONSERVE_MASS) then
                if (options%parameters%use_terrain_difference) then  ! 
                !! use the ratio between hi-res and lo-res grid deformation (i.e. due to 'addtional' terrain) for speedup
                    call mass_conservative_acceleration(domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, domain%zfr_u, domain%zfr_v)
                else    
                    call mass_conservative_acceleration(domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, domain%zr_u, domain%zr_v)
                endif   
            elseif (options%physics%windtype==kITERATIVE_WINDS) then
                call iterative_winds(domain, options, update_in=.True.)

            endif
            ! use horizontal divergence (convergence) to calculate vertical convergence (divergence)

            call balance_uvw(domain% u %meta_data%dqdt_3d,      &
                             domain% v %meta_data%dqdt_3d,      &
                             domain% w %meta_data%dqdt_3d,      &
                             domain%jacobian_u, domain%jacobian_v, domain%jacobian_w,         &
                             domain%advection_dz, domain%dx,    &
                             domain%jacobian, domain%density%data_3d, domain%smooth_height, options)
            write(*,*) "Just in Updaate of winds"

        endif


    end subroutine update_winds
    
    subroutine iterative_winds(domain, options, update_in)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        logical, optional, intent(in) :: update_in
        
        ! interal parameters
        real, allocatable, dimension(:,:,:) :: div, dial_weights, ADJ,ADJ_coef, U_cor, V_cor, current_u, current_v, current_w
        real, allocatable, dimension(:,:)   :: wind_layer
        real    :: corr_factor
        integer :: it, k, j, i, ims, ime, jms, jme, kms, kme, wind_k
        logical :: update
        
        update=.False.
        if (present(update_in)) update=update_in

        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)

        !If we are doing an update, we need to swap meta data into data_3d fields so it can be exchanged while balancing
        !First, we save a copy of the current data_3d so that we can substitute it back in later
        if (update) then
             current_u = domain%u%data_3d
             current_v = domain%v%data_3d
             current_w = domain%w%data_3d
             
             domain%u%data_3d = domain%u%meta_data%dqdt_3d
             domain%v%data_3d = domain%v%meta_data%dqdt_3d
             domain%w%data_3d = domain%w%meta_data%dqdt_3d
        endif
        
        !Do an initial exchange to make sure the U and V grids are similar for calculating w
        call domain%u%exchange_u()
        call domain%v%exchange_v()

        if (options%wind%Dial) then
            allocate(dial_weights(ims:ime,kms:kme,jms:jme))
            !do k = kms,kme
            !    if (k <= 10) then
            !        dial_weights(:,k,:) = 1.0
            !    else if (k < kme-10) then
            !        dial_weights(:,k,:) = (1 + (10/(kme-10)) - k/(kme-10))
            !    else 
            !        dial_weights(:,k,:) = (10/(kme-10))
            !    endif
            !enddo
            
            !Dial of 100% (all vertical) at Fr 1.25, 0% (all horizontal) at Fr <= 0.75
            dial_weights = max((domain%froude-0.75),0.0)/1.25
            
            !Sanity-bounding to min of 0, max of 1
            dial_weights = max(min(dial_weights,1.0),0.0)
                        
            if (this_image()==1) call io_write("dial_weights.nc","data",dial_weights)
            if (this_image()==1) call io_write("froude.nc","froude",domain%froude)
            if (this_image()==1) call io_write("Dick.nc","Ri",domain%Ri)
            
            !First call bal_uvw to generate an initial-guess for vertical winds
            call balance_uvw(domain%u%data_3d,domain%v%data_3d,domain%w%data_3d,domain%jacobian_u,domain%jacobian_v,domain%jacobian_w, & 
                             domain%advection_dz,domain%dx,domain%jacobian,domain%density%data_3d,domain%smooth_height,options,vert_weight=dial_weights)
        else
            !First call bal_uvw to generate an initial-guess for vertical winds
            call balance_uvw(domain%u%data_3d,domain%v%data_3d,domain%w%data_3d,domain%jacobian_u,domain%jacobian_v,domain%jacobian_w, & 
                        domain%advection_dz,domain%dx,domain%jacobian,domain%density%data_3d,domain%smooth_height,options)
        endif

        allocate(div(ims:ime,kms:kme,jms:jme))
        allocate(ADJ_coef(ims:ime,kms:kme,jms:jme))
        allocate(ADJ(ims:ime,kms:kme,jms:jme))
        allocate(U_cor(ims:ime,kms:kme,jms:jme))
        allocate(V_cor(ims:ime,kms:kme,jms:jme))
        allocate(wind_layer(ims:ime,jms:jme))
        ! Calculate and apply correction to w-winds
        wind_k = kme-2
        do k = kms,kme-2
            if (sum(domain%advection_dz(ims,1:k,jms)) > domain%smooth_height) then
                wind_k = k
                exit
            endif
        enddo

        !Store wind-layer to be subtracted, since we will be modifying vertical winds     
        wind_layer = domain%w%data_3d(:,wind_k,:)

        !Compute relative correction factors for U and V based on input speeds
        U_cor = ABS(domain%u%data_3d(ims:ime,:,jms:jme))/ &
                (ABS(domain%u%data_3d(ims:ime,:,jms:jme))+ABS(domain%v%data_3d(ims:ime,:,jms:jme)))
                
        do k = kms,kme
            !corr_factor = ((sum(domain%advection_dz(ims,1:k,jms)))/domain%smooth_height)
            corr_factor = (k*1.0)/wind_k
            corr_factor = min(corr_factor,1.0)
            do i = ims,ime
                do j = jms,jme
                    !domain%w%data_3d(i,k,j) = domain%w%data_3d(i,k,j) - corr_factor * wind_layer(i,j)
                    
                    !if ( (domain%u%data_3d(i,k,j)+domain%v%data_3d(i,k,j)) == 0) U_cor(i,k,j) = 0.5
                enddo
            enddo
            ! Compute this now, since it wont change in the loop
            ADJ_coef(:,k,:) = -(4*domain%jacobian(:,k,:))/domain%dx 
        enddo
        

        !domain%w%data_3d(:,kms,:) = 0
        if (this_image()==4) call io_write("ideal_w_grid_before.nc","w_grid",domain%w%data_3d)
        !do k = kms,kme
        !    if (kms==1) then
        !        domain%w%data_3d(:,k,:) = domain%w%data_3d(:,k,:) !( -2*(U_cor(:,k,:) + V_cor(:,k,:))/domain%jacobian(:,k,:) ) - 0
        !    else
        !        domain%w%data_3d(:,k,:) = ( -2*(U_cor(:,k,:) + V_cor(:,k,:))/domain%jacobian(:,k,:) ) - domain%w%data_3d(:,k-1,:)
        !    endif
        !enddo
        
        U_cor    = (domain%u%data_3d(ims+1:ime+1, :, :) * domain%dzdx(ims+1:ime+1,:,jms:jme) + &
                    domain%u%data_3d(ims:ime, :, :) * domain%dzdx(ims:ime,:,jms:jme))/2
                ! compute the V * dz/dy component of vertical motion
        V_cor    = (domain%v%data_3d(:, :, jms+1:jme+1) * domain%dzdy(ims:ime,:,jms+1:jme+1) + &
                    domain%v%data_3d(:, :, jms:jme) * domain%dzdy(ims:ime,:,jms:jme))/2
        
        !domain%w%data_3d = -( U_cor + V_cor)/domain%jacobian
        
        if (this_image()==2) call io_write("ideal_w_grid_after.nc","w_grid",domain%w%data_3d)
        
        !U_cor = ABS(sin(domain%dzdx(ims:ime,kms:kme,jms:jme))) / &
        !        ( ABS(sin(domain%dzdx(ims:ime,kms:kme,jms:jme))) + ABS(sin(domain%dzdy(ims:ime,kms:kme,jms:jme))) )
        
        !V_cor = ABS(domain%dzdx)!0.5
        U_cor = 0.5 !+ ( ABS(sin(domain%dzdx(ims:ime,kms:kme,jms:jme))) - ABS(sin(domain%dzdy(ims:ime,kms:kme,jms:jme))) )/2
        V_cor = 1 - U_cor 
        
        if (this_image()==2) call io_write("U_cor.nc","U_cor",U_cor)
        if (this_image()==2) call io_write("V_cor.nc","V_cor",V_cor)
        
        ! Now, fixing w-winds, iterate over U/V to reduce divergence with new w-winds
        do it = 0,options%parameters%wind_iterations
        
            !U_cor    = (domain%u%data_3d(ims+1:ime+1, :, :) * domain%dzdx(ims+1:ime+1,:,jms:jme) + &
            !        domain%u%data_3d(ims:ime, :, :) * domain%dzdx(ims:ime,:,jms:jme))/2
                ! compute the V * dz/dy component of vertical motion
            !V_cor    = (domain%v%data_3d(:, :, jms+1:jme+1) * domain%dzdy(ims:ime,:,jms+1:jme+1) + &
            !        domain%v%data_3d(:, :, jms:jme) * domain%dzdy(ims:ime,:,jms:jme))/2
        
            !domain%w%data_3d = -( U_cor + V_cor)/domain%jacobian
        
        
            !Compute divergence in new wind field
            call calc_divergence(div,domain%u%data_3d,domain%v%data_3d,domain%w%data_3d,domain%jacobian_u,domain%jacobian_v,domain%jacobian_w, &
            domain%advection_dz,domain%dx,domain%jacobian,domain%density%data_3d,domain%smooth_height,options)

            !Compute adjustment based on divergence
            ADJ = div/ADJ_coef
            
            !Distribute divergence among the U and V fields                                                                    
            domain%u%data_3d(ims+2:ime,:,jms+1:jme-1) = domain%u%data_3d(ims+2:ime,:,jms+1:jme-1) + &
                                                        (ADJ(ims+1:ime-1,:,jms+1:jme-1) * U_cor(ims+2:ime,:,jms+1:jme-1))
                                                        
            domain%u%data_3d(ims+2:ime,:,jms+1:jme-1) = domain%u%data_3d(ims+2:ime,:,jms+1:jme-1) - &
                                                        (ADJ(ims+2:ime,:,jms+1:jme-1) * U_cor(ims+2:ime,:,jms+1:jme-1))
                                                        
            domain%v%data_3d(ims+1:ime-1,:,jms+2:jme) = domain%v%data_3d(ims+1:ime-1,:,jms+2:jme) + &
                                                        (ADJ(ims+1:ime-1,:,jms+1:jme-1) * V_cor(ims+1:ime-1,:,jms+2:jme))
                                                        
            domain%v%data_3d(ims+1:ime-1,:,jms+2:jme) = domain%v%data_3d(ims+1:ime-1,:,jms+2:jme) - &
                                                        (ADJ(ims+1:ime-1,:,jms+2:jme) * V_cor(ims+1:ime-1,:,jms+2:jme))
            
            !U_cor = ABS(domain%u%data_3d(ims:ime,:,jms:jme))/ &
            !    (ABS(domain%u%data_3d(ims:ime,:,jms:jme))+ABS(domain%v%data_3d(ims:ime,:,jms:jme)))
            !U_cor = min(max(U_cor,0.0),1.0)
            !V_cor = 1.0 - U_cor 
        
            call domain%u%exchange_u()
            call domain%v%exchange_v()
            
        enddo
        
        

        !If an update loop, swap meta_data and data_3d fields back
        if (update) then
            domain%u%meta_data%dqdt_3d = domain%u%data_3d
            domain%v%meta_data%dqdt_3d = domain%v%data_3d
            domain%w%meta_data%dqdt_3d = domain%w%data_3d
            
            domain%u%data_3d = current_u
            domain%v%data_3d = current_v
            domain%w%data_3d = current_w
        endif
        
    end subroutine iterative_winds

    subroutine mass_conservative_acceleration(u, v, u_accel, v_accel)
        implicit none
        real, intent(inout) :: u(:,:,:)
        real, intent(inout) :: v(:,:,:)
        real, intent(in)    :: u_accel(:,:,:)
        real, intent(in)    :: v_accel(:,:,:)

        u = u / u_accel
        v = v / v_accel

    end subroutine mass_conservative_acceleration

    !>------------------------------------------------------------
    !! Setup initial fields (i.e. grid relative rotation fields)
    !!
    !!------------------------------------------------------------
    subroutine init_winds(domain,options)
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options

        integer :: i, j, ims, ime, jms, jme, kms, kme
        integer :: starti, endi
        double precision :: dist, dlat, dlon

        real, allocatable :: temporary_2d(:,:)

        call allocate_winds(domain)

        if (options%parameters%fixed_dz_advection) then
            do i=domain%grid%kms, domain%grid%kme
                domain%advection_dz(:,i,:) = options%parameters%dz_levels(i)
            enddo
        else
            domain%advection_dz = domain%dz_interface%data_3d
        endif


    end subroutine init_winds

    !>------------------------------------------------------------
    !! Allocate memory used in various wind related routines
    !!
    !!------------------------------------------------------------
    subroutine allocate_winds(domain)
        type(domain_t), intent(inout) :: domain
        integer :: ims, ime, jms, jme, kms, kme

        ims = lbound(domain%latitude%data_2d, 1)
        ime = ubound(domain%latitude%data_2d, 1)
        jms = lbound(domain%latitude%data_2d, 2)
        jme = ubound(domain%latitude%data_2d, 2)
        kms = lbound(domain%w%data_3d, 2)
        kme = ubound(domain%w%data_3d, 2)

        if (.not.allocated(domain%advection_dz)) then
            allocate(domain%advection_dz(ims:ime,kms:kme,jms:jme))
        endif

        ! note w is special cased because it does not have a forcing variable, so it is not necessarily allocated automatically
        if (.not.associated(domain%w%meta_data%dqdt_3d)) then
            allocate(domain%w%meta_data%dqdt_3d(ims:ime,kms:kme,jms:jme))
            domain%w%meta_data%dqdt_3d = 0
        endif

        ! if (.not.allocated(domain%dzdx)) then
        !     allocate(domain%dzdx(nx-1,ny))
        ! endif
        ! if (.not.allocated(domain%dzdy)) then
        !     allocate(domain%dzdy(nx,ny-1))
        ! endif

    end subroutine allocate_winds
    
    subroutine update_stability(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        real, allocatable, dimension(:,:,:) :: wind_speed, temp_froude
        
        integer :: k, j, i, n, ims, ime, jms, jme, kms, kme
        real :: z_top, z_bot, th_top, th_bot, stability
        integer :: ymin, ymax, xmin, xmax, n_smoothing_passes, nsmooth_gridcells
        
        n_smoothing_passes = 4
        nsmooth_gridcells = int(1000 / domain%dx)
        
        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)

        !If it is our first time calculating Fr, allocate and populate froude_terrain array
        if (.not.allocated(domain%froude)) then
            allocate(domain%froude(ims:ime,kms:kme,jms:jme))            
            allocate(domain%froude_terrain(ims:ime,jms:jme))

            call compute_terrain_blocking_heights(domain)
        endif
       
        if (.not.allocated(domain%Ri)) allocate(domain%Ri(ims:ime,kms:kme,jms:jme))
       
        allocate(wind_speed(ims:ime,kms:kme,jms:jme))
        !allocate(temp_froude(ims:ime,kms:kme,jms:jme))       
        
        wind_speed = sqrt( ((domain%u%data_3d(ims+1:ime+1,:,:) + domain%u%data_3d(ims:ime,:,:))/2)**2 + &
                           ((domain%v%data_3d(:,:,jms+1:jme+1) + domain%v%data_3d(:,:,jms:jme))/2)**2 )
        
        !Since we will loop up to nz-1, we set all Fr to 0.1, which will leave the upper layer as very stable
        domain%froude = 0.1
        
        !Since we will loop up to nz-1, we set all Ri here to 10
        domain%Ri = 10.0
        
        do i = ims,ime
            do j = jms,jme
                do k = kms,kme-1
                    th_bot = domain%potential_temperature%data_3d(i,k,j)
                    th_top = domain%potential_temperature%data_3d(i,k+1,j)
                    z_bot  = domain%z%data_3d(i,k,j)
                    z_top  = domain%z%data_3d(i,k+1,j)
                    stability = calc_dry_stability(th_top, th_bot, z_top, z_bot) 
                    
                    domain%Ri(i,k,j) =  calc_Ri(stability, wind_speed(i,k,j), domain%advection_dz(i,k,j) )

                    stability = sqrt(max(stability, 0.))
                    domain%froude(i,k,j) = calc_froude(stability, domain%froude_terrain(i,j), wind_speed(i,k,j))
                enddo
            enddo
        enddo

        !do n = 1,n_smoothing_passes
        !    do j=jms,jme
        !        ymin = max(j-nsmooth_gridcells, jms)
        !        ymax = min(j+nsmooth_gridcells, jme)
        !        do i=ims,ime
        !            xmin = max(i-nsmooth_gridcells, ims)
        !            xmax = min(i+nsmooth_gridcells, ime)
        !            do k=kms,kme
        !                !write(*,*) "temp_f:  ", sum(temp_froude(xmin:xmax,k,ymin:ymax))
        !                !write(*,*) "num_sum:  ", ((xmax-xmin+1) * (ymax-ymin+1))
        !                domain%froude(i,k,j) = sum(temp_froude(xmin:xmax,k,ymin:ymax)) / ((xmax-xmin+1) * (ymax-ymin+1))
        !            enddo
        !        enddo
        !    enddo

        !    if (n/=n_smoothing_passes) then
        !        temp_froude = domain%froude
        !    endif
        !enddo

    end subroutine update_stability

    !>-----------------------------------------
    !> Compute a smoothed terrain varience field for use in Froude number calculation
    !>
    !------------------------------------------
    subroutine compute_terrain_blocking_heights(domain) !froude_terrain, terrain)
        implicit none
        type(domain_t), intent(inout) :: domain

        integer :: nx, ny, x, y
        integer :: xs,xe, ys,ye, n, np
        integer :: window_size, smooth_window, n_smoothing_passes
        real, allocatable :: temp_terrain(:,:)

        n_smoothing_passes = 4
        window_size   = int(max(5000.0/domain%dx,1.0)) !Compute Froude-terrain as the diff in min and max over a 5000m search window
        smooth_window = 2

        nx = size(domain%global_terrain,1)
        ny = size(domain%global_terrain,2)
        
        allocate(temp_terrain(nx,ny))

        ! first compute lightly smoothed terrain
        !do y=1,ny
        !    ys = max( y - smooth_window, 1)
        !    ye = min( y + smooth_window, ny)
        !    do x=1,nx
        !        xs = max( x - smooth_window, 1)
        !        xe = min( x + smooth_window, nx)
        !        n = (xe-xs+1) * (ye-ys+1)
        !        temp_terrain(x,y) = sum(domain%global_terrain(xs:xe,ys:ye)) / n
        !    enddo
        !enddo

        ! then compute the range of terrain (max-min) in a given window
        do y=domain%jts-1,domain%jte+1
            ys = max( y - window_size, 1)
            ye = min( y + window_size, ny)
            do x=domain%its-1,domain%ite+1
                xs = max( x - window_size, 1)
                xe = min( x + window_size, nx)
                domain%froude_terrain(x,y) = maxval(domain%global_terrain(xs:xe,ys:ye)) - &
                                             minval(domain%global_terrain(xs:xe,ys:ye))
            enddo
        enddo
        ! call io_write("initial_terrain_delta.nc","data",temp_terrain)

        ! finally smooth that terrain delta field a few times as well

        !do np=1,n_smoothing_passes
        !    do y=1,ny
        !        ys = max( y - smooth_window, 1)
        !        ye = min( y + smooth_window, ny)
        !        do x=1,nx
        !            xs = max( x - smooth_window, 1)
        !            xe = min( x + smooth_window, nx)
        !            n = (xe-xs+1) * (ye-ys+1)
        !            froude_terrain(x,y) = sum(temp_terrain(xs:xe,ys:ye)) / n
        !        enddo
        !    enddo
        !    if (np /= n_smoothing_passes) then
        !        temp_terrain = froude_terrain
        !    endif
        !enddo
        call io_write("terrain_blocking.nc","data",domain%froude_terrain)

    end subroutine compute_terrain_blocking_heights


    !>------------------------------------------------------------
    !! Provides a routine to deallocate memory allocated in allocate_winds
    !!
    !!------------------------------------------------------------
    ! subroutine finalize_winds(domain)
    !     type(domain_t), intent(inout) :: domain
    !
    !     if (allocated(domain%sintheta)) then
    !         deallocate(domain%sintheta)
    !     endif
    !     if (allocated(domain%costheta)) then
    !         deallocate(domain%costheta)
    !     endif
    !     if (allocated(domain%dzdx)) then
    !         deallocate(domain%dzdx)
    !     endif
    !     if (allocated(domain%dzdy)) then
    !         deallocate(domain%dzdy)
    !     endif
    !
    ! end subroutine finalize_winds
end module wind
