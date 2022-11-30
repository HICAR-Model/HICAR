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
    use wind_iterative,      only : calc_iter_winds
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
    public:: balance_uvw, update_winds, init_winds, calc_w_real, wind_var_request

    real, parameter::deg2rad=0.017453293 !2*pi/360
    real, parameter :: rad2deg=57.2957779371
    real, parameter :: DEFAULT_FR_L = 1000.0
contains



        subroutine wind_linear_var_request(options)
            implicit none
            type(options_t), intent(inout) :: options

            ! List the variables that are required to be allocated for the linear wind solution
            call options%alloc_vars( &
                            [kVARS%nsquared,    kVARS%potential_temperature,   kVARS%exner,            &
                             kVARS%water_vapor, kVARS%cloud_water,             kVARS%rain_in_air,      &
                             kVARS%u,           kVARS%v,                       kVARS%w,                &
                             kVARS%dz ])

            ! List the variables that are required to be advected
            ! call options%advect_vars( &
            !               [, &
            !                , &
            !                ] )

            ! List the variables that are required for restarts with the linear wind solution
            call options%restart_vars( &
                            [kVARS%nsquared,    kVARS%potential_temperature,                           &
                             kVARS%water_vapor, kVARS%cloud_water,             kVARS%rain_in_air,      &
                             kVARS%u,           kVARS%v,                       kVARS%w,                &
                             kVARS%dz ])

        end subroutine

        subroutine wind_var_request(options)
            implicit none
            type(options_t), intent(inout) :: options

            if (options%physics%windtype == kWIND_LINEAR) then
                call wind_linear_var_request(options)
            endif
            if (options%physics%windtype == kLINEAR_ITERATIVE_WINDS) then
                call wind_linear_var_request(options)
            endif

            if (options%physics%windtype == kITERATIVE_WINDS) then
                call options%restart_vars([kVARS%w_real])
            endif
        end subroutine wind_var_request




    !>------------------------------------------------------------
    !! Forces u,v, and w fields to balance
    !!       du/dx + dv/dy = dw/dz
    !!
    !! Starts by setting w out of the ground=0 then works through layers
    !!
    !!------------------------------------------------------------

    subroutine balance_uvw(domain,options,update_in)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        logical, optional, intent(in) :: update_in
        
        real, allocatable, dimension(:,:,:) :: divergence
        integer :: ims, ime, jms, jme, kms, kme, i_s, i_e, j_s, j_e, k
        logical :: update
        
        
        associate(dx         => domain%dx,                         &
                  rho        => domain%density%data_3d,            &
                  dz         => domain%advection_dz,               &
                  jaco_u     => domain%jacobian_u,                 &
                  jaco_v     => domain%jacobian_v,                 &
                  jaco_w     => domain%jacobian_w)

        ims = domain%ims
        ime = domain%ime
        kms = domain%kms
        kme = domain%kme
        jms = domain%jms
        jme = domain%jme
        
        i_s = domain%its-1
        i_e = domain%ite+1
        j_s = domain%jts-1
        j_e = domain%jte+1
        
        !Safety checks
        if (i_s<domain%ids) i_s = domain%ids
        if (i_e>domain%ide) i_e = domain%ide
        if (j_s<domain%jds) j_s = domain%jds
        if (j_e>domain%jde) j_e = domain%jde

        
        update = .False.
        if(present(update_in)) update=update_in
        
        allocate(divergence(ims:ime,kms:kme,jms:jme))

        if (update) then
            call calc_divergence(divergence,domain%u%meta_data%dqdt_3d,domain%v%meta_data%dqdt_3d,domain%w%meta_data%dqdt_3d, &
                                 jaco_u,jaco_v,jaco_w,dz,dx,rho,options,horz_only=.True.)
            call calc_w(domain%w%meta_data%dqdt_3d,divergence,dz,jaco_w,rho,ims,ime,kms,kme,jms,jme,i_s,i_e,j_s,j_e,&
                        options%parameters%advect_density)
        else
            call calc_divergence(divergence,domain%u%data_3d,domain%v%data_3d,domain%w%data_3d, &
                                 jaco_u,jaco_v,jaco_w,dz,dx,rho,options,horz_only=.True.)
            call calc_w(domain%w%data_3d,divergence,dz,jaco_w,rho,ims,ime,kms,kme,jms,jme,i_s,i_e,j_s,j_e,options%parameters%advect_density)
        endif
            
        end associate

        !call domain%w%exchange(update)
        
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

        
    end subroutine balance_uvw

    subroutine calc_w(w,div,dz,jaco_w,rho,ims,ime,kms,kme,jms,jme,i_s,i_e,j_s,j_e,adv_den)
        real,    intent(inout)                                   :: w(ims:ime,kms:kme,jms:jme)
        real,    dimension(ims:ime,kms:kme,jms:jme), intent(in)  :: div, dz, jaco_w, rho
        logical, intent(in)    :: adv_den
        integer, intent(in)    :: ims, ime, kms, kme, jms, jme, i_s, i_e, j_s, j_e
        
        real, allocatable, dimension(:,:,:) :: rho_i
        integer k

        allocate(rho_i(ims:ime,kms:kme-1,jms:jme))
        
        rho_i(:,kms:kme-1,:) = ( rho(:,kms:kme-1,:)*dz(:,kms:kme-1,:) + rho(:,kms+1:kme,:)*dz(:,kms+1:kme,:) ) / (dz(:,kms:kme-1,:)+dz(:,kms+1:kme,:))
        
        w = 0
        
        do k = kms,kme
            if (adv_den) then
                if (k==kms) then
                    w(i_s:i_e,k,j_s:j_e) = 0 - div(i_s:i_e,k,j_s:j_e) * dz(i_s:i_e,k,j_s:j_e) &
                                                / (jaco_w(i_s:i_e,k,j_s:j_e) * rho_i(i_s:i_e,k,j_s:j_e) )
                elseif (k==kme) then
                    w(i_s:i_e,k,j_s:j_e) = ((w(i_s:i_e,k-1,j_s:j_e) * rho_i(i_s:i_e,k-1,j_s:j_e) &
                                                * jaco_w(i_s:i_e,k-1,j_s:j_e)) - div(i_s:i_e,k,j_s:j_e) * &
                                                dz(i_s:i_e,k,j_s:j_e)) / (jaco_w(i_s:i_e,k,j_s:j_e) * rho(i_s:i_e,k,j_s:j_e))
                else
                    w(i_s:i_e,k,j_s:j_e) = ( (w(i_s:i_e,k-1,j_s:j_e) * rho_i(i_s:i_e,k-1,j_s:j_e) * &
                                jaco_w(i_s:i_e,k-1,j_s:j_e)) - div(i_s:i_e,k,j_s:j_e) * dz(i_s:i_e,k,j_s:j_e)) / &
                                (jaco_w(i_s:i_e,k,j_s:j_e) *  rho_i(i_s:i_e,k,j_s:j_e) )
                endif
            else
                if (k==kms) then
                    w(i_s:i_e,k,j_s:j_e) = (0 - div(i_s:i_e,k,j_s:j_e) * dz(i_s:i_e,k,j_s:j_e)) / (jaco_w(i_s:i_e,k,j_s:j_e) )
                else 
                    w(i_s:i_e,k,j_s:j_e) = (w(i_s:i_e,k-1,j_s:j_e) * jaco_w(i_s:i_e,k-1,j_s:j_e) - &
                                                div(i_s:i_e,k,j_s:j_e) * dz(i_s:i_e,k,j_s:j_e))/ (jaco_w(i_s:i_e,k,j_s:j_e) )
                end if
            end if
        end do

    end subroutine

    subroutine calc_divergence(div, u, v, w, jaco_u, jaco_v, jaco_w, dz, dx, rho, options, horz_only)
        implicit none
        real,           intent(inout) :: div(:,:,:)
        real,           intent(in)    :: u(:,:,:), v(:,:,:), w(:,:,:), dz(:,:,:), jaco_u(:,:,:), jaco_v(:,:,:), jaco_w(:,:,:), rho(:,:,:)

        real,           intent(in)    :: dx
        type(options_t),intent(in)    :: options
        logical, optional, intent(in) :: horz_only
        
        real, allocatable, dimension(:,:,:) :: diff_U, diff_V, u_met, v_met, w_met, rho_i
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
        allocate(rho_i(ims:ime,kms:kme-1,jms:jme))

        !Multiplication of U/V by metric terms, converting jacobian to staggered-grid where possible, otherwise making assumption of
        !Constant jacobian at edges
        
        if (options%parameters%advect_density) then
            u_met(ims+1:ime,:,:) = u(ims+1:ime,:,:) * jaco_u(ims+1:ime,:,:) * (rho(ims:ime-1,:,:) + rho(ims+1:ime,:,:))/2
            u_met(ims,:,:) = u(ims,:,:) * jaco_u(ims,:,:) * (1.5*rho(ims,:,:) - 0.5*rho(ims+1,:,:))
            u_met(ime+1,:,:) = u(ime+1,:,:) * jaco_u(ime+1,:,:) * (1.5*rho(ime,:,:) - 0.5*rho(ime-1,:,:))

            v_met(:,:,jms+1:jme) = v(:,:,jms+1:jme) * jaco_v(:,:,jms+1:jme) * (rho(:,:,jms:jme-1) + rho(:,:,jms+1:jme))/2
            v_met(:,:,jms) = v(:,:,jms) * jaco_v(:,:,jms) * (1.5*rho(:,:,jms) - 0.5*rho(:,:,jms+1))
            v_met(:,:,jme+1) = v(:,:,jme+1) * jaco_v(:,:,jme+1) * (1.5*rho(:,:,jme) - 0.5*rho(:,:,jme-1))

            rho_i(:,kms:kme-1,:) = ( rho(:,kms:kme-1,:)*dz(:,kms:kme-1,:) + rho(:,kms+1:kme,:)*dz(:,kms+1:kme,:) ) / (dz(:,kms:kme-1,:)+dz(:,kms+1:kme,:))

        else
            u_met = u * jaco_u
            v_met = v * jaco_v
        end if

        diff_U = u_met(ims+1:ime+1, :, jms:jme) - u_met(ims:ime, :, jms:jme)
        diff_V = v_met(ims:ime, :, jms+1:jme+1) - v_met(ims:ime, :, jms:jme)

        div(ims:ime,kms:kme,jms:jme) = (diff_U+diff_V) /(dx)

        if (.NOT.(horz)) then
            if (options%parameters%advect_density) then
                w_met(:,kme,:) = w(:,kme,:) * jaco_w(:,kme,:) * rho(:,kme,:)
                w_met(:,kms:kme-1,:) = w(:,kms:kme-1,:) * jaco_w(:,kms:kme-1,:) * rho_i
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
        endif

    end subroutine calc_divergence
    


    !>------------------------------------------------------------
    !! Correct for a grid that is locally rotated with respect to EW,NS
    !!
    !! Assumes forcing winds are EW, NS relative, not grid relative.
    !!
    !!------------------------------------------------------------
    subroutine make_winds_grid_relative(u, v, sintheta, costheta, ims, ime, kms, kme, jms, jme)
        real, intent(inout)             :: u(ims:ime+1,kms:kme,jms:jme), v(ims:ime,kms:kme,jms:jme+1)
        double precision, intent(in)    :: sintheta(ims:ime,jms:jme), costheta(ims:ime,jms:jme)
        integer, intent(in)             :: ims, ime, kms, kme, jms, jme
        real, dimension(:,:,:), allocatable :: u_vstag,v_ustag
        real, dimension(:,:), allocatable :: costheta_ustag,sintheta_ustag,costheta_vstag,sintheta_vstag

        integer :: k
        
        allocate(v_ustag(ims+1:ime,kms:kme,jms:jme))
        allocate(u_vstag(ims:ime,kms:kme,jms+1:jme))
        
        allocate(costheta_ustag(ims+1:ime,jms:jme))
        allocate(sintheta_ustag(ims+1:ime,jms:jme))
        allocate(costheta_vstag(ims:ime,jms+1:jme))
        allocate(sintheta_vstag(ims:ime,jms+1:jme))

        
        v_ustag = (v(ims:ime-1,:,jms:jme)+v(ims+1:ime,:,jms:jme)+v(ims:ime-1,:,jms+1:jme+1)+v(ims+1:ime,:,jms+1:jme+1))/4
        u_vstag = (u(ims:ime,:,jms:jme-1)+u(ims:ime,:,jms+1:jme)+u(ims+1:ime+1,:,jms:jme-1)+u(ims+1:ime+1,:,jms+1:jme))/4
        
        costheta_ustag = (costheta(ims+1:ime,jms:jme)+costheta(ims:ime-1,jms:jme))/2
        sintheta_ustag = (sintheta(ims+1:ime,jms:jme)+sintheta(ims:ime-1,jms:jme))/2
        
        costheta_vstag = (costheta(ims:ime,jms+1:jme)+costheta(ims:ime,jms:jme-1))/2
        sintheta_vstag = (sintheta(ims:ime,jms+1:jme)+sintheta(ims:ime,jms:jme-1))/2

        do k = kms,kme
            u(ims+1:ime,k,:) = u(ims+1:ime,k,:) * costheta_ustag - v_ustag(:,k,:) * sintheta_ustag
            v(:,k,jms+1:jme) = v(:,k,jms+1:jme) * costheta_vstag + u_vstag(:,k,:) * sintheta_vstag

            !u(ims,k,:)       = u(ims,k,:) * costheta_ustag(ims+1,:) + v_ustag(ims+1,k,:) * sintheta_ustag(ims+1,:)
            !u(ime+1,k,:)     = u(ime,k,:) * costheta_ustag(ime,:) + v_ustag(ime,k,:) * sintheta_ustag(ime,:)
        
            !v(:,k,jms)       = v(:,k,jms) * costheta_vstag(:,jms+1) + u_vstag(:,k,jms+1) * sintheta_vstag(:,jms+1)
            !v(:,k,jme+1)     = v(:,k,jme) * costheta_vstag(:,jme) + u_vstag(:,k,jme) * sintheta_vstag(:,jme)
        enddo
        
        deallocate(v_ustag,u_vstag,costheta_ustag,sintheta_ustag,costheta_vstag,sintheta_vstag)

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

        real, allocatable, dimension(:,:,:) :: alpha, div
        integer :: nx, ny, nz, i, j
        integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte


        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)
        
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte
        
        if (.not.associated(domain%w%meta_data%dqdt_3d)) then

            call init_winds(domain, options)
            !call initialize_blocking(domain, options)
            call update_stability(domain)

            ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
            call make_winds_grid_relative(domain%u%data_3d, domain%v%data_3d, domain%sintheta, domain%costheta,ims, ime, kms, kme, jms, jme)

            call domain%u%exchange_x()
            call domain%v%exchange_y()

            ! flow blocking parameterization
            ! if (options%block_options%block_flow) then
            !     call add_blocked_flow(domain, options)
            ! endif

            if (options%wind%Sx) then
                call apply_Sx(domain%Sx,domain%TPI,domain%u%data_3d, domain%v%data_3d, domain%w%data_3d,domain%Ri,domain%dzdx,domain%dzdy,domain%z%data_3d)
            endif 

            ! linear winds
            if (options%physics%windtype==kWIND_LINEAR) then
                call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%parameters%advect_density, update=.False.)
            ! simple acceleration over topography
            elseif (options%physics%windtype==kCONSERVE_MASS) then
                if (options%parameters%use_terrain_difference) then  !
                !! use the ratio between hi-res and lo-res grid deformation (i.e. due to 'additional' terrain) for speedup
                    call mass_conservative_acceleration(domain%u%data_3d, domain%v%data_3d, domain%zfr_u, domain%zfr_v)
                else
                    call mass_conservative_acceleration(domain%u%data_3d, domain%v%data_3d, domain%zr_u, domain%zr_v)
                endif
            elseif (options%physics%windtype==kOBRIEN_WINDS) then
                call Obrien_winds(domain, options, its, ite, jts, jte)
            elseif (options%physics%windtype==kITERATIVE_WINDS) then
            
                allocate(alpha(ims:ime,kms:kme,jms:jme))
                allocate(div(ims:ime,kms:kme,jms:jme))
                
                call calc_alpha(alpha, domain%froude, ims, ime,  kms, kme, jms, jme, &
                            domain%ids, domain%ide, domain%jds, domain%jde)

                !Call this, passing 0 for w_grid, to get vertical components of vertical motion
                call calc_w_real(domain% u %data_3d,      &
                             domain% v %data_3d,      &
                             domain%w%data_3d*0.0,      &
                             domain%w%data_3d,      &
                             domain%dzdx_u, domain%dzdy_v,    &
                             domain%jacobian_w,ims,ime,kms,kme,jms,jme,its, ite, jts, jte)

                !If we have not read in W_real from forcing, set target w_real to 0.0. This minimizes vertical motion in solution
                if (options%parameters%wvar=="") domain%w_real%data_3d = 0.0  
                domain%w%data_3d = (domain%w_real%data_3d-domain%w%data_3d)/domain%jacobian
                call calc_divergence(div,domain%u%data_3d,domain%v%data_3d,domain%w%data_3d, &
                                domain%jacobian_u, domain%jacobian_v,domain%jacobian_w,domain%advection_dz,domain%dx, &
                                domain%density%data_3d,options,horz_only=.False.)

                call calc_iter_winds(domain,alpha,div,options%parameters%advect_density)


            elseif (options%physics%windtype==kLINEAR_ITERATIVE_WINDS) then
                call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%parameters%advect_density, update=.False.)
                call iterative_winds(domain, options)

            endif
            ! else assumes even flow over the mountains

            ! use horizontal divergence (convergence) to calculate vertical convergence (divergence)

            call balance_uvw(domain,options)
            
            call calc_w_real(domain% u %data_3d,      &
                             domain% v %data_3d,      &
                             domain% w %data_3d,      &
                             domain% w_real %data_3d,      &
                             domain%dzdx_u, domain%dzdy_v,    &
                             domain%jacobian_w,ims,ime,kms,kme,jms,jme,its, ite, jts, jte)
        else

            call update_stability(domain)

            ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
            call make_winds_grid_relative(domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, domain%sintheta, domain%costheta,ims, ime, kms, kme, jms, jme)
            
            call domain%u%exchange_x(do_metadata=.True.)
            call domain%v%exchange_y(do_metadata=.True.)
            
            if (options%wind%Sx) then
                call apply_Sx(domain%Sx,domain%TPI,domain%u%meta_data%dqdt_3d,domain%v%meta_data%dqdt_3d, domain%w%meta_data%dqdt_3d,domain%Ri,domain%dzdx,domain%dzdy,domain%z%data_3d)
            endif 

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
            elseif (options%physics%windtype==kOBRIEN_WINDS) then
                call Obrien_winds(domain, options, its, ite, jts, jte, update_in=.True.)
            elseif (options%physics%windtype==kITERATIVE_WINDS) then
            
                allocate(alpha(ims:ime,kms:kme,jms:jme))
                allocate(div(ims:ime,kms:kme,jms:jme))

                call calc_alpha(alpha, domain%froude, ims, ime,  kms, kme, jms, jme, &
                            domain%ids, domain%ide, domain%jds, domain%jde)

                !Call this, passing 0 for w_grid, to get vertical components of vertical motion
                call calc_w_real(domain% u %meta_data%dqdt_3d,      &
                             domain% v %meta_data%dqdt_3d,      &
                             domain%w%meta_data%dqdt_3d*0.0,      &
                             domain%w%meta_data%dqdt_3d,      &
                             domain%dzdx_u, domain%dzdy_v,    &
                             domain%jacobian_w,ims,ime,kms,kme,jms,jme,its, ite, jts, jte)
                             
                !If we have not read in W_real from forcing, set target w_real to 0.0. This minimizes vertical motion in solution
                if (options%parameters%wvar=="") then
                    domain%w%meta_data%dqdt_3d = (0.0-domain%w%meta_data%dqdt_3d)/domain%jacobian
                else
                    domain%w%meta_data%dqdt_3d = (domain%w_real%dqdt_3d-domain%w%meta_data%dqdt_3d)/domain%jacobian
                endif

                call calc_divergence(div,domain%u%meta_data%dqdt_3d,domain%v%meta_data%dqdt_3d,domain%w%meta_data%dqdt_3d, &
                                domain%jacobian_u, domain%jacobian_v,domain%jacobian_w,domain%advection_dz,domain%dx, &
                                domain%density%data_3d,options,horz_only=.False.)
                                
                call calc_iter_winds(domain,alpha,div,options%parameters%advect_density,update_in=.True.)

            elseif (options%physics%windtype==kLINEAR_ITERATIVE_WINDS) then
                call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%parameters%advect_density, update=.True.)
                call iterative_winds(domain, options, update_in=.True.)
            endif

            call balance_uvw(domain,options,update_in=.True.)
                             
            call calc_w_real(domain% u %meta_data%dqdt_3d,      &
                             domain% v %meta_data%dqdt_3d,      &
                             domain% w %meta_data%dqdt_3d,      &
                             domain% w_real %dqdt_3d,           &
                             domain%dzdx_u, domain%dzdy_v,    &
                             domain%jacobian_w,ims,ime,kms,kme,jms,jme,its, ite, jts, jte)


        endif

    end subroutine update_winds
    
    subroutine calc_alpha(alpha, froude, ims, ime, kms, kme, jms, jme, ids, ide, jds, jde)
        implicit none
        integer, intent(in)    :: ims, ime, kms, kme, jms, jme, ids, ide, jds, jde
        real,    intent(in)    :: froude(ims:ime,kms:kme,jms:jme)
        real,    intent(inout) :: alpha(ims:ime,kms:kme,jms:jme)
        
        alpha = 1.0
        !Following Moussiopoulos, et al. (1988). Bounding low Fr to avoid /0 error and negative Fr
        alpha = 1.0 - 0.5*(1./max(froude**4,0.00001))*(sqrt(1.0+4.0*max(froude**4,0.00001)) - 1.0) 
        alpha = sqrt(max(alpha,0.0))
        alpha = min(max(alpha,0.1),1.0)

        ! Ensure that there are no sharp transitions in alpha at boundary, 
        ! which can leak boundary effects into model (very high w_grid values result)
        if (jms==jds) alpha(:,:,jms) = alpha(:,:,jms+1)
        if (ims==ids) alpha(ims,:,:) = alpha(ims+1,:,:)
        if (jme==jde) alpha(:,:,jme) = alpha(:,:,jme-1)
        if (ime==ide) alpha(ime,:,:) = alpha(ime-1,:,:)
    end subroutine calc_alpha
    
    subroutine calc_w_real(u,v,w_grid,w_real,dzdx,dzdy,jaco_w,ims,ime,kms,kme,jms,jme,its, ite, jts, jte)

        implicit none
        real, intent(in), dimension(ims:ime,kms:kme,jms:jme)   :: w_grid,jaco_w
        real, intent(in), dimension(ims:ime+1,kms:kme,jms:jme) :: u,dzdx
        real, intent(in), dimension(ims:ime,kms:kme,jms:jme+1) :: v,dzdy
        real, intent(inout)                                    :: w_real(ims:ime,kms:kme,jms:jme)
        
        real, allocatable :: lastw(:,:)
        real, allocatable :: currw(:,:)
        real, allocatable :: uw(:,:)
        real, allocatable :: vw(:,:)
        integer :: z, ims, ime, kms, kme, jms, jme, its, ite, jts, jte

        if (.not.allocated(lastw)) then
            allocate( lastw( its:ite  , jts:jte  ))
            allocate( currw( its:ite  , jts:jte  ))
            allocate(    uw( its:ite+1, jts:jte  ))
            allocate(    vw( its:ite  , jts:jte+1))
        endif
        
        !calculate the real vertical motions (including U*dzdx + V*dzdy)
        lastw = 0
        do z = kms, kme
            
            ! ! if(options%parameters%use_terrain_difference) then
            !                 ! compute the U * dz/dx component of vertical motion
            !     uw    = u(its:ime,   z, jts:jts) * domain%delta_dzdx(:,z,jts:jts)
            !     ! compute the V * dz/dy component of vertical motion
            !     vw    = v(its:its, z, jts:jme  ) * domain%delta_dzdy(its:its,z,:)
            ! else    
                ! compute the U * dz/dx component of vertical motion
                uw    = u(its:ite+1,   z, jts:jte) * dzdx(its:ite+1,z,jts:jte)
                ! compute the V * dz/dy component of vertical motion
                vw    = v(its:ite, z,   jts:jte+1) * dzdy(its:ite,z,jts:jte+1)
            ! endif    
            ! ! convert the W grid relative motion to m/s
            ! currw = w(its:its, z, jts:jts) * dz_interface(its:its, z, jts:jts) / domain%dx

            ! the W grid relative motion
            currw = w_grid(its:ite, z, jts:jte) * jaco_w(its:ite, z, jts:jte)

            ! if (options%physics%convection>0) then
            !     currw = currw + domain%w_cu(2:nx-1,z,2:ny-1) * domain%dz_inter(2:nx-1,z,2:ny-1) / domain%dx
            ! endif
            
            ! compute the real vertical velocity of air by combining the different components onto the mass grid
            ! includes vertical interpolation between w_z-1/2 and w_z+1/2
            w_real(its:ite, z, jts:jte) = (uw(its:ite,:) + uw(its+1:ite+1,:))*0.5 &
                                                 +(vw(:,jts:jte) + vw(:,jts+1:jte+1))*0.5 &
                                                 +(lastw + currw) * 0.5
            lastw = currw ! could avoid this memcopy cost using pointers or a single manual loop unroll
        end do
                
    end subroutine calc_w_real
    
    subroutine Obrien_winds(domain, options, its, ite, jts, jte, update_in)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        logical, optional, intent(in) :: update_in

        ! interal parameters
        real, allocatable, dimension(:,:,:) :: div, ADJ,ADJ_coef, current_u, current_v, current_w
        real    :: corr_factor
        integer :: it, k, j, i, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, wind_k
        logical :: update

        update=.False.
        if (present(update_in)) update=update_in

        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)
        
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte
        
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
        call domain%u%exchange_x()
        call domain%v%exchange_y()

        !First call bal_uvw to generate an initial-guess for vertical winds
        call balance_uvw(domain,options)

        allocate(div(ims:ime,kms:kme,jms:jme))
        allocate(ADJ_coef(ims:ime,kms:kme,jms:jme))
        allocate(ADJ(ims:ime,kms:kme,jms:jme))
        
        ! Calculate and apply correction to w-winds
        wind_k = kme
        ! previously this code was solving for 0 vertical motion at the flat z height instead of the top boundary.
        ! left in for now as it could be useful to implement something similar in the future.
        ! however, this was also creating weird artifacts above the flat z height that need to be fixed if re-implementing.
        ! do k = kms,kme
        !     if (sum(domain%advection_dz(ims,1:k,jms)) > domain%smooth_height) then
        !         wind_k = k
        !         exit
        !     endif
        ! enddo
        ! domain%smooth_height = sum(domain%advection_dz(ims,:,jms))
        !Compute relative correction factors for U and V based on input speeds
        U_cor = ABS(domain%u%data_3d(ims:ime,:,jms:jme))/ &
                (ABS(domain%u%data_3d(ims:ime,:,jms:jme))+ABS(domain%v%data_3d(ims:ime,:,jms:jme)))

        do i = ims,ime
            do j = jms,jme
                domain%smooth_height = sum(domain%advection_dz(i,:,j)) !
                do k = kms,kme
                    corr_factor = ((sum(domain%advection_dz(i,1:k,j)))/domain%smooth_height)
                    corr_factor = min(corr_factor,1.0)
                    domain%w%data_3d(i,k,j) = domain%w%data_3d(i,k,j) - corr_factor * (domain%w%data_3d(i,wind_k,j))

                    !if ( (domain%u%data_3d(i,k,j)+domain%v%data_3d(i,k,j)) == 0) U_cor(i,k,j) = 0.5
                enddo
            enddo
        enddo

        do k = kms,kme
            ! Compute this now, since it wont change in the loop
            ADJ_coef(:,k,:) = -2/domain%dx
        enddo

        ! Now, fixing w-winds, iterate over U/V to reduce divergence with new w-winds
        do it = 0,options%parameters%wind_iterations
            !Compute divergence in new wind field
            call calc_divergence(div,domain%u%data_3d,domain%v%data_3d,domain%w%data_3d, &
                                domain%jacobian_u, domain%jacobian_v,domain%jacobian_w,domain%advection_dz,domain%dx, &
                                domain%density%data_3d,options)
                                
            !Because we want the true divergence, and the above function does not divide by the jacobian in the end, divide here
            div = div/domain%jacobian

            !Compute adjustment based on divergence
            ADJ = div/ADJ_coef

            !Distribute divergence among the U and V fields
            domain%u%data_3d(its+1:ite+1,:,jts:jte) = domain%u%data_3d(its+1:ite+1,:,jts:jte) + &
                                                        (ADJ(its:ite,:,jts:jte) * 0.5)

            domain%u%data_3d(its+1:ite+1,:,jts:jte) = domain%u%data_3d(its+1:ite+1,:,jts:jte) - &
                                                        (ADJ(its+1:ite+1,:,jts:jte) * 0.5)

            domain%v%data_3d(its:ite,:,jts+1:jte+1) = domain%v%data_3d(its:ite,:,jts+1:jte+1) + &
                                                        (ADJ(its:ite,:,jts:jte) * 0.5)

            domain%v%data_3d(its:ite,:,jts+1:jte+1) = domain%v%data_3d(its:ite,:,jts+1:jte+1) - &
                                                        (ADJ(its:ite,:,jts+1:jte+1) * 0.5)
            call domain%u%exchange_x()
            call domain%v%exchange_y()

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

    end subroutine Obrien_winds

    
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

        do i=domain%grid%kms, domain%grid%kme
            domain%advection_dz(:,i,:) = options%parameters%dz_levels(i)
        enddo
        
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

        real, allocatable, dimension(:,:,:) :: wind_speed, temp_froude, u_m, v_m, u_shear, v_shear, winddir, stability
        integer,  allocatable, dimension(:,:,:) :: dir_indices
        
        integer :: k, j, i, n, ims, ime, jms, jme, kms, kme, ob_k
        real :: z_top, z_bot, th_top, th_bot, obstacle_height
        integer :: ymin, ymax, xmin, xmax, n_smoothing_passes, nsmooth_gridcells
        

        n_smoothing_passes = 5
        nsmooth_gridcells = 20 !int(500 / domain%dx)
        
        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)

        !If it is our first time calculating Fr, allocate and populate froude_terrain array
        if (.not.allocated(domain%froude)) then
            allocate(domain%froude(ims:ime,kms:kme,jms:jme))            
            allocate(domain%froude_terrain(1:72,ims:ime,kms:kme,jms:jme))

            call compute_terrain_blocking_heights(domain)
        endif
       
        if (.not.allocated(domain%Ri)) allocate(domain%Ri(ims:ime,kms:kme,jms:jme))
       
        allocate(u_m(ims:ime,kms:kme,jms:jme))
        allocate(v_m(ims:ime,kms:kme,jms:jme))
        allocate(u_shear(ims:ime,kms:kme,jms:jme))
        allocate(v_shear(ims:ime,kms:kme,jms:jme))
        allocate(winddir(ims:ime,kms:kme,jms:jme))
        allocate(dir_indices(ims:ime,kms:kme,jms:jme))
        allocate(wind_speed(ims:ime,kms:kme,jms:jme))
        allocate(temp_froude(ims:ime,kms:kme,jms:jme))       
        allocate(stability(ims:ime,kms:kme,jms:jme))  
        
        u_m = (domain%u%data_3d(ims:ime,:,:) + domain%u%data_3d(ims+1:ime+1,:,:))/2
        v_m = (domain%v%data_3d(:,:,jms:jme) + domain%v%data_3d(:,:,jms+1:jme+1))/2
        
        u_shear(:,kms,:) = u_m(:,kms+4,:)
        u_shear(:,kms+1:kme,:) = u_m(:,kms+1:kme,:) - u_m(:,kms:kme-1,:)
        v_shear(:,kms,:) = v_m(:,kms+4,:)
        v_shear(:,kms+1:kme,:) = v_m(:,kms+1:kme,:) - v_m(:,kms:kme-1,:)
        
        !Compute wind direction for each cell on mass grid
        winddir = atan2(-u_m,-v_m)*rad2deg
        where(winddir < 0.0) winddir = winddir+360
        where(winddir == 360.0) winddir = 0.0
        dir_indices = int(winddir/5)+1

        !Build grid of Sx values based on wind direction at that cell
        do i = ims, ime
            do j = jms, jme
                do k=kms, kme
                    temp_froude(i,k,j) = domain%froude_terrain(dir_indices(i,k,j),i,k,j)
                enddo
            end do
        end do
        
        wind_speed = sqrt( (u_m)**2 + (v_m)**2 )
        
        !Since we will loop up to nz-1, we set all Fr to 0.1, which will leave the upper layer as very stable
        domain%froude = 0.1
        
        !Since we will loop up to nz-1, we set all Ri here to 10
        domain%Ri = 10.0
        
        do i = ims,ime
            do j = jms,jme
                do k = kms,kme-1
                    th_bot = domain%potential_temperature%data_3d(i,kms,j)
                    th_top = domain%potential_temperature%data_3d(i,kms+4,j)
                    z_bot  = domain%z%data_3d(i,kms,j)
                    z_top  = domain%z%data_3d(i,kms+4,j)
                    stability(i,k,j) = calc_dry_stability(th_top, th_bot, z_top, z_bot) 
                    
                    domain%Ri(i,k,j) =  calc_Ri(stability(i,k,j), u_shear(i,kms,j), v_shear(i,kms,j), (z_top-z_bot))
                    
                    th_bot = domain%potential_temperature%data_3d(i,k,j)
                    th_top = domain%potential_temperature%data_3d(i,k+1,j)
                    z_bot  = domain%z%data_3d(i,k,j)
                    z_top  = domain%z%data_3d(i,k+1,j)
                    stability(i,k,j) = calc_dry_stability(th_top, th_bot, z_top, z_bot) 
                    
                    ! If we have an upwind obstacle, use the obstacle height to calculate a bulk Froude Number over the column
                    ! If there is nothing blocking, we calculate a local bulk Froude Number, using the local th and z indexed 
                    ! above
                    if (.not.(temp_froude(i,k,j) == DEFAULT_FR_L)) then
                    !The height of the obstacle is calculated from the blocking terrain height (z_obst-z_loc+1000)
                        obstacle_height = temp_froude(i,k,j)-DEFAULT_FR_L+z_bot

                        do ob_k = k+1,kme
                            if (domain%z%data_3d(i,ob_k,j) > obstacle_height) exit
                        enddo
                        ob_k = min(ob_k,kme)
                        th_top = domain%potential_temperature%data_3d(i,ob_k,j)
                        z_top  = domain%z%data_3d(i,ob_k,j)

                        stability(i,k,j) = calc_dry_stability(th_top, th_bot, z_top, z_bot) 

                    endif

                    !Above function calculates N^2, but froude number wants just N
                    stability(i,k,j) = sqrt(max(stability(i,k,j), 0.0))
                    domain%froude(i,k,j) = calc_froude(stability(i,k,j), temp_froude(i,k,j), wind_speed(i,k,j))
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
                        !write(*,*) "num_sum:  ", ((xmax-xmin+1) * (ymax-ymin+1))
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
        real, allocatable    ::  azm(:,:), temp_ft_array(:,:,:,:)
        integer, allocatable :: azm_indices(:,:)
        integer           :: i, j, k, kms, kme, ang, i_s, j_s, i_start_buffer, i_end_buffer, j_start_buffer, j_end_buffer
        integer           :: rear_ang, fore_ang, test_ang, rear_ang_diff, fore_ang_diff, ang_diff, k_max, window_rear, window_fore, window_width
        integer :: nx, ny, x, y
        integer :: xs,xe, ys,ye, n, np, search_max
        real, allocatable :: temp_terrain(:,:), f_terrain(:,:)
        real              :: pt_height, temp_ft, maxFTVal

        search_max = int(max(4000.0/domain%dx,1.0))

        nx = size(domain%global_terrain,1)
        ny = size(domain%global_terrain,2)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        
        allocate(temp_ft_array( 1:72, domain%grid2d%ims:domain%grid2d%ime, kms:kme, domain%grid2d%jms:domain%grid2d%jme ))
        allocate(f_terrain(nx,ny))
        
        
        temp_ft_array = -100000.0
        
        allocate(azm( 2*search_max+1, 2*search_max+1 ))
        allocate(azm_indices( 2*search_max+1, 2*search_max+1 ))

        azm = 0
                
        !Setup azm so that it is looking along wind direction (i.e. negatives here)
        do i = 1, 2*search_max+1
            do j = 1, 2*search_max+1                
                azm(i,j)  = atan2(-1.0*(i-(search_max+1)),-1.0*(j-(search_max+1)))
            end do
        end do
        
        !convert azm to deg
        azm = azm*rad2deg
        where(azm < 0) azm = 360+azm
        where(azm >= 360.0) azm=0.0
        azm_indices = int(azm/5)+1

        ! then compute the range of terrain (max-min) in a given window
        do i=domain%grid2d%ims, domain%grid2d%ime
            do j=domain%grid2d%jms, domain%grid2d%jme
                do k=kms,kme
                    if (k == kms) then
                        pt_height = domain%global_terrain(i,j)
                    else if (k > kms) then
                        pt_height = pt_height + domain%global_dz_interface(i,k,j)
                    end if
                    
                    
                    ! Check to use buffers to avoid searching out of grid
                    i_start_buffer = -min(0,i-(search_max+1))
                    i_end_buffer = min(0,domain%grid2d%ide-(i+search_max))
                
                    j_start_buffer = -min(0,j-(search_max+1))
                    j_end_buffer = min(0,domain%grid2d%jde-(j+search_max))
                
                    do i_s = 1+i_start_buffer, (search_max*2+1)+i_end_buffer
                        do j_s = 1+j_start_buffer, (search_max*2+1)+j_end_buffer
                        
                            temp_ft = DEFAULT_FR_L + domain%global_terrain(i+(i_s-(search_max+1)),j+(j_s-(search_max+1))) - pt_height
                            
                            if (temp_ft > temp_ft_array(azm_indices(i_s,j_s),i,k,j)) then
                                                        
                                !Only save scale length if it is greater than the default -- otherwise copy that over
                                if (temp_ft > DEFAULT_FR_L) then
                                    temp_ft_array(azm_indices(i_s,j_s),i,k,j) = temp_ft
                                else
                                    temp_ft_array(azm_indices(i_s,j_s),i,k,j) = DEFAULT_FR_L
                                end if
                            end if
                        enddo
                    enddo

                    !After finding Fr-Terrain in each absolute direction around grid cell, 
                    !Pick max for each 20º window and perform interpolation to other directions if necesarry
                    
                    rear_ang = 1 
                    fore_ang = 1
                    
                    if (.not.( all((temp_ft_array(:,i,k,j) <= -100000.0)) )) then
                    
                        !Perform 20º window max search
                        window_width = 2
                        do ang = 1, 72
                            window_rear = ang-window_width
                            window_fore = ang+window_width
                        
                            if (ang <= window_width) then
                                window_rear = 72-(window_width-ang)
                                
                                maxFTVal = maxval(temp_ft_array(window_rear:72,i,k,j))

                                if (maxval(temp_ft_array(1:window_fore,i,k,j)) > maxFTVal) then
                                    maxFTVal = maxval(temp_ft_array(1:window_fore,i,k,j))
                                end if
                                
                            else if ( ang >= (72-(window_width-1)) ) then
                                window_fore = window_width-(72-ang)
                                
                                maxFTVal = maxval(temp_ft_array(window_rear:72,i,k,j))

                                if (maxval(temp_ft_array(1:window_fore,i,k,j)) > maxFTVal) then
                                    maxFTVal = maxval(temp_ft_array(1:window_fore,i,k,j))
                                end if
                            else
                                maxFTVal = maxval(temp_ft_array(window_rear:window_fore,i,k,j))
                            end if
                            domain%froude_terrain(ang,i,k,j) = maxFTVal
                        end do                    
                    
                        do ang = 1, 72
                            !Determine indices for interpolation
                            if ( (ang==fore_ang) ) then
                                !Update indices for interpolated Fr-Terrain's
                                rear_ang = ang
                            
                                fore_ang = ang+1
                                if (fore_ang > 72) fore_ang = 1
                                
                                do while (domain%froude_terrain(fore_ang,i,k,j) <= -100000.0)
                                    fore_ang = fore_ang+1
                                    if (fore_ang > 72) fore_ang = 1
                                end do
                            
                            end if
                            
                            if (ang==1) then
                                rear_ang = 72
                                do while(domain%froude_terrain(rear_ang,i,k,j) <= -100000.0)
                                    rear_ang = rear_ang-1
                                end do
                            end if
                    
                            !If we did not calculate Fr-Terrain for a given direction
                            if (domain%froude_terrain(ang,i,k,j) == -100000.0) then
                                !Weight the two surrounding Fr-Terrain values based on our angular-distance to them
                                rear_ang_diff = ang-rear_ang
                                fore_ang_diff = fore_ang-ang
                                ang_diff = fore_ang-rear_ang
                        
                                !Handle wrap-around case
                                if (ang > fore_ang) then
                                    fore_ang_diff = fore_ang+(72-ang)
                                    ang_diff = fore_ang+(72-rear_ang)
                                end if
                        
                                !Interpolation, linearly-weighted by angular-distance from values
                                domain%froude_terrain(ang,i,k,j) = (domain%froude_terrain(rear_ang,i,k,j)*fore_ang_diff + &
                                                    domain%froude_terrain(fore_ang,i,k,j)*rear_ang_diff)/ang_diff

                            end if
                        end do

                    else
                        !IF we only have -100000 for all entries, set to dz
                        domain%froude_terrain(:,i,k,j) = DEFAULT_FR_L
                    end if
                enddo

            enddo
        enddo
                                                               
        if (domain%jts==(domain%jds+1)) domain%froude_terrain(:,:,:,domain%grid2d%jms) = &
                                        domain%froude_terrain(:,:,:,domain%grid2d%jms+1)
                        
        if (domain%its==(domain%ids+1)) domain%froude_terrain(:,domain%grid2d%ims,:,:) = &
                                        domain%froude_terrain(:,domain%grid2d%ims+1,:,:)

        if (domain%jte==(domain%jde-1)) domain%froude_terrain(:,:,:,domain%grid2d%jme) = &
                                        domain%froude_terrain(:,:,:,domain%grid2d%jme-1)

        if (domain%ite==(domain%ide-1)) domain%froude_terrain(:,domain%grid2d%ime,:,:) = &
                                        domain%froude_terrain(:,domain%grid2d%ime-1,:,:)
                                 

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
