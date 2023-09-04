!>------------------------------------------------------------
!! Module to pre-compute terrain descriptors for thermal winds
!! and then apply a thermal wind correction based on Mahrt 1982
!! and Scire and Robe 2000
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!------------------------------------------------------------
module wind_thermal
    use data_structures
    use domain_interface,  only : domain_t
    use options_interface, only : options_t
    use io_routines,       only : io_write
    use mod_wrf_constants, only : wrf_cp, epsilon, piconst
    implicit none
    private
    public :: apply_thermal_winds, init_thermal_winds
    real, parameter :: thrd=1.0/3.0
    real, parameter :: L_e = 1000.0 ! Scale length for downslope flows (m)
    real, parameter :: Max_flow_height = 400 ! Maximum height for thermal flows above the surface (m)
    real, allocatable :: level_height(:,:,:)
    
    integer, save   :: therm_k_max = 0          
    
    integer, save :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer, save :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer, save :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions

contains
    
    subroutine init_thermal_winds(domain, options)
        implicit none
        type(domain_t),  intent(in) :: domain
        type(options_t), intent(in) :: options

        real :: z_mean
        integer :: k

        ids = domain%grid%ids
        ide = domain%grid%ide
        jds = domain%grid%jds
        jde = domain%grid%jde
        kds = domain%grid%kds
        kde = domain%grid%kde
        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        kms = domain%grid%kms
        kme = domain%grid%kme
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte
        kts = domain%grid%kts
        kte = domain%grid%kte

        !Max flow height
        if (therm_k_max==0) then
            do k = kms,kme
                z_mean = SUM(options%parameters%dz_levels(1:k))
                if (z_mean > Max_flow_height .and. therm_k_max==0) therm_k_max = max(2,k-1)
            enddo
            allocate( level_height(ims:ime,kms:therm_k_max,jms:jme))
            do k = kms,therm_k_max
                level_height(:,k,:) = domain%z%data_3d(:,k,:)-domain%z_interface%data_3d(ims:ime,kms,jms:jme)
            enddo
        endif

    end subroutine
    
    subroutine apply_thermal_winds(hfss, rho, temperature, r_dist, v_dist, r_drop, u, v, dzdx, dzdy, z_0) !v_dist, r_drop, v_drop,
        implicit none
        real, intent(in), dimension(ims:ime,kms:kme,jms:jme)      :: rho, temperature, dzdx, dzdy
        real, intent(in), dimension(ims:ime,jms:jme)              :: hfss, r_dist, v_dist, r_drop, z_0
        real, intent(inout), dimension(ims:ime+1,kms:kme,jms:jme) :: u
        real, intent(inout), dimension(ims:ime,kms:kme,jms:jme+1) :: v

        real, dimension(ims:ime,jms:jme)   :: x_norm, y_norm, slope, slp_wnd_spd, spd, C_d, upslope_peak_h, upslope_max_h, downslope_peak_h, downslope_max_h, up_z_scale, dwn_z_scale
        real, dimension(ims:ime,kms:therm_k_max,jms:jme) :: therm_U_corr, therm_V_corr

        integer ::  i, j, k

        do j = jms,jme
            do i = ims,ime
                x_norm(i,j) = 0
                y_norm(i,j) = 0
                slp_wnd_spd(i,j) = 0
            enddo
        enddo
        
        do j = jms,jme
            do k = kms,therm_k_max
                do i = ims,ime
                    therm_U_corr(i,k,j) = 0
                    therm_V_corr(i,k,j) = 0
                enddo
            enddo
        enddo
                
        !Draf coef for upslope or downslope flows -- taken (roughly) from Mahrt 1982
        spd = sqrt((u(ims+1:ime+1,kms,:)*0.5+u(ims:ime,kms,:)*0.5)**2 + (v(:,kms,jms+1:jme+1)*0.5+v(:,kms,jms:jme)*0.5)**2)
        C_d = 0.41/(log(level_height(:,kms,:)/z_0)*spd+epsilon)
        
        !Calculate horizontal normal vector components of terrain
        slope = sqrt( dzdx(ims:ime,kms,jms:jme)**2+dzdy(ims:ime,kms,jms:jme)**2)
        where(.not.(slope == 0.)) x_norm = dzdx(ims:ime,kms,jms:jme)/slope
        where(.not.(slope == 0.)) y_norm = dzdy(ims:ime,kms,jms:jme)/slope
        slope = atan(slope)
        
        !Don't allow for crazy large slope flows/valley flows just because we end up far from a ridge
        !r_dist = max(r_dist,10000)
        
        ! Solve for maximum wind perturbation from slope flows
        ! Equation from Scire and Robe 2000
        where(hfss <= 0)
            slp_wnd_spd = (hfss*9.81*r_dist*0.5*cos(slope))/(rho(:,kms,:)*temperature(:,kms,:)*wrf_cp*C_d*2+epsilon)
            slp_wnd_spd = sign(abs(slp_wnd_spd)**(thrd),slp_wnd_spd)
            slp_wnd_spd = slp_wnd_spd*(1 - exp(-r_dist/L_e))**(thrd)
        else where(hfss > 0)
            slp_wnd_spd = (hfss*9.81*v_dist*0.5*cos(slope))/(rho(:,kms,:)*temperature(:,kms,:)*wrf_cp*(0.5+C_d)+epsilon)
            slp_wnd_spd = sign(abs(slp_wnd_spd)**(thrd),slp_wnd_spd)
        end where
        
        
        !Dummy bounding of corrections, for safety
        slp_wnd_spd = min(max(slp_wnd_spd,-10.0),10.0)
        
        ! Determine vertical structure/application of slope flow correction
        
        
        ! LESs of upslope flows showed a peak in wind speed at around 30-50m for all slope angles and weakly stable conditions
        upslope_peak_h = 40.0
        ! Increases to a height of 2-10xm depending on slope angle. Increased slope angle reduces the vertical extent of slope flows
        upslope_max_h = 2*upslope_peak_h + upslope_peak_h*6*(1-min(1.0,slope/(piconst/4))**2)

        !Vertical scaling for downslope flows
        ! Downslope flows tend to peak in wind speed at around 5% of ridge drop, with the flow extending several (~3) times above this height due to entrainment
        downslope_peak_h = 0.05*r_drop
        ! Atmospheric stability and slope angle influence the depth of this flow
        downslope_max_h  = 0.15*r_drop

        if (therm_k_max >= kms) then
            do k=kms,therm_k_max
                up_z_scale = 0
                dwn_z_scale = 0
                ! From the surface to the height of the maximum disturbance, winds increase roughly exponentially. 
                ! From the heigt of maximum to the zero height, winds decrease roughly linearly
                where (level_height(:,k,:) <= upslope_peak_h)
                    up_z_scale = max(1.0,log(level_height(:,k,:)/z_0))/(log(upslope_peak_h/z_0)+epsilon)
                end where

                where (level_height(:,k,:) > upslope_peak_h .and. level_height(:,k,:) < upslope_max_h )
                    up_z_scale = 1.0 - (level_height(:,k,:)-upslope_peak_h)/(upslope_max_h-upslope_peak_h+epsilon)
                end where

                where (level_height(:,k,:) <= downslope_peak_h)
                    dwn_z_scale = max(1.0,log(level_height(:,k,:)/z_0))/(log(downslope_peak_h/z_0)+epsilon)
                end where

                where (level_height(:,k,:) > downslope_peak_h .and. level_height(:,k,:) < downslope_max_h )
                    dwn_z_scale = 1.0 - (level_height(:,k,:)-downslope_peak_h)/(downslope_max_h-downslope_peak_h+epsilon)
                end where

                !if (k==2) then
                !    call io_write("dwn_z_scale.nc", "dwn_z_scale", dwn_z_scale(:,:) )
                !    call io_write("up_z_scale.nc", "up_z_scale", up_z_scale(:,:) )
                !endif
            
                where((hfss < 0)) ! .and. (level_height(:,k,:) < r_drop*0.05 ))
                    therm_U_corr(:,k,:) = slp_wnd_spd*x_norm*dwn_z_scale
                    therm_V_corr(:,k,:) = slp_wnd_spd*y_norm*dwn_z_scale
                end where
                where((hfss > 0)) ! .and. (level_height(:,k,:) < r_drop*0.05 ))
                    therm_U_corr(:,k,:) = slp_wnd_spd*x_norm*up_z_scale
                    therm_V_corr(:,k,:) = slp_wnd_spd*y_norm*up_z_scale
                end where
            enddo
        endif
        !Finally, apply TPI and Sx corrections, staggering corrections to U/V grids
        u(ims+1:ime,kms:therm_k_max,:) = u(ims+1:ime,kms:therm_k_max,:) + ( (therm_U_corr(ims+1:ime,:,:) + therm_U_corr(ims:ime-1,:,:))/2 )
        u(ims,kms:therm_k_max,:) = u(ims,kms:therm_k_max,:) + therm_U_corr(ims,:,:)
        u(ime+1,kms:therm_k_max,:) = u(ime+1,kms:therm_k_max,:) + therm_U_corr(ime,:,:)

        v(:,kms:therm_k_max,jms+1:jme) = v(:,kms:therm_k_max,jms+1:jme) + ( (therm_V_corr(:,:,jms+1:jme) + therm_V_corr(:,:,jms:jme-1))/2 )
        v(:,kms:therm_k_max,jms) = v(:,kms:therm_k_max,jms) + therm_V_corr(:,:,jms)
        v(:,kms:therm_k_max,jme+1) = v(:,kms:therm_k_max,jme+1) + therm_V_corr(:,:,jme)

    end subroutine apply_thermal_winds
    
    
end module wind_thermal
