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
    use mod_wrf_constants, only : wrf_cp, epsilon
    implicit none
    private
    public :: apply_thermal_winds, init_thermal_winds
    real, parameter :: thrd=1.0/3.0
    real, parameter :: L_e = 1000.0 ! Scale length for downslope flows (m)
    real, parameter :: Max_flow_height = 200 ! Maximum height for thermal flows above the surface (m)
    integer, save   :: therm_k_max = 0          
    
    integer, save :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer, save :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer, save :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions

contains
    
    subroutine init_thermal_winds(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        
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
                z_mean = SUM(domain%global_z_interface(:,k,:)-domain%global_z_interface(:,kms,:))/SIZE(domain%global_z_interface(:,k,:))
                if (z_mean > Max_flow_height .and. therm_k_max==0) therm_k_max = max(2,k-1)
            enddo
        endif

    end subroutine
    
    subroutine apply_thermal_winds(hfss, rho, temperature, r_dist, v_dist, r_drop, u, v, dzdx, dzdy, z, z_0) !v_dist, r_drop, v_drop,
        implicit none
        real, intent(in), dimension(ims:ime,kms:kme,jms:jme)      :: rho, temperature, dzdx, dzdy, z
        real, intent(in), dimension(ims:ime,jms:jme)              :: hfss, r_dist, v_dist, r_drop, z_0
        real, intent(inout), dimension(ims:ime+1,kms:kme,jms:jme) :: u
        real, intent(inout), dimension(ims:ime,kms:kme,jms:jme+1) :: v

        real, dimension(ims:ime,jms:jme)   :: x_norm, y_norm, slope, slp_wnd_spd, spd, C_d
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
        C_d = 0.41/(log(z(:,kms,:)/z_0)*spd+epsilon)
        
        !Calculate horizontal normal vector components of terrain
        slope = sqrt( dzdx(ims:ime,kms,jms:jme)**2+dzdy(ims:ime,kms,jms:jme)**2)
        where(.not.(slope == 0.)) x_norm = dzdx(ims:ime,kms,jms:jme)/slope
        where(.not.(slope == 0.)) y_norm = dzdy(ims:ime,kms,jms:jme)/slope
        slope = atan(slope)
        
        !Don't allow for crazy large slope flows/valley flows just because we end up far from a ridge
        !r_dist = max(r_dist,10000)
        
        ! Equation from Scire and Robe 2000
        where(hfss <= 0)
            slp_wnd_spd = (hfss*9.81*r_dist*sin(slope))/(rho(:,kms,:)*temperature(:,kms,:)*wrf_cp*C_d*2+epsilon)
            slp_wnd_spd = sign(abs(slp_wnd_spd)**(thrd),slp_wnd_spd)
            slp_wnd_spd = slp_wnd_spd*(1 - exp(-r_dist/L_e))**(thrd)
        else where(hfss > 0)
            slp_wnd_spd = (hfss*9.81*v_dist*sin(slope))/(rho(:,kms,:)*temperature(:,kms,:)*wrf_cp*(0.5+C_d)+epsilon)
            slp_wnd_spd = sign(abs(slp_wnd_spd)**(thrd),slp_wnd_spd)
        end where
        
        !Dummy bounding of corrections, for safety
        slp_wnd_spd = min(max(slp_wnd_spd,-10.0),10.0)
        
        ! If we are at the first layer, apply correction
        therm_U_corr(:,kms,:) = slp_wnd_spd*x_norm
        therm_V_corr(:,kms,:) = slp_wnd_spd*y_norm

        ! If we have katabatic flows (hfss < 0) and we are below the max flow height for this point, apply correction
        if (therm_k_max > kms) then
            do k=kms+1,therm_k_max
                where((hfss < 0) .and. (z(:,k,:) < r_drop*0.05 ))
                    therm_U_corr(:,k,:) = slp_wnd_spd*x_norm
                    therm_V_corr(:,k,:) = slp_wnd_spd*y_norm
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
