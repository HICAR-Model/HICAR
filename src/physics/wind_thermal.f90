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
    use mod_wrf_constants, only : cp, epsilon, piconst
    implicit none
    private
    public :: apply_thermal_winds, init_thermal_winds
    !real, parameter :: thrd=1.0/3.0
    !real, parameter :: L_e = 1000.0 ! Scale length for downslope flows (m)
    real, parameter :: Pr = 2.            !assumed constant Prandl number of 2, from Grisogono et al., 2014
    real, parameter :: Pr_sqrt = sqrt(Pr) !assumed constant Prandl number of 2, from Grisogono et al., 2014
    real, parameter :: Max_flow_height = 400 ! Maximum height for thermal flows above the surface (m)
    
    real, parameter :: K_const = 0.30 !m2 s−1, from (https://doi.org/10.1175/1520-0469(2001)058<3349:KFASFG>2.0.CO;2)
    real, parameter :: h_const = 30.0 !m, from (https://doi.org/10.1175/1520-0469(2001)058<3349:KFASFG>2.0.CO;2)
    
    real, allocatable :: level_height(:,:,:)
    integer, allocatable :: k_200m(:,:), k_1200m(:,:) !indices at which we are roughly 200m or 1200m above the surface, respectively 
    integer, save   :: therm_k_max = 0          
        
    integer, save :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer, save :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer, save :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions
    
    integer, save :: i_s, i_e, j_s, j_e

contains
    
    subroutine init_thermal_winds(domain, options)
        implicit none
        type(domain_t),  intent(in) :: domain
        type(options_t), intent(in) :: options

        real :: z_mean
        integer :: i, k, j

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
        
        i_s = ims
        i_e = ime
        j_s = jms
        j_e = jme

        if (ims==ids) i_s = its
        if (ime==ide) i_e = ite
        if (jms==jds) j_s = jts
        if (jme==jde) j_e = jte

        allocate(k_200m(ims:ime,jms:jme))
        allocate(k_1200m(ims:ime,jms:jme))
        k_200m = 0
        k_1200m = 0

        do j = jms,jme
            do i = ims,ime
                do k = kms,kme
                    if ( (k_200m(i,j)==0) .and. ((domain%z%data_3d(i,k,j)-domain%z_interface%data_3d(i,kms,j)) >= 200.0) ) k_200m(i,j) = k
                    if ( (k_1200m(i,j)==0) .and. ((domain%z%data_3d(i,k,j)-domain%z_interface%data_3d(i,kms,j)) >= 1200.0) ) then
                        k_1200m(i,j) = k
                        exit
                    endif
                enddo
            enddo
        enddo

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
    
    subroutine apply_thermal_winds(tskin, exner, pot, z, dz_mass, r_dist, v_dist, r_drop, u, v, dzdx, dzdy, K_h) !v_dist, r_drop, v_drop,
        implicit none
        real, intent(in), dimension(ims:ime,kms:kme,jms:jme)      :: exner, pot, z, dz_mass, dzdx, dzdy
        real, intent(in), dimension(ims:ime,jms:jme)              :: tskin, r_dist, v_dist, r_drop
        real, intent(inout), dimension(ims:ime+1,kms:kme,jms:jme) :: u
        real, intent(inout), dimension(ims:ime,kms:kme,jms:jme+1) :: v
        real, optional, intent(in), dimension(ims:ime,kms:kme,jms:jme)      :: K_h

        real, dimension(ims:ime,jms:jme)   :: x_norm, y_norm, slope!, spd, C_d, upslope_peak_h, upslope_max_h, downslope_peak_h, downslope_max_h, up_z_scale, dwn_z_scale
        real, dimension(ims:ime,kms:therm_k_max,jms:jme) :: therm_U_corr, therm_V_corr

        integer ::  i, j, k
        real    :: C_1, C_2, C_3, C_4, C_5
        real    :: K_int, slp_wnd_spd, gamma, theta_0, C, N, h_p, sig_0, mu, small, Iz, u_A, u_0, u_1, K_sca

        do j = jms,jme
            do i = ims,ime
                x_norm(i,j) = 0
                y_norm(i,j) = 0
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
        !spd = sqrt((u(ims+1:ime+1,kms,:)*0.5+u(ims:ime,kms,:)*0.5)**2 + (v(:,kms,jms+1:jme+1)*0.5+v(:,kms,jms:jme)*0.5)**2)
        !C_d = 0.41/(log(level_height(:,kms,:)/z_0)*spd+epsilon)
        
        !Calculate horizontal normal vector components of terrain
        slope = sqrt( dzdx(ims:ime,kms,jms:jme)**2+dzdy(ims:ime,kms,jms:jme)**2)
        where(.not.(slope == 0.)) x_norm = dzdx(ims:ime,kms,jms:jme)/slope
        where(.not.(slope == 0.)) y_norm = dzdy(ims:ime,kms,jms:jme)/slope
        slope = atan(slope)
        
        
        C_1 = -1/3.
        C_2 = 2/15.
        C_3 = 1/30.
        C_4 = -1/30.
        C_5 = -1/10.
        
        !If we are on the domain boundary, then tskin is never updated, since LSM only gives values for the process tile grid cells
        do j = j_s,j_e
            do i = i_s,i_e
                theta_0 = pot(i,k_200m(i,j),j)  !background potential temperature
                gamma = (pot(i,k_1200m(i,j),j)-theta_0)/(z(i,k_1200m(i,j),j)-z(i,k_200m(i,j),j)) !potential temperature lapse rate (K/m)
                !Method only valid for statically stable conditions (plus, mathematically, gamma is essentially always positive where used (sqrts))
                gamma = max(0.0001,gamma)
                C = tskin(i,j)/exner(i,kms,j)-theta_0 !surface potential temperature anomoly

                N = sqrt(gamma*9.81/theta_0)
                sig_0 = N*sin(slope(i,j))/Pr_sqrt
                mu = sqrt(9.81/(theta_0*gamma*Pr+epsilon))
                
                K_int = 0.0
                do k = kms,therm_k_max
                    !if (present(K_h)) then
                    !    K_sca = K_h(i,k,j)+0.0001
                    !else
                    K_sca = K_const*( level_height(i,k,j)/h_const)*exp(-0.5*( level_height(i,k,j)/h_const)**2 )+0.0001
                    !endif
                    
                    K_int = K_int + (1./sqrt(K_sca)) * (dz_mass(i,k,j)*cos(slope(i,j))) !Adjust z for slanted slope...
                    Iz =  sqrt(sig_0*0.5)*K_int

                    u_A = sqrt(sig_0*0.5)*(C**2)*(mu/(gamma+epsilon))/sqrt(K_sca)
                    
                    u_1 = u_A*exp(-Iz)  *(C_1*sin(Iz)  +C_2*cos(Iz)) + &
                          u_A*exp(-2*Iz)*(C_3*sin(2*Iz)+C_4*cos(2*Iz)+C_5)
                        
                    u_0 = -C*mu*exp(-Iz)*sin(Iz)
                    
                    !Determine small parameter
                    if (k==kms) then
                        h_p = sqrt(2.)/sqrt((N*sin(slope(i,j)))/(K_sca*Pr_sqrt))

                        !!For katabatic flows, small should be bounded up to 0.01, while for anabatic flows it should be bounded up to 0.06)
                        if (C < 0 ) then
                            small = min(0.01,(2*gamma*h_p/abs(C+epsilon)))
                        else
                            small = min(0.06,(2*gamma*h_p/abs(C+epsilon)))
                        endif
                    endif
                    
                    slp_wnd_spd = u_0 + small*u_1
                    ! slp_wnd_spd must be multiplied by -1 here to get correct direction, 
                    ! due to how x/y_norm are calculated
                    slp_wnd_spd = sign(min(abs(slp_wnd_spd),6.),-slp_wnd_spd) !dummy bound
                    
                    therm_U_corr(i,k,j) = slp_wnd_spd*x_norm(i,j)
                    therm_V_corr(i,k,j) = slp_wnd_spd*y_norm(i,j)
                enddo
            enddo
        enddo
        !
        
        
        !!Don't allow for crazy large slope flows/valley flows just because we end up far from a ridge
        !!r_dist = max(r_dist,10000)
        !
        !! Solve for maximum wind perturbation from slope flows
        !! Equation from Scire and Robe 2000
        !where(hfss <= 0)
        !    slp_wnd_spd = (hfss*9.81*r_dist*cos(slope))/(rho(:,kms,:)*temperature(:,kms,:)*cp*C_d*2+epsilon)
        !    slp_wnd_spd = sign(abs(slp_wnd_spd)**(thrd),slp_wnd_spd)
        !    slp_wnd_spd = slp_wnd_spd*(1 - exp(-r_dist/L_e))**(thrd)
        !else where(hfss > 0)
        !    slp_wnd_spd = (hfss*9.81*v_dist*cos(slope))/(rho(:,kms,:)*temperature(:,kms,:)*cp*(0.5+C_d)+epsilon)
        !    slp_wnd_spd = sign(abs(slp_wnd_spd)**(thrd),slp_wnd_spd)
        !end where
        !
        !
        !!Dummy bounding of corrections, for safety
        !slp_wnd_spd = min(max(slp_wnd_spd,-10.0),10.0)
        !
        !! Determine vertical structure/application of slope flow correction
        !
        !
        !! LESs of upslope flows showed a peak in wind speed at around 30-50m for all slope angles and weakly stable conditions
        !upslope_peak_h = 40.0
        !! Increases to a height of 2-10xm depending on slope angle. Increased slope angle reduces the vertical extent of slope flows
        !upslope_max_h = 2*upslope_peak_h + upslope_peak_h*6*(1-min(1.0,slope/(piconst/4))**2)
!
        !!Vertical scaling for downslope flows
        !! Downslope flows tend to peak in wind speed at around 5% of ridge drop, with the flow extending several (~3) times above this height due to entrainment
        !downslope_peak_h = 0.05*r_drop
        !! Atmospheric stability and slope angle influence the depth of this flow
        !downslope_max_h  = 0.15*r_drop
!
        !if (therm_k_max >= kms) then
        !    do k=kms,therm_k_max
        !        up_z_scale = 0
        !        dwn_z_scale = 0
        !        ! From the surface to the height of the maximum disturbance, winds increase roughly exponentially. 
        !        ! From the heigt of maximum to the zero height, winds decrease roughly linearly
        !        where (level_height(:,k,:) <= upslope_peak_h)
        !            up_z_scale = max(1.0,log(level_height(:,k,:)/z_0))/(log(upslope_peak_h/z_0)+epsilon)
        !        end where
!
        !        where (level_height(:,k,:) > upslope_peak_h .and. level_height(:,k,:) < upslope_max_h )
        !            up_z_scale = 1.0 - (level_height(:,k,:)-upslope_peak_h)/(upslope_max_h-upslope_peak_h+epsilon)
        !        end where
!
        !        where (level_height(:,k,:) <= downslope_peak_h)
        !            dwn_z_scale = max(1.0,log(level_height(:,k,:)/z_0))/(log(downslope_peak_h/z_0)+epsilon)
        !        end where
!
        !        where (level_height(:,k,:) > downslope_peak_h .and. level_height(:,k,:) < downslope_max_h )
        !            dwn_z_scale = 1.0 - (level_height(:,k,:)-downslope_peak_h)/(downslope_max_h-downslope_peak_h+epsilon)
        !        end where
!
        !        !if (k==2) then
        !        !    call io_write("dwn_z_scale.nc", "dwn_z_scale", dwn_z_scale(:,:) )
        !        !    call io_write("up_z_scale.nc", "up_z_scale", up_z_scale(:,:) )
        !        !endif
        !    
        !        where((hfss < 0)) ! .and. (level_height(:,k,:) < r_drop*0.05 ))
        !            therm_U_corr(:,k,:) = slp_wnd_spd*x_norm*dwn_z_scale
        !            therm_V_corr(:,k,:) = slp_wnd_spd*y_norm*dwn_z_scale
        !        end where
        !        where((hfss > 0)) ! .and. (level_height(:,k,:) < r_drop*0.05 ))
        !            therm_U_corr(:,k,:) = slp_wnd_spd*x_norm*up_z_scale
        !            therm_V_corr(:,k,:) = slp_wnd_spd*y_norm*up_z_scale
        !        end where
        !    enddo
        !endif
        !Finally, apply TPI and Sx corrections, staggering corrections to U/V grids
        u(i_s+1:i_e,kms:therm_k_max,:) = u(i_s+1:i_e,kms:therm_k_max,:) + &
            ( (therm_U_corr(i_s+1:i_e,kms:therm_k_max,:) + therm_U_corr(i_s:i_e-1,kms:therm_k_max,:))/2 )
        u(i_s,kms:therm_k_max,:) = u(i_s,kms:therm_k_max,:) + therm_U_corr(i_s,kms:therm_k_max,:)
        u(i_e+1,kms:therm_k_max,:) = u(i_e+1,kms:therm_k_max,:) + therm_U_corr(i_e,kms:therm_k_max,:)

        v(:,kms:therm_k_max,j_s+1:j_e) = v(:,kms:therm_k_max,j_s+1:j_e) + &
            ( (therm_V_corr(:,kms:therm_k_max,j_s+1:j_e) + therm_V_corr(:,kms:therm_k_max,j_s:j_e-1))/2 )
        v(:,kms:therm_k_max,j_s) = v(:,kms:therm_k_max,j_s) + therm_V_corr(:,kms:therm_k_max,j_s)
        v(:,kms:therm_k_max,j_e+1) = v(:,kms:therm_k_max,j_e+1) + therm_V_corr(:,kms:therm_k_max,j_e)

    end subroutine apply_thermal_winds
    
    
end module wind_thermal
