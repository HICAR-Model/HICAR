!> ----------------------------------------------------------------------------
!!  Standard advection scheme with variable order
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!! ----------------------------------------------------------------------------
module adv_std
    use data_structures
    use icar_constants
    use options_interface, only: options_t
    use domain_interface,  only: domain_t
    use adv_fluxcorr,      only: WRF_flux_corr
    implicit none
    private
    real,dimension(:,:,:),allocatable :: U_m, V_m, W_m, rho
    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte, i_s, i_e, j_s, j_e, horder, vorder
    real    :: dx

    ! For use advecting a (convective?) wind field
    ! real,dimension(:,:,:),allocatable :: U_4cu_u, V_4cu_u, W_4cu_u
    ! real,dimension(:,:,:),allocatable :: U_4cu_v, V_4cu_v, W_4cu_v

    public :: adv_std_init, adv_std_var_request, adv_std_advect3d, adv_fluxcorr_advect3d, adv_std_compute_wind

contains

    subroutine adv_std_init(domain,options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options


        ims = domain%ims
        ime = domain%ime
        jms = domain%jms
        jme = domain%jme
        kms = domain%kms
        kme = domain%kme
        its = domain%its
        ite = domain%ite
        jts = domain%jts
        jte = domain%jte
        
        !set order of advection
        horder = options%adv_options%h_order
        vorder = options%adv_options%v_order
        
        !Define bounds of advection computation. If using monotonic flux-limiter, it is necesarry to increase
        !advection bounds by 1. The necesarry extension of the halo is handeled in domain_object
        i_s = its
        i_e = ite
        j_s = jts
        j_e = jte
        
        if (options%adv_options%flux_corr==kFLUXCOR_MONO) then
            i_s = its - 1
            i_e = ite + 1
            j_s = jts - 1
            j_e = jte + 1
        endif
        
        ! if module level arrays are already allocated for some reason, deallocate them first
        if (allocated(U_m)) deallocate(U_m)
        if (allocated(V_m)) deallocate(V_m)
        if (allocated(W_m)) deallocate(W_m)
        !if (allocated(lastqv_m)) deallocate(lastqv_m)

        ! allocate the module level arrays
        allocate(U_m     (i_s:i_e+1,kms:kme,j_s:j_e  ))
        allocate(V_m     (i_s:i_e,  kms:kme,j_s:j_e+1))
        allocate(W_m     (i_s:i_e,  kms:kme,j_s:j_e  ))
        allocate(rho     (ims:ime,  kms:kme,jms:jme  ))
        !allocate(lastqv_m(ims:ime,  kms:kme,jms:jme  ))

        !     if (.not.allocated(U_4cu_u)) then
        !         allocate(U_4cu_u(nx,  nz, ny))
        !         U_4cu_u = 0
        !         allocate(V_4cu_u(nx+1,nz, ny-1))
        !         V_4cu_u = 0
        !         allocate(W_4cu_u(nx+1,nz, ny))
        !         W_4cu_u = 0
        !
        !         allocate(U_4cu_v(nx-1,nz, ny+1))
        !         U_4cu_v = 0
        !         allocate(V_4cu_v(nx,  nz, ny))
        !         V_4cu_v = 0
        !         allocate(W_4cu_v(nx,  nz, ny+1))
        !         W_4cu_v = 0
        !     endif
    end subroutine

    subroutine adv_std_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for adv4 advection
        call options%alloc_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface, kVARS%water_vapor])

        call options%advect_vars([kVARS%water_vapor, kVARS%potential_temperature])

        ! List the variables that are required for restarts with adv4 advection
        call options%restart_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface, kVARS%water_vapor])

    end subroutine

!     Note this routine has been manually inlined because the compiler didn't seem to optimize it well and made the array copies
!     the routine is left in for documentation.
!     subroutine flux2(l,r,U,nx,nz,ny,f)
!     !     Calculate the donor cell flux function
!     !     l = left gridcell scalar
!     !     r = right gridcell scalar
!     !     U = Courant number (u*dt/dx)
!     !
!     !     If U is positive, return l*U if U is negative return r*U
!     !     By using the mathematical form instead of the logical form,
!     !     we can run on the entire grid simultaneously, and avoid branches
!
!     !   arguments
!         implicit none
!         real, dimension(1:nx,1:nz,1:ny), intent(in) :: l,r,U
!         real, dimension(1:nx,1:nz,1:ny), intent(inout) :: f
!         integer,intent(in) :: ny,nz,nx
!         !   internal parameter
!         integer ::  err,i!,j,Ny,Nz,Nx
!         !   main code
!         f= ((U+ABS(U)) * l + (U-ABS(U)) * r)/2
!
!     end subroutine flux2

    subroutine flux3(q,flux_x,flux_z,flux_y,t_factor)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),   intent(in)       :: q
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e),   intent(inout)    :: flux_x
        real, dimension(i_s:i_e,  kms:kme,j_s:j_e+1), intent(inout)    :: flux_y
        real, dimension(i_s:i_e,  kms:kme+1,j_s:j_e), intent(inout)    :: flux_z
        real, intent(in) :: t_factor
        integer :: i, j, k
        real :: tmp, coef, u, q0, q1, q2, qn1, qn2, qn3
                      
                      
        if (horder==1) then
            do j = j_s,j_e
                do k = kms,kme
                    do i = i_s,i_e+1
                        flux_x(i,k,j)= ((U_m(i,k,j) + ABS(U_m(i,k,j))) * q(i-1,k,j) + (U_m(i,k,j) - ABS(U_m(i,k,j))) * q(i,k,j))  * 0.5 * t_factor
                    enddo
                enddo
            enddo
            do j = j_s,j_e+1
                do k = kms,kme
                    do i = i_s,i_e
                        flux_y(i,k,j)= ((V_m(i,k,j) + ABS(V_m(i,k,j))) * q(i,k,j-1) + (V_m(i,k,j) - ABS(V_m(i,k,j))) * q(i,k,j))  * 0.5 * t_factor
                    enddo
                enddo
            enddo
        else if (horder==3) then
            coef = (1./12)*t_factor
            !DIR$ UNROLL 5
            do j = j_s,j_e
                do k = kms,kme
                    do i = i_s,i_e+1
                        u = U_m(i,k,j)
                        q0  = q(i,k,j);   q1  = q(i+1,k,j)
                        qn1 = q(i-1,k,j); qn2 = q(i-2,k,j)
                        !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                        tmp = 7*(q0+qn1) - (q1+qn2)
                        tmp = u*tmp
                        !Application of 3rd order diffusive terms
                        tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
                        flux_x(i,k,j) = tmp*coef                
                    enddo
                enddo
            enddo
            !DIR$ UNROLL 5
            do j = j_s,j_e+1
                do k = kms,kme
                    do i = i_s,i_e
                        u = V_m(i,k,j)
                        q0  = q(i,k,j);   q1  = q(i,k,j+1)
                        qn1 = q(i,k,j-1); qn2 = q(i,k,j-2)
                        !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                        tmp = 7*(q0+qn1) - (q1+qn2)
                        tmp = u*tmp
                        !Application of 3rd order diffusive terms
                        tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
                        flux_y(i,k,j) = tmp*coef                
                    enddo
                enddo
            enddo
        else if (horder==5) then
            coef = (1./60)*t_factor
            !DIR$ UNROLL 5
            do j = j_s,j_e
                do k = kms,kme
                    do i = i_s,i_e+1
                        u = U_m(i,k,j)
                        q0  = q(i,k,j);   q1  = q(i+1,k,j); q2  = q(i+2,k,j)
                        qn1 = q(i-1,k,j); qn2 = q(i-2,k,j); qn3 = q(i-3,k,j)

                        !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                        tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                        tmp = u*tmp
                        !Application of 5th order diffusive terms
                        tmp = tmp - abs(u) * (10*(q0-qn1) - 5*(q1-qn2) + (q2-qn3))
                        flux_x(i,k,j) = tmp*coef
                    enddo
                enddo
            enddo
            !DIR$ UNROLL 5
            do j = j_s,j_e+1
                do k = kms,kme
                    do i = i_s,i_e
                        u = V_m(i,k,j)
                        q0  = q(i,k,j);   q1  = q(i,k,j+1); q2  = q(i,k,j+2)
                        qn1 = q(i,k,j-1); qn2 = q(i,k,j-2); qn3 = q(i,k,j-3)
                        !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                        tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                        tmp = u*tmp
                        !Application of 5th order diffusive terms
                        tmp = tmp - abs(u) * (10*(q0-qn1) -  5*(q1-qn2) + (q2-qn3))
                        flux_y(i,k,j) = tmp*coef
                    enddo
                enddo
            enddo
        endif
        
        if (vorder==1) then
            do j = j_s,j_e
                do k = kms+1,kme
                    do i = i_s,i_e
                        flux_z(i,k,j) = ((W_m(i,k-1,j) + ABS(W_m(i,k-1,j))) * q(i,k-1,j) + &
                                     (W_m(i,k-1,j) - ABS(W_m(i,k-1,j))) * q(i,k,j))  * 0.5  * t_factor
                    enddo
                enddo
            enddo
        else if (vorder==3) then
            coef = (1./12)*t_factor
            do j = j_s,j_e
                do k = kms+2,kme-1
                    do i = i_s,i_e
                        u = W_m(i,k-1,j)
                        q0  = q(i,k,j);   q1  = q(i,k+1,j)
                        qn1 = q(i,k-1,j); qn2 = q(i,k-2,j)
                        !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                        tmp = 7*(q0+qn1) - (q1+qn2)
                        tmp = u*tmp
                        !Application of 3rd order diffusive terms
                        tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                        flux_z(i,k,j) = tmp*coef               
                    enddo
                enddo
            enddo
            do j = j_s,j_e
                do i = i_s,i_e
                    u = W_m(i,kms,j)
                    !Do simple upwind for the cells who's stencil does not allow higher-order
                    flux_z(i,kms+1,j) = ((u + ABS(u)) * q(i,kms,j) + (u - ABS(u)) * q(i,kms+1,j))  * 0.5  * t_factor
                    u = W_m(i,kme-1,j)
                    flux_z(i,kme,j) = ((u + ABS(u)) * q(i,kme-1,j) + (u - ABS(u)) * q(i,kme,j))  * 0.5  * t_factor
                enddo
            enddo
        else if (vorder==5) then
            coef = (1./60)*t_factor
            do j = j_s,j_e
                do k = kms+3,kme-2
                    do i = i_s,i_e
                        u = W_m(i,k-1,j)
                        q0  = q(i,k,j);   q1  = q(i,k+1,j);  q2 = q(i,k+2,j)
                        qn1 = q(i,k-1,j); qn2 = q(i,k-2,j); qn3 = q(i,k-3,j)
                        !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                        tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                        tmp = u*tmp
                        !Application of 5th order diffusive terms
                        tmp = tmp - abs(u) * (10*(q0-qn1) -  5*(q1-qn2) + (q2-qn3))
                        flux_z(i,k,j) = tmp*coef
                    enddo
                enddo
            enddo
            coef = (1./12)*t_factor
            do j = j_s,j_e
                do i = i_s,i_e
                    u = W_m(i,kms+1,j)
                    q0  = q(i,kms+2,j);   q1  = q(i,kms+3,j)
                    qn1 = q(i,kms-1,j);   qn2 = q(i,kms,j)
                    
                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                    flux_z(i,kms+2,j) = tmp*coef               
                    
                    u = W_m(i,kme-2,j)
                    q0  = q(i,kme-1,j);   q1  = q(i,kme,j)
                    qn1 = q(i,kme-2,j);   qn2 = q(i,kme-3,j)

                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                    flux_z(i,kme-1,j) = tmp*coef               

                    u = W_m(i,kms,j)
                    !Do simple upwind for the cells who's stencil does not allow higher-order
                    flux_z(i,kms+1,j) = ((u + ABS(u)) * q(i,kms,j) + &
                                         (u - ABS(u)) * q(i,kms+1,j))  * 0.5 * t_factor
                    u = W_m(i,kme-1,j)
                    flux_z(i,kme,j) = ((u + ABS(u)) * q(i,kme-1,j) + &
                                       (u - ABS(u)) * q(i,kme,j))  * 0.5 * t_factor
                enddo
            enddo
        endif
                                                          
        !Handle top and bottom boundaries for z here
        do j = j_s,j_e
            do i = i_s,i_e
                flux_z(i,kms,j) = 0
                flux_z(i,kme+1,j) = q(i,kme,j) * W_m(i,kme,j) * t_factor
            enddo
        enddo
    end subroutine flux3

    subroutine flux3_w_up(q,flux_x,flux_z,flux_y,flux_x_up,flux_z_up,flux_y_up)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),   intent(in)       :: q
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e),   intent(inout)    :: flux_x, flux_x_up
        real, dimension(i_s:i_e,  kms:kme,j_s:j_e+1), intent(inout)    :: flux_y, flux_y_up
        real, dimension(i_s:i_e,  kms:kme+1,j_s:j_e), intent(inout)    :: flux_z, flux_z_up
        integer :: i, j, k
        real :: tmp, coef, u, q0, q1, q2, qn1, qn2, qn3
                      
        if (horder==3) then
            coef = (1./12)
            !DIR$ UNROLL 5
            do j = j_s,j_e
                do k = kms,kme
                    do i = i_s,i_e+1
                        u = U_m(i,k,j)
                        q0  = q(i,k,j);   q1  = q(i+1,k,j)
                        qn1 = q(i-1,k,j); qn2 = q(i-2,k,j)
                        !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                        tmp = 7*(q0+qn1) - (q1+qn2)
                        tmp = u*tmp
                        !Application of 3rd order diffusive terms
                        tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
                        !Calculation of Upwind fluxes
                        flux_x_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                        !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                        flux_x(i,k,j) = tmp*coef
                    enddo
                enddo
            enddo
            !DIR$ UNROLL 5
            do j = j_s,j_e+1
                do k = kms,kme
                    do i = i_s,i_e
                        u = V_m(i,k,j)
                        q0  = q(i,k,j);   q1  = q(i,k,j+1)
                        qn1 = q(i,k,j-1); qn2 = q(i,k,j-2)
                        !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                        tmp = 7*(q0+qn1) - (q1+qn2)
                        tmp = u*tmp
                        !Application of 3rd order diffusive terms
                        tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
                        !Calculation of Upwind fluxes
                        flux_y_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                        !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                        flux_y(i,k,j) = tmp*coef         
                    enddo
                enddo
            enddo
        else if (horder==5) then
            coef = (1./60)
            !DIR$ UNROLL 5
            do j = j_s,j_e
                do k = kms,kme
                    do i = i_s,i_e+1
                        u = U_m(i,k,j)
                        q0  = q(i,k,j);   q1  = q(i+1,k,j); q2  = q(i+2,k,j)
                        qn1 = q(i-1,k,j); qn2 = q(i-2,k,j); qn3 = q(i-3,k,j)

                        !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                        tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                        tmp = u*tmp
                        !Application of 5th order diffusive terms
                        tmp = tmp - abs(u) * (10*(q0-qn1) - 5*(q1-qn2) + (q2-qn3))
                        !Calculation of Upwind fluxes
                        flux_x_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                        !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                        flux_x(i,k,j) = tmp*coef
                    enddo
                enddo
            enddo
            !DIR$ UNROLL 5
            do j = j_s,j_e+1
                do k = kms,kme
                    do i = i_s,i_e
                        u = V_m(i,k,j)
                        q0  = q(i,k,j);   q1  = q(i,k,j+1); q2  = q(i,k,j+2)
                        qn1 = q(i,k,j-1); qn2 = q(i,k,j-2); qn3 = q(i,k,j-3)
                        !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                        tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                        tmp = u*tmp
                        !Application of 5th order diffusive terms
                        tmp = tmp - abs(u) * (10*(q0-qn1) -  5*(q1-qn2) + (q2-qn3))
                        !Calculation of Upwind fluxes
                        flux_y_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                        !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                        flux_y(i,k,j) = tmp*coef
                    enddo
                enddo
            enddo
        endif
        
        if (vorder==3) then
            coef = (1./12)
            do j = j_s,j_e
                do k = kms+2,kme-1
                    do i = i_s,i_e
                        u = W_m(i,k-1,j)
                        q0  = q(i,k,j);   q1  = q(i,k+1,j)
                        qn1 = q(i,k-1,j); qn2 = q(i,k-2,j)
                        !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                        tmp = 7*(q0+qn1) - (q1+qn2)
                        tmp = u*tmp
                        !Application of 3rd order diffusive terms
                        tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                        !Calculation of Upwind fluxes
                        flux_z_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                        !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                        flux_z(i,k,j) = tmp*coef         
                    enddo
                enddo
            enddo
            do j = j_s,j_e
                do i = i_s,i_e
                    u = W_m(i,kms,j)
                    !Do simple upwind for the cells who's stencil does not allow higher-order
                    flux_z(i,kms+1,j) = ((u + ABS(u)) * q(i,kms,j) + (u - ABS(u)) * q(i,kms+1,j))  * 0.5  
                    flux_z_up(i,kms+1,j) = flux_z(i,kms+1,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
                    
                    u = W_m(i,kme-1,j)
                    flux_z(i,kme,j) = ((u + ABS(u)) * q(i,kme-1,j) + (u - ABS(u)) * q(i,kme,j))  * 0.5  
                    flux_z_up(i,kme,j) = flux_z(i,kme,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
                enddo
            enddo
        else if (vorder==5) then
            coef = (1./60)
            do j = j_s,j_e
                do k = kms+3,kme-2
                    do i = i_s,i_e
                        u = W_m(i,k-1,j)
                        q0  = q(i,k,j);   q1  = q(i,k+1,j);  q2 = q(i,k+2,j)
                        qn1 = q(i,k-1,j); qn2 = q(i,k-2,j); qn3 = q(i,k-3,j)
                        !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                        tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                        tmp = u*tmp
                        !Application of 5th order diffusive terms
                        tmp = tmp - abs(u) * (10*(q0-qn1) -  5*(q1-qn2) + (q2-qn3))
                        !Calculation of Upwind fluxes
                        flux_z_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                        !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                        flux_z(i,k,j) = tmp*coef
                    enddo
                enddo
            enddo
            coef = (1./12)
            do j = j_s,j_e
                do i = i_s,i_e
                    u = W_m(i,kms+1,j)
                    q0  = q(i,kms+2,j);   q1  = q(i,kms+3,j)
                    qn1 = q(i,kms-1,j);   qn2 = q(i,kms,j)
                    
                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                    flux_z(i,kms+2,j) = tmp*coef               
                    flux_z_up(i,kms+2,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step

                    u = W_m(i,kme-2,j)
                    q0  = q(i,kme-1,j);   q1  = q(i,kme,j)
                    qn1 = q(i,kme-2,j);   qn2 = q(i,kme-3,j)

                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                    flux_z(i,kme-1,j) = tmp*coef               
                    flux_z_up(i,kme-1,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step

                    u = W_m(i,kms,j)
                    !Do simple upwind for the cells who's stencil does not allow higher-order
                    flux_z(i,kms+1,j) = ((u + ABS(u)) * q(i,kms,j) + &
                                         (u - ABS(u)) * q(i,kms+1,j))  * 0.5
                    flux_z_up(i,kms+1,j) = flux_z(i,kms+1,j) * 0.5 ! additional "0.5" since we only want half of the upwind step

                    u = W_m(i,kme-1,j)
                    flux_z(i,kme,j) = ((u + ABS(u)) * q(i,kme-1,j) + &
                                       (u - ABS(u)) * q(i,kme,j))  * 0.5
                    flux_z_up(i,kme,j) = flux_z(i,kme,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
                enddo
            enddo
        endif
                                                          
        !Handle top and bottom boundaries for z here
        do j = j_s,j_e
            do i = i_s,i_e
                flux_z(i,kms,j) = 0
                flux_z_up(i,kms,j) = 0
                flux_z(i,kme+1,j) = q(i,kme,j) * W_m(i,kme,j)
                flux_z_up(i,kme+1,j) = flux_z(i,kme+1,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
            enddo
        enddo
    end subroutine flux3_w_up


    subroutine adv_std_advect3d(qfluxes,qold,dz,jaco,t_factor_in)
        ! !DIR$ INLINEALWAYS adv_std_advect3d
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(inout)   :: qfluxes
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: qold
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: dz
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: jaco
        real, optional,                              intent(in)      :: t_factor_in
        ! interal parameters
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e)   :: flux_x
        real, dimension(i_s:i_e,  kms:kme,j_s:j_e+1) :: flux_y
        real, dimension(i_s:i_e,  kms:kme+1,j_s:j_e) :: flux_z
        real    :: t_factor
        integer :: i, k, j
        
        
        !Initialize t_factor, which is used during RK time stepping to scale the time step
        t_factor = 1.0
        if (present(t_factor_in)) t_factor = t_factor_in
        
        call flux3(qfluxes,flux_x,flux_z,flux_y,t_factor)
        
        do j = jts,jte
            do k = kms,kme
                do i = its,ite
                    ! perform advection, from difference terms
                    qfluxes(i,k,j)  = qold(i,k,j)  - &
                                  ((flux_x(i+1,k,j) - flux_x(i,k,j))  + &
                                   (flux_y(i,k,j+1) - flux_y(i,k,j))  + &
                                   (flux_z(i,k+1,j) - flux_z(i,k,j))  / &
                                   dz(i,k,j))                         / &
                                   (jaco(i,k,j)*rho(i,k,j))
                enddo
            enddo
        enddo
        
    end subroutine adv_std_advect3d
    
    
    subroutine adv_fluxcorr_advect3d(qfluxes,qold,dz,jaco)
        ! !DIR$ INLINEALWAYS adv_std_advect3d
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(inout)   :: qfluxes
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: qold
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: dz
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: jaco
        ! interal parameters
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e)   :: flux_x, flux_x_up
        real, dimension(i_s:i_e,  kms:kme,j_s:j_e+1) :: flux_y, flux_y_up
        real, dimension(i_s:i_e,  kms:kme+1,j_s:j_e) :: flux_z, flux_z_up

        integer :: i, k, j
        
                
        call flux3_w_up(qfluxes,flux_x,flux_z,flux_y,flux_x_up,flux_z_up,flux_y_up)

        call WRF_flux_corr(qold,U_m,V_m,W_m,flux_x,flux_z,flux_y,flux_x_up,flux_z_up,flux_y_up,jaco,dz,rho)
        
        do j = jts,jte
            do k = kms,kme
                do i = its,ite
                    ! perform advection, from difference terms
                    qfluxes(i,k,j)  = qold(i,k,j)  - &
                                  ((flux_x(i+1,k,j) - flux_x(i,k,j))  + &
                                   (flux_y(i,k,j+1) - flux_y(i,k,j))  + &
                                   (flux_z(i,k+1,j) - flux_z(i,k,j))  / &
                                   dz(i,k,j))                         / &
                                   (jaco(i,k,j)*rho(i,k,j))
                enddo
            enddo
        enddo
        
    end subroutine adv_fluxcorr_advect3d


    ! subroutine setup_cu_winds(domain, options, dt)
    !     implicit none
    !     type(domain_type),  intent(in) :: domain
    !     type(options_type), intent(in) :: options
    !     real,               intent(in) :: dt
    !
    !     real    :: dx
    !     integer :: nx,nz,ny
    !
    !     dx = domain%dx
    !     nx = size(domain%dz,1)
    !     nz = size(domain%dz,2)
    !     ny = size(domain%dz,3)
    !
    !
    !     U_4cu_u           =  (dt/dx) * (domain%u(1:nx,:,:)      + domain%u(2:nx+1,:,:)) / 2
    !     V_4cu_u(2:nx,:,:) =  (dt/dx) * (domain%v(1:nx-1,:,2:ny) + domain%v(2:nx,:,2:ny)) / 2
    !     W_4cu_u(2:nx,:,:) =  (dt/dx) * (domain%w(1:nx-1,:,:)    + domain%w(2:nx,:,:)) / 2
    !     call rebalance_cu_winds(U_4cu_u, V_4cu_u, W_4cu_u)
    !
    !     U_4cu_v(:,:,2:ny) =  (dt/dx) * (domain%u(2:nx,:,1:ny-1) + domain%u(2:nx,:,2:ny)) / 2
    !     V_4cu_v           =  (dt/dx) * (domain%v(:,:,1:ny)      + domain%v(:,:,2:ny+1)) / 2
    !     W_4cu_v(:,:,2:ny) =  (dt/dx) * (domain%w(:,:,1:ny-1)    + domain%w(:,:,2:ny)) / 2
    !     call rebalance_cu_winds(U_4cu_v, V_4cu_v, W_4cu_v)
    !
    ! end subroutine setup_cu_winds

    ! subroutine rebalance_cu_winds(u,v,w)
    !     implicit none
    !     ! u, v, w 3D east-west, south-north, and up-down winds repsectively
    !     ! note for this code, u is [nx-1,nz,ny] and v is [nx,nz,ny-1]
    !     real, dimension(:,:,:), intent(inout) :: u, v, w
    !
    !     real, allocatable, dimension(:,:) :: divergence, du, dv
    !     integer :: i,nx,ny,nz
    !
    !     nx = size(w,1)
    !     nz = size(w,2)
    !     ny = size(w,3)
    !
    !     allocate(divergence(nx-2,ny-2))
    !     allocate(du(nx-2,ny-2))
    !     allocate(dv(nx-2,ny-2))
    !
    !     do i=1,nz
    !         ! calculate horizontal divergence
    !         dv = v(2:nx-1,i,2:ny-1) - v(2:nx-1,i,1:ny-2)
    !         du = u(2:nx-1,i,2:ny-1) - u(1:nx-2,i,2:ny-1)
    !         divergence = du + dv
    !         ! Then calculate w to balance
    !         if (i==1) then
    !             ! if this is the first model level start from 0 at the ground
    !             w(2:nx-1,i,2:ny-1) = 0 - divergence
    !         else
    !             ! else calculate w as a change from w at the level below
    !             w(2:nx-1,i,2:ny-1) = w(2:nx-1,i-1,2:ny-1) - divergence
    !         endif
    !     enddo
    !
    ! end subroutine rebalance_cu_winds
    !
    ! subroutine advect_cu_winds(domain, options, dt)
    !     implicit none
    !     type(domain_type),  intent(inout) :: domain
    !     type(options_type), intent(in)    :: options
    !     real,               intent(in)    :: dt
    !
    !     integer :: nx,nz,ny
    !
    !     nx = size(domain%dz,1)
    !     nz = size(domain%dz,2)
    !     ny = size(domain%dz,3)
    !
    !     ! first put the background u,v,w winds on a staggered grid with respect to the u grid
    !     ! then advect the u winds
    !     if (options%advect_density) then
    !         print*, "ERROR: Density advection not enabled when using convective winds"
    !         print*, "   Requires update to wind.f90 balance_uvw and advect.f90 (at least)"
    !         stop
    !     endif
    !
    !     call setup_cu_winds(domain, options, dt)
    !
    !     ! set the top boundary condition for CU winds to 0 to prevent artifacts coming in from the "top"
    !     domain%u_cu(:,nz,:) = 0
    !     domain%v_cu(:,nz,:) = 0
    !
    !     call advect3d(domain%u_cu, U_4cu_u,V_4cu_u,W_4cu_u, domain%rho, domain%dz_inter, nx+1,nz,ny, options)
    !     call advect3d(domain%v_cu, U_4cu_v,V_4cu_v,W_4cu_v, domain%rho, domain%dz_inter, nx,nz,ny+1, options)
    !
    ! end subroutine advect_cu_winds


    subroutine test_divergence(dz)
        implicit none
        real, intent(in) :: dz(ims:ime,kms:kme,jms:jme)

        real, allocatable :: du(:,:), dv(:,:), dw(:,:)
        integer :: i,j,k

        allocate(du(i_s:i_e,j_s:j_e))
        allocate(dv(i_s:i_e,j_s:j_e))
        allocate(dw(i_s:i_e,j_s:j_e))

        do i=i_s,i_e
            do j=j_s,j_e
                do k=kms,kme
                    du(i,j) = (U_m(i+1,k,j)-U_m(i,k,j))
                    dv(i,j) = (V_m(i,k,j+1)-V_m(i,k,j))
                    if (k==kms) then
                        dw(i,j) = (W_m(i,k,j))/dz(i,k,j)
                    else
                        dw(i,j) = (W_m(i,k,j)-W_m(i,k-1,j))/dz(i,k,j)
                    endif
                    if (abs(du(i,j) + dv(i,j) + dw(i,j)) > 1e-3) then
                        print*, this_image(), i,k,j , abs(du(i,j) + dv(i,j) + dw(i,j))
                        print*, "Winds are not balanced on entry to advect"
                        !error stop
                    endif
                enddo
            enddo
        enddo

    end subroutine test_divergence

    subroutine adv_std_compute_wind(domain, options, dt)
        implicit none

        type(options_t),    intent(in)  :: options
        type(domain_t),  intent(inout) :: domain
        real,intent(in)::dt
        
        integer :: i, j, k
        
        ! if (options%physics%convection > 0) then
            ! print*, "Advection of convective winds not enabled in ICAR >=1.5 yet"
            ! stop
            ! U_m = (domain%u_cu(2:nx,:,:) + domain%u(2:nx,:,:)) * (dt/dx)
            ! V_m = (domain%v_cu(:,:,2:ny) + domain%v(:,:,2:ny)) * (dt/dx)
            ! W_m = (domain%w_cu + domain%w)                     * (dt/dx)
            ! call rebalance_cu_winds(U_m,V_m,W_m)
        ! else
             ! Divide only U and V by dx. This minimizes the number of operations per advection step. W cannot be divided by dz,
             ! since non-uniform dz spacing does not allow for the same spacing to be assumed on either side of a k+1/2 interface,
             ! as is required for the adv4 scheme.
            dx = domain%dx
            
            if (options%parameters%advect_density) then
                do i = ims,ime
                    do j = jms,jme
                        do k = kms,kme
                            rho(i,k,j) = domain%density%data_3d(i,k,j)  
                        enddo
                    enddo
                enddo
            else
                do i = ims,ime
                    do j = jms,jme
                        do k = kms,kme
                            rho(i,k,j) = 1
                        enddo
                    enddo
                enddo
            endif
            do j = j_s,j_e
                do k = kms,kme
                    do i = i_s,i_e+1
                        U_m(i,k,j) = domain%u%data_3d(i,k,j) * dt * (rho(i,k,j)+rho(i-1,k,j))*0.5 * &
                                domain%jacobian_u(i,k,j) / domain%dx
                    enddo
                enddo
            enddo
            do j = j_s,j_e+1
                do k = kms,kme
                    do i = i_s,i_e
                        V_m(i,k,j) = domain%v%data_3d(i,k,j) * dt * (rho(i,k,j)+rho(i,k,j-1))*0.5 * &
                                domain%jacobian_v(i,k,j) / domain%dx
                    enddo
                enddo
            enddo
            do j = j_s,j_e
                do k = kms,kme-1
                    do i = i_s,i_e
                        W_m(i,k,j) = domain%w%data_3d(i,k,j) * dt * domain%jacobian_w(i,k,j) * &
                                                ( rho(i,k,j)*domain%advection_dz(i,k+1,j) + &
                                                  rho(i,k+1,j)*domain%advection_dz(i,k,j) ) / &
                                                 (domain%advection_dz(i,k,j)+domain%advection_dz(i,k+1,j))
                    enddo
                enddo
            enddo
            do j = j_s,j_e
                do i = i_s,i_e
                    W_m(i,kme,j) = domain%w%data_3d(i,kme,j) * dt * domain%jacobian_w(i,kme,j) * rho(i,kme,j)
                enddo
            enddo
    end subroutine adv_std_compute_wind

    subroutine setup_advection_dz(domain, options)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        integer :: i

        if (.not.allocated(domain%advection_dz)) then
            allocate(domain%advection_dz(ims:ime,kms:kme,jms:jme))
        else
            return
        endif

        do i=kms,kme
            domain%advection_dz(:,i,:) = options%parameters%dz_levels(i)
        enddo

    end subroutine setup_advection_dz


end module adv_std
