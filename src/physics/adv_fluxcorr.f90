!> ----------------------------------------------------------------------------
!!  A collection of flux correction schemes for advection
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!! ----------------------------------------------------------------------------
module adv_fluxcorr
    use data_structures
    use domain_interface,  only: domain_t

    implicit none
    private
    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte

    public :: WRF_flux_corr, init_fluxcorr

contains

    subroutine init_fluxcorr(domain)
        implicit none
        type(domain_t), intent(in) :: domain
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
    end subroutine init_fluxcorr

    subroutine WRF_flux_corr(q,u,v,w,flux_x,flux_z,flux_y,flux_x_up,flux_z_up,flux_y_up,jaco,dz,rho)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in) :: q, jaco, dz, rho
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+1),  intent(in) :: w
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+1),  intent(in) :: u
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+2),  intent(in) :: v
        
        real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+1),intent(inout)          :: flux_x, flux_x_up
        real, dimension(its-1:ite+1,kms:kme,  jts-1:jte+2),intent(inout)          :: flux_y, flux_y_up
        real, dimension(its-1:ite+1,kms:kme+1,jts-1:jte+1),intent(inout)    :: flux_z, flux_z_up
        
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+1)   :: scale_in, scale_out
        real :: temp, flux_in, flux_out, dz_t_i, jaco_rho_t_i, fx, fx1, fy, fy1, fz, fz1, qmax, qmin
        integer :: i, j ,k
        
        !Initialize some internal variables
        do j = jts-1,jte+1
            do k = kms,kme
                scale_in(:,k,j) = 1.0
                scale_out(:,k,j) = 1.0
            enddo
        enddo
        ! Get upwind fluxes
        call upwind_flux3(q,u,v,w,flux_x_up,flux_z_up,flux_y_up,jaco,dz,rho)


        ! Next compute max and min possible fluxes
        do j = jts-1,jte+1
            do k = kms,kme
                !dir$ safe_address
                do i = its-1,ite+1
                
                   qmax = q(i,k,j)
                   qmin = q(i,k,j)

                    if (u(i,k,j) > 0) then
                        qmax = max(q(i-1,k,j),qmax)
                        qmin = min(q(i-1,k,j),qmin)
                    else if (u(i+1,k,j) < 0) then
                        qmax = max(q(i+1,k,j),qmax)
                        qmin = min(q(i+1,k,j),qmin)
                    endif

                    if (v(i,k,j) > 0) then
                        qmax = max(q(i,k,j-1),qmax)
                        qmin = min(q(i,k,j-1),qmin)
                    else if (v(i,k,j+1) < 0) then
                        qmax = max(q(i,k,j+1),qmax)
                        qmin = min(q(i,k,j+1),qmin)
                    endif
                    
                    if (w(i,k,j) < 0 .and. k < kme) then
                        qmax = max(q(i,k+1,j),qmax)
                        qmin = min(q(i,k+1,j),qmin)
                    else if (w(i,k,j) > 0 .and. k > kms) then
                        qmax = max(q(i,k-1,j),qmax)
                        qmin = min(q(i,k-1,j),qmin)
                    endif

                    !Store reused variables to minimize memory accesses
                    dz_t_i   = 1./dz(i,k,j)
                    jaco_rho_t_i = 1./(jaco(i,k,j)*rho(i,k,j))
                    fx = flux_x(i,k,j)-flux_x_up(i,k,j); fx1 = flux_x(i+1,k,j)-flux_x_up(i+1,k,j)
                    fy = flux_y(i,k,j)-flux_y_up(i,k,j); fy1 = flux_y(i,k,j+1)-flux_y_up(i,k,j+1)
                    fz = flux_z(i,k,j)-flux_z_up(i,k,j); fz1 = flux_z(i,k+1,j)-flux_z_up(i,k+1,j)
                    
                    !Now that the respective flux arrays should be loaded into the cache, update flux_x/y/z with - flux_x/y/z_up
                    !This is needed later on
                    flux_x(i,k,j) = fx; flux_y(i,k,j) = fy; flux_z(i,k,j) = fz; 
                    
                    !Compute concentration if upwind only was used
                    temp  = q(i,k,j) - ((flux_x_up(i+1,k,j) - flux_x_up(i,k,j)) + &
                                        (flux_y_up(i,k,j+1) - flux_y_up(i,k,j)) + &
                                        (flux_z_up(i,k+1,j) - flux_z_up(i,k,j)) * &
                                         dz_t_i)*jaco_rho_t_i

        
                    flux_in = -((min(0.,fx1) - max(0.,fx)) + &
                                (min(0.,fy1) - max(0.,fy)) + &
                                (min(0.,fz1) - max(0.,fz)) * &
                                 dz_t_i)*jaco_rho_t_i

                    flux_out = ((max(0.,fx1) - min(0.,fx)) + &
                                (max(0.,fy1) - min(0.,fy)) + &
                                (max(0.,fz1) - min(0.,fz)) * &
                                 dz_t_i)*jaco_rho_t_i

                    if (flux_in  > (qmax-temp))  then
                        scale_in(i,k,j) = max(0.,(qmax-temp)/(flux_in  + 1.E-15))
                    endif
                    if (flux_out > (temp-qmin)) then
                        scale_out(i,k,j) = max(0.,(temp-qmin)/(flux_out+ 1.E-15))
                    endif
                enddo
            enddo
        enddo

        !Apply upwind fluxes to calculated higher order fluxes -- these were not already done within the above loop
        flux_x(ite+2,kms:kme,jts-1:jte+1) = flux_x(ite+2,kms:kme,jts-1:jte+1) - flux_x_up(ite+2,kms:kme,jts-1:jte+1)
        flux_y(its-1:ite+1,kms:kme,jte+2) = flux_y(its-1:ite+1,kms:kme,jte+2) - flux_y_up(its-1:ite+1,kms:kme,jte+2)
        flux_z(its-1:ite+1,kme+1,jts-1:jte+1) = flux_z(its-1:ite+1,kme+1,jts-1:jte+1) - flux_z_up(its-1:ite+1,kme+1,jts-1:jte+1)

        do j = jts-1,jte+1
            do k = kms,kme
                do i = its,ite+1
                    if (flux_x(i,k,j) > 0) then
                        flux_x(i,k,j) = min(scale_in(i,k,j),scale_out(i-1,k,j))*flux_x(i,k,j)
                    else
                        flux_x(i,k,j) = min(scale_out(i,k,j),scale_in(i-1,k,j))*flux_x(i,k,j)
                    endif
                enddo
            enddo
        enddo
        do j = jts,jte+1
            do k = kms,kme
                do i = its-1,ite+1
                    if (flux_y(i,k,j) > 0) then
                        flux_y(i,k,j) = min(scale_in(i,k,j),scale_out(i,k,j-1))*flux_y(i,k,j)
                    else
                        flux_y(i,k,j) = min(scale_out(i,k,j),scale_in(i,k,j-1))*flux_y(i,k,j)
                    endif
                enddo
            enddo
        enddo
        do j = jts-1,jte+1
            do k = kms+1,kme
                do i = its-1,ite+1
                    if (flux_z(i,k,j) > 0) then
                        flux_z(i,k,j) = min(scale_in(i,k,j),scale_out(i,k-1,j))*flux_z(i,k,j)
                    else
                        flux_z(i,k,j) = min(scale_out(i,k,j),scale_in(i,k-1,j))*flux_z(i,k,j)
                    endif
                enddo
            enddo
        enddo

        !Finally, compute the output fluxes
        do j = jts-1,jte+1
            do k = kms,kme
                do i = its-1,ite+1
                    flux_x(i,k,j) = flux_x(i,k,j) + flux_x_up(i,k,j)
                    flux_y(i,k,j) = flux_y(i,k,j) + flux_y_up(i,k,j)
                    flux_z(i,k,j) = flux_z(i,k,j) + flux_z_up(i,k,j)
                enddo
            enddo
        enddo
        flux_x(ite+2,kms:kme,jts-1:jte+1) = flux_x(ite+2,kms:kme,jts-1:jte+1) + flux_x_up(ite+2,kms:kme,jts-1:jte+1)
        flux_y(its-1:ite+1,kms:kme,jte+2) = flux_y(its-1:ite+1,kms:kme,jte+2) + flux_y_up(its-1:ite+1,kms:kme,jte+2)
        flux_z(its-1:ite+1,kme+1,jts-1:jte+1) = flux_z(its-1:ite+1,kme+1,jts-1:jte+1) + flux_z_up(its-1:ite+1,kme+1,jts-1:jte+1)

    end subroutine WRF_flux_corr

    subroutine upwind_flux3(q,u,v,w,flux_x,flux_z,flux_y,jaco,dz,rho)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in)      :: q, jaco, dz, rho
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+1),  intent(in)    :: u
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+2),  intent(in)    :: v
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+1),  intent(in) :: w

        real, dimension(its-1:ite+2,kms:kme,jts-1:jte+1),intent(inout)          :: flux_x
        real, dimension(its-1:ite+1,kms:kme,jts-1:jte+2),intent(inout)          :: flux_y
        real, dimension(its-1:ite+1,kms:kme+1,jts-1:jte+1),intent(inout)    :: flux_z
        
        real, dimension(ims:ime,  kms:kme,jms:jme) :: dumb_q
        integer :: i, j, k
        real    :: tmp
        !When using RK3, we may have a time step derived using a CFL constraint larger than 1
        !This means that our upwind advection here may be in violation of the CFL criterion,
        !since this does not take place within the RK3 scheme. Since the possible CFL constraints
        !under RK3 with the available advection orders are all < 2, we can simply do 2 upwind steps

        !Upwind fluxes for first (half of a )step already calculated in advection flux function -- apply them first here

        !Update intermediate concentration
        do j = jms,jme
            do k = kms,kme
                dumb_q(:,k,j) = q(:,k,j)
            enddo
        enddo
        do j = jts-1,jte+1
            do k = kms,kme
                do i = its-1,ite+1
                    dumb_q(i,k,j)  = q(i,k,j) - ((flux_x(i+1,k,j) - flux_x(i,k,j)) + &
                                                 (flux_y(i,k,j+1) - flux_y(i,k,j)) + &
                                                 (flux_z(i,k+1,j) - flux_z(i,k,j)) / &
                                                 dz(i,k,j))/(jaco(i,k,j)*rho(i,k,j))
                enddo
            enddo
        enddo
        
        !Now compute upwind fluxes after second step
        do j = jts-1,jte+1
            do k = kms,kme
                do i = its-1,ite+2
                    tmp = u(i,k,j)
                    flux_x(i,k,j)= flux_x(i,k,j) + 0.5*((tmp + ABS(tmp)) * dumb_q(i-1,k,j) + (tmp - ABS(tmp)) * dumb_q(i,k,j)) * 0.5
                enddo
            enddo
        enddo
        do j = jts-1,jte+2
            do k = kms,kme
                do i = its-1,ite+1
                    tmp = v(i,k,j)
                    flux_y(i,k,j)= flux_y(i,k,j) + 0.5*((tmp + ABS(tmp)) * dumb_q(i,k,j-1) + (tmp - ABS(tmp)) * dumb_q(i,k,j)) * 0.5
                enddo
            enddo
        enddo
        do j = jts-1,jte+1
            do k = kms+1,kme
                do i = its-1,ite+1
                    tmp = w(i,k-1,j)
                    flux_z(i,k,j) = flux_z(i,k,j) + 0.5*((tmp + ABS(tmp)) * dumb_q(i,k-1,j) + &
                                                         (tmp - ABS(tmp)) * dumb_q(i,k,j)) * 0.5
                enddo
            enddo
        enddo
        !Handle top and bottom boundaries for z here
        do j = jts-1,jte+1
            do i = its-1,ite+1
                flux_z(i,kme+1,j) = flux_z(i,kme+1,j) + 0.5*dumb_q(i,kme,j) * w(i,kme,j)
            enddo
        enddo
                                        
    end subroutine upwind_flux3
end module adv_fluxcorr