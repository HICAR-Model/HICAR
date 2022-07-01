!> ----------------------------------------------------------------------------
!!  The MPDATA advection scheme
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module adv_mpdata
    use data_structures
    use options_interface, only: options_t
    use domain_interface,  only: domain_t

    implicit none
    private
    real,dimension(:,:,:),allocatable::U_m,V_m,W_m, rho
    integer :: order, ims, ime, jms, jme, kms, kme, its,ite,jts,jte
        
    public:: mpdata, mpdata_init, mpdata_advect3d, mpdata_compute_wind
    

contains

    subroutine flux1(l,r,U,f)
    !     Calculate the donor cell flux function
    !     l = left gridcell scalar 
    !     r = right gridcell scalar
    !     U = Courant number (u*dt/dx)
    !     
    !     If U is positive, return l*U if U is negative return r*U
    !     By using the mathematical form instead of the logical form, 
    !     we can run on the entire grid simultaneously, and avoid branches

    !   arguments
        implicit none
        real, dimension(:), intent(in) :: l,r,U
        real, dimension(:), intent(inout) :: f
        
        !   main code
        f= ((U+ABS(U)) * l + (U-ABS(U)) * r)/2

    end subroutine flux1
    
    subroutine flux3(q,u,v,w,flux_x,flux_z,flux_y)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in) :: q, w
        real, dimension(ims+1:ime,  kms:kme,jms:jme),  intent(in) :: u
        real, dimension(ims:ime,  kms:kme,jms+1:jme),  intent(in) :: v
        
        real, dimension(ims+1:ime,kms:kme,jms:jme),intent(inout)          :: flux_x
        real, dimension(ims:ime,kms:kme,jms+1:jme),intent(inout)          :: flux_y
        real, dimension(ims+1:ime-1,kms:kme+1,jms+1:jme-1),intent(inout)  :: flux_z

        flux_x= ((u(ims+1:ime,:,:) + ABS(u(ims+1:ime,:,:))) * q(ims:ime-1,:,:)    + &
            (u(ims+1:ime,:,:) - ABS(u(ims+1:ime,:,:))) * q(ims+1:ime,:,:))   / 2

        flux_y= ((v(:,:,jms+1:jme)    + ABS(v(:,:,jms+1:jme)))    * q(:,:,jms:jme-1) + &
            (v(:,:,jms+1:jme)    - ABS(v(:,:,jms+1:jme)))    * q(:,:,jms+1:jme))  / 2

        flux_z(:,kms+1:kme,:) = ((w(ims+1:ime-1,kms:kme-1,jms+1:jme-1) + ABS(w(ims+1:ime-1,kms:kme-1,jms+1:jme-1))) &
                                        * q(ims+1:ime-1,kms:kme-1,jms+1:jme-1) + &
            (w(ims+1:ime-1,kms:kme-1,jms+1:jme-1) - ABS(w(ims+1:ime-1,kms:kme-1,jms+1:jme-1))) &
                                         * q(ims+1:ime-1,kms+1:kme,jms+1:jme-1))  / 2
                                         
        !Handle top and bottom boundaries for z here
        flux_z(:,kms,:) = 0
        flux_z(:,kme+1,:) = q(ims+1:ime-1,kme,jms+1:jme-1) * w(ims+1:ime-1,kme,jms+1:jme-1)

                                         
    end subroutine flux3
    !-------------------------------------------------------------------------------------------------
    !Function expects u,v,w to all be multiplied by dt, and u and v to be divided by dx
    !-------------------------------------------------------------------------------------------------
        subroutine upwind_advection(qin, u, v, w, q, dz, jaco)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in) :: qin, w 
        real, dimension(ims+1:ime,  kms:kme,jms:jme),  intent(in) :: u
        real, dimension(ims:ime,  kms:kme,jms+1:jme),  intent(in) :: v
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(inout) :: q
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in) :: jaco, dz
        
        ! interal parameters
        integer :: i
        real, dimension(ims+1:ime,kms:kme,jms:jme)    :: flux_x
        real, dimension(ims:ime,kms:kme,jms+1:jme)    :: flux_y
        real, dimension(ims+1:ime-1,kms:kme+1,jms+1:jme-1)  :: flux_z
        
        ! !$omp parallel shared(qin,q,u,v,w) firstprivate(nx,ny,nz) private(i,f1,f3,f4,f5)
        ! !$omp do schedule(static)
        !do i=jms,jme
        !    q(:,:,i)=qin(:,:,i)
        !enddo
        
        q = qin
        
        ! !$omp end do
        ! !$omp barrier
        ! !$omp do schedule(static)
            ! by manually inlining the flux2 call we should remove extra array copies that the compiler doesn't remove.
            ! equivalent flux2 calls are left in for reference (commented) to restore recall that f1,f3,f4... arrays should be 3D : n x m x 1

        call flux3(qin,u,v,w,flux_x,flux_z,flux_y)

        ! perform horizontal advection, from difference terms
        q(its:ite,:,jts:jte)  = q(its:ite,:,jts:jte)  - &
                                   ((flux_x(its+1:ite+1,:,jts:jte) - flux_x(its:ite,:,jts:jte)) + &
                                   (flux_y(its:ite,:,jts+1:jte+1) - flux_y(its:ite,:,jts:jte))) &
                                   / (jaco(its:ite,:,jts:jte)*rho(its:ite,:,jts:jte))                      
               ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
               ! add fluxes to middle layers
        q(its:ite,:,jts:jte) = q(its:ite,:,jts:jte)  &
                                   - (flux_z(its:ite,kms+1:kme+1,jts:jte) - flux_z(its:ite,kms:kme,jts:jte)) &
                                   / (dz(its:ite,:,jts:jte)*jaco(its:ite,:,jts:jte)*rho(its:ite,:,jts:jte))

        ! !$omp end do
        ! !$omp end parallel
        
    end subroutine upwind_advection

    subroutine mpdata_fluxes(q,u,v,w,u2,v2,w2,G)
        implicit none
        real, dimension(ims:ime,kms:kme,jms:jme),   intent(in) :: q,w
        real, dimension(ims+1:ime,kms:kme,jms:jme), intent(in) :: u
        real, dimension(ims:ime,kms:kme,jms+1:jme), intent(in) :: v
        real, dimension(ims+1:ime,kms:kme,jms:jme), intent(out) :: u2
        real, dimension(ims:ime,kms:kme,jms+1:jme), intent(out) :: v2
        real, dimension(ims:ime,kms:kme,jms:jme),   intent(out) :: w2
        real, dimension(ims:ime,kms:kme,jms:jme),   intent(in) :: G  ! Notation kept from SMOLARKIEWICZ AND MARGOLIN 1998
        
        real, dimension(ims+1:ime) :: rx, lx, denomx
        real, dimension(ims:ime) :: r, l, denom, edge_v, edge_q
        integer :: i, j
        
        u2 = 0
        v2 = 0
        w2 = 0
        ! This might run faster if tiled over x and y to be more cache friendly. 
        ! !$omp parallel shared(q,u,v,w,u2,v2,w2) firstprivate(nx,ny,nz) &
        ! !$omp private(i,j, rx,lx,r,l, denomx,denom)
        ! !$omp do schedule(static)
        do i=jms,jme
            do j=kms,kme
                ! -----------------------
                ! First compute the U component
                ! -----------------------
                if (i > 0) then !((i>jms).and.(i<jme)) then
                    rx=q(ims+1:ime,j,i)
                    lx=q(ims:ime-1,j,i)
                    ! In MPDATA papers (r-l)/(r+l) is usually refered to as "A"
                    ! compute the denomenator first so we can check that it is not zero
                    denomx=(rx + lx + 1e-10)
                    !where(denomx==0) denomx=1e-10
                    ! U2 is the diffusive pseudo-velocity
                    u2(:,j,i) = abs(u(:,j,i))*(1-abs(u(:,j,i))/(0.5* ( G(ims+1:ime,j,i) + G(ims:ime-1,j,i) )))
                    u2(:,j,i) = u2(:,j,i) * (rx-lx) / denomx
                    
                    edge_v = 0
                    edge_q = 0
                    ! UxV terms
                    if ( (i > jms) .and. (i < jme) ) then
                        edge_q(ims+1:ime) = (q(ims+1:ime,j,i+1) - q(ims+1:ime,j,i-1) + &
                             q(ims:ime-1,j,i+1) - q(ims:ime-1,j,i-1)) / &
                            (q(ims+1:ime,j,i+1) + q(ims+1:ime,j,i-1) + &
                             q(ims:ime-1,j,i+1) + q(ims:ime-1,j,i-1)+ 1e-10)
                        edge_v(ims+1:ime) = (1/4.0)*(v(ims+1:ime,j,i) + v(ims+1:ime,j,i+1) + v(ims:ime-1,j,i) + v(ims:ime-1,j,i+1))
                        u2(:,j,i) = u2(:,j,i) - u(:,j,i)*edge_v(ims+1:ime)*edge_q(ims+1:ime)/( G(ims+1:ime,j,i) + G(ims:ime-1,j,i) )
                    endif
                    
                    edge_v = 0
                    edge_q = 0
                    ! UxW terms
                    if ( (j > kms) .and. (j < kme) ) then
                        edge_q(ims+1:ime) = (q(ims+1:ime,j+1,i) - q(ims+1:ime,j-1,i) + &
                             q(ims:ime-1,j+1,i) - q(ims:ime-1,j-1,i)) / &
                            (q(ims+1:ime,j+1,i) + q(ims+1:ime,j-1,i) + &
                             q(ims:ime-1,j+1,i) + q(ims:ime-1,j-1,i)+ 1e-10)
                        edge_v(ims+1:ime) = (1/4.0)*(w(ims+1:ime,j,i) + w(ims+1:ime,j-1,i) + w(ims:ime-1,j,i) + w(ims:ime-1,j-1,i))
                        u2(:,j,i) = u2(:,j,i) - u(:,j,i)*edge_v(ims+1:ime)*edge_q(ims+1:ime)/( G(ims+1:ime,j,i) + G(ims:ime-1,j,i) )
                    endif
                endif
                

                ! next compute the V and W components
                if (i>jms) then
                    ! -----------------------
                    ! compute the V component
                    ! -----------------------
                    r=q(:,j,i)
                    l=q(:,j,i-1)
                    ! In MPDATA papers A = (r-l)/(r+l)
                    ! compute the denomenator first so we can check that it is not zero
                    denom=(r + l + 1e-10)
                    !where(denom==0) denom=1e-10
                    ! U2 is the diffusive pseudo-velocity
                    v2(:,j,i) = abs(v(:,j,i))*(1-abs(v(:,j,i))/(0.5* ( G(:,j,i) + G(:,j,i-1) )))
                    v2(:,j,i) = v2(:,j,i) * (r-l) / denom
                    
                    edge_v = 0
                    edge_q = 0
                    !VxU terms
                    edge_q(ims+1:ime-1) = (q(ims+2:ime,j,i-1) - q(ims:ime-2,j,i) + &
                             q(ims+2:ime,j,i) - q(ims:ime-2,j,i-1)) / &
                            (q(ims+2:ime,j,i) + q(ims+2:ime,j,i-1) + &
                             q(ims:ime-2,j,i) + q(ims:ime-2,j,i-1)+ 1e-10)
                    edge_v(ims+1:ime-1) = (1/4.0)* &
                                        (u(ims+2:ime,j,i) + u(ims+2:ime,j,i-1) + u(ims+1:ime-1,j,i) + u(ims+1:ime-1,j,i-1))
                    v2(:,j,i) = v2(:,j,i) - v(:,j,i)*edge_v*edge_q/( G(:,j,i) + G(:,j,i-1) )
                    
                    edge_v = 0
                    edge_q = 0
                    ! VxW terms
                    if ( (j > kms) .and. (j < kme) ) then
                        edge_q = (q(:,j+1,i-1) - q(:,j-1,i) + &
                             q(:,j+1,i) - q(:,j-1,i-1)) / &
                            (q(:,j+1,i-1) + q(:,j-1,i) + &
                             q(:,j+1,i) + q(:,j-1,i-1)+ 1e-10)
                        edge_v = (1/4.0)*(w(:,j,i) + w(:,j-1,i) + w(:,j,i-1) + w(:,j-1,i-1))
                        v2(:,j,i) = v2(:,j,i) - v(:,j,i)*edge_v*edge_q/( G(:,j,i) + G(:,j,i-1) )
                    endif
                endif
                
                
                ! -----------------------
                ! compute the w component
                ! -----------------------
                if (j==kme) then
                    w2(:,j,i)=0
                else
                    r=q(:,j+1,i)
                    l=q(:,j,i)
                    ! In MPDATA papers A = (r-l)/(r+l)
                    ! compute the denomenator first so we can check that it is not zero
                    denom=(r + l + 1e-10)
                    !where(denom==0) denom=1e-10
                    ! U2 is the diffusive pseudo-velocity
                    w2(:,j,i) = abs(w(:,j,i))*(1-abs(w(:,j,i))/(0.5* ( G(:,j+1,i) + G(:,j,i) )))
                    w2(:,j,i) = w2(:,j,i) * (r-l) / denom
                    
                    edge_v = 0
                    edge_q = 0
                    !WxU terms
                    edge_q(ims+1:ime-1) = (q(ims+2:ime,j+1,i) - q(ims:ime-2,j,i) + &
                             q(ims+2:ime,j,i) - q(ims:ime-2,j+1,i)) / &
                            (q(ims+2:ime,j,i) + q(ims+2:ime,j+1,i) + &
                             q(ims:ime-2,j,i) + q(ims:ime-2,j+1,i)+ 1e-10)
                    edge_v(ims+1:ime-1) = (1/4.0)* &
                                        (u(ims+2:ime,j,i) + u(ims+2:ime,j+1,i) + u(ims+1:ime-1,j,i) + u(ims+1:ime-1,j+1,i))
                    w2(:,j,i) = w2(:,j,i) - w(:,j,i)*edge_v*edge_q/( G(:,j+1,i) + G(:,j,i) )
                    
                    edge_v = 0
                    edge_q = 0
                    ! WxV terms
                    if ( (i > jms) .and. (i < jme) ) then
                        edge_q = (q(:,j+1,i+1) - q(:,j,i-1) + &
                             q(:,j,i+1) - q(:,j+1,i-1)) / &
                            (q(:,j,i+1) + q(:,j+1,i-1) + &
                             q(:,j+1,i+1) + q(:,j,i-1)+ 1e-10)
                        edge_v = (1/4.0)*(v(:,j,i) + v(:,j+1,i) + v(:,j,i+1) + v(:,j+1,i+1))
                        w2(:,j,i) = w2(:,j,i) - w(:,j,i)*edge_v*edge_q/( G(:,j+1,i) + G(:,j,i) )
                    endif
                endif
            end do
        end do
        ! !$omp end do
        ! !$omp end parallel
        
    end subroutine mpdata_fluxes
    
    subroutine netpos(f_pos,flux_x,flux_z,flux_y,dz,dx,dt,G_in)
        implicit none
        real, dimension(ims+1:ime-1,kms:kme,jms+1:jme-1), intent(inout) :: f_pos
        real, dimension(ims+1:ime,kms:kme,jms:jme),         intent(in)  :: flux_x
        real, dimension(ims:ime,kms:kme,jms+1:jme),         intent(in)  :: flux_y
        real, dimension(ims+1:ime-1,kms:kme+1,jms+1:jme-1), intent(in)  :: flux_z
        real,dimension(ims:ime,kms:kme,jms:jme),  intent(in)            :: dz, G_in

        real,    intent(in) :: dt, dx    
    
    
        f_pos =         (max(0.,flux_x(ims+1:ime-1,:,jms+1:jme-1))       + max(0.,abs(flux_x(ims+2:ime,:,jms+1:jme-1))))             * dt / &
                                                        (dx*G_in(ims+1:ime-1,:,jms+1:jme-1))
        f_pos = f_pos + (max(0.,flux_z(ims+1:ime-1,kms:kme,jms+1:jme-1)) + max(0.,abs(flux_z(ims+1:ime-1,kms+1:kme+1,jms+1:jme-1)))) * dt / &
                                                        (dz(ims+1:ime-1,:,jms+1:jme-1)*G_in(ims+1:ime-1,:,jms+1:jme-1))
        f_pos = f_pos + (max(0.,flux_y(ims+1:ime-1,:,jms+1:jme-1))       + max(0.,abs(flux_y(ims+1:ime-1,:,jms+2:jme))))             * dt / &
                                                        (dx*G_in(ims+1:ime-1,:,jms+1:jme-1))
    
    end subroutine netpos
    
    subroutine netneg(f_neg,flux_x,flux_z,flux_y,dz,dx,dt,G_in)
        implicit none
        real, dimension(ims+1:ime-1,kms:kme,jms+1:jme-1), intent(inout) :: f_neg
        real, dimension(ims+1:ime,kms:kme,jms:jme),         intent(in)  :: flux_x
        real, dimension(ims:ime,kms:kme,jms+1:jme),         intent(in)  :: flux_y
        real, dimension(ims+1:ime-1,kms:kme+1,jms+1:jme-1), intent(in)  :: flux_z
        real,dimension(ims:ime,kms:kme,jms:jme),  intent(in)            :: dz, G_in

        real,    intent(in) :: dt, dx    
    
    
        f_neg =         (max(0.,-flux_x(ims+1:ime-1,:,jms+1:jme-1))       + max(0.,abs(-flux_x(ims+2:ime,:,jms+1:jme-1))))             * dt / &
                                                        (dx*G_in(ims+1:ime-1,:,jms+1:jme-1))
        f_neg = f_neg + (max(0.,-flux_z(ims+1:ime-1,kms:kme,jms+1:jme-1)) + max(0.,abs(-flux_z(ims+1:ime-1,kms+1:kme+1,jms+1:jme-1)))) * dt / &
                                                        (dz(ims+1:ime-1,:,jms+1:jme-1)*G_in(ims+1:ime-1,:,jms+1:jme-1))
        f_neg = f_neg + (max(0.,-flux_y(ims+1:ime-1,:,jms+1:jme-1))       + max(0.,abs(-flux_y(ims+1:ime-1,:,jms+2:jme))))             * dt / &
                                                        (dx*G_in(ims+1:ime-1,:,jms+1:jme-1))
    
    end subroutine netneg    

    subroutine flux_limiter(q, q2, u,v,w, dz,dx,dt,G_in)
        implicit none
        real,dimension(ims:ime,kms:kme,jms:jme),  intent(in)    :: q, q2, dz, G_in
        real,dimension(ims+1:ime,kms:kme,jms:jme),intent(inout) :: u
        real,dimension(ims:ime,kms:kme,jms+1:jme),intent(inout) :: v
        real,dimension(ims:ime,kms:kme,jms:jme),  intent(inout) :: w
        real,    intent(in) :: dt, dx
        
        integer :: i,j,k,n,n_s
        real, dimension(:), pointer :: q1, U2, l, f_p, f_n
        real, dimension(ims+1:ime,kms:kme,jms:jme)    :: flux_x
        real, dimension(ims:ime,kms:kme,jms+1:jme)    :: flux_y
        real, dimension(ims+1:ime-1,kms:kme+1,jms+1:jme-1)  :: flux_z
        real, dimension(ims+1:ime-1,kms:kme,jms+1:jme-1)    :: f_pos
        real, dimension(ims+1:ime-1,kms:kme,jms+1:jme-1)    :: f_neg
        
        ! q1 = q after applying previous iteration advection
        ! l  = q before applying previous iteration
        ! U2 is the anti-diffusion pseudo-velocity
        ! f is the first pass calculation of MPDATA fluxes
        real, dimension(ims:ime),     target :: q1x,lx
        real, dimension(ims+1:ime),   target :: U2x
        real, dimension(ims+1:ime-1), target :: f_px,f_nx
        real, dimension(jms:jme),     target :: q1y,ly
        real, dimension(jms+1:jme),   target :: U2y
        real, dimension(jms+1:jme-1), target :: f_py,f_ny
        real, dimension(kms:kme),     target :: f_pz,f_nz, U2z, q1z,lz
        logical :: flux_is_w
        
        real :: qmax_i,qmin_i,qmax_i2,qmin_i2
        real :: beta_in_i, beta_out_i, beta_in_i2, beta_out_i2
        real :: fin_i, fout_i, fin_i2, fout_i2
        
        
        call flux3(q,u,v,w,flux_x,flux_z,flux_y)
        
        call netpos(f_pos,flux_x,flux_z,flux_y,dz,dx,dt,G_in)
        call netneg(f_neg,flux_x,flux_z,flux_y,dz,dx,dt,G_in)

        ! NOTE: before inclusion of FCT_core the following variables must be setup: 
        ! q1 and l (l=q0)
        ! !$omp parallel shared(q2,q,u,v,w,dz,G_in) firstprivate(ims,ime,kms,kme,jms,jme,dx,dt) default(private)
        ! !$omp do schedule(static)
        do j=jms+1,jme-1
            flux_is_w=.False.
            q1=>q1x
            l =>lx
            U2=>U2x
            n=ime-1
            n_s=ims+2
            f_p =>f_px
            f_n =>f_nx
            do k=kms,kme
                ! setup u
                q1=q2(:,k,j)
                U2=u(:,k,j)
                l =q(:,k,j)
                f_p = f_pos(:,k,j)
                f_n = f_neg(:,k,j)

                include "adv_mpdata_FCT_core.f90"
                u(:,k,j)=U2
            end do
            
            !K indices for z are a bit offset, such that the staggering notation
            !used in FCT core still works for w
            q1=>q1z
            l =>lz
            U2=>U2z
            n=kme
            n_s=kms+1
            f_p =>f_pz
            f_n =>f_nz
            flux_is_w=.True.
            do k=ims+1,ime-1
                ! setup w
                q1=q2(k,:,j)
                U2(kms+1:kme)=w(k,kms:kme-1,j)
                l =q(k,:,j)
                f_p = f_pos(k,:,j)
                f_n = f_neg(k,:,j)
                
                include "adv_mpdata_FCT_core.f90"
                w(k,kms:kme-1,j)=U2(kms+1:kme)
                w(k,kme,j)=0
            end do
            
        end do
        ! !$omp end do
        
        flux_is_w=.False.
        q1=>q1y
        l =>ly
        U2=>U2y
        n=jme-1
        n_s=jms+2
        f_p =>f_py
        f_n =>f_ny
        ! NOTE: This it typically not the correct order for the loop variables
        ! but in this case it permits parallelization over a larger number (nx instead of nz)
        ! and because all data are copied from an oddly spaced grid regardless, it *probably* doesn't slow it down
        ! I'd like to re-write the v-flux delimiter to operate on all x simulataneously at some point...
        ! !$omp do
        do j=ims+1,ime-1
            do k=kms,kme
                q1=q2(j,k,:)
                U2=v(j,k,:)
                l =q(j,k,:)
                f_p = f_pos(j,k,:)
                f_n = f_neg(j,k,:)
                
                include "adv_mpdata_FCT_core.f90"
                v(j,k,:)=U2
            end do
        end do
        ! !$omp end do
        ! !$omp end parallel
        
    end subroutine flux_limiter

    subroutine mpdata_advect3d(q,jaco,dz,dx,dt,options)
        implicit none
        real,dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: q
        real,dimension(ims:ime,kms:kme,jms:jme), intent(in) :: jaco,dz
        real, intent(in)    :: dt, dx
        type(options_t), intent(in)::options

        ! used for intermediate values in the mpdata calculation
        real,dimension(ims+1:ime,kms:kme,jms:jme) :: u2
        real,dimension(ims:ime,kms:kme,jms+1:jme) :: v2
        real,dimension(ims:ime,kms:kme,jms:jme)   :: w2
        real,dimension(ims:ime,kms:kme,jms:jme)   :: q2

        integer :: iord, i
        
        call upwind_advection(q, U_m, V_m, W_m, q2,dz,jaco)
        
        do iord=1,options%adv_options%mpdata_order
            if (iord>1) then
                ! Due to an uneven vertical grid spacing, W_m is not normalized by dz, because this causes problems in
                ! advection code. However, pseudo-velocity expects, all velocities normalized by time and dx/z
                ! so do that here before passing to pseudo-velocity calculations
                call mpdata_fluxes(q2, U_m, V_m, (W_m/dz), u2,v2,w2,(jaco*rho))
                ! and un-normalize, since upwind advection scheme includes dz
                ! Since pseudo-velocities cannot be gaurenteed to be non-divergent, we assume worst-case and multiply by 0.5 to
                ! ensure stability (from Smolarkiewicz 1984, after Eq. 24)
                if (options%adv_options%flux_corrected_transport) then
                    !These should be pure velocities, no dt/dx normalization
                    u2 = u2*dx/dt
                    v2 = v2*dx/dt
                    w2 = w2*dz/dt
                    call flux_limiter(q, q2, u2,v2,w2, dz,dx,dt,(jaco*rho))
                    !Advection expects all normalized, but w not divided by dz
                    u2 = u2*dt/dx
                    v2 = v2*dt/dx
                    w2 = w2*dt
                else
                    w2 = w2*dz
                endif
                call upwind_advection(q2, u2,v2,w2, q,dz,jaco)
            endif
            
            !  
            if (iord/=options%adv_options%mpdata_order) then
                if (iord>1) then
                    ! !$omp parallel shared(q,q2) firstprivate(ny) private(i)
                    ! !$omp do schedule(static)
                    do i=jms,jme
                        q2(:,:,i)=q(:,:,i)
                    enddo
                    ! !$omp end do
                    ! !$omp end parallel
                endif
            else 
                if (iord==1) then
                    ! !$omp parallel shared(q,q2) firstprivate(ny) private(i)
                    ! !$omp do schedule(static)
                    do i=jms,jme
                        q(:,:,i)=q2(:,:,i)
                    enddo
                    ! !$omp end do
                    ! !$omp end parallel
                endif
                
            endif
        end do
        

    end subroutine mpdata_advect3d
    
    subroutine mpdata_init(domain,options)
        type(domain_t), intent(in) :: domain
        type(options_t), intent(in) :: options
        
        ! originally used to permit the order of dimensions in advection to be rotated
        order    = 0
        
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
        
!       if this if the first time we are called, we need to allocate the module level arrays
        if (.not.allocated(U_m)) then
            allocate(U_m     (ims+1:ime,kms:kme,jms:jme  ))
            allocate(V_m     (ims:ime,  kms:kme,jms+1:jme))
            allocate(W_m     (ims:ime,  kms:kme,jms:jme  ))
        endif
        
        allocate(rho(ims:ime,  kms:kme,jms:jme  ))
        
    end subroutine mpdata_init
    
    
    subroutine test_divergence(dz)
        implicit none
        real, intent(in) :: dz(ims:ime,kms:kme,jms:jme)

        real, allocatable :: du(:,:), dv(:,:), dw(:,:)
        integer :: i,j,k

        allocate(du(ims+1:ime-1,jms+1:jme-1))
        allocate(dv(ims+1:ime-1,jms+1:jme-1))
        allocate(dw(ims+1:ime-1,jms+1:jme-1))

        do i=ims+1,ime-1
            do j=jms+1,jme-1
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
    
    subroutine mpdata_compute_wind(domain,options,dt)
        implicit none
        type(domain_t),intent(inout)::domain
        type(options_t), intent(in)::options
        real,intent(in)::dt
        
        real::dx
        integer :: i

        dx=domain%dx


        rho = 1
        if (options%parameters%advect_density) rho = domain%density%data_3d
        
        U_m = domain%u%data_3d(ims+1:ime,:,:) * dt * (rho(ims+1:ime,:,:)+rho(ims:ime-1,:,:))*0.5 * &
                    domain%jacobian_u(ims+1:ime,:,:) / dx
        V_m = domain%v%data_3d(:,:,jms+1:jme) * dt * (rho(:,:,jms+1:jme)+rho(:,:,jms:jme-1))*0.5 * &
                    domain%jacobian_v(:,:,jms+1:jme) / dx
        W_m(:,kms:kme-1,:) = domain%w%data_3d(:,kms:kme-1,:) * dt * domain%jacobian_w(:,kms:kme-1,:) * &
                    (rho(:,kms+1:kme,:)+rho(:,kms:kme-1,:)) * 0.5
        W_m(:,kme,:) = domain%w%data_3d(:,kme,:) * dt * domain%jacobian_w(:,kme,:) * rho(:,kme,:)

        if (options%parameters%debug) then
            call test_divergence(domain%advection_dz)
        endif
    end subroutine mpdata_compute_wind
        
!   primary entry point, advect all scalars in domain
    subroutine mpdata(domain,options,dt)
        implicit none
        type(domain_t),intent(inout)::domain
        type(options_t), intent(in)::options
        real,intent(in)::dt
        
        real::dx
        integer :: i
        
        dx=domain%dx

        if (.not.allocated(domain%advection_dz)) then
            allocate(domain%advection_dz(ims:ime,kms:kme,jms:jme))
            do i=kms,kme
                domain%advection_dz(:,i,:) = options%parameters%dz_levels(i)
            enddo
        endif
        
!       if this if the first time we are called, we need to allocate the module level arrays
        if (.not.allocated(U_m)) then
            allocate(U_m     (ims+1:ime,kms:kme,jms:jme  ))
            allocate(V_m     (ims:ime,  kms:kme,jms+1:jme))
            allocate(W_m     (ims:ime,  kms:kme,jms:jme  ))
        endif
        
        allocate(rho(ims:ime,  kms:kme,jms:jme  ))
        rho = 1
        if (options%parameters%advect_density) rho = domain%density%data_3d
        
        U_m = domain%u%data_3d(ims+1:ime,:,:) * dt * (rho(ims+1:ime,:,:)+rho(ims:ime-1,:,:))*0.5 * &
                    domain%jacobian_u(ims+1:ime,:,:) / dx
        V_m = domain%v%data_3d(:,:,jms+1:jme) * dt * (rho(:,:,jms+1:jme)+rho(:,:,jms:jme-1))*0.5 * &
                    domain%jacobian_v(:,:,jms+1:jme) / dx
        W_m(:,kms:kme-1,:) = domain%w%data_3d(:,kms:kme-1,:) * dt * domain%jacobian_w(:,kms:kme-1,:) * &
                    (rho(:,kms+1:kme,:)+rho(:,kms:kme-1,:)) * 0.5
        W_m(:,kme,:) = domain%w%data_3d(:,kme,:) * dt * domain%jacobian_w(:,kme,:) * rho(:,kme,:)

        if (options%parameters%debug) then
            call test_divergence(domain%advection_dz)
        endif

        !if (options%vars_to_advect(kVARS%water_vapor)>0)                  call mpdata_advect3d(domain%water_vapor%data_3d,             domain%jacobian, domain%advection_dz, dx,dt,options)
        !if (options%vars_to_advect(kVARS%cloud_water)>0)                  call mpdata_advect3d(domain%cloud_water_mass%data_3d,        domain%jacobian, domain%advection_dz, dx,dt,options)
        !if (options%vars_to_advect(kVARS%rain_in_air)>0)                  call mpdata_advect3d(domain%rain_mass%data_3d,               domain%jacobian, domain%advection_dz, dx,dt,options)
        !if (options%vars_to_advect(kVARS%snow_in_air)>0)                  call mpdata_advect3d(domain%snow_mass%data_3d,               domain%jacobian, domain%advection_dz, dx,dt,options)
        !if (options%vars_to_advect(kVARS%potential_temperature)>0)        call mpdata_advect3d(domain%potential_temperature%data_3d, domain%jacobian, domain%advection_dz, dx,dt,options)
        !if (options%vars_to_advect(kVARS%cloud_ice)>0)                    call mpdata_advect3d(domain%cloud_ice_mass%data_3d,          domain%jacobian, domain%advection_dz, dx,dt,options)
        !if (options%vars_to_advect(kVARS%graupel_in_air)>0)               call mpdata_advect3d(domain%graupel_mass%data_3d,            domain%jacobian, domain%advection_dz, dx,dt,options)
        !if (options%vars_to_advect(kVARS%ice_number_concentration)>0)     call mpdata_advect3d(domain%cloud_ice_number%data_3d,        domain%jacobian, domain%advection_dz, dx,dt,options)
        !if (options%vars_to_advect(kVARS%rain_number_concentration)>0)    call mpdata_advect3d(domain%rain_number%data_3d,             domain%jacobian, domain%advection_dz, dx,dt,options)
        !if (options%vars_to_advect(kVARS%snow_number_concentration)>0)    call mpdata_advect3d(domain%snow_number%data_3d,             domain%jacobian, domain%advection_dz, dx,dt,options)
        !if (options%vars_to_advect(kVARS%graupel_number_concentration)>0) call mpdata_advect3d(domain%graupel_number%data_3d,          domain%jacobian, domain%advection_dz, dx,dt,options)

    end subroutine mpdata
end module adv_mpdata
