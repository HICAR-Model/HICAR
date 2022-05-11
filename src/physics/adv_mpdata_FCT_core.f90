!>---------------
!! MP-DATA flux correction core code
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!---------------
! For Reference, this is the core 1D implementation of mpdata included earlier in the code from adv_mpdata_core.f90
!
! It is stored in a separate file so that the same code can
! be used for advect_u, advect_v, and advect_w after l, r, etc. are set up
! 
!     ! First Calculate the standard upwind advection
!     call flux1(l,r,U2,f)
!
!     ! q1 is a temporary array to store the data while working with it here.
!     q1(1)=l(1)
!     q1(n)=r(n-1)
!     q1(2:n-1) = l(2:n-1) + (f(:n-2) - f(2:n-1))
!
!     ! This is the MPDATA diffusion correction term for 1D flow
!     ! U is defined on the mass grid for the pseudo-velocities?
!     ! left and right pseudo (diffusive) velocities
!
!     ! we will copy the q1 data into r to potentially minimize aliasing problems
!     ! for the compiler, and improve memory alignment for vectorization
!     r  = q1(2:n)
!     ! l  = q1(1:n-1) ! no need to copy these data over again
!
!     ! In MPDATA papers (r-l)/(r+l) is usually refered to as "A"
!     ! compute the denomenator first so we can check that it is not zero
!     denom=(r + q1(1:n-1))
!     where(denom==0) denom=1e-10
!     ! U2 is the diffusive pseudo-velocity
!     U2 = abs(U2) - U2**2
!     U2 = U2 * (r-q1(1:n-1)) / denom
!
!     ! now calculate the MPDATA flux term
!     call flux1(q1(1:n-1),r,U2,f)

! Fluxes are added to the original scalar field in the advect_u and advect_v subroutines

! This is the Flux Corrected Transport option described in : 
! Smolarkiewicz and Grabowski (1990) J. of Comp. Phys. v86 p355-375

! for now at least this is in a loop instead of vectorized.  I'm not sure how easy this would be to vectorize. 
    do i=n_s,n
        ! first find the min and max values allowable in the final field based on the initial (stored in l) and upwind (q1) fields
        ! min and max are taken from the grid cells on either side of the flux cell wall to be corrected
        if (i==n_s) then
            ! l still equals q0
            if (flux_is_w) then
                qmax_i=max(q1(i-1),q1(i),l(i-1),l(i))
                qmin_i=min(q1(i-1),q1(i),l(i-1),l(i))            
            else
                qmax_i=max(q1(i-2),q1(i-1),q1(i),l(i-2),l(i-1),l(i))
                qmin_i=min(q1(i-2),q1(i-1),q1(i),l(i-2),l(i-1),l(i))
            endif

        else
            ! l still equals q0
            qmax_i=qmax_i2
            qmin_i=qmin_i2
        endif
        if ( (flux_is_w) .and. (i==n) ) then
            qmax_i2=max(q1(i-1),q1(i),l(i-1),l(i))
            qmin_i2=min(q1(i-1),q1(i),l(i-1),l(i))
        else
            qmax_i2=max(q1(i-1),q1(i),q1(i+1),l(i-1),l(i),l(i+1))
            qmin_i2=min(q1(i-1),q1(i),q1(i+1),l(i-1),l(i),l(i+1))
        endif
        
            ! for the boundary, q1(i+1)==q0(i+1), l is only 1:n-1
        !    qmax_i=qmax_i2
        !    qmin_i=qmin_i2
        !    qmax_i2=max(q1(i-1),q1(i),l(i-1))
        !    qmin_i2=min(q1(i-1),q1(i),l(i-1))
        !endif

        fin_i = f_p(i-1)
        fout_i = f_n(i-1)
        
        fin_i2 = f_p(i)
        fout_i2 = f_n(i)
        
        ! if wind is left to right we limit based on flow out of the left cell and into the right cell
        if (U2(i)>0) then
            beta_out_i = (q1(i-1)-qmin_i) / (fout_i+1e-15)
            beta_in_i2 = (qmax_i2-q1(i)) / (fin_i2+1e-15)
            
            U2(i) = min(1.,beta_in_i2, beta_out_i) * U2(i)
            
        ! if wind is right to left we limit based on flow out of the right cell and into the left cell
        elseif (U2(i)<0) then
            beta_in_i = (qmax_i-q1(i-1)) / (fin_i+1e-15)
            beta_out_i2 = (q1(i)-qmin_i2) / (fout_i2+1e-15)
            
            U2(i) = min(1.,beta_in_i, beta_out_i2) * U2(i)
        endif
    end do
    
