!>------------------------------------------------------------
!! Module to solve for a 3D wind field following mass-conservation
!! and reducing differennces between initial and final wind field.
!! Solver requires use of PETSc software package, and so is 
!! separated from the rest of wind core here to allow for dependancy
!! on PETSc library.
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!------------------------------------------------------------
!#include <petsc/finclude/petscdmda.h90>


module wind_iterative
    !include 'petsc/finclude/petscksp.h'
    !include 'petsc/finclude/petscdm.h'
    !include 'petsc/finclude/petscdmda.h'
    
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

!#include <petsc/finclude/petscdmstag.h>

    !use exchangeable_interface,   only : exchangeable_t
    use domain_interface,  only : domain_t
    !use options_interface, only : options_t
    !use grid_interface,    only : grid_t
    use mod_wrf_constants, only : epsilon
    use petscksp
    use petscdm
    use petscdmda
!    use petscdmstag
    use io_routines,          only : io_write, io_read
    use array_utilities,      only : smooth_array_3d

    implicit none
    private
    public:: init_iter_winds, calc_iter_winds, finalize_iter_winds
    real, parameter::deg2rad=0.017453293 !2*pi/360
    real, parameter :: rad2deg=57.2957779371
    real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, H_coef, I_coef, &
                                           J_coef, K_coef, L_coef, M_coef, N_coef, O_coef, tmp_x, tmp_y
    real    :: dx
    real, allocatable, dimension(:,:,:)  :: div, dz_if, jaco, dzhatdxdz, dzhatdydz, dzhatdzz, dzdxz, dzdyz, d2zdx2, d2zdy2, dzdx, dzdy, sigma, alpha
    real, allocatable, dimension(:,:)    :: w_0, dzdx_surf, dzdy_surf
    integer, allocatable :: xl(:), yl(:)
    integer              :: hs, i_s, i_e, k_s, k_e, j_s, j_e, ims, ime, kms, kme, jms, jme
contains



    !>------------------------------------------------------------
    !! Forces u,v, and w fields to balance
    !!       du/dx + dv/dy = dw/dz
    !!
    !! Starts by setting w out of the ground=0 then works through layers
    !!
    !!------------------------------------------------------------
    subroutine calc_iter_winds(domain,alpha_in,div_in,w_surf,adv_den,update_in)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, intent(in), dimension(domain%grid%ims:domain%grid%ime, &
                                    domain%grid%kms:domain%grid%kme, &
                                    domain%grid%jms:domain%grid%jme) :: alpha_in, div_in, w_surf
        logical, intent(in) :: adv_den
        logical, optional, intent(in) :: update_in

        PetscScalar,pointer :: lambda(:,:,:)
        logical             :: update

        integer k !, i_s, i_e, k_s, k_e, j_s, j_e
        
        PetscErrorCode ierr
        KSP            ksp
        PC             pc
        DM             da
        Vec            x, localX
        PetscInt       one, two, x_size, iteration
        PetscReal      norm, conv_tol
        KSPConvergedReason reason
        

        update=.False.
        if (present(update_in)) update=update_in
                
        call calc_RHS(domain%u%meta_data%dqdt_3d,domain%v%meta_data%dqdt_3d,domain%w%meta_data%dqdt_3d, &
                                domain%jacobian_u, domain%jacobian_v,domain%jacobian_w,domain%jacobian, &
                                domain%advection_dz,domain%dx, domain%dzdx, domain%dzdy, domain%density%data_3d,adv_den)

        !Initialize div to be the initial divergence of the input wind field
        div = div_in(ims:ime,kms:kme,jms:jme) 
        one = 1
        two = 2

        
        alpha = alpha_in(i_s:i_e,k_s:k_e,j_s:j_e)
        w_0 = w_surf(i_s:i_e,k_s,j_s:j_e)
        
        if (.not.(allocated(A_coef))) then
            call initialize_coefs(domain)
        else
            call update_coefs(domain)
        endif
                                                                
        call KSPCreate(domain%IO_comms,ksp,ierr)
        conv_tol = 1e-21
        call KSPSetFromOptions(ksp,ierr)

        !call KSPSetTolerances(ksp,conv_tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        call KSPSetType(ksp,KSPBCGSL,ierr);
        !call KSPGetPC(ksp, pc, ierr)
        !call PCSetType(pc, PCBJACOBI, ierr)
        
        call DMDACreate3d(domain%IO_comms,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, &
                          domain%ide,(domain%kde+2),domain%jde,domain%grid%ximages,one,domain%grid%yimages,one,two, &
                          xl, PETSC_NULL_INTEGER,yl,da,ierr)
        
        call DMSetFromOptions(da,ierr)
        call DMSetUp(da,ierr)
        
        call KSPSetDM(ksp,da,ierr)
        call KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,0,ierr)
        call KSPSetComputeRHS(ksp,ComputeRHS,0,ierr)
        call KSPSetComputeOperators(ksp,ComputeMatrix,0,ierr)

        call DMCreateLocalVector(da,localX,ierr)
        call KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)
        
        call KSPGetSolution(ksp,x,ierr)
        call KSPGetConvergedReason(ksp, reason, ierr)
        call KSPGetIterationNumber(ksp, iteration, ierr)
        if(this_image()==1) write(*,*) 'Solved PETSc after ',iteration,' iterations, reason: ',reason
        
        !Subset global solution x to local grid so that we can access ghost-points
        call DMGlobalToLocalBegin(da,x,INSERT_VALUES,localX,ierr)
        call DMGlobalToLocalEnd(da,x,INSERT_VALUES,localX,ierr)

        call DMDAVecGetArrayF90(da,localX,lambda, ierr)
        if ( this_image()==1 .and. reason < 0) then
            write (*,*) 'PETSc ERROR: convergence failed after ',iteration,' iterations'
            write (*,*) 'PETSc KSP Error code:  ',reason
            !stop
        endif
        call calc_updated_winds(domain, lambda, update, adv_den)
        call DMDAVecRestoreArrayF90(da,localX,lambda, ierr)

        !Exchange u and v, since the outer points are not updated in above function
        call domain%u%exchange_x(update)
        call domain%v%exchange_y(update)
        !call smooth_array_3d( domain%u%meta_data%dqdt_3d, windowsize = 1, ydim = 3)
        !call smooth_array_3d( domain%v%meta_data%dqdt_3d, windowsize = 1, ydim = 3)
        
        call VecDestroy(localX,ierr)
        call DMDestroy(da,ierr)
        call KSPDestroy(ksp,ierr)
        !call finalize_iter_winds()
                
    end subroutine calc_iter_winds
        
    subroutine calc_updated_winds(domain,lambda,update,adv_den) !u, v, w, jaco_u,jaco_v,jaco_w,u_dzdx,v_dzdy,lambda, ids, ide, jds, jde)
        type(domain_t), intent(inout) :: domain
        !real, intent(inout), dimension(:,:,:)  :: u,v,w
        !real, intent(in), dimension(:,:,:)     :: jaco_u,jaco_v,jaco_w, u_dzdx, v_dzdy
        PetscScalar, intent(in), pointer       :: lambda(:,:,:)
        logical,     intent(in)                :: update
        logical,     intent(in)                :: adv_den


        real, allocatable, dimension(:,:,:)    :: u_dlambdz, v_dlambdz, dlambdadz, dlambdadx, dlambdady, rho, rho_u, rho_v, u_cor, v_cor, tmp
        integer :: k, i_start, i_end, j_start, j_end !i_s, i_e, k_s, k_e, j_s, j_e, ids, ide, jds, jde
                        
        !i_s+hs, unless we are on global boundary, then i_s
        i_start = i_s+1
        if (i_s==domain%grid%ids) i_start = i_s
        
        !i_e, unless we are on global boundary, then i_e+1
        i_end = i_e
        if (i_e==domain%grid%ide) i_end = i_e+1
        
        !j_s+hs, unless we are on global boundary, then j_s
        j_start = j_s+1
        if (j_s==domain%grid%jds) j_start = j_s
        
        !j_e, unless we are on global boundary, then j_e+1
        j_end = j_e
        if (j_e==domain%grid%jde) j_end = j_e+1

        allocate(dlambdadz(i_s:i_e,k_s:k_e,j_s:j_e))
        
        allocate(tmp(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(dlambdadx(i_start:i_end,k_s:k_e,j_s:j_e))
        allocate(dlambdady(i_s:i_e,k_s:k_e,j_start:j_end))
        allocate(u_cor(i_s:i_end,k_s:k_e,j_s:j_e))
        allocate(v_cor(i_s:i_e,k_s:k_e,j_s:j_end))

        allocate(u_dlambdz(i_start:i_end,k_s:k_e,j_s:j_e))
        allocate(v_dlambdz(i_s:i_e,k_s:k_e,j_start:j_end))
        
        allocate(rho(domain%ims:domain%ime,k_s:k_e,domain%jms:domain%jme))
        allocate(rho_u(i_start:i_end,k_s:k_e,j_s:j_e))
        allocate(rho_v(i_s:i_e,k_s:k_e,j_start:j_end))

        rho = 1.0
        rho_u = 1.0
        rho_v = 1.0
        
        if (adv_den) rho(domain%ims:domain%ime,:,domain%jms:domain%jme)=domain%density%data_3d(domain%ims:domain%ime,:,domain%jms:domain%jme)
        
        if (i_s==domain%grid%ids .and. i_e==domain%grid%ide) then
            rho_u(i_start+1:i_end-1,:,j_s:j_e) = 0.5*(rho(i_start+1:i_end-1,:,j_s:j_e) + rho(i_start:i_end-2,:,j_s:j_e))
            rho_u(i_end,:,j_s:j_e) = 1.5*rho(i_end-1,:,j_s:j_e) - 0.5*rho(i_end-2,:,j_s:j_e)
            rho_u(i_start,:,j_s:j_e) = 1.5*rho(i_start,:,j_s:j_e) - 0.5*rho(i_start+1,:,j_s:j_e)
        else if (i_s==domain%grid%ids) then
            rho_u(i_start+1:i_end,:,j_s:j_e) = 0.5*(rho(i_start+1:i_end,:,j_s:j_e) + rho(i_start:i_end-1,:,j_s:j_e))
            rho_u(i_start,:,j_s:j_e) = 1.5*rho(i_start,:,j_s:j_e) - 0.5*rho(i_start+1,:,j_s:j_e)
        else if (i_e==domain%grid%ide) then
            rho_u(i_start:i_end-1,:,j_s:j_e) = 0.5*(rho(i_start:i_end-1,:,j_s:j_e) + rho(i_start-1:i_end-2,:,j_s:j_e))
            rho_u(i_end,:,j_s:j_e) = 1.5*rho(i_end-1,:,j_s:j_e) - 0.5*rho(i_end-2,:,j_s:j_e)
        else
            rho_u(i_start:i_end,:,j_s:j_e) = 0.5*(rho(i_start:i_end,:,j_s:j_e) + rho(i_start-1:i_end-1,:,j_s:j_e))
        endif
        
        if (j_s==domain%grid%jds .and. j_e==domain%grid%jde) then
            rho_v(i_s:i_e,:,j_start+1:j_end-1) = 0.5*(rho(i_s:i_e,:,j_start+1:j_end-1) + rho(i_s:i_e,:,j_start:j_end-2))
            rho_v(i_s:i_e,:,j_start) = 1.5*rho(i_s:i_e,:,j_start) - 0.5*rho(i_s:i_e,:,j_start+1)
            rho_v(i_s:i_e,:,j_end) = 1.5*rho(i_s:i_e,:,j_end-1) - 0.5*rho(i_s:i_e,:,j_end-2)
        else if (j_s==domain%grid%jds) then
            rho_v(i_s:i_e,:,j_start+1:j_end) = 0.5*(rho(i_s:i_e,:,j_start+1:j_end) + rho(i_s:i_e,:,j_start:j_end-1))
            rho_v(i_s:i_e,:,j_start) = 1.5*rho(i_s:i_e,:,j_start) - 0.5*rho(i_s:i_e,:,j_start+1)
        else if (j_e==domain%grid%jde) then
            rho_v(i_s:i_e,:,j_start:j_end-1) = 0.5*(rho(i_s:i_e,:,j_start:j_end-1) + rho(i_s:i_e,:,j_start-1:j_end-2))
            rho_v(i_s:i_e,:,j_end) = 1.5*rho(i_s:i_e,:,j_end-1) - 0.5*rho(i_s:i_e,:,j_end-2)
        else
            rho_v(i_s:i_e,:,j_start:j_end) = 0.5*(rho(i_s:i_e,:,j_start:j_end) + rho(i_s:i_e,:,j_start-1:j_end-1))
        endif
        
        !PETSc arrays are zero-indexed        

        !divide dz differennces by dz. Note that dz will be horizontally constant
        do k=k_s,k_e
            dlambdadz(:,k,:) = (lambda(i_s-1:i_e-1,k+1,j_s-1:j_e-1)*sigma(i_s,k,j_s)**2) - (sigma(i_s,k,j_s)**2 - 1)*lambda(i_s-1:i_e-1,k,j_s-1:j_e-1) - lambda(i_s-1:i_e-1,k-1,j_s-1:j_e-1)
            dlambdadz(:,k,:) = dlambdadz(:,k,:)/(dz_if(i_s,k+1,j_s)*(sigma(i_s,k,j_s)+sigma(i_s,k,j_s)**2))
            !dlambdadz(:,k,:) = (lambda(i_s-1:i_e-1,k+1,j_s-1:j_e-1) - lambda(i_s-1:i_e-1,k,j_s-1:j_e-1))/dz_if(i_s,k+1,j_s)
        enddo
        
        !dlambdadz(:,k_s,:) = -(lambda(i_s-1:i_e-1,k_s+2,j_s-1:j_e-1)*sigma(i_s,k_s+1,j_s)**2) + &
        !                    lambda(i_s-1:i_e-1,k_s+1,j_s-1:j_e-1)*(sigma(i_s,k_s+1,j_s)+1)**2 - lambda(i_s-1:i_e-1,k_s,j_s-1:j_e-1)*(2*sigma(i_s,k_s+1,j_s)+1)        
        !dlambdadz(:,k_s,:) = dlambdadz(:,k_s,:)/(dz_if(i_s,k_s+1,j_s)*(sigma(i_s,k_s+1,j_s)+1))
        
        dlambdadx(i_s+1:i_e,k_s:k_e,j_s:j_e) = (lambda(i_s:i_e-1,k_s:k_e,j_s-1:j_e-1) - lambda(i_s-1:i_e-2,k_s:k_e,j_s-1:j_e-1))/dx
        dlambdady(i_s:i_e,k_s:k_e,j_s+1:j_e) = (lambda(i_s-1:i_e-1,k_s:k_e,j_s:j_e-1) - lambda(i_s-1:i_e-1,k_s:k_e,j_s-1:j_e-2))/dx

        !stager dlambdadz to u grid
        u_dlambdz(i_s+1:i_e,k_s:k_e,j_s:j_e) = (dlambdadz(i_s+1:i_e,k_s:k_e,j_s:j_e) + &
                                                dlambdadz(i_s:i_e-1,k_s:k_e,j_s:j_e)) / 2
        !stager dlambdadz to v grid
        v_dlambdz(i_s:i_e,k_s:k_e,j_s+1:j_e) = (dlambdadz(i_s:i_e,k_s:k_e,j_s+1:j_e) + &
                                                dlambdadz(i_s:i_e,k_s:k_e,j_s:j_e-1)) / 2


        !stager dlambdadz to u grid
        !tmp(i_s+1:i_e,k_s:k_e,j_s:j_e) = max(min(domain%dzdx(i_s+1:i_e,k_s:k_e,j_s:j_e)/(domain%dzdx(i_s:i_e-1,k_s:k_e,j_s:j_e) + &
        !                                                                                   domain%dzdx(i_s+1:i_e,k_s:k_e,j_s:j_e)+0.0001),1.0),0.0)
        !u_dlambdz(i_s+1:i_e,k_s:k_e,j_s:j_e) = (dlambdadz(i_s+1:i_e,k_s:k_e,j_s:j_e)*tmp(i_s+1:i_e,k_s:k_e,j_s:j_e) + &
        !                                        dlambdadz(i_s:i_e-1,k_s:k_e,j_s:j_e)*(1-tmp(i_s+1:i_e,k_s:k_e,j_s:j_e)))
        !!!stager dlambdadz to v grid
        !tmp(i_s:i_e,k_s:k_e,j_s+1:j_e) = max(min(domain%dzdy(i_s:i_e,k_s:k_e,j_s+1:j_e)/(domain%dzdy(i_s:i_e,k_s:k_e,j_s:j_e-1) + &
        !                                                                                   domain%dzdy(i_s:i_e,k_s:k_e,j_s+1:j_e)+0.0001),1.0),0.0)
        !v_dlambdz(i_s:i_e,k_s:k_e,j_s+1:j_e) = (dlambdadz(i_s:i_e,k_s:k_e,j_s+1:j_e)*tmp(i_s:i_e,k_s:k_e,j_s+1:j_e) + &
        !                                        dlambdadz(i_s:i_e,k_s:k_e,j_s:j_e-1)*(1-tmp(i_s:i_e,k_s:k_e,j_s+1:j_e)))

        
        if (i_s==domain%grid%ids) then
            u_dlambdz(i_start,:,:) = 0*u_dlambdz(i_start+1,:,:)
            dlambdadx(i_start,:,:) = 0*dlambdadx(i_start+1,:,:)
        endif
        if (i_e==domain%grid%ide) then
            u_dlambdz(i_end,:,:) = 0*u_dlambdz(i_end-1,:,:)
            dlambdadx(i_end,:,:) = 0*dlambdadx(i_end-1,:,:)
        endif
        if (j_s==domain%grid%jds) then
            v_dlambdz(:,:,j_start) = 0*v_dlambdz(:,:,j_start+1)
            dlambdady(:,:,j_start) = 0*dlambdady(:,:,j_start+1)
        endif
        if (j_e==domain%grid%jde) then
            v_dlambdz(:,:,j_end) = 0*v_dlambdz(:,:,j_end-1)
            dlambdady(:,:,j_end) = 0*dlambdady(:,:,j_end-1)
        endif

        !domain%froude%data_3d(i_s:i_e,k_s:k_e,j_s:j_e) = d2zdx2(i_s:i_e,k_s:k_e,j_s:j_e)
        domain%Ri%data_3d(i_s:i_e,k_s:k_e,j_s:j_e) = lambda(i_s-1:i_e-1,k_s-1:k_e-1,j_s-1:j_e-1)
        
        if (update) then
            domain%u%meta_data%dqdt_3d(i_start:i_end,:,j_s:j_e) = domain%u%meta_data%dqdt_3d(i_start:i_end,:,j_s:j_e) + 0.5*(dlambdadx - &
                        domain%dzdx_u(i_start:i_end,:,j_s:j_e)*(u_dlambdz/domain%jacobian_u(i_start:i_end,:,j_s:j_e))) / &
                        (rho_u(i_start:i_end,:,j_s:j_e))!*domain%jacobian_u(i_start:i_end,:,j_s:j_e))
             
            domain%v%meta_data%dqdt_3d(i_s:i_e,:,j_start:j_end) = domain%v%meta_data%dqdt_3d(i_s:i_e,:,j_start:j_end) + 0.5*(dlambdady - &
                        domain%dzdy_v(i_s:i_e,:,j_start:j_end)*(v_dlambdz/domain%jacobian_v(i_s:i_e,:,j_start:j_end))) / &
                        (rho_v(i_s:i_e,:,j_start:j_end))!*domain%jacobian_v(i_s:i_e,:,j_start:j_end))
                        
            !domain%u%meta_data%dqdt_3d(i_start:i_end,:,j_s:j_e) = domain%u%meta_data%dqdt_3d(i_start:i_end,:,j_s:j_e) + u_cor(i_start:i_end,:,j_s:j_e)/rho_u(i_start:i_end,:,j_s:j_e)
            
            !domain%v%meta_data%dqdt_3d(i_s:i_e,:,j_start:j_end) = domain%v%meta_data%dqdt_3d(i_s:i_e,:,j_start:j_end) + v_cor(i_s:i_e,:,j_start:j_end)/rho_v(i_s:i_e,:,j_start:j_end)
            domain%w%meta_data%dqdt_3d(i_s:i_e,k_s:k_e,j_s:j_e) = domain%w%meta_data%dqdt_3d(i_s:i_e,k_s:k_e,j_s:j_e) + alpha(i_s:i_e,k_s:k_e,j_s:j_e)**2*0.5* &
                        ((lambda(i_s-1:i_e-1,k_s+1:k_e+1,j_s-1:j_e-1)-lambda(i_s-1:i_e-1,k_s:k_e,j_s-1:j_e-1))/(dz_if(i_s:i_e,k_s+1:k_e+1,j_s:j_e)*domain%jacobian_w(i_s:i_e,k_s:k_e,j_s:j_e)**2))
                        
            !domain%w%meta_data%dqdt_3d(i_s:i_e,k_s:k_e,j_s:j_e) = domain%w%meta_data%dqdt_3d(i_s:i_e,k_s:k_e,j_s:j_e) + alpha(i_s:i_e,k_s:k_e,j_s:j_e)**2*0.5* &
            !            (dlambdadz(i_s:i_e,k_s:k_e,j_s:j_e)/domain%jacobian_w(i_s:i_e,k_s:k_e,j_s:j_e))

        else
            domain%u%data_3d(i_start:i_end,:,j_s:j_e) = domain%u%data_3d(i_start:i_end,:,j_s:j_e) + 0.5*(dlambdadx - &
                        domain%dzdx_u(i_start:i_end,:,j_s:j_e)*(u_dlambdz/domain%jacobian_u(i_start:i_end,:,j_s:j_e))) / &
                        (rho_u(i_start:i_end,:,j_s:j_e))
            
            domain%v%data_3d(i_s:i_e,:,j_start:j_end) = domain%v%data_3d(i_s:i_e,:,j_start:j_end) + 0.5*(dlambdady - &
                        domain%dzdy_v(i_s:i_e,:,j_start:j_end)*(v_dlambdz/domain%jacobian_v(i_s:i_e,:,j_start:j_end))) / &
                        (rho_v(i_s:i_e,:,j_start:j_end))
        endif


    end subroutine calc_updated_winds
    
    
    subroutine calc_RHS(u, v, w, jaco_u, jaco_v, jaco_w, jaco, dz, dx, dzdx, dzdy, rho, adv_den)
        implicit none
        real, dimension(ims:ime,kms:kme,jms:jme),   intent(in)    :: w, dz, jaco_w, rho, jaco, dzdx, dzdy
        real, dimension(ims:ime+1,kms:kme,jms:jme), intent(in)    :: u, jaco_u
        real, dimension(ims:ime,kms:kme,jms:jme+1), intent(in)    :: v, jaco_v
        real,           intent(in)    :: dx
        logical,        intent(in)    :: adv_den
        
        real, dimension(ims:ime,kms:kme,jms:jme) :: diff_U, diff_V, w_met, dudz, dvdz, rho_i, tmp
        real, dimension(ims:ime+1,kms:kme,jms:jme) :: u_met
        real, dimension(ims:ime,kms:kme,jms:jme+1) :: v_met
        integer :: k
        
        !Multiplication of U/V by metric terms, converting jacobian to staggered-grid where possible, otherwise making assumption of
        !Constant jacobian at edges
        
        u_met = u * jaco_u
        v_met = v * jaco_v
        w_met = w * jaco_w
        
        if (adv_den) then
            u_met(ims+1:ime,:,:) = u_met(ims+1:ime,:,:) * (rho(ims:ime-1,:,:) + rho(ims+1:ime,:,:))/2
            u_met(ims,:,:) = u_met(ims,:,:) * (1.5*rho(ims,:,:) - 0.5*rho(ims+1,:,:))
            u_met(ime+1,:,:) = u_met(ime+1,:,:) * (1.5*rho(ime,:,:) - 0.5*rho(ime-1,:,:))

            v_met(:,:,jms+1:jme) = v_met(:,:,jms+1:jme) * (rho(:,:,jms:jme-1) + rho(:,:,jms+1:jme))/2
            v_met(:,:,jms) = v_met(:,:,jms) * (1.5*rho(:,:,jms) - 0.5*rho(:,:,jms+1))
            v_met(:,:,jme+1) = v_met(:,:,jme+1) * (1.5*rho(:,:,jme) - 0.5*rho(:,:,jme-1))

            rho_i(:,kms:kme-1,:) = ( rho(:,kms:kme-1,:)*dz(:,kms+1:kme,:) + rho(:,kms+1:kme,:)*dz(:,kms:kme-1,:) ) / (dz(:,kms:kme-1,:)+dz(:,kms+1:kme,:))
            rho_i(:,kme,:) = rho(:,kme,:)
            
            w_met = w * rho_i

            tmp(ims:ime,:,:) = (u(ims+1:ime+1,:,:)+u(ims:ime,:,:))*0.5
            tmp(:,kms:kme-1,:) = ( tmp(:,kms:kme-1,:)*dz(:,kms+1:kme,:) + tmp(:,kms+1:kme,:)*dz(:,kms:kme-1,:) ) / (dz(:,kms:kme-1,:)+dz(:,kms+1:kme,:))
            tmp(:,kme,:) = tmp(:,kme,:)
            
            dudz(:,kms+1:kme,:) = jaco_w(:,kms+1:kme,:)*tmp(:,kms+1:kme,:)*rho_i(:,kms+1:kme,:) - jaco_w(:,kms:kme-1,:)*tmp(:,kms:kme-1,:)*rho_i(:,kms:kme-1,:)
            dudz(:,kms,:) = jaco_w(:,kms,:)*tmp(:,kms,:)*rho_i(:,kms,:)
            
            tmp(:,:,jms:jme) = (v(:,:,jms+1:jme+1)+v(:,:,jms:jme))*0.5
            tmp(:,kms:kme-1,:) = ( tmp(:,kms:kme-1,:)*dz(:,kms+1:kme,:) + tmp(:,kms+1:kme,:)*dz(:,kms:kme-1,:) ) / (dz(:,kms:kme-1,:)+dz(:,kms+1:kme,:))
            tmp(:,kme,:) = tmp(:,kme,:)
            
            dvdz(:,kms+1:kme,:) = jaco_w(:,kms+1:kme,:)*tmp(:,kms+1:kme,:)*rho_i(:,kms+1:kme,:) - jaco_w(:,kms:kme-1,:)*tmp(:,kms:kme-1,:)*rho_i(:,kms:kme-1,:)
            dvdz(:,kms,:) = jaco_w(:,kms,:)*tmp(:,kms,:)*rho_i(:,kms,:)
        else
            tmp(ims:ime,:,:) = (u(ims+1:ime+1,:,:)+u(ims:ime,:,:))*0.5
            tmp(:,kms:kme-1,:) = ( tmp(:,kms:kme-1,:)*dz(:,kms+1:kme,:) + tmp(:,kms+1:kme,:)*dz(:,kms:kme-1,:) ) / (dz(:,kms:kme-1,:)+dz(:,kms+1:kme,:))
            tmp(:,kme,:) = tmp(:,kme,:)
            
            tmp = tmp*jaco_w
            dudz(:,kms+1:kme,:) = tmp(:,kms+1:kme,:) - tmp(:,kms:kme-1,:)
            dudz(:,kms,:) = tmp(:,kms,:)
            
            tmp(:,:,jms:jme) = (v(:,:,jms+1:jme+1)+v(:,:,jms:jme))*0.5
            tmp(:,kms:kme-1,:) = ( tmp(:,kms:kme-1,:)*dz(:,kms+1:kme,:) + tmp(:,kms+1:kme,:)*dz(:,kms:kme-1,:) ) / (dz(:,kms:kme-1,:)+dz(:,kms+1:kme,:))
            tmp(:,kme,:) = tmp(:,kme,:)
            
            tmp = tmp*jaco_w
            dvdz(:,kms+1:kme,:) = tmp(:,kms+1:kme,:) - tmp(:,kms:kme-1,:)
            dvdz(:,kms,:) = tmp(:,kms,:)
        end if

        diff_U = (u_met(ims+1:ime+1, :, jms:jme) - u_met(ims:ime, :, jms:jme))/dx
        diff_V = (v_met(ims:ime, :, jms+1:jme+1) - v_met(ims:ime, :, jms:jme))/dx
        
        div = diff_U+diff_V
        
        do k = kms,kme
            if (k == kms) then
                div(ims:ime, k, jms:jme) = div(ims:ime, k, jms:jme) + 0*w_met(ims:ime, k, jms:jme)/(dz(ims:ime, k, jms:jme))
            else
                div(ims:ime, k, jms:jme) = div(ims:ime, k, jms:jme) + 0*(w_met(ims:ime,k,jms:jme)-w_met(ims:ime,k-1,jms:jme))/(dz(ims:ime,k,jms:jme))
            endif
        enddo
        
        div = div - dzdx*dudz/(dz*jaco) - dzdy*dvdz/(dz*jaco)

    
    end subroutine calc_RHS
    

    subroutine ComputeRHS(ksp,vec_b,dummy,ierr)
        implicit none
        
        PetscErrorCode  ierr
        KSP ksp
        Vec vec_b, x
        integer dummy(*)
        DM             dm
        DMDALocalInfo       :: info(DMDA_LOCAL_INFO_SIZE)
        
        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs
        PetscScalar,pointer :: barray(:,:,:)

        call KSPGetDM(ksp,dm,ierr)
        
        !call DMStagGetCornersF90(dm,xs,zs,ys,xm,zm,ym,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER)
        !call DMStagGetGlobalSizesF90(dm,mx,mz,my)

        call DMDAGetLocalInfoF90(dm,info,ierr)
        !swap z and y in info calls since dims 2/3 are transposed for PETSc w.r.t. HICAR indexing
        mx = info(DMDA_LOCAL_INFO_MX)
        my = info(DMDA_LOCAL_INFO_MZ)
        mz = info(DMDA_LOCAL_INFO_MY)
        xm = info(DMDA_LOCAL_INFO_XM)
        ym = info(DMDA_LOCAL_INFO_ZM)
        zm = info(DMDA_LOCAL_INFO_YM)
        xs = info(DMDA_LOCAL_INFO_XS)
        ys = info(DMDA_LOCAL_INFO_ZS)
        zs = info(DMDA_LOCAL_INFO_YS)

        call DMDAVecGetArrayF90(dm,vec_b,barray, ierr)
                
        do j=ys,(ys+ym-1)
            do k=zs,(zs+zm-1)
                do i=xs,(xs+xm-1)
                    !For global boundary conditions
                    
                    ! Neumann BC at upper boundary (ghost point)
                    if (k.eq.mz-1 .and. i.ge.1 .and. i.le.mx-2 .and. &
                             j.ge.1 .and. j.le.my-2) then
                        barray(i,k,j) = 0.0
                    ! Neumann BC at lower boundary (ghost point)
                    else if (k.eq.0) then
                        barray(i,k,j) = 0*w_0(i+1,j+1) / &
                                (alpha(i+1,k+1,j+1)**2+dzdx(i+1,k+1,j+1)**2+dzdy(i+1,k+1,j+1)**2)
                    ! Equation for interior points and boundary points with Neumann BC
                    else if ((k.ge.2 .and. k.le.mz-2 .and. &
                             i.ge.1 .and. i.le.mx-2 .and. &
                             j.ge.1 .and. j.le.my-2) .or. &
                             (k.eq.1)) then
                        barray(i,k,j) = -2*div(i+1,k,j+1)!*jaco(i+1,k,j+1)
                    ! Dirlect BC for lateral boundaries
                    else
                        barray(i,k,j) = 0.0
                    endif
                enddo
            enddo
        enddo

        call DMDAVecRestoreArrayF90(dm,vec_b,barray, ierr)
    end subroutine ComputeRHS

    subroutine ComputeInitialGuess(ksp,vec_b,ctx,ierr)
        implicit none
        PetscErrorCode  ierr
        KSP ksp
        PetscInt ctx(*)
        Vec vec_b
        PetscScalar  i_guess

        i_guess = 0.0

        call VecSet(vec_b,i_guess,ierr)
    end subroutine ComputeInitialGuess

    subroutine ComputeMatrix(ksp,arr_A,arr_B,dummy,ierr)
        implicit none
        PetscErrorCode  ierr
        KSP ksp
        Mat arr_A,arr_B
        integer dummy(*)
        real denom


        DM             dm
        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,i1,i3,i7,i15
        PetscScalar    v(15)
        MatStencil     row(4),col(4,15),gnd_col(4,7),top_col(4,2)
        DMDALocalInfo       :: info(DMDA_LOCAL_INFO_SIZE)
        
        i1 = 1
        i3 = 2
        i7 = 7
        i15 = 15

        call KSPGetDM(ksp,dm,ierr)

        !call DMStagGetCornersF90(dm,xs,zs,ys,xm,zm,ym,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER)
        !call DMStagGetGlobalSizesF90(dm,mx,mz,my)
        call DMDAGetLocalInfoF90(dm,info,ierr)
        !swap z and y in info calls since dims 2/3 are transposed for PETSc w.r.t. HICAR indexing
        mx = info(DMDA_LOCAL_INFO_MX)
        my = info(DMDA_LOCAL_INFO_MZ)
        mz = info(DMDA_LOCAL_INFO_MY)
        xm = info(DMDA_LOCAL_INFO_XM)
        ym = info(DMDA_LOCAL_INFO_ZM)
        zm = info(DMDA_LOCAL_INFO_YM)
        xs = info(DMDA_LOCAL_INFO_XS)
        ys = info(DMDA_LOCAL_INFO_ZS)
        zs = info(DMDA_LOCAL_INFO_YS)
                
        do j=ys,(ys+ym-1)
        do k=zs,(zs+zm-1)
            do i=xs,(xs+xm-1)
                row(MatStencil_i) = i
                row(MatStencil_j) = k
                row(MatStencil_k) = j
                ! Neumann BC at ghost point rows
                if (k.eq.0 .or. (k.eq.mz-1 .and. i.ge.1 .and. i.le.mx-2 .and. &
                             j.ge.1 .and. j.le.my-2))  then
                    !at upper boundary (ghost point)
                    if (k.eq.mz-1) then
                        !k
                        !v(1) = sigma(i+1,k-1,j+1)**2/(dz_if(i+1,k,j+1)*(sigma(i+1,k-1,j+1)+sigma(i+1,k-1,j+1)**2))
                        !top_col(MatStencil_i,1) = i
                        !top_col(MatStencil_j,1) = k
                        !top_col(MatStencil_k,1) = j
                        !!k - 1
                        !v(2) = -(sigma(i+1,k-1,j+1)**2-1)/(dz_if(i+1,k,j+1)*(sigma(i+1,k-1,j+1)+sigma(i+1,k-1,j+1)**2))
                        !top_col(MatStencil_i,2) = i
                        !top_col(MatStencil_j,2) = k-1
                        !top_col(MatStencil_k,2) = j
                        !!k - 2
                        !v(3) = -1/(dz_if(i+1,k,j+1)*(sigma(i+1,k-1,j+1)+sigma(i+1,k-1,j+1)**2))
                        !top_col(MatStencil_i,3) = i
                        !top_col(MatStencil_j,3) = k-2
                        !top_col(MatStencil_k,3) = j
                        
                        v(1) = 1/dz_if(i+1,k,j+1)
                        top_col(MatStencil_i,1) = i
                        top_col(MatStencil_j,1) = k
                        top_col(MatStencil_k,1) = j
                        !k + 1
                        v(2) = -1/dz_if(i+1,k,j+1)
                        top_col(MatStencil_i,2) = i
                        top_col(MatStencil_j,2) = k-1
                        top_col(MatStencil_k,2) = j

                        
                        call MatSetValuesStencil(arr_B,i1,row,i3,top_col,v,INSERT_VALUES, ierr)
                    !lower boundary (ghost point)
                    else if (k.eq.0) then
                        denom = (dzdx_surf(i+1,j+1)**2 + dzdy_surf(i+1,j+1)**2 + &
                                        alpha(i+1,1,j+1)**2)!/(jaco(i+1,1,j+1))
                        !k
                        v(1) = -1/dz_if(i+1,k+1,j+1)
                        gnd_col(MatStencil_i,1) = i
                        gnd_col(MatStencil_j,1) = k
                        gnd_col(MatStencil_k,1) = j
                        !!k + 1
                        v(2) = 1/dz_if(i+1,k+1,j+1)
                        gnd_col(MatStencil_i,2) = i
                        gnd_col(MatStencil_j,2) = k+1
                        gnd_col(MatStencil_k,2) = j
                        !k + 2
                        v(3) = 0
                        gnd_col(MatStencil_i,3) = i
                        gnd_col(MatStencil_j,3) = k+2
                        gnd_col(MatStencil_k,3) = j

                        !k
                        !v(1) = -(2*sigma(i+1,k+1,j+1)+1)/(dz_if(i+1,k+1,j+1)*(sigma(i+1,k+1,j+1)+1))
                        !gnd_col(MatStencil_i,1) = i
                        !gnd_col(MatStencil_j,1) = k
                        !gnd_col(MatStencil_k,1) = j
                        !!k + 1
                        !v(2) = (sigma(i+1,k+1,j+1)+1)**2/(dz_if(i+1,k+1,j+1)*(sigma(i+1,k+1,j+1)+1))
                        !gnd_col(MatStencil_i,2) = i
                        !gnd_col(MatStencil_j,2) = k+1
                        !gnd_col(MatStencil_k,2) = j
                        !!k + 2
                        !v(3) = -(sigma(i+1,k+1,j+1)**2)/(dz_if(i+1,k+1,j+1)*(sigma(i+1,k+1,j+1)+1))
                        !gnd_col(MatStencil_i,3) = i
                        !gnd_col(MatStencil_j,3) = k+2
                        !gnd_col(MatStencil_k,3) = j
                        !
                        !!k + 1
                        !v(1) = -(sigma(i+1,k+1,j+1)**2-1)/(dz_if(i+1,k+2,j+1)*(sigma(i+1,k+1,j+1)+sigma(i+1,k+1,j+1)**2))
                        !gnd_col(MatStencil_i,1) = i
                        !gnd_col(MatStencil_j,1) = k+1
                        !gnd_col(MatStencil_k,1) = j
                        !!k
                        !v(2) = -1/(dz_if(i+1,k+2,j+1)*(sigma(i+1,k+1,j+1)+sigma(i+1,k+1,j+1)**2))
                        !gnd_col(MatStencil_i,2) = i
                        !gnd_col(MatStencil_j,2) = k
                        !gnd_col(MatStencil_k,2) = j
                        !!k + 2
                        !v(3) = sigma(i+1,k+1,j+1)**2/(dz_if(i+1,k+2,j+1)*(sigma(i+1,k+1,j+1)+sigma(i+1,k+1,j+1)**2))
                        !gnd_col(MatStencil_i,3) = i
                        !gnd_col(MatStencil_j,3) = k+2
                        !gnd_col(MatStencil_k,3) = j

                        
                        !Reminder: The equation for the lower BC has all of the lateral derivatives (dlambdadx, etc) NEGATIVE, hence opposite signs below
                        !If we are on left most border
                        if (i.eq.0) then
                            !i, k + 1
                            v(1) = v(1) + 3*dzdx_surf(i+1,j+1)/(denom*2*dx)
                            !i + 1, k + 1
                            v(4) = -4*dzdx_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,4) = i+1
                            gnd_col(MatStencil_j,4) = k
                            gnd_col(MatStencil_k,4) = j
                            !i + 2, k + 1
                            v(5) = dzdx_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,5) = i+2
                            gnd_col(MatStencil_j,5) = k
                            gnd_col(MatStencil_k,5) = j
                        !If we are on right most border
                        else if (i.eq.mx-1) then
                            !i, k + 1
                            v(1) = v(1) - 3*dzdx_surf(i+1,j+1)/(denom*2*dx)
                            !i - 1, k + 1
                            v(4) = 4*dzdx_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,4) = i-1
                            gnd_col(MatStencil_j,4) = k
                            gnd_col(MatStencil_k,4) = j
                            !i - 2, k + 1
                            v(5) = -dzdx_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,5) = i-2
                            gnd_col(MatStencil_j,5) = k
                            gnd_col(MatStencil_k,5) = j
                        else
                            !i - 1, k + 1
                            v(4) = dzdx_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,4) = i-1
                            gnd_col(MatStencil_j,4) = k
                            gnd_col(MatStencil_k,4) = j
                            !i + 1, k + 1
                            v(5) = -dzdx_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,5) = i+1
                            gnd_col(MatStencil_j,5) = k
                            gnd_col(MatStencil_k,5) = j
                        endif
                        
                        !!If we are on south most border
                        if (j.eq.0) then
                            !j, k + 1
                            v(1) = v(1) + 3*dzdy_surf(i+1,j+1)/(denom*2*dx)
                            !j + 1, k + 1
                            v(6) = -4*dzdy_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,6) = i
                            gnd_col(MatStencil_j,6) = k
                            gnd_col(MatStencil_k,6) = j+1
                            !j + 2, k + 1
                            v(7) = dzdy_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,7) = i
                            gnd_col(MatStencil_j,7) = k
                            gnd_col(MatStencil_k,7) = j+2
                        !If we are on north most border
                        else if (j.eq.my-1) then
                            !j, k + 1
                            v(1) = v(1) - 3*dzdy_surf(i+1,j+1)/(denom*2*dx)
                            !j + 1, k + 1
                            v(6) = 4*dzdy_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,6) = i
                            gnd_col(MatStencil_j,6) = k
                            gnd_col(MatStencil_k,6) = j-1
                            !j + 2, k + 1
                            v(7) = -dzdy_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,7) = i
                            gnd_col(MatStencil_j,7) = k
                            gnd_col(MatStencil_k,7) = j-2
                        else
                            !j - 1, k + 1
                            v(6) = dzdy_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,6) = i
                            gnd_col(MatStencil_j,6) = k
                            gnd_col(MatStencil_k,6) = j-1
                            !j + 1, k + 1
                            v(7) = -dzdy_surf(i+1,j+1)/(denom*2*dx)
                            gnd_col(MatStencil_i,7) = i
                            gnd_col(MatStencil_j,7) = k
                            gnd_col(MatStencil_k,7) = j+1
                        endif
                        call MatSetValuesStencil(arr_B,i1,row,i7,gnd_col,v,INSERT_VALUES, ierr)
                    endif
                ! interior points and boundary points with Neumann BC
                else if (   (k.ge.2 .and. k.le.mz-2 .and. &
                             i.ge.1 .and. i.le.mx-2 .and. &
                             j.ge.1 .and. j.le.my-2) .or. &
                             (k.eq.1)) then
                    v(:) = 0.0
                    
                    !Center
                    v(1) = A_coef(i+1,k,j+1)
                    col(MatStencil_i,1) = i
                    col(MatStencil_j,1) = k
                    col(MatStencil_k,1) = j
                    !k + 1
                    v(2) = B_coef(i+1,k,j+1)
                    col(MatStencil_i,2) = i
                    col(MatStencil_j,2) = k+1
                    col(MatStencil_k,2) = j
                    !k - 1
                    v(3) = C_coef(i+1,k,j+1)
                    col(MatStencil_i,3) = i
                    col(MatStencil_j,3) = k-1
                    col(MatStencil_k,3) = j
                    
                    if (j.eq.0) then
                        !j + 1, k - 1
                        v(8) = N_coef(i+1,k,j+1)
                        col(MatStencil_i,8) = i
                        col(MatStencil_j,8) = k-1
                        col(MatStencil_k,8) = j+1
                        !j + 1
                        v(6) = F_coef(i+1,k,j+1)
                        col(MatStencil_i,6) = i
                        col(MatStencil_j,6) = k
                        col(MatStencil_k,6) = j+1
                        !j + 1, k + 1
                        v(4) = L_coef(i+1,k,j+1)
                        col(MatStencil_i,4) = i
                        col(MatStencil_j,4) = k+1
                        col(MatStencil_k,4) = j+1
                        !j + 2, k - 1
                        v(9) = O_coef(i+1,k,j+1)
                        col(MatStencil_i,9) = i
                        col(MatStencil_j,9) = k-1
                        col(MatStencil_k,9) = j+2
                        !j + 2
                        v(7) = G_coef(i+1,k,j+1)
                        col(MatStencil_i,7) = i
                        col(MatStencil_j,7) = k
                        col(MatStencil_k,7) = j+2
                        !j + 2, k + 1
                        v(5) = M_coef(i+1,k,j+1)
                        col(MatStencil_i,5) = i
                        col(MatStencil_j,5) = k+1
                        col(MatStencil_k,5) = j+2
                    else if (j.eq.my-1) then
                        !j - 1, k - 1
                        v(9) = O_coef(i+1,k,j+1)
                        col(MatStencil_i,9) = i
                        col(MatStencil_j,9) = k-1
                        col(MatStencil_k,9) = j-1
                        !j - 1
                        v(7) = G_coef(i+1,k,j+1)
                        col(MatStencil_i,7) = i
                        col(MatStencil_j,7) = k
                        col(MatStencil_k,7) = j-1
                        !j - 1, k + 1
                        v(5) = M_coef(i+1,k,j+1)
                        col(MatStencil_i,5) = i
                        col(MatStencil_j,5) = k+1
                        col(MatStencil_k,5) = j-1
                        !j - 2, k - 1
                        v(8) = N_coef(i+1,k,j+1)
                        col(MatStencil_i,8) = i
                        col(MatStencil_j,8) = k-1
                        col(MatStencil_k,8) = j-2
                        !j - 2
                        v(6) = F_coef(i+1,k,j+1)
                        col(MatStencil_i,6) = i
                        col(MatStencil_j,6) = k
                        col(MatStencil_k,6) = j-2
                        !j - 2, k + 1
                        v(4) = L_coef(i+1,k,j+1)
                        col(MatStencil_i,4) = i
                        col(MatStencil_j,4) = k+1
                        col(MatStencil_k,4) = j-2
                    else
                        !j - 1, k - 1
                        v(9) = O_coef(i+1,k,j+1)
                        col(MatStencil_i,9) = i
                        col(MatStencil_j,9) = k-1
                        col(MatStencil_k,9) = j-1
                        !j - 1
                        v(7) = G_coef(i+1,k,j+1)
                        col(MatStencil_i,7) = i
                        col(MatStencil_j,7) = k
                        col(MatStencil_k,7) = j-1
                        !j - 1, k + 1
                        v(5) = M_coef(i+1,k,j+1)
                        col(MatStencil_i,5) = i
                        col(MatStencil_j,5) = k+1
                        col(MatStencil_k,5) = j-1
                        !j + 1, k + 1
                        v(4) = L_coef(i+1,k,j+1)
                        col(MatStencil_i,4) = i
                        col(MatStencil_j,4) = k+1
                        col(MatStencil_k,4) = j+1
                        !j + 1
                        v(6) = F_coef(i+1,k,j+1)
                        col(MatStencil_i,6) = i
                        col(MatStencil_j,6) = k
                        col(MatStencil_k,6) = j+1
                        !j + 1, k - 1
                        v(8) = N_coef(i+1,k,j+1)
                        col(MatStencil_i,8) = i
                        col(MatStencil_j,8) = k-1
                        col(MatStencil_k,8) = j+1
                    endif
                    
                    if (i.eq.0) then
                        !i + 1, k - 1
                        v(14) = J_coef(i+1,k,j+1)
                        col(MatStencil_i,14) = i+1
                        col(MatStencil_j,14) = k-1
                        col(MatStencil_k,14) = j
                        !i + 1
                        v(12) = D_coef(i+1,k,j+1)
                        col(MatStencil_i,12) = i+1
                        col(MatStencil_j,12) = k
                        col(MatStencil_k,12) = j
                        !i + 1, k + 1
                        v(10) = H_coef(i+1,k,j+1)
                        col(MatStencil_i,10) = i+1
                        col(MatStencil_j,10) = k+1
                        col(MatStencil_k,10) = j
                        !i + 2, k - 1
                        v(15) = K_coef(i+1,k,j+1)
                        col(MatStencil_i,15) = i+2
                        col(MatStencil_j,15) = k-1
                        col(MatStencil_k,15) = j
                        !i + 2
                        v(13) = E_coef(i+1,k,j+1)
                        col(MatStencil_i,13) = i+2
                        col(MatStencil_j,13) = k
                        col(MatStencil_k,13) = j
                        !i + 2, k + 1
                        v(11) = I_coef(i+1,k,j+1)
                        col(MatStencil_i,11) = i+2
                        col(MatStencil_j,11) = k+1
                        col(MatStencil_k,11) = j
                    else if (i.eq.mx-1) then
                        !i - 1, k - 1
                        v(15) = K_coef(i+1,k,j+1)
                        col(MatStencil_i,15) = i-1
                        col(MatStencil_j,15) = k-1
                        col(MatStencil_k,15) = j
                        !i - 1
                        v(13) = E_coef(i+1,k,j+1)
                        col(MatStencil_i,13) = i-1
                        col(MatStencil_j,13) = k
                        col(MatStencil_k,13) = j
                        !i - 1, k + 1
                        v(11) = I_coef(i+1,k,j+1)
                        col(MatStencil_i,11) = i-1
                        col(MatStencil_j,11) = k+1
                        col(MatStencil_k,11) = j
                        !i - 2, k - 1
                        v(14) = J_coef(i+1,k,j+1)
                        col(MatStencil_i,14) = i-2
                        col(MatStencil_j,14) = k-1
                        col(MatStencil_k,14) = j
                        !i - 2
                        v(12) = D_coef(i+1,k,j+1)
                        col(MatStencil_i,12) = i-2
                        col(MatStencil_j,12) = k
                        col(MatStencil_k,12) = j
                        !i - 2, k + 1
                        v(10) = H_coef(i+1,k,j+1)
                        col(MatStencil_i,10) = i-2
                        col(MatStencil_j,10) = k+1
                        col(MatStencil_k,10) = j
                    else 
                        !i - 1, k - 1
                        v(15) = K_coef(i+1,k,j+1)
                        col(MatStencil_i,15) = i-1
                        col(MatStencil_j,15) = k-1
                        col(MatStencil_k,15) = j
                        !i - 1
                        v(13) = E_coef(i+1,k,j+1)
                        col(MatStencil_i,13) = i-1
                        col(MatStencil_j,13) = k
                        col(MatStencil_k,13) = j
                        !i - 1, k + 1
                        v(11) = I_coef(i+1,k,j+1)
                        col(MatStencil_i,11) = i-1
                        col(MatStencil_j,11) = k+1
                        col(MatStencil_k,11) = j
                        !i + 1, k + 1
                        v(10) = H_coef(i+1,k,j+1)
                        col(MatStencil_i,10) = i+1
                        col(MatStencil_j,10) = k+1
                        col(MatStencil_k,10) = j
                        !i + 1
                        v(12) = D_coef(i+1,k,j+1)
                        col(MatStencil_i,12) = i+1
                        col(MatStencil_j,12) = k
                        col(MatStencil_k,12) = j
                        !i + 1, k - 1
                        v(14) = J_coef(i+1,k,j+1)
                        col(MatStencil_i,14) = i+1
                        col(MatStencil_j,14) = k-1
                        col(MatStencil_k,14) = j
                    endif
                    call MatSetValuesStencil(arr_B,i1,row,i15,col,v,INSERT_VALUES, ierr)
                ! Dirlect BC for lateral boundaries
                else
                    v(1) = 1.0
                    call MatSetValuesStencil(arr_B,i1,row,i1,row,v,INSERT_VALUES, ierr)
                endif
            enddo
        enddo
        enddo

        call MatAssemblyBegin(arr_B,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(arr_B,MAT_FINAL_ASSEMBLY,ierr)
        !if (arr_A .ne. arr_B) then
        ! call MatAssemblyBegin(arr_A,MAT_FINAL_ASSEMBLY,ierr)
        ! call MatAssemblyEnd(arr_A,MAT_FINAL_ASSEMBLY,ierr)
        !endif

    end subroutine ComputeMatrix
    
    subroutine initialize_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        
        real, allocatable, dimension(:,:,:) :: mixed_denom
        
        allocate(A_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(B_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(C_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(D_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(E_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(F_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(G_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(H_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(I_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(J_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(K_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(L_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(M_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(N_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(O_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        
        allocate(mixed_denom(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(tmp_x(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(tmp_y(i_s:i_e,k_s:k_e,j_s:j_e))

        A_coef = 1 !Because this corresponds to the centered node, must be non-zero, otherwise there is no solution
        B_coef = 0
        C_coef = 0
        D_coef = 0
        E_coef = 0
        F_coef = 0
        G_coef = 0
        H_coef = 0
        I_coef = 0
        J_coef = 0
        K_coef = 0
        L_coef = 0
        M_coef = 0
        N_coef = 0
        O_coef = 0

        mixed_denom = 2*dx*dz_if(:,k_s+1:k_e+1,:)*(sigma+sigma**2)


        D_coef = jaco/(dx**2) + (sigma**2 - 1)*(2*dzdx)/mixed_denom
        F_coef = jaco/(dx**2) + (sigma**2 - 1)*(2*dzdy)/mixed_denom

        E_coef = jaco/(dx**2) - (sigma**2 - 1)*(2*dzdx)/mixed_denom
        G_coef = jaco/(dx**2) - (sigma**2 - 1)*(2*dzdy)/mixed_denom

        H_coef = -(sigma**2)*(2*dzdx)/mixed_denom
                
        J_coef = (2*dzdx)/mixed_denom

        L_coef = -(sigma**2)*(2*dzdy)/mixed_denom
                
        N_coef = (2*dzdy)/mixed_denom

        !The following coeffs are just the opposites...
        I_coef = -H_coef
        K_coef = -J_coef
        M_coef = -L_coef
        O_coef = -N_coef
        
        tmp_x = J_coef
        tmp_y = N_coef

        if (domain%grid%ims==domain%grid%ids) then
            !become the "i+2" indices
            I_coef(i_s,:,j_s:j_e) = tmp_x(i_s,:,j_s:j_e)*(sigma(i_s,:,j_s:j_e)**2)
            E_coef(i_s,:,j_s:j_e) = jaco(i_s,:,j_s:j_e)/(dx**2) - tmp_x(i_s,:,j_s:j_e)*(sigma(i_s,:,j_s:j_e)**2-1)
            K_coef(i_s,:,j_s:j_e) = -tmp_x(i_s,:,j_s:j_e)
            !remain the "i+1" indices, but require modification
            H_coef(i_s,:,j_s:j_e) = -4*tmp_x(i_s,:,j_s:j_e)*(sigma(i_s,:,j_s:j_e)**2)
            D_coef(i_s,:,j_s:j_e) = -2*jaco(i_s,:,j_s:j_e)/(dx**2) + 4*tmp_x(i_s,:,j_s:j_e)*(sigma(i_s,:,j_s:j_e)**2-1)
            J_coef(i_s,:,j_s:j_e) = 4*tmp_x(i_s,:,j_s:j_e)
        endif
        
        if (domain%grid%ime==domain%grid%ide) then
            !become the "i-2" indices
            H_coef(i_e,:,j_s:j_e) = -tmp_x(i_e,:,j_s:j_e)*(sigma(i_e,:,j_s:j_e)**2)
            D_coef(i_e,:,j_s:j_e) = jaco(i_e,:,j_s:j_e)/(dx**2) + tmp_x(i_e,:,j_s:j_e)*(sigma(i_e,:,j_s:j_e)**2-1)
            J_coef(i_e,:,j_s:j_e) = tmp_x(i_e,:,j_s:j_e)
            !remain the "i-1" indices, but require modification
            I_coef(i_e,:,j_s:j_e) = 4*tmp_x(i_e,:,j_s:j_e)*(sigma(i_e,:,j_s:j_e)**2)
            E_coef(i_e,:,j_s:j_e) = -2*jaco(i_e,:,j_s:j_e)/(dx**2) - 4*tmp_x(i_e,:,j_s:j_e)*(sigma(i_e,:,j_s:j_e)**2-1)
            K_coef(i_e,:,j_s:j_e) = -4*tmp_x(i_e,:,j_s:j_e)
        endif

        if (domain%grid%jms==domain%grid%jds) then
            !become the "j+2" indices
            M_coef(i_s:i_e,:,j_s) = tmp_y(i_s:i_e,:,j_s)*(sigma(i_s:i_e,:,j_s)**2)
            G_coef(i_s:i_e,:,j_s) = jaco(i_s:i_e,:,j_s)/(dx**2) - tmp_y(i_s:i_e,:,j_s)*(sigma(i_s:i_e,:,j_s)**2-1)
            O_coef(i_s:i_e,:,j_s) = -tmp_y(i_s:i_e,:,j_s)
            !remain the "j+1" indices, but require modification
            L_coef(i_s:i_e,:,j_s) = -4*tmp_y(i_s:i_e,:,j_s)*(sigma(i_s:i_e,:,j_s)**2)
            F_coef(i_s:i_e,:,j_s) = -2*jaco(i_s:i_e,:,j_s)/(dx**2) + 4*tmp_y(i_s:i_e,:,j_s)*(sigma(i_s:i_e,:,j_s)**2-1)
            N_coef(i_s:i_e,:,j_s) = 4*tmp_y(i_s:i_e,:,j_s)
        endif
        
        if (domain%grid%jme==domain%grid%jde) then
            !become the "j-2" indices
            L_coef(i_s:i_e,:,j_e) = -tmp_y(i_s:i_e,:,j_e)*(sigma(i_s:i_e,:,j_e)**2)
            F_coef(i_s:i_e,:,j_e) = jaco(i_s:i_e,:,j_e)/(dx**2) + tmp_y(i_s:i_e,:,j_e)*(sigma(i_s:i_e,:,j_e)**2-1)
            N_coef(i_s:i_e,:,j_e) = tmp_y(i_s:i_e,:,j_e)
            !remain the "j-1" indices, but require modification
            M_coef(i_s:i_e,:,j_e) = 4*tmp_y(i_s:i_e,:,j_e)*(sigma(i_s:i_e,:,j_e)**2)
            G_coef(i_s:i_e,:,j_e) = -2*jaco(i_s:i_e,:,j_e)/(dx**2) - 4*tmp_y(i_s:i_e,:,j_e)*(sigma(i_s:i_e,:,j_e)**2-1)
            O_coef(i_s:i_e,:,j_e) = -4*tmp_y(i_s:i_e,:,j_e)
        endif

        !This function sets A, B, annd C coefficients
        call update_coefs(domain)

    end subroutine initialize_coefs
    
    
    !Update the coefs which change with time, i.e. those which depend on alpha
    subroutine update_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain        
        real, allocatable, dimension(:,:,:) :: mixed_denom, X_coef, Y_coef, dalphadz
        
        allocate(mixed_denom(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(X_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(Y_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(dalphadz(i_s:i_e,k_s:k_e,j_s:j_e))

        dalphadz(:,k_s+1:k_e-1,:) = (sigma(:,k_s+1:k_e-1,:)**2*(alpha(i_s:i_e,k_s+2:k_e,j_s:j_e)**2) - &
                (sigma(:,k_s+1:k_e-1,:)**2-1)*(alpha(i_s:i_e,k_s+1:k_e-1,j_s:j_e)**2)-(alpha(i_s:i_e,k_s:k_e-2,j_s:j_e)**2)) / &
                (dz_if(:,k_s+2:k_e,:)*(sigma(:,k_s+1:k_e-1,:)+sigma(:,k_s+1:k_e-1,:)**2))

        dalphadz(:,k_s,:) = -(sigma(:,k_s+1,:)**2*(alpha(i_s:i_e,k_s+2,j_s:j_e)**2))/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1)) + &
                          ((sigma(:,k_s+1,:)+1)*(alpha(i_s:i_e,k_s+1,j_s:j_e)**2))/dz_if(:,k_s+1,:) - &
                        ((2*sigma(:,k_s+1,:)+1)*(alpha(i_s:i_e,k_s,j_s:j_e)**2))/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1))

        dalphadz(:,k_e,:) = (alpha(i_s:i_e,k_e-2,j_s:j_e)**2)/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)*(1+sigma(:,k_e-1,:))) - &
                        ((1+sigma(:,k_e-1,:))*(alpha(i_s:i_e,k_e-1,j_s:j_e)**2))/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)) + &
                        ((2+sigma(:,k_e-1,:))*(alpha(i_s:i_e,k_e,j_s:j_e)**2))/(dz_if(:,k_e,:)*(1+sigma(:,k_e-1,:)))


        mixed_denom = dz_if(:,k_s+1:k_e+1,:)*(sigma+sigma**2)

        Y_coef = 2*(alpha**2 + dzdy**2 + dzdx**2) / (jaco*(dz_if(:,k_s+1:k_e+1,:)**2)*(sigma+sigma**2))
        
        X_coef = -d2zdx2 - d2zdy2 - (dzdx*dzhatdxdz+dzdy*dzhatdydz)*jaco + (dzdx*dzdxz+dzdy*dzdyz)/jaco + (dzdx**2 + dzdy**2 + alpha**2)*dzhatdzz
        
        B_coef = sigma * Y_coef + sigma**2*X_coef/mixed_denom 
        C_coef = Y_coef - X_coef/mixed_denom
        
        !B_coef = B_coef + sigma**2*X_coef/mixed_denom 
        !C_coef = C_coef - X_coef/mixed_denom
        
        A_coef = 0.0
        A_coef(i_s+1:i_e-1,:,j_s+1:j_e-1) = -(4*jaco(i_s+1:i_e-1,:,j_s+1:j_e-1)/(dx**2))
        !A_coef(i_s+1:i_e-1,:,j_s+1:j_e-1) = -((domain%jacobian(i_s+2:i_e,k_s:k_e,j_s+1:j_e-1) + &
        !                                       2*domain%jacobian(i_s+1:i_e-1,k_s:k_e,j_s+1:j_e-1) + &
        !                                       domain%jacobian(i_s:i_e-2,k_s:k_e,j_s+1:j_e-1))/(2*domain%dx**2)) &
        !                                    -((domain%jacobian(i_s+1:i_e-1,k_s:k_e,j_s+2:j_e) + &
        !                                       2*domain%jacobian(i_s+1:i_e-1,k_s:k_e,j_s+1:j_e-1) + &
        !                                       domain%jacobian(i_s+1:i_e-1,k_s:k_e,j_s:j_e-2))/(2*domain%dx**2))       

        if (domain%grid%ims==domain%grid%ids) then
            A_coef(i_s,:,j_s+1:j_e-1) = A_coef(i_s,:,j_s+1:j_e-1) - (  jaco(i_s,:,j_s+1:j_e-1)/(dx**2))
            B_coef(i_s,:,j_s+1:j_e-1) = B_coef(i_s,:,j_s+1:j_e-1) + 3*tmp_x(i_s,:,j_s+1:j_e-1)*sigma(i_s,:,j_s+1:j_e-1)**2
            C_coef(i_s,:,j_s+1:j_e-1) = C_coef(i_s,:,j_s+1:j_e-1) - 3*tmp_x(i_s,:,j_s+1:j_e-1)
        else
            A_coef(i_s,:,j_s+1:j_e-1) = A_coef(i_s,:,j_s+1:j_e-1) - (4*jaco(i_s,:,j_s+1:j_e-1)/(dx**2))
        endif
        
        if (domain%grid%ime==domain%grid%ide) then
            A_coef(i_e,:,j_s+1:j_e-1) = A_coef(i_e,:,j_s+1:j_e-1) - (  jaco(i_e,:,j_s+1:j_e-1)/(dx**2))
            B_coef(i_e,:,j_s+1:j_e-1) = B_coef(i_e,:,j_s+1:j_e-1) - 3*tmp_x(i_e,:,j_s+1:j_e-1)*sigma(i_e,:,j_s+1:j_e-1)**2
            C_coef(i_e,:,j_s+1:j_e-1) = C_coef(i_e,:,j_s+1:j_e-1) + 3*tmp_x(i_e,:,j_s+1:j_e-1)
        else
            A_coef(i_e,:,j_s+1:j_e-1) = A_coef(i_e,:,j_s+1:j_e-1)-(4*jaco(i_e,:,j_s+1:j_e-1)/(dx**2))
        endif

        if (domain%grid%jms==domain%grid%jds) then
            A_coef(i_s+1:i_e-1,:,j_s) = A_coef(i_s+1:i_e-1,:,j_s) - (  jaco(i_s+1:i_e-1,:,j_s)/(dx**2))
            B_coef(i_s+1:i_e-1,:,j_s) = B_coef(i_s+1:i_e-1,:,j_s) + 3*tmp_y(i_s+1:i_e-1,:,j_s)*sigma(i_s+1:i_e-1,:,j_s)**2
            C_coef(i_s+1:i_e-1,:,j_s) = C_coef(i_s+1:i_e-1,:,j_s) - 3*tmp_y(i_s+1:i_e-1,:,j_s)
        else
            A_coef(i_s+1:i_e-1,:,j_s) = A_coef(i_s+1:i_e-1,:,j_s)-(4*jaco(i_s+1:i_e-1,:,j_s)/(dx**2))
        endif
        
        if (domain%grid%jme==domain%grid%jde) then
            A_coef(i_s+1:i_e-1,:,j_e) = A_coef(i_s+1:i_e-1,:,j_e) - (  jaco(i_s+1:i_e-1,:,j_e)/(dx**2))
            B_coef(i_s+1:i_e-1,:,j_e) = B_coef(i_s+1:i_e-1,:,j_e) - 3*tmp_y(i_s+1:i_e-1,:,j_e)*sigma(i_s+1:i_e-1,:,j_e)**2
            C_coef(i_s+1:i_e-1,:,j_e) = C_coef(i_s+1:i_e-1,:,j_e) + 3*tmp_y(i_s+1:i_e-1,:,j_e)
        else
            A_coef(i_s+1:i_e-1,:,j_e) = A_coef(i_s+1:i_e-1,:,j_e)-(4*jaco(i_s+1:i_e-1,:,j_e)/(dx**2))
        endif
        
        
        !North-west corner
        if (domain%grid%ims==domain%grid%ids .and. domain%grid%jme==domain%grid%jde) then
            A_coef(i_s,:,j_e) = A_coef(i_s,:,j_e) + 2*(  jaco(i_s,:,j_e)/(dx**2))
            B_coef(i_s,:,j_e) = B_coef(i_s,:,j_e) + 3*tmp_x(i_s,:,j_e)*sigma(i_s,:,j_e)**2 - 3*tmp_y(i_s,:,j_e)*sigma(i_s,:,j_e)**2
            C_coef(i_s,:,j_e) = C_coef(i_s,:,j_e) - 3*tmp_x(i_s,:,j_e) + 3*tmp_y(i_s,:,j_e)
        else
            A_coef(i_s,:,j_e) = A_coef(i_s,:,j_e)-(4*jaco(i_s,:,j_e)/(dx**2))
        endif
        !North-east corner
        if (domain%grid%ime==domain%grid%ide .and. domain%grid%jme==domain%grid%jde) then
            A_coef(i_e,:,j_e) = A_coef(i_e,:,j_e) + 2*(  jaco(i_e,:,j_e)/(dx**2))
            B_coef(i_e,:,j_e) = B_coef(i_e,:,j_e) - 3*tmp_x(i_e,:,j_e)*sigma(i_e,:,j_e)**2 - 3*tmp_y(i_e,:,j_e)*sigma(i_e,:,j_e)**2
            C_coef(i_e,:,j_e) = C_coef(i_e,:,j_e) + 3*tmp_x(i_e,:,j_e) + 3*tmp_y(i_e,:,j_e)
        else
            A_coef(i_e,:,j_e) = A_coef(i_e,:,j_e)-(4*jaco(i_e,:,j_e)/(dx**2))
        endif
        !South-east corner
        if (domain%grid%ime==domain%grid%ide .and. domain%grid%jms==domain%grid%jds) then
            A_coef(i_e,:,j_s) = A_coef(i_e,:,j_s) + 2*(  jaco(i_e,:,j_s)/(dx**2))
            B_coef(i_e,:,j_s) = B_coef(i_e,:,j_s) - 3*tmp_x(i_e,:,j_s)*sigma(i_e,:,j_s)**2 + 3*tmp_y(i_e,:,j_s)*sigma(i_e,:,j_s)**2
            C_coef(i_e,:,j_s) = C_coef(i_e,:,j_s) + 3*tmp_x(i_e,:,j_s) - 3*tmp_y(i_e,:,j_s)
        else
            A_coef(i_e,:,j_s) = A_coef(i_e,:,j_s)-(4*jaco(i_e,:,j_s)/(dx**2))
        endif
        !South-west corner
        if (domain%grid%ims==domain%grid%ids .and. domain%grid%jms==domain%grid%jds) then
            A_coef(i_s,:,j_s) = A_coef(i_s,:,j_s) + 2*(  jaco(i_s,:,j_s)/(dx**2))
            B_coef(i_s,:,j_s) = B_coef(i_s,:,j_s) + 3*tmp_x(i_s,:,j_s)*sigma(i_s,:,j_s)**2 + 3*tmp_y(i_s,:,j_s)*sigma(i_s,:,j_s)**2
            C_coef(i_s,:,j_s) = C_coef(i_s,:,j_s) - 3*tmp_x(i_s,:,j_s) - 3*tmp_y(i_s,:,j_s)
        else
            A_coef(i_s,:,j_s) = A_coef(i_s,:,j_s)-(4*jaco(i_s,:,j_s)/(dx**2))
        endif
        
        A_coef = A_coef - B_coef - C_coef


    end subroutine


    subroutine finalize_iter_winds()
        implicit none

        PetscErrorCode ierr

        call PetscFinalize(ierr)
    end subroutine

    subroutine init_iter_winds(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        PetscErrorCode ierr

        
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        if (ierr .ne. 0) then
            print*,'Unable to initialize PETSc'
            stop
        endif 
        call init_module_vars(domain)
        if(this_image()==1) write(*,*) 'Initialized PETSc'
    end subroutine
    
    subroutine init_module_vars(domain)
        implicit none
        type(domain_t), intent(in) :: domain

        i_s = domain%its-1
        i_e = domain%ite+1
        k_s = domain%kts  
        k_e = domain%kte  
        j_s = domain%jts-1
        j_e = domain%jte+1
        
        ims = domain%ims
        ime = domain%ime
        kms = domain%kms  
        kme = domain%kme  
        jms = domain%jms
        jme = domain%jme

        
        !i_s+hs, unless we are on global boundary, then i_s
        if (domain%grid%ims==domain%grid%ids) i_s = domain%grid%ids
        
        !i_e, unless we are on global boundary, then i_e+1
        if (domain%grid%ime==domain%grid%ide) i_e = domain%grid%ide
        
        !j_s+hs, unless we are on global boundary, then j_s
        if (domain%grid%jms==domain%grid%jds) j_s = domain%grid%jds
        
        !j_e, unless we are on global boundary, then j_e+1
        if (domain%grid%jme==domain%grid%jde) j_e = domain%grid%jde

        hs = domain%grid%halo_size
        if (.not.(allocated(dzdx))) then
            !call MPI_INIT(ierr)
            allocate(dzdx(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(dzdy(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(jaco(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(dzhatdxdz(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(dzhatdydz(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(dzhatdzz(i_s:i_e,k_s:k_e,j_s:j_e))
            
            allocate(dzdxz(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(dzdyz(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(w_0(i_s:i_e,j_s:j_e))
            allocate(dzdx_surf(i_s:i_e,j_s:j_e))
            allocate(dzdy_surf(i_s:i_e,j_s:j_e))
            
            allocate(d2zdx2(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(d2zdy2(i_s:i_e,k_s:k_e,j_s:j_e))

            allocate(sigma(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(dz_if(i_s:i_e,k_s:k_e+1,j_s:j_e))
            allocate(alpha(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(div(ims:ime,kms:kme,jms:jme))
            allocate( xl( 1:domain%grid%ximages ))
            allocate( yl( 1:domain%grid%yimages ))
            xl = 0
            yl = 0
            
            dx = domain%dx
            dzdx  = domain%dzdx(i_s:i_e,k_s:k_e,j_s:j_e) 
            dzdyz = domain%dzdyz(i_s:i_e,k_s:k_e,j_s:j_e)
            dzdxz = domain%dzdxz(i_s:i_e,k_s:k_e,j_s:j_e) 
            dzdy  = domain%dzdy(i_s:i_e,k_s:k_e,j_s:j_e)
            jaco = domain%jacobian(i_s:i_e,k_s:k_e,j_s:j_e)
            
            
            !call smooth_array_3d( dzdx, windowsize = 2, ydim = 3)
            !call smooth_array_3d( dzdy, windowsize = 2, ydim = 3)
            !div = domain%jacobian !(i_s:i_e,k_s:k_e,j_s:j_e)
            !call smooth_array_3d( div, windowsize = 1, ydim = 3)
            !jaco = div(i_s:i_e,k_s:k_e,j_s:j_e)
            
            dz_if(:,k_s+1:k_e,:) = (domain%advection_dz(i_s:i_e,k_s+1:k_e,j_s:j_e) + &
                                   domain%advection_dz(i_s:i_e,k_s:k_e-1,j_s:j_e))/2
            dz_if(:,k_s,:)   = domain%advection_dz(i_s:i_e,k_s,j_s:j_e)
            dz_if(:,k_e+1,:) = domain%advection_dz(i_s:i_e,k_e,j_s:j_e)
            sigma = dz_if(:,k_s:k_e,:)/dz_if(:,k_s+1:k_e+1,:)
                          
            dzdx_surf = 0.1
            dzdy_surf = 0.1
            dzdx_surf(i_s+1:i_e-1,j_s:j_e) = (domain%neighbor_terrain(i_s+2:i_e,j_s:j_e)-domain%neighbor_terrain(i_s:i_e-2,j_s:j_e))/(2*dx)
            dzdy_surf(i_s:i_e,j_s+1:j_e-1) = (domain%neighbor_terrain(i_s:i_e,j_s+2:j_e)-domain%neighbor_terrain(i_s:i_e,j_s:j_e-2))/(2*dx)
                          
            dzdx_surf(i_s,j_s:j_e) = dzdx_surf(i_s+1,j_s:j_e)
            dzdx_surf(i_e,j_s:j_e) = dzdx_surf(i_e-1,j_s:j_e)
            dzdy_surf(i_s:i_e,j_s) = dzdy_surf(i_s:i_e,j_s+1)
            dzdy_surf(i_s:i_e,j_e) = dzdy_surf(i_s:i_e,j_e-1)

            dzhatdxdz(i_s:i_e,k_s:k_e,:) = &
                    ((1/domain%jacobian_u(i_s+1:i_e+1,k_s:k_e,j_s:j_e))-(1/domain%jacobian_u(i_s:i_e,k_s:k_e,j_s:j_e)))/(dx)
            dzhatdydz(:,k_s:k_e,j_s:j_e) = &
                    ((1/domain%jacobian_v(i_s:i_e,k_s:k_e,j_s+1:j_e+1))-(1/domain%jacobian_v(i_s:i_e,k_s:k_e,j_s:j_e)))/(dx)

            dzhatdxdz(i_s+1:i_e-1,k_s:k_e,:) = &
                    ((1/domain%jacobian(i_s+2:i_e,k_s:k_e,j_s:j_e))-(1/domain%jacobian(i_s:i_e-2,k_s:k_e,j_s:j_e)))/(2*dx)
            dzhatdydz(:,k_s:k_e,j_s+1:j_e-1) = &
                    ((1/domain%jacobian(i_s:i_e,k_s:k_e,j_s+2:j_e))-(1/domain%jacobian(i_s:i_e,k_s:k_e,j_s:j_e-2)))/(2*dx)


            d2zdx2(i_s+1:i_e-1,k_s:k_e,:) = (domain%z%data_3d(i_s+2:i_e,k_s:k_e,j_s:j_e) - &
                                                2*domain%z%data_3d(i_s+1:i_e-1,k_s:k_e,j_s:j_e) + &
                                                domain%z%data_3d(i_s:i_e-2,k_s:k_e,j_s:j_e))/(dx**2)
                                                
            d2zdy2(:,k_s:k_e,j_s+1:j_e-1) = (domain%z%data_3d(i_s:i_e,k_s:k_e,j_s+2:j_e) - &
                                                2*domain%z%data_3d(i_s:i_e,k_s:k_e,j_s+1:j_e-1) + &
                                                domain%z%data_3d(i_s:i_e,k_s:k_e,j_s:j_e-2))/(dx**2)
                                                
            d2zdx2(i_s:i_e,k_s:k_e,:) = (domain%dzdx_u(i_s+1:i_e+1,k_s:k_e,j_s:j_e) - &
                                                domain%dzdx_u(i_s:i_e,k_s:k_e,j_s:j_e))/(dx)
                                                
            d2zdy2(:,k_s:k_e,j_s:j_e) = (domain%dzdy_v(i_s:i_e,k_s:k_e,j_s+1:j_e+1) - &
                                                domain%dzdy_v(i_s:i_e,k_s:k_e,j_s:j_e))/(dx)

            !d2zdx2(i_s+2:i_e-2,k_s:k_e,:) = (-domain%z%data_3d(i_s+4:i_e,  k_s:k_e,j_s:j_e) &
            !                              +16*domain%z%data_3d(i_s+3:i_e-1,k_s:k_e,j_s:j_e) &
            !                              -30*domain%z%data_3d(i_s+2:i_e-2,k_s:k_e,j_s:j_e) &
            !                              +16*domain%z%data_3d(i_s+1:i_e-3,k_s:k_e,j_s:j_e) &
            !                                 -domain%z%data_3d(i_s  :i_e-4,k_s:k_e,j_s:j_e))/(12*dx**2)
            !d2zdy2(:,k_s:k_e,j_s+2:j_e-2) = (-domain%z%data_3d(i_s:i_e,k_s:k_e,j_s+4:j_e  ) &
            !                              +16*domain%z%data_3d(i_s:i_e,k_s:k_e,j_s+3:j_e-1) &
            !                              -30*domain%z%data_3d(i_s:i_e,k_s:k_e,j_s+2:j_e-2) &
            !                              +16*domain%z%data_3d(i_s:i_e,k_s:k_e,j_s+1:j_e-3) &
            !                                 -domain%z%data_3d(i_s:i_e,k_s:k_e,j_s  :j_e-4))/(12*dx**2)

            if (domain%grid%ims==domain%grid%ids) then
                dzhatdxdz(i_s,:,:) = (-(1/jaco(i_s+2,:,:)) + 4*(1/jaco(i_s+1,:,:)) - 3*(1/jaco(i_s,:,:)) )/(2*dx)
                !d2zdx2(i_s,:,:) = d2zdx2(i_s+1,:,:)
            else
                !dzhatdxdz(i_s,:,:) = ((1/domain%jacobian(i_s+1,k_s:k_e,j_s:j_e)) - &
                !                 (1/domain%jacobian(i_s-1,k_s:k_e,j_s:j_e)))/(2*dx)
                !d2zdx2(i_s,:,:) = (domain%z%data_3d(i_s+1,k_s:k_e,j_s:j_e) - &
                !                    2*domain%z%data_3d(i_s,k_s:k_e,j_s:j_e) + &
                !                    domain%z%data_3d(i_s-1,k_s:k_e,j_s:j_e))/(dx**2)
            endif
            
            if (domain%grid%ime==domain%grid%ide) then
                dzhatdxdz(i_e,:,:) = ((1/jaco(i_e-2,:,:)) - 4*(1/jaco(i_e-1,:,:)) + 3*(1/jaco(i_e,:,:)) )/(2*dx)
                !d2zdx2(i_e,:,:) = d2zdx2(i_e-1,:,:)
            else
                !dzhatdxdz(i_e,:,:) = ((1/domain%jacobian(i_e+1,k_s:k_e,j_s:j_e)) - &
                !                 (1/domain%jacobian(i_e-1,k_s:k_e,j_s:j_e)))/(2*dx)
                !d2zdx2(i_e,:,:) = (domain%z%data_3d(i_e+1,k_s:k_e,j_s:j_e) - &
                !                    2*domain%z%data_3d(i_e,k_s:k_e,j_s:j_e) + &
                !                    domain%z%data_3d(i_e-1,k_s:k_e,j_s:j_e))/(dx**2)
            endif
            
            if (domain%grid%jms==domain%grid%jds) then
                dzhatdydz(:,:,j_s) = (-(1/jaco(:,:,j_s+2)) + 4*(1/jaco(:,:,j_s+1)) - 3*(1/jaco(:,:,j_s)) )/(2*dx)
                !d2zdy2(:,:,j_s) = d2zdy2(:,:,j_s+1)
            else
                !dzhatdydz(:,:,j_s) = ((1/domain%jacobian(i_s:i_e,k_s:k_e,j_s+1)) - &
                !                 (1/domain%jacobian(i_s:i_e,k_s:k_e,j_s-1)))/(2*dx)
                !d2zdy2(:,:,j_s) = (domain%z%data_3d(i_s:i_e,k_s:k_e,j_s+1) - &
                !                    2*domain%z%data_3d(i_s:i_e,k_s:k_e,j_s) + &
                !                    domain%z%data_3d(i_s:i_e,k_s:k_e,j_s-1))/(dx**2)
            endif
            
            if (domain%grid%jme==domain%grid%jde) then
                dzhatdydz(:,:,j_e) = ((1/jaco(:,:,j_e-2)) - 4*(1/jaco(:,:,j_e-1)) + 3*(1/jaco(:,:,j_e)) )/(2*dx)
                !d2zdy2(:,:,j_e) = d2zdy2(:,:,j_e-1)
            else
                !dzhatdydz(:,:,j_e) = ((1/domain%jacobian(i_s:i_e,k_s:k_e,j_e+1)) - &
                !                 (1/domain%jacobian(i_s:i_e,k_s:k_e,j_e-1)))/(2*dx)
                !d2zdy2(:,:,j_e) = (domain%z%data_3d(i_s:i_e,k_s:k_e,j_e+1) - &
                !                    2*domain%z%data_3d(i_s:i_e,k_s:k_e,j_e) + &
                !                    domain%z%data_3d(i_s:i_e,k_s:k_e,j_e-1))/(dx**2)
            endif


            dzhatdzz(:,k_s+1:k_e-1,:) = (sigma(:,k_s+1:k_e-1,:)**2*(1/jaco(i_s:i_e,k_s+2:k_e,j_s:j_e)) - &
                    (sigma(:,k_s+1:k_e-1,:)**2-1)*(1/jaco(i_s:i_e,k_s+1:k_e-1,j_s:j_e))-(1/jaco(i_s:i_e,k_s:k_e-2,j_s:j_e))) / &
                    (dz_if(:,k_s+2:k_e,:)*(sigma(:,k_s+1:k_e-1,:)+sigma(:,k_s+1:k_e-1,:)**2))
                    
            dzhatdzz(:,k_s,:) = -(sigma(:,k_s+1,:)**2*(1/jaco(i_s:i_e,k_s+2,j_s:j_e)))/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1)) + &
                              ((sigma(:,k_s+1,:)+1)*(1/jaco(i_s:i_e,k_s+1,j_s:j_e)))/dz_if(:,k_s+1,:) - &
                            ((2*sigma(:,k_s+1,:)+1)*(1/jaco(i_s:i_e,k_s,j_s:j_e)))/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1))
                            
            dzhatdzz(:,k_e,:) = (1/jaco(i_s:i_e,k_e-2,j_s:j_e))/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)*(1+sigma(:,k_e-1,:))) - &
                            ((1+sigma(:,k_e-1,:))*(1/jaco(i_s:i_e,k_e-1,j_s:j_e)))/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)) + &
                            ((2+sigma(:,k_e-1,:))*(1/jaco(i_s:i_e,k_e,j_s:j_e)))/(dz_if(:,k_e,:)*(1+sigma(:,k_e-1,:)))
            
            dzhatdzz(:,k_s,:) = ((1/domain%jacobian_w(i_s:i_e,k_s,j_s:j_e)) - 1)/ &
                                       domain%advection_dz(i_s:i_e,k_s,j_s:j_e)
                                       
            dzhatdzz(:,k_s+1:k_e,:) = ((1/domain%jacobian_w(i_s:i_e,k_s+1:k_e,j_s:j_e)) - &
                                       (1/domain%jacobian_w(i_s:i_e,k_s:k_e-1,j_s:j_e)))/ &
                                       domain%advection_dz(i_s:i_e,k_s+1:k_e,j_s:j_e)
            
            !dzdxz(:,k_s+1:k_e-1,:) = (sigma(:,k_s+1:k_e-1,:)**2*dzdx(i_s:i_e,k_s+2:k_e,j_s:j_e) - &
            !        (sigma(:,k_s+1:k_e-1,:)**2-1)*dzdx(i_s:i_e,k_s+1:k_e-1,j_s:j_e)-dzdx(i_s:i_e,k_s:k_e-2,j_s:j_e)) / &
            !        (dz_if(:,k_s+2:k_e,:)*(sigma(:,k_s+1:k_e-1,:)+sigma(:,k_s+1:k_e-1,:)**2))
            !        
            !dzdxz(:,k_s,:) = -(sigma(:,k_s+1,:)**2*dzdx(i_s:i_e,k_s+2,j_s:j_e))/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1)) + &
            !                  ((sigma(:,k_s+1,:)+1)*dzdx(i_s:i_e,k_s+1,j_s:j_e))/dz_if(:,k_s+1,:) - &
            !                ((2*sigma(:,k_s+1,:)+1)*dzdx(i_s:i_e,k_s,j_s:j_e))/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1))
            !                
            !dzdxz(:,k_e,:) = dzdx(i_s:i_e,k_e-2,j_s:j_e)/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)*(1+sigma(:,k_e-1,:))) - &
            !                ((1+sigma(:,k_e-1,:))*dzdx(i_s:i_e,k_e-1,j_s:j_e))/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)) + &
            !                ((2+sigma(:,k_e-1,:))*dzdx(i_s:i_e,k_e,j_s:j_e))/(dz_if(:,k_e,:)*(1+sigma(:,k_e-1,:)))
!
            !        
            !dzdyz(:,k_s+1:k_e-1,:) = (sigma(:,k_s+1:k_e-1,:)**2*dzdy(i_s:i_e,k_s+2:k_e,j_s:j_e) - &
            !        (sigma(:,k_s+1:k_e-1,:)**2-1)*dzdy(i_s:i_e,k_s+1:k_e-1,j_s:j_e)-dzdy(i_s:i_e,k_s:k_e-2,j_s:j_e)) / &
            !        (dz_if(:,k_s+2:k_e,:)*(sigma(:,k_s+1:k_e-1,:)+sigma(:,k_s+1:k_e-1,:)**2))
            !
            !dzdyz(:,k_s,:) = -(sigma(:,k_s+1,:)**2*dzdy(i_s:i_e,k_s+2,j_s:j_e))/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1)) + &
            !                  ((sigma(:,k_s+1,:)+1)*dzdy(i_s:i_e,k_s+1,j_s:j_e))/dz_if(:,k_s+1,:) - &
            !                ((2*sigma(:,k_s+1,:)+1)*dzdy(i_s:i_e,k_s,j_s:j_e))/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1))
            !                
            !dzdyz(:,k_e,:) = dzdy(i_s:i_e,k_e-2,j_s:j_e)/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)*(1+sigma(:,k_e-1,:))) - &
            !                ((1+sigma(:,k_e-1,:))*dzdy(i_s:i_e,k_e-1,j_s:j_e))/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)) + &
            !                ((2+sigma(:,k_e-1,:))*dzdy(i_s:i_e,k_e,j_s:j_e))/(dz_if(:,k_e,:)*(1+sigma(:,k_e-1,:)))

            !dzhatdzz(:,k_s,:) = dzhatdzz(:,k_s+1,:)
            !dzhatdzz(:,k_e,:) = dzhatdzz(:,k_e-1,:)
            
            !dzdxz(:,k_s,:) = dzdxz(:,k_s+1,:)
            !dzdxz(:,k_e,:) = dzdxz(:,k_e-1,:)
            
            !dzdyz(:,k_s,:) = dzdyz(:,k_s+1,:)
            !dzdyz(:,k_e,:) = dzdyz(:,k_e-1,:)

            !Calculate how global grid is decomposed for DMDA
            !subtract halo size from boundaries of each cell to get x/y extent
            if (domain%grid%yimg == 1) xl(domain%grid%ximg) = domain%grid%nx-hs*2
            if (domain%grid%ximg == 1) yl(domain%grid%yimg) = domain%grid%ny-hs*2
        
            !Wait for all images to contribute their dimension            
            call CO_MAX(xl)
            call CO_MAX(yl)
            
            !Add points to xy-edges to accomodate ghost-points of DMDA grid
            !cells at boundaries have 1 extra for ghost-point, and should also be corrected
            !to have the hs which was falsely removed added back on
            xl(1) = xl(1)+hs
            xl(domain%grid%ximages) = xl(domain%grid%ximages)+hs

            yl(1) = yl(1)+hs
            yl(domain%grid%yimages) = yl(domain%grid%yimages)+hs
            
        endif

    
    end subroutine

end module wind_iterative
