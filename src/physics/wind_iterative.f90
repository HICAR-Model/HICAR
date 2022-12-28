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

    !use exchangeable_interface,   only : exchangeable_t
    use domain_interface,  only : domain_t
    !use options_interface, only : options_t
    !use grid_interface,    only : grid_t
    use petscksp
    use petscdm
    use petscdmda
    
    implicit none
    private
    public:: init_iter_winds, calc_iter_winds, finalize_iter_winds
    real, parameter::deg2rad=0.017453293 !2*pi/360
    real, parameter :: rad2deg=57.2957779371
    real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, H_coef, I_coef, &
                                           J_coef, K_coef, L_coef, M_coef, N_coef, O_coef
    real    :: dx
    real, allocatable, dimension(:,:,:)  :: div, dz_if, jaco, dzdx, dzdy, sigma, alpha
    integer, allocatable :: xl(:), yl(:)
    integer              :: hs, i_s, i_e, k_s, k_e, j_s, j_e
contains



    !>------------------------------------------------------------
    !! Forces u,v, and w fields to balance
    !!       du/dx + dv/dy = dw/dz
    !!
    !! Starts by setting w out of the ground=0 then works through layers
    !!
    !!------------------------------------------------------------
    subroutine calc_iter_winds(domain,alpha_in,div_in,adv_den,update_in)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, intent(in), dimension(domain%grid%ims:domain%grid%ime, &
                                    domain%grid%kms:domain%grid%kme, &
                                    domain%grid%jms:domain%grid%jme) :: alpha_in, div_in
        logical, intent(in) :: adv_den
        logical, optional, intent(in) :: update_in

        PetscScalar,pointer :: lambda(:,:,:)
        logical             :: update

        integer k !, i_s, i_e, k_s, k_e, j_s, j_e
        
        PetscErrorCode ierr
        KSP            ksp
        DM             da
        Vec            x, localX
        PetscInt       one, x_size
        PetscReal      norm, conv_tol

        update=.False.
        if (present(update_in)) update=update_in
                
                
        !Initialize div to be the initial divergence of the input wind field
        div = div_in(i_s:i_e,k_s:k_e,j_s:j_e) 
        one = 1
        
        alpha = alpha_in(i_s:i_e,k_s:k_e,j_s:j_e) 
        
        if (.not.(allocated(A_coef))) then
            call initialize_coefs(domain)
        else
            call update_coefs(domain)
        endif
                                
        !call init_iter_winds()
                                
        call KSPCreate(domain%IO_comms,ksp,ierr)
        conv_tol = 1e-4

        !call KSPSetTolerances(ksp,conv_tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        call KSPSetType(ksp,KSPBCGS,ierr);
        
        call DMDACreate3d(domain%IO_comms,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, &
                          (domain%ide+2),(domain%kde+2),(domain%jde+2),domain%grid%ximages,one,domain%grid%yimages,one,one, &
                          xl, PETSC_NULL_INTEGER,yl,da,ierr)
        
        call DMSetFromOptions(da,ierr)
        call DMSetUp(da,ierr)
        
        call KSPSetDM(ksp,da,ierr)
        call KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,0,ierr)
        call KSPSetComputeRHS(ksp,ComputeRHS,0,ierr)
        call KSPSetComputeOperators(ksp,ComputeMatrix,0,ierr)

        call KSPSetFromOptions(ksp,ierr)
        call DMCreateLocalVector(da,localX,ierr)
        call KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)
        
        call KSPGetSolution(ksp,x,ierr)
        if(this_image()==1) write(*,*) 'Solved PETSc'
        
        !Subset global solution x to local grid so that we can access ghost-points
        call DMGlobalToLocalBegin(da,x,INSERT_VALUES,localX,ierr)
        call DMGlobalToLocalEnd(da,x,INSERT_VALUES,localX,ierr)

        call DMDAVecGetArrayF90(da,localX,lambda, ierr)
        call calc_updated_winds(domain, lambda, update, adv_den)
        call DMDAVecRestoreArrayF90(da,localX,lambda, ierr)

        !Exchange u and v, since the outer points are not updated in above function
        call domain%u%exchange_x(update)
        call domain%v%exchange_y(update)
        
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


        real, allocatable, dimension(:,:,:)    :: u_dlambdz, v_dlambdz, u_temp, v_temp, lambda_too, rho, rho_u, rho_v
        integer k, i_start, i_end, j_start, j_end !i_s, i_e, k_s, k_e, j_s, j_e, ids, ide, jds, jde
                
        !i_s = domain%its-1
        !i_e = domain%ite+1
        !k_s = domain%kts  
        !k_e = domain%kte  
        !j_s = domain%jts-1
        !j_e = domain%jte+1


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

        allocate(u_temp(i_start:i_end,k_s-1:k_e+1,j_s:j_e))
        allocate(v_temp(i_s:i_e,k_s-1:k_e+1,j_start:j_end))

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
            rho_u(i_end,:,j_s:j_e) = rho(i_end-1,:,j_s:j_e)
            rho_u(i_start,:,j_s:j_e) = rho(i_start,:,j_s:j_e)
        else if (i_s==domain%grid%ids) then
            rho_u(i_start+1:i_end,:,j_s:j_e) = 0.5*(rho(i_start+1:i_end,:,j_s:j_e) + rho(i_start:i_end-1,:,j_s:j_e))
            rho_u(i_start,:,j_s:j_e) = rho(i_start,:,j_s:j_e)
        else if (i_e==domain%grid%ide) then
            rho_u(i_start:i_end-1,:,j_s:j_e) = 0.5*(rho(i_start:i_end-1,:,j_s:j_e) + rho(i_start-1:i_end-2,:,j_s:j_e))
            rho_u(i_end,:,j_s:j_e) = rho(i_end-1,:,j_s:j_e)
        else
            rho_u(i_start:i_end,:,j_s:j_e) = 0.5*(rho(i_start:i_end,:,j_s:j_e) + rho(i_start-1:i_end-1,:,j_s:j_e))
        endif
        
        if (j_s==domain%grid%jds .and. j_e==domain%grid%jde) then
            rho_v(i_s:i_e,:,j_start+1:j_end-1) = 0.5*(rho(i_s:i_e,:,j_start+1:j_end-1) + rho(i_s:i_e,:,j_start:j_end-2))
            rho_v(i_s:i_e,:,j_start) = rho(i_s:i_e,:,j_start)
            rho_v(i_s:i_e,:,j_end) = rho(i_s:i_e,:,j_end-1)
        else if (j_s==domain%grid%jds) then
            rho_v(i_s:i_e,:,j_start+1:j_end) = 0.5*(rho(i_s:i_e,:,j_start+1:j_end) + rho(i_s:i_e,:,j_start:j_end-1))
            rho_v(i_s:i_e,:,j_start) = rho(i_s:i_e,:,j_start)
        else if (j_e==domain%grid%jde) then
            rho_v(i_s:i_e,:,j_start:j_end-1) = 0.5*(rho(i_s:i_e,:,j_start:j_end-1) + rho(i_s:i_e,:,j_start-1:j_end-2))
            rho_v(i_s:i_e,:,j_end) = rho(i_s:i_e,:,j_end-1)
        else
            rho_v(i_s:i_e,:,j_start:j_end) = 0.5*(rho(i_s:i_e,:,j_start:j_end) + rho(i_s:i_e,:,j_start-1:j_end-1))
        endif
        
        !stager lambda to u grid
        u_temp = (lambda(i_start:i_end,k_s-1:k_e+1,j_s:j_e) + lambda(i_start-1:i_end-1,k_s-1:k_e+1,j_s:j_e)) / 2 

        !stager lambda to v grid
        v_temp = (lambda(i_s:i_e,k_s-1:k_e+1,j_start:j_end) + lambda(i_s:i_e,k_s-1:k_e+1,j_start-1:j_end-1)) / 2 

        !divide dz differennces by dz. Note that dz will be horizontally constant
        do k=k_s,k_e
            u_dlambdz(:,k,:) = (u_temp(:,k+1,:)*sigma(i_s,k,j_s)**2) - (sigma(i_s,k,j_s)**2 - 1)*u_temp(:,k,:) - u_temp(:,k-1,:)
            v_dlambdz(:,k,:) = (v_temp(:,k+1,:)*sigma(i_s,k,j_s)**2) - (sigma(i_s,k,j_s)**2 - 1)*v_temp(:,k,:) - v_temp(:,k-1,:)
        
            u_dlambdz(:,k,:) = u_dlambdz(:,k,:)/(dz_if(i_s,k+1,j_s)*(sigma(i_s,k,j_s)+sigma(i_s,k,j_s)**2))
            v_dlambdz(:,k,:) = v_dlambdz(:,k,:)/(dz_if(i_s,k+1,j_s)*(sigma(i_s,k,j_s)+sigma(i_s,k,j_s)**2))
        enddo
        
        u_dlambdz(:,k_s,:) = -(u_temp(:,k_s+2,:)*sigma(i_s,k_s+1,j_s)**2) + &
                            u_temp(:,k_s+1,:)*(sigma(i_s,k_s+1,j_s)+1)**2 - u_temp(:,k_s,:)*(2*sigma(i_s,k_s+1,j_s)+1)
        v_dlambdz(:,k_s,:) = -(v_temp(:,k_s+2,:)*sigma(i_s,k_s+1,j_s)**2) + &
                            v_temp(:,k_s+1,:)*(sigma(i_s,k_s+1,j_s)+1)**2 - v_temp(:,k_s,:)*(2*sigma(i_s,k_s+1,j_s)+1)
        
        u_dlambdz(:,k_s,:) = u_dlambdz(:,k_s,:)/(dz_if(i_s,k_s+1,j_s)*(sigma(i_s,k_s+1,j_s)+1))
        v_dlambdz(:,k_s,:) = v_dlambdz(:,k_s,:)/(dz_if(i_s,k_s+1,j_s)*(sigma(i_s,k_s+1,j_s)+1))
        
        !PETSc arrays are zero-indexed
        
        if (update) then
            domain%u%meta_data%dqdt_3d(i_start:i_end,:,j_s:j_e) = domain%u%meta_data%dqdt_3d(i_start:i_end,:,j_s:j_e) + &
                                                            0.5*((lambda(i_start:i_end,k_s:k_e,j_s:j_e) - &
                                                            lambda(i_start-1:i_end-1,k_s:k_e,j_s:j_e))/dx - &
            (1/domain%jacobian_u(i_start:i_end,:,j_s:j_e))*domain%dzdx_u(i_start:i_end,:,j_s:j_e)*(u_dlambdz))/rho_u(i_start:i_end,:,j_s:j_e)
            
            domain%v%meta_data%dqdt_3d(i_s:i_e,:,j_start:j_end) = domain%v%meta_data%dqdt_3d(i_s:i_e,:,j_start:j_end) + &
                                                            0.5*((lambda(i_s:i_e,k_s:k_e,j_start:j_end) - &
                                                            lambda(i_s:i_e,k_s:k_e,j_start-1:j_end-1))/dx - &
            (1/domain%jacobian_v(i_s:i_e,:,j_start:j_end))*domain%dzdy_v(i_s:i_e,:,j_start:j_end)*(v_dlambdz))/rho_v(i_s:i_e,:,j_start:j_end)
        else
            domain%u%data_3d(i_start:i_end,:,j_s:j_e) = domain%u%data_3d(i_start:i_end,:,j_s:j_e) + &
                                                            0.5*((lambda(i_start:i_end,k_s:k_e,j_s:j_e) - &
                                                            lambda(i_start-1:i_end-1,k_s:k_e,j_s:j_e))/dx - &
            (1/domain%jacobian_u(i_start:i_end,:,j_s:j_e))*domain%dzdx_u(i_start:i_end,:,j_s:j_e)*(u_dlambdz))/rho_u(i_start:i_end,:,j_s:j_e)
            
            domain%v%data_3d(i_s:i_e,:,j_start:j_end) = domain%v%data_3d(i_s:i_e,:,j_start:j_end) + &
                                                            0.5*((lambda(i_s:i_e,k_s:k_e,j_start:j_end) - &
                                                            lambda(i_s:i_e,k_s:k_e,j_start-1:j_end-1))/dx - &
            (1/domain%jacobian_v(i_s:i_e,:,j_start:j_end))*domain%dzdy_v(i_s:i_e,:,j_start:j_end)*(v_dlambdz))/rho_v(i_s:i_e,:,j_start:j_end)
        
        endif

    end subroutine calc_updated_winds

    subroutine ComputeRHS(ksp,vec_b,dummy,ierr)
        implicit none
        
        PetscErrorCode  ierr
        KSP ksp
        Vec vec_b, x
        integer dummy(*)
        DM             dm

        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs
        DMDALocalInfo       :: info(DMDA_LOCAL_INFO_SIZE)
        PetscScalar,pointer :: barray(:,:,:)

        call KSPGetDM(ksp,dm,ierr)
        
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
                    if (i.eq.0 .or. j.eq.0 .or. &
                        i.eq.mx-1 .or. j.eq.my-1 .or. k.eq.mz-1) then
                        barray(i,k,j) = 0.0
                    else if (k.eq.0) then
                        barray(i,k,j) = 0.0
                    else
                        barray(i,k,j) = -2*div(i,k,j)
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


        DM             da
        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,i1,i2,i10,i15
        DMDALocalInfo  info(DMDA_LOCAL_INFO_SIZE)
        PetscScalar    v(15)
        MatStencil     row(4),col(4,15),gnd_col(4,10),top_col(4,2)
        
        i1 = 1
        i2 = 2
        i10 = 10
        i15 = 15

        call KSPGetDM(ksp,da,ierr)

        call DMDAGetLocalInfoF90(da,info,ierr)
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
                if (i.eq.0 .or. j.eq.0 .or.k.eq.mz-1 .or. &
                    i.eq.mx-1 .or. j.eq.my-1) then
                    v(1) = 1.0
                    call MatSetValuesStencil(arr_B,i1,row,i1,row,v,INSERT_VALUES, ierr)
                else if (k.eq.0) then
                                
                    denom = 2*(1./alpha(i,1,j)**2 + dzdx(i,1,j)**2 + &
                                          dzdy(i,1,j)**2)/(jaco(i,1,j))
                    !k
                    v(1) = - 1
                    gnd_col(MatStencil_i,1) = i
                    gnd_col(MatStencil_j,1) = k
                    gnd_col(MatStencil_k,1) = j
                    !k + 1
                    v(2) = 1
                    gnd_col(MatStencil_i,2) = i
                    gnd_col(MatStencil_j,2) = k+1
                    gnd_col(MatStencil_k,2) = j
                    !i - 1
                    v(3) = dz_if(i,1,j)*dzdx(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,3) = i-1
                    gnd_col(MatStencil_j,3) = k+1
                    gnd_col(MatStencil_k,3) = j
                    !i - 1
                    v(4) = dz_if(i,1,j)*dzdx(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,4) = i-1
                    gnd_col(MatStencil_j,4) = k
                    gnd_col(MatStencil_k,4) = j
                    !i + 1
                    v(5) = -dz_if(i,1,j)*dzdx(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,5) = i+1
                    gnd_col(MatStencil_j,5) = k+1
                    gnd_col(MatStencil_k,5) = j
                    !i + 1
                    v(6) = -dz_if(i,1,j)*dzdx(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,6) = i+1
                    gnd_col(MatStencil_j,6) = k
                    gnd_col(MatStencil_k,6) = j
                    !j - 1
                    v(7) = dz_if(i,1,j)*dzdy(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,7) = i
                    gnd_col(MatStencil_j,7) = k+1
                    gnd_col(MatStencil_k,7) = j-1
                    !j - 1
                    v(8) = dz_if(i,1,j)*dzdy(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,8) = i
                    gnd_col(MatStencil_j,8) = k
                    gnd_col(MatStencil_k,8) = j-1
                    !j + 1
                    v(9) = -dz_if(i,1,j)*dzdy(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,9) = i
                    gnd_col(MatStencil_j,9) = k+1
                    gnd_col(MatStencil_k,9) = j+1
                    !j + 1
                    v(10) = -dz_if(i,1,j)*dzdy(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,10) = i
                    gnd_col(MatStencil_j,10) = k
                    gnd_col(MatStencil_k,10) = j+1
                    call MatSetValuesStencil(arr_B,i1,row,i10,gnd_col,v,INSERT_VALUES, ierr)
                else
                    !j - 1, k - 1
                    v(1) = O_coef(i,k,j)
                    col(MatStencil_i,1) = i
                    col(MatStencil_j,1) = k-1
                    col(MatStencil_k,1) = j-1
                    !i - 1, k - 1
                    v(2) = K_coef(i,k,j)
                    col(MatStencil_i,2) = i-1
                    col(MatStencil_j,2) = k-1
                    col(MatStencil_k,2) = j
                    !j - 1, k + 1
                    v(3) = M_coef(i,k,j)
                    col(MatStencil_i,3) = i
                    col(MatStencil_j,3) = k+1
                    col(MatStencil_k,3) = j-1
                    !i - 1, k + 1
                    v(4) = I_coef(i,k,j)
                    col(MatStencil_i,4) = i-1
                    col(MatStencil_j,4) = k+1
                    col(MatStencil_k,4) = j
                    !k - 1
                    v(5) = C_coef(i,k,j)
                    col(MatStencil_i,5) = i
                    col(MatStencil_j,5) = k-1
                    col(MatStencil_k,5) = j
                    !j - 1
                    v(6) = G_coef(i,k,j)
                    col(MatStencil_i,6) = i
                    col(MatStencil_j,6) = k
                    col(MatStencil_k,6) = j-1
                    !i - 1
                    v(7) = E_coef(i,k,j)
                    col(MatStencil_i,7) = i-1
                    col(MatStencil_j,7) = k
                    col(MatStencil_k,7) = j
                    !Center
                    v(8) = A_coef(i,k,j)
                    col(MatStencil_i,8) = i
                    col(MatStencil_j,8) = k
                    col(MatStencil_k,8) = j
                    !i + 1
                    v(9) = D_coef(i,k,j)
                    col(MatStencil_i,9) = i+1
                    col(MatStencil_j,9) = k
                    col(MatStencil_k,9) = j
                    !j + 1
                    v(10) = F_coef(i,k,j)
                    col(MatStencil_i,10) = i
                    col(MatStencil_j,10) = k
                    col(MatStencil_k,10) = j+1
                    !k + 1
                    v(11) = B_coef(i,k,j)
                    col(MatStencil_i,11) = i
                    col(MatStencil_j,11) = k+1
                    col(MatStencil_k,11) = j
                    !i + 1, k + 1
                    v(12) = H_coef(i,k,j)
                    col(MatStencil_i,12) = i+1
                    col(MatStencil_j,12) = k+1
                    col(MatStencil_k,12) = j
                    !j + 1, k + 1
                    v(13) = L_coef(i,k,j)
                    col(MatStencil_i,13) = i
                    col(MatStencil_j,13) = k+1
                    col(MatStencil_k,13) = j+1
                    !i + 1, k - 1
                    v(14) = J_coef(i,k,j)
                    col(MatStencil_i,14) = i+1
                    col(MatStencil_j,14) = k-1
                    col(MatStencil_k,14) = j
                    !j + 1, k - 1
                    v(15) = N_coef(i,k,j)
                    col(MatStencil_i,15) = i
                    col(MatStencil_j,15) = k-1
                    col(MatStencil_k,15) = j+1
                    call MatSetValuesStencil(arr_B,i1,row,i15,col,v,INSERT_VALUES, ierr)
                endif
            enddo
        enddo
        enddo

        call MatAssemblyBegin(arr_B,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(arr_B,MAT_FINAL_ASSEMBLY,ierr)
        if (arr_A .ne. arr_B) then
         call MatAssemblyBegin(arr_A,MAT_FINAL_ASSEMBLY,ierr)
         call MatAssemblyEnd(arr_A,MAT_FINAL_ASSEMBLY,ierr)
        endif

    end subroutine ComputeMatrix
    
    
    subroutine initialize_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        
        real, allocatable, dimension(:,:,:) :: mixed_denom
        !integer i_s, i_e, k_s, k_e, j_s, j_e
        
        !i_s = domain%its-1
        !i_e = domain%ite+1
        !k_s = domain%kts  
        !k_e = domain%kte  
        !j_s = domain%jts-1
        !j_e = domain%jte+1
        
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

        mixed_denom = 2*domain%dx*dz_if(:,k_s+1:k_e+1,:)*(sigma+sigma**2)

        !This function sets A, B, annd C coefficients
        call update_coefs(domain)
                
        D_coef(i_s:i_e-1,k_s:k_e,j_s:j_e) = (domain%jacobian(i_s+1:i_e,k_s:k_e,j_s:j_e) + &
            domain%jacobian(i_s:i_e-1,k_s:k_e,j_s:j_e))/(2*domain%dx**2) + &
            (sigma(i_s:i_e-1,:,:)**2 - 1)*(dzdx(i_s:i_e-1,:,:)+domain%dzdx_u(i_s+1:i_e,:,j_s:j_e))/&
            mixed_denom(i_s:i_e-1,k_s:k_e,j_s:j_e)
        E_coef(i_s+1:i_e,k_s:k_e,j_s:j_e) = (domain%jacobian(i_s+1:i_e,k_s:k_e,j_s:j_e) + &
            domain%jacobian(i_s:i_e-1,k_s:k_e,j_s:j_e))/(2*domain%dx**2) - &
            (sigma(i_s+1:i_e,:,:)**2 - 1)*(dzdx(i_s+1:i_e,:,:)+domain%dzdx_u(i_s+1:i_e,:,j_s:j_e))/&
            mixed_denom(i_s+1:i_e,k_s:k_e,j_s:j_e)
        F_coef(i_s:i_e,k_s:k_e,j_s:j_e-1) = (domain%jacobian(i_s:i_e,k_s:k_e,j_s+1:j_e) + &
            domain%jacobian(i_s:i_e,k_s:k_e,j_s:j_e-1))/(2*domain%dx**2) + &
            (sigma(:,:,j_s:j_e-1)**2 - 1)*(dzdy(:,:,j_s:j_e-1)+domain%dzdy_v(i_s:i_e,:,j_s+1:j_e))/&
            mixed_denom(i_s:i_e,k_s:k_e,j_s:j_e-1)
        G_coef(i_s:i_e,k_s:k_e,j_s+1:j_e) = (domain%jacobian(i_s:i_e,k_s:k_e,j_s+1:j_e) + &
            domain%jacobian(i_s:i_e,k_s:k_e,j_s:j_e-1))/(2*domain%dx**2) - &
            (sigma(:,:,j_s+1:j_e)**2 - 1)*(dzdy(:,:,j_s+1:j_e)+domain%dzdy_v(i_s:i_e,:,j_s+1:j_e))/&
            mixed_denom(i_s:i_e,k_s:k_e,j_s+1:j_e)

        D_coef(i_e,k_s:k_e,j_s:j_e) = D_coef(i_e-1,k_s:k_e,j_s:j_e)
        E_coef(i_s,k_s:k_e,j_s:j_e) = E_coef(i_s+1,k_s:k_e,j_s:j_e)
        F_coef(i_s:i_e,k_s:k_e,j_e) = F_coef(i_s:i_e,k_s:k_e,j_e-1)
        G_coef(i_s:i_e,k_s:k_e,j_s) = G_coef(i_s:i_e,k_s:k_e,j_s+1)

        H_coef(i_s:i_e-1,k_s:k_e-1,:) = &
                -(sigma(i_s:i_e-1,k_s:k_e-1,:)**2)* &
                (dzdx(i_s:i_e-1,k_s+1:k_e,:)+domain%dzdx_u(i_s+1:i_e,k_s:k_e-1,j_s:j_e))&
                /mixed_denom(i_s:i_e-1,k_s:k_e-1,:)
        I_coef(i_s+1:i_e,k_s:k_e-1,:) = &
                (sigma(i_s+1:i_e,k_s:k_e-1,:)**2)* &
                (dzdx(i_s+1:i_e,k_s+1:k_e,:)+&
                domain%dzdx_u(i_s+1:i_e,k_s:k_e-1,j_s:j_e))/mixed_denom(i_s+1:i_e,k_s:k_e-1,:)
        J_coef(i_s:i_e-1,k_s+1:k_e,:) = &
                (dzdx(i_s:i_e-1,k_s:k_e-1,:)+domain%dzdx_u(i_s+1:i_e,k_s+1:k_e,j_s:j_e))&
                /mixed_denom(i_s:i_e-1,k_s+1:k_e,:)
        K_coef(i_s+1:i_e,k_s+1:k_e,:) = &
                -(dzdx(i_s+1:i_e,k_s:k_e-1,:)+&
                domain%dzdx_u(i_s+1:i_e,k_s+1:k_e,j_s:j_e))/mixed_denom(i_s+1:i_e,k_s+1:k_e,:)

        L_coef(:,k_s:k_e-1,j_s:j_e-1) = &
                -(sigma(:,k_s:k_e-1,j_s:j_e-1)**2)* &
                (dzdy(:,k_s+1:k_e,j_s:j_e-1)+domain%dzdy_v(i_s:i_e,k_s:k_e-1,j_s+1:j_e))/mixed_denom(:,k_s:k_e-1,j_s:j_e-1)
        M_coef(:,k_s:k_e-1,j_s+1:j_e) = &
                (sigma(:,k_s:k_e-1,j_s+1:j_e)**2)* &
                (dzdy(:,k_s+1:k_e,j_s+1:j_e)+domain%dzdy_v(i_s:i_e,k_s:k_e-1,j_s+1:j_e))/mixed_denom(:,k_s:k_e-1,j_s+1:j_e)
        N_coef(:,k_s+1:k_e,j_s:j_e-1) = &
                (dzdy(:,k_s:k_e-1,j_s:j_e-1)+domain%dzdy_v(i_s:i_e,k_s+1:k_e,j_s+1:j_e))/mixed_denom(:,k_s+1:k_e,j_s:j_e-1)
        O_coef(:,k_s+1:k_e,j_s+1:j_e) = &
                -(dzdy(:,k_s:k_e-1,j_s+1:j_e)+domain%dzdy_v(i_s:i_e,k_s+1:k_e,j_s+1:j_e))/&
                mixed_denom(:,k_s+1:k_e,j_s+1:j_e)

        J_coef(i_s:i_e-1,k_s,:) = (dzdx(i_s:i_e-1,k_s,:)+domain%dzdx_u(i_s+1:i_e,k_s,j_s:j_e))&
                                    /mixed_denom(i_s:i_e-1,k_s,:)
        K_coef(i_s+1:i_e,k_s,:) = -(dzdx(i_s+1:i_e,k_s,:)+domain%dzdx_u(i_s+1:i_e,k_s,j_s:j_e))&
                                    /mixed_denom(i_s+1:i_e,k_s,:)
        N_coef(:,k_s,j_s:j_e-1) = (dzdy(:,k_s,j_s:j_e-1)+domain%dzdy_v(i_s:i_e,k_s,j_s+1:j_e))&
                                    /mixed_denom(:,k_s,j_s:j_e-1)
        O_coef(:,k_s,j_s+1:j_e) = -(dzdy(:,k_s,j_s+1:j_e)+domain%dzdy_v(i_s:i_e,k_s,j_s+1:j_e))&
                                    /mixed_denom(:,k_s,j_s+1:j_e)

        H_coef(i_e,:,:) = H_coef(i_e-1,:,:)
        I_coef(i_s,:,:) = I_coef(i_s+1,:,:)
        J_coef(i_e,:,:) = J_coef(i_e-1,:,:)
        K_coef(i_s,:,:) = K_coef(i_s+1,:,:)

        L_coef(:,:,j_e) = L_coef(:,:,j_e-1)
        M_coef(:,:,j_s) = M_coef(:,:,j_s+1)
        N_coef(:,:,j_e) = N_coef(:,:,j_e-1)
        O_coef(:,:,j_s) = O_coef(:,:,j_s+1)

    end subroutine initialize_coefs
    
    
    !Update the coefs which change with time, i.e. those which depend on alpha
    subroutine update_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain        
        real, allocatable, dimension(:,:,:) :: mixed_denom
        !integer i_s, i_e, k_s, k_e, j_s, j_e
        
        !i_s = domain%its-1
        !i_e = domain%ite+1
        !k_s = domain%kts  
        !k_e = domain%kte  
        !j_s = domain%jts-1
        !j_e = domain%jte+1
        
        allocate(mixed_denom(i_s:i_e,k_s:k_e,j_s:j_e))
        mixed_denom = 2*domain%dx*dz_if(:,k_s+1:k_e+1,:)*(sigma+sigma**2)


        B_coef(:,k_s:k_e-1,:) = sigma(:,k_s:k_e-1,:) * &
                              ( (1./alpha(i_s:i_e,k_s:k_e-1,j_s:j_e)**2 + dzdy(i_s:i_e,k_s:k_e-1,j_s:j_e)**2 + &
                              dzdx(i_s:i_e,k_s:k_e-1,j_s:j_e)**2) * (1./domain%jacobian(i_s:i_e,k_s:k_e-1,j_s:j_e)) + &
                              (1./alpha(i_s:i_e,k_s+1:k_e,j_s:j_e)**2 + dzdy(i_s:i_e,k_s+1:k_e,j_s:j_e)**2 + &
                              dzdx(i_s:i_e,k_s+1:k_e,j_s:j_e)**2) * (1./domain%jacobian(i_s:i_e,k_s+1:k_e,j_s:j_e))) / &
                          ((sigma(:,k_s:k_e-1,:)+sigma(:,k_s:k_e-1,:)**2)*dz_if(:,k_s+1:k_e,:)**2)
                          
                          
        C_coef(:,k_s+1:k_e,:) = ( (1./alpha(i_s:i_e,k_s:k_e-1,j_s:j_e)**2 + dzdy(i_s:i_e,k_s:k_e-1,j_s:j_e)**2 + &
                              dzdx(i_s:i_e,k_s:k_e-1,j_s:j_e)**2) * (1./domain%jacobian(i_s:i_e,k_s:k_e-1,j_s:j_e)) + &
                              (1./alpha(i_s:i_e,k_s+1:k_e,j_s:j_e)**2 + dzdy(i_s:i_e,k_s+1:k_e,j_s:j_e)**2 + &
                              dzdx(i_s:i_e,k_s+1:k_e,j_s:j_e)**2) * (1./domain%jacobian(i_s:i_e,k_s+1:k_e,j_s:j_e))) / &
                          ((sigma(:,k_s+1:k_e,:)+sigma(:,k_s+1:k_e,:)**2)*dz_if(:,k_s+2:k_e+1,:)**2)
                
        C_coef(:,k_s,:) = ( (1./alpha(i_s:i_e,k_s,j_s:j_e)**2 + dzdy(i_s:i_e,k_s,j_s:j_e)**2 + &
                              dzdx(i_s:i_e,k_s,j_s:j_e)**2) * (1./domain%jacobian(i_s:i_e,k_s,j_s:j_e)) + &
                              (1./alpha(i_s:i_e,k_s,j_s:j_e)**2 + dzdy(i_s:i_e,k_s,j_s:j_e)**2 + &
                              dzdx(i_s:i_e,k_s,j_s:j_e)**2) * (1./domain%jacobian(i_s:i_e,k_s,j_s:j_e))) / &
                          ((sigma(:,k_s,:)+sigma(:,k_s,:)**2)*dz_if(:,k_s+1,:)**2)
                          
        B_coef(i_s+1:i_e-1,:,:) = B_coef(i_s+1:i_e-1,:,:) - sigma(i_s+1:i_e-1,:,:)**2 * &
                                    (domain%dzdx_u(i_s+2:i_e,:,j_s:j_e)-domain%dzdx_u(i_s+1:i_e-1,:,j_s:j_e))/(mixed_denom(i_s+1:i_e-1,:,:))
                          
        B_coef(:,:,j_s+1:j_e-1) = B_coef(:,:,j_s+1:j_e-1) - sigma(:,:,j_s+1:j_e-1)**2 * &
                                    (domain%dzdy_v(i_s:i_e,:,j_s+2:j_e)-domain%dzdy_v(i_s:i_e,:,j_s+1:j_e-1))/(mixed_denom(:,:,j_s+1:j_e-1))
                          
        
        C_coef(i_s+1:i_e-1,:,:) = C_coef(i_s+1:i_e-1,:,:) + &
                                    (domain%dzdx_u(i_s+2:i_e,:,j_s:j_e)-domain%dzdx_u(i_s+1:i_e-1,:,j_s:j_e))/(mixed_denom(i_s+1:i_e-1,:,:))
                          
        C_coef(:,:,j_s+1:j_e-1) = C_coef(:,:,j_s+1:j_e-1) + &
                                    (domain%dzdy_v(i_s:i_e,:,j_s+2:j_e)-domain%dzdy_v(i_s:i_e,:,j_s+1:j_e-1))/(mixed_denom(:,:,j_s+1:j_e-1))
                                                    
        B_coef(i_s,:,:) = B_coef(i_s+1,:,:)
        B_coef(i_e,:,:) = B_coef(i_e-1,:,:)
        B_coef(:,:,j_s) = B_coef(:,:,j_s+1)
        B_coef(:,:,j_e) = B_coef(:,:,j_e-1)
        
        C_coef(i_s,:,:) = C_coef(i_s+1,:,:)
        C_coef(i_e,:,:) = C_coef(i_e-1,:,:)
        C_coef(:,:,j_s) = C_coef(:,:,j_s+1)
        C_coef(:,:,j_e) = C_coef(:,:,j_e-1)
        
        A_coef(i_s+1:i_e-1,k_s:k_e,j_s+1:j_e-1) = -((domain%jacobian(i_s+2:i_e,k_s:k_e,j_s+1:j_e-1) + &
                                               2*domain%jacobian(i_s+1:i_e-1,k_s:k_e,j_s+1:j_e-1) + &
                                               domain%jacobian(i_s:i_e-2,k_s:k_e,j_s+1:j_e-1))/(2*domain%dx**2)) &
                                            -((domain%jacobian(i_s+1:i_e-1,k_s:k_e,j_s+2:j_e) + &
                                               2*domain%jacobian(i_s+1:i_e-1,k_s:k_e,j_s+1:j_e-1) + &
                                               domain%jacobian(i_s+1:i_e-1,k_s:k_e,j_s:j_e-2))/(2*domain%dx**2)) - &
                                               B_coef(i_s+1:i_e-1,k_s:k_e,j_s+1:j_e-1) - C_coef(i_s+1:i_e-1,k_s:k_e,j_s+1:j_e-1)
                                               
        A_coef(i_s,:,j_s+1:j_e-1) = A_coef(i_s+1,:,j_s+1:j_e-1) 
        A_coef(i_e,:,j_s+1:j_e-1) = A_coef(i_e-1,:,j_s+1:j_e-1) 
        A_coef(i_s:i_e,:,j_s) = A_coef(i_s:i_e,:,j_s+1)
        A_coef(i_s:i_e,:,j_e) = A_coef(i_s:i_e,:,j_e-1)
                    
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
            allocate(sigma(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(dz_if(i_s:i_e,k_s:k_e+1,j_s:j_e))
            allocate(alpha(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(div(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate( xl( 1:domain%grid%ximages ))
            allocate( yl( 1:domain%grid%yimages ))
            xl = 0
            yl = 0
            
            dx = domain%dx
            dzdx = domain%dzdx(i_s:i_e,k_s:k_e,j_s:j_e) 
            dzdy = domain%dzdy(i_s:i_e,k_s:k_e,j_s:j_e)
            jaco = domain%jacobian(i_s:i_e,k_s:k_e,j_s:j_e)
            
            dz_if(:,k_s,:) = domain%advection_dz(i_s:i_e,k_s,j_s:j_e)
            dz_if(:,k_s+1:k_e,:) = (domain%advection_dz(i_s:i_e,k_s+1:k_e,j_s:j_e) + &
                                   domain%advection_dz(i_s:i_e,k_s:k_e-1,j_s:j_e))/2
            dz_if(:,k_e+1,:) = domain%advection_dz(i_s:i_e,k_e,j_s:j_e)
            sigma = dz_if(:,k_s:k_e,:)/dz_if(:,k_s+1:k_e+1,:)
            
                    
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
            xl(1) = xl(1)+1+hs
            xl(domain%grid%ximages) = xl(domain%grid%ximages)+1+hs

            yl(1) = yl(1)+1+hs
            yl(domain%grid%yimages) = yl(domain%grid%yimages)+1+hs
            
        endif

    
    end subroutine

end module wind_iterative