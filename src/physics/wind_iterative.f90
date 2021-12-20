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
    public::calc_iter_winds
    real, parameter::deg2rad=0.017453293 !2*pi/360
    real, parameter :: rad2deg=57.2957779371
    real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, H_coef, I_coef, &
                                           J_coef, K_coef, L_coef, M_coef, N_coef, O_coef
    real    :: dx
    real, allocatable, dimension(:,:,:)  :: div, dz_if, jaco, dzdx, dzdy, sigma, alpha
    integer, allocatable :: xl(:), yl(:)
contains



    !>------------------------------------------------------------
    !! Forces u,v, and w fields to balance
    !!       du/dx + dv/dy = dw/dz
    !!
    !! Starts by setting w out of the ground=0 then works through layers
    !!
    !!------------------------------------------------------------
    subroutine calc_iter_winds(domain,alpha_in,div_in,update_in)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, intent(in)              :: alpha_in(:,:,:), div_in(:,:,:)
        logical, optional, intent(in) :: update_in
        
        PetscScalar,pointer :: lambda(:,:,:)
        logical             :: update
        integer k, ims, ime, kms, kme, jms, jme
        
        PetscErrorCode ierr
        KSP            ksp
        PetscReal      norm, conv_tol
        DM             da
        Vec            x, localX
        PetscInt       one, x_size
        
        update=.False.
        if (present(update_in)) update=update_in
                
        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)

        if (.not.(allocated(dzdx))) then
            allocate(dzdx(ims:ime,kms:kme,jms:jme))
            allocate(dzdy(ims:ime,kms:kme,jms:jme))
            allocate(jaco(ims:ime,kms:kme,jms:jme))
            allocate(sigma(ims:ime,kms:kme,jms:jme))
            allocate(dz_if(ims:ime,kms:kme+1,jms:jme))
            allocate(alpha(ims:ime,kms:kme,jms:jme))
            allocate(div(ims:ime,kms:kme,jms:jme))
            allocate( xl( 1:domain%grid%ximages ))
            allocate( yl( 1:domain%grid%yimages ))
            xl = 0
            yl = 0
            
            dx = domain%dx
            dzdx = domain%dzdx 
            dzdy = domain%dzdy
            jaco = domain%jacobian
            
            dz_if(:,kms,:) = domain%advection_dz(:,kms,:)
            dz_if(:,kms+1:kme,:) = (domain%advection_dz(:,kms+1:kme,:)/2) + (domain%advection_dz(:,kms:kme-1,:)/2)
            dz_if(:,kme+1,:) = domain%advection_dz(:,kme,:)
            sigma = dz_if(:,kms:kme,:)/dz_if(:,kms+1:kme+1,:)
            
                    
            !Calculate how global grid is decomposed for DMDA
            if (domain%grid%yimg == 1) xl(domain%grid%ximg) = domain%grid%nx-2
            if (domain%grid%ximg == 1) yl(domain%grid%yimg) = domain%grid%ny-2
        
            !Wait for all images to contribute their dimension
            sync all
            
            call CO_MAX(xl)
            call CO_MAX(yl)
            
            !Add points to xy-edges to accomodate ghost-points of DMDA grid
            xl(1) = xl(1)+2
            xl(domain%grid%ximages) = xl(domain%grid%ximages)+2

            yl(1) = yl(1)+2
            yl(domain%grid%yimages) = yl(domain%grid%yimages)+2
            
        endif
                
        !Initialize div to be the initial divergence of the input wind field
        div=div_in
        one = 1
        
        alpha = alpha_in
        
        if (.not.(allocated(A_coef))) then
            call initialize_coefs(domain)
        else
            call update_coefs(domain)
        endif

        
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        if (ierr .ne. 0) then
            print*,'Unable to initialize PETSc'
            stop
        endif 
        if(this_image()==1) write(*,*) 'Initialized PETSc'
        
        call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
        conv_tol = 1e-4

        
        !call KSPSetTolerances(ksp,conv_tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        call KSPSetType(ksp,KSPBCGS,ierr);
        call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, &
                          (domain%ide+2),(domain%kde+2),(domain%jde+2),domain%grid%ximages,one,domain%grid%yimages,one,one, &
                          xl, PETSC_NULL_INTEGER,yl,da,ierr)
        
        call DMSetFromOptions(da,ierr)
        call DMSetUp(da,ierr)
        
        call KSPSetDM(ksp,da,ierr)
        call KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,0,ierr)
        call KSPSetComputeRHS(ksp,ComputeRHS,0,ierr)
        call KSPSetComputeOperators(ksp,ComputeMatrix,0,ierr)

        call KSPSetFromOptions(ksp,ierr)
        call DMCreateGlobalVector(da,x,ierr)
        call DMCreateLocalVector(da,localX,ierr)
        call KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)
        
        call KSPGetSolution(ksp,x,ierr)
        if(this_image()==1) write(*,*) 'Solved PETSc'
        
        !Subset global solution x to local grid so that we can access ghost-points
        call DMGlobalToLocalBegin(da,x,INSERT_VALUES,localX,ierr)
        call DMGlobalToLocalEnd(da,x,INSERT_VALUES,localX,ierr)

        call DMDAVecGetArrayF90(da,localX,lambda, ierr)
        
        call calc_updated_winds(domain, lambda, update)
        
        !Exchange u and v, since the outer points are not updated in above function
        if (update) then 
            call domain%u%exchange_u_metadata()
            call domain%v%exchange_v_metadata()
        else
            call domain%u%exchange_u()
            call domain%v%exchange_v()
        endif
        
        call DMDAVecRestoreArrayF90(da,x,lambda, ierr)
        call DMDestroy(da,ierr)
        call KSPDestroy(ksp,ierr)
        call PetscFinalize(ierr)

    end subroutine calc_iter_winds
    
    subroutine calc_updated_winds(domain,lambda,update) !u, v, w, jaco_u,jaco_v,jaco_w,u_dzdx,v_dzdy,lambda, ids, ide, jds, jde)
        type(domain_t), intent(inout) :: domain
        !real, intent(inout), dimension(:,:,:)  :: u,v,w
        !real, intent(in), dimension(:,:,:)     :: jaco_u,jaco_v,jaco_w, u_dzdx, v_dzdy
        PetscScalar, intent(in), pointer       :: lambda(:,:,:)
        logical,     intent(in)                :: update

        real, allocatable, dimension(:,:,:)    :: u_dlambdz, v_dlambdz, u_temp, v_temp, lambda_too
        integer k, i_start, i_end, j_start, j_end, ims, ime, kms, kme, jms, jme, ids, ide, jds, jde
                
        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)


        !ims+1, unless we are on global boundary, then ims
        i_start = ims+1
        if (ims==domain%grid%ids) i_start = ims
        
        !ime, unless we are on global boundary, then ime+1
        i_end = ime
        if (ime==domain%grid%ide) i_end = ime+1
        
        !jms+1, unless we are on global boundary, then jms
        j_start = jms+1
        if (jms==domain%grid%jds) j_start = jms
        
        !jme, unless we are on global boundary, then jme+1
        j_end = jme
        if (jme==domain%grid%jde) j_end = jme+1

        allocate(u_temp(i_start:i_end,kms-1:kme+1,jms:jme))
        allocate(v_temp(ims:ime,kms-1:kme+1,j_start:j_end))

        allocate(u_dlambdz(i_start:i_end,kms:kme,jms:jme))
        allocate(v_dlambdz(ims:ime,kms:kme,j_start:j_end))
        
        !stager lambda to u grid
        u_temp = (lambda(i_start:i_end,kms-1:kme+1,jms:jme) + lambda(i_start-1:i_end-1,kms-1:kme+1,jms:jme)) / 2 

        !stager lambda to v grid
        v_temp = (lambda(ims:ime,kms-1:kme+1,j_start:j_end) + lambda(ims:ime,kms-1:kme+1,j_start-1:j_end-1)) / 2 

        !divide dz differennces by dz. Note that dz will be horizontally constant
        do k=kms,kme
            u_dlambdz(:,k,:) = (u_temp(:,k+1,:)*sigma(ims,k,jms)**2) - (sigma(ims,k,jms)**2 - 1)*u_temp(:,k,:) - u_temp(:,k-1,:)
            v_dlambdz(:,k,:) = (v_temp(:,k+1,:)*sigma(ims,k,jms)**2) - (sigma(ims,k,jms)**2 - 1)*v_temp(:,k,:) - v_temp(:,k-1,:)
        
            u_dlambdz(:,k,:) = u_dlambdz(:,k,:)/(dz_if(ims,k+1,jms)*(sigma(ims,k,jms)+sigma(ims,k,jms)**2))
            v_dlambdz(:,k,:) = v_dlambdz(:,k,:)/(dz_if(ims,k+1,jms)*(sigma(ims,k,jms)+sigma(ims,k,jms)**2))
        enddo
        
        u_dlambdz(:,kms,:) = -(u_temp(:,kms+2,:)*sigma(ims,kms+1,jms)**2) + &
                            u_temp(:,kms+1,:)*(sigma(ims,kms+1,jms)+1)**2 - u_temp(:,kms,:)*(2*sigma(ims,kms+1,jms)+1)
        v_dlambdz(:,kms,:) = -(v_temp(:,kms+2,:)*sigma(ims,kms+1,jms)**2) + &
                            v_temp(:,kms+1,:)*(sigma(ims,kms+1,jms)+1)**2 - v_temp(:,kms,:)*(2*sigma(ims,kms+1,jms)+1)
        
        u_dlambdz(:,kms,:) = u_dlambdz(:,kms,:)/(dz_if(ims,kms+1,jms)*(sigma(ims,kms+1,jms)+1))
        v_dlambdz(:,kms,:) = v_dlambdz(:,kms,:)/(dz_if(ims,kms+1,jms)*(sigma(ims,kms+1,jms)+1))
        
        !PETSc arrays are zero-indexed
        
        if (update) then
            domain%u%meta_data%dqdt_3d(i_start:i_end,:,:) = domain%u%meta_data%dqdt_3d(i_start:i_end,:,:) + &
                                                            0.5*((lambda(i_start:i_end,kms:kme,jms:jme) - &
                                                            lambda(i_start-1:i_end-1,kms:kme,jms:jme))/dx - &
                                      (1/domain%jacobian_u(i_start:i_end,:,:))*domain%dzdx_u(i_start:i_end,:,:)*(u_dlambdz))
            domain%v%meta_data%dqdt_3d(:,:,j_start:j_end) = domain%v%meta_data%dqdt_3d(:,:,j_start:j_end) + &
                                                            0.5*((lambda(ims:ime,kms:kme,j_start:j_end) - &
                                                            lambda(ims:ime,kms:kme,j_start-1:j_end-1))/dx - &
                                      (1/domain%jacobian_v(:,:,j_start:j_end))*domain%dzdy_v(:,:,j_start:j_end)*(v_dlambdz))
        else
            domain%u%data_3d(i_start:i_end,:,:) = domain%u%data_3d(i_start:i_end,:,:) + &
                                                            0.5*((lambda(i_start:i_end,kms:kme,jms:jme) - &
                                                            lambda(i_start-1:i_end-1,kms:kme,jms:jme))/dx - &
                                      (1/domain%jacobian_u(i_start:i_end,:,:))*domain%dzdx_u(i_start:i_end,:,:)*(u_dlambdz))
            domain%v%data_3d(:,:,j_start:j_end) = domain%v%data_3d(:,:,j_start:j_end) + &
                                                            0.5*((lambda(ims:ime,kms:kme,j_start:j_end) - &
                                                            lambda(ims:ime,kms:kme,j_start-1:j_end-1))/dx - &
                                      (1/domain%jacobian_v(:,:,j_start:j_end))*domain%dzdy_v(:,:,j_start:j_end)*(v_dlambdz))
        
        endif

    end subroutine calc_updated_winds

    subroutine ComputeRHS(ksp,vec_b,dummy,ierr)
        implicit none
        
        PetscErrorCode  ierr
        KSP ksp
        Vec vec_b, x
        integer dummy(*)
        DM             dm

        PetscInt       i,j,k,HICAR_i,HICAR_j,HICAR_k,mx,my,mz,xm,ym,zm,xs,ys,zs
        DMDALocalInfo       :: info(DMDA_LOCAL_INFO_SIZE)
        PetscScalar,pointer :: barray(:,:,:), lambda(:,:,:)
        real                :: dlambdx, dlambdy

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
                    HICAR_i = i
                    HICAR_j = j
                    HICAR_k = k
                    !For global boundary conditions
                    if (i.eq.0 .or. j.eq.0 .or. &
                        i.eq.mx-1 .or. j.eq.my-1 .or. k.eq.mz-1) then
                        barray(i,k,j) = 0.0
                    else if (k.eq.0) then
                        barray(i,k,j) = 0.0
                    else
                        barray(i,k,j) = -2*div(HICAR_i,HICAR_k,HICAR_j)
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
        integer ims, ime, kms, kme, jms, jme
        
        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)
        
        allocate(A_coef(ims:ime,kms:kme,jms:jme))
        allocate(B_coef(ims:ime,kms:kme,jms:jme))
        allocate(C_coef(ims:ime,kms:kme,jms:jme))
        allocate(D_coef(ims:ime,kms:kme,jms:jme))
        allocate(E_coef(ims:ime,kms:kme,jms:jme))
        allocate(F_coef(ims:ime,kms:kme,jms:jme))
        allocate(G_coef(ims:ime,kms:kme,jms:jme))
        allocate(H_coef(ims:ime,kms:kme,jms:jme))
        allocate(I_coef(ims:ime,kms:kme,jms:jme))
        allocate(J_coef(ims:ime,kms:kme,jms:jme))
        allocate(K_coef(ims:ime,kms:kme,jms:jme))
        allocate(L_coef(ims:ime,kms:kme,jms:jme))
        allocate(M_coef(ims:ime,kms:kme,jms:jme))
        allocate(N_coef(ims:ime,kms:kme,jms:jme))
        allocate(O_coef(ims:ime,kms:kme,jms:jme))
        
        allocate(mixed_denom(ims:ime,kms:kme,jms:jme))
        
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

        mixed_denom = 2*domain%dx*dz_if(:,kms+1:kme+1,:)*(sigma+sigma**2)

        !This function sets A, B, annd C coefficients
        call update_coefs(domain)
                
        D_coef(ims:ime-1,kms:kme,jms:jme) = (domain%jacobian(ims+1:ime,kms:kme,jms:jme) + &
            domain%jacobian(ims:ime-1,kms:kme,jms:jme))/(2*domain%dx**2) + &
            (sigma(ims:ime-1,:,:)**2 - 1)*(dzdx(ims:ime-1,:,:)+domain%dzdx_u(ims+1:ime,:,:))/&
            mixed_denom(ims:ime-1,kms:kme,jms:jme)
        E_coef(ims+1:ime,kms:kme,jms:jme) = (domain%jacobian(ims+1:ime,kms:kme,jms:jme) + &
            domain%jacobian(ims:ime-1,kms:kme,jms:jme))/(2*domain%dx**2) - &
            (sigma(ims+1:ime,:,:)**2 - 1)*(dzdx(ims+1:ime,:,:)+domain%dzdx_u(ims+1:ime,:,:))/&
            mixed_denom(ims+1:ime,kms:kme,jms:jme)
        F_coef(ims:ime,kms:kme,jms:jme-1) = (domain%jacobian(ims:ime,kms:kme,jms+1:jme) + &
            domain%jacobian(ims:ime,kms:kme,jms:jme-1))/(2*domain%dx**2) + &
            (sigma(:,:,jms:jme-1)**2 - 1)*(dzdy(:,:,jms:jme-1)+domain%dzdy_v(:,:,jms+1:jme))/&
            mixed_denom(ims:ime,kms:kme,jms:jme-1)
        G_coef(ims:ime,kms:kme,jms+1:jme) = (domain%jacobian(ims:ime,kms:kme,jms+1:jme) + &
            domain%jacobian(ims:ime,kms:kme,jms:jme-1))/(2*domain%dx**2) - &
            (sigma(:,:,jms+1:jme)**2 - 1)*(dzdy(:,:,jms+1:jme)+domain%dzdy_v(:,:,jms+1:jme))/&
            mixed_denom(ims:ime,kms:kme,jms+1:jme)
            
        D_coef(ime,kms:kme,jms:jme) = D_coef(ime-1,kms:kme,jms:jme)
        E_coef(ims,kms:kme,jms:jme) = E_coef(ims+1,kms:kme,jms:jme)
        F_coef(ims:ime,kms:kme,jme) = F_coef(ims:ime,kms:kme,jme-1)
        G_coef(ims:ime,kms:kme,jms) = G_coef(ims:ime,kms:kme,jms+1)
        
        H_coef(ims:ime-1,kms:kme-1,:) = &
                -(sigma(ims:ime-1,kms:kme-1,:)**2)* &
                (dzdx(ims:ime-1,kms+1:kme,:)+domain%dzdx_u(ims+1:ime,kms:kme-1,:))/mixed_denom(ims:ime-1,kms:kme-1,:)
        I_coef(ims+1:ime,kms:kme-1,:) = &
                (sigma(ims+1:ime,kms:kme-1,:)**2)* &
                (dzdx(ims+1:ime,kms+1:kme,:)+domain%dzdx_u(ims+1:ime,kms:kme-1,:))/mixed_denom(ims+1:ime,kms:kme-1,:)
        J_coef(ims:ime-1,kms+1:kme,:) = &
                (dzdx(ims:ime-1,kms:kme-1,:)+domain%dzdx_u(ims+1:ime,kms+1:kme,:))/mixed_denom(ims:ime-1,kms+1:kme,:)
        K_coef(ims+1:ime,kms+1:kme,:) = &
                -(dzdx(ims+1:ime,kms:kme-1,:)+domain%dzdx_u(ims+1:ime,kms+1:kme,:))/mixed_denom(ims+1:ime,kms+1:kme,:)
        
        L_coef(:,kms:kme-1,jms:jme-1) = &
                -(sigma(:,kms:kme-1,jms:jme-1)**2)* &
                (dzdy(:,kms+1:kme,jms:jme-1)+domain%dzdy_v(:,kms:kme-1,jms+1:jme))/mixed_denom(:,kms:kme-1,jms:jme-1)
        M_coef(:,kms:kme-1,jms+1:jme) = &
                (sigma(:,kms:kme-1,jms+1:jme)**2)* &
                (dzdy(:,kms+1:kme,jms+1:jme)+domain%dzdy_v(:,kms:kme-1,jms+1:jme))/mixed_denom(:,kms:kme-1,jms+1:jme)
        N_coef(:,kms+1:kme,jms:jme-1) = &
                (dzdy(:,kms:kme-1,jms:jme-1)+domain%dzdy_v(:,kms+1:kme,jms+1:jme))/mixed_denom(:,kms+1:kme,jms:jme-1)
        O_coef(:,kms+1:kme,jms+1:jme) = &
                -(dzdy(:,kms:kme-1,jms+1:jme)+domain%dzdy_v(:,kms+1:kme,jms+1:jme))/mixed_denom(:,kms+1:kme,jms+1:jme)

        J_coef(ims:ime-1,kms,:) = (dzdx(ims:ime-1,kms,:)+domain%dzdx_u(ims+1:ime,kms,:))/mixed_denom(ims:ime-1,kms,:)
        K_coef(ims+1:ime,kms,:) = -(dzdx(ims+1:ime,kms,:)+domain%dzdx_u(ims+1:ime,kms,:))/mixed_denom(ims+1:ime,kms,:)
        N_coef(:,kms,jms:jme-1) = (dzdy(:,kms,jms:jme-1)+domain%dzdy_v(:,kms,jms+1:jme))/mixed_denom(:,kms,jms:jme-1)
        O_coef(:,kms,jms+1:jme) = -(dzdy(:,kms,jms+1:jme)+domain%dzdy_v(:,kms,jms+1:jme))/mixed_denom(:,kms,jms+1:jme)

        H_coef(ime,:,:) = H_coef(ime-1,:,:)
        I_coef(ims,:,:) = I_coef(ims+1,:,:)
        J_coef(ime,:,:) = J_coef(ime-1,:,:)
        K_coef(ims,:,:) = K_coef(ims+1,:,:)

        L_coef(:,:,jme) = L_coef(:,:,jme-1)
        M_coef(:,:,jms) = M_coef(:,:,jms+1)
        N_coef(:,:,jme) = N_coef(:,:,jme-1)
        O_coef(:,:,jms) = O_coef(:,:,jms+1)


    end subroutine initialize_coefs
    
    
    !Update the coefs which change with time, i.e. those which depend on alpha
    subroutine update_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain        
        real, allocatable, dimension(:,:,:) :: mixed_denom
        integer ims, ime, kms, kme, jms, jme
        
        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)
        
        allocate(mixed_denom(ims:ime,kms:kme,jms:jme))
        mixed_denom = 2*domain%dx*dz_if(:,kms+1:kme+1,:)*(sigma+sigma**2)


        B_coef(:,kms:kme-1,:) = sigma(:,kms:kme-1,:) * &
                              ( (1./alpha(ims:ime,kms:kme-1,jms:jme)**2 + dzdy(ims:ime,kms:kme-1,jms:jme)**2 + &
                              dzdx(ims:ime,kms:kme-1,jms:jme)**2) * (1./domain%jacobian(ims:ime,kms:kme-1,jms:jme)) + &
                              (1./alpha(ims:ime,kms+1:kme,jms:jme)**2 + dzdy(ims:ime,kms+1:kme,jms:jme)**2 + &
                              dzdx(ims:ime,kms+1:kme,jms:jme)**2) * (1./domain%jacobian(ims:ime,kms+1:kme,jms:jme))) / &
                          ((sigma(:,kms:kme-1,:)+sigma(:,kms:kme-1,:)**2)*dz_if(:,kms+1:kme,:)**2)
                          
                          
        C_coef(:,kms+1:kme,:) = ( (1./alpha(ims:ime,kms:kme-1,jms:jme)**2 + dzdy(ims:ime,kms:kme-1,jms:jme)**2 + &
                              dzdx(ims:ime,kms:kme-1,jms:jme)**2) * (1./domain%jacobian(ims:ime,kms:kme-1,jms:jme)) + &
                              (1./alpha(ims:ime,kms+1:kme,jms:jme)**2 + dzdy(ims:ime,kms+1:kme,jms:jme)**2 + &
                              dzdx(ims:ime,kms+1:kme,jms:jme)**2) * (1./domain%jacobian(ims:ime,kms+1:kme,jms:jme))) / &
                          ((sigma(:,kms+1:kme,:)+sigma(:,kms+1:kme,:)**2)*dz_if(:,kms+2:kme+1,:)**2)
                
        C_coef(:,kms,:) = ( (1./alpha(ims:ime,kms,jms:jme)**2 + dzdy(ims:ime,kms,jms:jme)**2 + &
                              dzdx(ims:ime,kms,jms:jme)**2) * (1./domain%jacobian(ims:ime,kms,jms:jme)) + &
                              (1./alpha(ims:ime,kms,jms:jme)**2 + dzdy(ims:ime,kms,jms:jme)**2 + &
                              dzdx(ims:ime,kms,jms:jme)**2) * (1./domain%jacobian(ims:ime,kms,jms:jme))) / &
                          ((sigma(:,kms,:)+sigma(:,kms,:)**2)*dz_if(:,kms+1,:)**2)
                          
        B_coef(ims+1:ime-1,:,:) = B_coef(ims+1:ime-1,:,:) - sigma(ims+1:ime-1,:,:)**2 * &
                                    (domain%dzdx_u(ims+2:ime,:,:)-domain%dzdx_u(ims+1:ime-1,:,:))/(mixed_denom(ims+1:ime-1,:,:))
                          
        B_coef(:,:,jms+1:jme-1) = B_coef(:,:,jms+1:jme-1) - sigma(:,:,jms+1:jme-1)**2 * &
                                    (domain%dzdy_v(:,:,jms+2:jme)-domain%dzdy_v(:,:,jms+1:jme-1))/(mixed_denom(:,:,jms+1:jme-1))
                          
        
        C_coef(ims+1:ime-1,:,:) = C_coef(ims+1:ime-1,:,:) + &
                                    (domain%dzdx_u(ims+2:ime,:,:)-domain%dzdx_u(ims+1:ime-1,:,:))/(mixed_denom(ims+1:ime-1,:,:))
                          
        C_coef(:,:,jms+1:jme-1) = C_coef(:,:,jms+1:jme-1) + &
                                    (domain%dzdy_v(:,:,jms+2:jme)-domain%dzdy_v(:,:,jms+1:jme-1))/(mixed_denom(:,:,jms+1:jme-1))
                                                    
        B_coef(ims,:,:) = B_coef(ims+1,:,:)
        B_coef(ime,:,:) = B_coef(ime-1,:,:)
        B_coef(:,:,jms) = B_coef(:,:,jms+1)
        B_coef(:,:,jme) = B_coef(:,:,jme-1)
        
        C_coef(ims,:,:) = C_coef(ims+1,:,:)
        C_coef(ime,:,:) = C_coef(ime-1,:,:)
        C_coef(:,:,jms) = C_coef(:,:,jms+1)
        C_coef(:,:,jme) = C_coef(:,:,jme-1)
        
        A_coef(ims+1:ime-1,kms:kme,jms+1:jme-1) = -((domain%jacobian(ims+2:ime,kms:kme,jms+1:jme-1) + &
                                               2*domain%jacobian(ims+1:ime-1,kms:kme,jms+1:jme-1) + &
                                               domain%jacobian(ims:ime-2,kms:kme,jms+1:jme-1))/(2*domain%dx**2)) &
                                            -((domain%jacobian(ims+1:ime-1,kms:kme,jms+2:jme) + &
                                               2*domain%jacobian(ims+1:ime-1,kms:kme,jms+1:jme-1) + &
                                               domain%jacobian(ims+1:ime-1,kms:kme,jms:jme-2))/(2*domain%dx**2)) - &
                                               B_coef(ims+1:ime-1,kms:kme,jms+1:jme-1) - C_coef(ims+1:ime-1,kms:kme,jms+1:jme-1)
                                               
        A_coef(ims,:,jms+1:jme-1) = A_coef(ims+1,:,jms+1:jme-1) 
        A_coef(ime,:,jms+1:jme-1) = A_coef(ime-1,:,jms+1:jme-1) 
        A_coef(ims:ime,:,jms) = A_coef(ims:ime,:,jms+1)
        A_coef(ims:ime,:,jme) = A_coef(ims:ime,:,jme-1)
                    
    end subroutine


end module wind_iterative