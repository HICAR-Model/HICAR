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
    use io_routines, only : io_read, io_write
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
    real, allocatable, dimension(:,:,:) :: div, a_i, b_i, c_i, d_i, e_i, f_i, g_i, h_i, i_i, j_i, k_i, l_i, m_i, n_i
    real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, H_coef, I_coef, &
                                           J_coef, K_coef, L_coef, M_coef, N_coef, O_coef
    integer :: ims, ime, jms, jme, kms, kme
    real    :: dx
    real, allocatable, dimension(:,:,:)  :: dz, jaco, dzdx, dzdy, dzdx2, dzdy2, alpha
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
        real, allocatable  ::  temp_z(:,:,:)
        logical             :: update
        integer k
        
        PetscErrorCode ierr
        KSP            ksp
        PetscReal      norm
        DM             da
        Vec            x
        PetscInt       one
        
        update=.False.
        if (present(update_in)) update=update_in
        
        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)
        
        if (.not.(allocated(dz))) then
            allocate(dz(ims:ime,kms:kme,jms:jme))
            allocate(dzdx(ims:ime,kms:kme,jms:jme))
            allocate(dzdy(ims:ime,kms:kme,jms:jme))
            allocate(dzdx2(ims:ime,kms:kme,jms:jme))
            allocate(dzdy2(ims:ime,kms:kme,jms:jme))
            allocate(jaco(ims:ime,kms:kme,jms:jme))
            allocate(alpha(ims:ime,kms:kme,jms:jme))
            allocate(temp_z(ims:ime,kms:kme,jms:jme))
            
            dx = domain%dx
            dz = domain%advection_dz
            ! Construct dzdx/dzdy using centered differences, on the mass grid. This differes from the domain objects which are 
            ! already staggered to u/v respectively
            
            
            !temp_z(:,1,:) = domain%global_terrain + (domain%advection_dz(1,1,1)/2)*domain%global_jacobian(:,1,:)
        
            !do k=2,kme
            !    temp_z(:,k,:) = temp_z(:,k-1,:) + (((domain%advection_dz(1,k,1)) / 2)*domain%global_jacobian(:,k,:)) + &
            !                                      (((domain%advection_dz(1,k-1,1)) / 2)*domain%global_jacobian(:,k-1,:))
            !enddo
        
            temp_z = domain%z%data_3d
            
            !dzdx2
            dzdx2(ims+1:ime-1,kms:kme,jms:jme) = (temp_z(ims+2:ime,kms:kme,jms:jme) - &
                                                  2*temp_z(ims+1:ime-1,kms:kme,jms:jme) + &
                                                  temp_z(ims:ime-2,kms:kme,jms:jme))/(dx**2)
                                                 
            !dzdx2 = domain%dzdx(ims:ime,kms:kme,jms:jme)**2
            dzdx2(ims,kms:kme,jms:jme) = (-temp_z(ims+3,kms:kme,jms:jme) + 4*temp_z(ims+2,kms:kme,jms:jme) &
                                     -5*temp_z(ims+1,kms:kme,jms:jme) + 2*temp_z(ims,kms:kme,jms:jme)) / (dx**2)

            dzdx2(ime,kms:kme,jms:jme) = (-temp_z(ime-3,kms:kme,jms:jme) + 4*temp_z(ime-2,kms:kme,jms:jme) &
                                     -5*temp_z(ime-1,kms:kme,jms:jme) + 2*temp_z(ime,kms:kme,jms:jme)) / (dx**2)
            
            !dzdy2
            dzdy2(ims:ime,kms:kme,jms+1:jme-1) = (temp_z(ims:ime,kms:kme,jms+2:jme) - &
                                                  2*temp_z(ims:ime,kms:kme,jms+1:jme-1) + &
                                                  temp_z(ims:ime,kms:kme,jms:jme-2))/(dx**2)
                                     
            !dzdy2 = domain%dzdy(ims:ime,kms:kme,jms:jme)**2
            dzdy2(ims:ime,kms:kme,jms) = (-temp_z(ims:ime,kms:kme,jms+3) + 4*temp_z(ims:ime,kms:kme,jms+2) &
                                     -5*temp_z(ims:ime,kms:kme,jms+1) + 2*temp_z(ims:ime,kms:kme,jms)) / (dx**2)

            dzdy2(ims:ime,kms:kme,jme) = (-temp_z(ims:ime,kms:kme,jme-3) + 4*temp_z(ims:ime,kms:kme,jme-2) &
                                     -5*temp_z(ims:ime,kms:kme,jme-1) + 2*temp_z(ims:ime,kms:kme,jme)) / (dx**2)            
            
            
            !dzdx
            dzdx(ims+1:ime-1,kms:kme,jms:jme) = (temp_z(ims+2:ime,kms:kme,jms:jme) - &
                                                 temp_z(ims:ime-2,kms:kme,jms:jme))/(2*dx)
            dzdx(ims,kms:kme,jms:jme) = (-3*temp_z(ims,kms:kme,jms:jme) + &
                                          4*temp_z(ims+1,kms:kme,jms:jme) - temp_z(ims+2,kms:kme,jms:jme)) / (2*dx)
                                          
            dzdx(ime,kms:kme,jms:jme) = (3*temp_z(ime,kms:kme,jms:jme) - &
                                         4*temp_z(ime-1,kms:kme,jms:jme) + temp_z(ime-2,kms:kme,jms:jme)) / (2*dx)
                     
            !dzdy
            dzdy(ims:ime,kms:kme,jms+1:jme-1) = (temp_z(ims:ime,kms:kme,jms+2:jme) - &
                                     temp_z(ims:ime,kms:kme,jms:jme-2))/(2*dx)
            dzdy(ims:ime,kms:kme,jms) = (-3*temp_z(ims:ime,kms:kme,jms) + &
                                          4*temp_z(ims:ime,kms:kme,jms+1) - temp_z(ims:ime,kms:kme,jms+2)) / (2*dx)
            dzdy(ims:ime,kms:kme,jme) = (3*temp_z(ims:ime,kms:kme,jme) - &
                                         4*temp_z(ims:ime,kms:kme,jme-1) + temp_z(ims:ime,kms:kme,jme-2)) / (2*dx)
            
            
            !dzdx = domain%dzdx(ims:ime,kms:kme,jms:jme)
            !dzdy = domain%dzdy(ims:ime,kms:kme,jms:jme)
            
                                    
            jaco = domain%jacobian
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
        write(*,*) 'initialized PETSc'
        
        call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
        call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, &
                          domain%ide+2,domain%kde+2,domain%jde+2,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,one,one, &
                          PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)
        call DMSetFromOptions(da,ierr)
        call DMSetUp(da,ierr)
        
        call KSPSetDM(ksp,da,ierr)
        call KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,0,ierr)
        call KSPSetComputeRHS(ksp,ComputeRHS,0,ierr)
        call KSPSetComputeOperators(ksp,ComputeMatrix,0,ierr)

        call KSPSetFromOptions(ksp,ierr)
        call KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)
        call KSPGetSolution(ksp,x,ierr)
        write(*,*) 'solved KSP'
        !call KSPGetRhs(ksp,b,ierr)
        !call VecDuplicate(b,r,ierr)
        !call KSPGetOperators(ksp,A,NULL,ierr)

        !call MatMult(A,x,r,ierr)
        !call VecAXPY(r,-1.0,b,ierr)
        !call VecNorm(r,NORM_2,&norm,ierr)
        !call PetscPrintf(PETSC_COMM_WORLD,"Residual norm %g\n",(double)norm)

        !call VecDestroy(&r,ierr)
        call DMDAVecGetArrayF90(da,x,lambda, ierr)

        !call reconstruct_Lambda(x,lambda,domain%ide,domain%kde,domain%jde)
        
        if (update) then
            call calc_updated_winds(domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, domain%w%meta_data%dqdt_3d, &
                                    domain%jacobian_u, domain%jacobian_v, domain%jacobian_w, domain%dzdx, domain%dzdy, lambda)
        else
            call calc_updated_winds(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, &
                                    domain%jacobian_u, domain%jacobian_v, domain%jacobian_w, domain%dzdx, domain%dzdy, lambda)
        endif
        
        call DMDAVecRestoreArrayF90(da,x,lambda, ierr)
        call DMDestroy(da,ierr)
        call KSPDestroy(ksp,ierr)
        call PetscFinalize(ierr)

    end subroutine calc_iter_winds
    
    subroutine calc_updated_winds(u, v, w, jaco_u,jaco_v,jaco_w,u_dzdx,v_dzdy,lambda)
        real, intent(inout), dimension(:,:,:)  :: u,v,w
        real, intent(inout), dimension(:,:,:)     :: jaco_u,jaco_v,jaco_w, u_dzdx, v_dzdy
        PetscScalar, intent(in), pointer       :: lambda(:,:,:)

        real, allocatable, dimension(:,:,:)    :: u_dlambdz, v_dlambdz, u_temp, v_temp, dummy_lambda
        integer k
        
        allocate(dummy_lambda(ims-1:ime+1,kms-1:kme+1,jms-1:jme+1))
        dummy_lambda = lambda
        call io_write("lambda.nc","lambda",dummy_lambda)
        
        allocate(u_temp(ims:ime+1,kms-1:kme+1,jms:jme))
        allocate(v_temp(ims:ime,kms-1:kme+1,jms:jme+1))

        allocate(u_dlambdz(ims:ime+1,kms:kme,jms:jme))
        allocate(v_dlambdz(ims:ime,kms:kme,jms:jme+1))
        
        !stager lambda to u grid
        u_temp = (lambda(ims-1:ime,kms-1:kme+1,jms:jme) + lambda(ims:ime+1,kms-1:kme+1,jms:jme)) / 2 
        !stager lambda to z interfaces on u grid
        u_temp(:,kms:kme+1,:) = (u_temp(:,kms-1:kme,:) + u_temp(:,kms:kme+1,:)) / 2 
        !difference
        u_dlambdz = (u_temp(:,kms+1:kme+1,:) - u_temp(:,kms:kme,:))

        !stager lambda to v grid
        v_temp = (lambda(ims:ime,kms-1:kme+1,jms-1:jme) + lambda(ims:ime,kms-1:kme+1,jms:jme+1)) / 2 
        !stager lambda to z interfaces on v grid
        v_temp(:,kms:kme+1,:) = (v_temp(:,kms-1:kme,:) + v_temp(:,kms:kme+1,:)) / 2 
        !difference
        v_dlambdz = (v_temp(:,kms+1:kme+1,:) - v_temp(:,kms:kme,:))
        
        
        !divide dz differennces by dz. Note that dz will be horizontally constant
        do k=kms,kme
            u_dlambdz(:,k,:) = u_dlambdz(:,k,:)/dz(1,k,1)
            v_dlambdz(:,k,:) = v_dlambdz(:,k,:)/dz(1,k,1)
        enddo

        u_dzdx(ims:ime,:,:) = dzdx
        v_dzdy(:,:,jms:jme) = dzdy
        
        !PETSc arrays are zero-indexed
        u = u - (lambda(ims:ime+1,kms:kme,jms:jme)-lambda(ims-1:ime,kms:kme,jms:jme))/dx + &
                (1/jaco_u)*u_dzdx*(u_dlambdz)
        v = v - (lambda(ims:ime,kms:kme,jms:jme+1)-lambda(ims:ime,kms:kme,jms-1:jme))/dx + &
                (1/jaco_v)*v_dzdy*(v_dlambdz)
        !w = w + (lambda(ims:ime,kms-1:kme-1,jms:jme)-lambda(ims:ime,kms:kme,jms:jme))/dz


    end subroutine calc_updated_winds

    subroutine ComputeRHS(ksp,vec_b,dummy,ierr)!,void *ctx)
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
        call KSPGetSolution(ksp,x,ierr)
        
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
        call DMDAVecGetArrayF90(dm,x,lambda, ierr)
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
                    !For global boundary conditions
                    else if (k.eq.0) then
                        if (i.eq.xs) then
                            dlambdx = (-3*lambda(i,k,j) + 4*lambda(i+1,k,j) - lambda(i+2,k,j)) / (2*dx)
                        else if (i.eq.(xs+xm-1)) then
                            dlambdx = (3*lambda(i,k,j) - 4*lambda(i-1,k,j) + lambda(i-2,k,j)) / (2*dx)
                        else
                            dlambdx = (lambda(i+1,k,j) - lambda(i-1,k,j)) / (2*dx)
                        endif
                        
                        if (j.eq.ys) then
                            dlambdy = (-3*lambda(i,k,j) + 4*lambda(i,k,j+1) - lambda(i,k,j+2)) / (2*dx)
                        else if (j.eq.(ys+ym-1)) then
                            dlambdy = (3*lambda(i,k,j) - 4*lambda(i,k,j-1) + lambda(i,k,j-2)) / (2*dx)
                        else
                            dlambdy = (lambda(i,k,j+1) - lambda(i,k,j-1)) / (2*dx)
                        endif
                        
                        barray(i,k,j) = - (dz(HICAR_i,HICAR_k+1,HICAR_j)*jaco(HICAR_i,HICAR_k+1,HICAR_j)* &
                                          ( dzdx(HICAR_i,HICAR_k+1,HICAR_j)*dlambdx + dzdy(HICAR_i,HICAR_k+1,HICAR_j)*dlambdy)) / &
                                          (alpha(HICAR_i,HICAR_k+1,HICAR_j)**2 + dzdx2(HICAR_i,HICAR_k+1,HICAR_j) + &
                                          dzdy2(HICAR_i,HICAR_k+1,HICAR_j))
                    else
                        barray(i,k,j) = -div(HICAR_i,HICAR_k,HICAR_j)
                    endif
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayF90(dm,x,lambda, ierr)
        call DMDAVecRestoreArrayF90(dm,vec_b,barray, ierr)
    end subroutine ComputeRHS

    subroutine ComputeInitialGuess(ksp,vec_b,ctx,ierr)!,void *ctx)
        implicit none
        PetscErrorCode  ierr
        KSP ksp
        PetscInt ctx(*)
        Vec vec_b
        PetscScalar  i_guess

        i_guess = 0.0

        call VecSet(vec_b,i_guess,ierr)
    end subroutine ComputeInitialGuess

    subroutine ComputeMatrix(ksp,arr_A,arr_B,dummy,ierr)!,void *ctx)
        implicit none
        PetscErrorCode  ierr
        KSP ksp
        Mat arr_A,arr_B
        integer dummy(*)


        DM             da
        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,i1,i2,i15
        DMDALocalInfo  info(DMDA_LOCAL_INFO_SIZE)
        PetscScalar    v(15)
        MatStencil     row(4),col(4,15),gnd_col(4,2)
        
        i1 = 1
        i2 = 2
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
                if (i.eq.0 .or. j.eq.0 .or. &
                    i.eq.mx-1 .or. j.eq.my-1 .or. k.eq.mz-1) then
                    v(1) = 1.0
                    call MatSetValuesStencil(arr_B,i1,row,i1,row,v,INSERT_VALUES, ierr)
                else if (k.eq.0) then
                    !center
                    v(1) = 1.0
                    gnd_col(MatStencil_i,1) = i
                    gnd_col(MatStencil_j,1) = k
                    gnd_col(MatStencil_k,1) = j
                    !k + 1
                    v(2) = -1.0
                    gnd_col(MatStencil_i,2) = i
                    gnd_col(MatStencil_j,2) = k+1
                    gnd_col(MatStencil_k,2) = j
                    call MatSetValuesStencil(arr_B,i1,row,i2,gnd_col,v,INSERT_VALUES, ierr)                
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
        
        
        allocate(a_i(ims:ime,kms:kme,jms:jme))
        allocate(b_i(ims:ime,kms:kme,jms:jme))
        allocate(c_i(ims:ime,kms:kme,jms:jme))
        allocate(d_i(ims:ime,kms:kme,jms:jme))
        allocate(e_i(ims:ime,kms:kme,jms:jme))
        allocate(f_i(ims:ime,kms:kme,jms:jme))
        allocate(g_i(ims:ime,kms:kme,jms:jme))
        allocate(h_i(ims:ime,kms:kme,jms:jme))
        allocate(i_i(ims:ime,kms:kme,jms:jme))
        allocate(j_i(ims:ime,kms:kme,jms:jme))
        allocate(k_i(ims:ime,kms:kme,jms:jme))
        allocate(l_i(ims:ime,kms:kme,jms:jme))
        allocate(m_i(ims:ime,kms:kme,jms:jme))
        allocate(n_i(ims:ime,kms:kme,jms:jme))
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
        
        a_i = 0
        b_i = 0
        c_i = 0
        d_i = 0
        e_i = 0
        f_i = 0
        g_i = 0
        h_i = 0
        i_i = 0
        j_i = 0
        k_i = 0
        l_i = 0
        m_i = 0
        n_i = 0
        A_coef = 1
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
        
        a_i(ims+1:ime,kms:kme,jms:jme) = dzdx(ims+1:ime,kms:kme,jms:jme) + dzdx(ims:ime-1,kms:kme,jms:jme)
        b_i(ims:ime-1,kms:kme,jms:jme) = dzdx(ims:ime-1,kms:kme,jms:jme) + dzdx(ims+1:ime,kms:kme,jms:jme)
        c_i(ims+1:ime-1,kms:kme,jms:jme) = dzdx(ims+2:ime,kms:kme,jms:jme) - dzdx(ims:ime-2,kms:kme,jms:jme)
        d_i(ims:ime,kms:kme,jms+1:jme) = dzdy(ims:ime,kms:kme,jms+1:jme) + dzdy(ims:ime,kms:kme,jms:jme-1)
        e_i(ims:ime,kms:kme,jms:jme-1) = dzdy(ims:ime,kms:kme,jms:jme-1) + dzdy(ims:ime,kms:kme,jms+1:jme)
        f_i(ims:ime,kms:kme,jms+1:jme-1) = dzdy(ims:ime,kms:kme,jms+2:jme) - dzdy(ims:ime,kms:kme,jms:jme-2)
        g_i(ims:ime,kms+1:kme,jms:jme) = dzdx(ims:ime,kms+1:kme,jms:jme) + dzdx(ims:ime,kms:kme-1,jms:jme)
        h_i(ims:ime,kms:kme-1,jms:jme) = dzdx(ims:ime,kms:kme-1,jms:jme) + dzdx(ims:ime,kms+1:kme,jms:jme)
        i_i(ims:ime,kms+1:kme-1,jms:jme) = dzdx(ims:ime,kms+2:kme,jms:jme) - dzdx(ims:ime,kms:kme-2,jms:jme)
        j_i(ims:ime,kms+1:kme,jms:jme) = dzdy(ims:ime,kms+1:kme,jms:jme) + dzdy(ims:ime,kms:kme-1,jms:jme)
        k_i(ims:ime,kms:kme-1,jms:jme) = dzdy(ims:ime,kms:kme-1,jms:jme) + dzdy(ims:ime,kms+1:kme,jms:jme)
        l_i(ims:ime,kms+1:kme-1,jms:jme) = dzdy(ims:ime,kms+2:kme,jms:jme) - dzdy(ims:ime,kms:kme-2,jms:jme)
        
        m_i(:,kms:kme-1,:) = (alpha(ims:ime,kms:kme-1,jms:jme)**2 + alpha(ims:ime,kms+1:kme,jms:jme)**2 + &
                          dzdy2(ims:ime,kms:kme-1,jms:jme) + dzdy2(ims:ime,kms+1:kme,jms:jme) + &
                          dzdx2(ims:ime,kms:kme-1,jms:jme) + dzdx2(ims:ime,kms+1:kme,jms:jme)) / &
                          ((domain%jacobian(ims:ime,kms:kme-1,jms:jme) + domain%jacobian(ims:ime,kms+1:kme,jms:jme)) * 4 * &
                          domain%advection_dz(ims:ime,kms:kme-1,jms:jme) * &
                          (domain%advection_dz(ims:ime,kms+1:kme,jms:jme)/2 + domain%advection_dz(ims:ime,kms:kme-1,jms:jme)/2 ))
        n_i(:,kms+1:kme,:) = (alpha(ims:ime,kms:kme-1,jms:jme)**2 + alpha(ims:ime,kms+1:kme,jms:jme)**2 + &
                          dzdy2(ims:ime,kms:kme-1,jms:jme) + dzdy2(ims:ime,kms+1:kme,jms:jme) + &
                          dzdx2(ims:ime,kms:kme-1,jms:jme) + dzdx2(ims:ime,kms+1:kme,jms:jme)) / &
                          ((domain%jacobian(ims:ime,kms:kme-1,jms:jme) + domain%jacobian(ims:ime,kms+1:kme,jms:jme)) * 4 * &
                          domain%advection_dz(ims:ime,kms+1:kme,jms:jme) * &
                          (domain%advection_dz(ims:ime,kms+1:kme,jms:jme)/2 + domain%advection_dz(ims:ime,kms:kme-1,jms:jme)/2 ))
        
        A_coef(ims+1:ime-1,kms:kme,jms+1:jme-1) = -((domain%jacobian(ims+2:ime,kms:kme,jms+1:jme-1) + &
                                               2*domain%jacobian(ims+1:ime-1,kms:kme,jms+1:jme-1) + &
                                               domain%jacobian(ims:ime-2,kms:kme,jms+1:jme-1))/(2*domain%dx**2)) &
                                            -((domain%jacobian(ims+1:ime-1,kms:kme,jms+2:jme) + &
                                               2*domain%jacobian(ims+1:ime-1,kms:kme,jms+1:jme-1) + &
                                               domain%jacobian(ims+1:ime-1,kms:kme,jms:jme-2))/(2*domain%dx**2)) - &
                                               m_i(ims+1:ime-1,kms:kme,jms+1:jme-1) - n_i(ims+1:ime-1,kms:kme,jms+1:jme-1)
        B_coef = m_i - ((c_i+f_i)/(8*domain%dx*domain%advection_dz))
        C_coef = n_i + ((c_i+f_i)/(8*domain%dx*domain%advection_dz))
        
        D_coef(ims:ime-1,kms:kme,jms:jme) = (domain%jacobian(ims+1:ime,kms:kme,jms:jme) + &
            domain%jacobian(ims:ime-1,kms:kme,jms:jme))/(2*domain%dx**2) - &
            i_i(ims:ime-1,kms:kme,jms:jme)/(8*domain%dx*domain%advection_dz(ims:ime-1,kms:kme,jms:jme))
        E_coef(ims+1:ime,kms:kme,jms:jme) = (domain%jacobian(ims+1:ime,kms:kme,jms:jme) + &
            domain%jacobian(ims:ime-1,kms:kme,jms:jme))/(2*domain%dx**2) + &
            i_i(ims+1:ime,kms:kme,jms:jme)/(8*domain%dx*domain%advection_dz(ims+1:ime,kms:kme,jms:jme))
        F_coef(ims:ime,kms:kme,jms:jme-1) = (domain%jacobian(ims:ime,kms:kme,jms+1:jme) + &
            domain%jacobian(ims:ime,kms:kme,jms:jme-1))/(2*domain%dx**2) - &
            l_i(ims:ime,kms:kme,jms:jme-1)/(8*domain%dx*domain%advection_dz(ims:ime,kms:kme,jms:jme-1))
        G_coef(ims:ime,kms:kme,jms+1:jme) = (domain%jacobian(ims:ime,kms:kme,jms+1:jme) + &
            domain%jacobian(ims:ime,kms:kme,jms:jme-1))/(2*domain%dx**2) + &
            l_i(ims:ime,kms:kme,jms+1:jme)/(8*domain%dx*domain%advection_dz(ims:ime,kms:kme,jms+1:jme))
            
        H_coef = -(b_i+h_i)/(8*domain%dx*domain%advection_dz)
        I_coef = (a_i+h_i)/(8*domain%dx*domain%advection_dz)
        J_coef = (b_i+g_i)/(8*domain%dx*domain%advection_dz)
        K_coef = -(a_i+g_i)/(8*domain%dx*domain%advection_dz)
        L_coef = -(e_i+k_i)/(8*domain%dx*domain%advection_dz)
        M_coef = (d_i+k_i)/(8*domain%dx*domain%advection_dz)
        N_coef = (e_i+j_i)/(8*domain%dx*domain%advection_dz)
        O_coef = -(d_i+j_i)/(8*domain%dx*domain%advection_dz)
        
    end subroutine initialize_coefs
    
    !Update the coefs which change with time, i.e. those which depend on alpha
    subroutine update_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        
        m_i(:,kms:kme-1,:) = (alpha(ims:ime,kms:kme-1,jms:jme)**2 + alpha(ims:ime,kms+1:kme,jms:jme)**2 + &
                          dzdy2(ims:ime,kms:kme-1,jms:jme) + dzdy2(ims:ime,kms+1:kme,jms:jme) + &
                          dzdx2(ims:ime,kms:kme-1,jms:jme) + dzdx2(ims:ime,kms+1:kme,jms:jme)) / &
                          ((domain%jacobian(ims:ime,kms:kme-1,jms:jme) + domain%jacobian(ims:ime,kms+1:kme,jms:jme)) * 4 * &
                          domain%advection_dz(ims:ime,kms:kme-1,jms:jme) * &
                          (domain%advection_dz(ims:ime,kms+1:kme,jms:jme)/2 + domain%advection_dz(ims:ime,kms:kme-1,jms:jme)/2 ))
        n_i(:,kms+1:kme,:) = (alpha(ims:ime,kms:kme-1,jms:jme)**2 + alpha(ims:ime,kms+1:kme,jms:jme)**2 + &
                          dzdy2(ims:ime,kms:kme-1,jms:jme) + dzdy2(ims:ime,kms+1:kme,jms:jme) + &
                          dzdx2(ims:ime,kms:kme-1,jms:jme) + dzdx2(ims:ime,kms+1:kme,jms:jme)) / &
                          ((domain%jacobian(ims:ime,kms:kme-1,jms:jme) + domain%jacobian(ims:ime,kms+1:kme,jms:jme)) * 4 * &
                          domain%advection_dz(ims:ime,kms+1:kme,jms:jme) * &
                          (domain%advection_dz(ims:ime,kms+1:kme,jms:jme)/2 + domain%advection_dz(ims:ime,kms:kme-1,jms:jme)/2 ))
        
        A_coef(ims+1:ime-1,kms:kme,jms+1:jme-1) = -((domain%jacobian(ims+2:ime,kms:kme,jms+1:jme-1) + &
                                               2*domain%jacobian(ims+1:ime-1,kms:kme,jms+1:jme-1) + &
                                               domain%jacobian(ims:ime-2,kms:kme,jms+1:jme-1))/(2*domain%dx**2)) &
                                            -((domain%jacobian(ims+1:ime-1,kms:kme,jms+2:jme) + &
                                               2*domain%jacobian(ims+1:ime-1,kms:kme,jms+1:jme-1) + &
                                               domain%jacobian(ims+1:ime-1,kms:kme,jms:jme-2))/(2*domain%dx**2)) - &
                                               m_i(ims+1:ime-1,kms:kme,jms+1:jme-1) - n_i(ims+1:ime-1,kms:kme,jms+1:jme-1)
        B_coef = m_i - ((c_i+f_i)/(8*domain%dx*domain%advection_dz))
        C_coef = n_i + ((c_i+f_i)/(8*domain%dx*domain%advection_dz))
                    
    end subroutine

end module wind_iterative