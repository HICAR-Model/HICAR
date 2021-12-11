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
    real, allocatable, dimension(:,:)   :: bot_winds
    real, allocatable, dimension(:,:,:) :: div, sigma, dz_if
    real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, H_coef, I_coef, &
                                           J_coef, K_coef, L_coef, M_coef, N_coef, O_coef
    integer :: ims, ime, jms, jme, kms, kme
    real    :: dx
    real, allocatable, dimension(:,:,:)  :: dz, jaco, dzdx, dzdy, dzdx2, dzdy2, d2zdx2, d2zdy2, alpha
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
        PetscReal      norm, conv_tol
        DM             da
        Vec            x
        PetscInt       one, x_size
        
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
            allocate(d2zdx2(ims:ime,kms:kme,jms:jme))
            allocate(d2zdy2(ims:ime,kms:kme,jms:jme))
            allocate(jaco(ims:ime,kms:kme,jms:jme))
            allocate(dz_if(ims:ime,kms:kme+1,jms:jme))
            allocate(sigma(ims:ime,kms:kme,jms:jme))
            allocate(alpha(ims:ime,kms:kme,jms:jme))
            allocate(temp_z(ims:ime,kms:kme,jms:jme))
            allocate(div(ims:ime,kms:kme,jms:jme))
            allocate(bot_winds(ims:ime,jms:jme))

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
            d2zdx2(ims+1:ime-1,kms:kme,jms:jme) = (temp_z(ims+2:ime,kms:kme,jms:jme) - &
                                                  2*temp_z(ims+1:ime-1,kms:kme,jms:jme) + &
                                                  temp_z(ims:ime-2,kms:kme,jms:jme))/(dx**2)
                                                 
            !dzdx2 = domain%dzdx(ims:ime,kms:kme,jms:jme)**2
            d2zdx2(ims,kms:kme,jms:jme) = (-temp_z(ims+3,kms:kme,jms:jme) + 4*temp_z(ims+2,kms:kme,jms:jme) &
                                     -5*temp_z(ims+1,kms:kme,jms:jme) + 2*temp_z(ims,kms:kme,jms:jme)) / (dx**2)

            d2zdx2(ime,kms:kme,jms:jme) = (-temp_z(ime-3,kms:kme,jms:jme) + 4*temp_z(ime-2,kms:kme,jms:jme) &
                                     -5*temp_z(ime-1,kms:kme,jms:jme) + 2*temp_z(ime,kms:kme,jms:jme)) / (dx**2)
            
            !dzdy2
            d2zdy2(ims:ime,kms:kme,jms+1:jme-1) = (temp_z(ims:ime,kms:kme,jms+2:jme) - &
                                                  2*temp_z(ims:ime,kms:kme,jms+1:jme-1) + &
                                                  temp_z(ims:ime,kms:kme,jms:jme-2))/(dx**2)
                                     
            !dzdy2 = domain%dzdy(ims:ime,kms:kme,jms:jme)**2
            d2zdy2(ims:ime,kms:kme,jms) = (-temp_z(ims:ime,kms:kme,jms+3) + 4*temp_z(ims:ime,kms:kme,jms+2) &
                                     -5*temp_z(ims:ime,kms:kme,jms+1) + 2*temp_z(ims:ime,kms:kme,jms)) / (dx**2)

            d2zdy2(ims:ime,kms:kme,jme) = (-temp_z(ims:ime,kms:kme,jme-3) + 4*temp_z(ims:ime,kms:kme,jme-2) &
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
            
            dzdx = domain%dzdx !(domain%dzdx(ims+1:ime+1,kms:kme,jms:jme) + domain%dzdx(ims:ime,kms:kme,jms:jme))*0.5
            dzdy = domain%dzdy !(domain%dzdy(ims:ime,kms:kme,jms+1:jme+1) + domain%dzdy(ims:ime,kms:kme,jms:jme))*0.5
            dzdx2 = dzdx**2
            dzdy2 = dzdy**2
                                    
            jaco = domain%jacobian
        endif
                
                
        !Initialize div to be the initial divergence of the input wind field
        div=div_in
        one = 1
        
        alpha = alpha_in
        if (update) then
            bot_winds = domain%jacobian(:,kms,:)**2 * &
                        (domain%u%meta_data%dqdt_3d(ims+1:ime+1,kms,:) + domain%u%meta_data%dqdt_3d(ims:ime,kms,:)) * 0 * dzdx(:,kms,:) + &
                        (domain%v%meta_data%dqdt_3d(:,kms,jms+1:jme+1) + domain%v%meta_data%dqdt_3d(:,kms,jms:jme)) * 0 * dzdy(:,kms,:) - domain%w%meta_data%dqdt_3d(:,kms,:)
        else
            bot_winds = (1./domain%jacobian(:,kms,:)) * &
                        ((domain%u%data_3d(ims+1:ime+1,kms,:)+domain%u%data_3d(ims:ime,kms,:))*0*dzdx(:,kms,:) + &
                        (domain%v%data_3d(:,kms,jms+1:jme+1)+domain%v%data_3d(:,kms,jms:jme))*0*dzdy(:,kms,:) - &
                        domain%w%data_3d(:,kms,:))
        endif
        
        call io_write('jacobian.nc',"jacobian",domain%jacobian)
        call io_write("div.nc","div",div)
        
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
        conv_tol = 1e-4
        call KSPSetTolerances(ksp,conv_tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        !call KSPSetType(ksp,KSPBCGS,ierr);
        call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DMDA_STENCIL_BOX, &
                          domain%ide+2,domain%kde+2,domain%jde+2,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,one,one, &
                          PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)
                          

        !call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, &
        !                  domain%ide+2,domain%kde+2,domain%jde+2,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,one,one, &
        !                  PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)
        call DMSetFromOptions(da,ierr)
        call DMSetUp(da,ierr)
        
        call KSPSetDM(ksp,da,ierr)
        call KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,0,ierr)
        call KSPSetComputeRHS(ksp,ComputeRHS,0,ierr)
        call KSPSetComputeOperators(ksp,ComputeMatrix,0,ierr)

        call KSPSetFromOptions(ksp,ierr)
        call DMCreateGlobalVector(da,x,ierr)
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

        write(*,*) 'lambda is: ', lbound(lambda,1)
        write(*,*) 'lambda ie: ', ubound(lambda,1)
        write(*,*) 'lambda ks: ', lbound(lambda,2)
        write(*,*) 'lambda ke: ', ubound(lambda,2)
        write(*,*) 'lambda js: ', lbound(lambda,3)
        write(*,*) 'lambda je: ', ubound(lambda,3)
        
        !call reconstruct_Lambda(x,lambda,domain%ide,domain%kde,domain%jde)
        
        if (update) then
            call calc_updated_winds(domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, domain%w%meta_data%dqdt_3d, &
                                  domain%jacobian_u, domain%jacobian_v, domain%jacobian_w, domain%dzdx_u, domain%dzdy_v, lambda)
        else
            call calc_updated_winds(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, &
                                  domain%jacobian_u, domain%jacobian_v, domain%jacobian_w, domain%dzdx_u, domain%dzdy_v, lambda)
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
        
        
        !l_ims = lbound(lambda,1)
        !l_ime = ubound(lambda,1)
        !l_kms = lbound(lambda,2)
        !l_kme = ubound(lambda,2)
        !l_jms = lbound(lambda,3)
        !l_jme = ubound(lambda,3)
        
        write(*,*) 'lambda is: ', lbound(lambda,1)
        write(*,*) 'lambda ie: ', ubound(lambda,1)
        write(*,*) 'lambda ks: ', lbound(lambda,2)
        write(*,*) 'lambda ke: ', ubound(lambda,2)
        write(*,*) 'lambda js: ', lbound(lambda,3)
        write(*,*) 'lambda je: ', ubound(lambda,3)
        
        allocate(dummy_lambda(ims:ime,kms-1:kme+1,jms-1:jme+1))
        dummy_lambda = lambda!(ims-1:ime+1,kms-1:kme+1,jms-1:jme+1)
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
        
                
        !stager lambda to u grid
        u_temp = (lambda(ims-1:ime,kms-1:kme+1,jms:jme) + lambda(ims:ime+1,kms-1:kme+1,jms:jme)) / 2 

        !stager lambda to v grid
        v_temp = (lambda(ims:ime,kms-1:kme+1,jms-1:jme) + lambda(ims:ime,kms-1:kme+1,jms:jme+1)) / 2 

        !divide dz differennces by dz. Note that dz will be horizontally constant
        do k=kms,kme
            u_dlambdz(:,k,:) = (u_temp(:,k+1,:)*sigma(1,k,1)**2) - (sigma(1,k,1)**2 - 1)*u_temp(:,k,:) - u_temp(:,k-1,:)
            v_dlambdz(:,k,:) = (v_temp(:,k+1,:)*sigma(1,k,1)**2) - (sigma(1,k,1)**2 - 1)*v_temp(:,k,:) - v_temp(:,k-1,:)
        
            u_dlambdz(:,k,:) = u_dlambdz(:,k,:)/(dz_if(1,k+1,1)*(sigma(1,k,1)+sigma(1,k,1)**2))
            v_dlambdz(:,k,:) = v_dlambdz(:,k,:)/(dz_if(1,k+1,1)*(sigma(1,k,1)+sigma(1,k,1)**2))
        enddo
        
        call io_write("u_dlambdz_raw.nc","u_dlambdz",(u_dlambdz))
        call io_write("v_dlambdz_raw.nc","v_dlambdz",(v_dlambdz))
        
        call io_write("u_dlambdz.nc","u_dlambdz",0.5*-(1/jaco_u)*u_dzdx*(u_dlambdz))
        call io_write("v_dlambdz.nc","v_dlambdz",0.5*-(1/jaco_v)*v_dzdy*(v_dlambdz))

        call io_write("dzdx.nc","dzdx",dzdx)
        call io_write("dzdx_u.nc","dzdx_u",u_dzdx)

        !dlambdx=dlambdx/dx
        !dlambdy=dlambdy/dx

        !u_dzdx(ims:ime,:,:) = dzdx
        !v_dzdy(:,:,jms:jme) = dzdy
        
        
        !PETSc arrays are zero-indexed
        u = u + 0.5*((lambda(ims:ime+1,kms:kme,jms:jme)-lambda(ims-1:ime,kms:kme,jms:jme))/dx - &
                (1/jaco_u)*u_dzdx*(u_dlambdz))
        v = v + 0.5*((lambda(ims:ime,kms:kme,jms:jme+1)-lambda(ims:ime,kms:kme,jms-1:jme))/dx - &
                (1/jaco_v)*v_dzdy*(v_dlambdz))
                
        !u(ims:ime,:,:) = u(ims:ime,:,:) - (1/jaco)*dzdx(ims:ime,:,:)*(u_dlambdz(ims:ime,:,:))
        !v(:,:,jms:jme) = v(:,:,jms:jme) - (1/jaco)*dzdy(:,:,jms:jme)*(v_dlambdz(:,:,jms:jme))
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
        
        write(*,*) 'image #: ',this_image()
        write(*,*) 'ims: ',ims
        write(*,*) 'ime: ',ime
        write(*,*) 'kms: ',kms
        write(*,*) 'kme: ',kme
        write(*,*) 'jms: ',jms
        write(*,*) 'jme: ',jme
        write(*,*) 'xs: ',xs
        write(*,*) 'xm: ',(xs+xm-1)
        write(*,*) 'zs: ',zs
        write(*,*) 'zm: ',(zs+zm-1)
        write(*,*) 'ys: ',ys
        write(*,*) 'ym: ',(ys+ym-1)

        
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
                        barray(i,k,j) = -0.0*bot_winds(i,j)*dz_if(i,1,j)/(1./alpha(i,1,j)**2 + dzdx2(i,1,j) + dzdy2(i,1,j))
                    else
                        barray(i,k,j) = -2*div(HICAR_i,HICAR_k,HICAR_j)!/jaco(HICAR_i,HICAR_k,HICAR_j)
                    endif
                enddo
            enddo
        enddo
        write(*,*) 'min of jack: ',minval(jaco)
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
        real denom


        DM             da
        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,i1,i2,i6,i15
        DMDALocalInfo  info(DMDA_LOCAL_INFO_SIZE)
        PetscScalar    v(15)
        MatStencil     row(4),col(4,15),gnd_col(4,6),top_col(4,2)
        
        i1 = 1
        i2 = 2
        i6 = 6
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
                else if (k.eq.mz-1) then
                    !k
                    v(1) = 1.0
                    top_col(MatStencil_i,1) = i
                    top_col(MatStencil_j,1) = k
                    top_col(MatStencil_k,1) = j
                    !k - 1
                    v(2) = -1.0
                    top_col(MatStencil_i,2) = i
                    top_col(MatStencil_j,2) = k-1
                    top_col(MatStencil_k,2) = j        
                    call MatSetValuesStencil(arr_B,i1,row,i2,top_col,v,INSERT_VALUES, ierr)
                !else if (k.eq.1) then
                ! 
                !    denom = (alpha(i,1,j)**2 + dzdx2(i,1,j) + &
                !                          dzdy2(i,1,j))
                !    !k - 1
                !    v(1) = 0!- (2*sigma(i,2,j) + 1.0)/(sigma(i,2,j) + 1.0)
                !    gnd_col(MatStencil_i,1) = i
                !    gnd_col(MatStencil_j,1) = k-1
                !    gnd_col(MatStencil_k,1) = j
                !    !k
                !    v(2) = -(sigma(i,1,j)**2 - 1)/(sigma(i,1,j) + sigma(i,1,j)**2)
                !    gnd_col(MatStencil_i,2) = i
                !    gnd_col(MatStencil_j,2) = k
                !    gnd_col(MatStencil_k,2) = j
                !    !k + 1
                !    v(3) = sigma(i,1,j)**2/(sigma(i,1,j) + sigma(i,1,j)**2)
                !    gnd_col(MatStencil_i,3) = i
                !    gnd_col(MatStencil_j,3) = k+1
                !    gnd_col(MatStencil_k,3) = j      
                !    !i - 1
                !    v(4) = dz_if(i,2,j)*jaco(i,1,j)*dzdx(i,1,j)/(denom*2*dx)
                !    gnd_col(MatStencil_i,4) = i-1
                !    gnd_col(MatStencil_j,4) = k
                !    gnd_col(MatStencil_k,4) = j
                !    !i + 1
                !    v(5) = -dz_if(i,2,j)*jaco(i,1,j)*dzdx(i,1,j)/(denom*2*dx)
                !    gnd_col(MatStencil_i,5) = i+1
                !    gnd_col(MatStencil_j,5) = k
                !    gnd_col(MatStencil_k,5) = j
                !    !j - 1
                !    v(6) = dz_if(i,2,j)*jaco(i,1,j)*dzdy(i,1,j)/(denom*2*dx)
                !    gnd_col(MatStencil_i,6) = i
                !    gnd_col(MatStencil_j,6) = k
                !    gnd_col(MatStencil_k,6) = j-1
                !    !j + 1
                !    v(7) = -dz_if(i,2,j)*jaco(i,1,j)*dzdy(i,1,j)/(denom*2*dx)
                !    gnd_col(MatStencil_i,7) = i
                !    gnd_col(MatStencil_j,7) = k
                !    gnd_col(MatStencil_k,7) = j+1
                !    call MatSetValuesStencil(arr_B,i1,row,i7,gnd_col,v,INSERT_VALUES, ierr)  
                else if (k.eq.0) then
                                
                    denom = (1./alpha(i,1,j)**2 + dzdx2(i,1,j) + &
                                          dzdy2(i,1,j))/(jaco(i,1,j))
                    !k
                    v(1) = - 1!(2*sigma(i,2,j) + 1.0)/(sigma(i,2,j) + 1.0)
                    gnd_col(MatStencil_i,1) = i
                    gnd_col(MatStencil_j,1) = k
                    gnd_col(MatStencil_k,1) = j
                    !k + 1
                    v(2) = 1!(sigma(i,2,j) + 1.0)
                    gnd_col(MatStencil_i,2) = i
                    gnd_col(MatStencil_j,2) = k+1
                    gnd_col(MatStencil_k,2) = j
                    !i - 1
                    v(3) = dz_if(i,1,j)*dzdx(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,3) = i-1
                    gnd_col(MatStencil_j,3) = k
                    gnd_col(MatStencil_k,3) = j
                    !i + 1
                    v(4) = -dz_if(i,1,j)*dzdx(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,4) = i+1
                    gnd_col(MatStencil_j,4) = k
                    gnd_col(MatStencil_k,4) = j
                    !j - 1
                    v(5) = dz_if(i,1,j)*dzdy(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,5) = i
                    gnd_col(MatStencil_j,5) = k
                    gnd_col(MatStencil_k,5) = j-1
                    !j + 1
                    v(6) = -dz_if(i,1,j)*dzdy(i,1,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,6) = i
                    gnd_col(MatStencil_j,6) = k
                    gnd_col(MatStencil_k,6) = j+1
                    call MatSetValuesStencil(arr_B,i1,row,i6,gnd_col,v,INSERT_VALUES, ierr)
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

        dz_if(:,kms,:) = domain%advection_dz(:,kms,:)
        dz_if(:,kms+1:kme,:) = (domain%advection_dz(:,kms+1:kme,:)/2) + (domain%advection_dz(:,kms:kme-1,:)/2)
        dz_if(:,kme+1,:) = domain%advection_dz(:,kme,:)
        sigma = dz_if(:,kms:kme,:)/dz_if(:,kms+1:kme+1,:)
 
        mixed_denom = 2*domain%dx*dz_if(:,kms+1:kme+1,:)*(sigma+sigma**2)

        B_coef(:,kms:kme-1,:) = sigma(:,kms:kme-1,:) * &
                              ( (1./alpha(ims:ime,kms:kme-1,jms:jme)**2 + dzdy2(ims:ime,kms:kme-1,jms:jme) + &
                              dzdx2(ims:ime,kms:kme-1,jms:jme)) * (1./domain%jacobian(ims:ime,kms:kme-1,jms:jme)) + &
                              (1./alpha(ims:ime,kms+1:kme,jms:jme)**2 + dzdy2(ims:ime,kms+1:kme,jms:jme) + &
                              dzdx2(ims:ime,kms+1:kme,jms:jme)) * (1./domain%jacobian(ims:ime,kms+1:kme,jms:jme))) / &
                          (2*(sigma(:,kms:kme-1,:)+sigma(:,kms:kme-1,:)**2)*dz_if(:,kms+1:kme,:)**2)
                          
                          
        C_coef(:,kms+1:kme,:) = ( (1./alpha(ims:ime,kms:kme-1,jms:jme)**2 + dzdy2(ims:ime,kms:kme-1,jms:jme) + &
                              dzdx2(ims:ime,kms:kme-1,jms:jme)) * (1./domain%jacobian(ims:ime,kms:kme-1,jms:jme)) + &
                              (1./alpha(ims:ime,kms+1:kme,jms:jme)**2 + dzdy2(ims:ime,kms+1:kme,jms:jme) + &
                              dzdx2(ims:ime,kms+1:kme,jms:jme)) * (1./domain%jacobian(ims:ime,kms+1:kme,jms:jme))) / &
                          (2*(sigma(:,kms+1:kme,:)+sigma(:,kms+1:kme,:)**2)*dz_if(:,kms+2:kme+1,:)**2)
                
        C_coef(:,kms,:) = ( (1./alpha(ims:ime,kms,jms:jme)**2 + dzdy2(ims:ime,kms,jms:jme) + &
                              dzdx2(ims:ime,kms,jms:jme)) * (1./domain%jacobian(ims:ime,kms,jms:jme)) + &
                              (1./alpha(ims:ime,kms,jms:jme)**2 + dzdy2(ims:ime,kms,jms:jme) + &
                              dzdx2(ims:ime,kms,jms:jme)) * (1./domain%jacobian(ims:ime,kms,jms:jme))) / &
                          (2*(sigma(:,kms,:)+sigma(:,kms,:)**2)*dz_if(:,kms+1,:)**2)
                          
                          
        !B_coef = B_coef - ( d2zdx2 + d2zdy2 - 2*(dzdx2+dzdy2)/domain%jacobian) * sigma**2 / &
        !                ((sigma+sigma**2)*dz_if(:,kms+1:kme+1,:))
                        
        !C_coef = C_coef + ( d2zdx2 + d2zdy2 - 2*(dzdx2+dzdy2)/domain%jacobian) / &
        !                ((sigma+sigma**2)*dz_if(:,kms+1:kme+1,:))
                        
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

        D_coef(ims:ime-1,kms:kme,jms:jme) = (domain%jacobian(ims+1:ime,kms:kme,jms:jme) + &
            domain%jacobian(ims:ime-1,kms:kme,jms:jme))/(2*domain%dx**2) + &
            (sigma(ims:ime-1,:,:)**2 - 1)*(dzdx(ims+1:ime,:,:)+dzdx(ims:ime-1,:,:))/mixed_denom(ims:ime-1,kms:kme,jms:jme)
        E_coef(ims+1:ime,kms:kme,jms:jme) = (domain%jacobian(ims+1:ime,kms:kme,jms:jme) + &
            domain%jacobian(ims:ime-1,kms:kme,jms:jme))/(2*domain%dx**2) - &
            (sigma(ims+1:ime,:,:)**2 - 1)*(dzdx(ims+1:ime,:,:)+dzdx(ims:ime-1,:,:))/mixed_denom(ims+1:ime,kms:kme,jms:jme)
        F_coef(ims:ime,kms:kme,jms:jme-1) = (domain%jacobian(ims:ime,kms:kme,jms+1:jme) + &
            domain%jacobian(ims:ime,kms:kme,jms:jme-1))/(2*domain%dx**2) + &
            (sigma(:,:,jms:jme-1)**2 - 1)*(dzdy(:,:,jms+1:jme)+dzdy(:,:,jms:jme-1))/mixed_denom(ims:ime,kms:kme,jms:jme-1)
        G_coef(ims:ime,kms:kme,jms+1:jme) = (domain%jacobian(ims:ime,kms:kme,jms+1:jme) + &
            domain%jacobian(ims:ime,kms:kme,jms:jme-1))/(2*domain%dx**2) - &
            (sigma(:,:,jms+1:jme)**2 - 1)*(dzdy(:,:,jms+1:jme)+dzdy(:,:,jms:jme-1))/mixed_denom(ims:ime,kms:kme,jms+1:jme)
            
       
        H_coef(ims:ime-1,kms:kme-1,:) = &
                -(sigma(ims:ime-1,kms:kme-1,:)**2)* &
                (dzdx(ims:ime-1,kms+1:kme,:)+dzdx(ims+1:ime,kms:kme-1,:))/mixed_denom(ims:ime-1,kms:kme-1,:)
        I_coef(ims+1:ime,kms:kme-1,:) = &
                (sigma(ims+1:ime,kms:kme-1,:)**2)* &
                (dzdx(ims+1:ime,kms+1:kme,:)+dzdx(ims:ime-1,kms:kme-1,:))/mixed_denom(ims+1:ime,kms:kme-1,:)
        J_coef(ims:ime-1,kms+1:kme,:) = &
                (dzdx(ims:ime-1,kms:kme-1,:)+dzdx(ims+1:ime,kms+1:kme,:))/mixed_denom(ims:ime-1,kms+1:kme,:)
        K_coef(ims+1:ime,kms+1:kme,:) = &
                -(dzdx(ims+1:ime,kms:kme-1,:)+dzdx(ims:ime-1,kms+1:kme,:))/mixed_denom(ims+1:ime,kms+1:kme,:)
        
        L_coef(:,kms:kme-1,jms:jme-1) = &
                -(sigma(:,kms:kme-1,jms:jme-1)**2)* &
                (dzdy(:,kms+1:kme,jms:jme-1)+dzdy(:,kms:kme-1,jms+1:jme))/mixed_denom(:,kms:kme-1,jms:jme-1)
        M_coef(:,kms:kme-1,jms+1:jme) = &
                (sigma(:,kms:kme-1,jms+1:jme)**2)* &
                (dzdy(:,kms+1:kme,jms+1:jme)+dzdy(:,kms:kme-1,jms:jme-1))/mixed_denom(:,kms:kme-1,jms+1:jme)
        N_coef(:,kms+1:kme,jms:jme-1) = &
                (dzdy(:,kms:kme-1,jms:jme-1)+dzdy(:,kms+1:kme,jms+1:jme))/mixed_denom(:,kms+1:kme,jms:jme-1)
        O_coef(:,kms+1:kme,jms+1:jme) = &
                -(dzdy(:,kms:kme-1,jms+1:jme)+dzdy(:,kms+1:kme,jms:jme-1))/mixed_denom(:,kms+1:kme,jms+1:jme)

        J_coef(ims:ime-1,kms,:) = (dzdx(ims:ime-1,kms,:)+dzdx(ims+1:ime,kms,:))/mixed_denom(ims:ime-1,kms,:)
        K_coef(ims+1:ime,kms,:) = -(dzdx(ims+1:ime,kms,:)+dzdx(ims:ime-1,kms,:))/mixed_denom(ims+1:ime,kms,:)
        N_coef(:,kms,jms:jme-1) = (dzdy(:,kms,jms:jme-1)+dzdy(:,kms,jms+1:jme))/mixed_denom(:,kms,jms:jme-1)
        O_coef(:,kms,jms+1:jme) = -(dzdy(:,kms,jms+1:jme)+dzdy(:,kms,jms:jme-1))/mixed_denom(:,kms,jms+1:jme)

        !H_coef = 0
        !I_coef = 0
        !J_coef = 0
        !K_coef = 0
        !L_coef = 0
        !M_coef = 0
        !N_coef = 0
        !O_coef = 0
        
    end subroutine initialize_coefs
    
    
    !Update the coefs which change with time, i.e. those which depend on alpha
    subroutine update_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain        

        B_coef(:,kms:kme-1,:) = sigma(:,kms:kme-1,:) * &
                              ( (1./alpha(ims:ime,kms:kme-1,jms:jme)**2 + dzdy2(ims:ime,kms:kme-1,jms:jme) + &
                              dzdx2(ims:ime,kms:kme-1,jms:jme)) * (1./domain%jacobian(ims:ime,kms:kme-1,jms:jme)) + &
                              (1./alpha(ims:ime,kms+1:kme,jms:jme)**2 + dzdy2(ims:ime,kms+1:kme,jms:jme) + &
                              dzdx2(ims:ime,kms+1:kme,jms:jme)) * (1./domain%jacobian(ims:ime,kms+1:kme,jms:jme))) / &
                          (2*(sigma(:,kms:kme-1,:)+sigma(:,kms:kme-1,:)**2)*dz_if(:,kms+1:kme,:)**2)
                          
                          
        C_coef(:,kms+1:kme,:) = ( (1./alpha(ims:ime,kms:kme-1,jms:jme)**2 + dzdy2(ims:ime,kms:kme-1,jms:jme) + &
                              dzdx2(ims:ime,kms:kme-1,jms:jme)) * (1./domain%jacobian(ims:ime,kms:kme-1,jms:jme)) + &
                              (1./alpha(ims:ime,kms+1:kme,jms:jme)**2 + dzdy2(ims:ime,kms+1:kme,jms:jme) + &
                              dzdx2(ims:ime,kms+1:kme,jms:jme)) * (1./domain%jacobian(ims:ime,kms+1:kme,jms:jme))) / &
                          (2*(sigma(:,kms+1:kme,:)+sigma(:,kms+1:kme,:)**2)*dz_if(:,kms+2:kme+1,:)**2)
                
        C_coef(:,kms,:) = ( (1./alpha(ims:ime,kms,jms:jme)**2 + dzdy2(ims:ime,kms,jms:jme) + &
                              dzdx2(ims:ime,kms,jms:jme)) * (1./domain%jacobian(ims:ime,kms,jms:jme)) + &
                              (1./alpha(ims:ime,kms,jms:jme)**2 + dzdy2(ims:ime,kms,jms:jme) + &
                              dzdx2(ims:ime,kms,jms:jme)) * (1./domain%jacobian(ims:ime,kms,jms:jme))) / &
                          (2*(sigma(:,kms,:)+sigma(:,kms,:)**2)*dz_if(:,kms+1,:)**2)
                          
        
        A_coef(ims+1:ime-1,kms:kme,jms+1:jme-1) = -((domain%jacobian(ims+2:ime,kms:kme,jms+1:jme-1) + &
                                               2*domain%jacobian(ims+1:ime-1,kms:kme,jms+1:jme-1) + &
                                               domain%jacobian(ims:ime-2,kms:kme,jms+1:jme-1))/(2*domain%dx**2)) &
                                            -((domain%jacobian(ims+1:ime-1,kms:kme,jms+2:jme) + &
                                               2*domain%jacobian(ims+1:ime-1,kms:kme,jms+1:jme-1) + &
                                               domain%jacobian(ims+1:ime-1,kms:kme,jms:jme-2))/(2*domain%dx**2)) - &
                                               B_coef(ims+1:ime-1,kms:kme,jms+1:jme-1) - C_coef(ims+1:ime-1,kms:kme,jms+1:jme-1)
        
                    
    end subroutine


end module wind_iterative