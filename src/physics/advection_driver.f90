!> ----------------------------------------------------------------------------
!!  Driver to call different advection schemes
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module advection
    use data_structures
    use icar_constants
    use adv_std,                    only : adv_std_init, adv_std_var_request, adv_std_advect3d, adv_fluxcorr_advect3d, adv_std_compute_wind
    use adv_mpdata,                 only : mpdata_init, mpdata_advect3d, mpdata_compute_wind
    use adv_fluxcorr,               only : init_fluxcorr
    ! use debug_module,               only : domain_fix
    use options_interface,          only: options_t
    use domain_interface,           only: domain_t
    use variable_dict_interface,  only : var_dict_t
    use variable_interface,       only : variable_t

    implicit none
    private
    real, allocatable :: temp(:,:,:)
    integer :: ims, ime, kms, kme, jms, jme
    public :: advect, adv_init, adv_var_request
contains

    subroutine adv_init(domain,options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in) :: options

        if (this_image()==1) write(*,*) ""
        if (this_image()==1) write(*,*) "Initializing Advection"
        if (options%physics%advection==kADV_STD) then
            if (this_image()==1) write(*,*) "    Standard"
            call adv_std_init(domain,options)
        else if(options%physics%advection==kADV_MPDATA) then
            if (this_image()==1) write(*,*) "    MP-DATA"
            call mpdata_init(domain,options)
        endif
        
        ims = domain%ims; ime = domain%ime
        kms = domain%kms; kme = domain%kme
        jms = domain%jms; jme = domain%jme
        
        if (options%adv_options%flux_corr > 0) call init_fluxcorr(domain)
        !Allocate storage variable for temp-quantities
        if (options%time_options%RK3) then
            allocate(temp(ims:ime,kms:kme,jms:jme))
        endif

    end subroutine adv_init

    subroutine adv_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        !if (options%physics%advection==kADV_UPWIND) then
        !    call upwind_var_request(options)
        if (options%physics%advection==kADV_MPDATA) then
            call adv_std_var_request(options)
        else
            call adv_std_var_request(options)
        endif
    end subroutine
    
    subroutine advect(domain, options, dt)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,intent(in) :: dt
        type(variable_t) :: var_to_advect
        integer :: j, k
        ! integer :: nx, nz, ny
        !
        ! nx=size(domain%p,1)
        ! nz=size(domain%p,2)
        ! ny=size(domain%p,3)
        !
        ! if (.not.allocated(domain%tend%qv_adv)) then
        !     allocate(domain%tend%qv_adv(nx,nz,ny))
        !     domain%tend%qv_adv=0
        ! endif
        

        if (options%physics%advection==kADV_STD) then
            call adv_std_compute_wind(domain,options,dt)
        else if(options%physics%advection==kADV_MPDATA) then
            call mpdata_compute_wind(domain,options,dt)
        endif
        
        if (options%vars_to_advect(kVARS%water_vapor)>0) call adv_var(domain%water_vapor%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)
        if (options%vars_to_advect(kVARS%potential_temperature)>0) call adv_var(domain%potential_temperature%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options) 
        if (options%vars_to_advect(kVARS%cloud_water)>0) call adv_var(domain%cloud_water_mass%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)                  
        if (options%vars_to_advect(kVARS%rain_in_air)>0) call adv_var(domain%rain_mass%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)                    
        if (options%vars_to_advect(kVARS%snow_in_air)>0) call adv_var(domain%snow_mass%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)                    
        if (options%vars_to_advect(kVARS%cloud_ice)>0) call adv_var(domain%cloud_ice_mass%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)                      
        if (options%vars_to_advect(kVARS%graupel_in_air)>0) call adv_var(domain%graupel_mass%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)                 
        if (options%vars_to_advect(kVARS%ice_number_concentration)>0)  call adv_var(domain%cloud_ice_number%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)       
        if (options%vars_to_advect(kVARS%rain_number_concentration)>0) call adv_var(domain%rain_number%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)      
        if (options%vars_to_advect(kVARS%snow_number_concentration)>0) call adv_var(domain%snow_number%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)      
        if (options%vars_to_advect(kVARS%graupel_number_concentration)>0) call adv_var(domain%graupel_number%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   
        if (options%vars_to_advect(kVARS%ice1_a)>0) call adv_var(domain%ice1_a%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   
        if (options%vars_to_advect(kVARS%ice1_c)>0) call adv_var(domain%ice1_c%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   
        if (options%vars_to_advect(kVARS%ice2_mass)>0) call adv_var(domain%ice2_mass%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   
        if (options%vars_to_advect(kVARS%ice2_number)>0) call adv_var(domain%ice2_number%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   
        if (options%vars_to_advect(kVARS%ice2_a)>0) call adv_var(domain%ice2_a%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   
        if (options%vars_to_advect(kVARS%ice2_c)>0) call adv_var(domain%ice2_c%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   
        if (options%vars_to_advect(kVARS%ice3_mass)>0) call adv_var(domain%ice3_mass%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   
        if (options%vars_to_advect(kVARS%ice3_number)>0) call adv_var(domain%ice3_number%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   
        if (options%vars_to_advect(kVARS%ice3_a)>0) call adv_var(domain%ice3_a%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   
        if (options%vars_to_advect(kVARS%ice3_c)>0) call adv_var(domain%ice3_c%data_3d, domain%advection_dz, domain%jacobian, &
        domain%density%data_3d, options)   

    end subroutine advect

    subroutine adv_var(var, dz, jaco, rho, options)
        implicit none
        real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: var
        real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: dz, jaco, rho
        type(options_t),intent(in)    :: options
        
        integer :: j, k

        if (options%time_options%RK3) then
            if (options%physics%advection==kADV_STD) then

                !Initial advection-tendency calculations
                do j = jms,jme
                    do k = kms,kme
                        temp(:,k,j)  = var(:,k,j) 
                    enddo
                enddo
                call adv_std_advect3d(temp,var, dz, jaco,t_factor_in=0.333)
                call adv_std_advect3d(temp,var, dz, jaco,t_factor_in=0.5)

                !final advection call with tendency-fluxes
                call adv_fluxcorr_advect3d(temp,var, dz, jaco)

                do j = jms,jme
                    do k = kms,kme
                        var(:,k,j) = temp(:,k,j)
                    enddo
                enddo
            else if(options%physics%advection==kADV_MPDATA) then
                ! Not yet implemented (is it compatable w/ RK3?)
            endif
        else
            if (options%physics%advection==kADV_STD) then
                call adv_std_advect3d(var,var, dz, jaco)
            else if(options%physics%advection==kADV_MPDATA) then                                    
                call mpdata_advect3d(var, rho, jaco, dz, options)
            endif
        endif
    end subroutine adv_var

end module advection
