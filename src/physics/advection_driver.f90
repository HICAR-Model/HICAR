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
    use adv_upwind,                 only : upwind, upwind_init, upwind_var_request, upwind_advect3d, upwind_compute_wind
    use adv_4th,                    only : adv4, adv4_init, adv4_var_request, adv4_advect3d, adv4_compute_wind
    use adv_mpdata,                 only : mpdata, mpdata_init, mpdata_advect3d, mpdata_compute_wind
    use adv_fluxcorr,               only : init_fluxcorr
    ! use debug_module,               only : domain_fix
    use options_interface,          only: options_t
    use domain_interface,           only: domain_t
    use variable_dict_interface,  only : var_dict_t
    use variable_interface,       only : variable_t

    implicit none
    private
    type(var_dict_t) :: var_list


    public :: advect, adv_init, adv_var_request
contains

    subroutine adv_init(domain,options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options

        if (this_image()==1) write(*,*) ""
        if (this_image()==1) write(*,*) "Initializing Advection"
        if (options%physics%advection==kADV_UPWIND) then
            if (this_image()==1) write(*,*) "    Upwind"
            call upwind_init(domain,options)
        else if(options%physics%advection==kADV_MPDATA) then
            if (this_image()==1) write(*,*) "    MP-DATA"
            call mpdata_init(domain,options)
        else if(options%physics%advection==kADV_4TH) then
            if (this_image()==1) write(*,*) "    4th-order Centered Difference"
            call adv4_init(domain,options)
        endif
        
        if (options%adv_options%flux_corr > 0) call init_fluxcorr(domain)
        
        if (options%vars_to_advect(kVARS%water_vapor)>0) call var_list%add_var('qv', domain%water_vapor%meta_data)
        
        !if (options%vars_to_advect(kVARS%potential_temperature)>0) call var_list%add_var('theta', domain%potential_temperature%meta_data) 
        !if (options%vars_to_advect(kVARS%cloud_water)>0) call var_list%add_var('qc', domain%cloud_water_mass%meta_data)                  
        !if (options%vars_to_advect(kVARS%rain_in_air)>0) call var_list%add_var('qr', domain%rain_mass%meta_data)                    
        !if (options%vars_to_advect(kVARS%snow_in_air)>0) call var_list%add_var('qs', domain%snow_mass%meta_data)                    
        !if (options%vars_to_advect(kVARS%cloud_ice)>0) call var_list%add_var('qi', domain%cloud_ice_mass%meta_data)                      
        !if (options%vars_to_advect(kVARS%graupel_in_air)>0) call var_list%add_var('qg', domain%graupel_mass%meta_data)                 
        !if (options%vars_to_advect(kVARS%ice_number_concentration)>0)  call var_list%add_var('ni', domain%cloud_ice_number%meta_data)       
        !if (options%vars_to_advect(kVARS%rain_number_concentration)>0) call var_list%add_var('nr', domain%rain_number%meta_data)      
        !if (options%vars_to_advect(kVARS%snow_number_concentration)>0) call var_list%add_var('ns', domain%snow_number%meta_data)      
        !if (options%vars_to_advect(kVARS%graupel_number_concentration)>0) call var_list%add_var('ng', domain%graupel_number%meta_data)   


    end subroutine adv_init

    subroutine adv_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        !if (options%physics%advection==kADV_UPWIND) then
        !    call upwind_var_request(options)
        if (options%physics%advection==kADV_MPDATA) then
            call upwind_var_request(options)
        else
            call upwind_var_request(options)
        endif
    end subroutine
    
    subroutine advect(domain, options, dt)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,intent(in) :: dt

        type(variable_t) :: var_to_advect
        real, allocatable :: temp(:,:,:)
        

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
        
        !Allocate storage variable for temp-quantities
        if (options%time_options%RK3) allocate(temp(domain%ims:domain%ime,domain%kms:domain%kme,domain%jms:domain%jme))

        if (options%physics%advection==kADV_UPWIND) then
            call upwind_compute_wind(domain,options,dt)
        else if(options%physics%advection==kADV_MPDATA) then
            call mpdata_compute_wind(domain,options,dt)
        else if(options%physics%advection==kADV_4TH) then
            call adv4_compute_wind(domain,options,dt)
        endif


        !Loop through all vars to advect
                
        ! make sure the dictionary is reset to point to the first variable
        call var_list%reset_iterator()

        ! Now iterate through the dictionary as long as there are more elements present
        do while (var_list%has_more_elements())
            ! get the next variable
            var_to_advect = var_list%next()
            if (options%time_options%RK3) then
            
                if (options%physics%advection==kADV_UPWIND) then
                
                    !Initial advection-tendency calculations
                    temp = var_to_advect%data_3d
                    call upwind_advect3d(temp,var_to_advect%data_3d, domain%advection_dz, domain%jacobian,t_factor_in=0.333)
                    call upwind_advect3d(temp,var_to_advect%data_3d, domain%advection_dz, domain%jacobian,t_factor_in=0.5)
                                            
                    !final advection call with tendency-fluxes
                    call upwind_advect3d(temp,var_to_advect%data_3d, domain%advection_dz, domain%jacobian)
                                            
                    var_to_advect%data_3d = temp
                    
                else if (options%physics%advection==kADV_4TH) then
                
                    !Initial advection-tendency calculations
                    temp = var_to_advect%data_3d
                    call adv4_advect3d(temp,var_to_advect%data_3d, domain%advection_dz, domain%jacobian,t_factor_in=0.333)
                    call adv4_advect3d(temp,var_to_advect%data_3d, domain%advection_dz, domain%jacobian,t_factor_in=0.5)

                    !final advection call with tendency-fluxes
                    call adv4_advect3d(temp,var_to_advect%data_3d, domain%advection_dz, &
                                       domain%jacobian,flux_corr=options%adv_options%flux_corr)
                    var_to_advect%data_3d = temp
                             
                else if(options%physics%advection==kADV_MPDATA) then
                
                    ! Not yet implemented (is it compatable w/ RK3?)
                endif
            else
                if (options%physics%advection==kADV_UPWIND) then
                    call upwind_advect3d(var_to_advect%data_3d,var_to_advect%data_3d, domain%advection_dz, domain%jacobian)
                else if(options%physics%advection==kADV_MPDATA) then
                    call mpdata_advect3d(var_to_advect%data_3d, domain%jacobian, domain%advection_dz, domain%dx,dt,options)
                else if (options%physics%advection==kADV_4TH) then
                    call adv4_advect3d(var_to_advect%data_3d,var_to_advect%data_3d, domain%advection_dz, domain%jacobian)
                    !Cheep flux correction -- to be replaced                    
                    where(var_to_advect%data_3d < 0) var_to_advect%data_3d = 0
                endif
            endif
        enddo

    end subroutine advect

end module advection
