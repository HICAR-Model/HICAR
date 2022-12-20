
!!
!!----------------------------------------------------------
module module_sf_FSMdrv
    use module_sf_noahdrv,   only : lsm_noah, lsm_noah_init
   !use module_lsm_basic,    only : lsm_basic
   !use module_lsm_simple,   only : lsm_simple, lsm_simple_init
    use module_water_simple, only : water_simple
    use mod_atm_utilities,   only : sat_mr
    use time_object,         only : Time_type
    use data_structures
    use options_interface,   only : options_t
    use domain_interface,    only : domain_t
    
    use FSM_interface , only:  FSM_SETUP,FSM_DRIVE,FSM_PHYSICS
    use FSM_interface , only:  Nx_HICAR, Ny_HICAR,lat_HICAR,lon_HICAR,terrain_HICAR,dx_HICAR
    use FSM_interface, only: &
      year,          &
      month,         &
      day,           &
      hour,          &
      dt,            &
      LW,            &
      Ps,            &
      Qa,            &
      Rf,            &
      Sdif,          &
      Sdir,          &
      Sf,            &
      Ta,            &
      Ua
    use FSM_interface, only: &
      Esrf_,         &
      Gsoil_,        &
      H_,            &
      LE_,           &
      Melt_,         &
      Rnet_,         &
      Roff_,         &
      snowdepth_,    &
      SWE_,          &  
      KH_,           &  
      meltflux_out_, &
      Sliq_out_ 
    use FSM_interface, only: &
      Tsrf,          &
      Tsnow,         &
      Sice,          &
      Sliq,          &
      Ds,            &
      fsnow,         &
      Nsnow,         &
      Tsoil,         &
      albs,          &
      theta
      
    implicit none

    private
    public :: lsm_FSM_init,lsm_FSM,Esrf_,H_,LE_,Tsrf,KH_
    public :: snowfall_sum, rainfall_sum, Roff_sum, meltflux_out_sum
    
    integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions
   
    real, allocatable :: &
        snowfall_sum(:,:),       & !aggregated per output interval
        rainfall_sum(:,:),       & !aggregated per output interval
        Roff_sum(:,:),           & !aggregated per output interval
        meltflux_out_sum(:,:)      !aggregated per output interval

contains
 
    subroutine lsm_FSM_init(domain,options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        integer :: i,j
        !!
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte
        !!       
        Nx_HICAR=ite-its+1
        Ny_HICAR=jte-jts+1
        !!
        allocate(lat_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(lon_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(terrain_HICAR(Nx_HICAR,Ny_HICAR))
        lat_HICAR=domain%latitude%data_2d(its:ite,jts:jte)
        lon_HICAR=domain%longitude%data_2d(its:ite,jts:jte)
        terrain_HICAR=domain%terrain%data_2d(its:ite,jts:jte)
        dx_HICAR=domain%dx
        !!
        allocate(Esrf_(Nx_HICAR,Ny_HICAR)); Esrf_=0.
        allocate(Gsoil_(Nx_HICAR,Ny_HICAR));Gsoil_=0.
        allocate(H_(Nx_HICAR,Ny_HICAR)); H_=0.
        allocate(LE_(Nx_HICAR,Ny_HICAR)); LE_=0.
        allocate(Melt_(Nx_HICAR,Ny_HICAR)); Melt_=0.
        allocate(Rnet_(Nx_HICAR,Ny_HICAR)); Rnet_=0.
        allocate(Roff_(Nx_HICAR,Ny_HICAR)); Roff_=0.
        allocate(snowdepth_(Nx_HICAR,Ny_HICAR)); snowdepth_=0.
        allocate(SWE_(Nx_HICAR,Ny_HICAR)); SWE_=0.0
        allocate(KH_(Nx_HICAR,Ny_HICAR)); KH_=0.0
        allocate(meltflux_out_(Nx_HICAR,Ny_HICAR)); meltflux_out_=0.
        allocate(Sliq_out_(Nx_HICAR,Ny_HICAR)); Sliq_out_=0.
        !!
        allocate(snowfall_sum(Nx_HICAR,Ny_HICAR)); snowfall_sum=0.
        allocate(rainfall_sum(Nx_HICAR,Ny_HICAR)); rainfall_sum=0.
        allocate(Roff_sum(Nx_HICAR,Ny_HICAR)); Roff_sum=0.
        allocate(meltflux_out_sum(Nx_HICAR,Ny_HICAR)); meltflux_out_sum=0.        
        !!
        call FSM_SETUP()
        !!        
        !! MJ added this block to read in while we use restart file:
        if (options%parameters%restart) then
            !! giving feedback to HICAR
            Tsrf = domain%skin_temperature%data_2d(its:ite,jts:jte)
            albs = domain%albs%data_2d(its:ite,jts:jte)
            fsnow = domain%fsnow%data_2d(its:ite,jts:jte)
            Nsnow = domain%Nsnow%data_2d(its:ite,jts:jte)                        
            !!
            do i=1,3
                Tsnow(i,:,:) = domain%Tsnow%data_3d(its:ite,i,jts:jte)
                Sice(i,:,:) = domain%Sice%data_3d(its:ite,i,jts:jte)
                Sliq(i,:,:) = domain%Sliq%data_3d(its:ite,i,jts:jte)
                Ds(i,:,:) = domain%Ds%data_3d(its:ite,i,jts:jte)
            enddo
            do i=1,4
                Tsoil(i,:,:) = domain%soil_temperature%data_3d(its:ite,i,jts:jte)
                theta(i,:,:) = domain%soil_water_content%data_3d(its:ite,i,jts:jte)
            enddo
        endif
        !!
        do j = 1, Ny_HICAR
            do i = 1, Nx_HICAR
                !if (this_image()==1) write(*,*) "  albsH, albsF  ",i, j, domain%albs%data_2d(i+its-1,j+jts-1), albs(i,j)
            end do
        end do
        !SYNC ALL
    end subroutine lsm_FSM_init

    subroutine lsm_FSM(domain,options,lsm_dt,current_rain,current_snow,windspd)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        !!       
        real, intent(in) :: &
            lsm_dt          ! Timestep (s) from HICAR for lsm
        real, dimension(:,:), intent(in) :: &   !! Note that: when the input is an array with whatever indexing either using (:,:) or (N,M), it always starts from index 1 in all dims
            current_rain,                 &! rainfall in kg m-2
            current_snow,                 &! snowfall in kg m-2
            windspd                        ! Wind speed (m/s)
        integer :: i,j
        !!
        !! giving the required input from HICAR to FSM
        year=real(domain%model_time%year)
        month=real(domain%model_time%month)
        day=real(domain%model_time%day)
        hour=real(domain%model_time%hour)
        dt=lsm_dt
        LW=domain%longwave%data_2d(its:ite,jts:jte)
        Ps=domain%surface_pressure%data_2d(its:ite,jts:jte)
        Rf=current_rain
        !!
        if (options%physics%radiation_downScaling==1) then 
            Sdir=domain%shortwave_direct%data_2d(its:ite,jts:jte)  !Sdir=domain%shortwave%data_2d(its:ite,jts:jte)
            Sdif=domain%shortwave_diffuse%data_2d(its:ite,jts:jte) !Sdif=0.0
        endif
        if (options%physics%radiation_downScaling==0) then 
            Sdir=domain%shortwave%data_2d(its:ite,jts:jte)
            Sdif=0.0
        endif
        !!
        if (options%parameters%factor_p_var == "") then 
            !if (this_image()==1) write(*,*) "facto_p is not read...FSM"
            Sf=current_snow
        endif
        if (options%parameters%factor_p_var /= "") then 
            !if (this_image()==1) write(*,*) "facto_p is read...FSM"
            Sf=current_snow*domain%factor_p%data_2d(its:ite,jts:jte)
        endif
        Ta= domain%temperature%data_3d(its:ite,domain%grid%kms,jts:jte)!domain%temperature_2m%data_2d(its:ite,jts:jte)
        Qa= domain%water_vapor%data_3d(its:ite,domain%grid%kms,jts:jte)!domain%humidity_2m%data_2d(its:ite,jts:jte)
        Ua=windspd
        !!  
        !! FSM processing      
        call FSM_DRIVE()
        call FSM_PHYSICS()
        !!
        !! reseting the water pixels....
        if (options%physics%watersurface==kWATER_SIMPLE) then
            do j=1,Ny_HICAR
                do i=1,Nx_HICAR
                    if (domain%land_mask(its+i-1,jts+j-1)==kLC_WATER) then
                        H_(i,j) = 0.0
                        LE_(i,j) = 0.0
                        SWE_(i,j) = 0.0
                        snowdepth_(i,j) = 0.0
                        Roff_(i,j) = 0.0
                        Roff_sum(i,j) = 0.0
                        meltflux_out_(i,j) = 0.0
                        meltflux_out_sum(i,j) = 0.0
                        Sliq_out_(i,j) = 0.0
                        !! skin temperature will be reset to sea surface temperature later in water scheme
                    endif
                end do
            end do
        endif
        !! giving feedback to HICAR
        domain%sensible_heat%data_2d(its:ite,jts:jte)=H_
        domain%latent_heat%data_2d(its:ite,jts:jte)=LE_
        domain%snow_water_equivalent%data_2d(its:ite,jts:jte)=SWE_
        !!
        domain%skin_temperature%data_2d(its:ite,jts:jte)=Tsrf
        domain%snowdepth%data_2d(its:ite,jts:jte)=snowdepth_
        domain%Sliq_out%data_2d(its:ite,jts:jte)=Sliq_out_
        !!
        do i=1,3
            domain%Tsnow%data_3d(its:ite,i,jts:jte) = Tsnow(i,:,:)
            domain%Sice%data_3d(its:ite,i,jts:jte) = Sice(i,:,:)
            domain%Sliq%data_3d(its:ite,i,jts:jte) = Sliq(i,:,:)
            domain%Ds%data_3d(its:ite,i,jts:jte) = Ds(i,:,:)
        enddo
        do i=1,4
            domain%soil_temperature%data_3d(its:ite,i,jts:jte) = Tsoil(i,:,:)
            domain%soil_water_content%data_3d(its:ite,i,jts:jte)=theta(i,:,:)
        enddo
        domain%fsnow%data_2d(its:ite,jts:jte)=fsnow
        domain%Nsnow%data_2d(its:ite,jts:jte)=Nsnow
        domain%albs%data_2d(its:ite,jts:jte)=albs
        !!
        snowfall_sum=snowfall_sum+Sf
        rainfall_sum=rainfall_sum+Rf
        Roff_sum=Roff_sum+Roff_
        meltflux_out_sum=meltflux_out_sum+meltflux_out_
        !!
        !SYNC ALL  
    end subroutine lsm_FSM
!!
end module module_sf_FSMdrv
