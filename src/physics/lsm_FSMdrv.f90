
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
      year,              &
      month,             &
      day,               &
      hour,              &
      dt,                &
      LW,                &
      Ps,                &
      Qa,                &
      Rf,                &
      Sdif,              &
      Sdir,              &
      Sf,                &
      Ta,                &
      Ua
    use FSM_interface, only: &
      Esrf_,       &
      Gsoil_,      &
      H_,          &
      LE_,         &
      Melt_,       &
      Rnet_,       &
      Roff_,       &
      snowdepth_,  &
      SWE_,        &  
      KH_,        &  
      meltflux_out_ 
    use FSM_interface, only: &
      Tsrf,               &
      Tsnow,          &
      Sice,               &
      Sliq,               &
      Ds,                 &
      fsnow,          &
     Nsnow,         &
     Tsoil,             &
     albs,              &
     theta
      
   !use land_surface,   only : windspd,current_snow, SNOWBL, current_rain,snow_bucket
   !use land_surface,   only : ids,ide,jds,jde,kds,kde ! Domain dimensions
   !use land_surface,   only : ims,ime,jms,jme,kms,kme ! Local Memory dimensions
   !use land_surface,   only : its,ite,jts,jte,kts,kte ! Processing Tile dimensions

    implicit none

    private
    public :: lsm_FSM_init,lsm_FSM,Esrf_,H_,LE_,Tsrf,KH_
    public :: snowfall_sum, rainfall_sum, Roff_sum, meltflux_out_sum
   !public :: windspd_FSM,current_rain_FSM,current_snow_FSM !,windspd_FSM,
    
    integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions
   
   real, allocatable :: &
     snowfall_sum(:,:),       &!aggregated per output interval
     rainfall_sum(:,:),       &!aggregated per output interval
     Roff_sum(:,:),           &!aggregated per output interval
     meltflux_out_sum(:,:)     !aggregated per output interval

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
        !!
        !allocate(windspd_FSM(Nx_HICAR,Ny_HICAR)); windspd_FSM=0.
        !allocate(current_rain_FSM(Nx_HICAR,Ny_HICAR)); current_rain_FSM=0.
        !allocate(current_snow_FSM(Nx_HICAR,Ny_HICAR)); current_snow_FSM=0.
        !!
        allocate(snowfall_sum(Nx_HICAR,Ny_HICAR)); snowfall_sum=0.
        allocate(rainfall_sum(Nx_HICAR,Ny_HICAR)); rainfall_sum=0.
        allocate(Roff_sum(Nx_HICAR,Ny_HICAR)); Roff_sum=0.
        allocate(meltflux_out_sum(Nx_HICAR,Ny_HICAR)); meltflux_out_sum=0.        
        !!
        !if (this_image()==1) write(*,*) "  its,ite = ", its,ite
        !!
        call FSM_SETUP()
        !!
        !if (this_image()==1) write(*,*) "  ********Model time = ", trim(domain%model_time%as_string())
        !if (this_image()==1) write(*,*) "  ********YEARc= ", domain%model_time%year
        !if (this_image()==1) write(*,*) "  ********MONTH = ", domain%model_time%month
        !if (this_image()==1) write(*,*) "  ********DAY = ", domain%model_time%day
        !if (this_image()==1) write(*,*) "  ********MINUTE = ", domain%model_time%minute
        !if (this_image()==1) write(*,*) "  ********SECOND = ", domain%model_time%second
        
        !! MJ added this to read while we use restart file:
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
        !! Note that when the arg is an array with whatever indexing either using (:,:) or (N,M) it always starts from index 1 in all dims
        real, dimension(:,:), intent(in) :: &   
          current_rain,                 &! rainfall in kg m-2
          current_snow,                 &! snowfall in kg m-2
          windspd                        ! Wind speed (m/s)
        !! Note that when the arg is an array with whatever indexing either using (:,:) or (N,M) it always starts from index 1 in all dims
        !real, intent(in) :: &      
        !current_rain(:,:),            &! rainfall in kg m-2
        !current_snow(:,:),            &! snowfall in kg m-2
        !windspd(:,:)                   ! Wind speed (m/s)
        integer :: i,j
        !! giving the required input from HICAR to FSM
        year=real(domain%model_time%year)
        month=real(domain%model_time%month)
        day=real(domain%model_time%day)
        hour=real(domain%model_time%hour)
        dt=lsm_dt
        LW=domain%longwave%data_2d(its:ite,jts:jte)
        Ps=domain%surface_pressure%data_2d(its:ite,jts:jte)
        Rf=current_rain
        Sdir=domain%shortwave%data_2d(its:ite,jts:jte)
        Sdif=0.0
        Sf=current_snow
        Ta=domain%temperature_2m%data_2d(its:ite,jts:jte)
        Qa=domain%humidity_2m%data_2d(its:ite,jts:jte)
        Ua=windspd
        !!  
!       do j = 1, Nx_HICAR
!           do i = 1, Ny_HICAR
!               if ( isnan(LW(i,j)) ) write(*,*),"img-LW",i,j,this_image(),LW(i,j)
!               if ( isnan(Ps(i,j)) ) write(*,*),"img-Ps",i,j,this_image(),Ps(i,j)
!               if ( isnan(Rf(i,j)) ) write(*,*),"img-Rf",i,j,this_image(),Rf(i,j)
!               if ( isnan(Sdir(i,j)) ) write(*,*),"img-Sdir",i,j,this_image(),Sdir(i,j)
!               if ( isnan(Sf(i,j)) ) write(*,*),"img-Sf",i,j,this_image(),Sf(i,j)
!               if ( isnan(Ta(i,j)) .or. Ta(i,j)>=300.) write(*,*),"img-Ta",i,j,this_image(),Ta(i,j)
!               if ( isnan(Qa(i,j)) ) write(*,*),"img-Qa",i,j,this_image(),Qa(i,j)
!               if ( isnan(Ua(i,j)) .or. Ua(i,j)>10) write(*,*),"img-Ua",i,j,this_image(),Ua(i,j)
!               if ( isnan(Tsrf(i,j)) .or. Tsrf(i,j)>=300.) write(*,*),"img-Tsrf",i,j,this_image(),Tsrf(i,j)
!               if ( isnan(H_(i,j)) ) write(*,*),"img-H_",i,j,this_image(),H_(i,j)
!               if ( isnan(LE_(i,j)) ) write(*,*),"img-LE_",i,j,this_image(),LE_(i,j)
!          end do
!       end do
        !!
!       do j = 1, Nx_HICAR
!           do i = 1, Ny_HICAR
!               if ( isnan(H_(i,j)) .or. abs(H_(i,j))>=100. ) write(*,*),"img-H_",i,j,this_image(),H_(i,j)
!           end do
!       end do

        !! FSM processing      
        call FSM_DRIVE()
        call FSM_PHYSICS()
        !!
        do j = 1, NY_HICAR
            do i = 1, Nx_HICAR
                !if ( isnan(H_(i,j)) .or. abs(H_(i,j))>300 ) write(*,*) "img-H222",i,j,this_image(), H_(i,j), Tsrf(i,j), Ta(i,j), Ua(i,j), KH_(i,j)
                !if ( isnan(Ua(i,j)) .or. abs(Ua(i,j))>20 .or. Ua(i,j)<0) write(*,*),"befDr-Ua",i,j,this_image(), Ua(i,j), Tsrf(i,j), Ta(i,j), KH_(i,j), H_(i,j), windspd(i,j)
            end do
        end do
!       if (this_image()==1) then
!           print*, "DRIVE_interface***************-------", Nx,size(Qa,1)
!           print*, "DRIVE_interface***************-------", Ny,size(Qa,2) 
!           print*, "DRIVE_interface***************-------", dt 
!           do j = 1, 5
!               do i = 1, 5
!                   write(*,*),"lat....", Ua(i,j),Tsrf(i,j),Ta(i,j),Qa(i,j)
!               end do
!           end do
!       endif
!       do j = 1, Ny_HICAR
!           do i = 1, Nx_HICAR
!               if (this_image()==4) write(*,*),"match ",i, j, H_(i,j), KH_(i,j)!,windspd(its+i-1,jts+j-1) 
!           end do
!       end do
!       do j = 1, Ny_HICAR
!           do i = 1, Nx_HICAR
!              if( Sdir(i,j)>5)write(*,*),"fuck ",i, j, Sdir(i,j) 
!           end do
!       end do

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
                        meltflux_out_(i,j) = 0.0
                        !! skin temperature will be reset to sea surface temperature later in water scheme
                    endif
                end do
            end do
        endif
        !! giving feedback to HICAR
        domain%sensible_heat%data_2d(its:ite,jts:jte)=H_
        domain%latent_heat%data_2d(its:ite,jts:jte)=LE_
        domain%snow_water_equivalent%data_2d(its:ite,jts:jte)=SWE_
        !domain%runoff_tstep%data_2d(its:ite,jts:jte)=Roff_
        !domain%meltflux_out_tstep%data_2d(its:ite,jts:jte)=meltflux_out_
        !!
        domain%skin_temperature%data_2d(its:ite,jts:jte)=Tsrf
        domain%snowdepth%data_2d(its:ite,jts:jte)=snowdepth_
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
