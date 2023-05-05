!-----------------------------------------------------------------------
! Factorial Snow Model interface 
!
!-----------------------------------------------------------------------
module FSM_interface

	use LANDUSE

    use MODCONF, only: SNTRAN,SNSLID

!-----------------------------------------------------------------------
! Meteorological driving variables
!-----------------------------------------------------------------------
    use DRIVING, only: &
      year,              &! Year
      month,             &! Month of year
      day,               &! Day of month
      hour,              &! Hour of day
      dt,                &! Timestep (s)
      LW,                &! Incoming longwave radiation (W/m2)
      Ps,                &! Surface pressure (Pa)
      Qa,                &! Specific humidity (kg/kg)
      Rf,                &! Rainfall rate (kg/m2/s)
      Sdif,              &! Diffuse shortwave radiation (W/m^2)
      Sdir,              &! Direct-beam shortwave radiation (W/m^2)
      Sf,                &! Snowfall rate (kg/m2/s)
      Sf24h,             &! Snowfall 24hr (kg/m2)
      Ta,                &! Air temperature (K)
      Ua,                &! Wind speed (m/s)
      Udir                ! Wind direction


!-----------------------------------------------------------------------
! Model state variables  
!-----------------------------------------------------------------------
    use STATE_VARIABLES, only: &
      Tsrf,              &! Surface skin temperature (K)
      Tsnow,             &! Snow layer temperatures (K)
      Sice,              &! Ice content of snow layers (kg/m^2)
      Sliq,              &! Liquid content of snow layers (kg/m^2)
      Ds,                &! Snow layer thicknesses (m)
      fsnow,             &! Surface skin temperature (K)
      Nsnow,             &! Number of snow layers
      Tsoil,             &! Soil layer temperatures (K)
      albs,              &! Snow albedo
      theta               ! Volumetric moisture content of soil layers
!-----------------------------------------------------------------------
! output to HICAR variables which re decaled in PHYSICS subroutine but not in MODULES 
!-----------------------------------------------------------------------
    use MODULES_interface
           
    use SNOWTRAN3D, only : SNOWTRAN3D_setup,SNOWTRAN3D_fluxes,SNOWTRAN3D_accum, STRAN_NORTH, STRAN_SOUTH, STRAN_EAST, STRAN_WEST, Utau, Utau_t, Ds_soft

    implicit none
    
    private
    !!
    public :: FSM_SETUP,FSM_DRIVE,FSM_PHYSICS,FSM_SNOWSLIDE, FSM_CUMULATE_SD, FSM_SNOWTRAN_SETUP, FSM_SNOWTRAN_FLUXES, FSM_SNOWTRAN_ACCUM
    
    public :: Nx_HICAR, Ny_HICAR,lat_HICAR,lon_HICAR,terrain_HICAR,dx_HICAR,slope_HICAR,shd_HICAR, STRAN_NORTH, STRAN_SOUTH, STRAN_EAST, STRAN_WEST, SNTRAN, SNSLID
    !!
    public ::            &
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
      Sf24h,             &
      Ta,                &
      Ua           
    public ::  	&
      firstit,           &
      Tsrf,              &
      Tsnow,             &
      Sice,              &
      Sliq,              &
      Ds,                &
      fsnow,             &
      Nsnow,             &
      Tsoil,             &
      albs,              &
      theta
    public ::  	&
      Esrf_,       &
      Gsoil_,      &
      H_,          &
      LE_,         &
      Melt_,       &
      Rnet_,       &
      Roff_,       &
      snowdepth_,  &
      SWE_,        &  
      KH_,         &  
      meltflux_out_, &
      Sliq_out_,     &
      dm_salt_,      &
      dm_susp_,      &
      dm_subl_,      &
      dm_slide_,     &
      Qsalt_u,       &
      Qsalt_v
      
    
    integer :: Nx_HICAR, Ny_HICAR    

    real, allocatable:: &
      lat_HICAR(:,:),   & ! lat from HICAR based WGS84
      lon_HICAR(:,:),   & ! lon from HICAR based WGS84
      slope_HICAR(:,:), & ! terrain slope
      shd_HICAR(:,:),   & ! snow holding depth
      terrain_HICAR(:,:)
 
    real :: & !!!!!!!!---->
      dx_HICAR            !   The grid spacing of the high-resolution data in HICAR
 
 contains
    
    !!!!!!!!!
    subroutine FSM_SETUP()
        implicit none
        call SETUP_interface()
    end subroutine FSM_SETUP
    !!
    subroutine FSM_DRIVE()
        implicit none
        call DRIVE_interface()
    end subroutine FSM_DRIVE    
    !!
     subroutine FSM_PHYSICS()
        implicit none  
        call PHYSICS_interface()
    end subroutine FSM_PHYSICS  
    
    subroutine FSM_CUMULATE_SD()
        implicit none
        
        call CUMULATE_SD_interface()
        
    end subroutine FSM_CUMULATE_SD


    subroutine FSM_SNOWSLIDE(aval,frame_in)
        implicit none
        
        logical, intent(in) :: &
          aval(Nx,Ny)
          
        logical, optional, intent(in) :: &
          frame_in
        
        real :: &
          snowdepth0(Nx,Ny), &
          Sice0(Nx,Ny), &
          dm_slide(Nx,Ny)
          
        logical :: FRAME
        
        if (SNSLID==1) then
            snowdepth0(:,:) = 0.
            Sice0(:,:) = 0.
            dm_slide(:,:) = 0.

            FRAME = .False.
            if (present(frame_in)) FRAME = frame_in
            
            call SNOWSLIDE_interface(snowdepth0,Sice0,dm_slide,FRAME,aval)
        
            ! Accumulation of new snow, calculation of snow cover fraction and relayering
            call SNOW_LAYERING(snowdepth0,Sice0)
        
            call CUMULATE_SNOWSLIDE_interface(dm_slide)
        endif
    end subroutine FSM_SNOWSLIDE

    subroutine FSM_SNOWTRAN_SETUP()
        implicit none
        
        if (SNTRAN==1) call SNOWTRAN3D_setup()
        
    end subroutine FSM_SNOWTRAN_SETUP

    subroutine FSM_SNOWTRAN_FLUXES(DIR)
        implicit none
        
        integer, intent(in) :: DIR
        
        if (SNTRAN==1) then
            if (present(DIR)) then
                call SNOWTRAN3D_fluxes(DIR)    
            else
                call SNOWTRAN3D_fluxes(STRAN_NORTH)    
                call SNOWTRAN3D_fluxes(STRAN_SOUTH)    
                call SNOWTRAN3D_fluxes(STRAN_EAST)    
                call SNOWTRAN3D_fluxes(STRAN_WEST)    
            endif
        endif 
        
    end subroutine FSM_SNOWTRAN_FLUXES
    
    subroutine FSM_SNOWTRAN_ACCUM()
        implicit none
        
        real :: &
          snowdepth0(Nx,Ny), &
          Sice0(Nx,Ny), &
          dm_salt(Nx,Ny), &
          dm_susp(Nx,Ny), &
          dm_subl(Nx,Ny), &
          dm_subgrid(Nx,Ny)

        if (SNTRAN==1) then
            snowdepth0(:,:) = 0.
            Sice0(:,:) = 0.
            dm_salt(:,:) = 0.
            dm_susp(:,:) = 0.
            dm_subl(:,:) = 0.
            dm_subgrid(:,:) = 0.

            call SNOWTRAN3D_accum(snowdepth0,Sice0,dm_salt,dm_susp,dm_subl,dm_subgrid)
        
            ! Accumulation of new snow, calculation of snow cover fraction and relayering
            call SNOW_LAYERING(snowdepth0,Sice0)
        
            call CUMULATE_SNOWTRAN3D_interface(dm_salt,dm_susp,dm_subl,dm_subgrid)
        endif
        
    end subroutine FSM_SNOWTRAN_ACCUM

end module FSM_interface
