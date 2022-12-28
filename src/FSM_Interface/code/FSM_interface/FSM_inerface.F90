!-----------------------------------------------------------------------
! Factorial Snow Model interface 
!
!-----------------------------------------------------------------------
module FSM_interface

	use LANDUSE

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
      Ua                  ! Wind speed (m/s)

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
    use MODULES_interface, only: &
      Esrf_,       &! Moisture flux from the surface (kg/m^2/s)
      Gsoil_,      &! Heat flux into soil (W/m^2)
      H_,          &! Sensible heat flux to the atmosphere (W/m^2)
      LE_,         &! Latent heat flux to the atmosphere (W/m^2)
      Melt_,       &! Surface melt rate (kg/m^2/s)
      Rnet_,       &! Net radiation (W/m^2)
      Roff_,       &! Total runoff (kg/m^2)
      snowdepth_,  &! Snow depth (m)
      SWE_,        &! Snow water equivalent (kg/m^2)	  
      KH_,        &! Snow water equivalent (kg/m^2)	  
      meltflux_out_, &! Runoff from snowmelt at base of snow (kg/m^2)
      Sliq_out_ ! Total LWC (kg/m^2)
           
    implicit none
    
    private
    !!
    public :: FSM_SETUP,FSM_DRIVE,FSM_PHYSICS,Nx_HICAR, Ny_HICAR,lat_HICAR,lon_HICAR,terrain_HICAR,dx_HICAR
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
      Sliq_out_
      
    
    integer :: Nx_HICAR, Ny_HICAR    

    real, allocatable:: &
      lat_HICAR(:,:),   & ! lat from HICAR based WGS84
      lon_HICAR(:,:),   & ! lon from HICAR based WGS84
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
end module FSM_interface
