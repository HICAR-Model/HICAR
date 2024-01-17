!-----------------------------------------------------------------------
! Physical constants
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! output to HICAR variables which re decaled in PHYSICS subroutine but not in MODULES 
!-----------------------------------------------------------------------
module MODULES_interface

!use GRID, only: &
!  Nsmax,         &! Maximum number of snow layers
!  Nsoil,         &! Number of soil layers
!  Nx,Ny           ! Grid dimensions

!real, allocatable :: &
!  Roff(:,:),       &! Total runoff (kg/m^2)
!  Ps(:,:),           &! Surface pressure (Pa)
!  Qa(:,:),           &! Specific humidity (kg/kg)
!  RH(:,:),           &! Relative humidity (%)
!  Rf(:,:),           &! Rainfall rate (kg/m^2/s)
!  Sf(:,:),           &! Snowfall rate (kg/m^2/s)
!  Sf24h(:,:),        &! Snowfall 24hr (kg/m^2)
!  Sdif(:,:),         &! Diffuse shortwave radiation (W/m^2)
!  Sdir(:,:),         &! Direct-beam shortwave radiation (W/m^2)
!  Ta(:,:),           &! Air temperature (K)
!  Tv(:,:),           &! Time-varying transmissivity for direct SWR (-)
!  Ua(:,:)             ! Wind speed (m/s)

!! Fluxes
!real :: &
!  Esrf_(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
!  G_(Nx,Ny),          &! Heat flux into surface (W/m^2)
!  Gsoil_(Nx,Ny),      &! Heat flux into soil (W/m^2)
!  H_(Nx,Ny),           &! Sensible heat flux to the atmosphere (W/m^2)
!  LE_(Nx,Ny),          &! Latent heat flux to the atmosphere (W/m^2)
!  Melt_(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
!  Rnet_(Nx,Ny),       &! Net radiation (W/m^2)
!  Roff_(Nx,Ny),       &! Total runoff (kg/m^2)
!  meltflux_out_(Nx,Ny) ! Runoff from snowmelt at base of snow (kg/m^2)

! Fluxes
real, allocatable :: &
  Esrf_(:,:),       &! Moisture flux from the surface (kg/m^2/s)
  Gsoil_(:,:),      &! Heat flux into soil (W/m^2)
  H_(:,:),           &! Sensible heat flux to the atmosphere (W/m^2)
  LE_(:,:),          &! Latent heat flux to the atmosphere (W/m^2)
  Melt_(:,:),       &! Surface melt rate (kg/m^2/s)
  Rnet_(:,:),       &! Net radiation (W/m^2)
  Roff_(:,:),       &! Total runoff (kg/m^2)
  snowdepth_(:,:),  &! Snow depth (m)
  SWE_(:,:),        & ! Snow water equivalent (kg/m^2)
  KH_(:,:),         &! Eddy diffusivity for heat to the atmosphere (m/s)  
  meltflux_out_(:,:), &! Runoff from snowmelt at base of snow (kg/m^2)
  Sliq_out_(:,:),   & ! Total LWC (kg/m^2)
  dm_salt_(:,:),    &! SWE change due to saltation (kg/m^2)
  dm_susp_(:,:),    &! SWE change due to suspension (kg/m^2)
  dm_subl_(:,:),    &! SWE change due to sublimation (kg/m^2)
  dm_slide_(:,:)     ! SWE change due to snow slides (kg/m^2)

!real, allocatable :: &
!  Qsalt_u(:,:),       &! Moisture flux from the surface (kg/m^2/s)
!  Qsalt_v(:,:)        ! Heat flux into soil (W/m^2)

end module MODULES_interface
