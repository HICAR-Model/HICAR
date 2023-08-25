!-----------------------------------------------------------------------
! Call physics subroutines
!-----------------------------------------------------------------------
subroutine PHYSICS_interface

use MODCONF, only: CANMOD, RADSBG, CHECKS

use GRID, only: &
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  Nt,                &! Number of iterations of physics within an input timestep
                      ! (Nt = 4 implies a 15 min timestep for an hourly input timestep)
  Nitr                ! Number of iterations in energy balance calulation

use STATE_VARIABLES, only : &
  firstit             ! First iteration identifier

implicit none

! Eddy diffusivities
real :: &
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHg(Nx,Ny),        &! Eddy diffusivity for heat from the ground (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny)          ! Eddy diffusivity for water from vegetation (m/s)

! Surface properties
real :: &
  alb(Nx,Ny),        &! Albedo
  asrf_out(Nx,Ny),   &! Surface albedo
  Ds1(Nx,Ny),        &! Surface layer thickness (m)
  gs1(Nx,Ny),        &! Surface moisture conductance (m/s)
  ks1(Nx,Ny),        &! Surface thermal conductivity (W/m/K)
  Ts1(Nx,Ny)          ! Surface layer temperature (K)

! Snow properties
real :: &
  ksnow(Nsmax,Nx,Ny)  ! Thermal conductivity of snow (W/m/K)

! Soil properties
real :: &
  csoil(Nsoil,Nx,Ny),&! Areal heat capacity of soil (J/K/m^2)
  ksoil(Nsoil,Nx,Ny)  ! Thermal conductivity of soil (W/m/K)

! Fluxes
real :: &
  Esrf(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
  Eveg(Nx,Ny),       &! Moisture flux from vegetation (kg/m^2/s)
  G(Nx,Ny),          &! Heat flux into surface (W/m^2)
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
  H(Nx,Ny),          &! Sensible heat flux to the atmosphere (W/m^2)
  Hsrf(Nx,Ny),       &! Sensible heat flux from the surface (W/m^2)
  intcpt(Nx,Ny),     &! Canopy interception (kg/m^2)
  LE(Nx,Ny),         &! Latent heat flux to the atmosphere (W/m^2)
  LEsrf(Nx,Ny),      &! Latent heat flux from the surface (W/m^2)
  LWsci(Nx,Ny),      &! Subcanopy incoming SWR (W/m^2)
  LWveg(Nx,Ny),      &! Subcanopy incoming LWR (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny),       &! Net radiation (W/m^2)
  Roff(Nx,Ny),       &! Total runoff (kg/m^2)
  meltflux_out(Nx,Ny),  &! Runoff from snowmelt at base of snow (kg/m^2)
  Sdirt(Nx,Ny),      &! Incoming direct beam radiation corrected for subgrid topography (W/m^2)
  Sdift(Nx,Ny),      &! Incoming diffuse radiation corrected for subgrid topography (W/m^2)
  LWt(Nx,Ny),        &! Incoming longwave radiation corrected for subgrid topography (W/m^2)
  Rsrf(Nx,Ny),       &! Net radiation absorbed by the surface (W/m^2)
  Sbsrf(Nx,Ny),      &! Sublimation from the snow surface (kg/m^2)
  Sbveg(Nx,Ny),      &! Sublimation from the vegetation (kg/m^2)
  SWsci(Nx,Ny),      &! Subcanopy incoming SWR (W/m^2)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny),      &! Net SW radiation absorbed by vegetation (W/m^2)
  Tveg0(Nx,Ny),      &! Vegetation temperature at start of timestep (K)
  Usc(Nx,Ny),        &! Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)
  unload(Nx,Ny)       ! Snow mass unloaded from canopy (kg/m^2)

integer :: & 
  tt,                &! Iteration counter for physics
  n                   ! Iteration counter in energy balance calulation

do tt = 1, Nt

  call RADIATIONNN(alb,SWsrf,SWveg,Sdirt,Sdift,asrf_out,SWsci)

  call THERMAL(csoil,Ds1,gs1,ks1,ksnow,ksoil,Ts1,Tveg0)

  do n = 1, Nitr
    call SFEXCH(gs1,KH,KHa,KHg,KHv,KWg,KWv,Usc)

    if (CANMOD == 1) then
      call EBALFOR(Ds1,KHa,KHg,KHv,KWg,KWv,ks1,SWsrf,SWveg,Ts1,Tveg0, &
                  Esrf,Eveg,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Rsrf,LWsci,LWveg)
    endif

    if (RADSBG == 1) then
      call EBALSRF_SBG(Ds1,KH,KHa,KHv,KWg,KWv,ks1,SWsrf,SWveg,Ts1, &
                  Esrf,Eveg,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Rsrf,LWt)
    else
      call EBALSRF(Ds1,KH,KHa,KHv,KWg,KWv,ks1,SWsrf,SWveg,Ts1, &
                  Esrf,Eveg,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Rsrf,LWt, &
                  LWsci,LWveg)
    endif  
    
  end do

  call CANOPY(Eveg,unload,intcpt,Sbveg)

  call SNOW(Esrf,G,ksnow,ksoil,Melt,unload,Gsoil,Roff,meltflux_out,Sbsrf)

  call SOIL(csoil,Gsoil,ksoil)

end do

call CUMULATE_interface(Roff,meltflux_out,Esrf,Gsoil,KH,LE,Melt,Rnet,H)

if (CHECKS /= 0) then
  call CHECK_SNOWPACK
end if

end subroutine PHYSICS_interface
