!-----------------------------------------------------------------------
! Cumulate fluxes
!-----------------------------------------------------------------------
subroutine CUMULATE_interface(Roff, meltflux_out,Sbsrf,Sdirt,Sdift,LWt,asrf_out,Melt, &
                    Esrf,Eveg,Gsoil,Hsrf,intcpt,KH,KHa,Khg,KHv,KWg,KWv,  &
                    LE,LEsrf,LWsci,LWveg,Rnet,Rsrf,Sbveg,H,Swsci,SWsrf,  &
                    SWveg,Usc,unload)
!MJ added-----------------------------------------------------------------
use MODULES_interface
!MJ added-----------------------------------------------------------------

use GRID, only: &
  Nsmax,         &! Maximum number of snow layers
  Nsoil,         &! Number of soil layers
  Nx,Ny           ! Grid dimensions

use STATE_VARIABLES, only: &
  albs,              &! Snow albedo
  fsnow,             &! Snow cover fraction
  Qcan,              &! Canopy air space humidity
  Ds,                &! Snow layer thicknesses (m)
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  snowdepthmin,      &! min Snow at time step of swemin (m)
  snowdepthmax,      &! max Snow at time step of swemax (m)
  swemin,            &! min swe during season (mm)
  swemax,            &! max swe during season (mm)  
  Sveg,              &! Canopy snow mass (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  Tsrf,              &! Surface skin temperature (K)
  Tsnow,             &
  Tsoil,             &
  Tveg                ! Vegetation temperature (K)

implicit none

real, intent(in) :: &
  asrf_out(Nx,Ny),   &! Surface albedo
  intcpt(Nx,Ny),     &! Canopy interception (kg/m^2)  
  Roff(Nx,Ny),       &! Total runoff (kg/m^2)
  meltflux_out(Nx,Ny),  &! Runoff from snowmelt at base of snow (kg/m^2)
  Esrf(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
  Eveg(Nx,Ny),       &! Moisture flux from vegetation (kg/m^2/s)$
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
  Hsrf(Nx,Ny),       &! Sensible heat flux from the surface (W/m^2)
  H(Nx,Ny),          &! Sensible heat flux to the atmosphere (W/m^2)
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHg(Nx,Ny),        &! Eddy diffusivity for heat from the ground (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny),        &! Eddy diffusivity for water from vegetation (m/s)
  LE(Nx,Ny),         &! Latent heat flux to the atmosphere (W/m^2)
  LEsrf(Nx,Ny),      &! Latent heat flux from the surface (W/m^2)
  LWsci(Nx,Ny),      &! Subcanopy incoming SWR (W/m^2)
  LWt(Nx,Ny),        &! Incoming longwave radiation corrected for subgrid topography (W/m^2)
  LWveg(Nx,Ny),      &! Net longwave radiation absorbed by vegetation (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny),       &! Net radiation (W/m^2)
  Rsrf(Nx,Ny),       &! Net radiation absorbed by the surface (W/m^2
  Sbsrf(Nx,Ny),      &! Sublimation from the snow surface (kg/m^2)
  Sbveg(Nx,Ny),      &! Sublimation from the vegetation (kg/m^2)
  Sdirt(Nx,Ny),      &! Incoming direct beam radiation corrected for subgrid topography (W/m^2)
  Sdift(Nx,Ny),      &! Incoming diffuse radiation corrected for subgrid topography (W/m^2)
  SWsci(Nx,Ny),      &! Subcanopy incoming SWR (W/m^2)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny),      &! Net SW radiation absorbed by vegetation (W/m^2)
  unload(Nx,Ny),     &! Snow mass unloaded from canopy (kg/m^2)
  Usc(Nx,Ny)          ! Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)               

  call CUMULATE_SD_interface()

!! MJ added-----------------------------------------------------------------
  Esrf_= Esrf
  Gsoil_ = Gsoil
  H_ = H
  LE_ = LE
  Melt_ = Melt
  Rnet_ = Rnet
  Roff_ = Roff
  meltflux_out_ = meltflux_out
  KH_=KH

end subroutine CUMULATE_interface


subroutine CUMULATE_SD_interface()

!MJ added-----------------------------------------------------------------
use MODULES_interface
!MJ added-----------------------------------------------------------------

use GRID, only: &
  Nx,Ny           ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq                ! Liquid content of snow layers (kg/m^2)

implicit none

real :: &
  Sliq_out(Nx,Ny),   &!
  snowdepth(Nx,Ny),  &! Snow depth (m)
  SWE(Nx,Ny)          ! Snow water equivalent (kg/m^2)

integer :: i,j,k

! BC just in case, these sums should be performed only until Nsnow.
do j = 1,Ny
  do i = 1,Nx
    Sliq_out(i,j) = sum(Sliq(:,i,j))
    snowdepth(i,j) = sum(Ds(:,i,j))
    SWE(i,j) = sum(Sice(:,i,j)) + sum(Sliq(:,i,j))
  end do
end do

  snowdepth_=snowdepth
  SWE_=SWE
  Sliq_out_=Sliq_out

end subroutine CUMULATE_SD_interface

subroutine CUMULATE_SNOWTRAN3D_interface(dm_salt,dm_susp,dm_subl,dm_subgrid)

!MJ added-----------------------------------------------------------------
use MODULES_interface
!MJ added-----------------------------------------------------------------

use GRID, only: &
  Nx,Ny           ! Grid dimensions

implicit none

real, intent(in) :: &
  dm_salt(Nx,Ny),   &! saltation flux per timestep (kg/m^2)
  dm_susp(Nx,Ny),   &! sublimation per timestep (kg/m^2)
  dm_subl(Nx,Ny),   &! blowing snow sublimation per timestep (kg/m^2)
  dm_subgrid(Nx,Ny)  ! subgrid redistribution per timestep (kg/m^2)

  call CUMULATE_SD_interface()

  dm_salt_=dm_salt
  dm_susp_=dm_susp
  dm_subl_=dm_subl
  dm_subgrid_=dm_subgrid

end subroutine CUMULATE_SNOWTRAN3D_interface

subroutine CUMULATE_SNOWSLIDE_interface(dm_slide)

!MJ added-----------------------------------------------------------------
use MODULES_interface
!MJ added-----------------------------------------------------------------

use GRID, only: &
  Nx,Ny           ! Grid dimensions

implicit none

real, intent(in) :: &
  dm_slide(Nx,Ny)    ! avalanched redistribution per timestep (kg/m^2)

  call CUMULATE_SD_interface()

  dm_slide_=dm_slide

end subroutine CUMULATE_SNOWSLIDE_interface
