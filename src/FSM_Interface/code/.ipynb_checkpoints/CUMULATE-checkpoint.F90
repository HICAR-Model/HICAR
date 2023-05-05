!-----------------------------------------------------------------------
! Cumulate fluxes
!-----------------------------------------------------------------------
subroutine CUMULATE(Roff, meltflux_out,Sbsrf,Sdirt,Sdift,LWt,asrf_out,Melt, &
                    Esrf,Eveg,Gsoil,Hsrf,intcpt,KH,KHa,Khg,KHv,KWg,KWv,  &
                    LE,LEsrf,LWsci,LWveg,Rnet,Rsrf,Sbveg,H,Swsci,SWsrf,  &
                    SWveg,Usc,unload,dm_salt,dm_susp,dm_subl,dm_subgrid,dm_slide)
                    
use MODE_WRITE, only: WRITE_2D

use DRIVING, only: &
  year,              &! Year
  month,             &! Month of year
  day,               &! Day of month
  hour                ! Hour of the day
  
use MODCONF, only: FOR_HN

use MODOUTPUT, only: &
  WRITE_DIAG_TABLE,  &
  WRITE_STATE_TABLE

use GRID, only: &
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
  Tveg,              &! Vegetation temperature (K)
  dm_tot_subl,       &! Cumulated SWE change due to sublimation (kg/m^2)
  dm_tot_trans,      &! Cumulated transported SWE (kg/m^2)
  dm_tot_slide        ! Cumulated SWE change due to snow slides (kg/m^2)
  
use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

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
  Usc(Nx,Ny),        &! Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)               
  dm_salt(Nx,Ny),    &! SWE change due to saltation (kg/m^2)
  dm_susp(Nx,Ny),    &! SWE change due to suspension (kg/m^2)
  dm_subl(Nx,Ny),    &! SWE change due to sublimation (kg/m^2)
  dm_subgrid(Nx,Ny), &! SWE change due to subgrid redistribution (kg/m^2)
  dm_slide(Nx,Ny)     ! SWE change due to snow slides (kg/m^2)

real :: &
  Sliq_out(Nx,Ny),   &!
  snowdepth(Nx,Ny),  &! Snow depth (m)
  SWE(Nx,Ny)          ! Snow water equivalent (kg/m^2)

integer :: i,j,where,ii

! BC just in case, these sums should be performed only until Nsnow.
do j = 1,Ny
do i = 1,Nx

  if (isnan(dem(i,j))) goto 1 ! Exclude points outside of the domain

  Sliq_out(i,j) = sum(Sliq(:,i,j))
  snowdepth(i,j) = sum(Ds(:,i,j)) * fsnow(i,j)
  SWE(i,j) = sum(Sice(:,i,j)) + sum(Sliq(:,i,j))
  
  1 continue 

end do
end do

inquire(unit=1401, pos=where)
write(1401,pos=where) year
inquire(unit=1402, pos=where)
write(1402,pos=where) month
inquire(unit=1403, pos=where)
write(1403,pos=where) day
inquire(unit=1404, pos=where)
write(1404,pos=where) hour

ii = 1
if (write_diag_table(ii)) call WRITE_2D(Roff,1404 + ii)        ! Total runoff, snow and bare soil (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(snowdepth, 1404 + ii)  ! Snow depth (m)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(SWE, 1404 + ii)        ! Snow water equivalent (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Sliq_out, 1404 + ii)   ! Sliq_out(i,j) = Sliq(1,i,j) + Sliq(2,i,j) + Sliq(3,i,j)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Sdirt,1404 + ii)       ! Incoming direct beam radiation corrected for subgrid topography (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Sdift, 1404 + ii)      ! Incoming diffuse radiation corrected for subgrid topography (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(LWt, 1404 + ii)        ! Incoming longwave radiation corrected for subgrid topography (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(meltflux_out, 1404 + ii)  ! Runoff from snowmelt at base of snow (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Sbsrf, 1404 + ii)      ! Snow sublimation rate (kg/m^2/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(asrf_out, 1404 + ii)   ! Surface albedo
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Melt, 1404 + ii)       ! Surface melt rate (kg/m^2/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Esrf, 1404 + ii)       ! Moisture flux from the surface (kg/m^2/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Eveg, 1404 + ii)       ! Moisture flux from vegetation (kg/m^2/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Gsoil, 1404 + ii)      ! Heat flux into soil (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Hsrf, 1404 + ii)       ! Sensible heat flux from the surface (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(intcpt, 1404 + ii)     ! Canopy interception (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KH, 1404 + ii)         ! Eddy diffusivity for heat to the atmosphere (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KHa, 1404 + ii)        ! Eddy diffusivity for heat from the canopy air space (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KHg, 1404 + ii)        ! Eddy diffusivity for heat from the ground (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KHv, 1404 + ii)        ! Eddy diffusivity for heat from vegetation (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KWg, 1404 + ii)        ! Eddy diffusivity for water from the ground (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KWv, 1404 + ii)        ! Eddy diffusivity for water from vegetation (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(LE, 1404 + ii)         ! Latent heat flux to the atmosphere (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(LEsrf, 1404 + ii)      ! Latent heat flux from the surface (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(LWsci, 1404 + ii)      ! Subcanopy incoming longwave radiation (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(LWveg, 1404 + ii)      ! Net longwave radiation absorbed by vegetation (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Rnet, 1404 + ii)       ! Net radiation (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Rsrf, 1404 + ii)       ! Net radiation at surface (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Sbveg, 1404 + ii)      ! Sublimation from vegetation (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(H, 1404 + ii)          ! Sensible heat flux to the atmosphere (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(SWsci, 1404 + ii)      ! Subcanopy incoming shortwave radiation (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(SWsrf, 1404 + ii)      ! Net SW radiation absorbed by the surface (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(SWveg, 1404 + ii)      ! Net SW radiation absorbed by vegetation (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Usc, 1404 + ii)        ! Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(unload, 1404 + ii)     ! Snow mass unloaded from canopy (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(dm_salt, 1404 + ii)    ! SWE change due to saltation (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(dm_susp, 1404 + ii)    ! SWE change due to suspension (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(dm_subl, 1404 + ii)    ! SWE change due to sublimation (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(dm_subgrid, 1404 + ii) ! SWE change due to subgrid redistribution (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(dm_slide, 1404 + ii)   ! SWE change due to snow slides (kg/m^2)

! If necessary, write state vars into results files:
ii = 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(albs,1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(fsnow, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(Tsrf, 1500 + ii)
ii = ii + 1
! if (CANMOD==0) then BC this is commented on purpose. DO NOT introduce this if.
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(snowdepthmin, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(snowdepthmax, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(swemin, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(swemax, 1500 + ii)
ii = ii + 1
! else BC this is commented on purpose. DO NOT introduce this if.
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(Qcan, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(Sveg, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(Tcan, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(Tveg, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(dm_tot_subl, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(dm_tot_trans, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(dm_tot_slide, 1500 + ii)
! endif BC this is commented on purpose. DO NOT introduce this if.

! write 19h-states for the subsequent HN model
if (FOR_HN) then
  if (hour > 18.5 .and. hour < 19.5) then
    write(1227) Ds
    write(1228) Tsrf
    write(1229) Tsnow
    write(1230) Tsoil

    close(1227)
    close(1228)
    close(1229)
    close(1230)
  endif
endif

end subroutine CUMULATE
