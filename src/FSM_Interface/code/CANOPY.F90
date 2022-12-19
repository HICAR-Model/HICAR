!-----------------------------------------------------------------------
! Mass balance of canopy snow
!-----------------------------------------------------------------------
subroutine CANOPY(Eveg,unload,intcpt,Sbveg)

use MODTILE, only: tthresh 

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use DRIVING, only: &
  dt,                &! Timestep (s)
  Sf                  ! Snowfall rate (kg/m2/s)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  tcnc,              &! Canopy unloading time scale for cold snow (s)
  tcnm,              &! Canopy unloading time scale for melting snow (s)
  psf,               &! Scaling factor for solid precipitation (within forest stand, at min fveg)
  psr                 ! Range of solid precipitation (within forest stand, spread min-max CC)
  
use PARAMMAPS, only: &
  scap                ! Canopy snow capacity (kg/m^2)

use STATE_VARIABLES, only: &
  Sveg,              &! Canopy snow mass (kg/m^2)
  Tveg                ! Vegetation temperature (K)

use LANDUSE, only: &
  fveg,              &! Canopy cover fraction
  pmultf,            &! Precip multiplier applied to open-area snowfall 
  tilefrac            ! Grid cell tile fraction

implicit none

real, intent(in) :: &
  Eveg(Nx,Ny)         ! Moisture flux from vegetation (kg/m^2/s)

real, intent(out) :: &
  intcpt(Nx,Ny),     &! Canopy interception (kg/m^2)
  Sbveg(Nx,Ny),      &! Sublimation from the vegetation (kg/m^2)
  unload(Nx,Ny)       ! Snow mass unloaded from canopy (kg/m^2)

real :: &
  Evegs,             &! Canopy snow sublimation rate (kg/m^2/s)
  tunl                ! Canopy snow unloading timescale (s)

integer :: & 
  i,j                 ! Grid coordinates

do j = 1, Ny
do i = 1, Nx 

  unload(i,j) = 0
  intcpt(i,j) = 0
  Sbveg(i,j)  = 0

  if (tilefrac(i,j) < tthresh) goto 1 ! exclude points outside tile of interest
  
  if (fveg(i,j) > epsilon(fveg(i,j))) then
    ! rescale precipitation to correct back precip multiplier applied to open area 
    Sf(i,j) = pmultf(i,j)*Sf(i,j)

    ! interception
    intcpt(i,j) = (scap(i,j) - Sveg(i,j))*(1 - exp(-fveg(i,j)*Sf(i,j)*dt/scap(i,j)))
    Sveg(i,j) = Sveg(i,j) + intcpt(i,j)
    Sf(i,j) = Sf(i,j) - intcpt(i,j)/dt 
    Sf(i,j) = (psf - psr*fveg(i,j))*Sf(i,j) ! including preferential deposition in canopy gaps; might have to be revisited to ensure mass conservation, potentially integrate with pmultf

    ! sublimation
    Evegs = 0
    if (Sveg(i,j) > epsilon(Sveg(i,j)) .or. Tveg(i,j) < Tm) Evegs = Eveg(i,j)
    Sveg(i,j) = Sveg(i,j) - Evegs*dt
    Sbveg(i,j) = Evegs*dt
    if (Sveg(i,j) < 0) Sbveg(i,j) =  Sbveg(i,j) + Sveg(i,j)
    Sveg(i,j) = max(Sveg(i,j), 0.)

    ! unloading
    tunl = tcnc
    if (Tveg(i,j) >= Tm) tunl = tcnm
    tunl = max(tunl, dt)
    unload(i,j) = Sveg(i,j)*dt/tunl
    Sveg(i,j) = Sveg(i,j) - unload(i,j)

  end if
  
  1 continue ! exclude points
  
end do
end do

end subroutine CANOPY
