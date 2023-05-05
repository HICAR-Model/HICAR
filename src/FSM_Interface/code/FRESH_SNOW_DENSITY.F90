!-----------------------------------------------------------------------
! Fresh snow density
!-----------------------------------------------------------------------
subroutine FRESH_SNOW_DENSITY(rhonew,i,j)

use MODCONF, only: DENSTY

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use DRIVING, only: &
  Ta,                &! Air temperature (K)
  Ua                  ! Wind speed (m/s)

use PARAMETERS, only: &
  rho0,              &! Fixed snow density (kg/m^3)
  rhob,              &! Temperature factor in fresh snow density (kg/m^3/K)
  rhoc,              &! Wind factor in fresh snow density (kg s^0.5/m^3.5)
  rhof,              &! Fresh snow density (kg/m^3)
  rhos_min            ! Minimum snow density (kg/m^3)
  
integer, intent(in) :: &
  i,j                 ! Point counters

real, intent(out) :: &
  rhonew              ! Density of new snow (kg/m^3)

real :: &
  t_decompaction      ! Decompaction time for fresh snow density (h)

if (DENSTY == 0) then
  rhonew = rho0
else
  ! Initial formulation
  !rhonew = max(rhof + rhob*(Ta(i,j) - Tm) + rhoc*Ua(i,j)**0.5, 50.)
  ! New formulation with decompaction
  rhonew = rhof + rhob*(Ta(i,j) - Tm) + rhoc*Ua(i,j)**0.5
  if (dem(i,j) <= 1000) then
    t_decompaction = 24
  else if (dem(i,j) > 2000) then
    t_decompaction = 0
  else
    t_decompaction = 24 + (dem(i,j) - 1000) / (2000 - 1000) * (0 - 24)
  end if
  rhonew = 300 + (rhonew - 300)*exp(t_decompaction/100)
  rhonew = max(rhonew, rhos_min)
endif

end subroutine FRESH_SNOW_DENSITY
