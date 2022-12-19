!-----------------------------------------------------------------------
! Write out state variables at end of run
!-----------------------------------------------------------------------
subroutine DUMP

use MODCONF, only: CANMOD,SNFRAC

use MODE_WRITE, only: WRITE_2D

use STATE_VARIABLES, only : &
  albs,              &! Snow albedo
  Ds,                &! Snow layer thicknesses (m)
  fsnow,             &! Snow cover fraction
  Nsnow,             &! Number of snow layers 
  Qcan,              &! Canopy air space humidity
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  snowdepthmin,      &! min Snow at time step of swemin (m)
  snowdepthmax,      &! max Snow at time step of swemax (m)
  snowdepthhist,     &! history of Snow depth during last 14 days (m)  
  swemin,            &! min swe during season (mm)
  swemax,            &! max swe during season (mm)  
  swehist,           &! history of Snow depth during last 14 days (kg/m^2)
  Sveg,              &! Canopy snow mass (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  Tsnow,             &! Snow layer temperatures (K)
  Tsoil,             &! Soil layer temperatures (K)
  Tsrf,              &! Surface skin temperature (K)
  Tveg                ! Vegetation temperature (K)

implicit none

! Write into state files.
write(1201) albs
write(1202) Ds
write(1203) fsnow
write(1204) Nsnow
write(1206) Sice
write(1207) Sliq
write(1216) Tsrf
write(1219) Tsnow
write(1220) Tsoil
if (SNFRAC == 0 .or. SNFRAC == 2) then
  write(1210) snowdepthmax
endif
if (SNFRAC == 0) then
  write(1209) snowdepthmin
  write(1211) snowdepthhist
  write(1213) swemin
  write(1214) swemax
  write(1215) swehist
endif
if (CANMOD == 1) then
  write(1223) Qcan
  write(1224) Sveg
  write(1225) Tcan
  write(1226) Tveg
endif

end subroutine DUMP
