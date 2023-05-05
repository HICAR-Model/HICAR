!-----------------------------------------------------------------------
! Read point driving data
!-----------------------------------------------------------------------
subroutine DRIVE(EoR)

use MODCONF, only: CANMOD,SNTRAN

use CONSTANTS, only: &
  eps,               &! Ratio of molecular weights of water and dry air
  e0,                &! Saturation vapour pressure at Tm (Pa)
  Tm                  ! Melting point (K)

!! ATTENTION: Sf and Rf are initially hourly precipitation (from OSHD matlab scripts), converted here to a precipitation rate
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
  Tv,                &! Time-varying transmissivity for direct SWR (-)
  RH,                &! Relative humidity (%)
  Ua,                &! Wind speed (m/s)
  Udir                ! Wind direction (degrees, clockwise from N)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

implicit none

logical, intent(out) :: &
  EoR                 ! End-of-run flag

real*4 :: &
  es(Nx,Ny),      &! Saturation vapour pressure (Pa)
  Tc(Nx,Ny)        ! Temperature (C)

integer :: i,j,where,eastatus 


! FSM driving data
inquire(unit=800, pos=where)
read(800,pos=where,IOSTAT=eastatus) year
inquire(unit=801, pos=where)
read(801,pos=where,IOSTAT=eastatus) month
inquire(unit=802, pos=where)
read(802,pos=where,IOSTAT=eastatus) day
inquire(unit=803, pos=where)
read(803,pos=where,IOSTAT=eastatus) hour
inquire(unit=804, pos=where)
read(804,pos=where,IOSTAT=eastatus) Sdir
inquire(unit=805, pos=where)
read(805,pos=where,IOSTAT=eastatus) Sdif
inquire(unit=806, pos=where)
read(806,pos=where,IOSTAT=eastatus) LW
inquire(unit=807, pos=where)
read(807,pos=where,IOSTAT=eastatus) Sf
inquire(unit=808, pos=where)
read(808,pos=where,IOSTAT=eastatus) Rf
inquire(unit=809, pos=where)
read(809,pos=where,IOSTAT=eastatus) Ta
inquire(unit=810, pos=where)
read(810,pos=where,IOSTAT=eastatus) RH
inquire(unit=811, pos=where)
read(811,pos=where,IOSTAT=eastatus) Ua
inquire(unit=812, pos=where)
read(812,pos=where,IOSTAT=eastatus) Ps
inquire(unit=813, pos=where)
read(813,pos=where,IOSTAT=eastatus) Sf24h
if (CANMOD == 1) then
  inquire(unit=814, pos=where)
  read(814,pos=where,IOSTAT=eastatus) Tv
endif
if (SNTRAN == 1) then
  inquire(unit=815, pos=where)
  read(815,pos=where,IOSTAT=eastatus) Udir
endif

Ua = max(Ua, 0.1)

do i = 1,Nx
do j = 1,Ny

  if (isnan(dem(i,j))) goto 1 ! Exclude points outside of the domain

  ! Convert precipitation rates from kg/m2/h to kg/m2/s
  Sf(i,j) = Sf(i,j) / 3600
  Rf(i,j) = Rf(i,j) / 3600
  Tc(i,j) = Ta(i,j) - Tm
  es(i,j) = e0*exp(17.5043*Tc(i,j)/(241.3 + Tc(i,j)))
  Qa(i,j) = (RH(i,j)/100)*eps*es(i,j)/Ps(i,j)
  
  ! Fix for Udir: in WindNinja outputs, if Ua=0, then Udir=NaN, which creates an error in SNOWTRAN3D
  if (SNTRAN == 1) then
    if (isnan(Udir(i,j))) then
      Udir(i,j) = 0.0
    end if
  end if
  
  1 continue 

end do
end do

! End of driving data file
!1 EoR = .true.
if (eastatus < 0) then
 EoR = .true. ! End of file
endif

return

end subroutine DRIVE
