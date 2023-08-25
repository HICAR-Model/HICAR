!-----------------------------------------------------------------------
! Read point driving data
!-----------------------------------------------------------------------
subroutine DRIVE_interface()

use STATE_VARIABLES, only: Tsrf

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

implicit none

!logical, intent(out) :: &
!  EoR                 ! End-of-run flag

real*4 :: &
  es(Nx,Ny),      &! Saturation vapour pressure (Pa)
  Tc(Nx,Ny)        ! Temperature (C)

integer :: i,j!,where,eastatus 

!Sdif(:,:)=0.0
!if (CANMOD == 1) then
  !inquire(unit=814, pos=where)
  !read(814,pos=where,IOSTAT=eastatus) ((Tv(i,j),j=1,Ny),i=1,Nx)
!endif

Ua = max(Ua, 0.1)

do i = 1,Nx
 do j = 1,Ny
  Sf(i,j) = Sf(i,j)/dt !! MJ: Please note that in OSHD the precipitaion is mm or kg per one hour and it is already rate. It uses Sf(i,j)/dt only to convert it to kg/m2 per second.
                       !! MJ: However, in HICAR we always have Sf as kg/m2 during lsm_dt and needs to be divided by lsm_dt to be rate. 
  Rf(i,j) = Rf(i,j)/dt
  Tc(i,j) = Ta(i,j) - Tm
  es(i,j) = e0*exp(17.5043*Tc(i,j)/(241.3 + Tc(i,j)))
  RH(i,j) = 100*Qa(i,j)*Ps(i,j)/(eps*es(i,j))
  
  !Qa(i,j) = (RH(i,j)/100)*eps*es(i,j)/Ps(i,j)
  
  ! Fix for Udir: in WindNinja outputs, if Ua=0, then Udir=NaN, which creates an error in SNOWTRAN3D
  if (SNTRAN == 1) then
    if (isnan(Udir(i,j))) then
      Udir(i,j) = 0.0
    end if
  end if

 end do
end do

! End of driving data file
!1 EoR = .true.
!if (eastatus < 0) then
! EoR = .true. ! End of file
!endif

!return

end subroutine DRIVE_interface
