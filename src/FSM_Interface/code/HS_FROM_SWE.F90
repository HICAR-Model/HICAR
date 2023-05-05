!-----------------------------------------------------------------------
! Equivalent snow depth of surface SWE
!-----------------------------------------------------------------------
subroutine HS_FROM_SWE(swe,hs,i,j)

! This subroutine calculate the snow depth of a given SWE at the top of the snowpack (e.g. to be eroded)

use PARAMETERS, only : &
  rhos_min,          &! Minimum snow density (kg/m^3)
  rhos_max            ! Maximum snow density (kg/m^3)

use PARAM_SNOWTRAN3D, only: &
  rho_snow            ! Constant snow density (kg/m^3)
  
use STATE_VARIABLES, only: &
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Ds                  ! Snow layer thicknesses (m)

implicit none

real, intent(in) :: &
  swe                 ! SWE on the top of the snowpack (kg/m^2)

real, intent(out) :: &
  hs                  ! Snow depth on the top of the snowpack (kg/m^2)

integer, intent(in) :: &
  i,j                 ! Point counters

real :: &
  rho_avg,           &! Bulk snow density (kg/m^3)
  weight,            &! Layer weight
  Ds_tmp,            &! Temporary snow thickness variable (m)
  swe_tmp,           &! Temporary SWE variable (kg/m^2)
  swe_layer,         &! Layer SWE (kg/m^2)
  swe_tot             ! Total snowpack SWE (kg/m^2)

integer :: &
  k                   ! Snow layer counters

! Initialization
rho_avg = rho_snow
hs = 0.0

swe_tot = sum(Sice(:,i,j)+Sliq(:,i,j))

if (swe > epsilon(swe) .and. swe <= swe_tot) then

  k = 1
  swe_tmp = 0.0
  Ds_tmp = 0.0

  do while (k <= Nsnow(i,j) .and. swe_tmp < swe)
    swe_layer = Sice(k,i,j) + Sliq(k,i,j)
    if (swe - swe_tmp > swe_layer) then
      Ds_tmp = Ds_tmp + Ds(k,i,j)
      swe_tmp = swe_tmp + swe_layer
    else
      weight = (swe-swe_tmp)/swe_layer
      Ds_tmp = Ds_tmp + weight*Ds(k,i,j)
      swe_tmp = swe
    end if
    k = k+1
  end do

  hs = Ds_tmp * fsnow(i,j)

  if (hs > epsilon(hs)) then
  
    rho_avg = swe / hs
    
  else if (hs < - epsilon(hs)) then
  
    write(*,*) 'WARNING HS_FROM_SWE: hs < 0. CHECK/DEBUG!', i, j
    write(*,*) 'hs', hs
    
  end if

else if (swe < -epsilon(swe)) then

  write(*,*) 'WARNING HS_FROM_SWE: swe < 0. CHECK/DEBUG!', i, j
  write(*,*) 'swe', swe
  
else if (swe > swe_tot) then

  write(*,*) 'WARNING HS_FROM_SWE: swe > swetot. CHECK/DEBUG!'
  write(*,*) 'swe', swe, 'swe_tot', swe_tot

end if

if ((rho_avg < rhos_min - 0.5 .or. rho_avg > rhos_max + 0.5) .and. Ds_tmp > 0.001) then
  write(*,*) 'WARNING HS_FROM_SWE: invalid density.', rho_avg, i, j
  write(*,*) 'Ds: ', Ds(:,i,j)
  write(*,*) 'SWE: ', Sice(:,i,j) + Sliq(:,i,j)
  write(*,*) 'rho: ', (Sice(1:Nsnow(i,j),i,j) + Sliq(1:Nsnow(i,j),i,j)) / Ds(1:Nsnow(i,j),i,j) / fsnow(i,j)
end if

end subroutine HS_FROM_SWE
