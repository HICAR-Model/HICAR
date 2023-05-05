!-----------------------------------------------------------------------
! Equivalent SWE of surface HS
!-----------------------------------------------------------------------
subroutine SWE_FROM_HS(hs,swe,i,j)

! This subroutine calculate the average SWE of a given snow depth at the top of the snowpack (e.g. to be eroded)

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
  hs                  ! Snow depth on the top of the snowpack (m)

real, intent(out) :: &
  swe                 ! SWE on the top of the snowpack (kg/m^2)

integer, intent(in) :: &
  i,j                 ! Point counters

real :: &
  dDs,               &! Snow thickness on the top of the snowpack (m)
  rho_avg,           &! Bulk snow density (kg/m^3)
  weight,            &! Layer weight
  Ds_tmp,            &! Temporary snow thickness variable (m)
  swe_tmp,           &! Temporary SWE variable (kg/m^2)
  swe_layer,         &! Layer SWE (kg/m^2)
  snowthickness       ! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)

integer :: &
  k                   ! Snow layer counters

! Initialization
rho_avg = rho_snow
swe = 0.0

if (fsnow(i,j) > epsilon(fsnow)) then
  
  dDs = hs / fsnow(i,j)
  snowthickness = sum(Ds(:,i,j))

  if (dDs > epsilon(dDs) .and. dDs <= snowthickness) then

    k = 1
    swe_tmp = 0.0
    Ds_tmp = 0.0

    do while (k <= Nsnow(i,j) .and. Ds_tmp < dDs)
      swe_layer = Sice(k,i,j) + Sliq(k,i,j)
      if (dDs - Ds_tmp > Ds(k,i,j)) then
        Ds_tmp = Ds_tmp + Ds(k,i,j)
        swe_tmp = swe_tmp + swe_layer
      else
        weight = (dDs-Ds_tmp)/Ds(k,i,j)
        swe_tmp = swe_tmp + weight*swe_layer
        Ds_tmp = dDs
      end if
      k = k+1
    end do

    swe = swe_tmp
    rho_avg = swe / hs

  else if (dDs < -epsilon(dDs)) then

    write(*,*) 'WARNING SWE_FROM_HS: dDs < 0. CHECK/DEBUG!', i, j
    write(*,*) 'dDs', dDs
    
  else if (dDs > snowthickness) then

    write(*,*) 'WARNING SWE_FROM_HS: dDs > snowthickness. CHECK/DEBUG!', i, j
    write(*,*) 'dDs', dDs, 'snowthickness', snowthickness

  end if

end if

if ((rho_avg < rhos_min - 0.5 .or. rho_avg > rhos_max + 0.5) .and. Ds_tmp > 0.001) then
  write(*,*) 'WARNING SWE_FROM_HS: invalid density.', rho_avg, i, j
  write(*,*) 'Ds: ', Ds(:,i,j)
  write(*,*) 'SWE: ', Sice(:,i,j) + Sliq(:,i,j)
  write(*,*) 'rho: ', (Sice(1:Nsnow(i,j),i,j) + Sliq(1:Nsnow(i,j),i,j)) / Ds(1:Nsnow(i,j),i,j) / fsnow(i,j)
end if

end subroutine SWE_FROM_HS
