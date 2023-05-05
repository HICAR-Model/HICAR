!-----------------------------------------------------------------------
! Snow ablation at the top of the snowpack
!-----------------------------------------------------------------------
subroutine SNOW_ABLATION(dhs,dswe,i,j)

! This subroutine erodes a snow depth dhs corresponding to mass dswe at the top of the snowpack, and reduce the number of layers if necessary.

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use GRID, only: &
  Ds_min,            &! Minimum possible snow layer thickness (m)
  Nsmax               ! Maximum number of snow layers

use STATE_VARIABLES, only: &
  rgrn,              &! Snow layer grain radius (m)
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Ds,                &! Snow layer thicknesses (m)
  histowet,          &! Historical variable for past wetting of a layer (0-1)
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Tsnow               ! Snow layer temperatures (K)

implicit none

real, intent(in) :: &
  dhs,               &! HS decrease (m)
  dswe                ! SWE decrease (kg/m^2)

integer, intent(in) :: &
  i,j                 ! Point counters

real :: &
  dDs,               &! Snow thickness on the top of the snowpack (m)
  Ds_tmp,            &! Temporary variable for cumulated layers thickness (m)
  swe_tmp,           &! Temporary variable for cumulated layers SWE (kg/m^2)
  swe_layer,         &! Local layer SWE (kg/m^2)
  swe_layer_new       ! New local layer SWE (kg/m^2)

integer :: &
  l,k                  ! Snow layer counters

if (fsnow(i,j) > epsilon(fsnow)) then

  k = 0
  swe_tmp = 0.0
  Ds_tmp = 0.0
  dDs = dhs / fsnow(i,j)

  if (dDs > epsilon(dDs) .and. dswe > epsilon(dswe)) then

    do while (k <= Nsnow(i,j) .and. Ds_tmp < dDs)
      k = k + 1
      swe_layer = Sice(k,i,j) + Sliq(k,i,j)
      Ds_tmp = Ds_tmp + Ds(k,i,j)
      swe_tmp = swe_tmp + swe_layer
    end do

    Nsnow(i,j) = Nsnow(i,j) - k + 1

    if (Nsnow(i,j) >= 1) then
      Ds(1,i,j) = Ds_tmp - dDs
      swe_layer_new = swe_tmp - dswe
      if (swe_layer > epsilon(swe_layer)) then
        Sice(1,i,j) = Sice(k,i,j) * (swe_layer_new / swe_layer)
        Sliq(1,i,j) = Sliq(k,i,j) * (swe_layer_new / swe_layer)
      else
        Sice(1,i,j) = swe_layer_new
        Sliq(1,i,j) = 0.0
      end if
      histowet(1,i,j) = histowet(k,i,j)      
      rgrn(1,i,j) = rgrn(k,i,j)
      Tsnow(1,i,j) = Tsnow(k,i,j)
    end if

    if (Nsnow(i,j) >= 2) then
      do l = 2, Nsnow(i,j)
        Ds(l,i,j) = Ds(l+k-1,i,j)
        Sice(l,i,j) = Sice(l+k-1,i,j)
        Sliq(l,i,j) = Sliq(l+k-1,i,j)
        histowet(l,i,j) = histowet(l+k-1,i,j)      
        rgrn(l,i,j) = rgrn(l+k-1,i,j)
        Tsnow(l,i,j) = Tsnow(l+k-1,i,j)
      end do
    end if

    do l = Nsnow(i,j)+1,Nsmax
      Ds(l,i,j) = 0.0
      Sice(l,i,j) = 0.0
      Sliq(l,i,j) = 0.0
      histowet(l,i,j) = 0.0
      rgrn(l,i,j) = 0.0
      Tsnow(l,i,j) = Tm
    end do

    ! If the top layer gets too thin, aggregate it with the next one, if it exists
    if (Nsnow(i,j) >= 2 .and. Ds(1,i,j) < Ds_min) then
      Ds(1,i,j) = Ds(1,i,j) + Ds(2,i,j)
      Sice(1,i,j) = Sice(1,i,j) + Sice(2,i,j)
      Sliq(1,i,j) = Sliq(1,i,j) + Sliq(2,i,j)
      histowet(1,i,j) = histowet(2,i,j)      
      rgrn(1,i,j) = rgrn(2,i,j)
      Tsnow(1,i,j) = Tsnow(2,i,j)
      Nsnow(i,j) = Nsnow(i,j) - 1
      if (Nsnow(i,j) >= 2) then
        do l = 2, Nsnow(i,j)
          Ds(l,i,j) = Ds(l+1,i,j)
          Sice(l,i,j) = Sice(l+1,i,j)
          Sliq(l,i,j) = Sliq(l+1,i,j)
          histowet(l,i,j) = histowet(l+1,i,j)      
          rgrn(l,i,j) = rgrn(l+1,i,j)
          Tsnow(l,i,j) = Tsnow(l+1,i,j)
        end do
      end if
      do l = Nsnow(i,j)+1,Nsmax
        Ds(l,i,j) = 0.0
        Sice(l,i,j) = 0.0
        Sliq(l,i,j) = 0.0
        histowet(l,i,j) = 0.0
        rgrn(l,i,j) = 0.0
        Tsnow(l,i,j) = Tm
      end do
    end if

  else
  
    write(*,*) 'WARNING SNOW_ABLATION, dDs', dDs, 'dswe', dswe
  
  end if

end if

end subroutine SNOW_ABLATION
