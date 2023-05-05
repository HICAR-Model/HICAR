!-----------------------------------------------------------------------
! Lateral redistribution of snow through gravity
! Implementation of the SnowSlide model (Bernhardt and Schulz, 2010)
! in the FSM model by Louis Quéno
! ATTENTION: y is the W->E axis, while x is the S->N
! South of (i,j): (i-1,j)
! North of (i,j): (i+1,j)
! West of (i,j): (i,j-1)
! East of (i,j): (i,j+1)
!-----------------------------------------------------------------------
subroutine SNOWSLIDE(snowdepth0,Sice0,dm_slide)

use GRID, only: &
  Nx,Ny                     ! Grid dimensions

use STATE_VARIABLES, only: &
  firstit,                 &! First iteration identifier
  fsnow,                   &! Snow cover fraction 
  Ds,                      &! Snow layer thicknesses (m)
  dm_tot_slide,            &! Cumulated SWE change due to snow slides (kg/m^2)
  index_grid_dem_sorted     ! Location (i,j) of sorted grid points

use LANDUSE, only: &
  dem,                     &! Terrain elevation (m)
  slope,                   &! Slope (deg)
  Shd                       ! Snow holding depth (m)

use PARAM_SNOWSLIDE, only: &
  dyn_ratio,               &! Dynamic snow holding depth ratio
  rho_deposit,             &! Constant snow avalanche deposit density (kg/m^3)
  slope_min,               &! Minimum slope for snow slide occurrence (deg)
  Shd_min                   ! Minimum snow holding depth (m)

implicit none

real, intent(inout) :: &
  snowdepth0(Nx,Ny),       &! Snow depth of deposited snow (m)
  Sice0(Nx,Ny),            &! Ice content of deposited snow (kg/m^2)
  dm_slide(Nx,Ny)           ! SWE change due to snow slides (kg/m^2)

integer :: &
  !is_sorted,               &! Boolean indicating if the DEM vector is sorted by decreasing elevations
  i,j,                     &! Point counters
  n                         ! Vector point counter

logical :: &
  snow_depo(Nx,Ny)          ! Boolean to identify pixels where snow is deposited

real :: &
  dswe,                    &! Mass of snow transported to a pixel (kg/m2)
  elev,                    &! Elevation with snowdepth of grid point n (m)
  elev_S,                  &! Elevation with snow depth of South grid point (m)
  elev_N,                  &! Elevation with snow depth of North grid point (m)
  elev_W,                  &! Elevation with snow depth of West grid point (m)
  elev_E,                  &! Elevation with snow depth of East grid point (m)
  elev_SW,                 &! Elevation with snow depth of South-West grid point (m)
  elev_SE,                 &! Elevation with snow depth of South-East grid point (m)
  elev_NW,                 &! Elevation with snow depth of North-West grid point (m)
  elev_NE,                 &! Elevation with snow depth of North-East grid point (m)
  w_S,                     &! Weight of South grid point for slide
  w_N,                     &! Weight of North grid point for slide
  w_W,                     &! Weight of West grid point for slide
  w_E,                     &! Weight of East grid point for slide
  w_SW,                    &! Weight of South-West grid point for slide
  w_SE,                    &! Weight of South-East grid point for slide
  w_NW,                    &! Weight of North-West grid point for slide
  w_NE,                    &! Weight of North-East grid point for slide
  wt,                      &! Ratio of fresh avalanche deposit available for new slide
  delev_tot,               &! Total elevation difference with lower pixels (m)
  snowdepth_available,     &! Available snow depth for slide (m)
  snowdepth_available2,    &! Available snow depth for slide after substracting fresh deposit(m)
  swe_available,           &! Available SWE for slide (kg/m2)
  swe_available2,          &! Available SWE for slide after substracting fresh deposit (kg/m2)
  snowdepth_updated         ! Local snow depth updated in the snow slide loop (m)

real :: &
  Shd_corr(Nx,Ny),         &! Snow holding depth, thresholded (m)
  snowdepth(Nx,Ny)          ! Snow depth averaged over the grid cell before slides (m)

! Initialize snow_depo array
snow_depo(:,:) = .FALSE.

! Compute snowdepth before snow slides, and threshold the snow holding depth with Shd_min
do j = 1, Ny
  do i = 1, Nx
    snowdepth(i,j) = sum(Ds(:,i,j)) * fsnow(i,j)
    Shd_corr(i,j) = max(Shd(i,j),Shd_min)
  end do
end do

! If first iteration of the model, we need to sort the DEM anyway
if (firstit == 1) then

  ! Sort grid points by decreasing elevation in a one-dimension vector
  call SORT_DEM(snowdepth)

!else
!
!  ! Check if the DEM vector is still sorted by decreasing elevation, despite snowdepth changes
!  is_sorted = 0
!  call CHECK_SORTED_DEM(snowdepth,is_sorted)
!
!  ! If the vector is not sorted anymore, we need to sort it again
!  if (is_sorted == 0) then
!    call SORT_DEM(snowdepth)
!  end if

end if

! Treating grid points from the highest to the lowest
do n = 1, Nx*Ny

  i = index_grid_dem_sorted(n,1)
  j = index_grid_dem_sorted(n,2)

  ! Start slide processes only if slope higher than the defined minimum. 
  ! If it is a pixel receiving avalanche snow, no slope threshold.
  if (slope(i,j) >= slope_min .or. snow_depo(i,j)) then

    ! Update snowdepth in case snow has been transported to this pixel earlier in the loop
    snowdepth_updated = sum(Ds(:,i,j)) * fsnow(i,j) + snowdepth0(i,j)

    ! Local elevation accounting for updated snowdepth
    elev = dem(i,j) + snowdepth_updated

    ! If an avalanche is occurring, the snow holding depth is reduced
    ! to micic the dynamic effect
    if (snow_depo(i,j)) then
      Shd_corr(i,j) = Shd_corr(i,j) * dyn_ratio
      Shd_corr(i,j) = max(Shd_corr(i,j),Shd_min)
    end if

    ! Compute the depth of snow available for avalanche transport
    snowdepth_available = max(0.0, snowdepth_updated - Shd_corr(i,j))

    ! if snowdepth_available == 0 then there is no snow available to slide.

    ! There is a snow excess at the pixel.
    if (snowdepth_available > epsilon(snowdepth_available)) then

      if (i == 1 .or. i == Nx .or. j == 1 .or. j == Ny) then
        ! Case 1: edge pixels. Excess snow is dumped out of the domain
        ! to avoid accumulation artefacts. It is not routed to domain pixels.

        ! Move first snow coming from fresh avalanche deposit
        if (snowdepth0(i,j) - snowdepth_available > epsilon(snowdepth0)) then

          wt = snowdepth_available / snowdepth0(i,j)
          snowdepth0(i,j) = snowdepth0(i,j) - snowdepth_available
          swe_available = wt * Sice0(i,j)
          Sice0(i,j) = Sice0(i,j) - swe_available

        else if (snowdepth_available - snowdepth0(i,j) > epsilon(snowdepth0)) then

          snowdepth_available2 = snowdepth_available - snowdepth0(i,j)
          snowdepth0(i,j) = 0.0
          swe_available = Sice0(i,j)
          Sice0(i,j) = 0.0

          ! Compute the mass of snow available for transport
          call SWE_FROM_HS(snowdepth_available2,swe_available2,i,j)

          call SNOW_ABLATION(snowdepth_available2,swe_available2,i,j)

          swe_available = swe_available + swe_available2

        else ! snowdepth_available == snowdepth0(i,j)

          snowdepth0(i,j) = 0.0
          swe_available = Sice0(i,j)
          Sice0(i,j) = 0.0

        end if

        dm_slide(i,j) = dm_slide(i,j) - swe_available
        dm_tot_slide(i,j) = dm_tot_slide(i,j) - swe_available

      else
        ! Case 2: inner domain pixels. Excess snow is routed to lower pixels if there are.
        ! Mass transfers are weighted by elevation (snowdepth included) differences.

        ! South of (i,j): (i-1,j)
        elev_S = dem(i-1,j) + sum(Ds(:,i-1,j)) * fsnow(i-1,j) + snowdepth0(i-1,j)
        w_S = max(0.0, elev - elev_S)
        ! North of (i,j): (i+1,j)
        elev_N = dem(i+1,j) + sum(Ds(:,i+1,j)) * fsnow(i+1,j) + snowdepth0(i+1,j)
        w_N = max(0.0, elev - elev_N)
        ! West of (i,j): (i,j-1)
        elev_W = dem(i,j-1) + sum(Ds(:,i,j-1)) * fsnow(i,j-1) + snowdepth0(i,j-1)
        w_W = max(0.0, elev - elev_W)
        ! East of (i,j): (i,j+1)
        elev_E = dem(i,j+1) + sum(Ds(:,i,j+1)) * fsnow(i,j+1) + snowdepth0(i,j+1)
        w_E = max(0.0, elev - elev_E)
        ! South-West of (i,j): (i-1,j-1)
        elev_SW = dem(i-1,j-1) + sum(Ds(:,i-1,j-1)) * fsnow(i-1,j-1) + snowdepth0(i-1,j-1)
        w_SW = max(0.0, elev - elev_SW)
        ! South-East of (i,j): (i-1,j+1)
        elev_SE = dem(i-1,j+1) + sum(Ds(:,i-1,j+1)) * fsnow(i-1,j+1) + snowdepth0(i-1,j+1)
        w_SE = max(0.0, elev - elev_SE)
        ! North-West of (i,j): (i+1,j-1)
        elev_NW = dem(i+1,j-1) + sum(Ds(:,i+1,j-1)) * fsnow(i+1,j-1) + snowdepth0(i+1,j-1)
        w_NW = max(0.0, elev - elev_NW)
        ! North-East of (i,j): (i+1,j+1)
        elev_NE = dem(i+1,j+1) + sum(Ds(:,i+1,j+1)) * fsnow(i+1,j+1) + snowdepth0(i+1,j+1)
        w_NE = max(0.0, elev - elev_NE)

        delev_tot = w_S + w_N + w_W + w_E + w_SW + w_SE + w_NW + w_NE
        ! if delev_tot == 0 then it's a "sink", no avalanche possible.

        if (delev_tot > epsilon(delev_tot)) then

          w_S = w_S / delev_tot
          w_N = w_N / delev_tot
          w_W = w_W / delev_tot
          w_E = w_E / delev_tot
          w_SW = w_SW / delev_tot
          w_SE = w_SE / delev_tot
          w_NW = w_NW / delev_tot
          w_NE = w_NE / delev_tot

          ! Move first snow coming from fresh avalanche deposit
          if (snowdepth0(i,j) - snowdepth_available > epsilon(snowdepth0)) then

            wt = snowdepth_available / snowdepth0(i,j)
            snowdepth0(i,j) = snowdepth0(i,j) - snowdepth_available
            swe_available = wt * Sice0(i,j)
            Sice0(i,j) = Sice0(i,j) - swe_available

          else if (snowdepth_available - snowdepth0(i,j) > epsilon(snowdepth0)) then

            snowdepth_available2 = snowdepth_available - snowdepth0(i,j)
            snowdepth0(i,j) = 0.0
            swe_available = Sice0(i,j)
            Sice0(i,j) = 0.0

            ! Compute the mass of snow available for transport
            call SWE_FROM_HS(snowdepth_available2,swe_available2,i,j)

            call SNOW_ABLATION(snowdepth_available2,swe_available2,i,j)

            swe_available = swe_available + swe_available2

          else ! snowdepth_available == snowdepth0(i,j)

            snowdepth0(i,j) = 0.0
            swe_available = Sice0(i,j)
            Sice0(i,j) = 0.0

          end if

          dm_slide(i,j) = dm_slide(i,j) - swe_available
          dm_tot_slide(i,j) = dm_tot_slide(i,j) - swe_available

          ! Transport to neighbour pixels, weighted by elevation difference
          ! Higher pixels have a weight of 0

          ! South
          if (w_S > epsilon(w_S)) then
            dswe = w_S * swe_available
            Sice0(i-1,j) = Sice0(i-1,j) + dswe
            snowdepth0(i-1,j) = snowdepth0(i-1,j) + dswe / rho_deposit
            snow_depo(i-1,j) = .TRUE.
            dm_slide(i-1,j) = dm_slide(i-1,j) + dswe
            dm_tot_slide(i-1,j) = dm_tot_slide(i-1,j) + dswe
          end if

          ! North
          if (w_N > epsilon(w_N)) then
            dswe = w_N * swe_available
            Sice0(i+1,j) = Sice0(i+1,j) + dswe
            snowdepth0(i+1,j) = snowdepth0(i+1,j) + dswe / rho_deposit
            snow_depo(i+1,j) = .TRUE.
            dm_slide(i+1,j) = dm_slide(i+1,j) + dswe
            dm_tot_slide(i+1,j) = dm_tot_slide(i+1,j) + dswe
          end if

          ! West
          if (w_W > epsilon(w_W)) then
            dswe = w_W * swe_available
            Sice0(i,j-1) = Sice0(i,j-1) + dswe
            snowdepth0(i,j-1) = snowdepth0(i,j-1) + dswe / rho_deposit
            snow_depo(i,j-1) = .TRUE.
            dm_slide(i,j-1) = dm_slide(i,j-1) + dswe
            dm_tot_slide(i,j-1) = dm_tot_slide(i,j-1) + dswe
          end if

          ! East
          if (w_E > epsilon(w_E)) then
            dswe = w_E * swe_available
            Sice0(i,j+1) = Sice0(i,j+1) + dswe
            snowdepth0(i,j+1) = snowdepth0(i,j+1) + dswe / rho_deposit
            snow_depo(i,j+1) = .TRUE.
            dm_slide(i,j+1) = dm_slide(i,j+1) + dswe
            dm_tot_slide(i,j+1) = dm_tot_slide(i,j+1) + dswe
          end if

          ! South-West
          if (w_SW > epsilon(w_SW)) then
            dswe = w_SW * swe_available
            Sice0(i-1,j-1) = Sice0(i-1,j-1) + dswe
            snowdepth0(i-1,j-1) = snowdepth0(i-1,j-1) + dswe / rho_deposit
            snow_depo(i-1,j-1) = .TRUE.
            dm_slide(i-1,j-1) = dm_slide(i-1,j-1) + dswe
            dm_tot_slide(i-1,j-1) = dm_tot_slide(i-1,j-1) + dswe
          end if

          ! South-East
          if (w_SE > epsilon(w_SE)) then
            dswe = w_SE * swe_available
            Sice0(i-1,j+1) = Sice0(i-1,j+1) + dswe
            snowdepth0(i-1,j+1) = snowdepth0(i-1,j+1) + dswe / rho_deposit
            snow_depo(i-1,j+1) = .TRUE.
            dm_slide(i-1,j+1) = dm_slide(i-1,j+1) + dswe
            dm_tot_slide(i-1,j+1) = dm_tot_slide(i-1,j+1) + dswe
          end if

          ! North-West
          if (w_NW > epsilon(w_NW)) then
            dswe = w_NW * swe_available
            Sice0(i+1,j-1) = Sice0(i+1,j-1) + dswe
            snowdepth0(i+1,j-1) = snowdepth0(i+1,j-1) + dswe / rho_deposit
            snow_depo(i+1,j-1) = .TRUE.
            dm_slide(i+1,j-1) = dm_slide(i+1,j-1) + dswe
            dm_tot_slide(i+1,j-1) = dm_tot_slide(i+1,j-1) + dswe
          end if

          ! North-East
          if (w_NE > epsilon(w_NE)) then
            dswe = w_NE * swe_available
            Sice0(i+1,j+1) = Sice0(i+1,j+1) + dswe
            snowdepth0(i+1,j+1) = snowdepth0(i+1,j+1) + dswe / rho_deposit
            snow_depo(i+1,j+1) = .TRUE.
            dm_slide(i+1,j+1) = dm_slide(i+1,j+1) + dswe
            dm_tot_slide(i+1,j+1) = dm_tot_slide(i+1,j+1) + dswe
          end if

        end if

      end if

    end if

  end if

end do

end subroutine SNOWSLIDE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CHECK_SORTED_DEM(snowdepth,is_sorted)

! This subroutine sorts the grid points by decreasing elevation in a one-dimension vector.

use GRID, only: &
  Nx,Ny                 ! Grid dimensions

use STATE_VARIABLES, only: &
  index_grid_dem_sorted ! Location (i,j) of sorted grid points

use LANDUSE, only: &
  dem                   ! Terrain elevation (m)

implicit none

real, intent(in) :: &
  snowdepth(Nx,Ny)      ! Snow depth before slides (m)

integer, intent(inout) :: &
  is_sorted             ! Boolean indicating if the DEM vector is sorted by decreasing elevations

integer :: &
  temp_sorted,         &! Temporary is_sorted
  i,j,i1,j1,           &! Grid point counters
  n                     ! Vector point counter

real :: &
  elev,                &! Local elevation accounting for snow depth at index n (m)
  elev_next             ! Local elevation accounting for snow depth at index n+1 (m)

write(*,*) '--CHECK_DEM'

n = 1
temp_sorted = 1
i1 = index_grid_dem_sorted(1,1)
j1 = index_grid_dem_sorted(1,2)
elev_next = dem(i1,j1) + snowdepth(i1,j1)

do while(n < Nx*Ny .and. temp_sorted == 1)

  i = index_grid_dem_sorted(n+1,1)
  j = index_grid_dem_sorted(n+1,2)
  elev = elev_next
  elev_next = dem(i,j) + snowdepth(i,j)
  if (elev_next > elev) then
    temp_sorted = 0
  end if
  n = n + 1

end do

is_sorted = temp_sorted

end subroutine CHECK_SORTED_DEM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SORT_DEM(snowdepth)

! This subroutine sorts the grid points by decreasing elevation in a one-dimension vector.

use GRID, only: &
  Nx,Ny                 ! Grid dimensions

use STATE_VARIABLES, only: &
  index_grid_dem_sorted ! Location (i,j) of sorted grid points

use LANDUSE, only: &
  dem                   ! Terrain elevation (m)

implicit none

real, intent(in) :: &
  snowdepth(Nx,Ny)      ! Snow depth before slides (m)

integer :: &
  iter,                &! Number of iterations
  index_first,         &! First boundary index of the bisection range
  index_last,          &! Last boundary index of the bisection range
  index_mid,           &! Middle index of the bisection range
  i,j,                 &! Grid point counters
  n                     ! Vector point counter

real :: &
  dem_with_snow_loc,   &! Local elevation accounting for snow depth at point i,j (m)
  snowdepth_loc         ! Local snow depth at point i,j (m)

real :: &
  dem_sorted(Nx*Ny)     ! Elevation with snow depth of grid points sorted by decreasing elevation (m)

write(*,*) '--SORT_DEM'

dem_sorted(:) = 0.0
iter = 0

! Fill sorted vectors
do j = 1, Ny
  do i = 1, Nx

    iter = iter + 1
    snowdepth_loc = snowdepth(i,j)
    dem_with_snow_loc = dem(i,j) !+ snowdepth_loc
    if (iter == 1) then
      dem_sorted(1) = dem_with_snow_loc
      index_grid_dem_sorted(1,1) = i
      index_grid_dem_sorted(1,2) = j
    else if (iter == 2) then
      if (dem_sorted(1) < dem_with_snow_loc) then
        dem_sorted(2) = dem_sorted(1)
        dem_sorted(1) = dem_with_snow_loc
        index_grid_dem_sorted(2,:) = index_grid_dem_sorted(1,:)
        index_grid_dem_sorted(1,1) = i
        index_grid_dem_sorted(1,2) = j
      else
        dem_sorted(2) = dem_with_snow_loc
        index_grid_dem_sorted(2,1) = i
        index_grid_dem_sorted(2,2) = j
      end if
    else
      ! Bisection method
      index_first = 1
      index_last = iter - 1
      do while (index_first < index_last)
        if (MODULO(index_last - index_first + 1,2) == 0) then ! even
          index_mid = index_first + ((index_last - index_first + 1) / 2) - 1
        else ! odd
          index_mid = index_first + ((index_last - index_first + 2) / 2) - 1
        end if
        if (dem_sorted(index_mid) < dem_with_snow_loc) then
          index_last = index_mid
        else
          index_first = index_mid + 1
        end if
      end do
      if (index_first == iter - 1 .and. dem_sorted(iter-1) > dem_with_snow_loc)  then
        dem_sorted(iter) = dem_with_snow_loc
        index_grid_dem_sorted(iter,1) = i
        index_grid_dem_sorted(iter,2) = j
      else
        do n = iter, index_first + 1, -1
          dem_sorted(n) = dem_sorted(n-1)
          index_grid_dem_sorted(n,:) = index_grid_dem_sorted(n-1,:)
        end do
        dem_sorted(index_first) = dem_with_snow_loc
        index_grid_dem_sorted(index_first,1) = i
        index_grid_dem_sorted(index_first,2) = j
      end if
    end if

  end do
end do

end subroutine SORT_DEM
