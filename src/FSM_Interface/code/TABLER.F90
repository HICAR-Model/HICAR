!-----------------------------------------------------------------------
! Tabler model to be initialized before SnowTran-3D
! References: Liston et al. (2007)
! Author: Glen E. Liston
! Adapted to FSM 2.0 by Louis Quéno
!-----------------------------------------------------------------------

subroutine TABLER(tabler_nn,tabler_ne,tabler_ee,tabler_se, &
                  tabler_ss,tabler_sw,tabler_ww,tabler_nw)


use CONSTANTS, only: &
  rho_wat                 ! Density of water (kg/m^3)

use GRID, only: &
  Nx,Ny                   ! Grid dimensions

use PARAMMAPS, only: &
  vegsnowd_xy             ! Vegetation snow holding capacity (m)

use STATE_VARIABLES, only: &
  firstit,               &! First iteration identifier
  Ds,                    &! Snow layer thicknesses (m)
  Sice,                  &! Ice content of snow layers (kg/m^2)
  Sliq                    ! Liquid content of snow layers (kg/m^2)

use LANDUSE, only: &
  Ld,                    &! Grid cell size (m)
  dem                     ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  subgrid_flag,          &! Flag for subgrid parameterization
  tp_scale                ! Tabler precipitation scale

implicit none

real, intent(inout) :: &
  tabler_nn(Nx,Ny),      &! Tabler surfaces NN
  tabler_ss(Nx,Ny),      &! Tabler surfaces SS
  tabler_ee(Nx,Ny),      &! Tabler surfaces EE
  tabler_ww(Nx,Ny),      &! Tabler surfaces WW
  tabler_ne(Nx,Ny),      &! Tabler surfaces NE
  tabler_se(Nx,Ny),      &! Tabler surfaces SE
  tabler_sw(Nx,Ny),      &! Tabler surfaces SW
  tabler_nw(Nx,Ny)        ! Tabler surfaces NW

integer :: &
  i,j                     ! Point counters
  
real :: &
  swe,                   &! Snow water equivalent of the snowpack (kg/m^2)
  snow_d,                &! Snow depth (m)
  rho_snow_grid,         &! Bulk snowpack density (kg/m^3)
  delta_EW,              &! Grid cell size in EW direction (m)
  delta_SN                ! Grid cell size in SN direction (m)

real :: &
  topo_tmp(Nx,Ny),       &! Temporary topography variable (m)
  snow_d_tabler(Nx,Ny)    ! Snow depth Tabler

real, allocatable :: &
  tabler_nn_orig(:,:),   &! Initialized Tabler surface NN
  tabler_ss_orig(:,:),   &! Initialized Tabler surface SS
  tabler_ee_orig(:,:),   &! Initialized Tabler surface EE
  tabler_ww_orig(:,:),   &! Initialized Tabler surface WW
  tabler_ne_orig(:,:),   &! Initialized Tabler surface NE
  tabler_se_orig(:,:),   &! Initialized Tabler surface SE
  tabler_sw_orig(:,:),   &! Initialized Tabler surface SW
  tabler_nw_orig(:,:)     ! Initialized Tabler surface NW
  
character(len=80) :: file_tabler_nn_orig
character(len=80) :: file_tabler_ne_orig
character(len=80) :: file_tabler_ee_orig
character(len=80) :: file_tabler_se_orig
character(len=80) :: file_tabler_ss_orig
character(len=80) :: file_tabler_sw_orig
character(len=80) :: file_tabler_ww_orig
character(len=80) :: file_tabler_nw_orig


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialization of delta_EW and delta_SN (must be > 1 m)
delta_EW = Ld(1,1)
delta_SN = Ld(1,1)

if (subgrid_flag==1.0) then

  file_tabler_nn_orig = 'tabler_nn_orig.bin'
  file_tabler_ne_orig = 'tabler_ne_orig.bin'
  file_tabler_ee_orig = 'tabler_ee_orig.bin'
  file_tabler_se_orig = 'tabler_se_orig.bin'
  file_tabler_ss_orig = 'tabler_ss_orig.bin'
  file_tabler_sw_orig = 'tabler_sw_orig.bin'
  file_tabler_ww_orig = 'tabler_ww_orig.bin'
  file_tabler_nw_orig = 'tabler_nw_orig.bin'
  
  open(1601, file = file_tabler_nn_orig, form='unformatted', access='stream', action='readwrite')
  open(1602, file = file_tabler_ne_orig, form='unformatted', access='stream', action='readwrite')
  open(1603, file = file_tabler_ee_orig, form='unformatted', access='stream', action='readwrite')
  open(1604, file = file_tabler_se_orig, form='unformatted', access='stream', action='readwrite')
  open(1605, file = file_tabler_ss_orig, form='unformatted', access='stream', action='readwrite')
  open(1606, file = file_tabler_sw_orig, form='unformatted', access='stream', action='readwrite')
  open(1607, file = file_tabler_ww_orig, form='unformatted', access='stream', action='readwrite')
  open(1608, file = file_tabler_nw_orig, form='unformatted', access='stream', action='readwrite') 
    
  allocate(tabler_nn_orig(Nx,Ny))
  allocate(tabler_ss_orig(Nx,Ny))
  allocate(tabler_ee_orig(Nx,Ny))
  allocate(tabler_ww_orig(Nx,Ny))
  allocate(tabler_ne_orig(Nx,Ny))
  allocate(tabler_se_orig(Nx,Ny))
  allocate(tabler_sw_orig(Nx,Ny))
  allocate(tabler_nw_orig(Nx,Ny))

  ! Perform some initialization steps that are unique to SnowTran-3D.
  if (firstit==1) then ! Initialization if first model iteration

    ! If this is the first time through, generate the Tabler snow
    !   accumulation surfaces for the land topography.
    call tabler_3d(dem,delta_EW,delta_SN, &
                   tabler_ww_orig,tabler_ee_orig,tabler_ss_orig, &
                   tabler_nn_orig,tabler_ne_orig, &
                   tabler_se_orig,tabler_sw_orig,tabler_nw_orig)

    ! As part of generating the Tabler surfaces, a -8888.0 has been
    !   used to identify the areas immediately upwind of any Tabler
    !   drift trap that is an erosion area where no snow is allowed
    !   to accumulate (see 'erosion_dist' in snowmodel.par).  Now
    !   take those areas and set them equal to the snow-holding depth
    !   to keep the snow relatively thin in those areas.
    do i = 1, Nx
      do j = 1, Ny
      
        if (isnan(dem(i,j))) goto 1 ! Exclude points outside of the domain

        tabler_nn_orig(i,j) = max(tabler_nn_orig(i,j),vegsnowd_xy(i,j))
        tabler_ne_orig(i,j) = max(tabler_ne_orig(i,j),vegsnowd_xy(i,j))
        tabler_ee_orig(i,j) = max(tabler_ee_orig(i,j),vegsnowd_xy(i,j))
        tabler_se_orig(i,j) = max(tabler_se_orig(i,j),vegsnowd_xy(i,j))
        tabler_ss_orig(i,j) = max(tabler_ss_orig(i,j),vegsnowd_xy(i,j))
        tabler_sw_orig(i,j) = max(tabler_sw_orig(i,j),vegsnowd_xy(i,j))
        tabler_ww_orig(i,j) = max(tabler_ww_orig(i,j),vegsnowd_xy(i,j))
        tabler_nw_orig(i,j) = max(tabler_nw_orig(i,j),vegsnowd_xy(i,j))
        
        1 continue  ! Exclude points outside of the domain
        
      end do
    end do
    
    ! Save initialized Tabler surfaces
    write(1601) tabler_nn_orig
    write(1602) tabler_ne_orig
    write(1603) tabler_ee_orig
    write(1604) tabler_se_orig
    write(1605) tabler_ss_orig
    write(1606) tabler_sw_orig
    write(1607) tabler_ww_orig
    write(1608) tabler_nw_orig
    
    close(1601)
    close(1602)
    close(1603)
    close(1604)
    close(1605)
    close(1606)
    close(1607)
    close(1608)
    
  else
  
    ! Open Tabler initialized surfaces
    read(1601) tabler_nn_orig
    read(1602) tabler_ne_orig
    read(1603) tabler_ee_orig
    read(1604) tabler_se_orig
    read(1605) tabler_ss_orig
    read(1606) tabler_sw_orig
    read(1607) tabler_ww_orig
    read(1608) tabler_nw_orig

    close(1601)
    close(1602)
    close(1603)
    close(1604)
    close(1605)
    close(1606)
    close(1607)
    close(1608)

  end if

  ! Generate the tabler surfaces at this time step, assuming the
  !   snow surface is the topographic surface.

  ! Take the snow depth coming out of SnowPack, and before any wind
  !   redistribution has been applied at this time step, and add it
  !   to dem.  Then use this to create the Tabler surface that
  !   will be used for this time step.  Also define this depth to be
  !   dependent on the SnowPack spatially distributed snow density,
  !   not the constant density used in SnowTran.
  do i = 1, Nx
    do j = 1, Ny
    
      if (isnan(dem(i,j))) goto 2 ! Exclude points outside of the domain

      swe = sum(Sice(:,i,j))+sum(Sliq(:,i,j))
      snow_d = sum(Ds(:,i,j))
      rho_snow_grid = swe / snow_d
      snow_d_tabler(i,j) = swe * rho_wat / rho_snow_grid
      topo_tmp(i,j) = tp_scale * snow_d_tabler(i,j) + dem(i,j)
      
      2 continue  ! Exclude points outside of the domain

    end do
  end do

  call tabler_3d(topo_tmp,delta_EW,delta_SN, &
                 tabler_ww_orig,tabler_ee_orig,tabler_ss_orig, &
                 tabler_nn_orig,tabler_ne_orig, &
                 tabler_se_orig,tabler_sw_orig,tabler_nw_orig)

  ! The Tabler surfaces that were just generated have had topo_tmp
  !   subtracted off of them, giving just the drift profiles with
  !   things like zero drift depth on ridges and windwards slopes.
  !   So, add the snow depth, prior to any wind redistribution, to
  !   these Tabler surfaces.  This will be the maximum snow depth
  !   allowed as part of the wind redistribution.
  do i = 1, Nx
    do j = 1, Ny
    
      if (isnan(dem(i,j))) goto 3 ! Exclude points outside of the domain
    
      tabler_ww(i,j) = tp_scale * snow_d_tabler(i,j) + tabler_ww_orig(i,j)
      tabler_ee(i,j) = tp_scale * snow_d_tabler(i,j) + tabler_ee_orig(i,j)
      tabler_ss(i,j) = tp_scale * snow_d_tabler(i,j) + tabler_ss_orig(i,j)
      tabler_nn(i,j) = tp_scale * snow_d_tabler(i,j) + tabler_nn_orig(i,j)
      tabler_ne(i,j) = tp_scale * snow_d_tabler(i,j) + tabler_ne_orig(i,j)
      tabler_se(i,j) = tp_scale * snow_d_tabler(i,j) + tabler_se_orig(i,j)
      tabler_sw(i,j) = tp_scale * snow_d_tabler(i,j) + tabler_sw_orig(i,j)
      tabler_nw(i,j) = tp_scale * snow_d_tabler(i,j) + tabler_nw_orig(i,j)
      
      3 continue  ! Exclude points outside of the domain
      
    end do
  end do

  ! As part of generating the Tabler surfaces, a -8888.0 has been
  !   used to identify the areas immediately upwind of any Tabler
  !   drift trap that is an erosion area where no snow is allowed
  !   to accumulate (see 'erosion_dist' in snowmodel.par).  Now
  !   take those areas and set them equal to the snow-holding depth
  !   to keep the snow relatively thin in those areas.
  do i = 1, Nx
    do j = 1, Ny
    
      if (isnan(dem(i,j))) goto 4 ! Exclude points outside of the domain
    
      tabler_nn(i,j) = max(tabler_nn(i,j),vegsnowd_xy(i,j))
      tabler_ne(i,j) = max(tabler_ne(i,j),vegsnowd_xy(i,j))
      tabler_ee(i,j) = max(tabler_ee(i,j),vegsnowd_xy(i,j))
      tabler_se(i,j) = max(tabler_se(i,j),vegsnowd_xy(i,j))
      tabler_ss(i,j) = max(tabler_ss(i,j),vegsnowd_xy(i,j))
      tabler_sw(i,j) = max(tabler_sw(i,j),vegsnowd_xy(i,j))
      tabler_ww(i,j) = max(tabler_ww(i,j),vegsnowd_xy(i,j))
      tabler_nw(i,j) = max(tabler_nw(i,j),vegsnowd_xy(i,j))
      
      4 continue  ! Exclude points outside of the domain
      
    end do
  end do

end if

end subroutine TABLER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tabler_3d(topo_land,delta_EW,delta_SN, &
                     tabler_ww,tabler_ee,tabler_ss,tabler_nn, &
                     tabler_ne,tabler_se,tabler_sw,tabler_nw)

! This subroutine uses Tabler (1975) to define equilibrium profiles
!   for the topographic drift catchments, for the case of winds
!   from the EAST, WEST, NORTH, and SOUTH, and anywhere inbetween.

use GRID, only: &
  Nx,Ny               ! Grid dimensions

implicit none

real, intent(in) :: &
  topo_land(Nx,Ny),  &! Topography (m)
  delta_EW,          &! Grid cell size in EW direction (m)
  delta_SN            ! Grid cell size in SN direction (m)

real, intent(inout) :: &
  tabler_nn(Nx,Ny),  &! Tabler surface NN
  tabler_ss(Nx,Ny),  &! Tabler surface SS
  tabler_ee(Nx,Ny),  &! Tabler surface EE
  tabler_ww(Nx,Ny),  &! Tabler surface WW
  tabler_ne(Nx,Ny),  &! Tabler surface NE
  tabler_se(Nx,Ny),  &! Tabler surface SE
  tabler_sw(Nx,Ny),  &! Tabler surface SW
  tabler_nw(Nx,Ny)    ! Tabler surface NW 

integer :: &
  irotate_flag        ! Flag for direction

! Here we generate maximum snow accumulation surfaces for n, ne, e,
!   se, s, sw, w, and nw winds.  I call these "tabler surfaces".  
!
! They are valid for the wind direction ranges: N=337.5-22.5,
!   NE=22.5-67.5, E=67.5-112.5, SE=112.5-157.5, S=157.5-202.5,
!   SW=202.5-247.5, W=247.5-292.5, and NW=292.5-337.5.
!
! These Tabler Surfaces define a "potential" snow surface that
!   represents the maximum possible snow-accumulation depth from winds
!   coming from these directions.
!
! Tabler, R. D., 1975: Predicting profiles of snowdrifts in
!   topographic catchments.  Proceedings of the 43rd Annual Western
!   Snow Conference, San Diego, California, 87-97.

! Consider N winds.
irotate_flag = 1
call tabler_n(topo_land,tabler_nn,delta_SN, &
              irotate_flag)

! Consider NE winds.
irotate_flag = 2
call tabler_e(topo_land,tabler_ne,1.41*delta_EW, &
              irotate_flag)

! Consider E winds.
irotate_flag = 1
call tabler_e(topo_land,tabler_ee,delta_EW, &
              irotate_flag)

! Consider SE winds.
irotate_flag = 2
call tabler_s(topo_land,tabler_se,1.41*delta_SN, &
              irotate_flag)

! Consider S winds.
irotate_flag = 1
call tabler_s(topo_land,tabler_ss,delta_SN, &
              irotate_flag)

! Consider SW winds.
irotate_flag = 2
call tabler_w(topo_land,tabler_sw,1.41*delta_EW, &
              irotate_flag)

! Consider W winds.
irotate_flag = 1
call tabler_w(topo_land,tabler_ww,delta_EW, &
              irotate_flag)

! Consider NW winds.
irotate_flag = 2
call tabler_n(topo_land,tabler_nw,1.41*delta_SN, &
              irotate_flag)

end subroutine tabler_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tabler_w(topo_land,tabler_ww,delta_EW, &
                    irotate_flag)

! This subroutine uses Tabler (1975) to define equilibrium profiles
!   for the topographic drift catchments, for the case of winds
!   from the WEST and NORTHWEST.

! Tabler, R. D., 1975: Predicting profiles of snowdrifts in
!   topographic catchments.  Proceedings of the 43rd Annual Western
!   Snow Conference, San Diego, California, 87-97.

use GRID, only: &
  Nx,Ny               ! Grid dimensions
  
use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  erosion_dist,      &! Erosion area
  slope_adjust        ! Slope adjustment factor

implicit none

real, intent(in) :: &
  topo_land(Nx,Ny),  &! Topography (m)
  delta_EW            ! Grid cell size in EW direction (m)

real, intent(inout) :: &
  tabler_ww(Nx,Ny)    ! Tabler surface WW

integer, intent(in) :: &
  irotate_flag        ! Flag

integer i,j,jstart,jend,II,nny,nnx,jj,jjj,maxlines
integer npts
real y,xmax_slope,dy,test,x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

real tabler(Ny)
real topo_line(Ny)
real drift_start_topo(Ny)

!integer (Ny_max_tabler=Ny*1000)
real topo_1m(Ny*1000)
real tabler_1m(Ny*1000)
real drift_start_topo_1m(Ny*1000)

! This program:
!   1) takes the coarse model topography and generates 1.0-m grid
!        increment topo lines;
!   2) uses those to generate the Tabler surface profiles; and
!   3) extracts the profiles at the coarse model grid cells.

! Fill the snow array with topo data.
do j = 1, Ny
  do i = 1, Nx
  
    if (isnan(dem(i,j))) goto 5 ! Exclude points outside of the domain
  
    tabler_ww(i,j) = topo_land(i,j)
    
    5 continue  ! Exclude points outside of the domain

  end do
end do

! Define the length of the j-looping, depending on whether there
!   is any rotation done to get 45 deg winds.
if (irotate_flag == 2) then
  nnx = Nx+Ny-1
else
  nnx = Nx
end if

! Required parameters.
if (irotate_flag == 2) then
  dy = 1.41
  test = amod(delta_EW/1.41,dy/1.41)
else
  dy = 1.0
  test = amod(delta_EW,dy)
end if

xmax_slope = -0.20

! This Tabler program has not been made general enough to deal
!   with delta_EW and delta_SN values that are not evenly divisible
!   by 1.0 (or 1.41 for diagonal profiles).
if (abs(test) > 1.0e10-5) then
  print *,'To generate the Tabler surfaces, delta_EW and delta_SN'
  print *,'  must be evenly divisible by 1.0 or 1.41.'
  print *,'  delta_EW = ',delta_EW
  stop
end if

! Define the number of 1-m grid cells in each model grid cell, and
!   calculate how many of these are in the entire Ny domain.
nny = nint(delta_EW/dy)
maxlines = (Ny - 1) * nny + 1

! Define the starting and ending points of the line we are going to
!   work with.  This is done to make sure we are never looking
!   outside the data for numbers to work with.
jstart = 1
if (irotate_flag == 2) then
  jend = maxlines - (32+11+10+11)
else
  jend = maxlines - (45+15+15+15)
end if

! Extract the line we are going to work with.
do i = 1, nnx

  if (irotate_flag == 2) then
    do j = 1, Ny
      II = i - j + 1
      drift_start_topo(j) = 0.0
      if (II <= 0) then
        tabler(j) = tabler(j-1)
        topo_line(j) = tabler(j-1)
      else if (II > Nx) then
        tabler(j) = tabler_ww(Nx,i-Nx+1)
        topo_line(j) = tabler_ww(Nx,i-Nx+1)
      else
        tabler(j) = tabler_ww(II,j)
        topo_line(j) = tabler_ww(II,j)
      end if
    end do
  else
    do j = 1, Ny
      tabler(j) = tabler_ww(i,j)
      topo_line(j) = topo_land(i,j)
      drift_start_topo(j) = 0.0
    end do
  end if

  ! To build the 1.0 m line, use linear interpolation between the
  !   model topo data.  Include the end point.
  do j = 1, Ny-1
    do jj = 1, nny
      jjj = (j - 1) * nny + jj
      x1 = 0.0
      x = real(jj - 1) * dy
      y2 = topo_line(j+1)
      y1 = topo_line(j)
      topo_1m(jjj) = y1 + ((y2 - y1)/delta_EW) * (x - x1)
    end do
  end do
  topo_1m((Ny - 1) * nny + 1) = topo_line(Ny)

  ! Use this topo array to be the starting point for generating the
  !   Tabler surfaces.
  do j = 1, maxlines
    tabler_1m(j) = topo_1m(j)
    drift_start_topo_1m(j) = 0.0
  end do

  ! Run the Tabler model.
  do j = jstart, jend
    if (irotate_flag == 2) then
      t1 = tabler_1m(j)
      t2 = tabler_1m(j+31)
      t3 = tabler_1m(j+31+11)
      t4 = tabler_1m(j+31+21)
      t5 = tabler_1m(j+31+32)

      x1 = (t2 - t1) / 45.0
      x2 = max((t3 - t2) / 15.0,xmax_slope)
      x3 = max((t4 - t3) / 15.0,xmax_slope)
      x4 = max((t5 - t4) / 15.0,xmax_slope)

      y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

      tabler_1m(j+32) = max(topo_1m(j+32), &
                        tabler_1m(j+31) + y * slope_adjust * dy)

      ! Create an erosion area 'erosion_dist' upwind of the tabler surface.
      npts = nint(erosion_dist/1.41)
      if (tabler_1m(j+32) /= topo_1m(j+32)) then
        do jj = 1, npts
          if (tabler_1m(j+32-jj) == topo_1m(j+32-jj)) drift_start_topo_1m(j+32-jj) = -8888.0
        end do
      end if
    else
      t1 = tabler_1m(j)
      t2 = tabler_1m(j+44)
      t3 = tabler_1m(j+44+15)
      t4 = tabler_1m(j+44+30)
      t5 = tabler_1m(j+44+45)

      x1 = (t2 - t1) / 45.0
      x2 = max((t3 - t2) / 15.0,xmax_slope)
      x3 = max((t4 - t3) / 15.0,xmax_slope)
      x4 = max((t5 - t4) / 15.0,xmax_slope)

      y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

      tabler_1m(j+45) = max(topo_1m(j+45), &
                        tabler_1m(j+44) + y * slope_adjust * dy)

      ! Create an erosion area 'erosion_dist' upwind of the tabler surface.
      npts = nint(erosion_dist/1.0)
      if (tabler_1m(j+45) /= topo_1m(j+45)) then
        do jj = 1, npts
          if (tabler_1m(j+45-jj) == topo_1m(j+45-jj)) drift_start_topo_1m(j+45-jj) = -8888.0
        end do
      end if
    end if
  end do

  ! Extract the profile at the model grid points.
  do j = 1, Ny
    jj = (j - 1) * nny + 1
    tabler(j) = tabler_1m(jj)
    drift_start_topo(j) = drift_start_topo_1m(jj)
  end do

  ! Use the 1-D arrays to fill in the 2-D tabler-surface array.
  do j = 1, Ny
    if (irotate_flag == 2) then
      II = i - j + 1
      if (II >= 1 .AND. II <= Nx) then
        tabler_ww(II,j) = tabler(j) + drift_start_topo(j)
      end if
    else
      tabler_ww(i,j) = tabler(j) + drift_start_topo(j)
    end if
  end do

end do

! Convert snow_traps back to actual snow depths instead of
!   depth plus topography.
do j = 1, Ny
  do i = 1, Nx
  
    if (isnan(dem(i,j))) goto 6 ! Exclude points outside of the domain
  
    tabler_ww(i,j) = tabler_ww(i,j) - topo_land(i,j)
    
    6 continue  ! Exclude points outside of the domain
    
  end do
end do

end subroutine tabler_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tabler_e(topo_land,tabler_ee,delta_EW, &
                    irotate_flag)

! This subroutine uses Tabler (1975) to define equilibrium profiles
!   for the topographic drift catchments, for the case of winds
!   from the EAST and SOUTHEAST.

! Tabler, R. D., 1975: Predicting profiles of snowdrifts in
!   topographic catchments.  Proceedings of the 43rd Annual Western
!   Snow Conference, San Diego, California, 87-97.

use GRID, only: &
  Nx,Ny               ! Grid dimensions
  
use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  erosion_dist,      &! Erosion area
  slope_adjust        ! Slope adjustment factor

implicit none

real, intent(in) :: &
  topo_land(Nx,Ny),  &! Topography (m)
  delta_EW            ! Grid cell size in EW direction (m)
  
real, intent(inout) :: &
  tabler_ee(Nx,Ny)    ! Tabler surface EE

integer, intent(in) :: &
  irotate_flag        ! Flag

integer i,j,jstart,jend,II,nny,nnx,jj,jjj,maxlines
integer npts
real y,xmax_slope,dy,test,x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

real tabler(Ny)
real topo_line(Ny)
real drift_start_topo(Ny)

!integer (Ny_max_tabler=Ny*1000)
real topo_1m(Ny*1000)
real tabler_1m(Ny*1000)
real drift_start_topo_1m(Ny*1000)

! This program:
!   1) takes the coarse model topography and generates 1.0-m grid
!        increment topo lines;
!   2) uses those to generate the Tabler surface profiles; and
!   3) extracts the profiles at the coarse model grid cells.

! Fill the snow array with topo data.
do j = 1, Ny
  do i = 1, Nx
  
    if (isnan(dem(i,j))) goto 7 ! Exclude points outside of the domain
  
    tabler_ee(i,j) = topo_land(i,j)
    
    7 continue  ! Exclude points outside of the domain
    
  end do
end do

! Define the length of the i-looping, depending on whether there
!   is any rotation done to get 45 deg winds.
if (irotate_flag == 2) then
  nnx = Nx+Ny-1
else
  nnx = Nx
end if

! Required parameters.
if (irotate_flag == 2) then
  dy = 1.41
  test = amod(delta_EW/1.41,dy/1.41)
else
  dy = 1.0
  test = amod(delta_EW,dy)
end if

xmax_slope = -0.20

! This Tabler program has not been made general enough to deal
!   with delta_EW and delta_SN values that are not evenly divisible
!   by 1.0 (or 1.41 for diagonal profiles).
if (abs(test) > 1.0e10-5) then
  print *,'To generate the Tabler surfaces, delta_EW and delta_SN'
  print *,'  must be evenly divisible by 1.0 or 1.41.'
  print *,'  delta_EW = ',delta_EW
  stop
end if

! Define the number of 1-m grid cells in each model grid cell, and
!   calculate how many of these are in the entire Ny domain.
nny = nint(delta_EW/dy)
maxlines = (Ny - 1) * nny + 1

! Define the starting and ending points of the line we are going to
!   work with.  This is done to make sure we are never looking
!   outside the data for numbers to work with.
jstart = maxlines
if (irotate_flag == 2) then
  jend = 1 + (32+11+10+11)
else
  jend = 1 + (45+15+15+15)
end if

! Extract the line we are going to work with.
do i = 1, nnx

  if (irotate_flag == 2) then
    do j = 1, Ny
      II = i - j + 1
      drift_start_topo(j) = 0.0
      if (II <= 0) then
        tabler(j) = tabler(j-1)
        topo_line(j) = tabler(j-1)
      else if (II > Nx) then
        tabler(j) = tabler_ee(Nx,i-Nx+1)
        topo_line(j) = tabler_ee(Nx,i-Nx+1)
      else
        tabler(j) = tabler_ee(II,j)
        topo_line(j) = tabler_ee(II,j)
      end if
    end do
  else
    do j = 1, Ny
      tabler(j) = tabler_ee(i,j)
      topo_line(j) = topo_land(i,j)
      drift_start_topo(j) = 0.0
    end do
  end if

  ! To build the 1.0 m line, use linear interpolation between the
  !   model topo data.  Include the end point.
  ! To build the 1.0 m line, use linear interpolation between the
  !   model topo data.  Include the end point.
  do j = 1, Ny-1
    do jj = 1, nny
      jjj = (j - 1) * nny + jj
      x1 = 0.0
      x = real(jj - 1) * dy
      y2 = topo_line(j+1)
      y1 = topo_line(j)
      topo_1m(jjj) = y1 + ((y2 - y1)/delta_EW) * (x - x1)
    end do
  end do
  topo_1m((Ny - 1) * nny + 1) = topo_line(Ny)

  ! Use this topo array to be the starting point for generating the
  !   Tabler surfaces.
  do j = 1, maxlines
    tabler_1m(j) = topo_1m(j)
    drift_start_topo_1m(j) = 0.0
  end do

  ! Run the Tabler model.
  do j = jstart, jend,-1
    if (irotate_flag == 2) then
      t1 = tabler_1m(j)
      t2 = tabler_1m(j-31)
      t3 = tabler_1m(j-31-11)
      t4 = tabler_1m(j-31-21)
      t5 = tabler_1m(j-31-32)

      x1 = (t2 - t1) / 45.0
      x2 = max((t3 - t2) / 15.0,xmax_slope)
      x3 = max((t4 - t3) / 15.0,xmax_slope)
      x4 = max((t5 - t4) / 15.0,xmax_slope)

      y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

      tabler_1m(j-32) = max(topo_1m(j-32), &
                        tabler_1m(j-31) + y * slope_adjust * dy)

      ! Create an erosion area 'erosion_dist' upwind of the tabler surface.
      npts = nint(erosion_dist/1.41)
      if (tabler_1m(j-32) /= topo_1m(j-32)) then
        do jj = 1, npts
          if (tabler_1m(j-32+jj) == topo_1m(j-32+jj)) drift_start_topo_1m(j-32+jj) = -8888.0
        end do
      end if
    else
      t1 = tabler_1m(j)
      t2 = tabler_1m(j-44)
      t3 = tabler_1m(j-44-15)
      t4 = tabler_1m(j-44-30)
      t5 = tabler_1m(j-44-45)

      x1 = (t2 - t1) / 45.0
      x2 = max((t3 - t2) / 15.0,xmax_slope)
      x3 = max((t4 - t3) / 15.0,xmax_slope)
      x4 = max((t5 - t4) / 15.0,xmax_slope)

      y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

      tabler_1m(j-45) = max(topo_1m(j-45), &
                        tabler_1m(j-44) + y * slope_adjust * dy)

      ! Create an erosion area 'erosion_dist' upwind of the tabler surface.
      npts = nint(erosion_dist/1.0)
      if (tabler_1m(j-45) /= topo_1m(j-45)) then
        do jj = 1, npts
          if (tabler_1m(j-45+jj) == topo_1m(j-45+jj)) drift_start_topo_1m(j-45+jj) = -8888.0
        end do
      end if
    end if
  end do

! Extract the profile at the model grid points.
  do j = 1, Ny
    jj = (j - 1) * nny + 1
    tabler(j) = tabler_1m(jj)
    drift_start_topo(j) = drift_start_topo_1m(jj)
  end do

! Use the 1-D arrays to fill in the 2-D tabler-surface array.
  do j = 1, Ny
    if (irotate_flag == 2) then
      II = i - j + 1
      if (II >= 1 .AND. II <= Nx) then
        tabler_ee(II,j) = tabler(j) + drift_start_topo(j)
      end if
    else
      tabler_ee(i,j) = tabler(j) + drift_start_topo(j)
    end if
  end do

end do

! Convert snow_traps back to actual snow depths instead of
!   depth plus topography.
do j = 1, Ny
  do i = 1, Nx
  
    if (isnan(dem(i,j))) goto 8 ! Exclude points outside of the domain
  
    tabler_ee(i,j) = tabler_ee(i,j) - topo_land(i,j)
    
    8 continue  ! Exclude points outside of the domain    
   
  end do
end do

end subroutine tabler_e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tabler_s(topo_land,tabler_ss,delta_SN, &
                    irotate_flag)

! This subroutine uses Tabler (1975) to define equilibrium profiles
!   for the topographic drift catchments, for the case of winds
!   from the SOUTH and SOUTHWEST.

! Tabler, R. D., 1975: Predicting profiles of snowdrifts in
!   topographic catchments.  Proceedings of the 43rd Annual Western
!   Snow Conference, San Diego, California, 87-97.

use GRID, only: &
  Nx,Ny               ! Grid dimensions
  
use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  erosion_dist,      &! Erosion area
  slope_adjust        ! Slope adjustment factor

implicit none

real, intent(in) :: &
  topo_land(Nx,Ny),  &! Topography (m)
  delta_SN            ! Grid cell size in SN direction (m)

real, intent(inout) :: &
  tabler_ss(Nx,Ny)    ! Tabler surface SS

integer, intent(in) :: &
  irotate_flag        ! Flag

integer i,j,istart,iend,JJ,nny,nnx,ii,iii,maxlines
integer npts
real y,xmax_slope,dx,test,x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

real tabler(Nx)
real topo_line(Nx)
real drift_start_topo(Nx)

!integer (Nx_max_tabler=Nx*1000)
real topo_1m(Nx*1000)
real tabler_1m(Nx*1000)
real drift_start_topo_1m(Nx*1000)

! This program:
!   1) takes the coarse model topography and generates 1.0-m grid
!        increment topo lines;
!   2) uses those to generate the Tabler surface profiles; and
!   3) extracts the profiles at the coarse model grid cells.

! Fill the snow array with topo data.
do j = 1, Ny
  do i = 1, Nx
  
    if (isnan(dem(i,j))) goto 9 ! Exclude points outside of the domain
  
    tabler_ss(i,j) = topo_land(i,j)
    
    9 continue  ! Exclude points outside of the domain   
    
  end do
end do

! Define the length of the j-looping, depending on whether there
!   is any rotation done to get 45 deg winds.
if (irotate_flag == 2) then
  nny = Nx+Ny-1
else
  nny = Ny
end if

! Required parameters.
if (irotate_flag == 2) then
  dx = 1.41
  test = amod(delta_SN/1.41,dx/1.41)
else
  dx = 1.0
  test = amod(delta_SN,dx)
end if

xmax_slope = -0.20

! This Tabler program has not been made general enough to deal
!   with delta_EW and delta_SN values that are not evenly divisible
!   by 1.0 (or 1.41 for diagonal profiles).
if (abs(test) > 1.0e10-5) then
  print *,'To generate the Tabler surfaces, delta_EW and delta_SN'
  print *,'  must be evenly divisible by 1.0 or 1.41.'
  print *,'  delta_SN = ',delta_SN
  stop
end if

! Define the number of 1-m grid cells in each model grid cell, and
!   calculate how many of these are in the entire Nx domain.
nnx = nint(delta_SN/dx)
maxlines = (Nx - 1) * nnx + 1

! Define the starting and ending points of the line we are going to
!   work with.  This is done to make sure we are never looking
!   outside the data for numbers to work with.
istart = 1
if (irotate_flag == 2) then
  iend = maxlines - (32+11+10+11)
else
  iend = maxlines - (45+15+15+15)
end if

! Extract the line we are going to work with.
do j = 1, nny

  if (irotate_flag == 2) then
    do i = 1, Nx
      JJ = j + i - Nx
      drift_start_topo(i) = 0.0
      if (JJ <= 0) then
        tabler(i) = tabler_ss(1-j+Nx,1)
        topo_line(i) = tabler_ss(1-j+Nx,1)
      else if (JJ > Ny) then
        tabler(i) = tabler(i-1)
        topo_line(i) = tabler(i-1)
      else
        tabler(i) = tabler_ss(i,JJ)
        topo_line(i) = tabler_ss(i,JJ)
      end if
    end do
  else
    do i = 1, Nx
      tabler(i) = tabler_ss(i,j)
      topo_line(i) = topo_land(i,j)
      drift_start_topo(i) = 0.0
    end do
  end if

  ! To build the 1.0 m line, use linear interpolation between the
  !   model topo data.  Include the end point.
  do i = 1, Nx-1
    do ii = 1, nnx
      iii = (i - 1) * nnx + ii
      x1 = 0.0
      x = real(ii - 1) * dx
      y2 = topo_line(i+1)
      y1 = topo_line(i)
      topo_1m(iii) = y1 + ((y2 - y1)/delta_SN) * (x - x1)
    end do
  end do
  topo_1m((Nx - 1) * nnx + 1) = topo_line(Nx)

  ! Use this topo array to be the starting point for generating the
  !   Tabler surfaces.
  do i = 1, maxlines
    tabler_1m(i) = topo_1m(i)
    drift_start_topo_1m(i) = 0.0
  end do

  ! Run the Tabler model.
  do i = istart, iend
    if (irotate_flag == 2) then
      t1 = tabler_1m(i)
      t2 = tabler_1m(i+31)
      t3 = tabler_1m(i+31+11)
      t4 = tabler_1m(i+31+21)
      t5 = tabler_1m(i+31+32)

      x1 = (t2 - t1) / 45.0
      x2 = max((t3 - t2) / 15.0,xmax_slope)
      x3 = max((t4 - t3) / 15.0,xmax_slope)
      x4 = max((t5 - t4) / 15.0,xmax_slope)

      y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

      tabler_1m(i+32) = max(topo_1m(i+32), &
                        tabler_1m(i+31) + y * slope_adjust * dx)

      ! Create an erosion area 'erosion_dist' upwind of the tabler surface.
      npts = nint(erosion_dist/1.41)
      if (tabler_1m(i+32) /= topo_1m(i+32)) then
        do ii = 1, npts
          if (tabler_1m(i+32-ii) == topo_1m(i+32-ii)) drift_start_topo_1m(i+32-ii) = -8888.0
        end do
      end if
    else
      t1 = tabler_1m(i)
      t2 = tabler_1m(i+44)
      t3 = tabler_1m(i+44+15)
      t4 = tabler_1m(i+44+30)
      t5 = tabler_1m(i+44+45)

      x1 = (t2 - t1) / 45.0
      x2 = max((t3 - t2) / 15.0,xmax_slope)
      x3 = max((t4 - t3) / 15.0,xmax_slope)
      x4 = max((t5 - t4) / 15.0,xmax_slope)

      y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

      tabler_1m(i+45) = max(topo_1m(i+45), &
                        tabler_1m(i+44) + y * slope_adjust * dx)

      ! Create an erosion area 'erosion_dist' upwind of the tabler surface.
      npts = nint(erosion_dist/1.0)
      if (tabler_1m(i+45) /= topo_1m(i+45)) then
        do ii = 1, npts
          if (tabler_1m(i+45-ii) == topo_1m(i+45-ii)) drift_start_topo_1m(i+45-ii) = -8888.0
        end do
      end if
    end if
  end do

  ! Extract the profile at the model grid points.
  do i = 1, Nx
    ii = (i - 1) * nnx + 1
    tabler(i) = tabler_1m(ii)
    drift_start_topo(i) = drift_start_topo_1m(ii)
  end do

  ! Use the 1-D arrays to fill in the 2-D tabler-surface array.
  do i = 1, Nx
    if (irotate_flag == 2) then
      JJ = j + i - Nx
      if (JJ >= 1 .AND. JJ <= Ny) then
        tabler_ss(i,JJ) = tabler(i) + drift_start_topo(i)
      end if
    else
      tabler_ss(i,j) = tabler(i) + drift_start_topo(i)
    end if
  end do

end do

! Convert snow_traps back to actual snow depths instead of
!   depth plus topography.
do j = 1, Ny
  do i = 1, Nx
  
    if (isnan(dem(i,j))) goto 10 ! Exclude points outside of the domain
  
    tabler_ss(i,j) = tabler_ss(i,j) - topo_land(i,j)
    
    10 continue  ! Exclude points outside of the domain   
    
  end do
end do
end subroutine tabler_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tabler_n(topo_land,tabler_nn,delta_SN, &
                    irotate_flag)

! This subroutine uses Tabler (1975) to define equilibrium profiles
!   for the topographic drift catchments, for the case of winds
!   from the NORTH and NORTHEAST.

! Tabler, R. D., 1975: Predicting profiles of snowdrifts in
!   topographic catchments.  Proceedings of the 43rd Annual Western
!   Snow Conference, San Diego, California, 87-97.

use GRID, only: &
  Nx,Ny               ! Grid dimensions
  
use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  erosion_dist,      &! Erosion area
  slope_adjust        ! Slope adjustment factor

implicit none

real, intent(in) :: &
  topo_land(Nx,Ny),  &! Topography (m)
  delta_SN            ! Grid cell size in SN direction (m)

real, intent(inout) :: &
  tabler_nn(Nx,Ny)    ! Tabler surface NN

integer, intent(in) :: &
  irotate_flag        ! Flag

integer i,j,istart,iend,JJ,nny,nnx,ii,iii,maxlines
integer npts
real y,xmax_slope,dx,test,x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

real tabler(Nx)
real topo_line(Nx)
real drift_start_topo(Nx)

!integer (Nx_max_tabler=Nx*1000)
real topo_1m(Nx*1000)
real tabler_1m(Nx*1000)
real drift_start_topo_1m(Nx*1000)

! This program:
!   1) takes the coarse model topography and generates 1.0-m grid
!        increment topo lines;
!   2) uses those to generate the Tabler surface profiles; and
!   3) extracts the profiles at the coarse model grid cells.

! Fill the snow array with topo data.
do j = 1, Ny
  do i = 1, Nx
  
    if (isnan(dem(i,j))) goto 11 ! Exclude points outside of the domain
    
    tabler_nn(i,j) = topo_land(i,j)
    
    11 continue  ! Exclude points outside of the domain
    
    end do
end do

! Define the length of the i-looping, depending on whether there
!   is any rotation done to get 45 deg winds.
if (irotate_flag == 2) then
  nny = Nx+Ny-1
else
  nny = Ny
end if

! Required parameters.
if (irotate_flag == 2) then
  dx = 1.41
  test = amod(delta_SN/1.41,dx/1.41)
else
  dx = 1.0
  test = amod(delta_SN,dx)
end if

xmax_slope = -0.20

! This Tabler program has not been made general enough to deal
!   with delta_EW and delta_SN values that are not evenly divisible
!   by 1.0 (or 1.41 for diagonal profiles).
if (abs(test) > 1.0e10-5) then
  print *,'To generate the Tabler surfaces, delta_EW and delta_SN'
  print *,'  must be evenly divisible by 1.0 or 1.41.'
  print *,'  delta_SN = ',delta_SN
  stop
end if

! Define the number of 1-m grid cells in each model grid cell, and
!   calculate how many of these are in the entire Nx domain.
nnx = nint(delta_SN/dx)
maxlines = (Nx - 1) * nnx + 1

! Define the starting and ending points of the line we are going to
!   work with.  This is done to make sure we are never looking
!   outside the data for numbers to work with.
istart = maxlines
if (irotate_flag == 2) then
  iend = 1 + (32+11+10+11)
else
  iend = 1 + (45+15+15+15)
end if

! Extract the line we are going to work with.
do j = 1, nny

  if (irotate_flag == 2) then
    do i = 1, Nx
      JJ = j + i - Nx
      drift_start_topo(i) = 0.0
      if (JJ <= 0) then
        tabler(i) = tabler_nn(1-j+Nx,1)
        topo_line(i) = tabler_nn(1-j+Nx,1)
      else if (JJ > Ny) then
        tabler(i) = tabler(i-1)
        topo_line(i) = tabler(i-1)
      else
        tabler(i) = tabler_nn(i,JJ)
        topo_line(i) = tabler_nn(i,JJ)
      end if
    end do
  else
    do i = 1, Nx
      tabler(i) = tabler_nn(i,j)
      topo_line(i) = topo_land(i,j)
      drift_start_topo(i) = 0.0
    end do
  end if

  ! To build the 1.0 m line, use linear interpolation between the
  !   model topo data.  Include the end point.
  do i = 1, Nx-1
    do ii = 1, nnx
      iii = (i - 1) * nnx + ii
      x1 = 0.0
      x = real(ii - 1) * dx
      y2 = topo_line(i+1)
      y1 = topo_line(i)
      topo_1m(iii) = y1 + ((y2 - y1)/delta_SN) * (x - x1)
    end do
  end do
  topo_1m((Nx - 1) * nnx + 1) = topo_line(Nx)

  ! Use this topo array to be the starting point for generating the
  !   Tabler surfaces.
  do i = 1, maxlines
    tabler_1m(i) = topo_1m(i)
    drift_start_topo_1m(i) = 0.0
  end do

  ! Run the Tabler model.
  do i = istart, iend,-1
    if (irotate_flag == 2) then
      t1 = tabler_1m(i)
      t2 = tabler_1m(i-31)
      t3 = tabler_1m(i-31-11)
      t4 = tabler_1m(i-31-21)
      t5 = tabler_1m(i-31-32)

      x1 = (t2 - t1) / 45.0
      x2 = max((t3 - t2) / 15.0,xmax_slope)
      x3 = max((t4 - t3) / 15.0,xmax_slope)
      x4 = max((t5 - t4) / 15.0,xmax_slope)

      y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

      tabler_1m(i-32) = max(topo_1m(i-32), &
                        tabler_1m(i-31) + y * slope_adjust * dx)

      ! Create an erosion area 'erosion_dist' upwind of the tabler surface.
      npts = nint(erosion_dist/1.41)
      if (tabler_1m(i-32) /= topo_1m(i-32)) then
        do ii = 1, npts
          if (tabler_1m(i-32+ii) == topo_1m(i-32+ii)) drift_start_topo_1m(i-32+ii) = -8888.0
        end do
      end if
    else
      t1 = tabler_1m(i)
      t2 = tabler_1m(i-44)
      t3 = tabler_1m(i-44-15)
      t4 = tabler_1m(i-44-30)
      t5 = tabler_1m(i-44-45)

      x1 = (t2 - t1) / 45.0
      x2 = max((t3 - t2) / 15.0,xmax_slope)
      x3 = max((t4 - t3) / 15.0,xmax_slope)
      x4 = max((t5 - t4) / 15.0,xmax_slope)

      y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

      tabler_1m(i-45) = max(topo_1m(i-45), &
                        tabler_1m(i-44) + y * slope_adjust * dx)

      ! Create an erosion area 'erosion_dist' upwind of the tabler surface.
      npts = nint(erosion_dist/1.0)
      if (tabler_1m(i-45) /= topo_1m(i-45)) then
        do ii = 1, npts
          if (tabler_1m(i-45+ii) == topo_1m(i-45+ii)) drift_start_topo_1m(i-45+ii) = -8888.0
        end do
      end if
    end if
  end do

  ! Extract the profile at the model grid points.
  do i = 1, Nx
    ii = (i - 1) * nnx + 1
    tabler(i) = tabler_1m(ii)
    drift_start_topo(i) = drift_start_topo_1m(ii)
  end do

  ! Use the 1-D arrays to fill in the 2-D tabler-surface array.
  do i = 1, Nx
    if (irotate_flag == 2) then
      JJ = j + i - Nx
      if (JJ >= 1 .AND. JJ <= Ny) then
        tabler_nn(i,JJ) = tabler(i) + drift_start_topo(i)
      end if
    else
      tabler_nn(i,j) = tabler(i) + drift_start_topo(i)
    end if
  end do

end do

! Convert snow_traps back to actual snow depths instead of
!   depth plus topography.
do j = 1, Ny
  do i = 1, Nx
  
    if (isnan(dem(i,j))) goto 12 ! Exclude points outside of the domain
    
    tabler_nn(i,j) = tabler_nn(i,j) - topo_land(i,j)
    
    12 continue  ! Exclude points outside of the domain  
    
  end do
end do

end subroutine tabler_n