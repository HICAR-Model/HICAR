!-----------------------------------------------------------------------
! Snow transport by wind - Liston's scheme adapted to FSM 2.0
! References: Liston and Sturm (1998), Liston et al. (2007)
! Author: Glen E. Liston
! Adapted to FSM 2.0 by Louis Quéno
! ATTENTION: y is the W->E axis, while x is the S->N
! South of (i,j): (i-1,j)
! North of (i,j): (i+1,j)
! West of (i,j): (i,j-1)
! East of (i,j): (i,j+1)
!-----------------------------------------------------------------------
subroutine SNOWTRAN3D(tabler_nn,tabler_ne,tabler_ee,tabler_se, &
                      tabler_ss,tabler_sw,tabler_ww,tabler_nw, &
                      snowdepth0,Sice0,dm_salt,dm_susp,dm_subl,dm_subgrid)

use CONSTANTS, only: &
  pi,                &! pi
  rho_wat,           &! Density of water (kg/m^3)
  Tm                  ! Melting point (K)

use DRIVING, only: &  
  Ua,                &! Wind speed (m/s)
  Udir                ! Wind direction (degrees, clockwise from N)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMMAPS, only: &
  vegsnowd_xy         ! Vegetation snow holding capacity (m)

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Tsnow               ! Snow layer temperatures (K)

use LANDUSE, only: &
  Ld,                &! Grid cell size (m)
  dem                 ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  Utau_t_flag,       &! Flag for variable threshold friction velocitiy (1.0) or constant (0.0)
  Utau_t_const,      &! Constant threshold friction velocity (m/s) used if Utau_t_flag = 0.0
  twolayer_flag,     &! Flag for soft/hard layers distinction
  rho_snow            ! Constant snow density (kg/m^3)

use CONSTANTS_SNOWTRAN3D, only: &
  wind_min            ! Minimum wind speed to compute snow transport (m/s)

implicit none

real, intent(in) :: &
  tabler_nn(Nx,Ny),  &! Tabler surfaces NN
  tabler_ss(Nx,Ny),  &! Tabler surfaces SS
  tabler_ee(Nx,Ny),  &! Tabler surfaces EE
  tabler_ww(Nx,Ny),  &! Tabler surfaces WW
  tabler_ne(Nx,Ny),  &! Tabler surfaces NE
  tabler_se(Nx,Ny),  &! Tabler surfaces SE
  tabler_sw(Nx,Ny),  &! Tabler surfaces SW
  tabler_nw(Nx,Ny)    ! Tabler surfaces NW

real, intent(inout) :: &
  snowdepth0(Nx,Ny), &! Snow depth of snowdrift accumulation, averaged over the grid cell (m)
  Sice0(Nx,Ny),      &! Ice content of snowdrift accumulation(kg/m^2)
  dm_salt(Nx,Ny),    &! SWE change due to saltation (kg/m^2)
  dm_susp(Nx,Ny),    &! SWE change due to suspension (kg/m^2)
  dm_subl(Nx,Ny),    &! SWE change due to sublimation (kg/m^2)
  dm_subgrid(Nx,Ny)   ! SWE change due to subgrid redistribution (kg/m^2)

integer :: &
  i,j                 ! Point counters

integer :: &
  index_ue(Nx,2*Ny+1), &! Wind index array E
  index_uw(Nx,2*Ny+1), &! Wind index array W
  index_vn(Ny,2*Nx+1), &! Wind index array N
  index_vs(Ny,2*Nx+1)   ! Wind index array S

real :: &
  bs_flag,           &! Blowing snow flag
  delta_WE,          &! Grid cell size in WE direction (m)
  delta_SN,          &! Grid cell size in SN direction (m)
  windspd_flag        ! Maximum wind speed on the domain (m/s)

real :: &
  Qsalt(Nx,Ny),             &! Saltation flux (kg/m/s)
  Qsalt_u(Nx,Ny),           &! x component of saltation flux (kg/m/s)
  Qsalt_v(Nx,Ny),           &! y component of saltation flux (kg/m/s)
  conc_salt(Nx,Ny),         &! Saltation-layer reference-level mass concentration (kg/m^3)
  Qsusp(Nx,Ny),             &! Suspension flux (kg/m/s)
  Qsusp_u(Nx,Ny),           &! x component of suspension flux (kg/m/s)
  Qsusp_v(Nx,Ny),           &! y component of suspension flux (kg/m/s)
  Qsubl(Nx,Ny),             &! Sublimation flux (kg/m^2/s)
  Ds_soft(Nx,Ny),           &! Soft snow thickness (m)
  snowthickness(Nx,Ny),     &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  uwind(Nx,Ny),             &! x component of wind speed (m/s)
  vwind(Nx,Ny),             &! y component of wind speed (m/s)
  Utau(Nx,Ny),              &! Friction velocity (m/s)
  Utau_t(Nx,Ny),            &! Threshold friction velocity (m/s)
  h_star(Nx,Ny),            &! Height of the saltation layer (m)
  veg_z0(Nx,Ny),            &! Vegetation roughness length (m)
  z_0(Nx,Ny)                 ! Surface roughness length (m)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialization of delta_WE and delta_SN
delta_WE = Ld(1,1)
delta_SN = Ld(1,1)

! Initialize Ds_soft and snowdepth0
Ds_soft(:,:) = 0.0
snowdepth0(:,:) = 0.0

! Initialize all fluxes at 0
Qsalt(:,:) = 0.0
Qsalt_u(:,:) = 0.0
Qsalt_v(:,:) = 0.0
dm_salt(:,:) = 0.0
conc_salt(:,:) = 0.0
Qsusp(:,:) = 0.0
Qsusp_u(:,:) = 0.0
Qsusp_v(:,:) = 0.0
dm_susp(:,:) = 0.0
Qsubl(:,:) = 0.0
dm_subl(:,:) = 0.0
dm_subgrid(:,:) = 0.0

! Initialization of uwind and vwind
uwind = Ua * cos(-Udir * pi/180.0 - pi/2.0)
vwind = Ua * sin(-Udir * pi/180.0 - pi/2.0)

! Initialization of maximum wind speed on the domain (m/s)
! Originally in micromet_code.f
! Initialization of snowthickness
windspd_flag = 0.0
do j = 1, Ny
  do i = 1, Nx

    if (isnan(dem(i,j))) goto 1 ! Exclude points outside of the domain

    windspd_flag = max(windspd_flag,Ua(i,j))
    snowthickness(i,j) = sum(Ds(:,i,j))

    1 continue ! Exclude points outside of the domain

  end do
end do

! Initialization of veg_z0 (copied from preprocess_code.f)
! Define the roughness lengths for each of the vegetation types.
! Note that this is pretty arbitrary, because these values are
! really only used when there is no blowing snow, and thus have
! no impact on the simulation except to provide a non-zero value
! for any such parts of the domain.
veg_z0(:,:) = 0.25 * vegsnowd_xy(:,:)

! If the two-layer scheme is turned on, update the thicknesses
! of the hard and soft layers.
if (twolayer_flag==1.0) then
  call compute_soft_snow(Ds_soft)
end if

! Update the threshold friction velocity.
if (Utau_t_flag==0.0) then

  Utau_t(:,:) = Utau_t_const

else if (Utau_t_flag==1.0) then

  call surface_snow(Utau_t)

end if

! Set the blowing snow flag to zero until it is clear that we will
! have blowing snow.
bs_flag = 0.0

! If the wind speed is lower that some threshold, then don't
! need to to any of the snow transport computations.
if (windspd_flag >= wind_min) then

  ! Get the wind direction indexing arrays for this particular
  ! wind event (time step).
  call getdirection(index_ue,index_uw,index_vn,index_vs,uwind,vwind)

  ! Solve for Utau and z_0 if snow is saltating, else solve assuming
  ! z_0 is known from snow depth and/or veg type, and solve for Utau.
  call solveUtau(Utau,z_0,h_star,snowthickness,veg_z0,bs_flag,Utau_t,Ds_soft)

  ! Update Ds_soft by removing too dense layers compared to U_tau
  call update_soft_snow(Utau,Ds_soft)

  ! If the blowing snow flag indicates wind transported snow
  ! somewhere within the domain (bs_flag = 1.0), run the saltation
  ! and suspension models.
  if (bs_flag == 1.0) then

    ! Solve for the saltation flux.
    call saltation(delta_WE,Utau,Utau_t,snowthickness,delta_SN, &
                   index_ue,index_uw,index_vn,index_vs,uwind,vwind,Ds_soft, &
                   Qsalt,Qsalt_u,Qsalt_v)

    ! Solve for the suspension flux.
    call suspension(Utau,z_0,h_star,Utau_t,uwind,vwind, &
                    Qsalt,conc_salt,Qsusp,Qsusp_u,Qsusp_v,Qsubl)

  end if

end if

! Compute the new snow depth due to accumulation from precipitation,
! saltation, and suspension, and the mass loss due to sublimation.
call accum(snowthickness, &
           delta_WE,delta_SN, &
           index_ue,index_uw,index_vn,index_vs, &
           bs_flag, &
           tabler_nn,tabler_ss,tabler_ee,tabler_ww, &
           tabler_ne,tabler_se,tabler_sw,tabler_nw, &
           uwind,vwind,Ds_soft, &
           Qsalt_u,Qsalt_v, &
           dm_salt, &
           Qsusp_u,Qsusp_v, &
           dm_susp, &
           Qsubl, &
           dm_subl, &
           dm_subgrid,&
           snowdepth0,Sice0)

end subroutine SNOWTRAN3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_soft_snow(Ds_soft)

use MODCONF, only: HISWET

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  histowet,          &! Historical variable for past wetting of a layer (0-1)
  Sliq                ! Liquid content of snow layers (kg/m^2)

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

implicit none

real, intent(inout) :: &
  Ds_soft(Nx,Ny)      ! Soft snow thickness (m)

integer :: &
  i,j,               &! Point counters
  k                   ! Layer counter

! Original SNOWTRAN3D determination of soft snow
if (HISWET == 0) then

  do i = 1, Nx
    do j = 1, Ny

      if (isnan(dem(i,j))) goto 2 ! Exclude points outside of the domain

      k = 1
      do while (k <= Nsnow(i,j) .AND. Sliq(k,i,j) < epsilon(Sliq))
        Ds_soft(i,j) = Ds_soft(i,j) + Ds(k,i,j)
        k = k + 1
      end do

      2 continue  ! Exclude points outside of the domain   

    end do
  end do

! New determination of soft snow based on layer wetting history and temperature
else ! HISWET == 1

  do i = 1, Nx
    do j = 1, Ny

      if (isnan(dem(i,j))) goto 3 ! Exclude points outside of the domain

      k = 1
      do while (k <= Nsnow(i,j) .AND. Sliq(k,i,j) < epsilon(Sliq) .AND. histowet(k,i,j) < 0.5)
        Ds_soft(i,j) = Ds_soft(i,j) + Ds(k,i,j)
        k = k + 1
      end do

      3 continue  ! Exclude points outside of the domain   

    end do
  end do

end if

end subroutine compute_soft_snow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_soft_snow(Utau,Ds_soft)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq                ! Liquid content of snow layers (kg/m^2)

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

implicit none

real, intent(in) :: &
  Utau(Nx,Ny)         ! Friction velocity (m/s)

real, intent(inout) :: &
  Ds_soft(Nx,Ny)      ! Soft snow thickness (m)

integer :: &
  i,j,               &! Point counters
  k                   ! Layer counter

logical :: &
  hard_layer_reached  ! A hard layer has been reached

real :: &
  Ds_soft_old,       &! Soft snow thickness before update (m)
  rho_layer,         &! Snow layer density
  Utau_t_layer        ! Threshold friction velocity of each layer (m/s)

do i = 1, Nx
  do j = 1, Ny

    if (isnan(dem(i,j))) goto 4 ! Exclude points outside of the domain

    Ds_soft_old = Ds_soft(i,j)
    Ds_soft(i,j) = 0.0
    hard_layer_reached = .FALSE.
    k = 1

    do while ((k <= Nsnow(i,j)) .AND. (Ds_soft_old - Ds_soft(i,j) > epsilon(Ds_soft)) .AND. (.NOT. hard_layer_reached))

      ! Calculate snow layer density
      rho_layer = (Sice(k,i,j) + Sliq(k,i,j)) / Ds(k,i,j) / fsnow(i,j)

      ! Calculate the snow threshold friction velocity.
      if (rho_layer <= 300.0) then
        Utau_t_layer = 0.10 * exp(0.003 * rho_layer)
      else
        Utau_t_layer = 0.005 * exp(0.013 * rho_layer)
      end if

      if (Utau(i,j) - Utau_t_layer > epsilon(Utau)) then
        Ds_soft(i,j) = Ds_soft(i,j) + Ds(k,i,j)
      else
        hard_layer_reached = .TRUE.
      end if

      k = k + 1

    end do

    4 continue  ! Exclude points outside of the domain   

  end do
end do

end subroutine update_soft_snow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine accum(snowthickness, &
                 delta_WE,delta_SN, &
                 index_ue,index_uw,index_vn,index_vs, &
                 bs_flag, &
                 tabler_nn,tabler_ss,tabler_ee,tabler_ww, &
                 tabler_ne,tabler_se,tabler_sw,tabler_nw, &
                 uwind,vwind,Ds_soft, &
                 Qsalt_u,Qsalt_v, &
                 dm_salt, &
                 Qsusp_u,Qsusp_v, &
                 dm_susp, &
                 Qsubl, &
                 dm_subl, &
                 dm_subgrid, &
                 snowdepth0,Sice0)

use CONSTANTS, only: &
  rho_wat                 ! Density of water (kg/m^3)

use DRIVING, only: &
  dt                      ! Timestep (s)

use GRID, only: &
  Nx,Ny                   ! Grid dimensions

use PARAMMAPS, only: &
  vegsnowd_xy             ! Vegetation snow holding capacity (m)

use STATE_VARIABLES, only: &
  Ds,                    &! Snow layer thicknesses (m)
  fsnow,                 &! Snow cover fraction 
  Sice,                  &! Ice content of snow layers (kg/m^2)
  Sliq,                  &! Liquid content of snow layers (kg/m^2)
  dm_tot_subl,           &! Cumulated SWE change due to sublimation (kg/m^2)
  dm_tot_trans            ! Cumulated transported SWE (kg/m^2)

use LANDUSE, only: &
  dem                     ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  subgrid_flag,          &! Flag for subgrid parameterization
  rho_snow                ! Constant snow density (kg/m^3)

implicit none

real, intent(in) :: &
  uwind(Nx,Ny),          &! WE component of wind speed (m/s)
  vwind(Nx,Ny),          &! SN component of wind speed (m/s)
  tabler_nn(Nx,Ny),      &! Tabler surfaces NN
  tabler_ss(Nx,Ny),      &! Tabler surfaces SS
  tabler_ee(Nx,Ny),      &! Tabler surfaces EE
  tabler_ww(Nx,Ny),      &! Tabler surfaces WW
  tabler_ne(Nx,Ny),      &! Tabler surfaces NE
  tabler_se(Nx,Ny),      &! Tabler surfaces SE
  tabler_sw(Nx,Ny),      &! Tabler surfaces SW
  tabler_nw(Nx,Ny)        ! Tabler surfaces NW

real, intent(in) :: &
  bs_flag,               &! Blowing snow flag
  delta_WE,              &! Grid cell size in WE direction (m)
  delta_SN                ! Grid cell size in SN direction (m)

real, intent(inout) :: &
  Qsalt_u(Nx,Ny),        &! x component of saltation flux (kg/m/s)
  Qsalt_v(Nx,Ny),        &! y component of saltation flux (kg/m/s)
  dm_salt(Nx,Ny),        &! SWE change due to saltation (kg/m^2)
  Qsusp_u(Nx,Ny),        &! x component of suspension flux (kg/m/s)
  Qsusp_v(Nx,Ny),        &! y component of suspension flux (kg/m/s)
  dm_susp(Nx,Ny),        &! SWE change due to suspension (kg/m^2)
  Qsubl(Nx,Ny),          &! Sublimation flux (kg/m^2/s)
  dm_subl(Nx,Ny),        &! SWE change due to sublimation (kg/m^2)
  dm_subgrid(Nx,Ny),     &! SWE change due to subgrid redistribution (kg/m^2)
  snowthickness(Nx,Ny),  &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  Ds_soft(Nx,Ny),        &! Soft snow thickness (m)
  snowdepth0(Nx,Ny),     &! Snow depth of snowdrift accumulation, averaged over the grid cell (m)
  Sice0(Nx,Ny)            ! Ice content of transported snow (kg/m^2)

integer, intent(in) :: &
  index_ue(Nx,2*Ny+1),   &! Wind index array E
  index_uw(Nx,2*Ny+1),   &! Wind index array W
  index_vn(Ny,2*Nx+1),   &! Wind index array N
  index_vs(Ny,2*Nx+1)     ! Wind index array S

integer :: &
  i,j                     ! Point counters

real :: &
  dm_subl_loss,          &! SWE loss due to sublimation (kg/m^2)
  dh_subl_loss,          &! Snow depth loss due to sublimation (m)
  dh_subgrid,            &! Snow depth change due to subgrid redistribution (m)
  Ds_hard,               &! Hard snow thickness (m)
  swe_loc,               &! Total snowpack SWE (kg/m^2)
  snowdmin                ! Minimum snow depth to allow transport (m)

real :: &
  snowdepth_tabler(Nx,Ny) ! Snow depth Tabler (m)

if (bs_flag == 1.0) then

  ! SALTATION
  call getnewdepth(delta_WE,delta_SN,Qsalt_u, &
                   Qsalt_v,dm_salt,index_ue,index_uw, &
                   index_vn,index_vs, &
                   snowthickness,Ds_soft,snowdepth0,Sice0)

  ! SUSPENSION
  call getnewdepth(delta_WE,delta_SN,Qsusp_u, &
                   Qsusp_v,dm_susp,index_ue,index_uw, &
                   index_vn,index_vs, &
                   snowthickness,Ds_soft,snowdepth0,Sice0)

  ! SUBLIMATION
  do i = 1, Nx
    do j = 1, Ny

      if (isnan(dem(i,j))) goto 5 ! Exclude points outside of the domain

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(vegsnowd_xy(i,j),Ds_hard)
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))

      ! Convert Qsubl to sublimated snow depth dh_subl_loss.
      dm_subl(i,j) = Qsubl(i,j) * dt * fsnow(i,j)
      dm_subl_loss = - dm_subl(i,j)
      dm_subl_loss = min(dm_subl_loss,swe_loc)
      if (dm_subl_loss > epsilon(dm_subl_loss)) then
        call HS_FROM_SWE(dm_subl_loss,dh_subl_loss,i,j)
      else
        dm_subl_loss = 0.0
        dh_subl_loss = 0.0
      end if

      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dh_subl_loss / fsnow(i,j) <= snowdmin) then
          dh_subl_loss = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dh_subl_loss > epsilon(dh_subl_loss)) then
            call SWE_FROM_HS(dh_subl_loss,dm_subl_loss,i,j)
            dm_subl(i,j) = - dm_subl_loss
          else
            dh_subl_loss = 0.0
            dm_subl(i,j) = 0.0
          end if
        end if
      else
        dm_subl(i,j) = 0.0
        dm_subl_loss = 0.0
        dh_subl_loss = 0.0
      end if

      if (dm_subl_loss > epsilon(dm_subl_loss) .and. dh_subl_loss > epsilon(dh_subl_loss)) then
        call SNOW_ABLATION(dh_subl_loss,dm_subl_loss,i,j)
      end if

      ! Update the snow layer thicknesses
      snowthickness(i,j) = sum(Ds(:,i,j))

      5 continue  ! Exclude points outside of the domain

    end do
  end do

  ! Run the subgrid parameterization to account for unrealistic
  ! snow accumulation spikes.
  if (subgrid_flag == 1.0) then

    ! Save a copy of the snow distribution to be used to calculate the
    ! snow distribution changes resulting from the subgrid
    ! redistribution.
    snowdepth_tabler(i,j) = snowthickness(i,j) * fsnow(i,j)

    ! Do the Tabler corrections
    call subgrid_1(snowdepth_tabler, &
                   index_ue,index_uw,index_vn,index_vs, &
                   tabler_nn,tabler_ss,tabler_ee,tabler_ww, &
                   tabler_ne,tabler_se,tabler_sw,tabler_nw,uwind,vwind)

    ! Calculate the snow depth resulting from the subgrid
    ! redistribution.
    do i = 1, Nx
      do j = 1, Ny

        if (isnan(dem(i,j))) goto 6 ! Exclude points outside of the domain

        dh_subgrid = snowdepth_tabler(i,j) - snowthickness(i,j) * fsnow(i,j)
        dm_subgrid(i,j) = dh_subgrid * rho_snow

        if (-dm_subgrid(i,j) > epsilon(-dm_subgrid) .and. -dh_subgrid > epsilon(-dh_subgrid)) then

          call SNOW_ABLATION(-dh_subgrid,-dm_subgrid(i,j),i,j)

        ! Net mass gain for this grid cell at this time step.
        else if (dm_subgrid(i,j) > epsilon(dm_subgrid) .and. dh_subgrid > epsilon(dh_subgrid)) then

          ! Add to the existing top layer.
          Sice0(i,j) = Sice0(i,j) + dm_subgrid(i,j)
          snowdepth0(i,j) = snowdepth0(i,j) + dh_subgrid

        end if

        ! Update the snow layer thicknesses
        snowthickness(i,j) = sum(Ds(:,i,j))

        6 continue  ! Exclude points outside of the domain

      end do
    end do

  else

    dm_subgrid(:,:) = 0.0

  end if

end if

! Update cumulated sublimation and transport
do i = 1, Nx
  do j = 1, Ny

    if (isnan(dem(i,j))) goto 7 ! Exclude points outside of the domain

    ! Fill summing arrays of the sublimation and transport quantities.
    dm_tot_subl(i,j) = dm_tot_subl(i,j) + dm_subl(i,j)
    dm_tot_trans(i,j) = dm_tot_trans(i,j) + dm_salt(i,j) + &
                        dm_susp(i,j) + dm_subgrid(i,j)

    7 continue  ! Exclude points outside of the domain

  end do
end do

end subroutine accum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine suspension(Utau,z_0,h_star,Utau_t,uwind,vwind, &
                      Qsalt,conc_salt,Qsusp,Qsusp_u,Qsusp_v,Qsubl)

use CONSTANTS, only: &
  vkman               ! Von Karman constant

use DRIVING, only: &
  Ta,                &! Air temperature (K)
  RH                  ! Relative humidity (%)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

use CONSTANTS_SNOWTRAN3D, only: &
  Up_const,          &! Constant coefficient for calculation of U_p
  dz_susp,           &! dz in suspension layer
  ztop_susp,         &! Height of the top of suspension layer (m)
  fall_vel,          &! Particle-settling velocity, s in L&S 1998 eq. 12 (m/s)
  Ur_const            ! Constant coefficient beta in L&S 1998, eq. 13

implicit none

real, intent(in) :: &
  Qsalt(Nx,Ny),      &! Saltation flux (kg/m/s)
  uwind(Nx,Ny),      &! x component of wind speed (m/s)
  vwind(Nx,Ny),      &! y component of wind speed (m/s)
  Utau_t(Nx,Ny),     &! Threshold friction velocity (m/s)
  Utau(Nx,Ny),       &! Friction velocity (m/s)
  z_0(Nx,Ny)          ! Surface roughness length (m)

real, intent(inout) :: &
  conc_salt(Nx,Ny),  &! Saltation-layer reference-level mass concentration (kg/m^3)
  Qsusp(Nx,Ny),      &! Suspension flux (kg/m/s)
  Qsusp_u(Nx,Ny),    &! x component of suspension flux (kg/m/s)
  Qsusp_v(Nx,Ny),    &! y component of suspension flux (kg/m/s)
  Qsubl(Nx,Ny)        ! Sublimation flux (kg/m^2/s)

real, intent(inout) :: &
  h_star(Nx,Ny)       ! Height of the saltation layer (m)

integer :: &
  i,j,               &! Point counters
  nzsteps,           &! Number of vertical steps
  iz                  ! Vertical steps counter

real :: &
  conc,              &! Concentration of the suspended snow at height z (kg/m^3)
  Utau_fallvel,      &! Auxiliary variable Utau/fall_vel
  prd,               &! Auxiliary variable for calculations
  U_p,               &! Horizontal particle velocity within the saltation layer (m/s)
  U_r,               &! Wind speed at height h_star (m/s)
  phistar_Cr,        &! Left member of L&S 1998 eq. 13
  V_susp,            &! Sublimation loss rate coefficient in the suspension layer (s^-1)
  V_salt,            &! Sublimation loss rate coefficient in the saltation layer (s^-1)
  z                   ! Height (m)

! Compute the mass concentration of suspended snow according to
! Kind (1992).

do i = 1, Nx
  do j = 1, Ny

    if (isnan(dem(i,j))) goto 8 ! Exclude points outside of the domain

    if (Qsalt(i,j) > epsilon(Qsalt)) then
      Utau_fallvel = Utau(i,j) / fall_vel
      if (h_star(i,j) == z_0(i,j)) h_star(i,j) = 2.0 * z_0(i,j)
      U_r = Utau(i,j)/vkman * log(h_star(i,j)/z_0(i,j))
      phistar_Cr = Utau(i,j)/U_r * Ur_const
      prd = phistar_Cr * Utau_fallvel
      U_p = Up_const * Utau_t(i,j)

      ! Compute the concentration in the saltation layer (kg/m^3).
      conc_salt(i,j) = Qsalt(i,j) / (h_star(i,j) * U_p)

      nzsteps = int((ztop_susp - h_star(i,j)) / dz_susp)

      Qsusp(i,j) = 0.0
      Qsubl(i,j) = 0.0

      do iz = 1, nzsteps
        z = h_star(i,j) + 0.5 * dz_susp + real(iz - 1) * dz_susp

        ! Compute the concentration of the suspended snow at height z.
        conc = conc_salt(i,j) * ((prd + 1.0) * &
               (z/h_star(i,j))**((-fall_vel)/(vkman*Utau(i,j))) - &
               prd)
        conc = max(conc,0.0)

        ! Only do the integration if the concentration is non-zero.
        if (conc > epsilon(conc)) then

          ! Compute the sublimation due to suspension.
          call getsublim(z,RH(i,j),Ta(i,j),Utau(i,j), &
                         z_0(i,j),V_susp,V_salt,Utau_t(i,j),1.0)

          ! Perform the quadrature (summation), without the constants.
          if (z == z_0(i,j)) z = 1.2 * z_0(i,j)
          Qsusp(i,j) = Qsusp(i,j) + conc * log(z/z_0(i,j)) * dz_susp
          Qsubl(i,j) = Qsubl(i,j) + conc * V_susp * dz_susp

        end if

      end do

      ! Finish the quadratures.
      ! Include the constants for Qsusp.
      Qsusp(i,j) = Utau(i,j) / vkman * Qsusp(i,j)

      ! Include the sublimation contribution due to saltation.
      z = h_star(i,j) / 2.0
      call getsublim(z,RH(i,j),Ta(i,j),Utau(i,j), &
                     z_0(i,j),V_susp,V_salt,Utau_t(i,j),0.0)

      Qsubl(i,j) = Qsubl(i,j) + V_salt * conc_salt(i,j) * h_star(i,j)

    else
      conc_salt(i,j) = 0.0
      Qsusp(i,j) = 0.0
      Qsubl(i,j) = 0.0
    end if

    8 continue  ! Exclude points outside of the domain

  end do
end do

! Separate the east-west and the north-south suspended transport
! components; the vector sum should equal Qsusp.
do i = 1, Nx
  do j = 1, Ny

    if (isnan(dem(i,j))) goto 9 ! Exclude points outside of the domain

    Qsusp_u(i,j) = Qsusp(i,j) * abs(uwind(i,j)) / &
                   sqrt(uwind(i,j)**2 + vwind(i,j)**2)
    Qsusp_v(i,j) = Qsusp(i,j) * abs(vwind(i,j)) / &
                   sqrt(uwind(i,j)**2 + vwind(i,j)**2)

    9 continue  ! Exclude points outside of the domain

  end do
end do

end subroutine suspension

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine saltation(delta_WE,Utau,Utau_t,snowthickness,delta_SN, &
                     index_ue,index_uw,index_vn,index_vs,uwind,vwind,Ds_soft, &
                     Qsalt,Qsalt_u,Qsalt_v)

use CONSTANTS, only: &
  grav,                 &! Acceleration due to gravity (m/s^2)
  rho_air                ! Density of air (kg/m^3)

use GRID, only: &
  Nx,Ny                  ! Grid dimensions

use PARAMMAPS, only: &
  vegsnowd_xy            ! Vegetation snow holding capacity (m)

use LANDUSE, only: &
  dem                    ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  bc_flag                ! Boundary condition flag

use CONSTANTS_SNOWTRAN3D, only: &
  fetch,                &! Equilibrium fetch distance (m)
  xmu                    ! Scaling constant for non-equilibrium saltation transport

implicit none

real, intent(in) :: &
  delta_WE,             &! Grid cell size in WE direction (m)
  delta_SN               ! Grid cell size in SN direction (m)

real, intent(in) :: &
  snowthickness(Nx,Ny), &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  uwind(Nx,Ny),         &! x component of wind speed (m/s)
  vwind(Nx,Ny),         &! y component of wind speed (m/s)
  Utau_t(Nx,Ny),        &! Threshold friction velocity (m/s)
  Utau(Nx,Ny),          &! Friction velocity (m/s)
  Ds_soft(Nx,Ny)         ! Soft snow thickness (m)

real, intent(inout) :: &
  Qsalt(Nx,Ny),         &! Saltation flux (kg/m/s)
  Qsalt_u(Nx,Ny),       &! x component of saltation flux (kg/m/s)
  Qsalt_v(Nx,Ny)         ! y component of saltation flux (kg/m/s)

integer, intent(in) :: &
  index_ue(Nx,2*Ny+1),    &! Wind index array E
  index_uw(Nx,2*Ny+1),    &! Wind index array W
  index_vn(Ny,2*Nx+1),    &! Wind index array N
  index_vs(Ny,2*Nx+1)      ! Wind index array S

integer :: &
  i,j,k,                &! Point counters
  istart,iend,          &! Point counters boundaries
  jstart,jend            ! Point counters boundaries

real :: &
  Qsalt_max(Nx,Ny),     &! Maximum possible saltation flux (kg/m/s)
  Qsalt_maxu(Nx,Ny),    &! x component of Qsalt_max (kg/m/s)
  Qsalt_maxv(Nx,Ny)      ! y component of Qsalt_max (kg/m/s)

real :: &
  blowby,               &! Blowby parameter (see description below)
  dUtau,                &! Utau difference (m/s)
  scale_EW,             &! Scaling coefficient for Eqn. 9 in L&S 1998
  scale_NS               ! Scaling coefficient for Eqn. 9 in L&S 1998

! The blowby parameter is implemented to account for the erosion
! of the tops of deep snow accumulations.  It corrects a
! deficiency in the du*/dx* < 0 formulation.  It is a number that
! should range from 0 to 1.0, and represents the fraction of the
! upwind saltation flux that is transfered farther downwind into
! the next grid cell.  So, the bigger the number, the less
! peaked the drift accumulation profile is.  blowby = 0.0 is the
! original model.  I am now using the Tabler surfaces to do the
! same kind of thing, so here I hard-code the parameter as in the
! original model.
blowby = 0.0

! Compute the maximum possible saltation flux, assuming that
! an abundance of snow is available at the surface.
do i = 1, Nx
  do j = 1, Ny

    if (isnan(dem(i,j))) goto 10 ! Exclude points outside of the domain

    ! For a given wind speed, find Qsalt_max.
    Qsalt_max(i,j) = 0.68 * rho_air / grav * &
                     Utau_t(i,j) / Utau(i,j) * (Utau(i,j)**2 - Utau_t(i,j)**2)
    Qsalt_max(i,j) = max(Qsalt_max(i,j),0.0)

    ! Now weight the max saltation flux for the u and v wind
    ! components, where the vector sum should equal Qsalt_max.
    Qsalt_maxu(i,j) = Qsalt_max(i,j) * abs(uwind(i,j)) / &
                      sqrt(uwind(i,j)**2 + vwind(i,j)**2)
    Qsalt_maxv(i,j) = Qsalt_max(i,j) * abs(vwind(i,j)) / &
                      sqrt(uwind(i,j)**2 + vwind(i,j)**2)

    10 continue  ! Exclude points outside of the domain

  end do
end do

! Define an upwind boundary condition.  If bc_flag = 1.0 then it is
! assumed that the inflow saltation flux has reached steady state.
! If bc_flag = 0.0 then the saltation flux is assumed to be zero.
! The boundary condition is implemented by initializing the arrays
! to Qsalt_max, and since upwind boundaries are not called in
! the Qsalt computation, they stay in effect for the future 
! accumulation/erosion computation.
if (bc_flag == 0.0) then
  do i = 1, Nx
    do j = 1, Ny

      if (isnan(dem(i,j))) goto 11 ! Exclude points outside of the domain

      ! Zero incoming flux at the boundaries.
      Qsalt_u(i,j) = 0.0
      Qsalt_v(i,j) = 0.0

      11 continue  ! Exclude points outside of the domain

    end do
  end do
else if (bc_flag == 1.0) then
  do i = 1, Nx
    do j = 1, Ny

      if (isnan(dem(i,j))) goto 12 ! Exclude points outside of the domain

      ! Steady-state (maximum) incoming flux at the boundaries.
      Qsalt_u(i,j) = Qsalt_maxu(i,j)
      Qsalt_v(i,j) = Qsalt_maxv(i,j)

      12 continue  ! Exclude points outside of the domain

    end do
  end do
end if

! Define the scaling coefficients for Eqn. 9 in L&S 1998. Don't
! let them be greater than 1.0 or you will make more snow than
! there was before.
scale_EW =  xmu * delta_WE / fetch
scale_EW = min(1.0,scale_EW)
scale_NS =  xmu * delta_SN / fetch
scale_NS = min(1.0,scale_NS)

! Consider WESTERLY winds.
do i = 1, Nx
  do k = 1, index_uw(i,1)
    jstart = index_uw(i,k*2)+1
    jend = index_uw(i,k*2+1)
    do j = jstart, jend
      dUtau = Utau(i,j) - Utau(i,j-1)
      if (dUtau >= epsilon(dUtau)) then
        Qsalt_u(i,j) = Qsalt_u(i,j-1) + scale_EW * &
                       (Qsalt_maxu(i,j) - Qsalt_u(i,j-1))
      else
!        Qsalt_u(i,j) = min(Qsalt_u(i,j-1),Qsalt_maxu(i,j))
        if (Qsalt_u(i,j-1) < Qsalt_maxu(i,j)) then
          Qsalt_u(i,j) = Qsalt_u(i,j-1)
        else
          Qsalt_u(i,j) = max(blowby*Qsalt_u(i,j-1),Qsalt_maxu(i,j))
        end if
      end if
    end do
  end do
end do

! Consider EASTERLY winds.
do i = 1, Nx
  do k = 1, index_ue(i,1)
    jend = index_ue(i,k*2)
    jstart = index_ue(i,k*2+1)-1
    do j = jstart, jend,-1
      dUtau = Utau(i,j) - Utau(i,j+1)
      if (dUtau >= epsilon(dUtau)) then
        Qsalt_u(i,j) = Qsalt_u(i,j+1) + scale_EW * &
                       (Qsalt_maxu(i,j) - Qsalt_u(i,j+1))
      else
!        Qsalt_u(i,j) = min(Qsalt_u(i,j+1),Qsalt_maxu(i,j))
        if (Qsalt_u(i,j+1) < Qsalt_maxu(i,j)) then
          Qsalt_u(i,j) = Qsalt_u(i,j+1)
        else
          Qsalt_u(i,j) = max(blowby*Qsalt_u(i,j+1),Qsalt_maxu(i,j))
        end if
      end if
    end do
  end do
end do

! Consider SOUTHERLY winds.
do j = 1, Ny
  do k = 1, index_vs(j,1)
    istart = index_vs(j,k*2)+1
    iend = index_vs(j,k*2+1)
    do i = istart, iend
      dUtau = Utau(i,j) - Utau(i-1,j)
      if (dUtau >= epsilon(dUtau)) then
        Qsalt_v(i,j) = Qsalt_v(i-1,j) + scale_NS * &
                       (Qsalt_maxv(i,j) - Qsalt_v(i-1,j))
      else
!        Qsalt_v(i,j) = min(Qsalt_v(i-1,j),Qsalt_maxv(i,j))
        if (Qsalt_v(i-1,j) < Qsalt_maxv(i,j)) then
          Qsalt_v(i,j) = Qsalt_v(i-1,j)
        else
          Qsalt_v(i,j) = max(blowby*Qsalt_v(i-1,j),Qsalt_maxv(i,j))
        end if
      end if
    end do
  end do
end do

! Consider NORTHERLY winds.
do j = 1, Ny
  do k = 1, index_vn(j,1)
    iend = index_vn(j,k*2)
    istart = index_vn(j,k*2+1)-1
    do i = istart, iend,-1
      dUtau = Utau(i,j) - Utau(i+1,j)
      if (dUtau >= epsilon(dUtau)) then
        Qsalt_v(i,j) = Qsalt_v(i+1,j) + scale_NS * &
                       (Qsalt_maxv(i,j) - Qsalt_v(i+1,j))
      else
!        Qsalt_v(i,j) = min(Qsalt_v(i+1,j),Qsalt_maxv(i,j))
        if (Qsalt_v(i+1,j) < Qsalt_maxv(i,j)) then
          Qsalt_v(i,j) = Qsalt_v(i+1,j)
        else
          Qsalt_v(i,j) = max(blowby*Qsalt_v(i+1,j),Qsalt_maxv(i,j))
        end if
      end if
    end do
  end do
end do

! Combine the u and v components to yield the total saltation flux
! at each grid cell.
do i = 1, Nx
  do j = 1, Ny

    if (isnan(dem(i,j))) goto 13 ! Exclude points outside of the domain

    Qsalt(i,j) = Qsalt_u(i,j) + Qsalt_v(i,j)

    13 continue  ! Exclude points outside of the domain

  end do
end do

! Adjust Qsalt to account for the availablity of snow for transport;
! taking into consideration whether there is snow on the ground,
! the holding depth of the vegetation, etc..
do i = 1, Nx
  do j = 1, Ny

    if (isnan(dem(i,j))) goto 14 ! Exclude points outside of the domain

    if (snowthickness(i,j) <= vegsnowd_xy(i,j)) then
      Qsalt(i,j) = 0.0
    end if
    if (Ds_soft(i,j) <= epsilon(Ds_soft)) then
      Qsalt(i,j) = 0.0
    end if

    14 continue  ! Exclude points outside of the domain

  end do
end do

end subroutine saltation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solveUtau(Utau,z_0,h_star,snowthickness,veg_z0,bs_flag,Utau_t,Ds_soft)

use CONSTANTS, only: &
  grav,                 &! Acceleration due to gravity (m/s^2)
  vkman                  ! Von Karman constant

use DRIVING, only: &
  Ua,                   &! Wind speed (m/s)
  zU                     ! Wind speed measurement height (m)

use GRID, only: &
  Nx,Ny                  ! Grid dimensions

use PARAMETERS, only : &
  z0sn                   ! Snow roughness length (m)

use PARAMMAPS, only: &
  vegsnowd_xy            ! Vegetation snow holding capacity (m)

use LANDUSE, only: &
  dem                    ! Terrain elevation (m)

use CONSTANTS_SNOWTRAN3D, only: &
  C_z,                  &! Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500
  h_const                ! Coefficient 1.6 in Liston and Sturm (1998) eq. 14 p. 501

implicit none

real, intent(in) :: &
  snowthickness(Nx,Ny), &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  Utau_t(Nx,Ny),        &! Threshold friction velocity (m/s)
  veg_z0(Nx,Ny),        &! Vegetation roughness length (m)
  Ds_soft(Nx,Ny)         ! Soft snow thickness (m)

real, intent(inout) :: &
  bs_flag                ! Blowing snow flag

real, intent(inout) :: &
  Utau(Nx,Ny)            ! Friction velocity (m/s)

real, intent(out) :: &  
  h_star(Nx,Ny),        &! Height of the saltation layer (m)
  z_0(Nx,Ny)             ! Surface roughness length (m)

integer :: &
  i,j                    ! Point counters

real :: &
  threshold,            &! Threshold for u* (m/s)
  threshold_flag,       &! Flag if u* is above threshold
  guess,                &! Initial guess for Utau
  sfrac,                &! Depth-fraction of vegetation covered by snow
  Utautmp,              &! Temporary Utau variable (m/s)
  windtmp,              &! Temporary wind speed variable (m/s)
  wind_max,             &! Maximum wind speed (m/s)
  z_0_tmp                ! Temporary z_0 variable (m)

! Initially set the blowing snow flag to no blowing snow
! (bs_flag = 0.0).  Then, if snow is found to blow in any
! domain grid cell, set the flag to on (bs_flag = 1.0).
bs_flag = 0.0

! Build the Utau array.
guess = 0.1
do i = 1, Nx
  do j = 1, Ny

    if (isnan(dem(i,j))) goto 15 ! Exclude points outside of the domain

    ! Determine whether snow is saltating (this influences how Utau
    ! and z_0 are computed).
    if (snowthickness(i,j) <= vegsnowd_xy(i,j)) then

      ! Saltation will not occur.
      sfrac = snowthickness(i,j) / max(vegsnowd_xy(i,j),veg_z0(i,j)) ! Eq. 3, LS 1998
      z_0(i,j) = sfrac * z0sn + (1.0 - sfrac) * veg_z0(i,j) ! Eq. 4, LS 1998
      z_0_tmp = min(0.25*zU,z_0(i,j))
      Utau(i,j) = Ua(i,j) * vkman / log(zU/z_0_tmp)
      h_star(i,j) = z_0(i,j) * h_const / C_z

    else if (Ds_soft(i,j) <= epsilon(Ds_soft)) then

      ! Saltation will not occur.
      z_0(i,j) = z0sn
      Utau(i,j) = Ua(i,j) * vkman / log(zU/z_0(i,j))
      h_star(i,j) = z_0(i,j)

    else

      ! Saltation may occur.  Test for that possibility by assuming that
      ! saltation is present, solving for Utau and z_0, and comparing
      ! whether Utau exceeds Utau_t.  If it does not, set z_0 to that
      ! of snow and recompute Utau.

      ! To help insure that the iteration converges, set the minimum
      ! wind speed to be 1.0 m/s, and the maximum wind speed to be
      ! 30 m/s at 10-m height.
      windtmp = max(1.0,Ua(i,j))
      wind_max = 30.0 * log(zU/z0sn)/log(10.0/z0sn)
      windtmp = min(windtmp,wind_max) 

      ! For u* over 0.6, use the relation z0 = 0.00734 u* - 0.0022,
      ! instead of Equation (5) in Liston and Sturm (1998).  Note that
      ! for windspeeds greater than about 35 m/s this will have to be
      ! modified for the solution algorithm to converge (because the
      ! roughness length will start to be higher than the obs height!).
      threshold = 0.6/vkman * log(zU/0.0022)
      if (windtmp <= threshold) then
        threshold_flag = 1.0
      else
        threshold_flag = 2.0
      end if

      call solve1(Utautmp,guess,windtmp, &
                  threshold_flag)

      if (Utautmp > Utau_t(i,j)) then

        ! We have saltation.
        Utau(i,j) = Utautmp
        z_0(i,j) = C_z * Utau(i,j)**2 / (2.0 * grav)
        h_star(i,j) = h_const * Utau(i,j)**2 / (2.0 * grav)
        bs_flag = 1.0

      else

        ! We do not have saltation, but the vegetation is covered by snow.
        ! Because we have determined that we do not have saltation, make
        ! sure Utau does not exceed Utau_t.
        z_0(i,j) = z0sn
        Utau(i,j) = Ua(i,j) * vkman / log(zU/z_0(i,j))
        Utau(i,j) = min(Utau(i,j),Utau_t(i,j))
        h_star(i,j) = z_0(i,j) * h_const / C_z

      end if
    end if

    15 continue  ! Exclude points outside of the domain

  end do
end do

end subroutine solveUtau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solve1(xnew,guess,windtmp,threshold_flag)

use CONSTANTS, only: &
  grav,              &! Acceleration due to gravity (m/s^2)
  vkman               ! Von Karman constant

use DRIVING, only: &
  zU                  ! Wind speed measurement height (m)

use CONSTANTS_SNOWTRAN3D, only: &
  C_z                 ! Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500

implicit none

real, intent(in) :: &
  guess,             &! Initial guess for Utau (m/s)
  windtmp,           &! Wind speed (m/s)
  threshold_flag      ! Flag if u* is above threshold

real, intent(inout) :: &  
  xnew                ! Solution for Utau (m/s)

integer :: &
  i,                 &! Point counter
  maxiter             ! Number of iterations

real :: &
  tol,               &! Convergence threshold
  old,               &! Temporary xnew
  fprime,            &! Temporary function 1
  funct               ! Temporary function 2

tol = 1.0e-3
maxiter = 20
old = guess

if (threshold_flag == 1.0) then

  do i = 1, maxiter
    fprime = - 1.0 + 2.0 / old * windtmp * vkman * &
            (log(zU) - log(C_z/(2.0*grav)) - 2.0*log(old))**(-2)
    funct = - old + windtmp * vkman * &
           (log(zU) - log(C_z/(2.0*grav)) - 2.0*log(old))**(-1)
    xnew = old - funct/fprime
    if (abs(xnew - old) < tol) return
    old = xnew
  end do

else if (threshold_flag == 2.0) then

  old = 0.6
  do i = 1, maxiter
    fprime = - 1.0 + windtmp * vkman * 0.00734 / (0.00734 * old - 0.0022) * &
            (log(zU) - log(0.00734 * old - 0.0022))**(-2)
    funct = - old + windtmp * vkman * &
           (log(zU) - log(0.00734 * old - 0.0022))**(-1)
    xnew = old - funct/fprime
    if (abs(xnew - old) < tol) return
    old = xnew
  end do

end if

write(*,*) 'max iteration exceeded when solving for Utau, Utau=',old

end subroutine solve1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getsublim(z,RH,Ta,Utau,z_0,V_susp,V_salt,Utau_t,flag)

use CONSTANTS, only: &
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  Ls,                &! Latent heat of sublimation (J/kg)
  pi,                &! pi
  Rair,              &! Gas constant for air (J/K/kg)
  Runi,              &! Universal gas constant (J/kmol/K)
  rho_ice,           &! Density of ice (kg/m^3)
  Tm,                &! Melting point (K)
  vkman,             &! Von Karman constant
  xM,                &! Molecular weight of water (kg/kmol)
  visc_air            ! Kinematic viscosity of air (m^2/s)

use DRIVING, only: &
  zRH                 ! Relative humidity measurement height (m)

implicit none

real, intent(in) :: &
  z,                 &! Height (m)
  RH,                &! Relative humidity (%)
  Ta,                &! Air temperature (K)
  Utau_t,            &! Threshold friction velocity (m/s)
  Utau,              &! Friction velocity (m/s)
  flag,              &! 1: sublimation due to suspension, 0: sublimation due to saltation
  z_0                 ! Surface roughness length (m)

real, intent(out) :: &
  V_susp,            &! Sublimation loss rate coefficient in the suspension layer (s^-1)
  V_salt              ! Sublimation loss rate coefficient in the saltation layer (s^-1)

real :: &
  D,                 &! Diffusivity of water vapor in the atmosphere (m^2/s)
  rho_sat,           &! Saturation density of water vapour (kg/m^3)
  rh_offset,         &! Relative humidity offset (fraction 0->1)
  sigma,             &! Atmospheric undersaturation of water vapour with respect to ice at height z
  alpha,             &! Coefficient alpha defined in L&S (1998), eq. A-5 p. 515 
  rbar_r,            &! Mean radius of snow particles at height z (m)
  xmbar,             &! Mean particle mass at height z (m)
  rbar,              &! Radius of a snow particle of mean particle mass xmbar at height z (m)
  u_z,               &! Wind speed at height z (m)
  x_r,               &! Fluctuating velocity component (m/s), defined in L&S (1998), eq. A-18 p. 516
  wbar,              &! Mean terminal fall velocity of suspended snow (m/s), L&S (1998), eq. A-17 p. 516
  V_r,               &! Ventilation velocity for turbulent suspension (m/s)
  V_rsalt,           &! Ventilation velocity for saltation (m/s)
  xN_r,              &! Particle Reynolds number Re
  xNu,               &! Nusselt number
  xSh,               &! Sherwood number
  tmp1,              &! Auxiliary variable 1
  tmp2,              &! Auxiliary variable 2
  top,               &! Top of eq. A-6, p. 516, L&S (1998)
  bottom              ! Bottom of eq. A-6, p. 516, L&S (1998)

D = 2.06e-5 * (Ta/273.)**(1.75)
! rho_sat = 0.622 * 10.0**(11.40 - 2353./Ta) / (Rair * Ta)
rho_sat = 0.622 / (Rair * Ta) * 610.78 * exp(21.875 * (Ta - Tm) / (Ta - 7.66))

! Assume that the rh varies according to a modification to 
! Pomeroy's humidity variation with height equation.
rh_offset = 1.0 - 0.027 * log(zRH)
sigma = (0.01 * RH - 1.0) * (rh_offset + 0.027 * log(z))
sigma = min(0.0,sigma)
sigma = max(-1.0,sigma)

alpha = 4.08 + 12.6 * z
rbar_r = 4.6e-5 * z**(-0.258)
xmbar = 4.0/3.0 * pi * rho_ice * rbar_r**3 * (1.0 + 3.0/alpha + 2.0/alpha**2)
rbar = ((3.0 * xmbar) / (4.0 * pi * rho_ice))**(0.33)
u_z = Utau/vkman * log(z/z_0)
x_r = 0.005 * u_z**(1.36)
wbar = 1.1e7 * rbar**(1.8)

if (flag == 1.0) then

  ! Compute the sublimation loss rate coefficient for the suspension
  ! layer.
  V_r = wbar + 3.0 * x_r * cos(pi/4.0)
  xN_r = 2.0 * rbar * V_r / visc_air
  xNu = 1.79 + 0.606 * xN_r**(0.5)
  xSh = xNu
  tmp1 = (Ls * xM)/(Runi * Ta) - 1.0
  tmp2 = hcon_air * Ta * xNu
  top = 2.0 * pi * rbar * sigma
  bottom = Ls/tmp2 * tmp1 + 1.0/(D * rho_sat * xSh)
  V_susp = (top/bottom)/xmbar
  V_salt = 0.0

else if (flag == 0.0) then

  ! Compute the sublimation loss rate coefficient for the saltation
  ! layer.
  V_rsalt = 0.68 * Utau + 2.3 * Utau_t
  xN_r = 2.0 * rbar * V_rsalt / visc_air
  xNu = 1.79 + 0.606 * xN_r**(0.5)
  xSh = xNu
  tmp1 = (Ls * xM)/(Runi * Ta) - 1.0
  tmp2 = hcon_air * Ta * xNu
  top = 2.0 * pi * rbar * sigma
  bottom = Ls/tmp2 * tmp1 + 1.0/(D * rho_sat * xSh)
  V_salt = (top/bottom)/xmbar
  V_susp = 0.0

end if

end subroutine getsublim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getdirection(index_ue,index_uw,index_vn,index_vs,uwind,vwind)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

implicit none

integer, intent(inout) :: &
  index_ue(Nx,2*Ny+1), &! Wind index array E
  index_uw(Nx,2*Ny+1), &! Wind index array W
  index_vn(Ny,2*Nx+1), &! Wind index array N
  index_vs(Ny,2*Nx+1)   ! Wind index array S

real, intent(in) :: &
  uwind(Nx,Ny),      &! x component of wind speed (m/s)
  vwind(Nx,Ny)        ! y component of wind speed (m/s)

integer :: &
  i,j,               &! Point counters
  npairs              ! ?

real :: &
  sign1,sign2         ! Sign indicators of wind

! Index whether the winds are blowing east or west.  The first
! column of the index array is the number of pairs of begining
! and ending array index of blocks of wind running in the same
! direction.

! Sweep looking for WESTERLY winds, looking for positive numbers.
do i = 1, Nx

  if (uwind(i,1) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 > 0.0) then
    npairs = 1
    index_uw(i,2) =  1
  else
    npairs = 0
  end if

  do j = 2, Ny

    if (isnan(dem(i,j))) goto 16 ! Exclude points outside of the domain

    if (uwind(i,j-1) <= 0.0) then
      sign1 = -1.0
    else
      sign1 = 1.0
    end if
    if (uwind(i,j) <= 0.0) then
      sign2 = -1.0
    else
      sign2 = 1.0
    end if

    if (sign2 /= sign1) then
      ! We have a sign change.
      if (sign2 > 0.0) then
        ! We have gone from negative to positive, indicating the start
        ! of a new positive group.
        npairs = npairs + 1
        index_uw(i,npairs*2) = j
      else
        ! We have gone from positive to negative, indicating the end of
        ! the group.
        index_uw(i,npairs*2+1) = j - 1
      end if
    end if

    16 continue  ! Exclude points outside of the domain

  end do

  if (uwind(i,Ny) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 > 0.0) then
    index_uw(i,npairs*2+1) = Ny
  end if
  index_uw(i,1) = npairs
end do

! Sweep looking for EASTERLY winds, looking for negative numbers.
do i = 1, Nx

  if (uwind(i,1) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 < 0.0) then
    npairs = 1
    index_ue(i,2) = 1
  else
    npairs = 0
  end if

  do j = 2, Ny

    if (isnan(dem(i,j))) goto 17 ! Exclude points outside of the domain

    if (uwind(i,j-1) <= 0.0) then
      sign1 = -1.0
    else
      sign1 = 1.0
    end if
    if (uwind(i,j) <= 0.0) then
      sign2 = -1.0
    else
      sign2 = 1.0
    end if

    if (sign2 /= sign1) then
      ! We have a sign change.
      if (sign2 < 0.0) then
        ! We have gone from positive to negative, indicating the start
        ! of a new negative group.
        npairs = npairs + 1
        index_ue(i,npairs*2) = j
      else
        ! We have gone from negative to positive, indicating the end of
        ! the group.
        index_ue(i,npairs*2+1) = j - 1
      end if
    end if

    17 continue  ! Exclude points outside of the domain

  end do

  if (uwind(i,Ny) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 < 0.0) then
    index_ue(i,npairs*2+1) = Ny
  end if
  index_ue(i,1) = npairs
end do

! Sweep looking for SOUTHERLY winds, looking for positive numbers.
do j = 1, Ny

  if (vwind(1,j) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 > 0.0) then
    npairs = 1
    index_vs(j,2) = 1
  else
    npairs = 0
  end if

  do i = 2, Nx

    if (isnan(dem(i,j))) goto 18 ! Exclude points outside of the domain

    if (vwind(i-1,j) <= 0.0) then
      sign1 = -1.0
    else
      sign1 = 1.0
    end if
    if (vwind(i,j) <= 0.0) then
      sign2 = -1.0
    else
      sign2 = 1.0
    end if

    if (sign2 /= sign1) then
      ! We have a sign change.
      if (sign2 > 0.0) then
        ! We have gone from negative to positive, indicating the start
        ! of a new positive group.
        npairs = npairs + 1
        index_vs(j,npairs*2) = i
      else
        ! We have gone from positive to negative, indicating the end of
        ! the group.
        index_vs(j,npairs*2+1) = i - 1
      end if
    end if

    18 continue  ! Exclude points outside of the domain

  end do

  if (vwind(Nx,j) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 > 0.0) then
    index_vs(j,npairs*2+1) = Nx
  end if
  index_vs(j,1) = npairs
end do

! Sweep looking for NORTHERLY winds, looking for negative numbers.
do j = 1, Ny

  if (vwind(1,j) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 < 0.0) then
    npairs = 1
    index_vn(j,2) = 1
  else
    npairs = 0
  end if

  do i = 2, Nx

    if (isnan(dem(i,j))) goto 19 ! Exclude points outside of the domain

    if (vwind(i-1,j) <= 0.0) then
      sign1 = -1.0
    else
      sign1 = 1.0
    end if
    if (vwind(i,j) <= 0.0) then
      sign2 = -1.0
    else
      sign2 = 1.0
    end if

    if (sign2 /= sign1) then
      ! We have a sign change.
      if (sign2 < 0.0) then
        ! We have gone from positive to negative, indicating the start
        ! of a new negative group.
        npairs = npairs + 1
        index_vn(j,npairs*2) = i
      else
        ! We have gone from negative to positive, indicating the end of
        ! the group.
        index_vn(j,npairs*2+1) = i - 1
      end if
    end if

    19 continue  ! Exclude points outside of the domain

  end do

  if (vwind(Nx,j) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 < 0.0) then
    index_vn(j,npairs*2+1) = Nx
  end if
  index_vn(j,1) = npairs
end do

end subroutine getdirection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getnewdepth(delta_WE,delta_SN,Qs_u,Qs_v, &
                       dm_s, &
                       index_ue,index_uw,index_vn,index_vs, &
                       snowthickness,Ds_soft,snowdepth0,Sice0)

use DRIVING, only: &
  dt                     ! Timestep (s)

use GRID, only: &
  Nx,Ny                  ! Grid dimensions

use PARAMMAPS, only: &
  vegsnowd_xy            ! Vegetation snow holding capacity (m)

use STATE_VARIABLES, only: &
  Ds,                   &! Snow layer thicknesses (m)
  fsnow,                &! Snow cover fraction 
  Sice,                 &! Ice content of snow layers (kg/m^2)
  Sliq                   ! Liquid content of snow layers (kg/m^2)

use LANDUSE, only: &
  dem                    ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  rho_snow               ! Constant snow density (kg/m^3)

implicit none

real, intent(in) :: &
  delta_WE,             &! Grid cell size in WE direction (m)
  delta_SN               ! Grid cell size in SN direction (m)

integer, intent(in) :: &
  index_ue(Nx,2*Ny+1),  &! Wind index array E
  index_uw(Nx,2*Ny+1),  &! Wind index array W
  index_vn(Ny,2*Nx+1),  &! Wind index array N
  index_vs(Ny,2*Nx+1)    ! Wind index array S

real, intent(inout) :: &
  snowthickness(Nx,Ny), &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  Qs_u(Nx,Ny),          &! x component of either Qsalt or Qsusp (kg/m/s)
  Qs_v(Nx,Ny),          &! y component of either Qsalt or Qsusp (kg/m/s)
  dm_s(Nx,Ny),          &! SWE change due to either saltation or suspension (kg/m^2)
  Ds_soft(Nx,Ny),       &! Soft snow thickness (m)
  snowdepth0(Nx,Ny),    &! Snow depth of snowdrift accumulation, averaged over the grid cell (m)
  Sice0(Nx,Ny)           ! Surface layer ice content for transported snow (kg/m^2)

real :: &
  dh_s_u(Nx,Ny),        &! x component of snow depth change due to either saltation or suspension (m)
  dh_s_v(Nx,Ny),        &! y component of snow depth change due to either saltation or suspension (m)
  dm_s_u(Nx,Ny),        &! x component of SWE change due to either saltation or suspension (kg/m^2)
  dm_s_v(Nx,Ny),        &! y component of SWE change due to either saltation or suspension (kg/m^2)
  dm_s_u_loss(Nx,Ny),   &! x component of SWE loss due to either saltation or suspension (kg/m^2)
  dm_s_v_loss(Nx,Ny),   &! y component of SWE loss due to either saltation or suspension (kg/m^2)
  dh_s_u_loss(Nx,Ny),   &! x component of snow depth loss due to either saltation or suspension (m)
  dh_s_v_loss(Nx,Ny),   &! y component of snow depth loss due to either saltation or suspension (m)
  dm_s_u_gain(Nx,Ny),   &! x component of SWE gain due to either saltation or suspension (kg/m^2)
  dm_s_v_gain(Nx,Ny)     ! y component of SWE gain due to either saltation or suspension (kg/m^2)

real :: &
  dh_s,                 &! Snow depth change due to either saltation or suspension (m)
  dh_s_u_gain,          &! x component of snow depth gain due to either saltation or suspension (m)
  dh_s_v_gain,          &! y component of snow depth gain due to either saltation or suspension (m)
  dh_s_gain,            &! Snow depth gain due to either saltation or suspension (m)
  dh_s_loss,            &! Snow depth loss due to either saltation or suspension (m)
  dm_s_loss,            &! SWE loss due to either saltation or suspension (kg/m^2)
  dm_s_gain,            &! SWE gain due to either saltation or suspension (kg/m^2)
  dm_s_u_loss_tmp,      &! Temporary variable for x component of SWE loss due to either saltation or suspension (kg/m^2)
  dm_s_v_loss_tmp,      &! Temporary variable for y component of SWE loss due to either saltation or suspension (kg/m^2)
  dh_s_u_loss_tmp,      &! Temporary variable for x component of snow depth loss due to either saltation or suspension (m)
  dh_s_v_loss_tmp,      &! Temporary variable for y component of snow depth loss due to either saltation or suspension (m)
  Ds_hard,              &! Hard snow thickness (m)
  eps,                  &! Epsilon
  weight_u,             &! x weight for dm_s
  weight_v,             &! y weight for dm_s
  snowdmin,             &! Minimum snow depth to allow transport (m)
  swe_loc                ! SWE in the snow covered part of the pixel (kg/m^2)

integer :: &
  i,j,k,                &! Point counters
  istart,iend,          &! Point counters boundaries
  jstart,jend            ! Point counters boundaries  

! Define an upwind boundary condition for saltation (here I have
! assumed that the transport is in equilibrium).
do i = 1, Nx
  do j = 1, Ny

    if (isnan(dem(i,j))) goto 20 ! Exclude points outside of the domain

    dh_s_u(i,j) = 0.0
    dh_s_v(i,j) = 0.0
    dm_s_u(i,j) = 0.0
    dm_s_v(i,j) = 0.0
    dm_s_u_gain(i,j) = 0.0
    dm_s_v_gain(i,j) = 0.0
    dm_s_u_loss(i,j) = 0.0
    dm_s_v_loss(i,j) = 0.0
    dh_s_u_loss(i,j) = 0.0
    dh_s_v_loss(i,j) = 0.0

    20 continue  ! Exclude points outside of the domain

  end do
end do

! Consider WESTERLY winds.
do i = 1, Nx
  do k = 1, index_uw(i,1)
    jstart = index_uw(i,k*2)+1
    jend = index_uw(i,k*2+1)
    do j = jstart, jend
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))
      dm_s_u_loss(i,j) = dt * Qs_u(i,j) * fsnow(i,j) / delta_WE
      dm_s_u_gain(i,j) = dt * Qs_u(i,j-1) * fsnow(i,j-1) / delta_WE
      dm_s_u_loss(i,j) = min(dm_s_u_loss(i,j),swe_loc) ! No need to adjust Qs_u here because if thresholded, it will be done in the next if loop anyway
      if (dm_s_u_loss(i,j) > epsilon(dm_s_u_loss)) then
        dm_s_u_loss_tmp = dm_s_u_loss(i,j)
        call HS_FROM_SWE(dm_s_u_loss_tmp,dh_s_u_loss_tmp,i,j)
        dh_s_u_loss(i,j) = dh_s_u_loss_tmp
      else
        dm_s_u_loss(i,j) = 0.0
        dh_s_u_loss(i,j) = 0.0
      end if
      dh_s_u_gain = dm_s_u_gain(i,j) / rho_snow
      dh_s_u(i,j) = dh_s_u_gain - dh_s_u_loss(i,j)
      dm_s_u(i,j) = dm_s_u_gain(i,j) - dm_s_u_loss(i,j)

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(vegsnowd_xy(i,j),Ds_hard)
      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dh_s_u_loss(i,j) / fsnow(i,j) <= snowdmin) then
          dh_s_u_loss(i,j) = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dh_s_u_loss(i,j) > epsilon(dh_s_u_loss)) then
            dh_s_u_loss_tmp = dh_s_u_loss(i,j)
            call SWE_FROM_HS(dh_s_u_loss_tmp,dm_s_u_loss_tmp,i,j)
            dm_s_u_loss(i,j) = dm_s_u_loss_tmp
          else
            dh_s_u_loss(i,j) = 0.0
            dm_s_u_loss(i,j) = 0.0
          end if
          dh_s_u(i,j) = dh_s_u_gain - dh_s_u_loss(i,j)
          dm_s_u(i,j) = dm_s_u_gain(i,j) - dm_s_u_loss(i,j)
          Qs_u(i,j) = Qs_u(i,j-1) - dm_s_u(i,j) * delta_WE / dt / fsnow(i,j)
        end if
      else
        Qs_u(i,j) = 0.0
        dh_s_u(i,j) = dh_s_u_gain
        dm_s_u_loss(i,j) = 0.0
        dh_s_u_loss(i,j) = 0.0
        dm_s_u(i,j) = dm_s_u_gain(i,j)
      end if
    end do
  end do
end do

! Consider EASTERLY winds.
do i = 1, Nx
  do k = 1, index_ue(i,1)
    jend = index_ue(i,k*2)
    jstart = index_ue(i,k*2+1)-1
    do j = jstart, jend,-1
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))
      dm_s_u_loss(i,j) = dt * Qs_u(i,j) * fsnow(i,j) / delta_WE
      dm_s_u_gain(i,j) = dt * Qs_u(i,j+1) * fsnow(i,j+1) / delta_WE
      dm_s_u_loss(i,j) = min(dm_s_u_loss(i,j),swe_loc) ! No need to adjust Qs_u here because if thresholded, it will be done in the next if loop anyway
      if (dm_s_u_loss(i,j) > epsilon(dm_s_u_loss)) then
        dm_s_u_loss_tmp = dm_s_u_loss(i,j)
        call HS_FROM_SWE(dm_s_u_loss_tmp,dh_s_u_loss_tmp,i,j)
        dh_s_u_loss(i,j) = dh_s_u_loss_tmp
      else
        dm_s_u_loss(i,j) = 0.0
        dh_s_u_loss(i,j) = 0.0
      end if
      dh_s_u_gain = dm_s_u_gain(i,j) / rho_snow
      dh_s_u(i,j) = dh_s_u_gain - dh_s_u_loss(i,j)
      dm_s_u(i,j) = dm_s_u_gain(i,j) - dm_s_u_loss(i,j)

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(vegsnowd_xy(i,j),Ds_hard)
      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dh_s_u_loss(i,j) / fsnow(i,j) <= snowdmin) then
          dh_s_u_loss(i,j) = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dh_s_u_loss(i,j) > epsilon(dh_s_u_loss)) then
            dh_s_u_loss_tmp = dh_s_u_loss(i,j)
            call SWE_FROM_HS(dh_s_u_loss_tmp,dm_s_u_loss_tmp,i,j)
            dm_s_u_loss(i,j) = dm_s_u_loss_tmp
          else
            dh_s_u_loss(i,j) = 0.0
            dm_s_u_loss(i,j) = 0.0
          end if
          dh_s_u(i,j) = dh_s_u_gain - dh_s_u_loss(i,j)
          dm_s_u(i,j) = dm_s_u_gain(i,j) - dm_s_u_loss(i,j)
          Qs_u(i,j) = Qs_u(i,j+1) - dm_s_u(i,j) * delta_WE / dt / fsnow(i,j)
        end if
      else
        Qs_u(i,j) = 0.0
        dh_s_u(i,j) = dh_s_u_gain
        dm_s_u_loss(i,j) = 0.0
        dh_s_u_loss(i,j) = 0.0
        dm_s_u(i,j) = dm_s_u_gain(i,j)
      end if
    end do
  end do
end do

! Consider SOUTHERLY winds.
do j = 1, Ny
  do k = 1, index_vs(j,1)
    istart = index_vs(j,k*2)+1
    iend = index_vs(j,k*2+1)
    do i = istart, iend
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))
      dm_s_v_loss(i,j) = dt * Qs_v(i,j) * fsnow(i,j) / delta_SN
      dm_s_v_gain(i,j) = dt * Qs_v(i-1,j) * fsnow(i-1,j) / delta_SN
      dm_s_v_loss(i,j) = min(dm_s_v_loss(i,j),swe_loc) ! No need to adjust Qs_v here because if thresholded, it will be done in the next if loop anyway
      if (dm_s_v_loss(i,j) > epsilon(dm_s_v_loss)) then
        dm_s_v_loss_tmp = dm_s_v_loss(i,j)
        call HS_FROM_SWE(dm_s_v_loss_tmp,dh_s_v_loss_tmp,i,j)
        dh_s_v_loss(i,j) = dh_s_v_loss_tmp
      else
        dm_s_v_loss(i,j) = 0.0
        dh_s_v_loss(i,j) = 0.0
      end if
      dh_s_v_gain = dm_s_v_gain(i,j) / rho_snow
      dh_s_v(i,j) = dh_s_v_gain - dh_s_v_loss(i,j)
      dm_s_v(i,j) = dm_s_v_gain(i,j) - dm_s_v_loss(i,j)

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(vegsnowd_xy(i,j),Ds_hard)
      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dh_s_v_loss(i,j) / fsnow(i,j) <= snowdmin) then
          dh_s_v_loss(i,j) = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dh_s_v_loss(i,j) > epsilon(dh_s_v_loss)) then
            dh_s_v_loss_tmp = dh_s_v_loss(i,j)
            call SWE_FROM_HS(dh_s_v_loss_tmp,dm_s_v_loss_tmp,i,j)
            dm_s_v_loss(i,j) = dm_s_v_loss_tmp
          else
            dh_s_v_loss(i,j) = 0.0
            dm_s_v_loss(i,j) = 0.0
          end if
          dh_s_v(i,j) = dh_s_v_gain - dh_s_v_loss(i,j)
          dm_s_v(i,j) = dm_s_v_gain(i,j) - dm_s_v_loss(i,j)
          Qs_v(i,j) = Qs_v(i-1,j) - dm_s_v(i,j) * delta_SN / dt / fsnow(i,j)
        end if
      else
        Qs_v(i,j) = 0.0
        dh_s_v(i,j) = dh_s_v_gain
        dm_s_v_loss(i,j) = 0.0
        dh_s_v_loss(i,j) = 0.0
        dm_s_v(i,j) = dm_s_v_gain(i,j)
      end if
    end do
  end do
end do

! Consider NORTHERLY winds.
do j = 1, Ny
  do k = 1, index_vn(j,1)
    iend = index_vn(j,k*2)
    istart = index_vn(j,k*2+1)-1
    do i = istart, iend,-1
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))
      dm_s_v_loss(i,j) = dt * Qs_v(i,j) * fsnow(i,j) / delta_SN
      dm_s_v_gain(i,j) = dt * Qs_v(i+1,j) * fsnow(i+1,j) / delta_SN
      dm_s_v_loss(i,j) = min(dm_s_v_loss(i,j),swe_loc) ! No need to adjust Qs_v here because if thresholded, it will be done in the next if loop anyway
      if (dm_s_v_loss(i,j) > epsilon(dm_s_v_loss)) then
        dm_s_v_loss_tmp = dm_s_v_loss(i,j)
        call HS_FROM_SWE(dm_s_v_loss_tmp,dh_s_v_loss_tmp,i,j)
        dh_s_v_loss(i,j) = dh_s_v_loss_tmp
      else
        dm_s_v_loss(i,j) = 0.0
        dh_s_v_loss(i,j) = 0.0
      end if
      dh_s_v_gain = dm_s_v_gain(i,j) / rho_snow
      dh_s_v(i,j) = dh_s_v_gain - dh_s_v_loss(i,j)
      dm_s_v(i,j) = dm_s_v_gain(i,j) - dm_s_v_loss(i,j)

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(vegsnowd_xy(i,j),Ds_hard)
      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dh_s_v_loss(i,j) / fsnow(i,j) <= snowdmin) then
          dh_s_v_loss(i,j) = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dh_s_v_loss(i,j) > epsilon(dh_s_v_loss)) then
            dh_s_v_loss_tmp = dh_s_v_loss(i,j)
            call SWE_FROM_HS(dh_s_v_loss_tmp,dm_s_v_loss_tmp,i,j)
            dm_s_v_loss(i,j) = dm_s_v_loss_tmp
          else
            dh_s_v_loss(i,j) = 0.0
            dm_s_v_loss(i,j) = 0.0
          end if
          dh_s_v(i,j) = dh_s_v_gain - dh_s_v_loss(i,j)
          dm_s_v(i,j) = dm_s_v_gain(i,j) - dm_s_v_loss(i,j)
          Qs_v(i,j) = Qs_v(i+1,j) - dm_s_v(i,j) * delta_SN / dt / fsnow(i,j)
        end if
      else
        Qs_v(i,j) = 0.0
        dh_s_v(i,j) = dh_s_v_gain
        dm_s_v_loss(i,j) = 0.0
        dh_s_v_loss(i,j) = 0.0
        dm_s_v(i,j) = dm_s_v_gain(i,j)
      end if
    end do
  end do
end do

! Update the snow depth changes due to saltation transport from the
! the east and west, and north and south.
eps = 1e-6
do i = 1, Nx
  do j = 1, Ny

    if (isnan(dem(i,j))) goto 21 ! Exclude points outside of the domain

    ! LQ: I can't understand the meaning of this weighting, I keep it as original SNOWTRAN3D though
    weight_u = abs(dh_s_u(i,j)) / &
               (abs(dh_s_u(i,j)) + abs(dh_s_v(i,j)) + eps)

    weight_v = abs(dh_s_v(i,j)) / &
               (abs(dh_s_u(i,j)) + abs(dh_s_v(i,j)) + eps)

    dm_s_u(i,j) = weight_u * dm_s_u(i,j)
    dm_s_v(i,j) = weight_v * dm_s_v(i,j)
    dm_s_u_gain(i,j) = weight_u * dm_s_u_gain(i,j)
    dm_s_v_gain(i,j) = weight_v * dm_s_v_gain(i,j)
    dm_s_u_loss(i,j) = weight_u * dm_s_u_loss(i,j)
    dm_s_v_loss(i,j) = weight_v * dm_s_v_loss(i,j)

    dm_s(i,j) = dm_s_u(i,j) + dm_s_v(i,j)
    dm_s_gain = dm_s_u_gain(i,j) + dm_s_v_gain(i,j)
    dm_s_loss = dm_s_u_loss(i,j) + dm_s_v_loss(i,j)
    dh_s_gain = dm_s_gain / rho_snow
    if (dm_s_loss > epsilon(dm_s_loss)) then
      call HS_FROM_SWE(dm_s_loss,dh_s_loss,i,j)
    else
      dm_s_loss = 0.0
      dh_s_loss = 0.0
    end if
    dh_s = dh_s_gain - dh_s_loss

    Ds_soft(i,j) = max(Ds_soft(i,j) - dh_s_loss,0.0)

    if (dm_s_loss > epsilon(dm_s_loss) .and. dh_s_loss > epsilon(dh_s_loss)) then
      call SNOW_ABLATION(dh_s_loss,dm_s_loss,i,j)
    end if

    ! Net mass gain for this grid cell at this time step.
    if (dm_s_gain > epsilon(dm_s_gain) .and. dh_s_gain > epsilon(dh_s_gain)) then

      ! Add to the existing top layer.
      Sice0(i,j) = Sice0(i,j) + dm_s_gain
      snowdepth0(i,j) = snowdepth0(i,j) + dh_s_gain

    end if

    ! Update the snow layer thicknesses
    snowthickness(i,j) = sum(Ds(:,i,j))

    21 continue  ! Exclude points outside of the domain

  end do
end do

end subroutine getnewdepth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine subgrid_1(snowdepth_tabler, &
                     index_ue,index_uw,index_vn,index_vs, &
                     tabler_nn,tabler_ss,tabler_ee,tabler_ww, &
                     tabler_ne,tabler_se,tabler_sw,tabler_nw,uwind,vwind)

! This subroutine forces SnowTran-3D's snow accumluation profiles
! to be bounded by the equilibrium topographic drift catchment
! profiles observed and modeled by Tabler (1975).

! Tabler, R. D., 1975: Predicting profiles of snowdrifts in
! topographic catchments.  Proceedings of the 43rd Annual Western
! Snow Conference, San Diego, California, 87-97.

use GRID, only: &
  Nx,Ny                   ! Grid dimensions

use LANDUSE, only: &
  dem                     ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  tabler_dir              ! Dominant wind direction

implicit none

real, intent(in) :: &
  uwind(Nx,Ny),          &! x component of wind speed (m/s)
  vwind(Nx,Ny),          &! y component of wind speed (m/s)
  tabler_nn(Nx,Ny),      &! Tabler surfaces NN
  tabler_ss(Nx,Ny),      &! Tabler surfaces SS
  tabler_ee(Nx,Ny),      &! Tabler surfaces EE
  tabler_ww(Nx,Ny),      &! Tabler surfaces WW
  tabler_ne(Nx,Ny),      &! Tabler surfaces NE
  tabler_se(Nx,Ny),      &! Tabler surfaces SE
  tabler_sw(Nx,Ny),      &! Tabler surfaces SW
  tabler_nw(Nx,Ny)        ! Tabler surfaces NW

real, intent(inout) :: &
  snowdepth_tabler(Nx,Ny) ! Snow depth (m)

integer, intent(in) :: &
  index_ue(Nx,2*Ny+1),   &! Wind index array E
  index_uw(Nx,2*Ny+1),   &! Wind index array W
  index_vn(Ny,2*Nx+1),   &! Wind index array N
  index_vs(Ny,2*Nx+1)     ! Wind index array S

integer :: &
  i,j,k,                 &! Point counters
  istart,iend,           &! Point counters boundaries
  jstart,jend             ! Point counters boundaries

real :: &
  tabler,                &! Tabler surface (m)
  snow_sfc,              &! Snow surface (m)
  snow_d_extra            ! Extra snow depth (m)

real :: &
  snow_d1(Nx,Ny),        &! Temporary snow depth variable (m)
  snow_d2(Nx,Ny),        &! Temporary snow depth variable (m)
  weight_u(Nx,Ny),       &! Weight for snow depth on x axis
  weight_v(Nx,Ny)         ! Weight for snow depth on y axis

! Create a copy of the incoming snow depth distribution.  Also define
! the u and v weighting functions.
do j = 1, Ny
  do i = 1, Nx

    if (isnan(dem(i,j))) goto 22 ! Exclude points outside of the domain

    snow_d1(i,j) = snowdepth_tabler(i,j)
    snow_d2(i,j) = snowdepth_tabler(i,j)

    weight_u(i,j) = abs(uwind(i,j)) / &
                    sqrt(uwind(i,j)**2 + vwind(i,j)**2)
    weight_v(i,j) = abs(vwind(i,j)) / &
                    sqrt(uwind(i,j)**2 + vwind(i,j)**2)

    22 continue  ! Exclude points outside of the domain

  end do
end do

! Consider WESTERLY winds.
do i = 1, Nx
  do k = 1, index_uw(i,1)
    jstart = index_uw(i,k*2)+1
    jend = index_uw(i,k*2+1)
    do j = jstart, jend

      if (tabler_dir > 337.5 .AND. tabler_dir <= 360.0 .OR. &
         tabler_dir >= 0.0 .AND. tabler_dir <= 22.5) then
        tabler = tabler_nn(i,j)
      else if (tabler_dir > 157.5 .AND. tabler_dir <= 202.5) then
        tabler = tabler_ss(i,j)
      else if (tabler_dir > 202.5 .AND. tabler_dir <= 247.5) then
        tabler = tabler_sw(i,j)
      else if (tabler_dir > 247.5 .AND. tabler_dir <= 292.5) then
        tabler = tabler_ww(i,j)
      else if (tabler_dir > 292.5 .AND. tabler_dir <= 337.5) then
        tabler = tabler_nw(i,j)
      end if

      snow_sfc = tabler

      if (snow_d1(i,j) > snow_sfc) then
        snow_d_extra = (snow_d1(i,j) - snow_sfc) * weight_u(i,j)
        snow_d1(i,j) = snow_d1(i,j) - snow_d_extra
        if (i < Nx) then
          snow_d1(i,j+1) = snow_d1(i,j+1) + snow_d_extra
        else
          snow_d1(i,j) = snow_d1(i,j)
        end if
      end if

    end do
  end do
end do

! Consider EASTERLY winds.
do i = 1, Nx
  do k = 1, index_ue(i,1)
    jend = index_ue(i,k*2)
    jstart = index_ue(i,k*2+1)-1
    do j = jstart, jend,-1

      if (tabler_dir > 337.5 .AND. tabler_dir <= 360.0 .OR. &
          tabler_dir >= 0.0 .AND. tabler_dir <= 22.5) then
        tabler = tabler_nn(i,j)
      else if (tabler_dir > 22.5 .AND. tabler_dir <= 67.5) then
        tabler = tabler_ne(i,j)
      else if (tabler_dir > 67.5 .AND. tabler_dir <= 112.5) then
        tabler = tabler_ee(i,j)
      else if (tabler_dir > 112.5 .AND. tabler_dir <= 157.5) then
        tabler = tabler_se(i,j)
      else if (tabler_dir > 157.5 .AND. tabler_dir <= 202.5) then
        tabler = tabler_ss(i,j)
      end if

      snow_sfc = tabler

      if (snow_d1(i,j) > snow_sfc) then
        snow_d_extra = (snow_d1(i,j) - snow_sfc) * weight_u(i,j)
        snow_d1(i,j) = snow_d1(i,j) - snow_d_extra
        if (i > 1) then
          snow_d1(i,j-1) = snow_d1(i,j-1) + snow_d_extra
        else
          snow_d1(i,j) = snow_d1(i,j)
        end if
      end if
    end do
  end do
end do

! Consider SOUTHERLY winds.
do j = 1, Ny
  do k = 1, index_vs(j,1)
    istart = index_vs(j,k*2)+1
    iend = index_vs(j,k*2+1)
    do i = istart, iend

      if (tabler_dir > 67.5 .AND. tabler_dir <= 112.5) then
        tabler = tabler_ee(i,j)
      else if (tabler_dir > 112.5 .AND. tabler_dir <= 157.5) then
        tabler = tabler_se(i,j)
      else if (tabler_dir > 157.5 .AND. tabler_dir <= 202.5) then
        tabler = tabler_ss(i,j)
      else if (tabler_dir > 202.5 .AND. tabler_dir <= 247.5) then
        tabler = tabler_sw(i,j)
      else if (tabler_dir > 247.5 .AND. tabler_dir <= 292.5) then
        tabler = tabler_ww(i,j)
      end if

      snow_sfc = tabler

      if (snow_d2(i,j) > snow_sfc) then
        snow_d_extra = (snow_d2(i,j) - snow_sfc) * weight_v(i,j)
        snow_d2(i,j) = snow_d2(i,j) - snow_d_extra
        if (j < Ny) then
          snow_d2(i+1,j) = snow_d2(i+1,j) + snow_d_extra
        else
          snow_d2(i,j) = snow_d2(i,j)
        end if
      end if
    end do
  end do
end do

! Consider NORTHERLY winds.
do j = 1, Ny
  do k = 1, index_vn(j,1)
    iend = index_vn(j,k*2)
    istart = index_vn(j,k*2+1)-1
    do i = istart, iend,-1

      if (tabler_dir > 337.5 .AND. tabler_dir <= 360.0 .OR. &
          tabler_dir >= 0.0 .AND. tabler_dir <= 22.5) then
        tabler = tabler_nn(i,j)
      else if (tabler_dir > 22.5 .AND. tabler_dir <= 67.5) then
        tabler = tabler_ne(i,j)
      else if (tabler_dir > 67.5 .AND. tabler_dir <= 112.5) then
        tabler = tabler_ee(i,j)
      else if (tabler_dir > 247.5 .AND. tabler_dir <= 292.5) then
        tabler = tabler_ww(i,j)
      else if (tabler_dir > 292.5 .AND. tabler_dir <= 337.5) then
        tabler = tabler_nw(i,j)
      end if

      snow_sfc = tabler

      if (snow_d2(i,j) > snow_sfc) then
        snow_d_extra = (snow_d2(i,j) - snow_sfc) * weight_v(i,j)
        snow_d2(i,j) = snow_d2(i,j) - snow_d_extra
        if (j > 1) then
          snow_d2(i-1,j) = snow_d2(i-1,j) + snow_d_extra
        else
          snow_d2(i,j) = snow_d2(i,j)
        end if
      end if

    end do
  end do
end do

! Update the snow depths resulting from these redistributions.
do j = 1, Ny
  do i = 1, Nx

    if (isnan(dem(i,j))) goto 23 ! Exclude points outside of the domain

    snowdepth_tabler(i,j) = snow_d1(i,j) * weight_u(i,j) + &
                  snow_d2(i,j) * weight_v(i,j)

    23 continue  ! Exclude points outside of the domain

  end do
end do

! Clean up the boundaries.  Make the boundary values equal to
! the values just inside the boundaries.
do i = 2, Nx-1
  snowdepth_tabler(i,1) = snowdepth_tabler(i,2)
  snowdepth_tabler(i,Ny) = snowdepth_tabler(i,Ny-1)
end do
do j = 1, Ny
  snowdepth_tabler(1,j) = snowdepth_tabler(2,j)
  snowdepth_tabler(Nx,j) = snowdepth_tabler(Nx-1,j)
end do

end subroutine subgrid_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine surface_snow(Utau_t)

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

use DRIVING, only: &
  Ua,                &! Wind speed (m/s)
  dt,                &! Timestep (s)
  zU                  ! Wind speed measurement height (m)

use PARAMETERS, only : &
  rhos_min,          &! Minimum snow density (kg/m^3)
  rhos_max,          &! Maximum snow density (kg/m^3)
  z0sn                ! Snow roughness length (m)

use PARAM_SNOWTRAN3D, only: &
  rho_snow            ! Constant snow density (kg/m^3)

use STATE_VARIABLES, only: &
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Ds,                &! Snow layer thicknesses (m)
  Tsnow               ! Snow layer temperatures (K)

implicit none

real, intent(inout) :: &
  Utau_t(Nx,Ny)       ! Threshold friction velocity (m/s)

integer :: &
  i,j                 ! Point counters

real :: &
  A1,                &! Coefficient A1 in Liston et al. (2007), eq. 17 p. 245 (m^-1)
  A2,                &! Coefficient A2 in Liston et al. (2007), eq. 17 p. 245 (m^3/kg)
  B,                 &! Coefficient B in Liston et al. (2007), eq. 17 p. 245 (K-1)
  U,                 &! Wind speed contribution to snow compaction
  rho_surf_snow,     &! Density of top snow layer  (kg/m^3)
  C,                 &! Density rate coefficient, C in Liston et al. (2007), eq. 17 p. 245
  alpha,             &! Coefficient E3 in Liston et al. (2007), eq. 18 p. 245 (m/s)
  windspd_2m          ! Windspeed at 2m height (m/s)

! Define the density rate coefficients.
C = 0.10

! Define alpha. 
alpha = 0.2

do j = 1, Ny
  do i = 1, Nx

    if (isnan(dem(i,j))) goto 24 ! Exclude points outside of the domain

    ! Calculate the 2-m wind speed.
    windspd_2m = Ua(i,j) * log(2.0/z0sn)/log(zU/z0sn)

    ! Initialize coefficients
    A1 = 0.0013
    A2 = 0.021
    B = 0.08

    ! Evolve the near-surface snow density under the influence of
    ! temperature and snow-transporting wind speeds.

    ! Update the snow density of the soft snow layer.  Eliminate the
    ! wind speed influence for speeds below 5 m/s, but account for it
    ! if speeds are >= 5 m/s.
    if (windspd_2m >= 5.0) then
      U = 5.0 + 15.0 * (1.0 - exp(-(alpha*(windspd_2m - 5.0))))
    else
      U = 1.0
    end if

    rho_surf_snow = rho_snow
    if (Nsnow(i,j) > 0 .and. Ds(1,i,j) > epsilon(Ds)) then
      rho_surf_snow = (Sice(1,i,j) + Sliq(1,i,j)) / Ds(1,i,j) / fsnow(i,j)
      rho_surf_snow = rho_surf_snow + dt * &
                      (C * A1 * U * rho_surf_snow * &
                      exp((- B)*(Tm-Tsnow(1,i,j))) * exp((- A2)*rho_surf_snow))
      ! Bound the calculated density.
      rho_surf_snow = min(rhos_max,rho_surf_snow)
      rho_surf_snow = max(rhos_min,rho_surf_snow)
      ! Update surface snow layer thickness after wind compaction
      Ds(1,i,j) = (Sice(1,i,j) + Sliq(1,i,j)) / rho_surf_snow / fsnow(i,j)
    end if

    ! Calculate the snow threshold friction velocity.
    if (rho_surf_snow <= 300.0) then
      Utau_t(i,j) = 0.10 * exp(0.003 * rho_surf_snow)
    else
      Utau_t(i,j) = 0.005 * exp(0.013 * rho_surf_snow)
    end if

    24 continue  ! Exclude points outside of the domain

  end do
end do

end subroutine surface_snow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!