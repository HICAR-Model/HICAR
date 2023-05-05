!-----------------------------------------------------------------------
! This routine performs three tasks:
! 1) Addition of new snow accumulation (snowfall or redistribution)
!    to the layers
! 2) Calculation of snow cover fraction and rescaling of layer 
!    thicknesses with the new fsnow
! 3) Relayering of the snowpack with the chosen method
!-----------------------------------------------------------------------
subroutine SNOW_LAYERING(snowdepth0,Sice0)

use MODCONF, only: CANMOD,HN_ON,SNOLAY

use CONSTANTS, only: &
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcap_wat,          &! Specific heat capacity of water (J/K/kg)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  Ta                  ! Air temperature (K)
  
use GRID, only: &
  Nsmax,             &! Maximum number of snow layers
  Nx,Ny,             &! Grid dimensions
  Dzsnow,            &! Minimum snow layer thicknesses (m)
  Ds_min,            &! Minimum possible snow layer thickness (m)
  Ds_surflay          ! Maximum thickness of surface fine snow layering (m)
  
use PARAMETERS, only: &
  rho0,              &! Fixed snow density (kg/m^3)
  rgr0,              &! Fresh snow grain radius (m)
  fthresh             ! Forest fraction required for forest tile to be considered

use STATE_VARIABLES, only: &
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Ds,                &! Snow layer thicknesses (m)
  histowet,          &! Historical variable for past wetting of a layer (0-1)
  rgrn,              &! Snow layer grain radius (m)
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Tsnow               ! Snow layer temperatures (K)

use LANDUSE, only: &
  dem,               &! Terrain elevation (m)
  forest              ! Grid cell forest fraction

implicit none

real, intent(in) :: &
  snowdepth0(Nx,Ny), &! Snow depth of new accumulation (m)
  Sice0(Nx,Ny)        ! Ice content of new snow accumulation(kg/m^2)

real :: &
  Ds0(Nx,Ny)          ! Layer thickness of new snow accumulation, i.e. on the snow-covered part (m)

! For original layering
real :: &
  csnow(Nsmax),      &! Areal heat capacity of snow (J/K/m^2)
  D(Nsmax),          &! Layer thickness before adjustment (m)
  E(Nsmax),          &! Energy contents before adjustment (J/m^2)
  hw(Nsmax),         &! histowet before adjustment (0-1)
  R(Nsmax),          &! Snow grain radii before adjustment (kg/m^2)
  S(Nsmax),          &! Ice contents before adjustment (kg/m^2)
  U(Nsmax),          &! Layer internal energy contents (J/m^2)
  W(Nsmax)            ! Liquid contents before adjustment (kg/m^2)

! For density-dependent layering
real :: &
  rho(Nsmax+1),          &! Layer density (kg/m^3)
  diff_rho(Nsmax),       &! Absolute difference of density with next layer (kg/m^3)
  csnow_loc(Nsmax+1),    &! Areal heat capacity of snow (J/K/m^2)
  Sice_loc(Nsmax+1),     &! Ice content of snow layers (kg/m^2)
  Sliq_loc(Nsmax+1),     &! Liquid content of snow layers (kg/m^2)
  Ds_loc(Nsmax+1),       &! Snow layer thicknesses (m)
  histowet_loc(Nsmax+1), &! Historical variable for past wetting of a layer (0-1)
  rgrn_loc(Nsmax+1),     &! Snow layer grain radius (m)
  U_loc(Nsmax+1),        &! Layer internal energy contents (J/m^2)
  Tsnow_loc(Nsmax+1)      ! Snow layer temperatures (K)
  
real :: &
  fold,              &! Previous snow cover fraction
  snowdepth,         &! Snow depth, i.e. scaled by fsnow (m)
  snowthickness,     &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  SWEtmp,            &! Temporary snow water equivalent for snow covered fraction calculation (kg/m^2)
  Ds_excess,         &! Temporary variable for snow thickness to transfer to lowest layer (m)
  Ds_old,            &! Temporary variable to save previous Ds (m)
  Sice_old,          &! Temporary variable to save previous Sice (kg/m^2)
  Sliq_old,          &! Temporary variable to save previous Sliq (kg/m^2)
  U_old,             &! Temporary variable to save previous U (J/m^2)
  rho_kmin,          &! Density of layer kmin (kg/m^3)
  rho_kup,           &! Density of layer kup (kg/m^3)
  rho_kdown,         &! Density of layer kdown (kg/m^3)
  Dtemp_surflay,     &! Temporary variable cumulating surface layer thicknesses (m)
  Tsnow0,            &! Falling snow temperature (K)
  dnew,              &! New snow layer thickness (m)
  wt                  ! Layer weighting

integer :: &
  i,j,               &! Point counters
  k,n,               &! Layer counters
  kmin,              &! Layer index with lowest thickness
  kmax,              &! Layer index with highest thickness
  kup,               &! Index of layer above kmin
  kdown,             &! Index of layer under kmin
  kmerge,            &! Index of layer to merge with the next one (minimum rho difference)
  k_surflay,         &! Layer index where Ds_surflay is exceeded
  knew,              &! New snow layer pointer
  kold,              &! Old snow layer pointer
  Nold,              &! Previous number of snow layers
  Nsnow_loc           ! Number of snow layers

! Initialize Ds0
Ds0(:,:) = 0.0

do j = 1, Ny
  do i = 1, Nx

    if (isnan(dem(i,j))) goto 1 ! Exclude points outside of the domain

    if (CANMOD == 1 .and. forest(i,j) < fthresh) goto 1 ! Exclude points that have no forest

    ! Decrease Nsnow if necessary (e.g. after melting)
    do while (Nsnow(i,j) > 0 .and. Ds(1,i,j) < epsilon(Ds))
      if (Nsnow(i,j) > 1) then
        do n = 1 , Nsnow(i,j) - 1
          Ds(n,i,j) = Ds(n+1,i,j)
          Sice(n,i,j) = Sice(n+1,i,j)
          Sliq(n,i,j) = Sliq(n+1,i,j)
          Tsnow(n,i,j) = Tsnow(n+1,i,j)
          rgrn(n,i,j) = rgrn(n+1,i,j)
          histowet(n,i,j) = histowet(n+1,i,j)
        end do
      end if
      Ds(Nsnow(i,j),i,j) = 0
      Sice(Nsnow(i,j),i,j) = 0
      Sliq(Nsnow(i,j),i,j) = 0
      Tsnow(Nsnow(i,j),i,j) = Tm
      rgrn(Nsnow(i,j),i,j) = 0
      histowet(Nsnow(i,j),i,j) = 0.0
      Nsnow(i,j) = Nsnow(i,j) - 1
    end do

    if (SNOLAY == 0) then
      if (Sice(1,i,j) + Sice0(i,j) > epsilon(Sice)) then
        rgrn(1,i,j) = (Sice(1,i,j)*rgrn(1,i,j) + Sice0(i,j)*rgr0) / (Sice(1,i,j) + Sice0(i,j))
      endif
      Sice(1,i,j) = Sice(1,i,j) + Sice0(i,j)
      ! histowet(1,i,j) stays the same as former surface layer
    end if
    snowdepth = sum(Ds(:,i,j)) * fsnow(i,j) + snowdepth0(i,j)

    ! Store previous snow cover fraction
    fold = fsnow(i,j)
    ! Updated Fractional Snow-Covered Area
    SWEtmp = sum(Sice(:,i,j) + Sliq(:,i,j))
    if (SNOLAY == 1) then
      SWEtmp = SWEtmp + Sice0(i,j)
    end if
    call SNOWCOVERFRACTION(snowdepth,SWEtmp,i,j)
    
    ! Rescale Ds with new snow cover fraction
    if (fsnow(i,j) > epsilon(fsnow)) then
      Ds0(i,j) = snowdepth0(i,j) / fsnow(i,j)
      ! Update surface layer thickness based on new fsnow
      if (SNOLAY == 0) then
        Ds(1,i,j) = Ds(1,i,j) * fold / fsnow(i,j) + Ds0(i,j)
      else ! SNOLAY == 1
        Ds(1,i,j) = Ds(1,i,j) * fold / fsnow(i,j)
      end if
    else
      Nsnow(i,j) = 0
      Ds(:,i,j) = 0
      Sice(:,i,j) = 0
      Sliq(:,i,j) = 0
      Tsnow(:,i,j) = Tm
      rgrn(:,i,j) = 0
      histowet(:,i,j) = 0.0
    end if
    if (Nsnow(i,j) > 1) then
      do k = 2, Nsnow(i,j)
        Ds(k,i,j) = Ds(k,i,j) * fold / fsnow(i,j)
      end do
    end if
    
    ! New snow temperature
    Tsnow0 = min(Ta(i,j), Tm)
    if (HN_ON) then
      Tsnow0 = max(Tsnow0, (Tm-40))
    end if
    
    if (SNOLAY == 0) then

      !-----------------------------------------------------------------------
      ! Original layering routine
      !-----------------------------------------------------------------------

      ! New snowpack
      if (Nsnow(i,j) == 0 .and. Sice(1,i,j) > epsilon(Sice(1,i,j))) then
        Nsnow(i,j) = 1
        Tsnow(1,i,j) = Tsnow0
        histowet(1,i,j) = 0
      end if

      ! Store state of old layers
      D(:) = Ds(:,i,j)
      R(:) = rgrn(:,i,j)
      S(:) = Sice(:,i,j)
      W(:) = Sliq(:,i,j)
      hw(:) = histowet(:,i,j)
      if (fsnow(i,j) > epsilon(fsnow)) then
        csnow(1) = (Sice(1,i,j)*hcap_ice + Sliq(1,i,j)*hcap_wat) / fsnow(i,j)
        E(1) = csnow(1)*(Tsnow(1,i,j) - Tm) + (Sice0(i,j) * hcap_ice / fsnow(i,j)) * (Tsnow0 - Tsnow(1,i,j)) ! Adjustment given that csnow(1) already includes the new snow
      else
        E(:) = 0
      end if
      if (Nsnow(i,j) > 1) then
        do k = 2, Nsnow(i,j)
          csnow(k) = (Sice(k,i,j)*hcap_ice + Sliq(k,i,j)*hcap_wat) / fsnow(i,j)
          E(k) = csnow(k)*(Tsnow(k,i,j) - Tm)
        end do
      end if
      Nold = Nsnow(i,j)

      ! Initialise new layers
      Ds(:,i,j) = 0
      rgrn(:,i,j) = 0
      Sice(:,i,j) = 0
      Sliq(:,i,j) = 0
      Tsnow(:,i,j) = Tm
      U(:) = 0
      Nsnow(i,j) = 0
      histowet(:,i,j) = 0

      if (fsnow(i,j) > epsilon(fsnow)) then  ! Existing or new snowpack

        ! Original layering routine
        ! Re-assign and count snow layers
        dnew = snowdepth / fsnow(i,j)
        Ds(1,i,j) = dnew
        k = 1
        if (Ds(1,i,j) > Dzsnow(1)) then 
          do k = 1, Nsmax
            Ds(k,i,j) = Dzsnow(k)
            dnew = dnew - Dzsnow(k)
            if (dnew <= Dzsnow(k) .or. k == Nsmax) then
              Ds(k,i,j) = Ds(k,i,j) + dnew
              exit
            end if
          end do
        end if
        Nsnow(i,j) = k

        ! Fill new layers from the top downwards
        knew = 1
        dnew = Ds(1,i,j)
        do kold = 1, Nold
          do
            if (D(kold) < dnew) then
              ! All snow from old layer partially fills new layer
              rgrn(knew,i,j) = rgrn(knew,i,j) + S(kold)*R(kold)
              Sice(knew,i,j) = Sice(knew,i,j) + S(kold)
              Sliq(knew,i,j) = Sliq(knew,i,j) + W(kold)
              U(knew) = U(knew) + E(kold)
              histowet(knew,i,j) = max(histowet(knew,i,j),hw(kold))
              dnew = dnew - D(kold)
              exit
            else
              ! Some snow from old layer fills new layer
              wt = dnew / D(kold)
              rgrn(knew,i,j) = rgrn(knew,i,j) + wt*S(kold)*R(kold)
              Sice(knew,i,j) = Sice(knew,i,j) + wt*S(kold) 
              Sliq(knew,i,j) = Sliq(knew,i,j) + wt*W(kold)
              U(knew) = U(knew) + wt*E(kold)
              histowet(knew,i,j) = max(histowet(knew,i,j),hw(kold))
              D(kold) = (1 - wt)*D(kold)
              E(kold) = (1 - wt)*E(kold)
              S(kold) = (1 - wt)*S(kold)
              W(kold) = (1 - wt)*W(kold)
              hw(kold) = hw(kold)
              knew = knew + 1
              if (knew > Nsnow(i,j)) exit
              dnew = Ds(knew,i,j)
            end if
          end do
        end do
        ! Diagnose snow layer temperatures
        do k = 1, Nsnow(i,j)
          csnow(k) = (Sice(k,i,j)*hcap_ice + Sliq(k,i,j)*hcap_wat) / fsnow(i,j)
          Tsnow(k,i,j) = Tm + U(k) / csnow(k)
          if (HN_ON) then
            Tsnow(k,i,j) = max(Tsnow(k,i,j), (Tm-40))
          endif
          rgrn(k,i,j) = rgrn(k,i,j) / Sice(k,i,j)
        end do

      end if ! Existing or new snowpack

    else ! SNOLAY == 1

      !-----------------------------------------------------------------------
      ! Density dependent snowpack layering
      ! Note that rgrn is not updated, since unused
      !-----------------------------------------------------------------------

      if (Ds0(i,j) > epsilon(Ds0)) then
        snowthickness = Ds0(i,j) + sum(Ds(:,i,j))
      else
        snowthickness = sum(Ds(:,i,j))
      end if
      rho(:) = rho0

      ! Step 0: Save state variables in local variables that can be up to Nsmax+1
      Sice_loc(1:Nsmax) = Sice(:,i,j)
      Sice_loc(Nsmax+1) = 0
      Sliq_loc(1:Nsmax) = Sliq(:,i,j)
      Sliq_loc(Nsmax+1) = 0
      Ds_loc(1:Nsmax) = Ds(:,i,j)
      Ds_loc(Nsmax+1) = 0
      histowet_loc(1:Nsmax) = histowet(:,i,j)
      histowet_loc(Nsmax+1) = 0
      rgrn_loc(1:Nsmax) = rgrn(:,i,j)
      rgrn_loc(Nsmax+1) = 0
      Tsnow_loc(1:Nsmax) = Tsnow(:,i,j)
      Tsnow_loc(Nsmax+1) = Tm
      Nsnow_loc = Nsnow(i,j)
      if (fsnow(i,j) > epsilon(fsnow)) then
        do k = 1, Nsmax
          csnow_loc(k) = (Sice_loc(k)*hcap_ice + Sliq_loc(k)*hcap_wat) / fsnow(i,j)
          U_loc(k) = csnow_loc(k)*(Tsnow_loc(k) - Tm)
        end do
        U_loc(Nsmax+1) = 0
      else
        U_loc(:) = 0
      end if

      ! Step 1: If there is fresh snow, add the top fresh snow layer and shift layer numbers
      if (Ds0(i,j) > epsilon(Ds0)) then
        Nsnow_loc = Nsnow_loc + 1
        if (Nsnow_loc > 1) then
          do k = 1 , Nsnow_loc-1
            Ds_loc(Nsnow_loc-k+1) = Ds_loc(Nsnow_loc-k)
            Sice_loc(Nsnow_loc-k+1) = Sice_loc(Nsnow_loc-k)
            Sliq_loc(Nsnow_loc-k+1) = Sliq_loc(Nsnow_loc-k)
            U_loc(Nsnow_loc-k+1) = U_loc(Nsnow_loc-k)
            histowet_loc(Nsnow_loc-k+1) = histowet_loc(Nsnow_loc-k)
          end do
        end if
        Ds_loc(1) = Ds0(i,j)
        Sice_loc(1) = Sice0(i,j)
        Sliq_loc(1) = 0
        Tsnow_loc(1) = Tsnow0
        if (fsnow(i,j) > epsilon(fsnow)) then
          csnow_loc(1) = (Sice_loc(1)*hcap_ice + Sliq_loc(1)*hcap_wat) / fsnow(i,j)
          U_loc(1) = csnow_loc(1)*(Tsnow_loc(1) - Tm)
        else
          U_loc(1) = 0
        end if
        rgrn_loc(1) = rgr0
        histowet_loc(1) = 0
      end if

      ! Step 2: Initialise new layers
      if (snowthickness < epsilon(snowthickness)) then
        Ds_loc(:) = 0
        Sice_loc(:) = 0
        Sliq_loc(:) = 0
        Nsnow_loc = 0
        U_loc(:) = 0
        histowet_loc(:) = 0
      end if
      rgrn_loc(:) = 0
      Tsnow_loc(:) = Tm

      if (snowthickness >= epsilon(snowthickness)) then

        ! Step 3: Restrict surface fine snow layering to Ds_surflay if this thickness is exceeded by the fresh snow addition
        if (Nsnow_loc > 1) then
          if (snowthickness - Ds_loc(Nsnow_loc) - Ds_surflay > epsilon(snowthickness)) then
            Dtemp_surflay = Ds_loc(1)
            k_surflay = 1
            do while (Ds_surflay - Dtemp_surflay > epsilon(Dtemp_surflay))
              k_surflay = k_surflay + 1
              Dtemp_surflay = Dtemp_surflay + Ds_loc(k_surflay)
            end do
            Ds_excess = Dtemp_surflay - Ds_surflay
            Ds_old = Ds_loc(k_surflay)
            Sice_old = Sice_loc(k_surflay)
            Sliq_old = Sliq_loc(k_surflay)
            U_old = U_loc(k_surflay)
            Ds_loc(k_surflay) = Ds_old - Ds_excess
            Sice_loc(k_surflay) = Sice_old * (Ds_loc(k_surflay) / Ds_old)
            Sliq_loc(k_surflay) = Sliq_old * (Ds_loc(k_surflay) / Ds_old)
            U_loc(k_surflay) = U_old * (Ds_loc(k_surflay) / Ds_old)
            ! histowet_loc(k_surflay) unchanged
            Ds_loc(k_surflay+1) = snowthickness - Ds_surflay
            histowet_loc(k_surflay+1) = (sum((Sice_loc(k_surflay+1:Nsnow_loc) + Sliq_loc(k_surflay+1:Nsnow_loc)) &
                                        * histowet_loc(k_surflay+1:Nsnow_loc)) &
                                        + (Sice_old + Sliq_old) * (1 - Ds_loc(k_surflay) / Ds_old) &
                                        * histowet_loc(k_surflay)) &
                                        / (sum(Sice_loc(k_surflay+1:Nsnow_loc) + Sliq_loc(k_surflay+1:Nsnow_loc)) &
                                        + (Sice_old + Sliq_old) * (1 - Ds_loc(k_surflay) / Ds_old))
            Sice_loc(k_surflay+1) = sum(Sice_loc(k_surflay+1:Nsnow_loc)) + Sice_old * (1 - Ds_loc(k_surflay) / Ds_old)
            Sliq_loc(k_surflay+1) = sum(Sliq_loc(k_surflay+1:Nsnow_loc)) + Sliq_old * (1 - Ds_loc(k_surflay) / Ds_old)
            U_loc(k_surflay+1) = sum(U_loc(k_surflay+1:Nsnow_loc)) + U_old * (1 - Ds_loc(k_surflay) / Ds_old)
            if (Nsnow_loc > k_surflay + 1) then
              do k = k_surflay + 2 , Nsnow_loc
                Ds_loc(k) = 0
                Sice_loc(k) = 0
                Sliq_loc(k) = 0
                U_loc(k) = 0
                histowet_loc(k) = 0
              end do      
              Nsnow_loc = k_surflay + 1
            end if
          end if
        end if

        ! Step 4: If one layer is too thin, merge it with the neighbouring layer of closest density
        do while (Ds_min - minval(Ds_loc(1:Nsnow_loc)) > epsilon(Ds_loc)  .and. Nsnow_loc > 1)
          kmin = minloc(Ds_loc(1:Nsnow_loc),DIM=1)
          if (kmin == 1) then
            Ds_loc(1) = Ds_loc(1) + Ds_loc(2)
            histowet_loc(1) = ((Sice_loc(1) + Sliq_loc(1)) * histowet_loc(1) + (Sice_loc(2) + Sliq_loc(2)) * histowet_loc(2)) / &
                              (Sice_loc(1) + Sliq_loc(1) + Sice_loc(2) + Sliq_loc(2))
            Sice_loc(1) = Sice_loc(1) + Sice_loc(2)
            Sliq_loc(1) = Sliq_loc(1) + Sliq_loc(2)
            U_loc(1) = U_loc(1) + U_loc(2)
            if (Nsnow_loc > 2) then
              do k = 2 , Nsnow_loc - 1
                Ds_loc(k) = Ds_loc(k+1)
                Sice_loc(k) = Sice_loc(k+1)
                Sliq_loc(k) = Sliq_loc(k+1)
                U_loc(k) = U_loc(k+1)
                histowet_loc(k) = histowet_loc(k+1)
              end do
            end if
            Ds_loc(Nsnow_loc) = 0
            Sice_loc(Nsnow_loc) = 0
            Sliq_loc(Nsnow_loc) = 0
            U_loc(Nsnow_loc) = 0
            histowet_loc(Nsnow_loc) = 0
            Nsnow_loc = Nsnow_loc - 1
          else if (kmin == Nsnow_loc) then
            Ds_loc(Nsnow_loc-1) = Ds_loc(Nsnow_loc-1) + Ds_loc(Nsnow_loc)
            histowet_loc(Nsnow_loc-1) = ((Sice_loc(Nsnow_loc-1) + Sliq_loc(Nsnow_loc-1)) * histowet_loc(Nsnow_loc-1) &
                                        + (Sice_loc(Nsnow_loc) + Sliq_loc(Nsnow_loc)) * histowet_loc(Nsnow_loc)) &
                                        / (Sice_loc(Nsnow_loc-1) + Sliq_loc(Nsnow_loc-1) &
                                        + Sice_loc(Nsnow_loc) + Sliq_loc(Nsnow_loc))
            Sice_loc(Nsnow_loc-1) = Sice_loc(Nsnow_loc-1) + Sice_loc(Nsnow_loc)
            Sliq_loc(Nsnow_loc-1) = Sliq_loc(Nsnow_loc-1) + Sliq_loc(Nsnow_loc)
            U_loc(Nsnow_loc-1) = U_loc(Nsnow_loc-1) + U_loc(Nsnow_loc)
            Ds_loc(Nsnow_loc) = 0
            Sice_loc(Nsnow_loc) = 0
            Sliq_loc(Nsnow_loc) = 0
            U_loc(Nsnow_loc) = 0
            histowet_loc(Nsnow_loc) = 0
            Nsnow_loc = Nsnow_loc - 1
          else
            kup = kmin - 1
            kdown = kmin + 1
            rho_kup = (Sice_loc(kup) + Sliq_loc(kup)) / Ds_loc(kup) / fsnow(i,j)
            rho_kdown = (Sice_loc(kdown) + Sliq_loc(kdown)) / Ds_loc(kdown) / fsnow(i,j)
            rho_kmin = (Sice_loc(kmin) + Sliq_loc(kmin)) / Ds_loc(kmin) / fsnow(i,j)
            if (abs(rho_kmin - rho_kup) < abs(rho_kmin - rho_kdown)) then
              ! Layer with closest density is up
              Ds_loc(kup) = Ds_loc(kup) + Ds_loc(kmin)
              histowet_loc(kup) = ((Sice_loc(kup) + Sliq_loc(kup)) * histowet_loc(kup) &
                                  + (Sice_loc(kmin) + Sliq_loc(kmin)) * histowet_loc(kmin)) &
                                  / (Sice_loc(kup) + Sliq_loc(kup) &
                                  + Sice_loc(kmin) + Sliq_loc(kmin))
              Sice_loc(kup) = Sice_loc(kup) + Sice_loc(kmin)
              Sliq_loc(kup) = Sliq_loc(kup) + Sliq_loc(kmin)
              U_loc(kup) = U_loc(kup) + U_loc(kmin)
              do k = kmin , Nsnow_loc - 1
                Ds_loc(k) = Ds_loc(k+1)
                Sice_loc(k) = Sice_loc(k+1)
                Sliq_loc(k) = Sliq_loc(k+1)
                U_loc(k) = U_loc(k+1)
                histowet_loc(k) = histowet_loc(k+1)
              end do
              Ds_loc(Nsnow_loc) = 0
              Sice_loc(Nsnow_loc) = 0
              Sliq_loc(Nsnow_loc) = 0
              U_loc(Nsnow_loc) = 0
              histowet_loc(Nsnow_loc) = 0
              Nsnow_loc = Nsnow_loc - 1
            else
              ! Layer with closest density is down
              Ds_loc(kmin) = Ds_loc(kmin) + Ds_loc(kdown)
              histowet_loc(kmin) = ((Sice_loc(kmin) + Sliq_loc(kmin)) * histowet_loc(kmin) &
                                   + (Sice_loc(kdown) + Sliq_loc(kdown)) * histowet_loc(kdown)) &
                                   / (Sice_loc(kmin) + Sliq_loc(kmin) &
                                   + Sice_loc(kdown) + Sliq_loc(kdown))
              Sice_loc(kmin) = Sice_loc(kmin) + Sice_loc(kdown)
              Sliq_loc(kmin) = Sliq_loc(kmin) + Sliq_loc(kdown)
              U_loc(kmin) = U_loc(kmin) + U_loc(kdown)
              if (kdown < Nsnow_loc) then
                do k = kdown , Nsnow_loc - 1
                  Ds_loc(k) = Ds_loc(k+1)
                  Sice_loc(k) = Sice_loc(k+1)
                  Sliq_loc(k) = Sliq_loc(k+1)
                  U_loc(k) = U_loc(k+1)
                  histowet_loc(k) = histowet_loc(k+1)
                end do
              end if
              Ds_loc(Nsnow_loc) = 0
              Sice_loc(Nsnow_loc) = 0
              Sliq_loc(Nsnow_loc) = 0
              U_loc(Nsnow_loc) = 0
              histowet_loc(Nsnow_loc) = 0
              Nsnow_loc = Nsnow_loc - 1
            end if
          end if
        end do

        ! Step 5: If too many layers, merge the two neighbours with closest density
        do while(Nsnow_loc > Nsmax)
          rho(1:Nsnow_loc) = (Sice_loc(1:Nsnow_loc) + Sliq_loc(1:Nsnow_loc)) / Ds_loc(1:Nsnow_loc) / fsnow(i,j)
          do k = 1 , Nsnow_loc - 1
            diff_rho(k) = abs(rho(k) - rho(k+1))
           end do
          kmerge = minloc(diff_rho(1:Nsnow_loc-1),DIM=1)
          Ds_loc(kmerge) = Ds_loc(kmerge) + Ds_loc(kmerge+1)
          histowet_loc(kmerge) = ((Sice_loc(kmerge) + Sliq_loc(kmerge)) * histowet_loc(kmerge) &
                                 + (Sice_loc(kmerge+1) + Sliq_loc(kmerge+1)) * histowet_loc(kmerge+1)) &
                                 / (Sice_loc(kmerge) + Sliq_loc(kmerge) &
                                 + Sice_loc(kmerge+1) + Sliq_loc(kmerge+1))
          Sice_loc(kmerge) = Sice_loc(kmerge) + Sice_loc(kmerge+1)
          Sliq_loc(kmerge) = Sliq_loc(kmerge) + Sliq_loc(kmerge+1)
          U_loc(kmerge) = U_loc(kmerge) + U_loc(kmerge+1)
          if (kmerge + 1 < Nsnow_loc) then
            do k = kmerge + 1 , Nsnow_loc - 1
              Ds_loc(k) = Ds_loc(k+1)
              Sice_loc(k) = Sice_loc(k+1)
              Sliq_loc(k) = Sliq_loc(k+1)
              U_loc(k) = U_loc(k+1)
              histowet_loc(k) = histowet_loc(k+1)
            end do
          end if
          Ds_loc(Nsnow_loc) = 0
          Sice_loc(Nsnow_loc) = 0
          Sliq_loc(Nsnow_loc) = 0
          U_loc(Nsnow_loc) = 0
          histowet_loc(Nsnow_loc) = 0
          Nsnow_loc = Nsnow_loc - 1
        end do

        ! Step 6: If more layers could be used, split the thickest ones
        if (Nsnow_loc > 0) then
          do while(Nsnow_loc < Nsmax)
            kmax = maxloc(Ds_loc(1:Nsnow_loc),DIM=1)
            2 continue
            if (Ds_loc(kmax) - 2.0 * Ds_min < epsilon(Ds_loc)) then
              goto 3
            end if
            if (kmax == Nsnow_loc) then
              if (Nsnow_loc > 1) then
                if (Ds_surflay - (sum(Ds_loc(1:kmax-1))) - Ds_min > epsilon(Ds_loc)) then
                  if (Ds_loc(kmax) - (Ds_surflay - (sum(Ds_loc(1:kmax-1)))) - Ds_min > epsilon(Ds_loc)) then
                    Nsnow_loc = Nsnow_loc + 1
                    wt = (Ds_surflay - (sum(Ds_loc(1:kmax-1)))) / Ds_loc(kmax) ! Ratio of layer taken away
                    Ds_loc(Nsnow_loc) = (1 - wt) * Ds_loc(kmax)
                    Sice_loc(Nsnow_loc) = (1 - wt) * Sice_loc(kmax)
                    Sliq_loc(Nsnow_loc) = (1 - wt) * Sliq_loc(kmax)
                    histowet_loc(Nsnow_loc) = histowet_loc(kmax)
                    U_loc(Nsnow_loc) = (1 - wt) * U_loc(kmax)
                    Ds_loc(kmax) = wt * Ds_loc(kmax)
                    Sice_loc(kmax) = wt * Sice_loc(kmax)
                    Sliq_loc(kmax) = wt * Sliq_loc(kmax)
                    histowet_loc(kmax) = histowet_loc(kmax)
                    U_loc(kmax) = wt * U_loc(kmax)
                  else ! more than Ds_min is missing in the "surface layers" 
                       ! and Ds_loc(kmax) > 2 * Ds_min anyway
                    Nsnow_loc = Nsnow_loc + 1
                    wt = Ds_min / Ds_loc(kmax) ! Ratio of layer taken away
                    Ds_loc(Nsnow_loc) = (1 - wt) * Ds_loc(kmax)
                    Sice_loc(Nsnow_loc) = (1 - wt) * Sice_loc(kmax)
                    Sliq_loc(Nsnow_loc) = (1 - wt) * Sliq_loc(kmax)
                    histowet_loc(Nsnow_loc) = histowet_loc(kmax)
                    U_loc(Nsnow_loc) = (1 - wt) * U_loc(kmax)
                    Ds_loc(kmax) = wt * Ds_loc(kmax)
                    Sice_loc(kmax) = wt * Sice_loc(kmax)
                    Sliq_loc(kmax) = wt * Sliq_loc(kmax)
                    histowet_loc(kmax) = histowet_loc(kmax)
                    U_loc(kmax) = wt * U_loc(kmax)
                  end if
                else
                  kmax = maxloc(Ds_loc(1:Nsnow_loc-1),DIM=1)
                  goto 2
                end if
              else ! Nsnow_loc = 1
                  Nsnow_loc = Nsnow_loc + 1
                  Ds_loc(Nsnow_loc) = Ds_loc(kmax) / 2
                  Sice_loc(Nsnow_loc) = Sice_loc(kmax) / 2
                  Sliq_loc(Nsnow_loc) = Sliq_loc(kmax) / 2
                  histowet_loc(Nsnow_loc) = histowet_loc(kmax)
                  U_loc(Nsnow_loc) = U_loc(kmax) / 2
                  Ds_loc(kmax) = Ds_loc(kmax) / 2
                  Sice_loc(kmax) = Sice_loc(kmax) / 2
                  Sliq_loc(kmax) = Sliq_loc(kmax) / 2
                  histowet_loc(kmax) = histowet_loc(kmax)
                  U_loc(kmax) = U_loc(kmax) / 2
              end if
            else
              Nsnow_loc = Nsnow_loc + 1
              do k = Nsnow_loc, kmax + 2, -1
                Ds_loc(k) = Ds_loc(k-1)
                Sice_loc(k) = Sice_loc(k-1)
                Sliq_loc(k) = Sliq_loc(k-1)
                histowet_loc(k) = histowet_loc(k-1)
                U_loc(k) = U_loc(k-1)
              end do
              Ds_loc(kmax+1) = Ds_loc(kmax) / 2
              Sice_loc(kmax+1) = Sice_loc(kmax) / 2
              Sliq_loc(kmax+1) = Sliq_loc(kmax) / 2
              histowet_loc(kmax+1) = histowet_loc(kmax)
              U_loc(kmax+1) = U_loc(kmax) / 2
              Ds_loc(kmax) = Ds_loc(kmax) / 2
              Sice_loc(kmax) = Sice_loc(kmax) / 2
              Sliq_loc(kmax) = Sliq_loc(kmax) / 2
              histowet_loc(kmax) = histowet_loc(kmax)
              U_loc(kmax) = U_loc(kmax) / 2
            end if
          end do
        end if
        3 continue

        ! Step 7: Diagnose snow layer temperatures
        do k = 1, Nsnow_loc
          if (fsnow(i,j) > epsilon(fsnow)) then
            csnow_loc(k) = (Sice_loc(k)*hcap_ice + Sliq_loc(k)*hcap_wat) / fsnow(i,j)
            Tsnow_loc(k) = Tm + U_loc(k) / csnow_loc(k)
          else
            Tsnow_loc(k) = Tm
          end if
          ! Bring back histowet to the [0,1] range if it is out by epsilon
          histowet_loc(k) = min(max(histowet_loc(k), 0.0), 1.0)
        end do

      end if

      ! Step 8: Copy local variables to state variables
      Sice(:,i,j) = Sice_loc(1:Nsmax)
      Sliq(:,i,j) = Sliq_loc(1:Nsmax)
      Ds(:,i,j) = Ds_loc(1:Nsmax)
      histowet(:,i,j) = histowet_loc(1:Nsmax)
      rgrn(:,i,j) = rgrn_loc(1:Nsmax)
      Tsnow(:,i,j) = Tsnow_loc(1:Nsmax)
      Nsnow(i,j) = Nsnow_loc

    end if

    1 continue

  end do
end do

end subroutine SNOW_LAYERING
