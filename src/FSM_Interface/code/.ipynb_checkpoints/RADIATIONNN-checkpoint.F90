!-----------------------------------------------------------------------
! Surface and canopy net shortwave radiation
!-----------------------------------------------------------------------
subroutine RADIATIONNN(alb,SWsrf,SWveg,Sdirt,Sdift,asrf_out,SWsci)

use MODCONF, only: CANMOD, RADSBG, ALBEDO

use CONSTANTS, only: &
  I0,                &! Solar constant (W/m^2)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  year,              &! Year
  month,             &! Month of year
  day,               &! Day of month
  hour,              &! Hour of day
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
  Sf,                &! Snowfall rate (kg/m2/s)
  Sf24h,             &! Snowfall 24hr (kg/m2)
  Ta,                &! Air temperature (K)
  Tv                  ! Time-varying transmissivity for dSWR

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  asmx,              &! Maximum albedo for fresh snow
  asmn,              &! Minimum albedo for melting snow
  avg0,              &! Snow-free vegetation albedo
  avgs,              &! Snow-covered vegetation albedo
  Salb,              &! Snowfall to refresh albedo (kg/m^2)
  Talb,              &! Albedo decay temperature threshold (C)
  tcld,              &! Cold snow albedo decay timescale (s)
  tmlt,              &! Melting snow albedo decay timescale (s)
  fsar,              &! Albedo adjustment range for canopy dependence 
  fthresh             ! Forest fraction required for forest tile to be considered 

use PARAMMAPS, only: &
  alb0,              &! Snow-free ground albedo
  fsky,              &! Sky view fraction
  scap,              &! Canopy snow capacity (kg/m^2)
  trcn                ! Canopy transmissivity

use STATE_VARIABLES, only: &
  albs,              &! Snow albedo
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  fsnow,             &! Snow cover fraction
  Sveg,              &! Snow mass on vegetation (kg/m^2)
  Tsnow,             &! Snow layer temperatures (K)
  Tsrf                ! Surface skin temperature (K)

use LANDUSE, only: &
  fsky_terr,         &! Terain sky view factor
  fveg,              &! Canopy cover fraction
  dem,               &! Terrain elevation (m)
  forest              ! Grid cell forest fraction

implicit none

real, intent(inout) :: &
  alb(Nx,Ny),        &! Albedo
  Sdirt(Nx,Ny),      &! Incoming direct beam radiation corrected for subgrid topography (W/m^2)
  Sdift(Nx,Ny),      &! Incoming diffuse radiation corrected for subgrid topography (W/m^2)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny)        ! Net SW radiation absorbed by vegetation (W/m^2)

integer :: &
  i,j                 ! Point counters

real, intent(out) :: &
  asrf_out(Nx,Ny),   &! Surface albedo
  SWsci(Nx,Ny)        ! Subcanopy incoming SWR (W/m^2)

real :: &
  alim,              &! Limiting snow albedo
  acan,              &! Canopy albedo
  asrf,              &! Surface albedo
  aveg,              &! Vegetation albedo
  fcans,             &! Canopy snowcover fraction
  rt,                &! Reciprocal timescale for albedo adjustment (1/s)
  SWEtmp,            &! Temporary snow water equivalent for snow covered fraction calculation (kg/m^2)
  Sdif_aux,          &! Auxiliary variable containing diffuse SW radiation relevant to canopy transfer
  tau,               &! Snow albedo decay timescale (s)
  tdif,              &! Canopy transmissivity for diffuse radiation
  tdir,              &! Canopy transmissivity for direct-beam radiation
  SWtopo_out,        &! Outgoing SW radiation corrected for subgrid topography (W/m^2)
  Sun_elev,          &! Solar elevation angle
  adm,               &! Melting snow albedo decay time (h)
  afs                 ! Fresh snow albedo
  
real, parameter :: &  ! Albedo model 2 (JIM) parameters
  adc = 3000,    &! Cold snow albedo decay time (h), original JIM: 3000, FSM: 1000
  Sfmin = 10      ! Minimum snowfall to refresh albedo (kg/m^2)

! Snow albedo
do j = 1, Ny
do i = 1, Nx

  if (isnan(dem(i,j))) goto 1 ! Exclude points outside of the domain

  if (CANMOD == 1 .and. forest(i,j) < fthresh) goto 1

  if (ALBEDO == 0) then
  ! Diagnostic
    albs(i,j) = asmn + (asmx - asmn)*(Tsrf(i,j) - Tm) / Talb
    if (albs(i,j) < min(asmx, asmn)) albs(i,j) = min(asmx, asmn)
    if (albs(i,j) > max(asmx, asmn)) albs(i,j) = max(asmx, asmn)

  elseif (ALBEDO == 1) then
  ! Prognostic
    tau = tcld
    if (Tsrf(i,j) >= Tm) tau = tmlt
    rt = 1/tau + Sf(i,j)/Salb
    alim = (asmn/tau + Sf(i,j)*asmx/Salb)/rt
    albs(i,j) = alim + (albs(i,j) - alim)*exp(-rt*dt)
    if (albs(i,j) < min(asmx, asmn)) albs(i,j) = min(asmx, asmn)
    if (albs(i,j) > max(asmx, asmn)) albs(i,j) = max(asmx, asmn)

  else ! ALBEDO == 2
  ! Prognostic, tuned, copied from JIM
    SWEtmp = sum(Sice(:,i,j) + Sliq(:,i,j))
    !! DECAY RATES
    ! Melting snow decay time (was fixed 100)
    if (month > 6 .AND. month < 10) then
      adm = 50
    else
      adm = 130
    end if
    ! New Snow albedo (was fixed 0.85)
    ! 11/2021 tuning: high elevation afs changed from 0.86 to 0.92
    if (dem(i,j) >= 2300) then
      afs = 0.92
    else if (dem(i,j) <= 1500) then
      afs = 0.80
    else
      afs = 0.92 + (2300 - dem(i,j)) / (2300 - 1500) * (0.80 - 0.92)
    end if
    if (Tsrf(i,j) >= Tm .AND. Tsnow(1,i,j) >= Tm) then  ! was only based on Tss
      albs(i,j) = (albs(i,j) - asmn)*exp(-(dt/3600)/adm) + asmn
    else
      albs(i,j) = albs(i,j) - (dt/3600)/adc
    end if
    if (SWEtmp < 75.0) then ! more stuff showing on and up through snow
      afs = afs * 0.80
    end if
    ! Reset to fresh snow albedo (wasn't originally available; only else term)
    if ((Sf(i,j) * dt) > 0.0 .AND. Sf24h(i,j) > Sfmin) then
      albs(i,j) = afs
    else
      albs(i,j) = albs(i,j) + (afs - albs(i,j))*Sf(i,j)*dt/Sfmin
    end if
    !! End Adjustments
      if (albs(i,j) > afs) albs(i,j) = afs
      if (albs(i,j) < asmn) albs(i,j) = asmn
  endif

  1 continue 
  
end do
end do

! Surface and canopy net shortwave radiation
do j = 1, Ny
do i = 1, Nx

  if (isnan(dem(i,j))) goto 2 ! Exclude points outside of the domain

  if (CANMOD == 1 .and. forest(i,j) < fthresh) goto 2 ! Exclude points that have no forest

  ! Surface albedo
  asrf = albs(i,j) 
  asrf = albs(i,j)*(1-fveg(i,j)*fsar)  ! simple snow albedo dependence on canopy density to account e.g. for litter
  if (fsnow(i,j) == 0) asrf = alb0(i,j)

! Partial snowcover on canopy
  fcans = 0
  if (scap(i,j) > epsilon(scap)) fcans = Sveg(i,j) / scap(i,j)
  aveg = (1 - fcans)*avg0 + fcans*avgs
  acan = fveg(i,j)*aveg
  ! Canopy surface albedo for computing terrain radiation over canopy
  alb(i,j) = fveg(i,j)*aveg + (1-fveg(i,j))*asrf
  
  ! Surface albedo is storef in asurf_out to write in results
  asrf_out(i,j) = alb(i,j)
  
  if (RADSBG == 1) then
    ! Call Subgrid parameterization for SW radiation to compute SWtopo,netto SWtn
    call SWRADTOPO(alb(i,j),Sdir(i,j),Sdif(i,j),SWsrf(i,j),Sdirt(i,j),Sdift(i,j),SWtopo_out,Sun_elev,year,month,day,hour,i,j)
  endif

  if (RADSBG == 0) then
    SWtopo_out =  alb(i,j)*(Sdir(i,j)+Sdif(i,j))
    Sdirt(i,j) = Sdir(i,j)
    Sdift(i,j) = Sdif(i,j) 
  endif

  ! Solar radiation trasmission + thermal emissions from surroundings
  if (CANMOD == 0) then
    SWveg(i,j) = 0
    SWsrf(i,j) = (1 - alb(i,j))*(Sdir(i,j)+Sdif(i,j))
    ! LW(i,j) = fsky(i,j)*LW(i,j) + (1 - fsky(i,j))*sb*Ta(i,j)**4
    LW(i,j) = fsky_terr(i,j)*LW(i,j) + (1 - fsky_terr(i,j))*sb*Ta(i,j)**4
  endif

  if (CANMOD == 1) then
    Sdif_aux = fsky(i,j)/fsky_terr(i,j)*Sdift(i,j)
    tdif = trcn(i,j)
    tdir = Tv(i,j) 

    ! Effective albedo and net radiation
    alb(i,j) = acan + (1 - acan)*asrf*tdif**2
    if (Sdift(i,j) + Sdirt(i,j) > epsilon(Sdift(i,j)))  & 
      alb(i,j) = (acan*(Sdif_aux+tdir*Sdirt(i,j)) + asrf*tdif*(tdif*Sdif_aux+tdir*Sdirt(i,j))) / &
                (Sdif_aux + Sdirt(i,j))
    SWsrf(i,j) = (1 - asrf)*(tdif*Sdif_aux + tdir*Sdirt(i,j))
    SWveg(i,j) = ((1-tdif)*(1-aveg)+tdif*asrf*(1-tdif))*Sdif_aux + &   
                  (tdir*fveg(i,j)*(1-aveg)+tdir*asrf*(1-tdif))*Sdirt(i,j)   ! local SWR absorption by vegetation correlates with local tdir  
    SWsci(i,j) = tdif*Sdif_aux + tdir*Sdirt(i,j)
    ! Terrain LWR if not calculated later
    if (fveg(i,j) == 0) LW(i,j) = fsky_terr(i,j)*LW(i,j) + (1 - fsky_terr(i,j))*sb*Ta(i,j)**4      
  endif 
            
  2 continue 
  
end do
end do

end subroutine RADIATIONNN
