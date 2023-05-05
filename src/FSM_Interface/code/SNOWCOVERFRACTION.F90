!--------------------------------------------------------------------------
! Compute Snow covered fraction with 
!      - parameterized standard deviation of snow depth at peak of winter (Helbig et al., 2015)
!      - tanh for fsnow as in Helbig et al., 2015
!      - with a snow depth history at each pixel in the season originally suggested from Jan Magnusson
! Author: Nora Helbig & Michael Schirmer & Jan Magnusson
!--------------------------------------------------------------------------

subroutine SNOWCOVERFRACTION(snowdepth,SWEtmp,i,j)

use MODCONF, only: SNFRAC

use PARAMETERS, only: &
  hfsn                        ! Snowcover fraction depth scale (m)

use STATE_VARIABLES, only: &
  fsnow,             &! Snow cover fraction 
  swehist,           &! history of Snow water equivalent during last 14 days (kg/m^2)
  swemin,            &! Minimum swe during the season (m)
  swemax,            &! Maximum swe during the season (m)
  snowdepthhist,     &! history of snow depth during last 14 days (m)
  snowdepthmin,      &! Minimum Snow depth at time step of swemin (m)
  snowdepthmax        ! Maximum Snow depth at time stemp of swemax(m)

use LANDUSE, only: &
  slopemu,           &! slope parameter 
  xi,                &! terrain correlation length
  Ld                  ! Model domain size (e.g. 1000m grid cell size)
  
use DRIVING, only: &
  hour                ! By convention, simulations start at 6am and the daily history of HS and SWE should be taken at this time.

implicit none

integer :: &
  iabsmin,iabsmax,irecentmin, &! Indices of extrema found in SWEbuffer and applied to snowdepthbuffer
  iloop                        ! do loop index

integer, intent(in) :: &
  i,j                        ! Point counters

real, intent(in) :: &
  snowdepth,                &! Snow depth (m)
  SWEtmp                     ! Temporary snow water equivalent for snow covered fraction calculation (kg/m^2)

real :: &
  sd_snowdepth0,             &! Parameterized standard deviation of snow depths in a grid cell cf. Helbig et al. (2014)
  sd_snowdepth1,             &! Standard deviation parameter cf. Helbig et al. (2014)
  sd_snowdepth2,             &! Standard deviation parameter cf. Helbig et al. (2014)
  sd_snowdepth3,             &! Standard deviation parameter cf. Helbig et al. (2014)
  snowdepth_threshold,       &! snowdepth threshold for setting fsnow to zero
  dsnowdepth,                &! change in Snow depth in past 14 days (m)
  dsnowdepthmax,             &! Maximum change in Snow depth in past 14 days (m)
  dsnowdepth_recent,         &! Recent change in Snow depth (m)
  SWEbuffer(15),             &! buffer with 14 day history of swe plus hourly current SWEtmp
  snowdepthbuffer(15),       &! buffer with 14 day history of snowdepth plus hourly current snowdepth
  diffSWEbuffer(14),         &! difference betwen entries in SWEbuffer
  snowdepthmin_buffer,       &! absolute minimum of snowdepth within the buffer
  snowdepthmax_buffer,       &! absolute maximum of snowdepth within the buffer
  snowdepthmin_recent,       &! recent minimum of snowdepth within the buffer
  fsnow_season,              &! seasonal snow covered fraction
  fsnow_nsnow,               &! new snow snow covered fraction
  fsnow_nsnow_recent,        &! recent new snow snow covered fraction
  coeff_vari,                &! coefficient of variation of snow depth
  sd_snowdepth0_dhs,         &! Parameterized standard deviation of new snow depths in a grid cell 
  sd_snowdepth0_dhs_recent    ! Parameterized standard deviation of recent new snow depths in a grid cell 
  
! lower snowdepth value which is taken from swe_threshold = 2 (in EKF) 
! converted using the standard_density = 350 (in EKF) to snowdepth_threshold = 0.005714286
snowdepth_threshold = 0.005714286

if (SNFRAC == 0) then
  ! reads SWE amounts from previous model runs, and calculate dswe in as difference between the minimum in the previous days to current swe
  ! This has three SCF regimes: 
  ! (1) scf_season calculated with swemax (for the topographic dependent standard deviation (Helbig et al., 2015)), swemin (for the nominator in the tanh function) 
  ! (2) scf_nsnow calculated with a SWE amounts of previous model runs and the current swe (SWEbuffer), swemin_buffer is the minimum in this buffer, 
  !     swemax_buffer is the maximum from the minimun to current. These values are used to calculate dswemax (for standard deviation, here Egli and Jonas 2009 is used 
  !     to account for new snow is behaving more like a flat field) and to calculate dswe (for the nominator in the tanh function). 
  ! (3) scf_nsnow_recent ensures that pixels with a SWE history have similar new snow SCF as pixels which were bare and not less. Here dswe_recent is used for both standard
  !     deviation (as well Egli and Jonas, 2009) and nominator, no dswemax_recent. Here the most recent minimum is searched in the SWEbuffer
  ! Snowfalls during melt may add variability. This is accounted for with Jan's idea of backwards re-calculating swemax with current swe. This swemax is the value as if
  !     this current swe would have been reached with melting only and no intermediate snow falls. cf. the 'ireset' part.
  ! The maximum of the three SCF regimes is used to calculate the final SCF. This is Jan's great idea to ensure that during the ablation period a new snow fall is adding significantly
  !     to the SCF, but allows to fall back to similar SCF values when this new snow has melted away. dswe used for scf_nsnow can be interpreted as if the last snowfall would fall on bare ground. 
  !     After melting away (i.e. scf_nsnow<scf_season), the SCF will switch from scf_nsnow to scf_season
  ! Important is that swemax, swemin is used without settling, i.e. on SWE and not HS. Also important is that Helbig's and Egli's standard deviations are calculated on HS. 

  ! calculate topo terms needed for standard deviation of snow depth (done)
  sd_snowdepth1 = exp(-1 / (Ld(i,j)/xi(i,j))**2)
  sd_snowdepth3 = slopemu(i,j)**0.309

  ! merge current SWEtmp with SWEtmp history from past 14 days
  SWEbuffer(1) = SWEtmp
  SWEbuffer(2:15) = swehist(i,j,:)
  snowdepthbuffer(1) = snowdepth
  snowdepthbuffer(2:15) = snowdepthhist(i,j,:)

  ! calculate snowdepthmin_buffer, snowdepthmax_buffer, snowdepthmin_recent 
  ! find indices of global min and max in SWEbuffer
  iabsmax = maxloc(SWEbuffer,DIM=1)
  iabsmin = minloc(SWEbuffer,DIM=1)

  ! find index of recent min in SWEbuffer
  ! calculate diff vector of SWEBuffer
  do iloop = 1, 14
    diffSWEbuffer = SWEbuffer(iloop+1)-SWEbuffer(iloop)
    if (diffSWEbuffer(iloop) > 0.5) then
      EXIT
    end if
  end do
  irecentmin = minloc(SWEbuffer(1:iloop),DIM=1)

  ! use indices to determine snowdepth amounts
  snowdepthmin_buffer = snowdepthbuffer(iabsmin)
  snowdepthmax_buffer = snowdepthbuffer(iabsmax)
  snowdepthmin_recent = snowdepthbuffer(irecentmin)

  ! Compute storage of new snow on old snow in snowdepthbuffer 
  dsnowdepth = snowdepth - snowdepthmin_buffer
  if (dsnowdepth < epsilon(dsnowdepth)) then
    dsnowdepth = 0
  end if

  ! compute dswemax in SWEbuffer 
  dsnowdepthmax = snowdepthmax_buffer - snowdepthmin_buffer
  if (dsnowdepthmax < epsilon(dsnowdepthmax)) then
    dsnowdepthmax = 0
  end if

  ! don't accept dsnowdepthmax to be larger then dsnowdepth, otherwise larger fnsnow values
  ! todo: think about doing this for snowdepthmin and snowdepthmax as well, and swemin and swemax
  if (dsnowdepthmax < dsnowdepth) then
    dsnowdepthmax = dsnowdepth
  end if

  ! Compute storage of recent new snow on old snow in SWEbuffer (done)
  dsnowdepth_recent = snowdepth - snowdepthmin_recent
  if (dsnowdepth_recent < epsilon(dsnowdepth_recent)) then
    dsnowdepth_recent = 0
  end if

  !!! state variables interpeting the whole SWEtmp history, not only the past 14 days in the buffer 
  ! Set swemax and swemin equal to zero if no snow, same with corresponding snow depth values 
  if (SWEtmp < epsilon(SWEtmp)) then
    swemax(i,j) = 0
    swemin(i,j) = 0
  end if
  if (snowdepth < epsilon(snowdepth)) then
    snowdepthmax(i,j) = 0
    snowdepthmin(i,j) = 0
  end if

  ! Set swemax and swemin equal to SWEtmp if maximum, store also snowdepthmax and snowdepthmin of those time steps 
  if (SWEtmp >= swemax(i,j)) then
    swemax(i,j)       = SWEtmp
    swemin(i,j)       = SWEtmp
  end if
  if (snowdepth >= snowdepthmax(i,j)) then
    snowdepthmax(i,j) = snowdepth
    snowdepthmin(i,j) = snowdepth
  end if

  ! Set swemin equal SWEtmp if smaller than swemin, same with corresponding snow depth value 
  if (SWEtmp < swemax(i,j) .and. SWEtmp  < swemin(i,j)) then
    swemin(i,j) = SWEtmp
  end if
  if (snowdepth < snowdepthmax(i,j) .and. snowdepth  < snowdepthmin(i,j)) then
    snowdepthmin(i,j) = snowdepth
  end if

  !!! calculating SCF
  ! Initial guess of snow covered fraction 
  fsnow_season       = 0
  fsnow_nsnow        = 0
  fsnow_nsnow_recent = 0

  !! seasonal scf, inserting snow depth in formulas of Helbig et al. and Egli and Jonas
  ! calculate standard deviation (done)
  sd_snowdepth2 = snowdepthmax(i,j)**0.549
  sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3
  ! set completely flat pixels to values determined by Luca (instead of 1 or 0) 
  if (.not.(slopemu(i,j) > epsilon(0.0))) then
    sd_snowdepth0 = snowdepthmax(i,j)**0.84
  end if
  ! calculate snow covered fraction 
  if (snowdepthmax(i,j) > epsilon(snowdepthmax(i,j))) then
    fsnow_season = tanh(1.3 * snowdepthmin(i,j) / sd_snowdepth0)
  end if

  ! calculate cv 
  coeff_vari = sd_snowdepth0 / snowdepthmax(i,j)

  !! scf based on dswe of last 14 days
  ! calculate standard deviation of dhs, taking Luca's formula (flat field approximation) 
  sd_snowdepth0_dhs = dsnowdepthmax**0.84
  ! calculate snow covered fraction of nsnow 
  if (dsnowdepthmax > epsilon(dsnowdepthmax)) then
    fsnow_nsnow = tanh(1.3 * dsnowdepth / sd_snowdepth0_dhs)
  end if

  !! scf based on dswe_recent since last minimum
  ! calculate standard deviation of dsnowdepth_recent, taking Luca's formula (flat field approximation)
  sd_snowdepth0_dhs_recent = dsnowdepth_recent**0.84
  ! calculate snow covered fraction of nsnow with recent dswe, converting SWEtmp into snow depth
  if (dsnowdepth_recent > epsilon(dsnowdepth_recent)) then
    fsnow_nsnow_recent = tanh(1.3 * dsnowdepth_recent / sd_snowdepth0_dhs_recent)
  end if

  ! take maximum between the two new snow scf, similar to taking the maximum of all three regimes at the end (done)
  fsnow_nsnow = max(fsnow_nsnow,fsnow_nsnow_recent)

  ! RESET PART OF THE CODE IS TEMPORARILY REMOVED - SOLUTION TO BE FOUND TO AVOID INSTABILITIES
  !    !! recalculate scf_season if new snow has melted after a snow fall to account for a higher CV
  !    ! If new snow has melted away, update parameters swemin and swemax of
  !    ! "seasonal snow" so that the scf trajectory continues along last
  !    ! scf-value given by scf_nsnow, added that it is really (SWE yesterday > SWE current + threshold of 2) melting, including the threshold for more than spurious differences
  !    if (fsnow_nsnow .NE. 0 .and. fsnow_nsnow < fsnow_season .and. SWEbuffer(2) > SWEtmp + 2) then
  !      swemin(i,j)       = SWEtmp
  !      snowdepthmin(i,j) = snowdepth
  !      rhomax = swemax(i,j)/snowdepthmax(i,j) ! rhomax should remain constant with time, i.e the modelleded density at timestep of swemax
  !      snowdepthmax(i,j) = 1.3 * snowdepthmin(i,j) / (coeff_vari * atanh(fsnow_season))
  !      swemax(i,j) = rhomax * snowdepthmax(i,j)
  !      ! re-calculate standard deviation with new snowdepthmax
  !      sd_snowdepth2 = snowdepthmax(i,j)**0.549
  !      sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3
  !      ! set completely flat pixels to values determined by Luca (instead of 1 or 0) 
  !      if (.not.(slopemu(i,j) > epsilon(0.0))) then
  !        sd_snowdepth0 = snowdepthmax(i,j)**0.84
  !      end if
  !      ! calculate snow covered fraction
  !      fsnow_season = tanh(1.3 * snowdepthmin(i,j) / !sd_snowdepth0)
  !      fsnow_nsnow = 0
  !    end if

  ! If snow amounts are too low, set snow covered fraction to zero 
  if (snowdepthmin(i,j) < snowdepth_threshold) then
    fsnow_season = 0
  end if
  if (dsnowdepth < snowdepth_threshold) then
    fsnow_nsnow = 0
  end if

  ! Use the largest of the two fsnow estimates
  fsnow(i,j) = max(fsnow_season,fsnow_nsnow)

else if (SNFRAC == 1) then
  ! HelbigHS
  ! calculate standard deviation
  sd_snowdepth2 = snowdepth**0.549
  sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3

  ! calculate snow covered fraction 
  fsnow(i,j) = tanh(1.3 * snowdepth / sd_snowdepth0)

  ! If snow amounts are too low, set snow covered fraction to zero
  if (snowdepth < snowdepth_threshold) then
    fsnow(i,j) = 0
  end if

else if (SNFRAC == 2) then
  ! HelbigHS0 !todo: think about switching this to HelbigSWEMAX to avoid settling dependent differenes between snowdepth and snowdepthmax 
  ! Set snowdepthmax equal to zero if no snow
  if (snowdepth == 0) then
    snowdepthmax(i,j) = 0
  end if

  ! Set snowdepthmax equal to snowdepth if maximum
  if (snowdepth > snowdepthmax(i,j)) then
    snowdepthmax(i,j) = snowdepth
  end if

  ! calculate standard deviation
  sd_snowdepth2 = snowdepthmax(i,j)**0.549
  sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3

  ! calculate snow covered fraction 
  fsnow(i,j) = tanh(1.3 * snowdepth / sd_snowdepth0)

  ! If snow amounts are too low, set snow covered fraction to zero
  if (snowdepth < snowdepth_threshold) then
    fsnow(i,j) = 0
  end if

else if (SNFRAC == 3) then
  ! Point model
  if (snowdepth > epsilon(snowdepth)) then
    fsnow(i,j) = 1
  else
    fsnow(i,j) = 0
  end if

else ! SNFRAC == 4
  !tanh model / original FSM  
  fsnow(i,j) = tanh(snowdepth/hfsn)
endif

if (snowdepth < epsilon(snowdepth)) then
  fsnow(i,j) = 0
else
  fsnow(i,j) = min(fsnow(i,j), 1.)
end if

! BC update history of SWE and hs only if they correspond to 6:00am values
if (hour > 4.5 .and. hour < 5.5) then
  SWEhist(i,j,:) = SWEbuffer(1:14)
  snowdepthhist(i,j,:) = snowdepthbuffer(1:14)
end if

end subroutine SNOWCOVERFRACTION
