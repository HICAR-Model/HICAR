!--------------------------------------------------------------------------
! Include subgrid topographical influences on surface shortwave energy balance terms
! Author: Nora Helbig
!--------------------------------------------------------------------------

subroutine SWRADTOPO(alball_n,Sdir_n,Sdif_n,SWtn,Sdirt_n,Sdift_n,SWtopo_out,sun_elev,year,month,day,hour,i,j)

use CONSTANTS, only : &
  pi               ! Pi

use LANDUSE, only : &
  fsky_terr,      &! Terrain sky view fraction
  slopemu,        &! slope parameter
  lat,            &! latitude of each grid cell (center?)
  lon              ! longitude of each grid cell (center?) 

implicit none

integer*4 :: &
  year,          &! Year
  jday,          &! Julian date of actual day
  jday1,         &! Julian date of 1 January of current year
  month,         &! Month of year
  day,           &! Day of month
  met_day_number  ! Day number of the year

real, intent(in) :: &
  alball_n,      &! Snow cover fraction weighted albedo
  Sdir_n,         &! Incoming flat field direct beam shortwave radiation (W/m^2)
  Sdif_n           ! Incoming flat field diffuse sky shortwave radiation (W/m^2)

real, intent(out) :: &
  SWtn            ! Net SW radiation accounting for subgrid topography (W/m^2)
  
real*4, intent(out) :: &
  sun_elev,            &! Solar elevation angle
  SWtopo_out,          &! Outgoing SW radiation corrected for subgrid topography (W/m^2)
  Sdirt_n,              &! Incoming direct beam radiation corrected for subgrid topography (W/m^2)
  Sdift_n                ! Incoming diffuse radiation corrected for subgrid topography (W/m^2)


real*4 :: &
  hour,                &! Hour of day
  day_angle,           &! Day_angle in radians
  day_angle_decl,      &! Day angle in radians
  spring_equi,         &! Spring equinox time in days from the beginning of the year
  met_day_number_equi, &! Day number starting from spring equinox
  solar_decl,          &! Solar declination
  solar_eq_time,       &! Equation of time (hours)
  solar_solar_time,    &! True solar time
  solar_hr_angle,      &! Hour angle in radians
  sun_zenith            ! Solar zenith angle

real*4 :: &
  SWtopo                ! Incident SW radiation accounting for subgrid topography (W/m^2)

! Double precision parameter to avoid rounding errors for small or very big numbers:
real*8 :: &
  lambda_mu,           &! Radiation parameter for subgrid topography cf. Loewe and Helbig (2012)
  lambda_mu1,          &! Radiation parameter for subgrid topography cf. Loewe and Helbig (2012)
  lambda_mu2,          &! Radiation parameter for subgrid topography cf. Loewe and Helbig (2012)
  lambda_mu3,          &! Radiation parameter for subgrid topography cf. Loewe and Helbig (2012)
  lambda_elev,         &! Radiation parameter for subgrid topography cf. Loewe and Helbig (2012)
  slopemul              ! slope parameter

integer, intent(in) :: &
  i,j                   ! Model indices

! Store slope parameter in a double precision parameter
slopemul = slopemu(i,j)

!! Incident radiation WITH subgrid topography
!! Radiation terrain parameter cf. Loewe and Helbig (2012):
lambda_mu1 = sqrt(pi/(2.*slopemul**2))
lambda_mu2 = exp(1./(2.*slopemul**2))
lambda_mu3 = erfc(1./sqrt(2.*slopemul**2))
lambda_mu = lambda_mu1 * lambda_mu2 * lambda_mu3

!! To eliminate lambda_mu=NaN for mu's almost zero (lower than 5° slope angle here)  
!! resulting in Inf*Inf*0=NaN we set lamba_mu=1 for slope angles
!! lower equal 1 degree 
!! since cos(mu=0) equals one as well
if (slopemul <= 0.07) then
  lambda_mu = 1.
  fsky_terr(i,j) = 1.
end if
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Derive sun elevation angle of COSMO data grid
!! All equations are as in Helbig et al., 2010    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Calculate the day number of the year which is always one on 1 January and 365 on 31 December
!! Julian date of current date
jday  = day - 32075 + 1461*(year + 4800 + (month - 14) / 12) / 4 +  &
        &367 * (month - 2 - ((month - 14) / 12) * 12) / 12 - 3 * ((year&
        & + 4900 + (month - 14) / 12) / 100) / 4
!! Julian date at 1 January of current year:
jday1 = 1   - 32075 + 1461*(year + 4800 + (1     - 14) / 12) / 4 +  &
        &367 * (1   - 2 - ((1     - 14) / 12) * 12) / 12 - 3 * ((year&
        & + 4900 + (1     - 14) / 12) / 100) / 4

met_day_number = jday - jday1 + 1

!! Day_angle in radians (Iqbal (1983) p. 3)
day_angle = (pi / 180.) * ( (met_day_number - 1.) * 360. / 365.2425 )

!! Spring equinox time in days from the beginning of the year (Bourges (1985))
!! spring_equi = 78.801 + 0.2422 * (YYYY - 1969) - (int) ( 0.25 * (YYYY - 1969) );
!! starting from spring equinox at the year 2000 (20.3. 8:30(UTC+1-> 79.3542)
spring_equi = 79.3542 + 0.2422 * (year - 2000) - INT( 0.25 * (year - 2000) )

!! Actual day number transformed that the time in days starts from spring equinox (Bourges (1985))
met_day_number_equi = INT( met_day_number - 0.5 - spring_equi )

if ( met_day_number_equi >= epsilon(0.) ) then
  met_day_number_equi = met_day_number_equi
else
  !! because the day_number is lower than zero
  !! for days before spring equinox in the year using the Gregorian calendar with
  !! a mean year length of 365.2425 days
  met_day_number_equi = 365.2425 + met_day_number_equi
end if

!! Day angle in radians analog to Iqbal (1983); with met.day_number_equi as day_number
!! starting from spring equinox
!! using the Gregorian calendar mean number of days in a year = 365.2425
day_angle_decl = pi/180. * ( ((met_day_number_equi - 1.) * 360. / 365.2425) )

!! Solar declination after Bourges (1985) with a Fourier series approximation
!! and day number starting from 20.3. = 0 ...
!! and errors smaller than 0.02 degree
!! solar declination in degree
solar_decl  = 0.3723&
           &+ 23.2567 * sin (      day_angle_decl )&
           &- 0.758   * cos (      day_angle_decl )&
           &+ 0.1149  * sin ( 2. * day_angle_decl )&
           &+ 0.3656  * cos ( 2. * day_angle_decl )&
           &- 0.1712  * sin ( 3. * day_angle_decl )&
           &+ 0.0201  * cos ( 3. * day_angle_decl )

!! Transformation to radians
solar_decl = pi/180. * ( solar_decl )

!! Equation of time in hours (Spencer (1971), Iqbal (1983), p.11)
!! 229.18 = 24*60/2PI (conversion in minutes)
solar_eq_time = 229.18 * (0.000075&
                       &+ 0.001868 * cos( day_angle )&
                       &- 0.032077 * sin( day_angle )&
                       &- 0.014615 * cos( 2. * day_angle )&
                       &- 0.040849 * sin( 2. * day_angle )) / 60.

!!! Local apparent time or true SOLAR TIME in hours (Oke (1987), p.340 or Iqbal (1983), p.12):
!!! Meteo input always has to be in local standard time (winter in Davos: UTC+1) !!!!

if (lon(i,j) >= 0) then
   solar_solar_time = hour&
                       &+ (4. / 60.) * (lon(i,j)&
                       &-INT( lon(i,j) / 15. + 0.5 ) * 15.)&
                       &+ solar_eq_time
else
   !! For the western hemisphere
   solar_solar_time = hour&
                      &+ (4. / 60.) * (lon(i,j)&
                      &- INT( lon(i,j) / 15. - 0.5 ) * 15.)&
                      &+ solar_eq_time
end if

!! hour angle in radians (Oke (1987), p.340 or Iqbal (1983), p.15)
solar_hr_angle = pi/180. * ( 15. * (12. - solar_solar_time) )


!! Solar elevation in radians, whereby solar zenith is the solar zenith angle in radians Oke (1987), p.339
sun_zenith = acos( sin( solar_decl ) * sin( pi/180. *( lat(i,j) ) )&
             &+ cos( solar_decl ) * cos( pi/180. *( lat(i,j) ) ) *&
             &cos( solar_hr_angle ) )

sun_elev = pi / 2. - sun_zenith
sun_elev = sun_elev

!! ATTENTION: Because Sdir includes the sin(sun_elev) we can neglect that in lambda_elev 
lambda_elev = erf((tan(sun_elev) / (slopemul * 0.3498))&
              &**0.4980) * lambda_mu

if (sun_elev <= 0) then
  sun_elev = 0
  lambda_elev = 0
end if

!! ATTENTION: Actually we have to weight every albedo with the corresponding fraction and rad-flux
!! for now this doesn't have an impact because rad_snowfree=rad_snowcovered 
!! but with vegetation it would make a difference since rad_vegetation equals transmissivitaet*rad_snowfree
!! Incident SW radiation taking into account subgrid topographic influences:
SWtopo = (1. + alball_n * (1. - fsky_terr(i,j))) * (Sdir_n * lambda_elev + Sdif_n * fsky_terr(i,j))

! Incoming radiation per unit area of the model grid cell:
Sdirt_n = Sdir_n * lambda_elev / lambda_mu
Sdift_n = (Sdif_n * fsky_terr(i,j) + alball_n * (1. - fsky_terr(i,j)) * (Sdir_n * lambda_elev &
                  + Sdif_n * fsky_terr(i,j))) / lambda_mu

!! Net SW radiation per unit area of the model grid cell
SWtopo_out = (SWtopo * fsky_terr(i,j) * alball_n) / lambda_mu
SWtn = Sdir_n + Sdif_n - SWtopo_out

if (sun_elev <= epsilon(0.0)) then
  SWtn = 0.
end if
 
end subroutine SWRADTOPO
