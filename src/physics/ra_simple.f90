!>----------------------------------------------------------
!! Very simple radiation code modeled after the description in Finch and Best(2004?)
!! of Reiff et al. 1984 shortwave and Idso and Jackson (1969) longwave.
!!
!! Clearsky Shortwave radiation is calculated as a function of day of year and time of day.
!! Cloudy Shortwave is calculated as clearsky SW * f(cloud cover) [0.25-1]
!!
!! Cloud cover is calculated as in Xu and Randal (1996) as f(surface_RH, qc+qs+qr)
!!
!! Clearsky Longwave radiation is calculated as f(Tair)
!! Cloudy longwave is scaled up by a f(cloud cover) [1-1.2]
!!
!! The entry point to the code is ra_simple.
!!
!! <pre>
!! Call tree graph :
!! ra_simple->
!!  [cloudfrac->],
!!  [shortwave->],
!!  [longwave->]
!!
!! High level routine descriptions / purpose
!!   ra_simple           - loops over X,Y grid cells, calls cloudfrac, shortwave,longwave on columns
!!   cloudfrac           - calculates the cloud fraction following Xu and Randall (1996)
!!   shortwave           - calculates shortwave at the surface following Reiff et al (1984)
!!   longwave            - calculates longwave at the surface following Idso and Jackson (1969)
!!
!! Driver inputs: p,th,pii,rho,qv,qc,qr,qs,rain,snow,dt,dz,nx,ny,nz
!!   p   = pressure                      - 3D - input  - Pa     - (nx,nz,ny)
!!   th  = potential temperature         - 3D - in/out - K      - (nx,nz,ny)
!!   pii = inverse exner function        - 3D - input  - []     - (nx,nz,ny)
!!   rho = air density                   - 3D - input  - kg/m^3 - (nx,nz,ny)
!!   qv  = specific humidity             - 3D - input  - kg/kg  - (nx,nz,ny)
!!   qc  = cloud water content           - 3D - input  - kg/kg  - (nx,nz,ny)
!!   qr  = rain water content            - 3D - input  - kg/kg  - (nx,nz,ny)
!!   qs  = snow water content            - 3D - input  - kg/kg  - (nx,nz,ny)
!!   swdown = shortwave down at surface  - 2D - output - W/m^2  - (nx,ny)
!!   lwdown = longwave down at surface   - 2D - output - W/m^2  - (nx,ny)
!!   dt = time step                      - 0D - input  - seconds    - scalar
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module module_ra_simple
    use time_object,        only : Time_type
    use mod_atm_utilities,  only : relative_humidity
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use data_structures
    ! use time
    use iso_fortran_env, only: real128 !!MJ added
    
    implicit none

    real, allocatable, dimension(:,:) :: cos_lat_m,sin_lat_m
    integer :: nrad_layers
    real, parameter :: So = 1367.0  ! Solar "constant" W/m^2
    real, parameter :: qcmin = 1e-6 ! Minimum "cloud" water content to affect radiation
    real    :: tzone !! MJ adedd,tzone is UTC Offset and 1 here for centeral Erupe
contains

    subroutine ra_simple_init(domain,options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in)    :: options

        associate(ims => domain%grid%ims,   &
                  ime => domain%grid%ime,   &
                  jms => domain%grid%jms,   &
                  jme => domain%grid%jme)

        allocate(cos_lat_m(ims:ime, jms:jme))
        allocate(sin_lat_m(ims:ime, jms:jme))

        cos_lat_m = cos(domain%latitude%data_2d / 360.0 * 2*pi)
        sin_lat_m = sin(domain%latitude%data_2d / 360.0 * 2*pi)

        nrad_layers = 5

        end associate
        
        tzone=options%rad_options%tzone  !! MJ adedd,tzone is UTC Offset and 1 here for centeral Erupe
    end subroutine ra_simple_init


    function shortwave(day_frac, cloud_cover, solar_elevation, ims,ime, its,ite)
        ! compute shortwave down at the surface based on solar elevation, fractional day of the year, and cloud fraction
        ! based on Reiff et al. (1984)
        implicit none
        real             :: shortwave       (ims:ime)
        real, intent(in) :: day_frac        (ims:ime)
        real, intent(in) :: cloud_cover     (ims:ime)
        real, intent(in) :: solar_elevation (ims:ime)
        integer, intent(in) :: ims,ime, its,ite

        real :: sin_solar_elev(its:ite)

        sin_solar_elev = sin(solar_elevation(its:ite))

        shortwave(its:ite) = So * (1 + 0.035 * cos(day_frac(its:ite) * 2*pi)) * sin_solar_elev * (0.48 + 0.29 * sin_solar_elev)

        ! note it is cloud_cover**3.4 in Reiff, but this makes almost no difference and integer powers are fast so could use **3
        shortwave(its:ite) = shortwave(its:ite) * (1 - (0.75 * (cloud_cover(its:ite)**3.4)) )

    end function shortwave

    function longwave(T_air, cloud_cover, ims,ime, its,ite)
        ! compute longwave down at the surface based on air temperature and cloud fraction
        ! based on Idso and Jackson (1969)
        implicit none
        real                :: longwave(ims:ime)
        real,    intent(in) :: T_air(ims:ime), cloud_cover(ims:ime)
        integer, intent(in) :: ims,ime, its,ite
        real :: effective_emissivity(its:ite)

        effective_emissivity = 1 - 0.261 * exp((-7.77e-4) * (273.16-T_air(its:ite))**2)

        longwave(its:ite) = effective_emissivity * stefan_boltzmann * T_air(its:ite)**4

        longwave(its:ite) = min(longwave(its:ite) * (1 + 0.2 * cloud_cover(its:ite)), 600.0)

    end function longwave

    function cloudfrac(rh, qc, ims,ime, its,ite)
        ! Calculate the cloud fraction based on cloud water content (qc) and relative humidity
        ! based on equations from Xu and Randal (1996)
        implicit none
        real                :: cloudfrac(ims:ime)
        real,    intent(in) :: rh(ims:ime), qc(ims:ime)
        integer, intent(in) :: ims, ime, its, ite

        real :: temporary(its:ite)

        cloudfrac = 0

        temporary = ((1 - rh(its:ite)) * qc(its:ite))**0.25
        where(temporary > 1) temporary=1
        where(temporary < 0.0001) temporary=0.0001

        cloudfrac(its:ite) = qc(its:ite) - qcmin
        where(cloudfrac < 5e-8) cloudfrac = 5e-8

        cloudfrac(its:ite) = (rh(its:ite)**0.25) * (1 - exp((-2000*(cloudfrac(its:ite))) / temporary))

        where(cloudfrac < 0) cloudfrac = 0
        where(cloudfrac > 1) cloudfrac = 1

    end function

    !! MJ corrected, as calc_solar_elevation has largley understimated the zenith angle in Switzerland
    !! MJ added: this is Tobias Jonas (TJ) scheme based on swr function in metDataWizard/PROCESS_COSMO_DATA_1E2E.m and also https://github.com/Tobias-Jonas-SLF/HPEval
    !! MJ: note that this works everywhere and may be checked by https://gml.noaa.gov/grad/solcalc/index.html
    !! MJ: the only parameter needs to be given is https://gml.noaa.gov/grad/solcalc/index.html UTC Offset here referred to tzone=1 for centeral Erupe. HACK: this should be given by use in the namelist file
    !! MJ: Julian_day is a large value, we need to use the real128 format when applying TJ scheme in HICAR.
    function calc_solar_elevation_corr(date, lon, j, ims,ime, jms,jme, its,ite, day_frac, solar_azimuth)
        implicit none
        real                       :: calc_solar_elevation_corr(ims:ime)
        type(Time_type),intent(in) :: date
        real,           intent(in) :: lon(ims:ime, jms:jme)
        integer,        intent(in) :: j
        integer,        intent(in) :: ims, ime, jms, jme
        integer,        intent(in) :: its, ite
        real,           intent(out):: day_frac(ims:ime)
        real, optional, intent(inout):: solar_azimuth(ims:ime)
        
        integer :: i
        real, dimension(ims:ime) :: declination, day_of_year, hour_angle
        integer :: year, month, day, hour, minute, second
        real(real128) :: timeofday, julian_day, julian_century!, tzone
        real(real128) :: geom_mean_long_sun_deg, geom_mean_anom_sun_deg, eccent_earth_orbit
        real(real128) :: sun_eq_of_ctr, sun_true_long_deg, sun_app_long_deg
        real(real128):: mean_obliq_ecliptic_deg, obliq_corr_deg, sun_declin_deg, var_y, eq_of_time_minutes, true_solar_time_min
        real(real128) :: hour_angle_deg, solar_zenith_angle_deg, solar_elev_angle_deg
        real(real128) :: lat_hr, lon_hr
        real(real128) :: approx_atm_refrac_deg, solar_elev_corr_atm_ref_deg, solar_azimuth_angle

        !!
        do i = its, ite
            day_of_year(i) = date%day_of_year(lon=lon(i,j))

            ! hour angle is 0 at noon
            hour_angle(i) = 2*pi* mod(day_of_year(i)+0.5, 1.0)

            day_frac(i) = date%year_fraction(lon=lon(i,j))
        end do
    
        !!
        year=date%year
        month=date%month
        day=date%day
        hour=date%hour
        minute=date%minute
        second=date%second
        !tzone=1.
        !if (this_image()==2) write(*,*),"tzone ", tzone
        !year=2016
        !month=3
        !day=3
        !hour=12
        !minute=0
        !second=0
        !tzone=-5.

        !lon_hr=-71.05
        !lat_hr=42.35        
        !!
        calc_solar_elevation_corr = 0
        if(present(solar_azimuth)) solar_azimuth = 0


        !!       
        timeofday        = real(hour)/24.+real(minute)/60./24.+real(second)/3600./24.
        julian_day       = date%date_to_jd(year,month,day,hour,minute,second)-tzone/24.
        julian_century   = (julian_day - 2451545) / 36525
        !!
        geom_mean_long_sun_deg = mod(280.46646 + julian_century * (36000.76983 + julian_century * 0.0003032),360.)
        geom_mean_anom_sun_deg = 357.52911 + julian_century * (35999.05029 - 0.0001537 * julian_century)
        eccent_earth_orbit = 0.016708634 - julian_century * (0.000042037 + 0.0000001267 * julian_century)
        !!
        sun_eq_of_ctr = sin(pi/180 *(geom_mean_anom_sun_deg)) * (1.914602 - julian_century * (0.004817 + 0.000014 * julian_century)) + sin(pi/180 *(2  * geom_mean_anom_sun_deg)) * ( 0.019993 - 0.000101 * julian_century) + sin(pi/180 *(3 * geom_mean_anom_sun_deg)) * 0.000289
        sun_true_long_deg = sun_eq_of_ctr + geom_mean_long_sun_deg
        sun_app_long_deg = sun_true_long_deg - 0.00569 - 0.00478 * sin(pi/180 *(125.04 - 1934.136 * julian_century))
        !!
        mean_obliq_ecliptic_deg = 23 + (26 + ((21.448 - julian_century * (46.815 + julian_century * (0.00059 - julian_century * 0.001813)))) / 60) / 60
        obliq_corr_deg = mean_obliq_ecliptic_deg + 0.00256  * cos(pi/180 *(125.04 - 1934.136 * julian_century))
        sun_declin_deg = 180./pi*(asin(sin(pi/180 *(obliq_corr_deg)) * sin(pi/180 *(sun_app_long_deg))))
        var_y = tan(pi/180 *(obliq_corr_deg / 2)) * tan(pi/180 *(obliq_corr_deg / 2))
        eq_of_time_minutes = 4 * 180./pi*(var_y  * sin(2 * pi/180. *(geom_mean_long_sun_deg)) - 2 * eccent_earth_orbit * sin(pi/180. *(geom_mean_anom_sun_deg)) + 4 * eccent_earth_orbit * var_y * sin(pi/180. *(geom_mean_anom_sun_deg)) * cos(2  * pi/180. *(geom_mean_long_sun_deg)) - 0.5 * var_y * var_y * sin(4 * pi/180. *(geom_mean_long_sun_deg)) - 1.25 * eccent_earth_orbit * eccent_earth_orbit * sin(2 * pi/180. *(geom_mean_anom_sun_deg)))
        !!
        do i = its, ite           
            !!
            lon_hr=lon(i,j)
            lat_hr=180./pi*asin(sin_lat_m(i,j))                   
            true_solar_time_min = mod(timeofday * 1440 + eq_of_time_minutes + 4 * lon_hr - 60. * tzone,1440.);
            !!
            if (true_solar_time_min /4 < 0) then
                hour_angle_deg=true_solar_time_min /4 + 180
            elseif (true_solar_time_min /4 >= 0) then 
                hour_angle_deg=true_solar_time_min /4 - 180
            endif
            !!
            solar_zenith_angle_deg = 180./pi*(acos(sin(pi/180. *(lat_hr)) * sin(pi/180. *(sun_declin_deg)) + cos(pi/180. *(lat_hr)) * cos(pi/180. *(sun_declin_deg)) * cos(pi/180. *(hour_angle_deg))))
            solar_elev_angle_deg = 90 - solar_zenith_angle_deg;

            !! calculate atmospheric diffraction dependent on solar elevation angle
            if (solar_elev_angle_deg > 85) then
               approx_atm_refrac_deg=0. 
            elseif (solar_elev_angle_deg > 5 .and. solar_elev_angle_deg <= 85) then
                approx_atm_refrac_deg = (58.1 / tan(pi/180. *(solar_elev_angle_deg)) - 0.07 / (tan(pi/180. *(solar_elev_angle_deg)))**3. + 0.000086 / (tan(pi/180. *(solar_elev_angle_deg)))**5.) / 3600 
            elseif (solar_elev_angle_deg > -0.757 .and. solar_elev_angle_deg <= 5) then 
                approx_atm_refrac_deg = (1735 + solar_elev_angle_deg * (-518.2 + solar_elev_angle_deg * (103.4 + solar_elev_angle_deg * (-12.79 + solar_elev_angle_deg * 0.711)))) / 3600
            elseif (solar_elev_angle_deg <= -0.757) then 
                approx_atm_refrac_deg = (-20.772 / tan(pi/180. *(solar_elev_angle_deg))) / 3600
            endif                       
            solar_elev_corr_atm_ref_deg = solar_elev_angle_deg + approx_atm_refrac_deg
            
            !! calculate solar azimuth angle depending on hour angle
            if (hour_angle_deg > 0) then
                solar_azimuth_angle = mod(floor((180./pi*(acos(((sin(pi/180.*(lat_hr)) * cos(pi/180.*(solar_zenith_angle_deg))) - sin(pi/180.*(sun_declin_deg))) / (cos(pi/180.*(lat_hr)) * sin(pi/180.*(solar_zenith_angle_deg))))) + 180)*100000)/100000,360);
            elseif (hour_angle_deg <= 0) then
                solar_azimuth_angle = mod(floor((540 - 180./pi*(acos(((sin(pi/180.*(lat_hr)) * cos(pi/180.*(solar_zenith_angle_deg))) - sin(pi/180.*(sun_declin_deg))) / (cos(pi/180.*(lat_hr)) * sin(pi/180.*(solar_zenith_angle_deg))))))*100000)/100000,360);      
            endif                       
            
            calc_solar_elevation_corr(i)=solar_elev_corr_atm_ref_deg*pi/180.
            if(present(solar_azimuth)) solar_azimuth(i)=solar_azimuth_angle*pi/180.
        end do

        where(calc_solar_elevation_corr<0) calc_solar_elevation_corr=0

!       if (this_image()==2 .and. j==jms+1) then
!           write(*,*), trim(date%as_string()),solar_elev_corr_atm_ref_deg, solar_azimuth_angle
!       endif
!       if (this_image()==2 .and. j==jms+1) then
!           write(*,*), trim(date%as_string()), julian_day
!       endif
!       if (this_image()==2) write(*,*),"julian_day ", trim(date%as_string()), julian_day
!       if (this_image()==2) write(*,*),"julian_century ", julian_century
!       if (this_image()==2) write(*,*),"geom_mean_long_sun_deg ", geom_mean_long_sun_deg
!       if (this_image()==2) write(*,*),"eccent_earth_orbit ", eccent_earth_orbit
!       if (this_image()==2) write(*,*),"elev ", calc_solar_elevation_corr(its)
!       stop
    end function calc_solar_elevation_corr

    function calc_solar_elevation(date, lon, j, ims,ime, jms,jme, its,ite, day_frac)
        implicit none
        real                       :: calc_solar_elevation(ims:ime)
        type(Time_type),intent(in) :: date
        real,           intent(in) :: lon(ims:ime, jms:jme)
        integer,        intent(in) :: j
        integer,        intent(in) :: ims, ime, jms, jme
        integer,        intent(in) :: its, ite
        real,           intent(out):: day_frac(ims:ime)

        integer :: i
        real, dimension(ims:ime) :: declination, day_of_year, hour_angle

        calc_solar_elevation = 0

        do i = its, ite
            day_of_year(i) = date%day_of_year(lon=lon(i,j))

            ! hour angle is 0 at noon
            hour_angle(i) = 2*pi* mod(day_of_year(i)+0.5, 1.0)

            day_frac(i) = date%year_fraction(lon=lon(i,j))
        end do

        ! fast approximation see : http://en.wikipedia.org/wiki/Position_of_the_Sun
        declination = (-0.4091) * cos(2.0*pi/365.0*(day_of_year+10))

        calc_solar_elevation(its:ite) = sin_lat_m(its:ite,j) * sin(declination(its:ite)) + &
                               cos_lat_m(its:ite,j) * cos(declination(its:ite)) * cos(hour_angle(its:ite))

        ! due to float precision errors, it is possible to exceed (-1 - 1) in which case asin will break
        where(calc_solar_elevation < -1)
            calc_solar_elevation = -1
        elsewhere(calc_solar_elevation > 1)
            calc_solar_elevation = 1
        endwhere

        calc_solar_elevation = asin(calc_solar_elevation)

        ! if the sun is below the horizon just set elevation to 0
        where(calc_solar_elevation<0) calc_solar_elevation=0
    end function calc_solar_elevation


    !! MJ added: based on https://solarsena.com/solar-azimuth-angle-calculator-solar-panels/
    function calc_solar_azimuth(date, lon, j, ims,ime, jms,jme, its,ite, day_frac, solar_elevation)
        implicit none
        real                       :: calc_solar_azimuth(ims:ime)
        type(Time_type),intent(in) :: date
        real,           intent(in) :: lon(ims:ime, jms:jme)
        integer,        intent(in) :: j
        integer,        intent(in) :: ims, ime, jms, jme
        integer,        intent(in) :: its, ite
        real,           intent(out):: day_frac(ims:ime)
        real,           intent(in):: solar_elevation(ims:ime)

        integer :: i
        real, dimension(ims:ime) :: declination, day_of_year, hour_angle

        calc_solar_azimuth = 0

        do i = its, ite
            day_of_year(i) = date%day_of_year(lon=lon(i,j))

            ! hour angle is 0 at noon
            hour_angle(i) = 2*pi* mod(day_of_year(i)+0.5, 1.0)

            day_frac(i) = date%year_fraction(lon=lon(i,j))
        end do

        ! fast approximation see : http://en.wikipedia.org/wiki/Position_of_the_Sun
        declination = (-0.4091) * cos(2.0*pi/365.0*(day_of_year+10))

        calc_solar_azimuth(its:ite) = ( cos_lat_m(its:ite,j) * sin(declination(its:ite)) - &
                               sin_lat_m(its:ite,j) * cos(declination(its:ite)) * cos(hour_angle(its:ite)) )/(1.e-16+cos(solar_elevation(its:ite)))

        ! due to float precision errors, it is possible to exceed (-1 - 1) in which case asin will break
        where(calc_solar_azimuth < -1)
            calc_solar_azimuth = -1
        elsewhere(calc_solar_azimuth > 1)
            calc_solar_azimuth = 1
        endwhere

        ! partitioning the answer based on the hour angle:
        where(hour_angle > pi)
            calc_solar_azimuth = acos(calc_solar_azimuth)
        elsewhere(calc_solar_azimuth <= pi)
            calc_solar_azimuth = 2*pi - acos(calc_solar_azimuth)
        endwhere
        
    end function calc_solar_azimuth
    
    subroutine ra_simple(theta, pii, qv, qc, qs, qr, p, swdown, lwdown, cloud_cover, lat, lon, date, options, dt, &
                ims, ime, jms, jme, kms, kme, &
                its, ite, jts, jte, kts, kte, &
                F_runlw)
        implicit none
        real, intent(inout):: theta (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: pii   (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: qv    (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: qc    (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: qs    (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: qr    (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: p     (ims:ime, kms:kme, jms:jme)
        real, intent(out)  :: swdown     (ims:ime, jms:jme)
        real, intent(out)  :: lwdown     (ims:ime, jms:jme)
        real, intent(out)  :: cloud_cover(ims:ime, jms:jme)
        real, intent(in)   :: lat        (ims:ime, jms:jme)
        real, intent(in)   :: lon        (ims:ime, jms:jme)
        type(Time_type),    intent(in) :: date
        type(options_t), intent(in) :: options
        real,               intent(in) :: dt
        integer,            intent(in) :: ims, ime, jms, jme, kms, kme
        integer,            intent(in) :: its, ite, jts, jte, kts, kte
        logical,            intent(in), optional :: F_runlw

        logical :: runlw
        real :: coolingrate
        integer :: j, k
        real, allocatable, dimension(:) :: rh, T_air, solar_elevation, hydrometeors, day_frac


        runlw = .True.
        if (present(F_runlw)) runlw = F_runlw

        !$omp parallel private(j,k,rh,T_air,solar_elevation,hydrometeors,day_frac,coolingrate) &
        !$omp shared(theta,pii,qv,p,qc,qs,qr,date,lon,cloud_cover,swdown,lwdown)                      &
        !$omp firstprivate(runlw, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        allocate(rh             (ims:ime))
        allocate(T_air          (ims:ime))
        allocate(solar_elevation(ims:ime))
        allocate(hydrometeors   (ims:ime))
        allocate(day_frac       (ims:ime))

        ! 1.5K/day radiative cooling rate (300 = W/m^2 at 270K)
        coolingrate = 1.5 * (dt / 86400.0) * stefan_boltzmann / 300.0

        !$omp do
        do j = jts, jte

            T_air = 0
            rh = 0
            do k = kts, kts + nrad_layers - 1
                T_air = T_air + (theta(:,k,j)*pii(:,k,j))
                rh    = rh    + relative_humidity((theta(:,k,j)*pii(:,k,j)), qv(:,k,j), p(:,k,j))
            enddo
            T_air = T_air / nrad_layers
            rh    = rh    / nrad_layers
            where(rh > 1) rh = 1

            hydrometeors = qc(:,kts,j) + qs(:,kts,j) + qr(:,kts,j)
            do k = kts+1, kte
                hydrometeors = hydrometeors + qc(:,k,j) + qs(:,k,j) + qr(:,k,j)
            end do
            where(hydrometeors<0) hydrometeors = 0

            !solar_elevation  = calc_solar_elevation(date, lon, j, ims,ime, jms,jme, its,ite, day_frac) !! MJ: it does not work in Erupoe
            solar_elevation  = calc_solar_elevation_corr(date, lon, j, ims,ime, jms,jme, its,ite, day_frac=day_frac)

            cloud_cover(:,j) = cloudfrac(rh, hydrometeors, ims,ime, its,ite)
            swdown(:,j)      = shortwave(day_frac, cloud_cover(:,j), solar_elevation, ims,ime, its,ite)
            if (runlw) then
                lwdown(:,j)      = longwave(T_air, cloud_cover(:,j), ims,ime, its,ite)

                ! apply a simple radiative cooling to the atmosphere
                theta(its:ite, kts:kte, j) = theta(its:ite, kts:kte, j) - (((theta(its:ite, kts:kte, j) * pii(its:ite, kts:kte, j)) ** 4) * coolingrate)
            endif
        end do
        !$omp end do

        deallocate(rh,T_air, solar_elevation, hydrometeors, day_frac)
        !$omp end parallel

    end subroutine ra_simple
end module module_ra_simple
