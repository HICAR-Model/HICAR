!-----------------------------------------------------------------------
! Perform a consistency check of snowpack state variables
!-----------------------------------------------------------------------
subroutine CHECK_SNOWPACK

use MODCONF, only: CHECKS

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use DRIVING, only: &
  hour                ! Hour of day

use PARAMETERS, only : &
  rhos_min,          &! Minimum snow density (kg/m^3)
  rhos_max            ! Maximum snow density (kg/m^3)

use GRID, only: &
  Nsmax,             &! Maximum number of snow layers
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Ds,                &! Snow layer thicknesses (m)
  histowet,          &! Historical variable for past wetting of a layer (0-1)
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Tsnow               ! Snow layer temperatures (K)

implicit none

integer :: &
  i,j,               &! Point counters
  k                   ! Snow layer counters

logical :: &
  critical,          &! Critical warning (0-1)
  warning             ! Warning (0-1)

real :: &
  rhos                ! SNow layer density (kg/m3)

character(len=80) :: &
  report              ! Warning text

do j = 1, Ny
  do i = 1, Nx

    warning = .false.
    critical = .false.

    if (Nsnow(i,j) < 0 &
        .or. Nsnow(i,j) > Nsmax) then

      warning = .true.
      critical = .true.
      
      report = 'Nsnow invalid'

    else if (size(Ds,1) > Nsmax &
             .or. size(Sice,1) > Nsmax &
             .or. size(Sliq,1) > Nsmax &
             .or. size(Tsnow,1) > Nsmax &
             .or. size(histowet,1) > Nsmax) then

      warning = .true.
      critical = .true.
      report = 'Invalid layer variable dimension'

    else if ((0.0 - fsnow(i,j)) > epsilon(fsnow) &
             .or. (fsnow(i,j) - 1.0) > epsilon(fsnow) &
             .or. isnan(fsnow(i,j))) then

      warning = .true.
      critical = .true.
      report = 'fsnow invalid'

    else if (fsnow(i,j) > epsilon(fsnow) &
             .and. (sum(Ds(:,i,j)) < epsilon(Ds) &
                    .or. sum(Sice(:,i,j)+Sliq(:,i,j)) < epsilon(Sice))) then

      warning = .true.
      critical = .true.
      report = 'fsnow inconsistent with layer variables'

    else if (fsnow(i,j) < epsilon(fsnow) &
             .and. (sum(Ds(:,i,j)) > epsilon(Ds) &
                    .or. sum(Sice(:,i,j)+Sliq(:,i,j)) > epsilon(Sice))) then

      warning = .true.
      critical = .true.
      report = 'fsnow inconsistent with layer variables'

    else

      do k = 1, Nsmax

        if ((0.0 - Ds(k,i,j)) > epsilon(Ds)  &
            .or. (0.0 - Sice(k,i,j)) > epsilon(Sice) &
            .or. (0.0 - Sliq(k,i,j)) > epsilon(Sliq) &
            .or. (Tsnow(k,i,j)-Tm) > epsilon(Tsnow) &
            .or. (0.0 - histowet(k,i,j)) > epsilon(histowet) &
            .or. (histowet(k,i,j) - 1.0) > epsilon(histowet) &
            .or. isnan(Ds(k,i,j)) &
            .or. isnan(Sice(k,i,j)) &
            .or. isnan(Sliq(k,i,j)) &
            .or. isnan(Tsnow(k,i,j)) &
            .or. isnan(histowet(k,i,j))) then

          warning = .true.
          critical = .true.
          report = 'Invalid layer variables'

        else if (Ds(k,i,j) > epsilon(Ds) &
                 .and. (Sice(k,i,j) + Sliq(k,i,j)) < epsilon(Sice)) then

          warning = .true.
          critical = .true.
          report = 'Inconsistent layer variables'

        else if (Ds(k,i,j) < epsilon(Ds) &
                 .and. (Sice(k,i,j) + Sliq(k,i,j)) > epsilon(Sice)) then

          warning = .true.
          critical = .true.
          report = 'Inconsistent layer variables'

        else if (Ds(k,i,j) > epsilon(Ds) &
                 .and. fsnow(i,j) > epsilon(fsnow)) then

          rhos = (Sice(k,i,j) + Sliq(k,i,j)) / Ds(k,i,j) / fsnow(i,j)

          if (isnan(rhos) &
              .or. rhos < rhos_min - 0.5 &
              .or. rhos > rhos_max + 0.5) then

            warning = .true.
            critical = .true.
            report = 'Invalid layer density'

          end if

        else if (k > Nsnow(i,j) &
                 .and. (Ds(k,i,j) > epsilon(Ds) &
                        .or. Sice(k,i,j) > epsilon(Sice) &
                        .or. Sliq(k,i,j) > epsilon(Sliq))) then

          warning = .true.
          critical = .true.
          report = 'Inconsistent layer variables'

        else if (Tsnow(k,i,j) < Tm - 100) then
        
          warning = .true.
          critical = .false.
          report = 'Snow temperature extremely low'
          
        end if

      end do

    end if

  if (warning) then
  
    write(*,*) '--CHECK_SNOWPACK'
    write(*,*) 'i',i,'j',j
    write(*,*) report
    write(*,*) 'hour', int(hour)
    write(*,*) 'Nsnow', Nsnow(i,j)
    write(*,*) 'fsnow', fsnow(i,j)
    write(*,*) 'Ds', Ds(:,i,j)
    write(*,*) 'Sice', Sice(:,i,j)
    write(*,*) 'Sliq', Sliq(:,i,j)
    write(*,*) 'Tsnow', Tsnow(:,i,j)
    write(*,*) 'histowet', histowet(:,i,j)
    write(*,*) 'rhos', (Sice(1:Nsnow(i,j),i,j)+Sliq(1:Nsnow(i,j),i,j))/ Ds(1:Nsnow(i,j),i,j) / fsnow(i,j)
    
    if (CHECKS == 2 .and. critical) then
      error stop
    end if
    
  end if  

  end do
end do

end subroutine CHECK_SNOWPACK
