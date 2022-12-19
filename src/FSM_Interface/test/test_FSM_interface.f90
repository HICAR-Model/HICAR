!-----------------------------------------------------------------------
! Factorial Snow Model
!
! Richard Essery
! School of GeoSciences
! University of Edinburgh
!-----------------------------------------------------------------------
program test_FSM_interface

use FSM_interface , only:  FSM_SETUP, FSM_DRIVE
use LANDUSE
use GRID

implicit none


integer :: & 
  i,j,               &! Point counters
  k,                 &! Level counter
  iresults_count
  
  
!call FSM_SETUP(10,10)
!print*, pmultf(1,1) 

!  do j = 1, Ny
!    do i = 1, Nx
!      print*,pmultf(i,j),i*i+j*j
!    end do
!  end do 

!call FSM_DRIVE()
!call FSM_DRIVE()
!call FSM_DRIVE()
!call FSM_DRIVE()

print*, pmultf(1,1) 

if (this_image()==1) write(*,*) "Reading restart data", NUM_IMAGES()
if (this_image()==2) write(*,*) "Reading restart data", NUM_IMAGES()

end program test_FSM_interface
