!-----------------------------------------------------------------------
! Update soil temperatures
!-----------------------------------------------------------------------
subroutine SOIL(csoil,Gsoil,ksoil)

use MODCONF, only: CANMOD

use DRIVING, only: &
  dt                  ! Timestep (s)

use GRID, only: &
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  fthresh             ! Forest fraction required for forest tile to be considered

use STATE_VARIABLES, only: &
  Tsoil               ! Soil layer temperatures (K)

use LANDUSE, only: &
  dem,               &! Terrain elevation (m)
  forest              ! Grid cell forest fraction
  
implicit none

real, intent(in) :: &
  csoil(Nsoil,Nx,Ny),&! Areal heat capacity of soil (J/K/m^2)
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
  ksoil(Nsoil,Nx,Ny)  ! Thermal conductivity of soil (W/m/K)

integer :: &
  i,j,               &! Point counters
  k                   ! Level counter

real :: &
  a(Nsoil),          &! Below-diagonal matrix elements
  b(Nsoil),          &! Diagonal matrix elements
  c(Nsoil),          &! Above-diagonal matrix elements
  dTs(Nsoil),        &! Temperature increments (k)
  Gs(Nsoil),         &! Thermal conductivity between layers (W/m^2/k)
  rhs(Nsoil)          ! Matrix equation rhs

do j = 1, Ny
do i = 1, Nx

  if (isnan(dem(i,j))) goto 1 ! Exclude points outside of the domain

  if (CANMOD == 1 .and. forest(i,j) < fthresh) goto 1 ! Exclude points that have no forest

  do k = 1, Nsoil - 1
    Gs(k) = 2 / (Dzsoil(k)/ksoil(k,i,j) + Dzsoil(k+1)/ksoil(k+1,i,j))
  end do
  a(1) = 0
  b(1) = csoil(1,i,j) + Gs(1)*dt
  c(1) = - Gs(1)*dt
  rhs(1) = (Gsoil(i,j) - Gs(1)*(Tsoil(1,i,j) - Tsoil(2,i,j)))*dt
  do k = 2, Nsoil - 1
    a(k) = c(k-1)
    b(k) = csoil(k,i,j) + (Gs(k-1) + Gs(k))*dt
    c(k) = - Gs(k)*dt
    rhs(k) = Gs(k-1)*(Tsoil(k-1,i,j) - Tsoil(k,i,j))*dt  &
             + Gs(k)*(Tsoil(k+1,i,j) - Tsoil(k,i,j))*dt 
  end do
  k = Nsoil
  Gs(k) = ksoil(k,i,j)/Dzsoil(k)
  a(k) = c(k-1)
  b(k) = csoil(k,i,j) + (Gs(k-1) + Gs(k))*dt
  c(k) = 0
  rhs(k) = Gs(k-1)*(Tsoil(k-1,i,j) - Tsoil(k,i,j))*dt
  call TRIDIAG(Nsoil,Nsoil,a,b,c,rhs,dTs)
  do k = 1, Nsoil
    Tsoil(k,i,j) = Tsoil(k,i,j) + dTs(k)
  end do

  1 continue
  
end do
end do

end subroutine SOIL
 