!-----------------------------------------------------------------------
! Cumulate fluxes
!-----------------------------------------------------------------------
subroutine CUMULATE_interface(Roff,meltflux_out,Esrf,Gsoil,KH,LE,Melt,Rnet,H)
!MJ added-----------------------------------------------------------------
use MODULES_interface, only : Roff_,meltflux_out_,Esrf_,Gsoil_,KH_,LE_,Melt_,Rnet_,H_
!MJ added-----------------------------------------------------------------

use GRID, only: &
  Nx,Ny           ! Grid dimensions

implicit none

real, intent(in) :: &
  Roff(Nx,Ny),       &! Total runoff (kg/m^2)
  meltflux_out(Nx,Ny),  &! Runoff from snowmelt at base of snow (kg/m^2)
  Esrf(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
  H(Nx,Ny),          &! Sensible heat flux to the atmosphere (W/m^2)
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  LE(Nx,Ny),         &! Latent heat flux to the atmosphere (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny)         ! Net radiation (W/m^2)

  !call CUMULATE_SD_interface()

!! MJ added-----------------------------------------------------------------
  Esrf_= Esrf
  Gsoil_ = Gsoil
  H_ = H
  LE_ = LE
  Melt_ = Melt
  Rnet_ = Rnet
  Roff_ = Roff
  meltflux_out_ = meltflux_out
  KH_=KH

end subroutine CUMULATE_interface


subroutine CUMULATE_SD_interface()

!MJ added-----------------------------------------------------------------
use MODULES_interface, only : snowdepth_, SWE_, Sliq_out_
!MJ added-----------------------------------------------------------------

use GRID, only: &
  Nx,Ny           ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Sice,              &! Ice content of snow layers (kg/m^2)
  fsnow,             &! Snow cover fraction
  Sliq                ! Liquid content of snow layers (kg/m^2)

implicit none

real :: &
  Sliq_out(Nx,Ny),   &!
  snowdepth(Nx,Ny),  &! Snow depth (m)
  SWE(Nx,Ny)          ! Snow water equivalent (kg/m^2)

integer :: i,j,k

! BC just in case, these sums should be performed only until Nsnow.
do j = 1,Ny
  do i = 1,Nx
    Sliq_out(i,j) = sum(Sliq(:,i,j))
    snowdepth(i,j) = sum(Ds(:,i,j)) * fsnow(i,j)
    SWE(i,j) = sum(Sice(:,i,j)) + sum(Sliq(:,i,j))
  end do
end do

  snowdepth_=snowdepth
  SWE_=SWE
  Sliq_out_=Sliq_out

end subroutine CUMULATE_SD_interface

subroutine CUMULATE_SNOWTRAN3D_interface(dm_salt,dm_susp,dm_subl,dm_subgrid)

!MJ added-----------------------------------------------------------------
use MODULES_interface, only : dm_salt_, dm_susp_, dm_subl_
!MJ added-----------------------------------------------------------------

use GRID, only: &
  Nx,Ny           ! Grid dimensions

implicit none

real, intent(in) :: &
  dm_salt(Nx,Ny),   &! saltation flux per timestep (kg/m^2)
  dm_susp(Nx,Ny),   &! sublimation per timestep (kg/m^2)
  dm_subl(Nx,Ny),   &! blowing snow sublimation per timestep (kg/m^2)
  dm_subgrid(Nx,Ny)  ! subgrid redistribution per timestep (kg/m^2)

  !call CUMULATE_SD_interface()

  dm_salt_=dm_salt
  dm_susp_=dm_susp
  dm_subl_=dm_subl
  !dm_subgrid_=dm_subgrid

end subroutine CUMULATE_SNOWTRAN3D_interface

subroutine CUMULATE_SNOWSLIDE_interface(dm_slide)

!MJ added-----------------------------------------------------------------
use MODULES_interface, only : dm_slide_
!MJ added-----------------------------------------------------------------

use GRID, only: &
  Nx,Ny           ! Grid dimensions

implicit none

real, intent(in) :: &
  dm_slide(Nx,Ny)    ! avalanched redistribution per timestep (kg/m^2)

  !call CUMULATE_SD_interface()

  dm_slide_=dm_slide

end subroutine CUMULATE_SNOWSLIDE_interface
