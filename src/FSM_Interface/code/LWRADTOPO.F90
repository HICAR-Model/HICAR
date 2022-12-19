!--------------------------------------------------------------------------
! Include subgrid topographical influences on surface longwave energy balance terms
! Author: Nora Helbig
!--------------------------------------------------------------------------
subroutine LWRADTOPO(LW_n,LWtn,LWt_n,lambda_mu,Tsall_n,em_w,i,j,trcn,Tveg)

use CONSTANTS, only : &
  pi,            &! Pi
  sb,            &! Stefan-Boltzmann constant (W/m^2/K^4)
  em_snow,       &! Emissivity snow for Stefan-Boltzmann
  em_soil         ! Emissivity soil for Stefan-Boltzmann 

use LANDUSE, only : &  
  fsky_terr,     &! Terrain sky view fraction
  slopemu         ! slope parameter 

implicit none

real, intent(in) :: &
  LW_n,           &! Incoming flat field longwave radiation (W/m^2)
  Tsall_n,        &! Snow cover fraction weighted surface skin temperature (K)
  em_w,           &! Snow cover fraction weighted surface emissivity
  trcn,           &! Canopy transmissivity
  Tveg             ! Vegetation temperature (K)

real, intent(out) :: &
  LWt_n,          &! Incoming longwave radiation corrected for subgrid topography (W/m^2)
  LWtn             ! Net LW radiation accounting for subgrid topography (W/m^2)

real*8, intent(out) :: &
  lambda_mu       ! Radiation parameter for subgrid topography cf. Loewe and Helbig (2012)
 
real*4 :: &
  LWtopo,        &! Outgoing LW radiation accounting for subgrid topography (W/m^2)
  LWter           ! Longwave terrain emission
 

! Double precision parameter to avoid rounding errors for small or very big numbers:
real*8 :: &
  lambda_mu1,    &! Radiation parameter for subgrid topography cf. Loewe and Helbig (2012)
  lambda_mu2,    &! Radiation parameter for subgrid topography cf. Loewe and Helbig (2012)
  lambda_mu3,    &! Radiation parameter for subgrid topography cf. Loewe and Helbig (2012)
  slopemul        ! slope parameter
 
integer, intent(in) :: &
  i,j             ! Model indices

slopemul = slopemu(i,j)

  
!! Equations are as in Helbig, 2009 and Loewe and Helbig, 2012    

!! Incident radiation WITH subgrid topography
!! Radiation terrain parameter cf. Loewe and Helbig (2012):
lambda_mu1 = sqrt(pi/(2.*slopemul**2))
lambda_mu2 = exp(1./(2.*slopemul**2))
lambda_mu3 = erfc(1./sqrt(2.*slopemul**2))
lambda_mu = lambda_mu1 * lambda_mu2 * lambda_mu3

!! To eliminate lambda_mu=NaN for mu almost 0 (lower 5° slope angle here)
!! resulting in Inf*Inf*0=NaN we set lamba_mu=1
!! since cos(mu=0) equals one as well
if (slopemul <= 0.07) then
  lambda_mu = 1.
  fsky_terr(i,j) = 1.
end if

!! Longwave terrain emission to include emissions from surrounding terrain
LWter = sb * em_w * Tsall_n**4

!! Longwave radiation (cf. Oliphant et al. 2003)
LWtopo = (LWter + (1-em_w) * fsky_terr(i,j) * LW_n) / (1 - (1 - fsky_terr(i,j)) * (1 - em_w)) 

!! Incoming longwave radiation per unit area of the model grid cell
LWt_n = LWtopo / lambda_mu

!! Outgoing LW radiation accounting for subgrid topography (W/m^2) per unit area of the model grid cell
LWtopo = LWtopo * fsky_terr(i,j) / lambda_mu

! Include canopy (new FSM version)
LWtn = trcn * LW_n + (1 - trcn) * sb * Tveg**4 - LWtopo

end subroutine LWRADTOPO
