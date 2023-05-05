!-----------------------------------------------------------------------
! Surface energy balance in open areas or zero-layer forest canopy model
!-----------------------------------------------------------------------
subroutine EBALSRF_SBG(Ds1,KH,KHa,KHv,KWg,KWv,ks1,SWsrf,SWveg,Ts1, &
                   Esrf,Eveg,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Rsrf,LWt)

use MODCONF, only: CANMOD

use CONSTANTS, only: &
  cp,                &! Specific heat capacity of air (J/K/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  Lv,                &! Latent heat of vapourisation (J/kg)
  Rair,              &! Gas constant for air (J/K/kg)
  Rwat,              &! Gas constant for water vapour (J/K/kg)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  em_snow,           &! Emissivity snow for Stefan-Boltzmann
  em_soil,           &! Emissivity soil for Stefan-Boltzmann 
  Tm                  ! Melting point (K)

use DRIVING, only: &
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Ta                  ! Air temperature (K)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  fthresh             ! Forest fraction required for forest tile to be considered

use PARAMMAPS, only: &
  trcn                ! Canopy transmissivity

use STATE_VARIABLES, only: &
  Sice,              &! Ice content of snow layers (kg/m^2)
  Tsrf,              &! Surface temperature (K)
  fsnow,             &! fractional snow cover
  Tveg                ! Vegetation temperature (K)

use LANDUSE, only : &
  dem,               &! Terrain elevation (m)
  forest,            &! Grid cell forest fraction
  fveg                ! Canopy cover fraction

implicit none

real, intent(in) :: &
  Ds1(Nx,Ny),        &! Surface layer thickness (m)
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny),        &! Eddy diffusivity for water from vegetation (m/s)
  ks1(Nx,Ny),        &! Surface layer thermal conductivity (W/m/K)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny),      &! Net SW radiation absorbed by vegetation (W/m^2)
  Ts1(Nx,Ny)          ! Surface layer temperature (K)

real, intent(out) :: &
  Esrf(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
  Eveg(Nx,Ny),       &! Moisture flux from vegetation (kg/m^2/s)
  G(Nx,Ny),          &! Heat flux into surface (W/m^2)
  H(Nx,Ny),          &! Sensible heat flux to the atmosphere (W/m^2)
  Hsrf(Nx,Ny),       &! Sensible heat flux from the surface (W/m^2)
  LE(Nx,Ny),         &! Latent heat flux to the atmosphere (W/m^2)
  LEsrf(Nx,Ny),      &! Latent heat flux from the surface (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny),       &! Net radiation (W/m^2)
  Rsrf(Nx,Ny),       &! Net radiation absorbed by the surface (W/m^2)
  LWt(Nx,Ny)          ! Incoming longwave radiation corrected for subgrid topography (W/m^2)

integer :: & 
  i,j                 ! Point counters

real :: &
  D,                 &! dQsat/dT (1/K)
  dE,                &! Change in surface moisture flux (kg/m^2/s)
  dG,                &! Change in surface heat flux (W/m^2)
  dH,                &! Change in sensible heat flux (W/m^2)
  dR,                &! Change in net radiation (W/m^2)
  dTs,               &! Change in surface skin temperatures (K)
  Lh,                &! Latent heat (J/kg)
  Qs,                &! Saturation humidity
  rho,               &! Air density (kg/m^3)
  em_w,              &! Snow cover fraction weighted surface emissivity
  LWsrf(Nx,Ny),      &! Net LW radiation absorbed by the surface (W/m^2)
  Tsall,             &! With fsnow weighted surface temperature
  Ssub                ! Mass of snow available for sublimation (kg/m^2)
  
real*8 :: &
  lambda_mu           ! Radiation parameter for subgrid topography cf. Loewe and Helbig (2012)

do j = 1, Ny
do i = 1, Nx

  if (isnan(dem(i,j))) goto 1 ! Exclude points outside of the domain

  if (CANMOD == 1 .and. forest(i,j) < fthresh) goto 1 ! Exclude points that have no forest

  if ((CANMOD == 1 .and. fveg(i,j) == 0) .or. CANMOD == 0) then  ! Surface energy balance in forests handled by subroutine EBALFOR
    Tveg(i,j) = Ta(i,j) 

    ! Saturation humidity and density of air
    call QSAT(Ps(i,j),Tsrf(i,j),Qs)
    Lh = Lv
    if (Tsrf(i,j) < Tm .or. Sice(1,i,j) > epsilon(Sice(1,i,j))) Lh = Ls
    D = Lh*Qs/(Rwat*Tsrf(i,j)**2)
    rho = Ps(i,j) / (Rair*Ta(i,j))

    ! Explicit fluxes
    Esrf(i,j) = rho*KWg(i,j)*(Qs - Qa(i,j))
    G(i,j) = 2*ks1(i,j)*(Tsrf(i,j) - Ts1(i,j))/Ds1(i,j)
    H(i,j) = cp*rho*KH(i,j)*(Tsrf(i,j) - Ta(i,j))
    LE(i,j) = Lh*Esrf(i,j)
    Melt(i,j) = 0
    ! Call Subgrid parameterization for LW radiation to compute LWtopo, net LW 
    ! There is currently no separate snow and bare soil surface skin temperature in FSM2. 
    ! As a quick fix we therefore set bare soil surface temperature to air temperature 
    ! in case of a snow cover
    if (fsnow(i,j) > epsilon(fsnow(i,j))) then
      Tsall = (((Tsrf(i,j)**4)*fsnow(i,j) &
              + (Ta(i,j)**4)*(1-fsnow(i,j)))**0.25)
       em_w = fsnow(i,j) * em_snow + (1-fsnow(i,j)) * em_soil
    else
      Tsall = Ta(i,j)
      em_w = em_soil  
    end if
    LWsrf(i,j) = 0
    ! Call Subgrid parameterization for LW radiation to compute LWtopo, netto LWtn
    call LWRADTOPO(LW(i,j),LWsrf(i,j),LWt(i,j),lambda_mu,Tsall,em_w,i,j,trcn(i,j),Tveg(i,j))

    Rnet(i,j) = SWsrf(i,j) + LWsrf(i,j)
    !Rnet(i,j) = SWsrf(i,j) + trcn(i,j)*LW(i,j) - sb*Tsrf(i,j)**4  &
    !                       + (1 - trcn(i,j))*sb*Tveg(i,j)**4

    ! Surface energy balance increments without melt
    dTs = (Rnet(i,j) - G(i,j) - H(i,j) - LE(i,j)) /  &
          (4*sb*Tsrf(i,j)**3 + 2*ks1(i,j)/Ds1(i,j) + rho*(cp*KH(i,j) + Lh*D*KWg(i,j)))
    dE = rho*KWg(i,j)*D*dTs
    dG = 2*ks1(i,j)*dTs/Ds1(i,j) 
    dH = cp*rho*KH(i,j)*dTs
    dR = -4*sb*Tsrf(i,j)**3*dTs

    ! Surface melting
    if (Tsrf(i,j) + dTs > Tm .and. Sice(1,i,j) > epsilon(Sice(1,i,j))) then
      Melt(i,j) = sum(Sice(:,i,j))/dt
      dTs = (Rnet(i,j) - G(i,j) - H(i,j) - LE(i,j) - Lf*Melt(i,j)) /  &
            (4*sb*Tsrf(i,j)**3 + 2*ks1(i,j)/Ds1(i,j) + rho*(cp*KH(i,j) + Ls*D*KWg(i,j)))
      dE = rho*KWg(i,j)*D*dTs
      dG = 2*ks1(i,j)*dTs/Ds1(i,j)
      dH = cp*rho*KH(i,j)*dTs
      dR = -4*sb*Tsrf(i,j)**3*dTs
      if (Tsrf(i,j) + dTs < Tm) then
        call QSAT(Ps(i,j),Tm,Qs)
        Esrf(i,j) = rho*KWg(i,j)*(Qs - Qa(i,j))  
        G(i,j) = 2*ks1(i,j)*(Tm - Ts1(i,j))/Ds1(i,j)
        H(i,j) = cp*rho*KH(i,j)*(Tm - Ta(i,j))
        LE(i,j) = Ls*Esrf(i,j)
        ! Call Subgrid parameterization for LW radiation to compute LWtopo, net LW 
        ! There is currently no separate snow and bare soil surface skin temperature in FSM2. 
        ! As a quick fix we therefore set bare soil surface temperature to air temperature 
        ! in case of a snow cover
        ! HERE: Snow surface is melting therefore Tsurf is set to Tm
        Tsall = (((Tm**4)*fsnow(i,j) + (Ta(i,j)**4)*(1-fsnow(i,j)))**0.25)
        ! Call Subgrid parameterization for LW radiation to compute LWtopo, netto LWtn
        call LWRADTOPO(LW(i,j),LWsrf(i,j),LWt(i,j),lambda_mu,Tsall,em_w,i,j,trcn(i,j),Tveg(i,j))
        Rnet(i,j) = SWsrf(i,j) + LWsrf(i,j)
        !Rnet(i,j) = SWsrf(i,j) + trcn(i,j)*LW(i,j) - sb*Tm**4  &
        !                       + (1 - trcn(i,j))*sb*Tveg(i,j)**4
        Melt(i,j) = (Rnet(i,j) - H(i,j) - LE(i,j) - G(i,j)) / Lf
        Melt(i,j) = max(Melt(i,j), 0.)
        dE = 0
        dG = 0
        dH = 0
        dR = 0
        dTs = Tm - Tsrf(i,j)
      end if

    end if

    ! Update surface temperature and fluxes
    Esrf(i,j) = Esrf(i,j) + dE
    G(i,j) = G(i,j) + dG
    H(i,j) = H(i,j) + dH
    LE(i,j) = Lh*Esrf(i,j)
    Rnet(i,j) = Rnet(i,j) + dR
    Tsrf(i,j) = Tsrf(i,j) + dTs

    ! Sublimation limited by amount of snow after melt
    Ssub = sum(Sice(:,i,j)) - Melt(i,j)*dt
    if (Ssub > epsilon(Ssub) .and. Esrf(i,j)*dt > Ssub) then
      Esrf(i,j) = Ssub / dt
      LE(i,j) = Ls*Esrf(i,j)
      H(i,j) = Rnet(i,j) - G(i,j) - LE(i,j) - Lf*Melt(i,j)
    end if
    Hsrf(i,j) = H(i,j)
    LEsrf(i,j) = LE(i,j)
    Rsrf(i,j) = Rnet(i,j)

    ! In case a pixel is forested we set LWin to LWt to account for topography effects:
    LW(i,j) = LWt(i,j)

    if (CANMOD == 0) then
      ! Add fluxes from canopy in zero-layer model
      Eveg(i,j) = 0
      if (fveg(i,j) > epsilon(fveg(i,j))) then
        Eveg(i,j) = - KWv(i,j)*Esrf(i,j) / (KHa(i,j) + KWv(i,j))
        H(i,j) = KHa(i,j)*H(i,j) / (KHa(i,j) + KHv(i,j))
        Lh = Ls
        if (Tveg(i,j) > Tm) Lh = Lv
        LE(i,j) = LE(i,j) + Lh*Eveg(i,j)
        Rnet(i,j) = Rnet(i,j) + SWveg(i,j) +  &
                    (1 - trcn(i,j))*(LW(i,j) + sb*Tsrf(i,j)**4 - 2*sb*Tveg(i,j)**4)
      end if
    end if
  end if

  1 continue 

end do
end do

end subroutine EBALSRF_SBG
