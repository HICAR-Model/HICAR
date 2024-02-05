!-----------------------------------------------------------------------
! Set parameters, initialize prognostic variables and write metadata
!-----------------------------------------------------------------------
subroutine SETUP_interface()

!MJ added-----------------------------------------------------------------
use FSM_interface, only: Nx_HICAR,Ny_HICAR,NNsmax_HICAR,lat_HICAR,lon_HICAR,terrain_HICAR,dx_HICAR,slope_HICAR,shd_HICAR
!MJ added-----------------------------------------------------------------

use MODCONF, only: CANMOD,DENSTY,ALBEDO,CANMOD,CONDCT,DENSTY,&
  EXCHNG,HYDROL,SNFRAC,RADSBG,ZOFFST,SNTRAN,SNSLID,SNOLAY,HISWET,CHECKS,HN_ON,FOR_HN
  
use MODOUTPUT, only: LIST_DIAG_RESULTS, LIST_STATE_RESULTS

use CONSTANTS

use DIAGNOSTICS, only : Nave

use DRIVING, only: &
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  RH,                &! Relative humidity (%)
  Rf,                &! Rainfall rate (kg/m2/s)
  Sf,                &! Snowfall rate (kg/m2/s)
  Sf24h,             &! Snowfall 24hr (kg/m2)
  Sdir,              &! Incoming direct beam radiation on flat,unobstructed surface (W/m2)
  Sdif,              &! Incoming diffuse radiation on flat,unobstructed (W/m2)
  Ta,                &! Air temperature (K)
  Tv,                &! Time-varying canopy transmissivity for dSWR (-)
  Ua,                &! Wind speed (m/s)
  Udir,              &! Wind direction (degrees, clockwise from N)
  zT,                &! Temperature measurement height (m)
  zU,                &! Wind speed measurement height (m)
  zRH                 ! Relative humidity measurement height (m)
  
use GRID

use PARAMETERS

use PARAMMAPS

use SOILPARAMS 

use STATE_VARIABLES

use LANDUSE

implicit none
 
integer :: & 
  i,j,               &! Point counters
  k,                 &! Level counter
  iresults_count
  
integer :: &
  NNsmax, NNsoil, NNx, NNy, NNt, nml_unit
  
integer :: &
  NALBEDO,NCANMOD,NCONDCT,NDENSTY,NEXCHNG,NHYDROL,NSNFRAC,NRADSBG,NZOFFST,&
  NSNTRAN,NSNSLID,NSNOLAY,NHISWET,NCHECKS
  
real :: &
  DDs_min,DDs_surflay,zzT,zzU,zzRH

real :: &
  hcon_min            ! Thermal conductivity of soil minerals (W/m/K)
  
real, allocatable :: &
  fsat(:),           &! Initial moisture content of soil layers as fractions of saturation
  Tprof(:)            ! Initial soil layer temperatures (K)

character(len=200) :: nlst_file

character(len=4), dimension(36) :: CLIST_DIAG_RESULTS

character(len=4), dimension(12) :: CLIST_STATE_RESULTS

logical :: LHN_ON  ! activate the HN model

logical :: LFOR_HN

logical :: lexist

!-1- !!!!!!!!!!!!!!!!!!!!  READ THE NAMELIST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
namelist  /nam_grid/    NNx,NNy,NNsmax,NNsoil,DDs_min,DDs_surflay
namelist  /nam_driving/ NNt,zzT,zzU,zzRH
namelist  /nam_modconf/ NALBEDO,NCANMOD,NCONDCT,NDENSTY,NEXCHNG,NHYDROL,NSNFRAC,NRADSBG,NZOFFST,&
                        NSNTRAN,NSNSLID,NSNOLAY,NHISWET,NCHECKS,LHN_ON,LFOR_HN
namelist /nam_results/ CLIST_DIAG_RESULTS, CLIST_STATE_RESULTS

! get namelist path from first command argument.
call getarg(1, nlst_file)
INQUIRE (FILE=nlst_file, EXIST=lexist)
if (lexist) then
  open(newunit=nml_unit, file = nlst_file)
else
  print*, '  namelist file: ', trim(nlst_file), ' does not exist'
  call exit(1)
endif

! Initialization of variables
NNsmax = 3
NNsoil = 4 ! Attention, 5 soil layers in JIM ... 
NNx = 1
NNy = 1
DDs_min = 0.02 ! Minimum possible snow layer thickness (m)
DDs_surflay = 0.5 ! Maximum thickness of surface fine snow layering (m)
read(nml_unit, nam_grid)
Nsoil = NNsoil

!!!!!!!!----> MJ: we receive it once in intializing for each image
Nx = Nx_HICAR !NNx
Ny = Ny_HICAR !NNy
Nsmax = NNsmax_HICAR

Ds_min = DDs_min
Ds_surflay = DDs_surflay

allocate(Dzsnow(Nsmax))
allocate(Dzsoil(Nsoil))
if (Nsmax == 3) then
  Dzsnow = (/0.1, 0.2, 0.4/)
else if (Nsmax == 6) then
  Dzsnow = (/0.05, 0.05, 0.1, 0.1, 0.2, 0.2/)
else
  print*, "Nsmax /= 3 or 6, Please specify your values for Dzsnow" 
  call exit()
endif
if (Nsoil == 4) then
  Dzsoil = (/0.1, 0.2, 0.4, 0.8/) ! Attention, 0.05, 0.1, 0.1, 0.5, 1.25 in JIM...
else
  print*, "Nsoil /= 4, Please specify your values for Dzsoil" 
  call exit()
endif

! Driving data
NNt = 1
zzT = 10
zzU = 10
zzRH = 10
read(nml_unit,nam_driving)
Nt = NNt
zT = zzT
zU = zzU
zRH = zzRH
dt = 1.0 / Nt

! Model configuration
! -1 for mandatory NLST arguments.
NALBEDO = -1
NCANMOD = -1
NCONDCT = -1
NDENSTY = -1
NEXCHNG = -1
NHYDROL = -1
NSNFRAC = -1
NRADSBG = -1
NZOFFST = -1
NSNTRAN = -1
NSNSLID = -1
NSNOLAY = -1
NHISWET = -1
NCHECKS = -1
LHN_ON = .FALSE.
LFOR_HN = .FALSE.
read(nml_unit, nam_modconf)
ALBEDO = NALBEDO
CANMOD = NCANMOD
CONDCT = NCONDCT
DENSTY = NDENSTY
EXCHNG = NEXCHNG
HYDROL = NHYDROL
SNFRAC = NSNFRAC
RADSBG = NRADSBG
ZOFFST = NZOFFST
SNTRAN = NSNTRAN
SNSLID = NSNSLID
SNOLAY = NSNOLAY
HISWET = NHISWET
CHECKS = NCHECKS
HN_ON = LHN_ON
FOR_HN = LFOR_HN
if (ALBEDO==-1 .or. CANMOD==-1 .or. CONDCT==-1 .or. DENSTY==-1 .or. EXCHNG==-1 &
    .or. HYDROL==-1 .or. SNFRAC==-1 .or. RADSBG==-1 .or. ZOFFST==-1 &
    .or. SNTRAN==-1 .or. SNSLID==-1 .or. SNOLAY==-1 .or. HISWET==-1 .or. CHECKS ==-1) then
  print*, 'model configuration error:\n please specify all the fields of MODCONF in the namelist (&nam_modconf)'
  call exit(1)
endif

! List of output variables
CLIST_DIAG_RESULTS(:)  = '    ' ! BC 13-char length
CLIST_STATE_RESULTS(:) = '    ' ! BC 13-char length
read(nml_unit, nam_results)
iresults_count = count(CLIST_DIAG_RESULTS /= '             ') ! BC 13-char length
LIST_DIAG_RESULTS(1:iresults_count) = CLIST_DIAG_RESULTS(1:iresults_count)
iresults_count = count(CLIST_STATE_RESULTS /= '             ') ! BC 13-char length
LIST_STATE_RESULTS(1:iresults_count) = CLIST_STATE_RESULTS(1:iresults_count)

! Outputs
Nave = 1 !24 Set to one, 24h averages created in OSHD Matlab wrapper

close(nml_unit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-2- !!!!!!!!!!!!!!!!!!!!    OPEN THE FILES   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open Files for OSHD states input and output and results output
!!!!!!!!----> call OPEN_FILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-3- !!!!!!!!!!!!!!!!!!!!  ALLOCATE VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(LW(Nx,Ny))
allocate(Ps(Nx,Ny))
allocate(Qa(Nx,Ny))
allocate(Rf(Nx,Ny))
allocate(Sf(Nx,Ny))
allocate(Sdir(Nx,Ny))
allocate(Sdif(Nx,Ny))
allocate(Ta(Nx,Ny))
allocate(Ua(Nx,Ny))
if (SNTRAN == 1) then
  allocate(Udir(Nx,Ny))
endif
allocate(RH(Nx,Ny))
allocate(Sf24h(Nx,Ny))
allocate(Tv(Nx,Ny))

! use Tv dummy in case of open simulations
if (CANMOD == 0) then
  Tv(:,:) = 1
endif

! Defaults for numerical solution parameters
Nitr = 4

! Defaults for canopy parameters
avg0 = 0.1
avgs = 0.4
cden = 0.004
cvai = 6.6
cveg = 20
Gcn1 = 0.5
Gcn2 = 0
gsnf = 0
kdif = 0.5
kveg = 1
rchd = 0.67
rchz = 0.2          
tcnc = 240
tcnm = 48

! Defaults for snow parameters
a_eta = 0.1
asmx = 0.8
asmn = 0.5
b_eta = 0.023
bstb = 5
bthr = 2
c_eta = 250
eta0 = 3.7e7
eta1 = 7.62237e6
hfsn = 0.1
kfix = 0.24
rho0 = 300
rhob = 6
rhoc = 26
rhof = 109
rhos_min = 50.0
rhos_max = 750.0
rcld = 300
rgr0 = 5e-5
rmlt = 500
Salb = 10
snda = 2.8e-6
Talb = -2
tcld = 1000
tmlt = 100
trho = 200
Wirr = 0.03
z0sn = 0.002

! Defaults for ground surface parameters
bstb = 5
gsat = 0.01

! Defaults for additional parameters required for forest snow process parametrization
adfs = 3
adfl = 2
fsar = 0.2
psf  = 1.1
psr  = 0.2
wcan = 2.5
zsub = 2
zgf = 5
zgr = 5
khcf = 3

! Default for tile parameters 
fthresh = 0.1

if (DENSTY == 0) then
  rhof = rho0
endif

! Surface data from defaults or namelists 
! Surface properties 
allocate(alb0(Nx,Ny))
allocate(fcly(Nx,Ny))
allocate(fsnd(Nx,Ny))
allocate(z0sf(Nx,Ny))
allocate(vegsnowd_xy(Nx,Ny))
alb0(:,:) = 0.2
fcly(:,:) = 0.3
fsnd(:,:) = 0.6
z0sf(:,:) = 0.2
vegsnowd_xy(:,:) = 0.1

! Canopy parameters
allocate(canh(Nx,Ny))
allocate(fsky(Nx,Ny))
allocate(fveg(Nx,Ny)) 
allocate(fves(Nx,Ny))
allocate(hcan(Nx,Ny))
allocate(lai(Nx,Ny))
allocate(scap(Nx,Ny))
allocate(trcn(Nx,Ny))
allocate(VAI(Nx,Ny))
allocate(vfhp(Nx,Ny))
canh(:,:) = undef
fsky(:,:) = undef
fveg(:,:) = undef
hcan(:,:) = undef
scap(:,:) = undef
trcn(:,:) = undef
VAI(:,:)  = undef

!Terrain properties
allocate(slopemu(Nx,Ny))
allocate(xi(Nx,Ny))
allocate(Ld(Nx,Ny))
allocate(lat(Nx,Ny))
allocate(lon(Nx,Ny))
allocate(dem(Nx,Ny))
allocate(slope(Nx,Ny))
allocate(Shd(Nx,Ny))
allocate(forest(Nx,Ny))
slopemu(:,:) = undef
xi(:,:) = undef
Ld(:,:) = undef
lat(:,:) = undef
lon(:,:) = undef
dem(:,:) = undef
slope(:,:) = undef
Shd(:,:) = undef
forest(:,:) = undef

! Derived soil parameters
allocate(b(Nx,Ny))
allocate(hcap_soil(Nx,Ny))
allocate(hcon_soil(Nx,Ny))
allocate(sathh(Nx,Ny))
allocate(Vsat(Nx,Ny))
allocate(Vcrit(Nx,Ny))
do j = 1, Ny
  do i = 1, Nx
    if (fcly(i,j) + fsnd(i,j) > 1) then
      fcly(i,j) = 1 - fsnd(i,j)
    endif
    b(i,j) = 3.1 + 15.7*fcly(i,j) - 0.3*fsnd(i,j)
    hcap_soil(i,j) = (2.128*fcly(i,j) + 2.385*fsnd(i,j))*1e6 / (fcly(i,j) + fsnd(i,j))
    sathh(i,j) = 10**(0.17 - 0.63*fcly(i,j) - 1.58*fsnd(i,j))
    Vsat(i,j) = 0.505 - 0.037*fcly(i,j) - 0.142*fsnd(i,j)
    Vcrit(i,j) = Vsat(i,j)*(sathh(i,j)/3.364)**(1/b(i,j))
    hcon_min = (hcon_clay**fcly(i,j)) * (hcon_sand**(1 - fcly(i,j)))
    hcon_soil(i,j) = (hcon_air**Vsat(i,j)) * (hcon_min**(1 - Vsat(i,j)))
  end do
end do

! Convert time scales from hours to seconds
dt = 3600*dt
tcnc = 3600*tcnc
tcnm = 3600*tcnm
tcld = 3600*tcld
tmlt = 3600*tmlt
trho = 3600*trho

! Allocate state variables
allocate(albs(Nx,Ny))
allocate(Ds(Nsmax,Nx,Ny))
allocate(Nsnow(Nx,Ny))
allocate(Qcan(Nx,Ny))
allocate(rgrn(Nsmax,Nx,Ny))
allocate(Sice(Nsmax,Nx,Ny))
allocate(Sliq(Nsmax,Nx,Ny))
allocate(Sveg(Nx,Ny))
allocate(Tcan(Nx,Ny))
allocate(theta(Nsoil,Nx,Ny))
allocate(Tsnow(Nsmax,Nx,Ny))
allocate(Tsoil(Nsoil,Nx,Ny))
allocate(Tsrf(Nx,Ny))
allocate(fsnow(Nx,Ny))
allocate(Tveg(Nx,Ny))
allocate(snowdepthmin(Nx,Ny))
allocate(snowdepthmax(Nx,Ny))
allocate(snowdepthhist(Nx,Ny,14))
allocate(swemin(Nx,Ny))
allocate(swemax(Nx,Ny))
allocate(swehist(Nx,Ny,14))
allocate(histowet(Nsmax,Nx,Ny))
allocate(fsky_terr(Nx,Ny))
allocate(dm_tot_subl(Nx,Ny))
allocate(dm_tot_trans(Nx,Ny))
allocate(dm_tot_slide(Nx,Ny))
allocate(index_grid_dem_sorted(Nx*Ny,2))

! Default initialization of state variables 
firstit      = 1;
albs(:,:)    = undef
Ds(:,:,:)    = undef
fsnow(:,:)   = undef
Nsnow(:,:)   = iundef
Qcan(:,:)    = undef
rgrn(:,:,:)  = undef !*GM watch out: rgrn currently not tracked
Sice(:,:,:)  = undef
Sliq(:,:,:)  = undef
Sveg(:,:)    = undef
Tcan(:,:)    = undef
Tsnow(:,:,:) = undef
Tsoil(:,:,:) = undef
Tveg(:,:)    = undef
snowdepthmin(:,:) = undef
snowdepthmax(:,:) = undef
snowdepthhist(:,:,:) = undef
swemin(:,:) = undef
swemax(:,:) = undef
swehist(:,:,:) = undef
histowet(:,:,:) = undef
fsky_terr(:,:) = undef
dm_tot_subl(:,:) = undef
dm_tot_trans(:,:) = undef
dm_tot_slide(:,:) = undef
index_grid_dem_sorted(:,:) = iundef

! Initial soil profiles from namelist
allocate(fsat(Nsoil))
allocate(Tprof(Nsoil))
fsat(:)  = 0.5
Tprof(:) = 285
do k = 1, Nsoil
  theta(k,:,:) = fsat(k)*Vsat(:,:)
  Tsoil(k,:,:) = Tprof(k)
end do
Tsrf(:,:) = Tsoil(1,:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-4- !!!!!!!!!!!!!!!!!!!! READ DRIVING/STATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! states relevant to both open and forest simulation
albs(:,:)      =0.85 !!!!!!!!---->read(1101) albs
Ds(:,:,:)      =0.0  !!!!!!!!---->read(1102) Ds
fsnow(:,:)     =0.0  !!!!!!!!---->read(1103) fsnow
Nsnow(:,:)     =0.0  !!!!!!!!---->read(1104) Nsnow
Sice(:,:,:)    =0.0  !!!!!!!!---->read(1106) Sice
Sliq(:,:,:)    =0.0  !!!!!!!!---->read(1107) Sliq
!!!!!!!!---->read(1116) Tsrf
Tsnow(:,:,:)   =273.15	!!!!!!!!---->read(1119) Tsnow
!!!!!!!!---->read(1120) Tsoil
fsky_terr(:,:) =1.0  !!!!!!!!---->read(1123) fsky_terr
lat(:,:)=lat_HICAR(:,:) !!!!!!!!---->read(1127) lat
lon(:,:)=lon_HICAR(:,:) !!!!!!!!---->read(1128) lon
dem(:,:)=terrain_HICAR(:,:) !!!!!!!!---->read(1129) dem

if (SNFRAC == 0 .or. SNFRAC == 2) then
  snowdepthmax(:,:) =0.0 !!!!!!!!---->read(1110) snowdepthmax
endif

if (SNFRAC == 0 ) then
  ! states specific to open runs
  snowdepthmin(:,:)    =0.0  !!!!!!!!---->read(1109) snowdepthmin
  snowdepthhist(:,:,:) =0.0  !!!!!!!!---->read(1111) snowdepthhist
  swemin(:,:)          =0.0  !!!!!!!!---->read(1113) swemin
  swemax(:,:)          =0.0  !!!!!!!!---->read(1114) swemax
  swehist(:,:,:)       =0.0  !!!!!!!!---->read(1115) swehist
  slopemu(:,:)         =1.0  !!!!!!!!---->slopemu(:,:)=slopemu_(:,:) !read(1124) slopemu
  xi(:,:)              =1.0  !!!!!!!!---->read(1125) xi
  Ld(:,:)              =dx_HICAR !!!!!!!!---->read(1126) Ld
endif

if (CANMOD == 0) then
  ! canopy properties (no canopy)
  VAI(:,:)  = 0
  hcan(:,:) = 0
  fsky(:,:) = 1
  trcn(:,:) = exp(-kdif*VAI(:,:))
  fveg(:,:) = 1 - exp(-kveg*VAI(:,:))
  fves(:,:) = 1 - exp(-kveg*VAI(:,:))
else ! CANMOD == 1
  ! states specific to forest runs
  Qcan(:,:)    =0.0   !!!!!!!!---->read(1130) Qcan
  Sveg(:,:)    =0.0   !!!!!!!!---->read(1131) Sveg
  Tcan(:,:)    =285   !!!!!!!!---->read(1132) Tcan
  Tveg(:,:)    =285   !!!!!!!!---->read(1133) Tveg
  fveg(:,:)    =0.0   !!!!!!!!---->read(1134) fveg
  hcan(:,:)    =0.0   !!!!!!!!---->read(1135) hcan
  lai(:,:)     =0.0   !!!!!!!!---->read(1136) lai
  vfhp(:,:)    =1.0   !!!!!!!!---->read(1137) vfhp

  ! derived canopy properties 
  VAI(:,:) = lai(:,:) 
  fves(:,:) = fveg(:,:)  !!! *GM to be revised 
  do j = 1, Ny
    do i = 1, Nx
      trcn(i,j) = 1-VAI(i,j)/5
      fsky(i,j) = vfhp(i,j)/trcn(i,j)
      if ( fsky(i,j) > 1 ) trcn(i,j) = vfhp(i,j)
      if ( fsky(i,j) > 1 ) fsky(i,j) = 1
    end do
  end do 
endif

! derived canopy parameters
canh(:,:) = 12500*VAI(:,:)
scap(:,:) = cvai*VAI(:,:)

if (SNTRAN == 1) then
  ! states specific to SNOWTRAN3D
  !!!!!! read(1140) dm_tot_subl
  !!!!!! read(1141) dm_tot_trans
  Ld(:,:)              =dx_HICAR !!!!!!!!---->read(1126) Ld
endif

if (SNSLID == 1) then
  ! states specific to SnowSlide
  slope(:,:)=slope_HICAR(:,:)
  Shd(:,:)=shd_HICAR(:,:)
  !!!!!! read(1142) dm_tot_slide
  !!!!!! read(1143) index_grid_dem_sorted
endif

end subroutine SETUP_interface
