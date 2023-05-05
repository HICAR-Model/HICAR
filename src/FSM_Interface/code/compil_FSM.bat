::-------------------------------------------------------------------------------------------
:: Flexible Snow Model DOS compilation script
::
:: Richard Essery
:: School of GeoSciences
:: University of Edinburgh
::
:: Edited by OSHD SLF (LQ, GM, BC) - 20210716
::-------------------------------------------------------------------------------------------
echo off
if %ComputerName%==STACHLERKOPF echo compilation on stachlerkopf denied, please compile on local machine
if %ComputerName%==STACHLERKOPF exit /b
if %ComputerName%==MADRISAHORN echo compilation on madrisahorn denied, please compile on local machine 
if %ComputerName%==MADRISAHORN exit /b
echo on
set mods= MODULES.F90 MODE_WRITE.F90
set routines= CANOPY.F90 CHECK_SNOWPACK.F90 CUMULATE.F90 ^
DRIVE.F90 DUMP.F90 EBALFOR.F90 EBALSRF_SBG.F90 EBALSRF.F90 FRESH_SNOW_DENSITY.F90 FSM2.F90 ^
HS_FROM_SWE.F90 LUDCMP.F90 LWRADTOPO.F90 OPEN_FILES.F90 PHYSICS.F90 QSAT.F90 RADIATION.F90 ^
SETUP.F90 SNOW.F90 SNOW_ABLATION.F90 SNOWCOVERFRACTION.F90 SNOWTRAN3D.F90 SNOWSLIDE.F90 ^
SNOW_LAYERING.F90 SOIL.F90 SFEXCH.F90 SWE_FROM_HS.F90 SWRADTOPO.F90 TABLER.F90 ^
THERMAL.F90 TRIDIAG.F90
set optim= %1
set profil= %2
gfortran %mods% %routines% %optim% %profil% -cpp -ffpe-trap=overflow -o FSM2
del *.mod
move FSM2.exe ..\FSM2.exe
pause