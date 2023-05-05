!-----------------------------------------------------------------------
! Open files for reading and writing
!-----------------------------------------------------------------------
subroutine OPEN_FILES
  
use MODCONF, only: CANMOD, SNFRAC, SNTRAN, SNSLID, FOR_HN

use MODOUTPUT, only: &
  LIST_DIAG_RESULTS, WRITE_DIAG_TABLE, WRITE_DIAG_VARS, &
  LIST_STATE_RESULTS, WRITE_STATE_TABLE, WRITE_STATE_VARS

implicit none   

integer :: i,j, iwrite_pos,iwrite_count,iwrite_search

! Driving data
character(len=80) :: met_file_year
character(len=80) :: met_file_month
character(len=80) :: met_file_day
character(len=80) :: met_file_hour
character(len=80) :: met_file_SWb
character(len=80) :: met_file_SWd
character(len=80) :: met_file_LW
character(len=80) :: met_file_Sf
character(len=80) :: met_file_Rf
character(len=80) :: met_file_Ta
character(len=80) :: met_file_RH
character(len=80) :: met_file_Ua
character(len=80) :: met_file_Ps
character(len=80) :: met_file_Sf24h
character(len=80) :: met_file_Tvt  
character(len=80) :: met_file_Udir

! Initial state variables
! States needed in all runs
character(len=80) :: states_in_firstit
character(len=80) :: states_in_albs
character(len=80) :: states_in_Ds
character(len=80) :: states_in_fsnow
character(len=80) :: states_in_Nsnow
character(len=80) :: states_in_Sice
character(len=80) :: states_in_Sliq
character(len=80) :: states_in_Tss
character(len=80) :: states_in_Tsnow
character(len=80) :: states_in_Tsoil
character(len=80) :: states_in_histowet
character(len=80) :: landuse_skyvf
character(len=80) :: landuse_lat
character(len=80) :: landuse_lon
character(len=80) :: landuse_dem
character(len=80) :: landuse_slope
character(len=80) :: landuse_shd

! States specific to open runs
character(len=80) :: states_in_snowdepthmin
character(len=80) :: states_in_snowdepthmax
character(len=80) :: states_in_snowdepthhist
character(len=80) :: states_in_swemin
character(len=80) :: states_in_swemax
character(len=80) :: states_in_swehist
character(len=80) :: landuse_slopemu
character(len=80) :: landuse_xi
character(len=80) :: landuse_Ld

! States specific to forest runs
character(len=80) :: states_in_Qcan
character(len=80) :: states_in_Sveg
character(len=80) :: states_in_Tcan
character(len=80) :: states_in_Tveg
character(len=80) :: landuse_fveg
character(len=80) :: landuse_hcan
character(len=80) :: landuse_lai
character(len=80) :: landuse_vfhp
character(len=80) :: landuse_forest

! States specific to SNOWTRAN3D runs
character(len=80) :: states_in_dm_tot_subl
character(len=80) :: states_in_dm_tot_trans

! States specific to SnowSlide runs
character(len=80) :: states_in_dm_tot_slide
character(len=80) :: states_in_igds

! Final state variables
character(len=80) :: states_out_firstit
character(len=80) :: states_out_albs
character(len=80) :: states_out_Ds
character(len=80) :: states_out_fsnow
character(len=80) :: states_out_Nsnow
character(len=80) :: states_out_Sice
character(len=80) :: states_out_Sliq
character(len=80) :: states_out_Tss
character(len=80) :: states_out_Tsnow
character(len=80) :: states_out_Tsoil
character(len=80) :: states_out_histowet

character(len=80) :: states_out_snowdepthmin
character(len=80) :: states_out_snowdepthmax
character(len=80) :: states_out_snowdepthhist
character(len=80) :: states_out_swemin
character(len=80) :: states_out_swemax
character(len=80) :: states_out_swehist

character(len=80) :: states_out_Qcan   
character(len=80) :: states_out_Sveg   
character(len=80) :: states_out_Tcan   
character(len=80) :: states_out_Tveg   

character(len=80) :: states_out_dm_tot_subl
character(len=80) :: states_out_dm_tot_trans

character(len=80) :: states_out_dm_tot_slide
character(len=80) :: states_out_igds

! State variables for the HN model
character(len=80) :: states_out_Ds19
character(len=80) :: states_out_Ts19
character(len=80) :: states_out_Tsnow19
character(len=80) :: states_out_Tsoil19

! Result data
character(len=80) :: res_year
character(len=80) :: res_month
character(len=80) :: res_day
character(len=80) :: res_hour

! Driving data    
met_file_year  = 'drive_year.bin'
met_file_month = 'drive_month.bin'
met_file_day   = 'drive_day.bin'
met_file_hour  = 'drive_hour.bin'
met_file_SWb   = 'drive_sdrx.bin'
met_file_SWd   = 'drive_sdfx.bin'
met_file_LW    = 'drive_lwrx.bin'
met_file_Sf    = 'drive_snfx.bin'
met_file_Rf    = 'drive_rnfx.bin'
met_file_Ta    = 'drive_taix.bin'
met_file_RH    = 'drive_rhux.bin'
met_file_Ua    = 'drive_wnsx.bin'
met_file_Ps    = 'drive_paix.bin'
met_file_Sf24h = 'drive_snfh.bin'
if (CANMOD == 1) then
  met_file_Tvt   = 'drive_stdx.bin' 
endif
if (SNTRAN == 1) then
  met_file_Udir  = 'drive_wndx.bin'
endif

open(800, file = met_file_year,  form='unformatted', access='stream', status='old')
open(801, file = met_file_month, form='unformatted', access='stream', status='old')
open(802, file = met_file_day,   form='unformatted', access='stream', status='old')
open(803, file = met_file_hour,  form='unformatted', access='stream', status='old')
open(804, file = met_file_SWb,   form='unformatted', access='stream', status='old')
open(805, file = met_file_SWd,   form='unformatted', access='stream', status='old')
open(806, file = met_file_LW,    form='unformatted', access='stream', status='old')
open(807, file = met_file_Sf,    form='unformatted', access='stream', status='old')
open(808, file = met_file_Rf,    form='unformatted', access='stream', status='old')
open(809, file = met_file_Ta,    form='unformatted', access='stream', status='old')
open(810, file = met_file_RH,    form='unformatted', access='stream', status='old')
open(811, file = met_file_Ua,    form='unformatted', access='stream', status='old')
open(812, file = met_file_Ps,    form='unformatted', access='stream', status='old')
open(813, file = met_file_Sf24h, form='unformatted', access='stream', status='old')
if (CANMOD == 1) then
  open(814, file = met_file_Tvt, form='unformatted', access='stream', status='old')
endif
if (SNTRAN == 1) then
  open(815, file = met_file_Udir, form='unformatted', access='stream', status='old')
endif

! Initial state variables

! Common states 
landuse_skyvf           = 'landuse_skyvf.bin'
landuse_lat             = 'landuse_lat.bin'
landuse_lon             = 'landuse_lon.bin'
landuse_dem             = 'landuse_dem.bin'
states_in_firstit       = 'states_in_init.bin'
states_in_albs          = 'states_in_alse.bin'
states_in_Ds            = 'states_in_hsnl.bin'
states_in_fsnow         = 'states_in_scfe.bin'
states_in_Nsnow         = 'states_in_nsne.bin'
states_in_Sice          = 'states_in_sicl.bin'
states_in_Sliq          = 'states_in_slql.bin'
states_in_Tss           = 'states_in_tsfe.bin'
states_in_Tsnow         = 'states_in_tsnl.bin'
states_in_Tsoil         = 'states_in_tsll.bin'
states_in_histowet      = 'states_in_hiwl.bin'

! Open simulations
if (SNFRAC == 0 .or. SNFRAC == 2) then
  states_in_snowdepthmax  = 'states_in_hsmx.bin'
endif

if (SNFRAC == 0) then
  states_in_snowdepthmin  = 'states_in_hsmn.bin'
  states_in_snowdepthhist = 'states_in_hshs.bin'
  states_in_swemin        = 'states_in_swmn.bin'
  states_in_swemax        = 'states_in_swmx.bin'
  states_in_swehist       = 'states_in_swhs.bin'
  landuse_slopemu         = 'landuse_slopemu.bin'
  landuse_xi              = 'landuse_xi.bin'
  landuse_Ld              = 'landuse_Ld.bin'
endif

if (CANMOD == 1) then
! Forest simulations
  states_in_Qcan          = 'states_in_qcan.bin'
  states_in_Sveg          = 'states_in_sveg.bin'
  states_in_Tcan          = 'states_in_tcan.bin'
  states_in_Tveg          = 'states_in_tveg.bin'
  landuse_fveg            = 'landuse_fveg.bin'
  landuse_hcan            = 'landuse_hcan.bin'
  landuse_lai             = 'landuse_lai.bin'
  landuse_vfhp            = 'landuse_vfhp.bin' 
  landuse_forest          = 'landuse_forest.bin' 
endif

if (SNTRAN == 1) then
! SNOWTRAN3D simulations
  states_in_dm_tot_subl   = 'states_in_sblt.bin'
  states_in_dm_tot_trans  = 'states_in_trat.bin'
  landuse_Ld              = 'landuse_Ld.bin'
endif

if (SNSLID == 1) then
! SnowSlide simulations
  landuse_slope           = 'landuse_slope.bin'
  landuse_shd             = 'landuse_shd.bin'
  states_in_dm_tot_slide  = 'states_in_sldt.bin'
  states_in_igds          = 'states_in_igds.bin'
endif

! NOTE / GM: the numbering of the files is not continuous, but it is currently consistent with JIM. 
! If we have no more interest in keeping this consistency, I would suggest to renumber the files to have continuous numbers 
open(1100, file = states_in_firstit,  form='unformatted', access='stream', status='old')
open(1101, file = states_in_albs,     form='unformatted', access='stream', status='old')
open(1102, file = states_in_Ds,       form='unformatted', access='stream', status='old')
open(1103, file = states_in_fsnow,    form='unformatted', access='stream', status='old')
open(1104, file = states_in_Nsnow,    form='unformatted', access='stream', status='old')
open(1106, file = states_in_Sice,     form='unformatted', access='stream', status='old')
open(1107, file = states_in_Sliq,     form='unformatted', access='stream', status='old')
open(1116, file = states_in_Tss,      form='unformatted', access='stream', status='old')
open(1119, file = states_in_Tsnow,    form='unformatted', access='stream', status='old')
open(1120, file = states_in_Tsoil,    form='unformatted', access='stream', status='old')
open(1121, file = states_in_histowet, form='unformatted', access='stream', status='old')
open(1123, file = landuse_skyvf,      form='unformatted', access='stream', status='old')
open(1127, file = landuse_lat,        form='unformatted', access='stream', status='old')
open(1128, file = landuse_lon,        form='unformatted', access='stream', status='old')
open(1129, file = landuse_dem,        form='unformatted', access='stream', status='old')

if (SNFRAC == 0 .or. SNFRAC == 2) then
  open(1110, file = states_in_snowdepthmax, form='unformatted', access='stream', status='old')
endif

if (SNFRAC == 0 ) then  
  open(1109, file = states_in_snowdepthmin,  form='unformatted', access='stream', status='old')
  open(1111, file = states_in_snowdepthhist, form='unformatted', access='stream', status='old')
  open(1113, file = states_in_swemin,        form='unformatted', access='stream', status='old')
  open(1114, file = states_in_swemax,        form='unformatted', access='stream', status='old')
  open(1115, file = states_in_swehist,       form='unformatted', access='stream', status='old')
  open(1124, file = landuse_slopemu,         form='unformatted', access='stream', status='old')
  open(1125, file = landuse_xi,              form='unformatted', access='stream', status='old')
  open(1126, file = landuse_Ld,              form='unformatted', access='stream', status='old')
endif

if (CANMOD == 1) then
  open(1130, file = states_in_Qcan, form='unformatted', access='stream', status='old') 
  open(1131, file = states_in_Sveg, form='unformatted', access='stream', status='old')
  open(1132, file = states_in_Tcan, form='unformatted', access='stream', status='old')
  open(1133, file = states_in_Tveg, form='unformatted', access='stream', status='old')
  open(1134, file = landuse_fveg,   form='unformatted', access='stream', status='old')
  open(1135, file = landuse_hcan,   form='unformatted', access='stream', status='old')
  open(1136, file = landuse_lai,    form='unformatted', access='stream', status='old')
  open(1137, file = landuse_vfhp,   form='unformatted', access='stream', status='old')
  open(1138, file = landuse_forest, form='unformatted', access='stream', status='old')
endif

if (SNTRAN == 1) then
  open(1140, file = states_in_dm_tot_subl, form='unformatted', access='stream', status='old') 
  open(1141, file = states_in_dm_tot_trans, form='unformatted', access='stream', status='old')
  open(1126, file = landuse_Ld,              form='unformatted', access='stream', status='old')
endif

if (SNSLID == 1) then
  open(1139, file = landuse_slope,          form='unformatted', access='stream', status='old')
  open(1122, file = landuse_shd,            form='unformatted', access='stream', status='old')
  open(1142, file = states_in_dm_tot_slide, form='unformatted', access='stream', status='old')
  open(1143, file = states_in_igds,         form='unformatted', access='stream', status='old')
endif

! Final state variables
states_out_firstit       = 'states_out_init.bin'
states_out_albs          = 'states_out_alse.bin'
states_out_Ds            = 'states_out_hsnl.bin'
states_out_fsnow         = 'states_out_scfe.bin'
states_out_Nsnow         = 'states_out_nsne.bin'
states_out_Sice          = 'states_out_sicl.bin'
states_out_Sliq          = 'states_out_slql.bin'
states_out_Tss           = 'states_out_tsfe.bin'
states_out_Tsnow         = 'states_out_tsnl.bin'
states_out_Tsoil         = 'states_out_tsll.bin'
states_out_histowet      = 'states_out_hiwl.bin'

if (SNFRAC == 0 .or. SNFRAC == 2) then
  states_out_snowdepthmax  = 'states_out_hsmx.bin'
endif

if (SNFRAC == 0 ) then
  states_out_snowdepthmin  = 'states_out_hsmn.bin'
  states_out_snowdepthhist = 'states_out_hshs.bin'
  states_out_swemin        = 'states_out_swmn.bin'
  states_out_swemax        = 'states_out_swmx.bin'
  states_out_swehist       = 'states_out_swhs.bin'
endif

if (CANMOD == 1) then
  states_out_Qcan          = 'states_out_qcan.bin'
  states_out_Sveg          = 'states_out_sveg.bin'
  states_out_Tcan          = 'states_out_tcan.bin'
  states_out_Tveg          = 'states_out_tveg.bin'
endif

if (SNTRAN == 1) then
  states_out_dm_tot_subl   = 'states_out_sblt.bin'
  states_out_dm_tot_trans  = 'states_out_trat.bin'
endif

if (SNSLID == 1) then
  states_out_dm_tot_slide  = 'states_out_sldt.bin'
  states_out_igds          = 'states_out_igds.bin'
endif

! states at 19h for subsequent hn runs.
states_out_Ds19          = 'states_out_hsn9.bin'
states_out_Ts19          = 'states_out_tsr9.bin'
states_out_Tsnow19       = 'states_out_tsn9.bin'
states_out_Tsoil19       = 'states_out_tsl9.bin'

open(1200, file = states_out_firstit,  form='unformatted', access='stream')
open(1201, file = states_out_albs,     form='unformatted', access='stream')
open(1202, file = states_out_Ds,       form='unformatted', access='stream')
open(1203, file = states_out_fsnow,    form='unformatted', access='stream')
open(1204, file = states_out_Nsnow,    form='unformatted', access='stream')
open(1206, file = states_out_Sice,     form='unformatted', access='stream')
open(1207, file = states_out_Sliq,     form='unformatted', access='stream')
open(1216, file = states_out_Tss,      form='unformatted', access='stream')
open(1219, file = states_out_Tsnow,    form='unformatted', access='stream')
open(1220, file = states_out_Tsoil,    form='unformatted', access='stream')
open(1221, file = states_out_histowet, form='unformatted', access='stream')

if (SNFRAC == 0 .or. SNFRAC == 2) then
  open(1210, file = states_out_snowdepthmax, form='unformatted', access='stream')
endif

if (SNFRAC == 0) then 
  open(1209, file = states_out_snowdepthmin,  form='unformatted', access='stream')
  open(1211, file = states_out_snowdepthhist, form='unformatted', access='stream')
  open(1213, file = states_out_swemin,        form='unformatted', access='stream')
  open(1214, file = states_out_swemax,        form='unformatted', access='stream')
  open(1215, file = states_out_swehist,       form='unformatted', access='stream')
endif

if (CANMOD == 1) then
  open(1223, file = states_out_Qcan, form='unformatted', access='stream')  
  open(1224, file = states_out_Sveg, form='unformatted', access='stream')
  open(1225, file = states_out_Tcan, form='unformatted', access='stream')
  open(1226, file = states_out_Tveg, form='unformatted', access='stream')
endif

if (SNTRAN == 1) then
  open(1240, file = states_out_dm_tot_subl,  form='unformatted', access='stream')  
  open(1241, file = states_out_dm_tot_trans, form='unformatted', access='stream')
endif

if (SNSLID == 1) then
  open(1242, file = states_out_dm_tot_slide, form='unformatted', access='stream')
  open(1243, file = states_out_igds,         form='unformatted', access='stream')
endif

! states at 19h for subsequent hn runs.
! BC yes you're right they should not be always written...
if (FOR_HN) then
  open(1227, file = states_out_Ds19,    form='unformatted', access='stream')  
  open(1228, file = states_out_Ts19,    form='unformatted', access='stream')
  open(1229, file = states_out_Tsnow19, form='unformatted', access='stream')
  open(1230, file = states_out_Tsoil19, form='unformatted', access='stream')
endif

! Result data
res_year        = 'res_year.bin'
res_month       = 'res_month.bin'
res_day         = 'res_day.bin'
res_hour        = 'res_hour.bin'
open(1401, file = res_year,  form='unformatted', access='stream')
open(1402, file = res_month, form='unformatted', access='stream')
open(1403, file = res_day,   form='unformatted', access='stream')
open(1404, file = res_hour,  form='unformatted', access='stream')

! Flexible output
!   1. check if the used asked to write variable into output
!   2. if so:
!     - set .true. in the table of vars.
!     - open the corresponding file

! diagnostic variables
iwrite_search = count(LIST_DIAG_RESULTS /= '            ')
WRITE_DIAG_TABLE(:) = .false.
iwrite_count = 1
do i=1, size(WRITE_DIAG_VARS)
  do  j=1,iwrite_search
    iwrite_pos = index(trim(WRITE_DIAG_VARS(i)), trim(LIST_DIAG_RESULTS(j)))
    if (iwrite_pos == 1) then
      WRITE_DIAG_TABLE(iwrite_count) = .true.
      open(1404 + i, file = 'res_'//trim(WRITE_DIAG_VARS(i))//'.bin',   form='unformatted', access='stream')
      exit
    endif
  end do
  iwrite_count =  iwrite_count + 1
end do

! state variables
iwrite_search = count(LIST_STATE_RESULTS /= '            ')
WRITE_STATE_TABLE(:) = .false.
iwrite_count = 1
do i=1, size(WRITE_STATE_VARS)
  do  j=1,iwrite_search
    iwrite_pos = index(trim(WRITE_STATE_VARS(i)), trim(LIST_STATE_RESULTS(j)))
    if (iwrite_pos == 1) then
      WRITE_STATE_TABLE(iwrite_count) = .true.
      open(1500 + i, file = 'res_'//trim(WRITE_STATE_VARS(i))//'.bin',   form='unformatted', access='stream')
      exit
    endif
  end do
  iwrite_count =  iwrite_count + 1
end do

end subroutine OPEN_FILES
