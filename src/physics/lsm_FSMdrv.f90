
!!
!!----------------------------------------------------------
module module_sf_FSMdrv
    use module_sf_noahdrv,   only : lsm_noah, lsm_noah_init
    use module_water_simple, only : water_simple
    use mod_atm_utilities,   only : sat_mr
    use time_object,         only : Time_type
    use data_structures
    use mod_wrf_constants,   only : piconst, XLS
    use icar_constants
    use options_interface,   only : options_t
    use domain_interface,    only : domain_t
    use io_routines,         only : io_write, io_read, io_add_attribute
    use FSM_interface , only:  FSM_SETUP,FSM_DRIVE,FSM_PHYSICS, FSM_SNOWSLIDE, FSM_SNOWSLIDE_END, FSM_CUMULATE_SD, FSM_SNOWTRAN_SETUP, FSM_SNOWTRAN_FLUXES, FSM_SNOWTRAN_ACCUM
    use FSM_interface , only:  Nx_HICAR, Ny_HICAR,lat_HICAR,lon_HICAR,terrain_HICAR,dx_HICAR,slope_HICAR,shd_HICAR
    use FSM_interface, only: &
      year,          &
      month,         &
      day,           &
      hour,          &
      dt,            &
      LW,            &
      Ps,            &
      Qa,            &
      Rf,            &
      Sdif,          &
      Sdir,          &
      Sf,            &
      Sf24h,         &
      Ta,            &
      Ua,            &
      Udir
    use FSM_interface, only: &
      Esrf_,         &
      Gsoil_,        &
      H_,            &
      LE_,           &
      Melt_,         &
      Rnet_,         &
      Roff_,         &
      snowdepth_,    &
      SWE_,          &  
      KH_,           &  
      meltflux_out_, &
      Sliq_out_,     &
      dm_salt_,      &
      dm_susp_,      &
      dm_subl_,      &
      dm_slide_!,     &
      !Qsalt_u,       &
      !Qsalt_v

    use FSM_interface, only: &
      firstit,       &
      Tsrf,          &
      Tsnow,         &
      Sice,          &
      Sliq,          &
      Ds,            &
      fsnow,         &
      Nsnow,         &
      Tsoil,         &
      albs,          &
      theta,         &
      z0sn
    use FSM_interface, only: SNTRAN, SNSLID
    
    implicit none

    private
    public :: sm_FSM_init,sm_FSM
    
    integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions
    
    type(Time_type) :: last_output
    real            :: last_snowslide
   
    real, allocatable :: &
        snowfall_sum(:,:),       & !aggregated per output interval
        rainfall_sum(:,:),       & !aggregated per output interval
        Roff_sum(:,:),           & !aggregated per output interval
        meltflux_out_sum(:,:)      !aggregated per output interval

contains
 
    subroutine sm_FSM_init(domain,options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        integer :: i,j
        !!
        
        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        
        !if (SNTRAN+SNSLID > 0) then
        its = domain%grid%its-1
        ite = domain%grid%ite+1
        jts = domain%grid%jts-1
        jte = domain%grid%jte+1
        !else
        !    its = domain%grid%its
        !    ite = domain%grid%ite
        !    jts = domain%grid%jts
        !    jte = domain%grid%jte
        !endif
        !!       
        
        Ny_HICAR=ite-its+1
        Nx_HICAR=jte-jts+1
        
        last_output = options%parameters%start_time
        last_snowslide = 4000
        !!
        allocate(lat_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(lon_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(terrain_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(slope_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(shd_HICAR(Nx_HICAR,Ny_HICAR))

        lat_HICAR=TRANSPOSE(domain%latitude%data_2d(its:ite,jts:jte))
        lon_HICAR=TRANSPOSE(domain%longitude%data_2d(its:ite,jts:jte))
        terrain_HICAR=TRANSPOSE(domain%terrain%data_2d(its:ite,jts:jte))
        if (associated(domain%slope_angle%data_2d)) then
            slope_HICAR=TRANSPOSE(domain%slope_angle%data_2d(its:ite,jts:jte))*180.0/piconst !Convert from radians to degrees
        else
            slope_HICAR = 0.
        endif
        if (associated(domain%shd%data_2d)) then
            shd_HICAR=TRANSPOSE(domain%shd%data_2d(its:ite,jts:jte))
        else
            shd_HICAR = 10000.
        endif
        
        dx_HICAR=domain%dx
        !!
        allocate(Esrf_(Nx_HICAR,Ny_HICAR)); Esrf_=0.
        allocate(Gsoil_(Nx_HICAR,Ny_HICAR));Gsoil_=0.
        allocate(H_(Nx_HICAR,Ny_HICAR)); H_=0.
        allocate(LE_(Nx_HICAR,Ny_HICAR)); LE_=0.
        allocate(Melt_(Nx_HICAR,Ny_HICAR)); Melt_=0.
        allocate(Rnet_(Nx_HICAR,Ny_HICAR)); Rnet_=0.
        allocate(Roff_(Nx_HICAR,Ny_HICAR)); Roff_=0.
        allocate(snowdepth_(Nx_HICAR,Ny_HICAR)); snowdepth_=0.
        allocate(SWE_(Nx_HICAR,Ny_HICAR)); SWE_=0.0
        allocate(KH_(Nx_HICAR,Ny_HICAR)); KH_=0.0
        allocate(meltflux_out_(Nx_HICAR,Ny_HICAR)); meltflux_out_=0.
        allocate(Sliq_out_(Nx_HICAR,Ny_HICAR)); Sliq_out_=0.
        allocate(dm_salt_(Nx_HICAR,Ny_HICAR)); dm_salt_=0.
        allocate(dm_susp_(Nx_HICAR,Ny_HICAR)); dm_susp_=0.
        allocate(dm_subl_(Nx_HICAR,Ny_HICAR)); dm_subl_=0.
        allocate(dm_slide_(Nx_HICAR,Ny_HICAR)); dm_slide_=0.
        !allocate(Qsalt_u(Nx_HICAR,Ny_HICAR)); Qsalt_u=0.
        !allocate(Qsalt_v(Nx_HICAR,Ny_HICAR)); Qsalt_v=0.
        !!
        allocate(snowfall_sum(Nx_HICAR,Ny_HICAR)); snowfall_sum=0.
        allocate(rainfall_sum(Nx_HICAR,Ny_HICAR)); rainfall_sum=0.
        allocate(Roff_sum(Nx_HICAR,Ny_HICAR)); Roff_sum=0.
        allocate(meltflux_out_sum(Nx_HICAR,Ny_HICAR)); meltflux_out_sum=0.   
        !!
        call FSM_SETUP()
        !!        
        
        !! MJ added this block to read in while we use restart file:
        if (options%parameters%restart) then
            !! giving feedback to HICAR
            Tsrf = TRANSPOSE(domain%skin_temperature%data_2d(its:ite,jts:jte))
            if (options%lsm_options%monthly_albedo) then
                albs = TRANSPOSE(domain%albedo%data_3d(its:ite, domain%model_time%month, jts:jte))
            else
                albs = TRANSPOSE(domain%albedo%data_3d(its:ite, 1, jts:jte))
            endif
            
            fsnow = TRANSPOSE(domain%fsnow%data_2d(its:ite,jts:jte))
            Nsnow = TRANSPOSE(domain%Nsnow%data_2d(its:ite,jts:jte))                        
            !!
            do i=1,3
                Tsnow(i,:,:) = TRANSPOSE(domain%snow_temperature%data_3d(its:ite,i,jts:jte))
                Sice(i,:,:) = TRANSPOSE(domain%Sice%data_3d(its:ite,i,jts:jte))
                Sliq(i,:,:) = TRANSPOSE(domain%Sliq%data_3d(its:ite,i,jts:jte))
                Ds(i,:,:) = TRANSPOSE(domain%Ds%data_3d(its:ite,i,jts:jte))
            enddo
            do i=1,4
                Tsoil(i,:,:) = TRANSPOSE(domain%soil_temperature%data_3d(its:ite,i,jts:jte))
                theta(i,:,:) = TRANSPOSE(domain%soil_water_content%data_3d(its:ite,i,jts:jte))
            enddo
        endif
        !!
        do j = 1, Ny_HICAR
            do i = 1, Nx_HICAR
                !if (this_image()==1) write(*,*) "  albsH, albsF  ",i, j, domain%albs%data_2d(i+its-1,j+jts-1), albs(i,j)
            end do
        end do
        !SYNC ALL
    end subroutine sm_FSM_init

    subroutine sm_FSM(domain,options,lsm_dt,current_rain,current_snow,windspd)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        !!       
        real, intent(in) :: &
            lsm_dt          ! Timestep (s) from HICAR for lsm
        real, dimension(ims:ime,jms:jme), intent(in) :: &   !! Note that: when the input is an array with whatever indexing either using (:,:) or (N,M), it always starts from index 1 in all dims
            current_rain,                 &! rainfall in kg m-2
            current_snow,                 &! snowfall in kg m-2
            windspd                        ! Wind speed (m/s)

        integer :: i,j,k, hj, hi, i_s, i_e, j_s, j_e
        real :: Delta_t
        real, dimension(its:ite,jts:jte) :: SWE_pre
        real, dimension(Nx_HICAR,Ny_HICAR) :: SD_0, Sice_0, SD_0_buff, Sice_0_buff, Qsalt_u, Qsalt_v
        logical, dimension(Nx_HICAR,Ny_HICAR) :: aval
        logical :: first_SLIDE, do_snowslide
        character(len=1024) :: filename
        
        !if (SNTRAN+SNSLID > 0) then
        do_snowslide = .False.
        if (SNSLID > 0) then
            do_snowslide = (last_snowslide > 3600.) !Hard code to just be called every hour...
        endif

        j_s = 2
        i_s = 2
        j_e = Nx_HICAR-1
        i_e = Ny_HICAR-1
        !else
        !    j_s = 1
        !    i_s = 1
        !    j_e = Ny_HICAR
        !    i_e = Nx_HICAR
        !endif
        
        !!
        !! giving the required input from HICAR to FSM
        year=real(domain%model_time%year)
        month=real(domain%model_time%month)
        day=real(domain%model_time%day)
        hour=real(domain%model_time%hour)
        dt=lsm_dt
        
        
        LW=TRANSPOSE(domain%longwave%data_2d(its:ite,jts:jte))
        Ps=TRANSPOSE(domain%surface_pressure%data_2d(its:ite,jts:jte))
        Rf=TRANSPOSE(current_rain(its:ite,jts:jte))
        !!
        Sdir=TRANSPOSE(domain%shortwave_direct%data_2d(its:ite,jts:jte))  !Sdir=domain%shortwave%data_2d(its:ite,jts:jte)
        Sdif=TRANSPOSE(domain%shortwave_diffuse%data_2d(its:ite,jts:jte)) !Sdif=0.0
        !!
        if (options%parameters%factor_p_var == "") then 
            !if (this_image()==1) write(*,*) "facto_p is not read...FSM"
            Sf=TRANSPOSE(current_snow(its:ite,jts:jte))
        endif
        if (options%parameters%factor_p_var /= "") then 
            !if (this_image()==1) write(*,*) "facto_p is read...FSM"
            Sf=TRANSPOSE(current_snow(its:ite,jts:jte)*domain%factor_p%data_2d(its:ite,jts:jte))
        endif
        !
        Sf24h=TRANSPOSE(domain%snowfall_tstep%data_2d(its:ite,jts:jte))
        !
        Ta= TRANSPOSE(domain%temperature%data_3d(its:ite,domain%grid%kms,jts:jte))!domain%temperature_2m%data_2d(its:ite,jts:jte)
        Qa= TRANSPOSE(domain%water_vapor%data_3d(its:ite,domain%grid%kms,jts:jte))!domain%humidity_2m%data_2d(its:ite,jts:jte)
        Ua=TRANSPOSE(windspd(its:ite,jts:jte))
        Udir = ATAN2(TRANSPOSE(domain%u_10m%data_2d(its:ite,jts:jte)/windspd(its:ite,jts:jte)), &
                     TRANSPOSE(domain%v_10m%data_2d(its:ite,jts:jte)/windspd(its:ite,jts:jte))) 
        Udir = (Udir * 180/piconst)
        where(Udir<0) Udir=Udir+360

        if ((domain%model_time%seconds() - dt <= last_output%seconds()) .and. &
            (domain%model_time%seconds()   >=    last_output%seconds())) then
            !If we are the first call since the last output, reset the per-output counters
            domain%dm_slide%data_2d = 0.
            domain%dm_salt%data_2d  = 0.
            domain%dm_susp%data_2d  = 0.
            domain%dm_subl%data_2d  = 0.
            last_output = last_output + options%io_options%output_dt
        endif
        
        SWE_pre = domain%snow_water_equivalent%data_2d(its:ite,jts:jte)
                !!  
        !! FSM processing      
        call FSM_DRIVE()
        call FSM_PHYSICS()
                        
        !Call Snowtran here -- must be done here and not in FSM since we need control over the parallelization of the routine
        if (SNTRAN > 0) then        
            call exch_FSM_state_vars(domain)
            
            call FSM_SNOWTRAN_SETUP(Qsalt_u,Qsalt_v)
            
            !First guess for fluxes
            call FSM_SNOWTRAN_FLUXES(Qsalt_u,Qsalt_v)
            !Exchange fluxes between processes
            call exch_SNTRAN_Qsalt(domain,Qsalt_u,Qsalt_v)
            !Recalculate fluxes with intermediate values from neighbors
            call FSM_SNOWTRAN_FLUXES(Qsalt_u,Qsalt_v)

            call FSM_SNOWTRAN_ACCUM(Qsalt_u,Qsalt_v)
        endif

        
        !Call Snowslide here -- must be done here and not in FSM since we need control over the parallelization of the routine
        if (do_snowslide .and. SNSLID > 0) then
            SD_0 = 0.0
            Sice_0 = 0.0
            SD_0_buff = 0.0
            Sice_0_buff = 0.0

            dm_slide_ = 0.0
            first_SLIDE = .True.
            aval = .False.
            !Snowslide needs the corner snow depth information from corner neighbor processes
            call exch_FSM_state_vars(domain,corners_in=.True.)
            
            do i=1,10
                call FSM_SNOWSLIDE(SD_0,Sice_0,SD_0_buff,Sice_0_buff,aval,first_SLIDE,dm_slide_)
                
                ! Copy interior buffer for exchange
                call exch_SLIDE_buffers(domain,SD_0_buff,Sice_0_buff)
                
                !Now, where the buffer is positive for halo cells, record that an avalanche was passed here
                aval = .False.
                where(SD_0_buff > SD_0) aval=.True.
                aval(j_s:j_e,i_s:i_e) = .False.
                
                !Must accumulate slide changes here, since we will loop over calls to snowslide
                domain%dm_slide%data_2d(domain%its:domain%ite,domain%jts:domain%jte) = &
                        domain%dm_slide%data_2d(domain%its:domain%ite,domain%jts:domain%jte) + TRANSPOSE(dm_slide_(j_s:j_e,i_s:i_e))

                ! trade over buffer to halo cells, in case avalanching snow was passed on previous iteration.
                SD_0(1,:) = SD_0_buff(1,:)
                SD_0(Nx_HICAR,:) = SD_0_buff(Nx_HICAR,:)
                SD_0(:,1) = SD_0_buff(:,1)
                SD_0(:,Ny_HICAR) = SD_0_buff(:,Ny_HICAR)

                Sice_0(1,:) = Sice_0_buff(1,:)
                Sice_0(Nx_HICAR,:) = Sice_0_buff(Nx_HICAR,:)
                Sice_0(:,1) = Sice_0_buff(:,1)
                Sice_0(:,Ny_HICAR) = Sice_0_buff(:,Ny_HICAR)


            enddo
            call FSM_SNOWSLIDE_END(SD_0,Sice_0)
            last_snowslide = 0
        endif
        
        
        if (SNSLID > 0) then
            last_snowslide = last_snowslide + lsm_dt
        endif

        !The End.
        if ( (SNTRAN > 0) .or. do_snowslide) call exch_FSM_state_vars(domain)
        call FSM_CUMULATE_SD()
        
        !! giving feedback to HICAR -- should only be done for snow-covered cells, or cells which were just snowed on
        do j=j_s,j_e
            do i=i_s,i_e
                hj = j-j_s+domain%jts
                hi = i-i_s+domain%its
                if ( SWE_(j,i) > 0 .and. .not.(options%physics%watersurface==kWATER_SIMPLE .and. domain%land_mask(hi,hj)==kLC_WATER)) then
                    domain%sensible_heat%data_2d(hi,hj)=H_(j,i)
                    domain%latent_heat%data_2d(hi,hj)=LE_(j,i)
                    domain%snow_water_equivalent%data_2d(hi,hj)=SWE_(j,i)
                    if (options%lsm_options%monthly_albedo) then
                        domain%albedo%data_3d(hi, domain%model_time%month, hj) = albs(j,i)
                    else
                        domain%albedo%data_3d(hi, 1, hj) = albs(j,i)
                    endif

                    !
                    domain%skin_temperature%data_2d(hi,hj)=Tsrf(j,i)
                    domain%snow_height%data_2d(hi,hj)=snowdepth_(j,i)
                    domain%Sliq_out%data_2d(hi,hj)=Sliq_out_(j,i)
                    domain%fsnow%data_2d(hi,hj)=fsnow(j,i)
                    domain%Nsnow%data_2d(hi,hj)=Nsnow(j,i)
                    !
                    domain%qfx%data_2d(hi,hj)=Esrf_(j,i)
                    domain%chs%data_2d(hi,hj)=KH_(j,i)
                    domain%roughness_z0%data_2d(hi,hj) = z0sn*fsnow(j,i)+(1-fsnow(j,i))*domain%roughness_z0%data_2d(hi,hj)
                    
                    do k=1,3
                        domain%snow_temperature%data_3d(hi,k,hj) = Tsnow(k,j,i)
                        domain%Sice%data_3d(hi,k,hj) = Sice(k,j,i)
                        domain%Sliq%data_3d(hi,k,hj) = Sliq(k,j,i)
                        domain%Ds%data_3d(hi,k,hj) = Ds(k,j,i)
                    enddo
                    do k=1,4
                        domain%soil_temperature%data_3d(hi,k,hj) = Tsoil(k,j,i)
                        domain%soil_water_content%data_3d(hi,k,hj)=theta(k,j,i)
                    enddo
                
                    if (SNTRAN>0) then
                        ! Convert to rate 1/s
                        domain%dm_salt%data_2d(hi,hj)=domain%dm_salt%data_2d(hi,hj) + dm_salt_(j,i)/dt 
                        domain%dm_susp%data_2d(hi,hj)= Qsalt_u(j,i)!domain%dm_susp%data_2d(hi,hj) + dm_susp_(j,i)/dt
                        domain%dm_subl%data_2d(hi,hj)= Qsalt_v(j,i)!domain%dm_subl%data_2d(hi,hj) + dm_subl_(j,i)/dt 
                        
                        !Add sublimated snow to latent heat flux. 
                        !Sometimes FSM returns NaN values for blowing snow sublimation, so mask those out here
                        if (abs(dm_subl_(j,i))>1) dm_subl_(j,i) = 0.0
                        domain%latent_heat%data_2d(hi,hj)   = domain%latent_heat%data_2d(hi,hj)   + (-dm_subl_(j,i))*XLS/dt 
                        domain%sensible_heat%data_2d(hi,hj) = domain%sensible_heat%data_2d(hi,hj) + ( dm_subl_(j,i))*XLS/dt 
                    endif
                endif
            enddo
        enddo
        
        !If we just melted out the snow in this step, set this variable. This is to ensure good hand-off from snow-covered to not with NoahMP
        where(SWE_pre(domain%its:domain%ite,domain%jts:domain%jte) > 0 .and. &
                domain%snow_water_equivalent%data_2d(domain%its:domain%ite,domain%jts:domain%jte)==0) 
            domain%ground_surf_temperature%data_2d(domain%its:domain%ite,domain%jts:domain%jte)=273.15
        end where
        
        ! Let FSM know that it has done the first iteration 
        if (firstit==1) firstit=0
        
        !These are FSM2 Diagnostic/output vars, so we can update them everywhere
        !domain%Nsnow%data_2d(its:ite,jts:jte)=Nsnow
        !!
        !snowfall_sum=snowfall_sum+Sf
        !rainfall_sum=rainfall_sum+Rf
        !Roff_sum=Roff_sum+Roff_
        !meltflux_out_sum=meltflux_out_sum+meltflux_out_
        

        !Delta_t=mod(domain%model_time%seconds(),options%io_options%out_dt)
        !if ( abs(options%io_options%out_dt-(Delta_t+dt)) <= 1.e-3 ) then
        !    if (this_image()==1) write(*,*) "resetting/aggregating vars e.g. runoff during t-1->t"!, Delta_t,Delta_t+dt
        !    !!
        !    domain%rainfall_tstep%data_2d(its:ite,jts:jte)=rainfall_sum
        !    domain%snowfall_tstep%data_2d(its:ite,jts:jte)=snowfall_sum
        !    domain%runoff_tstep%data_2d(its:ite,jts:jte)=Roff_sum
        !    domain%meltflux_out_tstep%data_2d(its:ite,jts:jte)=meltflux_out_sum
        !    !! reseting the container to zero for next output interval
        !    rainfall_sum = 0.
        !    snowfall_sum = 0.                
        !    Roff_sum = 0.
        !    meltflux_out_sum = 0.
        !endif
        
        !!
        !SYNC ALL  
    end subroutine sm_FSM
    
    
    subroutine exch_FSM_state_vars(domain,corners_in)
        implicit none
        
        type(domain_t), intent(inout) :: domain
        logical, optional, intent(in) :: corners_in
        
        integer :: i
        logical :: corners
        
        corners=.False.
        if (present(corners_in)) corners=corners_in
        
        domain%fsnow%data_2d(domain%its:domain%ite,domain%jts:domain%jte) = TRANSPOSE(fsnow(2:Nx_HICAR-1,2:Ny_HICAR-1))
        domain%Nsnow%data_2d(domain%its:domain%ite,domain%jts:domain%jte) = TRANSPOSE(Nsnow(2:Nx_HICAR-1,2:Ny_HICAR-1))                        
        !!
        do i=1,3
            domain%snow_temperature%data_3d(domain%its:domain%ite,i,domain%jts:domain%jte) = TRANSPOSE(Tsnow(i,2:Nx_HICAR-1,2:Ny_HICAR-1))
            domain%Sice%data_3d(domain%its:domain%ite,i,domain%jts:domain%jte) = TRANSPOSE(Sice(i,2:Nx_HICAR-1,2:Ny_HICAR-1))
            domain%Sliq%data_3d(domain%its:domain%ite,i,domain%jts:domain%jte) = TRANSPOSE(Sliq(i,2:Nx_HICAR-1,2:Ny_HICAR-1))
            domain%Ds%data_3d(domain%its:domain%ite,i,domain%jts:domain%jte) = TRANSPOSE(Ds(i,2:Nx_HICAR-1,2:Ny_HICAR-1))
        enddo


        call domain%halo_2d_exchange_batch()      
        call domain%halo_3d_exchange_batch(exch_var_only=.True.)      

        fsnow = TRANSPOSE(domain%fsnow%data_2d(its:ite,jts:jte))
        Nsnow = TRANSPOSE(domain%Nsnow%data_2d(its:ite,jts:jte))                        
        !!
        do i=1,3
            Tsnow(i,:,:) = TRANSPOSE(domain%snow_temperature%data_3d(its:ite,i,jts:jte))
            Sice(i,:,:) = TRANSPOSE(domain%Sice%data_3d(its:ite,i,jts:jte))
            Sliq(i,:,:) = TRANSPOSE(domain%Sliq%data_3d(its:ite,i,jts:jte))
            Ds(i,:,:) = TRANSPOSE(domain%Ds%data_3d(its:ite,i,jts:jte))
        enddo
        
        if (corners) then
            ! Need to handle corner exchanges necesarry for snowslide
            if (.not.(domain%south_boundary) .and. .not.(domain%west_boundary)) then
                !DIR$ PGAS DEFER_SYNC
                domain%east_in_2d(1,1,1)[domain%southwest_neighbor] = fsnow(2,2)
                !DIR$ PGAS DEFER_SYNC
                domain%east_in_3d(1,1,1:3,1)[domain%southwest_neighbor] = Ds(1:3,2,2)
            endif
            if (.not.(domain%north_boundary) .and. .not.(domain%west_boundary)) then
                !DIR$ PGAS DEFER_SYNC
                domain%south_in_2d(1,1,1)[domain%northwest_neighbor] = fsnow(Nx_HICAR-1,2)
                !DIR$ PGAS DEFER_SYNC
                domain%south_in_3d(1,1,1:3,1)[domain%northwest_neighbor] = Ds(1:3,Nx_HICAR-1,2)
            endif
            if (.not.(domain%south_boundary) .and. .not.(domain%east_boundary)) then
                !DIR$ PGAS DEFER_SYNC
                domain%north_in_2d(1,1,1)[domain%southeast_neighbor] = fsnow(2,Ny_HICAR-1)
                !DIR$ PGAS DEFER_SYNC
                domain%north_in_3d(1,1,1:3,1)[domain%southeast_neighbor] = Ds(1:3,2,Ny_HICAR-1)
            endif
            if (.not.(domain%north_boundary) .and. .not.(domain%east_boundary)) then
                !DIR$ PGAS DEFER_SYNC
                domain%west_in_2d(1,1,1)[domain%northeast_neighbor] = fsnow(Nx_HICAR-1,Ny_HICAR-1)
                !DIR$ PGAS DEFER_SYNC
                domain%west_in_3d(1,1,1:3,1)[domain%northeast_neighbor] = Ds(1:3,Nx_HICAR-1,Ny_HICAR-1)
            endif

            sync images ( domain%corner_neighbors )

            if (.not.(domain%south_boundary) .and. .not.(domain%west_boundary)) then
                fsnow(1,1) = domain%west_in_2d(1,1,1)
                Ds(1:3,1,1) = domain%west_in_3d(1,1,1:3,1)
            endif
            if (.not.(domain%north_boundary) .and. .not.(domain%west_boundary)) then
                fsnow(Nx_HICAR,1) = domain%north_in_2d(1,1,1)
                Ds(1:3,Nx_HICAR,1) = domain%north_in_3d(1,1,1:3,1)
            endif
            if (.not.(domain%south_boundary) .and. .not.(domain%east_boundary)) then
                fsnow(1,Ny_HICAR) = domain%south_in_2d(1,1,1)
                Ds(1:3,1,Ny_HICAR) = domain%south_in_3d(1,1,1:3,1)
            endif
            if (.not.(domain%north_boundary) .and. .not.(domain%east_boundary)) then
                fsnow(Nx_HICAR,Ny_HICAR) = domain%east_in_2d(1,1,1)
                Ds(1:3,Nx_HICAR,Ny_HICAR) = domain%east_in_3d(1,1,1:3,1)
            endif

        
        endif
        
    end subroutine exch_FSM_state_vars
    
    
    subroutine exch_SNTRAN_Qsalt(domain, Qsalt_u, Qsalt_v)
        implicit none
        
        type(domain_t), intent(inout) :: domain
        real, dimension(Nx_HICAR,Ny_HICAR), intent(inout) :: Qsalt_u, Qsalt_v

        if (.not.(domain%south_boundary)) then
            domain%south_buffer_2d(1,1:Ny_HICAR,1) = Qsalt_v(2,1:Ny_HICAR)
            !       !DIR$ PGAS DEFER_SYNC
            domain%north_in_2d(:,:,:)[domain%south_neighbor] = domain%south_buffer_2d(:,:,:)
        endif
        if (.not.(domain%north_boundary)) then
            domain%north_buffer_2d(1,1:Ny_HICAR,1) = Qsalt_v(Nx_HICAR-1,1:Ny_HICAR)
            !       !DIR$ PGAS DEFER_SYNC
            domain%south_in_2d(:,:,:)[domain%north_neighbor] = domain%north_buffer_2d(:,:,:)
        endif
        
        if (.not.(domain%east_boundary)) then
            domain%east_buffer_2d(1,1,1:Nx_HICAR) = Qsalt_u(1:Nx_HICAR,Ny_HICAR-1)
            !       !DIR$ PGAS DEFER_SYNC
            domain%west_in_2d(:,:,:)[domain%east_neighbor] = domain%east_buffer_2d(:,:,:)
        endif
        if (.not.(domain%west_boundary)) then
            domain%west_buffer_2d(1,1,1:Nx_HICAR) = Qsalt_u(1:Nx_HICAR,2)
            !       !DIR$ PGAS DEFER_SYNC
            domain%east_in_2d(:,:,:)[domain%west_neighbor] = domain%west_buffer_2d(:,:,:)
        endif

        sync images ( domain%neighbors )
        
        if (.not.(domain%south_boundary)) Qsalt_v(1,2:Ny_HICAR-1) = domain%south_in_2d(1,2:Ny_HICAR-1,1)
        if (.not.(domain%north_boundary)) Qsalt_v(Nx_HICAR,2:Ny_HICAR-1) = domain%north_in_2d(1,2:Ny_HICAR-1,1)
        if (.not.(domain%east_boundary)) Qsalt_u(2:Nx_HICAR-1,Ny_HICAR) = domain%east_in_2d(1,1,2:Nx_HICAR-1)
        if (.not.(domain%west_boundary)) Qsalt_u(2:Nx_HICAR-1,1) = domain%west_in_2d(1,1,2:Nx_HICAR-1)
        
        sync images ( domain%neighbors )

    end subroutine exch_SNTRAN_Qsalt

    subroutine exch_SLIDE_buffers(domain,SD_0,Sice_0)
        implicit none
        
        type(domain_t), intent(inout) :: domain
        
        real, dimension(Nx_HICAR,Ny_HICAR), intent(inout) :: SD_0, Sice_0

        if (.not.(domain%south_boundary)) then
            domain%south_buffer_2d(1,2:Ny_HICAR-1,1) = SD_0(2,2:Ny_HICAR-1)
            domain%south_buffer_2d(2,2:Ny_HICAR-1,1) = Sice_0(2,2:Ny_HICAR-1)
            !DIR$ PGAS DEFER_SYNC
            domain%north_in_2d(1:2,2:Ny_HICAR-1,1)[domain%south_neighbor] = domain%south_buffer_2d(1:2,2:Ny_HICAR-1,1)

            if (.not.(domain%east_boundary)) then
                !DIR$ PGAS DEFER_SYNC
                domain%north_in_2d(1:2,1,1)[domain%southeast_neighbor] = domain%south_buffer_2d(1:2,Ny_HICAR-1,1)
            endif
        endif
        if (.not.(domain%north_boundary)) then
            domain%north_buffer_2d(1,2:Ny_HICAR-1,1) = SD_0(Nx_HICAR-1,2:Ny_HICAR-1)
            domain%north_buffer_2d(2,2:Ny_HICAR-1,1) = Sice_0(Nx_HICAR-1,2:Ny_HICAR-1)
            !DIR$ PGAS DEFER_SYNC
            domain%south_in_2d(1:2,2:Ny_HICAR-1,1)[domain%north_neighbor] = domain%north_buffer_2d(1:2,2:Ny_HICAR-1,1)
            
            if (.not.(domain%west_boundary)) then
                !DIR$ PGAS DEFER_SYNC
                domain%south_in_2d(1:2,1,1)[domain%northwest_neighbor] = domain%north_buffer_2d(1:2,2,1)
            endif
        endif
        
        if (.not.(domain%east_boundary)) then
            domain%east_buffer_2d(1,1,2:Nx_HICAR-1) = SD_0(2:Nx_HICAR-1,Ny_HICAR-1)
            domain%east_buffer_2d(2,1,2:Nx_HICAR-1) = Sice_0(2:Nx_HICAR-1,Ny_HICAR-1)
            !DIR$ PGAS DEFER_SYNC
            domain%west_in_2d(1:2,1,2:Nx_HICAR-1)[domain%east_neighbor] = domain%east_buffer_2d(1:2,1,2:Nx_HICAR-1)

            if (.not.(domain%north_boundary)) then
                !DIR$ PGAS DEFER_SYNC
                domain%west_in_2d(1:2,1,1)[domain%northeast_neighbor] = domain%east_buffer_2d(1:2,1,Nx_HICAR-1)
            endif
        endif
        if (.not.(domain%west_boundary)) then
            domain%west_buffer_2d(1,1,2:Nx_HICAR-1) = SD_0(2:Nx_HICAR-1,2)
            domain%west_buffer_2d(2,1,2:Nx_HICAR-1) = Sice_0(2:Nx_HICAR-1,2)
            !DIR$ PGAS DEFER_SYNC
            domain%east_in_2d(1:2,1,2:Nx_HICAR-1)[domain%west_neighbor] = domain%west_buffer_2d(1:2,1,2:Nx_HICAR-1)

            if (.not.(domain%south_boundary)) then
                !DIR$ PGAS DEFER_SYNC
                domain%east_in_2d(1:2,1,1)[domain%southwest_neighbor] = domain%west_buffer_2d(1:2,1,2)
            endif
        endif

        sync images ( domain%neighbors )
        sync images ( domain%corner_neighbors )

        if (.not.(domain%south_boundary)) then
            SD_0(1,2:Ny_HICAR-1) = domain%south_in_2d(1,2:Ny_HICAR-1,1)
            Sice_0(1,2:Ny_HICAR-1) = domain%south_in_2d(2,2:Ny_HICAR-1,1)
            if (.not.(domain%east_boundary)) then
                SD_0(1,Ny_HICAR) = domain%south_in_2d(1,1,1)
                Sice_0(1,Ny_HICAR) = domain%south_in_2d(2,1,1)
            endif
        endif
        if (.not.(domain%north_boundary)) then
            SD_0(Nx_HICAR,2:Ny_HICAR-1) = domain%north_in_2d(1,2:Ny_HICAR-1,1)
            Sice_0(Nx_HICAR,2:Ny_HICAR-1) = domain%north_in_2d(2,2:Ny_HICAR-1,1)
            if (.not.(domain%west_boundary)) then
                SD_0(Nx_HICAR,1) = domain%north_in_2d(1,1,1)
                Sice_0(Nx_HICAR,1) = domain%north_in_2d(2,1,1)
            endif
        endif
        if (.not.(domain%east_boundary)) then
            SD_0(2:Nx_HICAR-1,Ny_HICAR) = domain%east_in_2d(1,1,2:Nx_HICAR-1)
            Sice_0(2:Nx_HICAR-1,Ny_HICAR) = domain%east_in_2d(2,1,2:Nx_HICAR-1)
            if (.not.(domain%north_boundary)) then
                SD_0(Nx_HICAR,Ny_HICAR) = domain%east_in_2d(1,1,1)
                Sice_0(Nx_HICAR,Ny_HICAR) = domain%east_in_2d(2,1,1)
            endif
        endif
        if (.not.(domain%west_boundary)) then
            SD_0(2:Nx_HICAR-1,1) = domain%west_in_2d(1,1,2:Nx_HICAR-1)
            Sice_0(2:Nx_HICAR-1,1) = domain%west_in_2d(2,1,2:Nx_HICAR-1)
            if (.not.(domain%south_boundary)) then
                SD_0(1,1) = domain%west_in_2d(1,1,1)
                Sice_0(1,1) = domain%west_in_2d(2,1,1)
            endif
        endif


    end subroutine exch_SLIDE_buffers


!!
end module module_sf_FSMdrv
