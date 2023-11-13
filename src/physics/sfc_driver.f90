!>----------------------------------------------------------
!! This module provides a wrapper to call various PBL models
!! It sets up variables specific to the physics package to be used including both
!!
!! The main entry point to the code is pbl(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!!  pbl_init->[ external initialization routines]
!!  pbl->[  external PBL routines]
!!  pbl_finalize
!!
!! High level routine descriptions / purpose
!!   pbl_init           - initializes physics package
!!   pbl                - sets up and calls main physics package
!!   pbl_finalize       - permits physics package cleanup (close files, deallocate memory)
!!
!! Inputs: domain, options, dt
!!      domain,options  = as defined in data_structures
!!      dt              = time step (seconds)
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module surface_layer
    use data_structures
    use domain_interface,   only : domain_t
    use options_interface,  only : options_t
    use module_sf_sfclayrev, only : sfclayrevinit, SFCLAYREV
    use mod_wrf_constants,   only : KARMAN, gravity, cp, R_d, rcp, EP_1, EP_2, SVPT0, SVP1, SVP2, SVP3, EOMEG, STBOLT, p1000mb, XLV
    use ieee_arithmetic ! for debugging


    implicit none
    real,allocatable, dimension(:,:)    ::  windspd, gz1oz0, t2d, th2d, regime, flhc, flqc, q2d, &
                                            rmol, mol, qgh, qsfc, cpm, mavail, zol

    private
    public :: sfc_var_request, sfc_init, sfc   !, sfc_finalize

    integer :: ids, ide, jds, jde, kds, kde,  &
               ims, ime, jms, jme, kms, kme,  &
               its, ite, jts, jte, kts, kte, j, k, i

    logical :: allowed_to_read, restart, flag_qi

contains

    subroutine sfc_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%surfacelayer == kSFC_MM5REV) then
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%temperature, kVARS%potential_temperature, kVARS%surface_pressure, &
                         kVARS%dz_interface, kVARS%pressure,  kVARS%skin_temperature, &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,                  &
                         kVARS%roughness_z0, kVARS%hpbl, kVARS%QFX,       &
                         kVARS%land_mask, kVARS%br, kVARS%ustar,                      &
                         kVARS%chs, kVARS%chs2, kVARS%cqs2,                           &
                         kVARS%u, kVARS%v, kVARS%psim, kVARS%psih, kVARS%fm, kVARS%fh])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%temperature, kVARS%potential_temperature, kVARS%surface_pressure, &
                         kVARS%dz_interface, kVARS%pressure,  kVARS%skin_temperature, &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,                  &
                         kVARS%roughness_z0, kVARS%hpbl, kVARS%QFX,       &
                         kVARS%chs, kVARS%chs2, kVARS%cqs2,               &
                         kVARS%u, kVARS%v ])

        endif
    end subroutine sfc_var_request


    subroutine sfc_init(domain,options)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(options_t),    intent(in)      :: options

        ids = domain%ids ; ide = domain%ide ; jds = domain%jds ; jde = domain%jde ; kds = domain%kds ; kde = domain%kde
        ims = domain%ims ; ime = domain%ime ; jms = domain%jms ; jme = domain%jme ; kms = domain%kms ; kme = domain%kme
        its = domain%its ; ite = domain%ite ; jts = domain%jts ; jte = domain%jte ; kts = domain%kts ; kte = domain%kte

        if (this_image()==1) write(*,*) "Initializing Surface Layer scheme"

        if (options%physics%surfacelayer==kSFC_MM5REV) then
            if (this_image()==1) write(*,*) "    Revised MM5"
            
            allocate(windspd(ims:ime, jms:jme))
            allocate(rmol(ims:ime, jms:jme))
            allocate(mol(ims:ime, jms:jme))
            allocate(zol(ims:ime, jms:jme))
            allocate(qgh(ims:ime, jms:jme))
            allocate(qsfc(ims:ime, jms:jme))
            allocate(cpm(ims:ime, jms:jme))
            allocate(flhc(ims:ime, jms:jme))
            allocate(flqc(ims:ime, jms:jme))
            allocate(regime(ims:ime, jms:jme))
            allocate(mavail(ims:ime, jms:jme))
            rmol = 0.0
            mol = 0.0
            zol = 0.0
            qgh = 0.0
            qsfc = 0.0
            cpm = 0.0
            flhc = 0.0
            flqc = 0.0
            mavail = 0.0

            allocate(t2d(ims:ime, jms:jme))
            allocate(th2d(ims:ime, jms:jme))
            allocate(q2d(ims:ime, jms:jme))
            allocate(gz1oz0(ims:ime, jms:jme))  !-- gz1oz0      log(z/z0) where z0 is roughness length
            gz1oz0 = log((domain%z%data_3d(:,kts,:) - domain%terrain%data_2d) / domain%roughness_z0%data_2d)
            
            
            call sfclayrevinit(ims,ime,jms,jme,                    &
                            its,ite,jts,jte,                       &
                            options%sfc_options%sbrlim)!,          &
!                            bathymetry_flag, shalwater_z0,      &
!                            shalwater_depth, water_depth,       &
!                            xland,LakeModel,lake_depth,lakemask )       
        endif
    end subroutine sfc_init

    subroutine sfc(domain, options, dt_in)
        implicit none
        type(domain_t),  intent(inout)  :: domain
        type(options_t), intent(in)     :: options
        real,            intent(in)     :: dt_in  !  =real(dt%seconds())

        if (options%physics%surfacelayer==kSFC_MM5REV) then
        
            windspd = sqrt(domain%u_mass%data_3d(ims:ime,kms,jms:jme)**2 + domain%v_mass%data_3d(ims:ime,kms,jms:jme)**2)
            where(windspd==0) windspd=1e-5            

            call sfclayrev(u3d=domain%u_mass%data_3d   & !-- u3d         3d u-velocity interpolated to theta points (m/s)
               ,v3d=domain%v_mass%data_3d              & !-- v3d         3d v-velocity interpolated to theta points (m/s)
               ,t3d=domain%temperature%data_3d         &
               ,qv3d=domain%water_vapor%data_3d        &
               ,p3d=domain%pressure%data_3d            & !-- p3d         3d pressure (pa)
               ,dz8w=domain%dz_interface%data_3d       & !-- dz8w        dz between full levels (m)
               ,cp=cp                                  &
               ,g=gravity                              &
               ,rovcp=rcp                              & ! rovcp = Rd/cp
               ,r=R_d                                  &  ! J/(kg K) specific gas constant for dry air
               ,xlv=XLV                                & !-- xlv         latent heat of vaporization (j/kg)
               ,psfc=domain%surface_pressure%data_2d   &
               ,chs=domain%chs%data_2d                 &
               ,chs2=domain%chs2%data_2d               &
               ,cqs2=domain%cqs2%data_2d               &
               ,cpm=cpm                                &
               ,mavail=mavail                          &
               ,regime=regime                          &
               ,psim=domain%psim%data_2d               &
               ,psih=domain%psih%data_2d               &
               ,fm=domain%fm%data_2d                   &
               ,fh=domain%fh%data_2d                   &
               ,flhc=flhc                              & !these are only used by 1 PBL scheme in WRF, which is not in ICAR, so we dont need to use them
               ,flqc=flqc                              & !these are only used by 1 PBL scheme in WRF, which is not in ICAR, so we dont need to use them
               ,mol=mol                                &
               ,rmol=rmol                              &
               ,qgh=qgh                                &
               ,qsfc=qsfc                              &
               ,u10=domain%u_10m%data_2d               &
               ,v10=domain%v_10m%data_2d               &
               ,znt=domain%roughness_z0%data_2d        & ! i/o -- znt		roughness length (m) (input only)
               ,ust=domain%ustar                       & ! i/o -- ust		u* in similarity theory (m/s)
               ,zol=zol                                & ! i/o -- zol		z/l height over monin-obukhov length - intent(inout) - but appears to not be used really?
               ,pblh=domain%hpbl%data_2d               & ! i/o -- hpbl	pbl height (m) - intent(inout)
               ,xland=real(domain%land_mask)           &
               ,hfx=domain%sensible_heat%data_2d       & !  HFX  - net upward heat flux at the surface (W/m^2)
               ,qfx=domain%qfx%data_2d                 & !  QFX  - net upward moisture flux at the surface (kg/m^2/s)
               ,lh=domain%latent_heat%data_2d          & !  LH  - net upward latent flux at the surface (W/m^2/s)
               ,th2=th2d                               &
               ,t2=t2d                                 &
               ,q2=q2d                                 &
               ,gz1oz0=gz1oz0                          &
               ,br=domain%br%data_2d                   &
               ,wspd=windspd                           & ! i/o -- wspd        wind speed at lowest model level (m/s)
               ,tsk=domain%skin_temperature%data_2d    &
               ,isfflx=options%sfc_options%isfflx      &
               ,dx=domain%dx                           &
               ,svp1=SVP1                              & !-- svp1        constant for saturation vapor pressure (kpa)
               ,svp2=SVP2                              & !-- svp2        constant for saturation vapor pressure (dimensionless)
               ,svp3=SVP3                              & !-- svp3        constant for saturation vapor pressure (k)
               ,svpt0=SVPT0                            & !-- svpt0       constant for saturation vapor pressure (k)
               ,ep1=EP_1                               & !-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
               ,ep2=EP_2                               & !-- ep2         constant for specific humidity calculation
               ,karman=KARMAN                          & !-- karman      von karman constant
               ,eomeg=EOMEG                            & 
               ,stbolt=STBOLT                          & 
               ,P1000mb=p1000mb                        &
               ,ids=ids, ide=ide, jds=jds, jde=jde     &
               ,kds=kds, kde=kde, ims=ims, ime=ime     &
               ,jms=jms, jme=jme, kms=kms, kme=kme     &
               ,its=its, ite=ite, jts=jts, jte=jte     &
               ,kts=kts, kte=kte                       &
               ,shalwater_z0=0                         &
               ,iz0tlnd=options%sfc_options%iz0tlnd    &
               ,isftcflx=options%sfc_options%isftcflx  &
               ,scm_force_flux=options%sfc_options%scm_force_flux)
                                             
                endif
    end subroutine sfc

!    subroutine sfc_finalize(options)
!        implicit none
!        type(options_t), intent(in) :: options
!
!        if (options%physics%surfacelayer==kSFC_MM5REV) then
!        endif
!
!    end subroutine sfc_finalize
end module surface_layer
