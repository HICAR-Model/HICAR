!>------------------------------------------------------------
!! Module to pre-compute and apply near-surface corrections to 
!! the domain wind field. This includes topographic 
!! exposure/sheltering via the Sx parameter (Winstral et al. 2017) 
!! and surface roughness.
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!------------------------------------------------------------
module wind_surf
    use data_structures
    use domain_interface,  only : domain_t
    use options_interface, only : options_t
    use io_routines,       only : io_write

    implicit none
    private
    public :: calc_Sx, calc_TPI, apply_Sx
    real, parameter :: deg2rad=0.017453293 !2*pi/360
    real, parameter :: rad2deg=57.2957779371
    !Values coming from COSMO-2 downscaling regression analysis per Winstral et al. 2017
    real, parameter, dimension(3,3) :: beta = reshape([-0.464, -0.236, 0.229, 0.155, 0.413, -0.055, 0.033, 0.031, 0.0],[3,3])
                                              
    
    
contains

    subroutine calc_TPI(domain)
        implicit none
        class(domain_t), intent(inout) :: domain
        
        real, allocatable    :: dist(:,:)
        integer           :: d_max, search_max, i, j, i_s, j_s, i_start_buffer, i_end_buffer, j_start_buffer, j_end_buffer
        integer           :: TPI_num
        real              :: search_height, TPI_sum
        
        !Use 2km per Winstral et al. 2017 paper
        d_max = 1000
        search_max = floor(d_max/domain%dx)
        
        allocate(domain%TPI(domain%grid2d% ims : domain%grid2d% ime, domain%grid2d% jms : domain%grid2d% jme))
        allocate(dist( 2*search_max+1, 2*search_max+1 ))
        
        do i = 1, 2*search_max+1
            do j = 1, 2*search_max+1
                dist(i,j) = sqrt(abs(i-(search_max+1.0))**2 + abs(j-(search_max+1.0))**2)
            end do
        end do
        
        !Convert distances to meters
        dist = dist*domain%dx
        
        !Now calc TPI
        do i=domain%grid2d%ims, domain%grid2d%ime
            do j=domain%grid2d%jms, domain%grid2d%jme
                TPI_num = 0
                TPI_sum = 0
                    
                ! Check to use buffers to avoid searching out of grid
                i_start_buffer = -min(0,i-(search_max+1))
                i_end_buffer = min(0,domain%grid2d%ide-(i+search_max))
                
                j_start_buffer = -min(0,j-(search_max+1))
                j_end_buffer = min(0,domain%grid2d%jde-(j+search_max))
                
                do i_s = 1+i_start_buffer, (search_max*2+1)+i_end_buffer
                    do j_s = 1+j_start_buffer, (search_max*2+1)+j_end_buffer 
                        if (dist(i_s,j_s) <= d_max .and. .not.(dist(i_s,j_s) == 0) ) then
                            
                            search_height = domain%global_terrain(i+(i_s-(search_max+1)),j+(j_s-(search_max+1)))
                        
                            TPI_sum = TPI_sum + search_height
                            TPI_num = TPI_num + 1
                        end if                               
                    end do
                end do
                domain%TPI(i,j) = domain%global_terrain(i,j) - TPI_sum/TPI_num
            end do
        end do
        
        
        !if ( this_image() == 1 ) then
        !    write (*,*) "Saving *_Sx.nc"
        !    !Save file
        !    call io_write(filename, "Sx", Sx(:,:,:,:) ) 
        !    call io_write("TPI_out.nc", "TPI", domain%TPI(:,:) ) 
        !endif
             
    end subroutine calc_TPI
    
    subroutine calc_Sx(domain, options, filename)
        implicit none
        class(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        character(len=*),   intent(in) :: filename
        
        real, allocatable    :: dist(:,:), azm(:,:), Sx_array_temp(:,:,:,:)
        integer, allocatable :: azm_indices(:,:), valid_ks(:)
        integer           :: d_max, search_max, i, j, k, ang, i_s, j_s, i_start_buffer, i_end_buffer, j_start_buffer, j_end_buffer
        integer           :: rear_ang, fore_ang, test_ang, rear_ang_diff, fore_ang_diff, ang_diff, k_max, window_rear, window_fore
        real              :: search_height, pt_height, h_diff, Sx_temp
        

        d_max = options%wind%Sx_dmax
        k_max = 30 !Max number of layers to compute Sx for
        search_max = floor(d_max/domain%dx)
        
        ! Initialize Sx and set to lowest possible (yet, un-physical) value of -90
        allocate(domain%Sx( 1:72, domain%grid2d%ims:domain%grid2d%ime, 1:k_max, domain%grid2d%jms:domain%grid2d%jme ))
        allocate(Sx_array_temp( 1:72, domain%grid2d%ims:domain%grid2d%ime, 1:k_max, domain%grid2d%jms:domain%grid2d%jme ))
        Sx_array_temp = -90.0
        
        !Pre-compute distances, using search_max to build grid
        allocate(dist( 2*search_max+1, 2*search_max+1 ))
        allocate(azm( 2*search_max+1, 2*search_max+1 ))
        allocate(azm_indices( 2*search_max+1, 2*search_max+1 ))
        azm = 0
                
        do i = 1, 2*search_max+1
            do j = 1, 2*search_max+1
                dist(i,j) = sqrt(abs(i-(search_max+1.0))**2 + abs(j-(search_max+1.0))**2)
                
                if (dist(i,j) > 0) azm(i,j)  = atan2(1.0*(i-(search_max+1)),1.0*(j-(search_max+1)))

            end do
        end do
        
        !Convert distances to meters
        dist = dist*domain%dx
        !convert azm to deg
        azm = azm*rad2deg
        where(azm < 0) azm = 360+azm
        where(azm >= 360.0) azm=0.0
        azm_indices = int(azm/5)+1
                
        !call unique_sort_ind(azm_indices,valid_ks)
        
        do i=domain%grid2d%ims, domain%grid2d%ime
            do j=domain%grid2d%jms, domain%grid2d%jme
                do k = 1, k_max

                    if (k == 1) then
                        pt_height = domain%global_terrain(i,j)
                    else if (k > 1) then
                        pt_height = pt_height + domain%global_dz_interface(i,k,j)
                    end if
                    
                    ! Check to use buffers to avoid searching out of grid
                    i_start_buffer = -min(0,i-(search_max+1))
                    i_end_buffer = min(0,domain%grid2d%ide-(i+search_max))
                
                    j_start_buffer = -min(0,j-(search_max+1))
                    j_end_buffer = min(0,domain%grid2d%jde-(j+search_max))
                
                    do i_s = 1+i_start_buffer, (search_max*2+1)+i_end_buffer
                        do j_s = 1+j_start_buffer, (search_max*2+1)+j_end_buffer
                            !Since dist is for a grid, and we want a search RADIUS, check that we are within d_max
                            if (dist(i_s,j_s) <= d_max .and. .not.(dist(i_s,j_s) == 0) ) then
                                !Calculate height difference
                                search_height = domain%global_terrain(i+(i_s-(search_max+1)),j+(j_s-(search_max+1)))
                                h_diff = search_height - pt_height

                                !Calculate Sx slope to search-cell
                                Sx_temp = atan(h_diff/dist(i_s,j_s))*rad2deg
                            
                                ! If new Sx is greater than existing Sx for a given search angle, replace
                                if (Sx_temp > Sx_array_temp(azm_indices(i_s,j_s),i,k,j)) then
                                    Sx_array_temp(azm_indices(i_s,j_s),i,k,j) = Sx_temp
                                end if
                            end if
                        end do
                    end do
                                    
                    !After finding Sx in each absolute direction around grid cell, 
                    !Pick max for each 30º window and perform interpolation to other directions if necesarry
                    
                    rear_ang = 1 
                    fore_ang = 1

                    do ang = 1, 72
                                    
                        !Determine indices for interpolation
                        if (ang==fore_ang) then
                            !Update indices for interpolated Sx's
                            rear_ang = ang
                            test_ang = ang+1
                            do while (Sx_array_temp(test_ang,i,k,j) <= -90.0)
                                test_ang = test_ang+1
                                if (test_ang > 72) exit 
                            end do
                            fore_ang = test_ang
                         
                            !Handle wrap-around case
                            if (fore_ang > 72) fore_ang = 1

                        end if
                        
                        !Perform 30º window max search
                        window_rear = ang-6
                        window_fore = ang+6
                        
                        if (ang <= 6) then
                            window_rear = 72-(6-ang)
                            domain%Sx(ang,i,k,j) = max(maxval(Sx_array_temp(window_rear:72,i,k,j)), &
                                                       maxval(Sx_array_temp(1:window_fore,i,k,j)))
                        else if (ang >= 67) then
                            window_fore = 6-(72-ang)
                            domain%Sx(ang,i,k,j) = max(maxval(Sx_array_temp(window_rear:72,i,k,j)), &
                                                       maxval(Sx_array_temp(1:window_fore,i,k,j)))
                        else
                            domain%Sx(ang,i,k,j) = maxval(Sx_array_temp(window_rear:window_fore,i,k,j))
                        end if
                    
                        !If we did not calculate Sx for a given direction
                        if (domain%Sx(ang,i,k,j) == -90.0) then
                            !Weight the two surrounding Sx values based on our angular-distance to them
                            rear_ang_diff = ang-rear_ang
                            fore_ang_diff = fore_ang-ang
                            ang_diff = fore_ang-rear_ang
                        
                            !Handle wrap-around case
                            if (ang > fore_ang) then
                                fore_ang_diff = fore_ang+(72-ang)
                                ang_diff = fore_ang+(72-rear_ang)
                            end if
                        
                            !Interpolation, linearly-weighted by angular-distance from values
                            domain%Sx(ang,i,k,j) = (Sx_array_temp(rear_ang,i,k,j)*fore_ang_diff + &
                                                    Sx_array_temp(fore_ang,i,k,j)*rear_ang_diff)/ang_diff

                        end if
                    end do
                end do
            end do
        end do
        
    
        !if ( this_image() == 1 ) then
        !    write (*,*) "Saving *_Sx.nc"
        !    !Save file
        !    call io_write(filename, "Sx", domain%Sx(:,:,:,:) ) 
        !    call io_write("TPI_out.nc", "TPI", domain%TPI(:,:) ) 
        !endif
        
        deallocate(dist)
        deallocate(Sx_array_temp)
        deallocate(azm)
        deallocate(azm_indices)
        
        !Sync images befoer exiting so that we don't try to read Sx while it is being written
        !sync all

    end subroutine calc_Sx
    
    subroutine unique_sort_ind(arr,us_arr)
        implicit none
        integer, allocatable, intent(in)    :: arr(:,:)
        integer, allocatable, intent(inout) :: us_arr(:)
        real, allocatable                :: unique(:)
        integer   :: i
        real      :: min_arr, max_arr
        
        allocate(unique(1:size(arr)))
        
        i = 0
        min_arr = minval(arr)-1
        max_arr = maxval(arr)
        
        do while (min_arr < max_arr)
            i = i+1
            min_arr = minval(arr, mask=(arr>min_arr))
            unique(i) = min_arr
        end do
        
        if (allocated(us_arr)) deallocate(us_arr)
        allocate(us_arr(1:i))
        us_arr = unique(1:i)
        
        deallocate(unique)
       
        
    end subroutine unique_sort_ind
    
    
    
    subroutine apply_Sx(Sx, TPI, u, v, w)
        implicit none
        real, intent(in)                       :: Sx(:,:,:,:), TPI(:,:)
        real, intent(inout),  dimension(:,:,:) :: u, v, w
        
        real, allocatable, dimension(:,:)   :: winddir, Sx_curr
        real, allocatable, dimension(:,:,:)   :: Sx_U, Sx_V

        integer ::  i, j, k, ims, ime, jms, jme, CV_u, CV_v
        real    :: TPI_corr, Sx_U_corr, Sx_V_corr
                
        ims = lbound(w,1)
        ime = ubound(w,1)
        jms = lbound(w,3)
        jme = ubound(w,3)
                
        allocate(Sx_curr(ims:ime,jms:jme))
        
        allocate(Sx_U(ims:ime+1,1:30,jms:jme))
        allocate(Sx_V(ims:ime,1:30,jms:jme+1))

        !Initialize staggered-Sx. This will keep the border values from being updated
        Sx_U = 0
        Sx_V = 0
        
        CV_u = 0
        CV_v = 0

        !Pick appropriate Sx for wind direction and interpolate Sx to staggered grids
        do k = 1, 30
            call pick_Sx(Sx(:,:,k,:), Sx_curr, u(:,k,:), v(:,k,:))
            Sx_U(ims+1:ime,k,:) = (Sx_curr(ims:ime-1,:) + Sx_curr(ims+1:ime,:))/2
            Sx_V(:,k,jms+1:jme) = (Sx_curr(:,jms:jme-1) + Sx_curr(:,jms+1:jme))/2
        end do

        !Loop through i,j
        do i = ims+1, ime
            do j = jms+1, jme
            
                !Loop through vertical column
                do k = 1, 30
                
                
                    !if ( (k == 1) .and. (TPI(i,j) > 100.0)) then !ridges
                        !u(i,k,j) = u(i,k,j) - u(i,k,j)*beta(1,1) + Sx_U(i,1,j)*beta(1,3) !+ CV_u*beta(1,2) 
                        !v(i,k,j) = v(i,k,j) - v(i,k,j)*beta(1,1) + Sx_V(i,1,j)*beta(1,3) !+ CV_v*beta(1,2) 
                    !else if (TPI(i,j) <= -100.0) then !valleys
                        !u(i,k,j) = u(i,k,j) - u(i,k,j)*beta(3,1) + Sx_U(i,1,j)*beta(3,3) !+ CV_u*beta(3,2) 
                        !v(i,k,j) = v(i,k,j) - v(i,k,j)*beta(3,1) + Sx_V(i,1,j)*beta(3,3) !+ CV_v*beta(3,2) 
                    !else !upper slopes
                        !u(i,k,j) = u(i,k,j) - u(i,k,j)*beta(2,1) + Sx_U(i,1,j)*beta(2,3) !+ CV_u*beta(2,2) 
                        !v(i,k,j) = v(i,k,j) - v(i,k,j)*beta(2,1) + Sx_V(i,1,j)*beta(2,3) !+ CV_v*beta(2,2) 
                    !end if
                
                    !Consider sheltering corrections
                    !If surface was sheltered and we are still sheltered: U-grid
                    if ((Sx_U(i,1,j) > 0) .and. (Sx_U(i,k,j) > 0)) then
                        !Sheltered correction
                        Sx_U_corr = Sx_U(i,k,j) / (Sx_U(i,k,j)+5)
                        !w(i,k,j) = -5
                    else 
                        Sx_U_corr = 0
                    end if
                    
                    !If surface was sheltered and we are still sheltered: V-grid
                    if ((Sx_V(i,1,j) > 0) .and. (Sx_V(i,k,j) > 0)) then
                        !Sheltered correction
                        Sx_V_corr = Sx_V(i,k,j) / (Sx_V(i,k,j)+5)
                        !w(i,k,j) = -5
                    else 
                        Sx_V_corr = 0
                    end if

                
                    !Compute TPI_corr
                    if (k <= 10) then
                        !Scale TPI correction and exposure with height so we smoothly merge with forcing data
                        TPI_corr = (TPI(i,j)/200) * ((11-k)/10.0) !This is setup such that ridges are >= 100, valleys <= -100
                        if (TPI_corr > 0.5) TPI_corr = 0.5
                        if (TPI_corr < -0.5) TPI_corr = -0.5
                        
                        !Consider exposure, only for first 10-levels
                        ! If Sx is negative, compute exposure correction
                        if (Sx_U(i,1,j) < 0) then
                            Sx_U_corr = (Sx_U(i,1,j)/30) * ((11-k)/10.0) !Exposure used to be /30
                            !w(i,k,j) = 5 * ((11-k)/10.0)
                        end if
                        
                        ! If Sx is negative, compute exposure correction
                        if (Sx_V(i,1,j) < 0) then
                            Sx_V_corr = (Sx_V(i,1,j)/30) * ((11-k)/10.0) !Exposure used to be /30
                            !w(i,k,j) = 5 * ((11-k)/10.0)
                        end if
                    else 
                        TPI_corr = 0
                    end if
                    
                       
                    !If we have no correction to apply, we are done with this cell
                    if ( Sx_U_corr == 0 .and. Sx_V_corr == 0 .and. TPI_corr == 0 ) exit
                    
                    !Finally, apply TPI and Sx corrections
                    u(i,k,j) = u(i,k,j) * (1 + TPI_corr - Sx_U_corr)
                    v(i,k,j) = v(i,k,j) * (1 + TPI_corr - Sx_V_corr)
                end do
            end do
        end do
        
        
        
        
        
        

        
        !Make mass-centered u/v at ground-level (domain%u_m/v_m are not initialized on spin-up)
        !do k = 1,ubound(Sx,3)
        !    Sx_curr = pick_Sx(Sx(:,:,k,:), u(:,k,:), v(:,k,:))

        !    !Interpolate Sx grid to staggered U/V grid for wind-field modification
        !    Sx_U(ims+1:ime,:) = (Sx_curr(ims:ime-1,:) + Sx_curr(ims+1:ime,:))/2
        !    Sx_V(:,jms+1:jme) = (Sx_curr(:,jms:jme-1) + Sx_curr(:,jms+1:jme))/2
        
        !    !If Sx is >= 20º, apply reduction of 50% U/V fields
        !    where(Sx_U >= 10.0) u(:,k,:) = u(:,k,:)*0.5
        !    where(Sx_V >= 10.0) v(:,k,:) = v(:,k,:)*0.5
        !end do
        
    end subroutine apply_Sx
    
    subroutine pick_Sx(Sx,Sx_curr,u,v)
        implicit none
        real, intent(in)                     :: Sx(:,:,:)
        real, intent(inout)                  :: Sx_curr(:,:)
        real, intent(in),  dimension(:,:)    :: u, v
                        
        real, allocatable, dimension(:,:)   :: winddir, u_m, v_m
        integer, allocatable                ::  dir_indices(:,:)
        integer ::  i, j, ims, ime, jms, jme
        
        ims = lbound(v,1)
        ime = ubound(v,1)
        jms = lbound(u,2)
        jme = ubound(u,2)
                
        u_m = (u(ims:ime,:) + u(ims+1:ime+1,:))/2
        v_m = (v(:,jms:jme) + v(:,jms+1:jme+1))/2
        
        !Compute wind direction for each cell on mass grid
        winddir = atan2(-u_m,-v_m)*rad2deg
        where(winddir < 0.0) winddir = winddir+360
        where(winddir == 360.0) winddir = 0.0
        dir_indices = int(winddir/5)+1
                
        !Build grid of Sx values based on wind direction at that cell
        do i = ims, ime
            do j = jms, jme
                Sx_curr(i,j) = Sx(dir_indices(i,j),i,j)
            end do
        end do
        

        
    end subroutine pick_Sx

end module wind_surf
