module array_utilities

    implicit none

    interface smooth_array
        module procedure smooth_array_2d, smooth_array_3d
    end interface

contains
    !>----------------------------------------------------------
    !! Generate an array with n_values elements spanning the range min to max
    !!
    !! This is similar to the matlab/numpy/etc linspace function
    !!
    !!----------------------------------------------------------
    subroutine linear_space(input_array, min_value, max_value, n_values)
        implicit none
        real,   intent(inout), allocatable :: input_array(:)
        real,   intent(in)  :: min_value, max_value
        integer,intent(in)  :: n_values

        integer :: i

        if (allocated(input_array)) then
            if (size(input_array)/=n_values) then
                deallocate(input_array)
            endif
        endif
        ! note, this can't be an else statement because the above block might change the result
        if (.not.allocated(input_array)) then
            allocate(input_array(n_values))
        endif

        do i=1,n_values
            input_array(i) = (i-1.0)/real(n_values-1.0) * (max_value - min_value) + min_value
        enddo

    end subroutine linear_space

    pure function check_array_dims(input_array, d1, d2, d3, d4, d5) result(passed)
        implicit none
        real,    intent(in) :: input_array(:,:,:,:,:)
        integer, intent(in) :: d1, d2, d3, d4, d5
        logical :: passed

        passed = .True.

        if (size(input_array,1)/=d1) then
            passed = .False.
            return
        endif

        if (size(input_array,2)/=d2) then
            passed = .False.
            return
        endif

        if (size(input_array,3)/=d3) then
            passed = .False.
            return
        endif

        if (size(input_array,4)/=d4) then
            passed = .False.
            return
        endif

        if (size(input_array,5)/=d5) then
            passed = .False.
            return
        endif

    end function check_array_dims

    !>----------------------------------------------------------
    !! Calculate the weights between the positions bestpos and nextpos
    !! based on the distance between match and indata(nextpos) (normalized by nextpos - bestpos)
    !! assumes indata is monotonically increasing,
    !! bestpos must be set prior to entry
    !! nextpos is calculated internally (either 1, bestpos+1, or n)
    !!
    !!----------------------------------------------------------
    function calc_weight(indata, bestpos, nextpos, match) result(weight)
        implicit none
        real :: weight
        real, dimension(:), intent(in) :: indata
        integer, intent(in) :: bestpos
        integer, intent(inout) :: nextpos
        real, intent(in) :: match

        integer :: n

        n=size(indata)

        if (match<indata(1)) then
            nextpos=1
            weight=1
        else
            if (bestpos==n) then
                nextpos=n
                weight=1
            else
                nextpos=bestpos+1
                weight=(indata(nextpos)-match) / (indata(nextpos) - indata(bestpos))
            endif
        endif

    end function

    !>------------------------------------------------------------
    !! Smooth an array (written for wind but will work for anything)
    !!
    !! Only smooths over the first (x) and second (y or z) or third (y or z) dimension
    !! ydim can be specified to allow working with (x,y,z) data or (x,z,y) data
    !! WARNING: this is a moderately complex setup to be efficient for the ydim=3 (typically large arrays, SLOW) case
    !! be careful when editing.
    !! For the complex case it pre-computes the sum of all columns for a given row,
    !! then to move from one column to the next it just has add the next column from the sums and subtracts the last one
    !! similarly, moving to the next row just means adding the next row to the sums, and subtracting the last one.
    !! Each point also has to be divided by N, but this decreases the compution from O(windowsize^2) to O(constant)
    !! Where O(constant) = 2 additions, 2 subtractions, and 1 divide regardless of windowsize!
    !!
    !! @param wind          3D array to be smoothed
    !! @param windowsize    size to smooth in both directions (i.e. the full window is this * 2 + 1)
    !! @param ysim          axis the y dimension is on (2 or 3) in the wind array
    !!
    !!------------------------------------------------------------
    subroutine smooth_array_3d(wind,windowsize,ydim)
        implicit none
        real, intent(inout), dimension(:,:,:):: wind    !> 3 dimensional wind field to be smoothed
        integer,intent(in)::windowsize                  !> halfwidth-1/2 of window to smooth over
                                                        ! Specified in grid cells, (+/- windowsize)
        integer,intent(in)::ydim                        !> the dimension to use for the y coordinate
                                                        ! It can be 2, or 3 (but not 1)
        real,allocatable,dimension(:,:,:)::inputwind    !> temporary array to store the input data in
        integer::i,j,k,nx,ny,nz,startx,endx,starty,endy ! various array indices/bounds
        ! intermediate sums to speed up the computation
        real,allocatable,dimension(:) :: rowsums,rowmeans
        real :: cursum
        integer :: cur_n,curcol,ncols,nrows

        ncols=windowsize*2+1
        nx=size(wind,1)
        ny=size(wind,2) !note, this is Z for the high-res domain (ydim=3)
        nz=size(wind,3) !note, this could be the Y or Z dimension depending on ydim

        allocate(inputwind(nx,ny,nz)) ! Can't be module level because nx,ny,nz could change between calls,

        inputwind=wind !make a copy so we always use the unsmoothed data when computing the smoothed data
        if (( ( (windowsize*2+1)>nz) .or. ((windowsize*2+1)>nx) ) .and. (ydim==3)) then
            write(*,*) "WARNING smoothing windowsize*2+1 is larger than nx or ny."
            write(*,*) "  This might lead to artifacts in the wind field especially near the borders."
            write(*,*) "  NX         = ", nx
            write(*,*) "  NY         = ", nz
            write(*,*) "  windowsize = ", windowsize
        endif

        !parallelize over a slower dimension (not the slowest because it is MUCH easier this way)
        ! as long as the inner loops (the array operations) are over the fastest dimension we are mostly OK
        ! $omp parallel firstprivate(windowsize,nx,ny,nz,ydim), &
        ! $omp private(i,j,k,startx,endx,starty,endy, rowsums,rowmeans,nrows,ncols,cursum), &
        ! $omp shared(wind,inputwind)
        allocate(rowsums(nx)) !this is only used when ydim=3, so nz is really ny
        allocate(rowmeans(nx)) !this is only used when ydim=3, so nz is really ny
        nrows=windowsize*2+1
        ncols=windowsize*2+1

        ! $omp do schedule(static)
        do j=1,ny

            ! so we pre-compute the sum over rows for each column in the current window
            if (ydim==3) then
                rowsums=inputwind(1:nx,j,1)*(windowsize+1)
                do i=1, min(windowsize,nz)
                    rowsums=rowsums+inputwind(1:nx,j,i)
                enddo
                if (windowsize > nz) then
                    rowsums = rowsums + inputwind(1:nx,j,nz) * (windowsize-nz)
                endif
            endif
            ! don't parallelize over this loop because it is much more efficient to be able to assume
            ! that you ran the previous row in serial for the slow (ydim=3) case
            do k=1,nz ! note this is y for ydim=3
                ! ydim=3 for the main model grid which is large and takes a long time
                if (ydim==3) then
                    ! if we are pinned to the top edge
                    if ((k-windowsize)<=1) then
                        starty=1
                        endy  =min(nz, k+windowsize)
                        rowsums=rowsums-inputwind(1:nx,j,starty)+inputwind(1:nx,j,endy)
                    ! if we are pinned to the bottom edge
                    else if ((k+windowsize)>nz) then
                        starty=max(2,k-windowsize)
                        endy  =nz
                        rowsums=rowsums-inputwind(1:nx,j,starty-1)+inputwind(1:nx,j,endy)
                    ! if we are in the middle (this is the most common)
                    else
                        starty=max(2,k-windowsize)
                        endy  =min(nz,k+windowsize)
                        rowsums=rowsums-inputwind(1:nx,j,starty-1)+inputwind(:,j,endy)
                    endif
                    rowmeans=rowsums/nrows
                    cursum=sum(rowmeans(1:windowsize))+rowmeans(1)*(windowsize+1)
                endif

                do i=1,nx
                    if (ydim==3) then
                        ! if we are pinned to the left edge
                        if ((i-windowsize)<=1) then
                            startx=1
                            endx  =min(nx,i+windowsize)
                            cursum=cursum-rowmeans(startx)+rowmeans(endx)
                        ! if we are pinned to the right edge
                        else if ((i+windowsize)>nx) then
                            startx=max(2,i-windowsize)
                            endx  =nx
                            cursum=cursum-rowmeans(startx-1)+rowmeans(endx)
                        ! if we are in the middle (this is the most common)
                        else
                            startx=max(2,i-windowsize)
                            endx  =min(nx,i+windowsize)
                            cursum=cursum-rowmeans(startx-1)+rowmeans(endx)
                        endif

                        wind(i,j,k)=cursum/ncols

                    else ! ydim==2
                        ! ydim=2 for the input data which is a small grid, thus cheap so we still use the slow method
                        !first find the current window bounds
                        startx=max(1, i-windowsize)
                        endx  =min(nx,i+windowsize)
                        starty=max(1, j-windowsize)
                        endy  =min(ny,j+windowsize)
                        ! then compute the mean within that window (sum/n)
                        ! note, artifacts near the borders (mentioned in ydim==3) don't affect this
                        ! because the borders *should* be well away from the domain
                        ! if the forcing data are not much larger than the model domain this *could* create issues
                        wind(i,j,k)=sum(inputwind(startx:endx,starty:endy,k)) &
                                    / ((endx-startx+1)*(endy-starty+1))
                    endif
                enddo
            enddo
        enddo
        ! $omp end do
        deallocate(rowmeans,rowsums)
        ! $omp end parallel

        deallocate(inputwind)
    end subroutine smooth_array_3d

    subroutine smooth_array_2d(input,windowsize)
        implicit none
        real, intent(inout), dimension(:,:):: input     !> 2 dimensional input field to be smoothed
        integer,intent(in)::windowsize                  !> halfwidth-1/2 of window to smooth over
                                                        ! Specified in grid cells, (+/- windowsize)
        real,allocatable,dimension(:,:)::inputtemp      !> temporary array to store the input data in
        integer::i,j,nx,ny,startx,endx,starty,endy      ! various array indices/bounds
        ! intermediate sums to speed up the computation
        real,allocatable,dimension(:) :: rowsums,rowmeans
        real :: cursum
        integer :: cur_n,curcol,ncols,nrows

        nx = size(input,1)
        ny = size(input,2)

        allocate(inputtemp(nx,ny))

        inputtemp = input !make a copy so we always use the unsmoothed data when computing the smoothed data
        if ((windowsize*2+1)>nx) then
            write(*,*) "WARNING can not operate if windowsize*2+1 is larger than nx"
            write(*,*) "NX         = ", nx
            write(*,*) "windowsize = ", windowsize
            stop
        endif

        !parallelize over a slower dimension (not the slowest because it is MUCH easier this way)
        ! as long as the inner loops (the array operations) are over the fastest dimension we are mostly OK
        !$omp parallel firstprivate(windowsize,nx,ny), &
        !$omp private(i,j,startx,endx,starty,endy, rowsums,rowmeans,nrows,ncols,cursum), &
        !$omp shared(input,inputtemp)
        allocate(rowsums(nx))
        allocate(rowmeans(nx))
        nrows = windowsize * 2 + 1
        ncols = windowsize * 2 + 1
        !$omp do schedule(static)
        do j=1,ny

            ! so we pre-compute the sum over rows for each column in the current window
            rowsums = 0
            do i = -1*windowsize, windowsize, 1
                starty = max(1,min(ny,j+i))
                rowsums = rowsums + inputtemp(1:nx,starty)
            enddo

            ! for a parallel algorithm (in 2D) we recompute the row sums for every y
            rowmeans = rowsums / nrows
            cursum = sum(rowmeans(1:windowsize)) + rowmeans(1) * (windowsize + 1)

            do i=1,nx
                ! if we are pinned to the left edge
                if ((i - windowsize)<=1) then
                    startx = 1
                    endx   = i + windowsize
                    cursum = cursum - rowmeans(startx) + rowmeans(endx)
                ! if we are pinned to the right edge
                else if ((i + windowsize) > nx) then
                    startx = i - windowsize
                    endx   = nx
                    cursum = cursum - rowmeans(startx-1) + rowmeans(endx)
                ! if we are in the middle (this is the most common)
                else
                    startx = i - windowsize
                    endx   = i + windowsize
                    cursum = cursum - rowmeans(startx-1) + rowmeans(endx)
                endif

                input(i,j) = cursum / ncols

            enddo
        enddo
        !$omp end do
        deallocate(rowmeans,rowsums)
        !$omp end parallel

        deallocate(inputtemp)
    end subroutine smooth_array_2d


end module array_utilities
