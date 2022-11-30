module time_io

    use data_structures
    use icar_constants
    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use string,             only : get_integer
    use io_routines,        only : io_read, io_read_attribute
    use iso_fortran_env, only: real64, real128

    implicit none

contains

    function find_timestep_in_file(filename, time_var, time, time_at_step, precision, error) result(step)
        implicit none
        character(len=*),  intent(in) :: filename
        character(len=*),  intent(in) :: time_var
        type(Time_type),   intent(in) :: time
        type(Time_type),   intent(inout), optional :: time_at_step
        type(time_delta_t),intent(in),    optional :: precision
        integer,           intent(inout), optional :: error
        integer :: step

        type(Time_type), allocatable :: times_in_file(:)
        type(time_delta_t) :: max_dt
        integer :: i,n
        logical :: found

        call max_dt%set(seconds=1.0)
        if (present(precision)) max_dt = precision
        if (present(error)) error=0

        ! read the times for all timesteps in the specified file
        call read_times(filename, time_var, times_in_file)

        step = -1
        found= .False.
        n    = size(times_in_file)
        ! loop through times looking for a time that matches the input time to within
        ! a specified maximum delta t
        do i = 1, n
            if (.not.found) then
                if (times_in_file(i)%equals(time, precision=max_dt)) then
                    step = i
                    found=.True.
                elseif (times_in_file(i) > time) then
                    step = i-1
                    found=.True.
                endif
            endif
        enddo

        if (step < 1) then
            if (present(error)) then
                error = 1
            else
                write(*,*) "ERROR: Unable to find requested date in file."
                write(*,*) "Filename: ",trim(filename)
                write(*,*) "  time  : ",trim(time%as_string())
                write(*,*) "First time in file : ", trim(times_in_file(1)%as_string())
                write(*,*) " Last time in file : ", trim(times_in_file(n)%as_string())
                stop "Unable to find date in file"
            endif
        endif

        if (present(time_at_step)) time_at_step = times_in_file(step)

        deallocate(times_in_file)

    end function find_timestep_in_file


    function time_gain_from_units(units) result(gain)
        implicit none
        character(len=*), intent(in) :: units
        real(real128) :: gain

        if ((units(1:4)=="days").or.(units(1:4)=="Days")) then
            gain = 1.0Q0
        else if ((units(1:4)=="hour").or.(units(1:4)=="Hour")) then
            gain = 24.0Q0
        else if ((units(1:3)=="min").or.(units(1:3)=="Min")) then
            gain = 1440.0Q0
        else if ((units(1:3)=="sec").or.(units(1:3)=="Sec")) then
            gain = 86400.0Q0
        else
            write(*,*) trim(units)
            stop "Error: unknown units"
        endif

    end function time_gain_from_units

    function year_from_units(units) result(year)
        implicit none
        character(len=*), intent(in) :: units
        integer :: year

        integer :: since_loc, year_loc

        since_loc = index(units,"since")

        year_loc = index(units(since_loc:)," ")
        year_loc = year_loc+since_loc

        year = get_integer(units(year_loc:year_loc+3))

    end function year_from_units

    function month_from_units(units) result(month)
        implicit none
        character(len=*), intent(in) :: units
        integer :: month

        integer :: since_loc, month_loc

        since_loc = index(units,"since")

        month_loc = index(units(since_loc:)," ") + 5
        month_loc = month_loc + since_loc

        month = get_integer(units(month_loc:month_loc+1))

    end function month_from_units

    function day_from_units(units) result(day)
        implicit none
        character(len=*), intent(in) :: units
        integer :: day

        integer :: since_loc, day_loc

        since_loc = index(units,"since")

        day_loc = index(units(since_loc:)," ") + 8
        day_loc = day_loc + since_loc

        day = get_integer(units(day_loc:day_loc+1))

    end function day_from_units

    function hour_from_units(units, error) result(hour)
        implicit none
        character(len=*), intent(in) :: units
        integer, intent(out), optional :: error
        integer :: hour

        integer :: since_loc, hour_loc
        if (present(error)) error = 0

        since_loc = index(units,"since")

        hour_loc = index(units(since_loc:)," ") + 11
        hour_loc = hour_loc + since_loc

        ! default return value if hours can't be read from the units attribute (e.g. they aren't present)
        hour = 0

        if( hour_loc+1 <= len(units) ) then
           if (trim(units(hour_loc:hour_loc+1)) /= "") then
               hour = get_integer(units(hour_loc:hour_loc+1))
           endif
        else
           if (present(error)) error = 1
        endif

    end function hour_from_units

    function find_timestep_in_filelist(filelist, time_var, time, filename, error) result(step)
        implicit none
        character(len=*),  intent(in) :: filelist(:)
        character(len=*),  intent(in) :: time_var
        type(Time_type),   intent(in) :: time
        character(len=*),  intent(inout), optional :: filename
        integer,           intent(inout), optional :: error
        integer :: step

        type(Time_type), allocatable :: times_in_file_up(:), times_in_file_down(:)
        type(time_delta_t) :: max_dt
        integer :: i, n_t_down, n_t_up, nup, ndown, nup_o, ndown_o
        logical :: found

        call max_dt%set(seconds=1.0)
        if (present(error)) error=0
        
        step = -1
        found = .False.
        
        nup = size(filelist)
        ndown=1
        ndown_o = ndown
        nup_o = nup

        call read_times(filelist(ndown), time_var, times_in_file_down)
        call read_times(filelist(nup), time_var, times_in_file_up)


        do while(.not.(found))
            ! read the times for all timesteps in the specified file
            if (.not.(ndown==ndown_o)) call read_times(filelist(ndown), time_var, times_in_file_down)

            n_t_down  = size(times_in_file_down)
            ! loop through times looking for a time that matches the input time to within
            ! a specified maximum delta t
            do i = 1, n_t_down
                if (.not.found) then
                    if (times_in_file_down(i)%equals(time, precision=max_dt)) then
                        step = i
                        found=.True.
                        filename = filelist(ndown)
                    elseif (i > 1) then
                        if ((times_in_file_down(i) > time) .and. (times_in_file_down(i-1) < time)) then
                            step = i-1
                            found=.True.
                            filename = filelist(ndown)
                        endif
                    endif
                endif
            enddo
                
            ! read the times for all timesteps in the specified file
            if (.not.(nup==nup_o)) call read_times(filelist(nup), time_var, times_in_file_up)

            n_t_up  = size(times_in_file_up)
            ! loop through times looking for a time that matches the input time to within
            ! a specified maximum delta t
            do i = 1, n_t_up
                if (.not.found) then
                    if (times_in_file_up(i)%equals(time, precision=max_dt)) then
                        step = i
                        found=.True.
                        filename = filelist(nup)
                    elseif (i > 1) then
                        if ((times_in_file_up(i) > time) .and. (times_in_file_up(i-1) < time)) then
                            step = i-1
                            found=.True.
                            filename = filelist(nup)
                        endif
                    endif
                endif
            enddo
            
            if (.not.found) then
                !Either it is bounded
                if (times_in_file_down(n_t_down) < time .and. times_in_file_up(1) > time) then
                    !If bounded, save these positions
                    ndown_o = ndown
                    nup_o = nup
                    !Advance ndown by half of distance between nup and ndown
                    ndown = ndown + nint((nup-ndown)/2.0)
                ! or we overshot (since on the bounded condition we advance up)
                elseif(times_in_file_down(1) > time .and. times_in_file_up(1) > time) then
                    !set upper bound to former lower bound
                    nup = ndown
                    times_in_file_up = times_in_file_down
                    nup_o = nup
                    !Set lower bound to last bounding lower bound
                    ndown = ndown_o
                    ndown_o = 0
                ! only other option is that the requested time is not included in the filelist
                else
                    if (present(error)) then
                        error = 1
                        found = .True.
                    else
                        write(*,*) "ERROR: Unable to find requested date in filelist."
                        write(*,*) "First filename: ",trim(filelist(1))
                        write(*,*) "  time  : ",trim(time%as_string())
                        stop "Unable to find date in file"
                    endif
                endif
            endif
            
            !Handle the limiting case
            if (ndown==nup) ndown=nup-1
        enddo
        
        deallocate(times_in_file_down)
        deallocate(times_in_file_up)

    end function find_timestep_in_filelist


    subroutine read_times(filename, varname, times, timezone_offset, curstep)
        implicit none
        character(len=*),   intent(in) :: filename, varname
        type(Time_type),    intent(inout), allocatable, dimension(:) :: times
        real(real128),      intent(in), optional :: timezone_offset
        integer,            intent(in), optional :: curstep

        real(real64),  allocatable, dimension(:) :: temp_times_64
        real(real128), allocatable, dimension(:) :: temp_times_128
        integer :: time_idx, error
        integer :: start_year, start_month, start_day, start_hour
        character(len=MAXSTRINGLENGTH) :: calendar, units
        real(real128) :: calendar_gain

        ! first read the time variable (presumebly a 1D real(real64) array)
        if (present(curstep)) then
            call io_read(trim(filename), trim(varname), temp_times_64, curstep=curstep)
        else
            call io_read(trim(filename), trim(varname), temp_times_64)
        endif

        ! attempt to read the calendar attribute from the time variable
        call io_read_attribute(trim(filename),"calendar", calendar, var_name=trim(varname), error=error)
        ! if time attribute it not present, set calendar to one specified in the config file
        if (error/=0) then
            if (this_image()==1) write(*,*) "WARNING: assuming standard/gregorian calendar for file "//trim(filename)
            calendar = "standard"
        endif

        ! attempt to read the units for this time variable
        call io_read_attribute(trim(filename), "units", units, var_name=trim(varname), error=error)

        ! if units attribute was present, then read information from it.
        if (error==0) then
            start_year    = year_from_units(units)
            start_month   = month_from_units(units)
            start_day     = day_from_units(units)
            start_hour    = hour_from_units(units)
            ! based off of the string "Days since" (or "seconds" or...)
            calendar_gain = time_gain_from_units(units)
        else

            stop "Time variable does not have units attribute"
        endif

        ! converts the input units to "days since ..."
        ! in case it is in units of e.g. "hours since" or "seconds since"
        allocate(temp_times_128(size(temp_times_64)))
        temp_times_128 = temp_times_64 / calendar_gain

        if (present(timezone_offset)) then
            temp_times_128 = temp_times_128 + timezone_offset / 24.0
        endif

        if (allocated(times)) deallocate(times)
        allocate(times(size(temp_times_128)))

        do time_idx = 1, size(temp_times_128,1)
            call times(time_idx)%init(calendar, start_year, start_month, start_day, start_hour)
            call times(time_idx)%set(days=temp_times_128(time_idx))
        end do

        deallocate(temp_times_64, temp_times_128)

    end subroutine read_times

    function get_output_time(time, units, round_seconds) result(output_time)
        implicit none
        type(Time_type),  intent(in) :: time
        character(len=*), intent(in), optional :: units
        logical,          intent(in), optional :: round_seconds

        type(Time_type) :: output_time
        type(time_delta_t) :: half_minute

        integer :: year, month, day, hour, minute, seconds
        integer :: year0, month0, day0, hour0, minute0, seconds0
        character(len=kMAX_NAME_LENGTH) :: use_units

        if (present(units)) then
            use_units = units
        else
            use_units = time%units()
        endif

        call time%date(year, month, day, hour, minute, seconds)
        year0 = year_from_units(use_units)
        month0 = month_from_units(use_units)
        day0 = day_from_units(use_units)
        hour0 = hour_from_units(use_units)
        minute0 = 0 ! minute_from_units(use_units)
        seconds0 = 0 ! seconds_from_units(use_units)

        call output_time%init(time%get_calendar(), year0, month0, day0, hour0)

        if (present(round_seconds)) then
            if (round_seconds) then
                if (seconds > 30) then
                    call output_time%set(year, month, day, hour, minute, seconds)
                    call half_minute%set(seconds=30)
                    output_time = output_time + half_minute
                    ! get a new date after adding 30 seconds
                    call output_time%date(year, month, day, hour, minute, seconds)
                    ! use that date after setting seconds to 0 this rounds the old date up by up to 30s
                    call output_time%set(year, month, day, hour, minute, 0)
                else
                    call output_time%set(year, month, day, hour, minute, 0)
                endif
            else
                call output_time%set(year, month, day, hour, minute, seconds)
            endif
        else
            call output_time%set(year, month, day, hour, minute, seconds)
        endif

    end function get_output_time

end module time_io
