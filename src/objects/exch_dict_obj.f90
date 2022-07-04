!>------------------------------------------------
!! Defines the implementation for a variable dictionary
!!
!! Dictionary can be accessed by string keys
!! Dictionary stores variable_t types
!!
!!------------------------------------------------
submodule(variable_dict_interface) variable_dict_implementation

    implicit none

    ! note that this dictionary is stored fairly inefficiently (both speed and space!)
    ! if this value gets larger one may need to think about re-writing to use a hash Table
    ! or similar
    integer, parameter :: kMIN_VAR_DICT_SIZE = 128

contains

    !>------------------------------------------------
    !! Initialize the variable dictionary
    !!
    !! Allocates the var_list array to a minimum size
    !! sets the max_vars and n_vars to the correct initial values
    !! then sets initialized==true
    !!
    !! Typically assumes lazy instantiation, e.g. not initialized until used
    !!
    !!------------------------------------------------
    module subroutine init(this)
        implicit none
        class(var_dict_t),   intent(inout)  :: this

        if (this%initialized) return

        if (allocated(this%var_list)) then
            stop "var_dict list is somehow allocated before var_dict was initialized"
        endif

        allocate(this%var_list(kMIN_VAR_DICT_SIZE))

        this%max_vars = kMIN_VAR_DICT_SIZE
        this%n_vars = 0

        this%initialized = .True.

    end subroutine


    !>------------------------------------------------
    !! Retrieve a variable for a given key
    !!
    !! @param varname : The key to look for in the dictionary
    !! @param err     : optional, if available this will be set to indicate if the key is not found
    !!
    !!------------------------------------------------
    module function get_var(this, varname, err) result(var_data)
        implicit none
        class(var_dict_t),   intent(in) :: this
        character(len=*),    intent(in) :: varname
        integer,             intent(out),  optional :: err

        type(variable_t)                :: var_data

        integer :: i

        if (present(err)) err = 0

        ! if this dictionary hasn't been initialized yet that is bad
        if (.not.this%initialized) then

            ! if the user requested error codes, simply return an error
            if (present(err)) then
                err = 2
                return
            endif

            ! otherwise stop execution
            stop "variable_dictionary read before being initialized"
        endif

        ! This is the core logic
        ! loop through the stored variables looking for the requested key
        do i = 1, this%n_vars

            ! if this key matches the supplied key, return the variable
            if (trim(this%var_list(i)%name) == trim(varname)) then
                var_data = this%var_list(i)%var
                return
            endif
        end do

        ! if we get here, we must not have found the key!
        if (present(err)) then
            err = 1
        else
            ! if the user did not request an error code, then we have to stop
            write(*,*) "Searching for : "//trim(varname)
            stop "ERROR: var not found in dictionary"
        endif

    end function

    !>------------------------------------------------
    !! Add a supplied variable to the dictionary using the key supplied
    !!
    !! @param varname       key to store for the var_data supplied
    !! @param var_data      variable to store associated with key
    !! @param domain_var    optional additional variable to associated with key
    !! @param save_state    optional store the input array data in a locally allocated array instead of just pointing to it
    !! @param err           optional value to store an error code in so execution doesn't have to stop
    !!
    !!------------------------------------------------
    module subroutine add_var(this, varname, var_data, err)
        implicit none
        class(var_dict_t),   intent(inout)  :: this
        character(len=*),    intent(in)     :: varname
        type(variable_t),    intent(in)     :: var_data
        integer,             intent(out),optional :: err

        integer :: ims,ime, jms,jme, kms,kme

        if (.not.this%initialized) then
            call this%init()
        endif

        if (this%n_vars==this%max_vars) then
            if (present(err)) then
                err = 3
                return
            endif
            write(*,*) this_image(), this%n_vars, this%max_vars
            stop "Ran out of space in var_dict"
        endif

        this%n_vars = this%n_vars + 1
        this%var_list(this%n_vars)%name = varname

        ! warning, this copies all information directly from the var_data into the dict, including the POINTER to the data
        this%var_list(this%n_vars)%var  = var_data
        this%var_list(this%n_vars)%var%name = varname


    end subroutine



    !>------------------------------------------------
    !! Reset the iteration sequence in the dictionary
    !!
    !!------------------------------------------------
    module subroutine reset_iterator(this)
        implicit none
        class(var_dict_t),   intent(inout)  :: this

        this%current_variable = 1
    end subroutine

    !>------------------------------------------------
    !! Test to see if there are more values in the iteration sequence
    !!
    !! @result boolean      returns true of there are values left to iterate over
    !!
    !!------------------------------------------------
    module function has_more_elements(this) result(boolean)
        implicit none
        class(var_dict_t),   intent(in) :: this
        logical :: boolean

        boolean = this%current_variable <= this%n_vars
    end function

    !>------------------------------------------------
    !! Get the next variable in the iteration
    !!
    !! @param name     optional, returns the name associated with the variable
    !! @param err      optional, returns an error if there are no variables left to iterate
    !!
    !!------------------------------------------------
    module function next(this, name, err) result(var_data)
        implicit none
        class(var_dict_t),   intent(inout)  :: this
        character(len=*),    intent(out),   optional :: name
        integer,             intent(out),   optional :: err
        type(variable_t)                    :: var_data

        if (present(err)) err = 0

        if (this%has_more_elements()) then
            ! get the variable data to return
            var_data = this%var_list(this%current_variable)%var

            ! if requested, get the variable name too
            if (present(name)) name = this%var_list(this%current_variable)%name

            ! step the iterator forward
            this%current_variable = this%current_variable + 1
        else
            if (present(err)) then
                err = 1
            else
                stop "Attempt to iterate past the end of the dictionary"
            endif
        endif

    end function


end submodule
