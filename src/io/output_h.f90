!>----------------------------------------------------------
!!  Define the interface for the output object
!!
!!  Output objects store all of the data and references to data necessary to write
!!  an output file.  This includes primarily internal netcdf related IDs.
!!  Output objects also store an array of variables to output.
!!  These variables maintain pointers to the data to be output as well as
!!  Metadata (e.g. dimension names, units, other attributes)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module output_interface
  use mpi
  use netcdf
  use icar_constants
  use variable_interface, only : variable_t
  use domain_interface,   only : domain_t
  use meta_data_interface,only : meta_data_t
  use time_object,        only : Time_type, THREESIXTY, GREGORIAN, NOCALENDAR, NOLEAP
  use options_interface,  only : options_t

  implicit none

  private
  public :: output_t

  !>----------------------------------------------------------
  !! Output type definition
  !!
  !!----------------------------------------------------------
  type, extends(meta_data_t) :: output_t

      ! all components are private and should only be modified through procedures
      private

      ! Store the variables to be written
      ! Note n_variables may be smaller then size(variables) so that it doesn't
      ! have to keep reallocating variables whenever something is added or removed
      integer, public :: n_vars = 0
      type(variable_t), public, allocatable :: variables(:)
      ! time variable , publicis stored outside of the variable list... probably need to think about that some
      type(variable_t) :: time

      ! store status of the object
      logical :: is_initialized = .false.
      logical :: creating = .false.
      logical :: block_checked = .false.
      logical :: is_blocked = .false.
      logical :: blocked_LL = .false.
      logical :: blocked_UR = .false.

      ! The filename of the netcdf file to write
      character(len=kMAX_FILE_LENGTH), public :: filename, base_file_name

      character(len=49)   :: file_date_format = '(I4,"-",I0.2,"-",I0.2,"_",I0.2,"-",I0.2,"-",I0.2)'

      ! the netcdf ID for an open file
      integer :: ncfile_id = -1

      ! number of dimensions in the file
      integer :: n_dims = 0

      integer :: start_3d(3), cnt_3d(3), cnt_2d(2)
      
      !same as above, but for block, if present
      integer :: start_3d_b(3), cnt_3d_b(3), cnt_2d_b(2)
      integer :: start_3d_b2(3), cnt_3d_b2(3), cnt_2d_b2(2)

      integer :: its, ite, kts, kte, jts, jte

      integer :: output_counter, output_count

      integer :: global_dim_len(3)

      ! list of netcdf dimension IDs
      integer :: dim_ids(kMAX_DIMENSIONS)
      ! name of the dimensions in the file
      character(len=kMAX_DIM_LENGTH) :: dimensions(kMAX_DIMENSIONS)

  contains

      procedure, public  :: add_to_output
      procedure, public  :: add_variables
      procedure, public  :: set_attrs
      procedure, public  :: save_file
      procedure, public  :: close_file

      procedure, public  :: init
      procedure, private :: increase_var_capacity
  end type

  interface

      !>----------------------------------------------------------
      !! Initialize the object (e.g. allocate the variables array)
      !!
      !!----------------------------------------------------------
      module subroutine init(this, domain, options, its, ite, kts, kte, jts, jte)
        implicit none
        class(output_t),  intent(inout)  :: this
        type(domain_t),   intent(in)     :: domain
        type(options_t),  intent(in)     :: options
        integer,          intent(in)     :: its, ite, kts, kte, jts, jte

      end subroutine

      !>----------------------------------------------------------
      !! Increase the size of the variables array if necessary
      !!
      !!----------------------------------------------------------
      module subroutine increase_var_capacity(this)
          implicit none
          class(output_t),   intent(inout)  :: this
      end subroutine

      !>----------------------------------------------------------
      !! Set the domain data structure to be used when writing
      !!
      !!----------------------------------------------------------
      module subroutine set_attrs(this, domain)
          implicit none
          class(output_t),  intent(inout)  :: this
          type(domain_t),   intent(in)     :: domain
      end subroutine

      !>----------------------------------------------------------
      !! Add a variable to the list of output variables
      !!
      !!----------------------------------------------------------
      module subroutine add_to_output(this, in_variable)
          implicit none
          class(output_t),   intent(inout)  :: this
          type(variable_t),  intent(in)     :: in_variable
      end subroutine

      !>----------------------------------------------------------
      !! Add multiple variables to the list of output variables given their integer constant index
      !!
      !!----------------------------------------------------------
      module subroutine add_variables(this, var_list)
          implicit none
          class(output_t), intent(inout)  :: this
          integer,         intent(in)     :: var_list(:)
      end subroutine

      !>----------------------------------------------------------
      !! Save a new timestep (time) to the file "filename"
      !!
      !!----------------------------------------------------------
      module subroutine save_file(this, time, par_comms)
          implicit none
          class(output_t),  intent(inout) :: this
          type(Time_type),  intent(in)    :: time
          integer,          intent(in)    :: par_comms
      end subroutine

    module subroutine close_file(this)
        implicit none
        class(output_t),   intent(inout)  :: this
    end subroutine

  end interface
end module
