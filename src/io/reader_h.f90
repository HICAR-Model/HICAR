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
module reader_interface
  use mpi
  use netcdf
  use icar_constants
  use options_interface,  only : options_t
  use variable_dict_interface,  only : var_dict_t
  use boundary_interface, only : boundary_t
  use meta_data_interface,only : meta_data_t
  use time_object,        only : Time_type

  implicit none

  private
  public :: reader_t

  !>----------------------------------------------------------
  !! Output type definition
  !!
  !!----------------------------------------------------------
  type reader_t

      ! all components are private and should only be modified through procedures
      private

      ! Store the variables to be written
      ! Note n_variables may be smaller then size(variables) so that it doesn't
      ! have to keep reallocating variables whenever something is added or removed
      integer, public :: n_vars = 0
      
      type(var_dict_t)    :: variables      ! a dictionary with all forcing data
      
      ! list of input files
      character (len=kMAX_FILE_LENGTH), allocatable :: file_list(:)
      character (len=kMAX_NAME_LENGTH)   :: time_var

      ! the netcdf ID for an open file
      integer :: ncfile_id
      
      integer :: its, ite, kts, kte, jts, jte

      integer :: curfile, curstep
  contains
      procedure, public :: init
      procedure, public :: read_next_step
      procedure, public :: close_file
  end type

  interface

    module subroutine init(this, its, ite, kts, kte, jts, jte, options)
        class(reader_t), intent(inout) :: this
        integer, intent(in) :: its, ite, kts, kte, jts, jte
        type(options_t), intent(in) :: options
    end subroutine

      !>----------------------------------------------------------
      !! Read the next timestep (time) from the input file list
      !!
      !!----------------------------------------------------------
    module subroutine read_next_step(this, buffer, par_comms)
          implicit none
          class(reader_t), intent(inout)   :: this
          real, allocatable, intent(inout) :: buffer(:,:,:,:)
          integer, intent(in)              :: par_comms
    end subroutine

    module subroutine close_file(this)
        implicit none
        class(reader_t),   intent(inout)  :: this
    end subroutine

  end interface
end module
