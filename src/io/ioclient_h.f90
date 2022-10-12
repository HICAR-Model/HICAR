!>----------------------------------------------------------
!!  Define the interface for the ioclient
!!
!!  The I/O client serves as a parallelized interface between the main program
!!  and the I/O routines. The I/O client handles the buffering and exchange
!!  of data between the processes of the main program (child processes)
!!  and the I/O processes (parent processes)
!! 
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!----------------------------------------------------------
module ioclient_interface
  use mpi
  use netcdf
  use icar_constants
  use variable_interface, only : variable_t
  use boundary_interface, only : boundary_t
  use domain_interface,   only : domain_t

!  use time_object,        only : Time_type, THREESIXTY, GREGORIAN, NOCALENDAR, NOLEAP

  implicit none

  private
  public :: ioclient_t

  !>----------------------------------------------------------
  !! Output type definition
  !!
  !!----------------------------------------------------------
  type :: ioclient_t

      ! all components are private and should only be modified through procedures
      private

      ! Store the variables to be written
      ! Note n_variables may be smaller then size(variables) so that it doesn't
      ! have to keep reallocating variables whenever something is added or removed
      integer, public :: n_input_variables, n_output_variables, IO_comms

      type(variable_t), public, allocatable :: variables(:)
      ! time variable , publicis stored outside of the variable list... probably need to think about that some
      
      ! store status of the object -- are we a parent or child process
      logical :: creating = .false.
            
      integer, public :: server
      
      integer, public ::  i_s_w, i_e_w, k_s_w, k_e_w, j_s_w, j_e_w, n_w, i_s_r, i_e_r, k_s_r, k_e_r, j_s_r, j_e_r, n_r
      integer :: restart_counter = 0
      integer :: output_counter = 0
      integer :: frames_per_outfile, restart_count

  contains

      procedure, public  :: push
      procedure, public  :: receive

      procedure, public  :: init
  end type

  interface

      !>----------------------------------------------------------
      !! Initialize the object (e.g. allocate the variables array)
      !!
      !!----------------------------------------------------------
      module subroutine init(this, domain, forcing)
          implicit none
          class(ioclient_t),  intent(inout)  :: this
          type(domain_t),     intent(in)     :: domain
          type(boundary_t),   intent(in)     :: forcing
      end subroutine

      !>----------------------------------------------------------
      !! Increase the size of the variables array if necessary
      !!
      !!----------------------------------------------------------
      module subroutine push(this, domain, write_buffer)
          implicit none
          class(ioclient_t),   intent(inout) :: this
          type(domain_t),   intent(inout)    :: domain
          real, intent(inout), allocatable   :: write_buffer(:,:,:,:)[:]
      end subroutine

      !>----------------------------------------------------------
      !! Set the domain data structure to be used when writing
      !!
      !!----------------------------------------------------------
      module subroutine receive(this, forcing, read_buffer)
          implicit none
          class(ioclient_t), intent(inout) :: this
          type(boundary_t), intent(inout)  :: forcing
          real, intent(in), allocatable    :: read_buffer(:,:,:,:)[:]
      end subroutine

      module subroutine close_files(this)
          class(ioclient_t), intent(inout) :: this
      end subroutine
      
  end interface
end module
