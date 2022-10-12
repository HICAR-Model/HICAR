!>----------------------------------------------------------
!!  Define the interface for the ioserver
!!
!!  The I/O server serves as a parallelized interface between the main program
!!  and the I/O routines. The I/O server handles the buffering and exchange
!!  of data between the processes of the main program (child processes)
!!  and the I/O processes (parent processes)
!! 
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!----------------------------------------------------------
module ioserver_interface
  use mpi
  use netcdf
  use icar_constants
  use variable_interface, only : variable_t
  use reader_interface,   only : reader_t
  use output_interface,   only : output_t
  use options_interface,  only : options_t
  use time_object,        only : Time_type
  use boundary_interface, only : boundary_t
  use domain_interface,   only : domain_t

!  use time_object,        only : Time_type, THREESIXTY, GREGORIAN, NOCALENDAR, NOLEAP

  implicit none

  private
  public :: ioserver_t

  !>----------------------------------------------------------
  !! Output type definition
  !!
  !!----------------------------------------------------------
  type :: ioserver_t

      ! all components are private and should only be modified through procedures
      private

      ! Store the variables to be written
      ! Note n_variables may be smaller then size(variables) so that it doesn't
      ! have to keep reallocating variables whenever something is added or removed
      integer, public :: n_input_variables, n_output_variables, n_children, IO_comms
      integer, public :: server_id = -1
      type(variable_t), public, allocatable :: variables(:)
      ! time variable , publicis stored outside of the variable list... probably need to think about that some
      type(output_t) :: outputer, rster
      type(reader_t) :: reader
      
      type(domain_t) :: output_domain

      
      real, pointer :: parent_write_buffer(:,:,:,:)

      ! store status of the object -- are we a parent or child process
      logical :: creating = .false.
      
      ! coarray-indices of child io processes, indexed according to the COMPUTE_TEAM
      integer, public, allocatable :: children(:)
      
      ! coarray-index of parent io process, indexed according to the IO_TEAM
      integer :: parent_id = 0
            
      ! indices of process extent. For child processes, this is just the tile extents, 
      ! except for edges, where we go to the domain extents
      !Convention is:
      ! (i/k/j) -- x, z, or y index
      ! (s/e) -- start or end index
      ! (r/w) -- index for the read or write buffer
      ! c -- defines this index as the index of a particular child, 
      ! where the index in this array corresponds to the child image at the same index in the children array

      integer, allocatable, dimension(:) :: isrc, ierc, ksrc, kerc, jsrc, jerc, iswc, iewc, kswc, kewc, jswc, jewc
      
      integer, public ::  i_s_w, i_e_w, k_s_w, k_e_w, j_s_w, j_e_w, n_w, i_s_r, i_e_r, k_s_r, k_e_r, j_s_r, j_e_r, n_r
      integer :: restart_counter = 0
      integer :: output_counter = 0
      integer :: frames_per_outfile, restart_count


      ! The filename of the netcdf file to write
      character(len=kMAX_FILE_LENGTH) :: output_file_name, restart_file_name

      ! number of dimensions in the file
      integer :: n_dims = 0

      ! list of netcdf dimension IDs
      integer :: dim_ids(kMAX_DIMENSIONS)
      ! name of the dimensions in the file
      character(len=kMAX_DIM_LENGTH) :: dimensions(kMAX_DIMENSIONS)

  contains

      procedure, public  :: write_file
      procedure, public  :: read_file
      procedure, public  :: close_files

      procedure, public  :: init
  end type

  interface

      !>----------------------------------------------------------
      !! Initialize the object (e.g. allocate the variables array)
      !!
      !!----------------------------------------------------------
      module subroutine init(this, domain, options, isrc, ierc, ksrc, kerc, jsrc, jerc, iswc, iewc, kswc, kewc, jswc, jewc)
          implicit none
          class(ioserver_t),   intent(inout)  :: this
          type(domain_t),   intent(in)        :: domain
          type(options_t), intent(in) :: options
          integer, allocatable, dimension(:), intent(in) :: isrc, ierc, ksrc, kerc, jsrc, jerc, iswc, iewc, kswc, kewc, jswc, jewc

      end subroutine

      !>----------------------------------------------------------
      !! Increase the size of the variables array if necessary
      !!
      !!----------------------------------------------------------
      module subroutine write_file(this, domain, time, write_buffer)
          implicit none
          class(ioserver_t),   intent(inout)  :: this
          type(domain_t),   intent(in)        :: domain
          type(Time_type),  intent(in)        :: time
          real, allocatable, intent(in)       :: write_buffer(:,:,:,:)[:]
      end subroutine

      !>----------------------------------------------------------
      !! Set the domain data structure to be used when writing
      !!
      !!----------------------------------------------------------
      module subroutine read_file(this, forcing, read_buffer)
          implicit none
          class(ioserver_t),  intent(inout)  :: this
          type(boundary_t),   intent(inout)  :: forcing
          real, allocatable, intent(inout)   :: read_buffer(:,:,:,:)[:]
      end subroutine

      module subroutine close_files(this)
          class(ioserver_t), intent(inout) :: this
      end subroutine
      
  end interface
end module
