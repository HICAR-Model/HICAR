submodule(exchangeable_interface) exchangeable_implementation

  implicit none

  ! these are all global for a given image, so we save them in the module instead of in the individual objects
  integer, save :: halo_size
  integer, save :: north_neighbor, south_neighbor
  integer, save :: east_neighbor, west_neighbor
  integer, save, allocatable :: neighbors(:)

contains

  module subroutine const(this, grid, metadata)
    class(exchangeable_t), intent(inout) :: this
    type(grid_t),          intent(in)    :: grid
    class(variable_t),     intent(in), optional :: metadata

    halo_size = grid%halo_size

    this%north_boundary = (grid%yimg == grid%yimages)
    this%south_boundary = (grid%yimg == 1)
    this%east_boundary  = (grid%ximg == grid%ximages)
    this%west_boundary  = (grid%ximg == 1)

    if (associated(this%local)) then
        deallocate(this%local)
        nullify(this%local)
    endif

    allocate(this%local(grid%ims:grid%ime, &
                        grid%kms:grid%kme, &
                        grid%jms:grid%jme))

    allocate( this%halo_south_in( grid%ns_halo_nx+halo_size*2, grid%halo_nz,   halo_size    )[*])
    allocate( this%halo_north_in( grid%ns_halo_nx+halo_size*2, grid%halo_nz,   halo_size    )[*])
    allocate( this%halo_east_in(    halo_size,     grid%halo_nz, grid%ew_halo_ny+halo_size*2)[*])
    allocate( this%halo_west_in(    halo_size,     grid%halo_nz, grid%ew_halo_ny+halo_size*2)[*])


    if (.not.allocated(neighbors)) call this%set_neighbors(grid)

    if (present(metadata)) call this%set_outputdata(metadata)

  end subroutine

  module subroutine set_neighbors(this, grid)
      class(exchangeable_t), intent(inout) :: this
      type(grid_t),          intent(in)    :: grid
      integer :: n_neighbors, current

      ! set up the neighbors array so we can sync with our neighbors when needed
      if (.not.allocated(neighbors)) then
        associate(me=>this_image())
          south_neighbor = me - grid%ximages
          north_neighbor = me + grid%ximages
          east_neighbor  = me + 1
          west_neighbor  = me - 1

          n_neighbors = merge(0,1,this%south_boundary)  &
                       +merge(0,1,this%north_boundary)  &
                       +merge(0,1,this%east_boundary)   &
                       +merge(0,1,this%west_boundary)
          n_neighbors = max(1, n_neighbors)

          allocate(neighbors(n_neighbors))

          current = 1
          if (.not. this%south_boundary) then
              neighbors(current) = south_neighbor
              current = current+1
          endif
          if (.not. this%north_boundary) then
              neighbors(current) = north_neighbor
              current = current+1
          endif
          if (.not. this%east_boundary) then
              neighbors(current) = east_neighbor
              current = current+1
          endif
          if (.not. this%west_boundary) then
              neighbors(current) = west_neighbor
              current = current+1
          endif
          ! if current = 1 then all of the boundaries were set, just store ourself as our "neighbor"
          if (current == 1) then
              neighbors(current) = me
          endif

        end associate
      endif

  end subroutine

  module subroutine set_outputdata(this, metadata)
    implicit none
    class(exchangeable_t), intent(inout)  :: this
    class(variable_t),     intent(in)     :: metadata

    this%meta_data = metadata

    this%meta_data%data_3d => this%local
    this%meta_data%three_d = .True.

    if (.not.allocated(this%meta_data%dim_len)) allocate(this%meta_data%dim_len(3))
    this%meta_data%dim_len(1) = size(this%local,1)
    this%meta_data%dim_len(2) = size(this%local,2)
    this%meta_data%dim_len(3) = size(this%local,3)

  end subroutine


  module subroutine send(this)
    class(exchangeable_t), intent(inout) :: this
    if (.not. this%north_boundary) call this%put_north
    if (.not. this%south_boundary) call this%put_south
    if (.not. this%east_boundary)  call this%put_east
    if (.not. this%west_boundary)  call this%put_west
  end subroutine

  module subroutine retrieve(this, no_sync)
    class(exchangeable_t), intent(inout) :: this
    logical,               intent(in),   optional :: no_sync

    if (.not. present(no_sync)) then
        sync images( neighbors )
    else
        if (.not. no_sync) then
            sync images( neighbors )
        endif
    endif

    if (.not. this%north_boundary) call this%retrieve_north_halo
    if (.not. this%south_boundary) call this%retrieve_south_halo
    if (.not. this%east_boundary) call this%retrieve_east_halo
    if (.not. this%west_boundary) call this%retrieve_west_halo
  end subroutine

  module subroutine exchange(this)
    class(exchangeable_t), intent(inout) :: this
    if (.not. this%north_boundary) call this%put_north
    if (.not. this%south_boundary) call this%put_south
    if (.not. this%east_boundary)  call this%put_east
    if (.not. this%west_boundary)  call this%put_west

    sync images( neighbors )

    if (.not. this%north_boundary) call this%retrieve_north_halo
    if (.not. this%south_boundary) call this%retrieve_south_halo
    if (.not. this%east_boundary)  call this%retrieve_east_halo
    if (.not. this%west_boundary)  call this%retrieve_west_halo
  end subroutine

  module subroutine put_north(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: n, nx
      n = ubound(this%local,3)
      nx = size(this%local,1)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        call assert( shape(this%halo_south_in(:nx,:,1:halo_size)[north_neighbor]) &
                     == shape(this%local(:,:,n-halo_size+1:n)),         &
                     "put_north: conformable halo_south_in and local " )
      end if

      !dir$ pgas defer_sync
      this%halo_south_in(1:nx,:,1:halo_size)[north_neighbor] = this%local(:,:,n-halo_size*2+1:n-halo_size)
  end subroutine

  module subroutine put_south(this)
      class(exchangeable_t), intent(inout) :: this

      integer :: start, nx
      start = lbound(this%local,3)
      nx = size(this%local,1)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        call assert( shape(this%halo_north_in(:nx,:,1:halo_size)[south_neighbor]) &
                     == shape(this%local(:,:,start:start+halo_size-1)), &
                     "put_south: conformable halo_north_in and local " )
      end if
      !dir$ pgas defer_sync
      this%halo_north_in(1:nx,:,1:halo_size)[south_neighbor] = this%local(:,:,start+halo_size:start+halo_size*2-1)
  end subroutine

  module subroutine retrieve_north_halo(this)
      class(exchangeable_t), intent(inout) :: this

      integer :: n, nx
      n = ubound(this%local,3)
      nx = size(this%local,1)

      this%local(:,:,n-halo_size+1:n) = this%halo_north_in(:nx,:,1:halo_size)
  end subroutine

  module subroutine retrieve_south_halo(this)
      class(exchangeable_t), intent(inout) :: this

      integer :: start, nx
      start = lbound(this%local,3)
      nx = size(this%local,1)

      this%local(:,:,start:start+halo_size-1) = this%halo_south_in(:nx,:,1:halo_size)
  end subroutine

  module subroutine put_east(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: n, ny
      n = ubound(this%local,1)
      ny = size(this%local,3)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        call assert( shape(this%halo_west_in(1:halo_size,:,:ny)[east_neighbor])       &
                     == shape(this%local(n-halo_size*2+1:n-halo_size,:,:)), &
                     "put_east: conformable halo_west_in and local " )
      end if

      !dir$ pgas defer_sync
      this%halo_west_in(1:halo_size,:,1:ny)[east_neighbor] = this%local(n-halo_size*2+1:n-halo_size,:,:)
  end subroutine

  module subroutine put_west(this)
      class(exchangeable_t), intent(inout) :: this

      integer :: start, ny
      start = lbound(this%local,1)
      ny = size(this%local,3)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        call assert( shape(this%halo_east_in(1:halo_size,:,:ny)[west_neighbor])               &
                     == shape(this%local(start+halo_size:start+halo_size*2-1,:,:)), &
                     "put_west: conformable halo_east_in and local " )
      end if

      !dir$ pgas defer_sync
      this%halo_east_in(1:halo_size,:,1:ny)[west_neighbor] = this%local(start+halo_size:start+halo_size*2-1,:,:)
  end subroutine

  module subroutine retrieve_east_halo(this)
      class(exchangeable_t), intent(inout) :: this

      integer :: n, ny
      n = ubound(this%local,1)
      ny = size(this%local,3)

      this%local(n-halo_size+1:n,:,:) = this%halo_east_in(1:halo_size,:,1:ny)
  end subroutine

  module subroutine retrieve_west_halo(this)
      class(exchangeable_t), intent(inout) :: this

      integer :: start, ny
      start = lbound(this%local,1)
      ny = size(this%local,3)

      this%local(start:start+halo_size-1,:,:) = this%halo_west_in(1:halo_size,:,1:ny)
  end subroutine

end submodule