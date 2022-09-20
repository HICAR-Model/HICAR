submodule(exchangeable_interface) exchangeable_implementation

  implicit none

  ! these are all global for a given image, so we save them in the module instead of in the individual objects
  integer, save :: halo_size
  integer, save :: north_neighbor, south_neighbor
  integer, save :: east_neighbor, west_neighbor
  integer, save, allocatable :: neighbors(:)

contains

  module subroutine const(this, grid, metadata, forcing_var)
    class(exchangeable_t),           intent(inout) :: this
    type(grid_t),                    intent(in)    :: grid
    type(variable_t),                intent(in),    optional :: metadata
    character(len=kMAX_NAME_LENGTH), intent(in),    optional :: forcing_var

    integer :: err

    halo_size = grid%halo_size


    this%north_boundary = (grid%yimg == grid%yimages)
    this%south_boundary = (grid%yimg == 1)
    this%east_boundary  = (grid%ximg == grid%ximages)
    this%west_boundary  = (grid%ximg == 1)

    if (associated(this%data_3d)) then
        deallocate(this%data_3d)
        nullify(this%data_3d)
    endif

    allocate(this%data_3d(grid%ims:grid%ime, &
                          grid%kms:grid%kme, &
                          grid%jms:grid%jme), stat=err)
    if (err /= 0) stop "exchangeable:dqdt_3d: Allocation request failed"
    this%data_3d = 0

    !Based on how grids are intialized, we can find out if this is a x/y staggered grid
    this%xe = 0
    this%ye = 0
    
    this%its = grid%its
    this%ite = grid%ite
    this%jts = grid%jts
    this%jte = grid%jte
    this%kts = grid%kts
    this%kte = grid%kte

    
    if ((grid%ns_halo_nx - (grid%nx_global / grid%ximages + 1)) > 0) then
        this%xe = 1
    else if ((grid%ew_halo_ny - (grid%ny_global / grid%yimages + 1)) > 0) then
        this%ye = 1
    endif

    allocate( this%halo_south_in( grid%ns_halo_nx+halo_size*2, grid%halo_nz, halo_size+this%ye   )[*])
    allocate( this%halo_north_in( grid%ns_halo_nx+halo_size*2, grid%halo_nz, halo_size           )[*])
    allocate( this%halo_east_in(  halo_size        ,  grid%halo_nz, grid%ew_halo_ny+halo_size*2  )[*])
    allocate( this%halo_west_in(  halo_size+this%xe,  grid%halo_nz, grid%ew_halo_ny+halo_size*2  )[*])


    if (.not.allocated(neighbors)) call this%set_neighbors(grid)

    if (present(metadata)) then
        this%meta_data = metadata
    else
        call this%meta_data%initialize(grid, forcing_var)
    endif
    
    this%meta_data%data_3d => this%data_3d
    this%meta_data%three_d = .True.
        
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
    type(variable_t),      intent(in),    optional :: metadata




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
  
  module subroutine exchange_y(this,do_metadata)
    class(exchangeable_t), intent(inout) :: this
    logical, optional, intent(in) :: do_metadata
    logical :: metadata
      
    metadata=.False.
    if (present(do_metadata)) metadata=do_metadata

    if (.not. this%north_boundary) call this%put_north(metadata)
    if (.not. this%south_boundary) call this%put_south(metadata)

    sync images( neighbors )
    
    if (.not. this%north_boundary) call this%retrieve_north_halo(metadata)
    if (.not. this%south_boundary) call this%retrieve_south_halo(metadata)
    
    if (.not. this%east_boundary)  call this%put_east(metadata)
    if (.not. this%west_boundary)  call this%put_west(metadata)

    sync images( neighbors )

    if (.not. this%east_boundary)  call this%retrieve_east_halo(metadata)
    if (.not. this%west_boundary)  call this%retrieve_west_halo(metadata)
    
  end subroutine
  
  module subroutine exchange_x(this,do_metadata)
    class(exchangeable_t), intent(inout) :: this
    logical, optional, intent(in) :: do_metadata
    logical :: metadata
      
    metadata=.False.
    if (present(do_metadata)) metadata=do_metadata

    if (.not. this%east_boundary)  call this%put_east(metadata)
    if (.not. this%west_boundary)  call this%put_west(metadata)

    sync images( neighbors )

    if (.not. this%east_boundary)  call this%retrieve_east_halo(metadata)
    if (.not. this%west_boundary)  call this%retrieve_west_halo(metadata)
    
    if (.not. this%north_boundary) call this%put_north(metadata)
    if (.not. this%south_boundary) call this%put_south(metadata)

    sync images( neighbors )
    
    if (.not. this%north_boundary) call this%retrieve_north_halo(metadata)
    if (.not. this%south_boundary) call this%retrieve_south_halo(metadata)
    
  end subroutine


  module subroutine exchange(this,do_metadata)
    class(exchangeable_t), intent(inout) :: this
    logical, optional, intent(in) :: do_metadata
    logical :: metadata
      
    metadata=.False.
    if (present(do_metadata)) metadata=do_metadata
    
    if (.not. this%north_boundary) call this%put_north(metadata)
    if (.not. this%south_boundary) call this%put_south(metadata)
    if (.not. this%east_boundary)  call this%put_east(metadata)
    if (.not. this%west_boundary)  call this%put_west(metadata)

    sync images( neighbors )

    if (.not. this%north_boundary) call this%retrieve_north_halo(metadata)
    if (.not. this%south_boundary) call this%retrieve_south_halo(metadata)
    if (.not. this%east_boundary)  call this%retrieve_east_halo(metadata)
    if (.not. this%west_boundary)  call this%retrieve_west_halo(metadata)
  end subroutine

  module subroutine put_north(this,do_metadata)
      class(exchangeable_t), intent(inout) :: this
      logical, optional, intent(in) :: do_metadata
      integer :: n, nx
      logical :: metadata
      
      metadata=.False.
      if (present(do_metadata)) metadata=do_metadata
      
      n = ubound(this%data_3d,3)
      nx = size(this%data_3d,1)
      
      if (metadata) then
          this%halo_south_in(1+halo_size:nx-halo_size,this%kts:this%kte,1:(halo_size+this%ye))[north_neighbor] = this%meta_data%dqdt_3d(this%its:this%ite,this%kts:this%kte,(n-halo_size*2+1-this%ye):(n-halo_size))
      else
          this%halo_south_in(1+halo_size:nx-halo_size,this%kts:this%kte,1:(halo_size+this%ye))[north_neighbor] = this%data_3d(this%its:this%ite,this%kts:this%kte,(n-halo_size*2+1-this%ye):(n-halo_size))
      endif
  end subroutine

  module subroutine put_south(this,do_metadata)
      class(exchangeable_t), intent(inout) :: this
      logical, optional, intent(in) :: do_metadata
      integer :: start, nx
      logical :: metadata
      
      metadata=.False.
      if (present(do_metadata)) metadata=do_metadata
      
      start = lbound(this%data_3d,3)
      nx = size(this%data_3d,1)
      
      if (metadata) then
          this%halo_north_in(1+halo_size:nx-halo_size,this%kts:this%kte,1:halo_size)[south_neighbor] = this%meta_data%dqdt_3d(this%its:this%ite,this%kts:this%kte,(start+halo_size+this%ye):(start+halo_size*2-1+this%ye))
      else
          this%halo_north_in(1+halo_size:nx-halo_size,this%kts:this%kte,1:halo_size)[south_neighbor] = this%data_3d(this%its:this%ite,this%kts:this%kte,(start+halo_size+this%ye):(start+halo_size*2-1+this%ye))
      endif
  end subroutine

  module subroutine retrieve_north_halo(this,do_metadata)
      class(exchangeable_t), intent(inout) :: this
      logical, optional, intent(in) :: do_metadata
      integer :: n, nx
      logical :: metadata
      
      metadata=.False.
      if (present(do_metadata)) metadata=do_metadata
      
      n = ubound(this%data_3d,3)
      nx = size(this%data_3d,1)
      
      if (metadata) then
          this%meta_data%dqdt_3d(this%its:this%ite,this%kts:this%kte,n-halo_size+1:n) = this%halo_north_in(1+halo_size:nx-halo_size,this%kts:this%kte,1:halo_size)
      else
          this%data_3d(this%its:this%ite,this%kts:this%kte,n-halo_size+1:n) = this%halo_north_in(1+halo_size:nx-halo_size,this%kts:this%kte,1:halo_size)
      endif
  end subroutine

  module subroutine retrieve_south_halo(this,do_metadata)
      class(exchangeable_t), intent(inout) :: this
      logical, optional, intent(in) :: do_metadata
      integer :: start, nx
      logical :: metadata
      
      metadata=.False.
      if (present(do_metadata)) metadata=do_metadata
      
      start = lbound(this%data_3d,3)
      nx = size(this%data_3d,1)
      
      if (metadata) then
          this%meta_data%dqdt_3d(this%its:this%ite,this%kts:this%kte,start:start+halo_size-1+this%ye) = this%halo_south_in(1+halo_size:nx-halo_size,this%kts:this%kte,1:(halo_size+this%ye))
      else
          this%data_3d(this%its:this%ite,this%kts:this%kte,start:start+halo_size-1+this%ye) = this%halo_south_in(1+halo_size:nx-halo_size,this%kts:this%kte,1:(halo_size+this%ye))
      endif
  end subroutine



  module subroutine put_east(this,do_metadata)
      class(exchangeable_t), intent(inout) :: this
      logical, optional, intent(in) :: do_metadata
      integer :: n, ny
      logical :: metadata
      
      metadata=.False.
      if (present(do_metadata)) metadata=do_metadata
      
      n = ubound(this%data_3d,1)
      ny = size(this%data_3d,3)
      
      if (metadata) then
          this%halo_west_in(1:(halo_size+this%xe),this%kts:this%kte,1+halo_size:ny-halo_size)[east_neighbor] = this%meta_data%dqdt_3d((n-halo_size*2+1-this%xe):(n-halo_size),this%kts:this%kte,this%jts:this%jte)
      else
          this%halo_west_in(1:(halo_size+this%xe),this%kts:this%kte,1+halo_size:ny-halo_size)[east_neighbor] = this%data_3d((n-halo_size*2+1-this%xe):(n-halo_size),this%kts:this%kte,this%jts:this%jte)
      endif
  end subroutine

  module subroutine put_west(this,do_metadata)
      class(exchangeable_t), intent(inout) :: this
      logical, optional, intent(in) :: do_metadata
      integer :: start, ny
      logical :: metadata
      
      metadata=.False.
      if (present(do_metadata)) metadata=do_metadata
      
      start = lbound(this%data_3d,1)
      ny = size(this%data_3d,3)
      
      if (metadata) then
          this%halo_east_in(1:halo_size,this%kts:this%kte,1+halo_size:ny-halo_size)[west_neighbor] = this%meta_data%dqdt_3d((start+halo_size+this%xe):(start+halo_size*2-1+this%xe),this%kts:this%kte,this%jts:this%jte)
      else
          this%halo_east_in(1:halo_size,this%kts:this%kte,1+halo_size:ny-halo_size)[west_neighbor] = this%data_3d((start+halo_size+this%xe):(start+halo_size*2-1+this%xe),this%kts:this%kte,this%jts:this%jte)
      endif
  end subroutine

  module subroutine retrieve_east_halo(this,do_metadata)
      class(exchangeable_t), intent(inout) :: this
      logical, optional, intent(in) :: do_metadata
      integer :: n, ny
      logical :: metadata
      
      metadata=.False.
      if (present(do_metadata)) metadata=do_metadata
      
      n = ubound(this%data_3d,1)
      ny = size(this%data_3d,3)
      
      if (metadata) then
          this%meta_data%dqdt_3d(n-halo_size+1:n,this%kts:this%kte,this%jts:this%jte) = this%halo_east_in(1:halo_size,this%kts:this%kte,1+halo_size:ny-halo_size)
      else
          this%data_3d(n-halo_size+1:n,this%kts:this%kte,this%jts:this%jte) = this%halo_east_in(1:halo_size,this%kts:this%kte,1+halo_size:ny-halo_size)
      endif
  end subroutine

  module subroutine retrieve_west_halo(this,do_metadata)
      class(exchangeable_t), intent(inout) :: this
      logical, optional, intent(in) :: do_metadata
      integer :: start, ny
      logical :: metadata
      
      metadata=.False.
      if (present(do_metadata)) metadata=do_metadata
      
      start = lbound(this%data_3d,1)
      ny = size(this%data_3d,3)
      
      if (metadata) then
          this%meta_data%dqdt_3d(start:start+halo_size-1+this%xe,this%kts:this%kte,this%jts:this%jte) = this%halo_west_in(1:(halo_size+this%xe),this%kts:this%kte,1+halo_size:ny-halo_size)
      else
          this%data_3d(start:start+halo_size-1+this%xe,this%kts:this%kte,this%jts:this%jte) = this%halo_west_in(1:(halo_size+this%xe),this%kts:this%kte,1+halo_size:ny-halo_size)
      endif
  end subroutine
end submodule
