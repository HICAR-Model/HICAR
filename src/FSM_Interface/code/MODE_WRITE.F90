module MODE_WRITE

  contains
  subroutine WRITE_2D(var,fileid)
    ! 
    ! BC created August 2021
    !  Write 2D variable var into fileid (append mode)
    ! @TODO: interface to handle non-real variables.
    implicit none
    real, dimension(:,:), intent(in):: var
    integer, intent(in) :: fileid
    integer :: pos
    inquire(unit = fileid, pos = pos)
    write(fileid, pos = pos) var
  end subroutine WRITE_2D

end module MODE_WRITE
