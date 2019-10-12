subroutine print_rec_matrix (iunit, n, a, label)

  ! print a n by n matrix to file (fortran unit = 'unit'),
  ! with a label stored in string 'label'

  use definitions

  implicit none

  integer, intent(in) :: iunit, n
  real(KREAL), intent(in) :: a(n,n)
  character*(*) :: label

  integer :: dbg, i
  character*100 :: fmt

  ! ---------------------------------------------------------------------------

  dbg = 0

  ! output format of matrix
  write(fmt,'(a,i5,a)') '(1x,', n, 'f20.8)'
  if (dbg>0) write (iunit,*) 'fmt =', fmt
  
  ! print matrix
  write(iunit,'(/1x,a/)') trim(label)
  do i = 1,n
    write (iunit, fmt) a(i,:)
  end do
  write (iunit,*)

end subroutine print_rec_matrix
