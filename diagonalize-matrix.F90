subroutine diagonalize_matrix (n, a)

  ! diagonalize a complex hermitean matrix with Lapack routines

  ! (c) Jochen Autschbach, SUNY Buffalo

  use definitions

  implicit none

  integer(KINT), intent(in)  :: n
  complex(KREAL), intent(inout) :: a(n,n)

  integer(KINT), dimension(:), allocatable :: iwork
  real((KREAL)), dimension(:), allocatable :: rwork, eval
  complex(KREAL), dimension(:), allocatable :: work
  real(KREAL), parameter :: small=1.0d-10, one=1.0d0, zero=0.0d0
  integer :: i, lda, lwork, lrwork, liwork, info, dbg
  character*1 :: jobz, uplo

  integer, parameter :: out = 6, err = 0

  ! ----------------------------------------------------------------------------
  
  dbg = 1 ! debug switch

  if (dbg>0) then
    call print_rec_matrix(out, n, real(a), &
      'Zeeman Hamiltonian REAL part')
    call print_rec_matrix(out, n, aimag(a), &
      'Zeeman Hamiltonian IMAG part')
  end if



  ! lapack routine dimensioning:

  liwork = 3 + 5 * n
  lrwork = 1 + 5 * n + 2 * n**2
  lwork = 2*n + n**2
  lda = n
  jobz = 'V'
  uplo = 'L'
  info = 0

  ! allocate arrays: 

  allocate(iwork(liwork))
  allocate(work(lwork))
  allocate(rwork(lrwork))

  allocate (eval(n))
      
  ! call lapack routine to diagonalize the matrix a:

  call zheevd( jobz, uplo, n, a, lda, eval, work, lwork, rwork, &
     &                   lrwork, iwork, liwork, info )

  if (info .ne. 0) then
    write (out,*) 'info =',info
    stop 'zheevd error. aborting'
  end if

  ! free up some workspace 

  deallocate(iwork,work,rwork)

  if (dbg>1) write(out,*) a(:,:)

  ! in order to compare with Mathematica et al., 
  ! let's fix the phases of the eigenvectors such that
  ! the first element is positive

  do i = 1, n
   if (abs(aimag(a(1,i))).gt.small) &
    stop 'eigenvector array has non-real first row'
   if (real(a(1,i)) < 0) a(:,i) = -a(:,i) 
  end do

  if (dbg>0) then
    call print_rec_matrix(out, n, real(a), &
      'Eigenvectors of Zeeman Hamiltonian REAL part')
    call print_rec_matrix(out, n, aimag(a), &
      'Eigenvectors of Zeeman Hamiltonian IMAG part')
  end if

  write (out,'(1x,a)') 'Eigenvalues of Zeeman Operator'
  write(out,'((1x,f20.8))') eval(:)

  ! clean up and exit

  deallocate(eval)

  return

end subroutine diagonalize_matrix
