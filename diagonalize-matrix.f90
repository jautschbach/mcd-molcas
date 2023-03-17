subroutine diagonalize_matrix (n, a)

  ! diagonalize a complex hermitean matrix with Lapack routines

  ! (c) Jochen Autschbach, SUNY Buffalo

  use definitions

  use constants_parameters

  use shared_variables

  implicit none

  integer(KINT), intent(in)  :: n
  complex(KREAL), intent(inout) :: a(n,n)

  integer(KINT), dimension(:), allocatable :: iwork
  real((KREAL)), dimension(:), allocatable :: rwork, eval
  complex(KREAL), dimension(:), allocatable :: work
  integer :: i, lda, lwork, lrwork, liwork, info
  character*1 :: jobz, uplo


  ! ----------------------------------------------------------------------------
  

  if (dbg>1) then
    write (out,*) 'hello from diagonalize_matrix'
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
   if (abs(aimag(a(1,i))).gt.tiny) &
    stop 'eigenvector array has non-real first row'
   if (real(a(1,i)) < 0) a(:,i) = -a(:,i) 
  end do

  if (dbg>1) then
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
