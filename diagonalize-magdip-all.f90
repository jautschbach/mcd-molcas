
subroutine diagonalize_magdip_all (idir)

  ! -------------------------------------------------------
  ! diagonalize the Zeeman Hamiltonian component idir
  ! within all set of degenerate states one after the other
  ! -------------------------------------------------------
  
  use definitions
  
  use constants_parameters

  use namelist_module
  
  use shared_variables
  
  implicit none

  integer(KINT), intent(in) :: idir

  integer(KINT) :: jdir, n, ndeg, nint, nfin
  
  complex(KREAL), dimension(:,:), allocatable :: eigv, tmpmat1, tmpmat2

  do n = 1,nlevels ! loop over all states considering their degeneracy

    ndeg = levels(n) ! number of degeneracy
    nint = accl(n) + 1 ! ref to the first state degenerated in the original list
    nfin = accl(n) + ndeg
    
    if (nint.eq.1) then
      if (ndeg.ne.degen) then
        stop 'behavior different than diagonalize-magdip-gs for the GS!'
      end if
    end if

    if (ndeg.gt.1) then ! if not, there is nothing to do

      if (dbg>0) then
        if (dbg>1) then
          write(out,'(1x,25(''.'')/,/1x,a,1x,i5/)') 'For level', n
        end if
        write(out,'(/1x,a/)') 'diagonalize-magdip-all is working :'
      end if ! dbg

      if ((nfin - nint + 1).ne.ndeg) then
        write(out,*) 'n, nint, nfin, ndeg', n, nint, nfin, ndeg
        stop 'error in diagonalize magdip'
      end if

      ! the next call uses the complex lapack routine zheevd.
      ! upon return, the input matrix is replaced with the
      ! eigenvectors.
      ! we then transform the blocks

      allocate(eigv(ndeg,ndeg))

      eigv = 0
      eigv(1:ndeg,1:ndeg) = magdip(nint:nfin,nint:nfin,idir)

      call diagonalize_matrix (ndeg, eigv)

      allocate(tmpmat1(ndeg, ndeg))
      tmpmat1 = transpose(conjg(eigv))


      ! transform the relevant blocks of the magnetic dipole matrix

      do jdir = 1,3
        allocate(tmpmat2(ndeg, nstates))    
        tmpmat2 = matmul(tmpmat1, magdip(nint:nfin,:,jdir))
        magdip(nint:nfin,:, jdir) = tmpmat2
        deallocate (tmpmat2)
        allocate (tmpmat2(nstates, ndeg))
        tmpmat2 = matmul(magdip(:, nint:nfin, jdir), eigv)
        magdip(:, nint:nfin, jdir) = tmpmat2
        deallocate (tmpmat2)
      end do

      if (dbg>0) then

        write (out,'(/1x,a/)') 'Thus'

        call print_rec_matrix(out, ndeg, real(magdip(nint:nfin,nint:nfin,idir)),&
          & 'Transformed Zeeman Hamiltonian REAL part')

        call print_rec_matrix(out, ndeg, aimag(magdip(nint:nfin,nint:nfin,idir)),&
          & 'Transformed Zeeman Hamiltonian IMAG part')

      end if ! dbg

      ! -------------------------------------------------------
      ! transform the dipole and quad matrix elements between the GS
      ! components and the ESs accordingly, store in *
      ! -------------------------------------------------------

      do jdir = 1,3
        allocate(tmpmat2(ndeg, nstates))         
        tmpmat2 = matmul(tmpmat1, eldip(nint:nfin,:,jdir))
        eldip(nint:nfin,:, jdir) = tmpmat2
        deallocate (tmpmat2)      
        allocate (tmpmat2(nstates, ndeg))      
        tmpmat2 = matmul(eldip(:, nint:nfin, jdir), eigv)
        eldip(:, nint:nfin, jdir) = tmpmat2
        deallocate (tmpmat2)
      end do ! jdir

      if (havequad) then
        do jdir = 1,6
          allocate(tmpmat2(ndeg, nstates))          
          tmpmat2 = matmul(tmpmat1, elquad(nint:nfin,:,jdir))
          elquad(nint:nfin,:, jdir) = tmpmat2
          deallocate (tmpmat2)      
          allocate (tmpmat2(nstates, ndeg))      
          tmpmat2 = matmul(elquad(:, nint:nfin, jdir), eigv)
          elquad(:, nint:nfin, jdir) = tmpmat2
          deallocate (tmpmat2)
        end do ! jdir
      end if ! havequad

      deallocate(tmpmat1)
  
      deallocate(eigv)

    end if ! end of the procdure for one degeneracy
  
  end do ! loop over all degenerated states
  
end subroutine diagonalize_magdip_all
