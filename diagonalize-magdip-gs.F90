
subroutine diagonalize_magdip_gs (idir)

  ! -------------------------------------------------
  ! diagonalize the Zeeman Hamiltonian component idir
  ! within the degenerate GS
  ! -------------------------------------------------
  
  use definitions
  
  use constants_parameters

  use namelist_module
  
  use shared_variables
  
  implicit none

  integer(KINT), intent(in) :: idir

  integer(KINT) :: jdir
  
  complex(KREAL), dimension(:,:), allocatable :: eigv, tmpmat1, tmpmat2
  
  allocate(eigv(degen,degen))
  
  ! the next call uses the complex lapack routine zheevd.
  ! upon return, the input matrix is replaced with the
  ! eigenvectors.
  ! we then transform the blocks 
  
  eigv = 0
  eigv(1:degen,1:degen) = magdip(1:degen,1:degen,idir)
  
  call diagonalize_matrix (degen, eigv)
  
  allocate(tmpmat1(degen, degen))
  tmpmat1 = transpose(conjg(eigv))
  
  
  ! transform the relevant blocks of the magnetic dipole matrix
  
  do jdir = 1,3
    allocate(tmpmat2(degen, nstates))    
    tmpmat2 = matmul(tmpmat1, magdip(1:degen,:,jdir))
    magdip(1:degen,:, jdir) = tmpmat2
    deallocate (tmpmat2)
    allocate (tmpmat2(nstates, degen))
    tmpmat2 = matmul(magdip(:, 1:degen, jdir), eigv)
    magdip(:, 1:degen, jdir) = tmpmat2(:, 1:degen)
    deallocate (tmpmat2)
  end do
  
  call print_rec_matrix(out, degen, real(magdip(1:degen,1:degen,idir)),&
    & 'Transformed Zeeman Hamiltonian REAL part')
  
  call print_rec_matrix(out, degen, aimag(magdip(1:degen,1:degen,idir)),&
    & 'Transformed Zeeman Hamiltonian IMAG part')
  
  ! -------------------------------------------------------
  ! transform the dipole and quad matrix elements between the GS
  ! components and the ESs accordingly, store in *
  ! -------------------------------------------------------
  
  do jdir = 1,3
    allocate(tmpmat2(degen, nstates))         
    tmpmat2 = matmul(tmpmat1, eldip(1:degen,:,jdir))
    eldip(1:degen,:, jdir) = tmpmat2
    deallocate (tmpmat2)      
    allocate (tmpmat2(nstates, degen))      
    tmpmat2 = matmul(eldip(:, 1:degen, jdir), eigv)
    eldip(:, 1:degen, jdir) = tmpmat2(:, 1:degen)
    deallocate (tmpmat2)
  end do ! jdir
  
  if (havequad) then
    do jdir = 1,6
      allocate(tmpmat2(degen, nstates))          
      tmpmat2 = matmul(tmpmat1, elquad(1:degen,:,jdir))
      elquad(1:degen,:, jdir) = tmpmat2
      deallocate (tmpmat2)      
      allocate (tmpmat2(nstates, degen))      
      tmpmat2 = matmul(elquad(:, 1:degen, jdir), eigv)
      elquad(:, 1:degen, jdir) = tmpmat2(:, 1:degen)
      deallocate (tmpmat2)
    end do ! jdir
  end if ! havequad
  
  deallocate(tmpmat1)
  
  deallocate(eigv)

!!$      if (dbg>1) then
!!$      allocate (tmpmat1(nstates, nstates))      
!!$      do jdir = 1,3
!!$        tmpmat1 = transpose(conjg(eldip(:,:,jdir)))
!!$        tmpmat1 = tmpmat1 - eldip(:,:,jdir)        
!!$        write (out,*) jdir, pack(tmpmat1(:,:), abs(tmpmat1(:,:))>tiny)        
!!$      end do ! jdir      
!!$      deallocate (tmpmat1)      
!!$    end if ! dbg
  
end subroutine diagonalize_magdip_gs
