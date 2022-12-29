
subroutine read_data_files

  ! this routine is part of J. Autschbach's set of programs to process
  ! Molcas data for the generation of various types of spectral
  ! intensities
  
  ! (c) 2019-2022 Jochen Autschbach, SUNY Buffalo

  use definitions
  
  use namelist_module
  
  use constants_parameters
  
  use shared_variables

  implicit none

  complex(KREAL), dimension(:,:,:), allocatable :: quadtmp

  integer(KINT) :: idir, i, j, idum, jdum

  ! ============================================================================

  
  ! ------------------------------------------------------
  ! read the electron spin matrices in a loop over x, y, z
  ! and add to magdip, if applicable
  ! ------------------------------------------------------

  if (havespin .or. haveang) then
    allocate (angmom(nstates,nstates,3))
    if (.not. allocated(magdip)) &
      stop 'read_data_files: magdip not allocated, though it should be'
  end if
 
  if (havespin) then

    angmom = 0
    
    do idir = 1, 3
      
      write (cs,'(a,i1,a)') 'spin-',idir,'.txt'
      
      open(unit=iu_d, file=trim(cs), status='old', iostat=ios)
      if (ios /= 0) then
        write (err,*) 'error: file '//trim(cs)//' does not exist'
        stop 'error termination'
      end if
      
      read(iu_d,*) cstemp

      ! by trial & error, we leared that the spin-y matrix must be
      ! processed the same way as angmon, and spin x,z are processed like
      ! electric dipoles and similar properties. 
      ! n.b. inner loop must be the row index
      do j = 1, nstates
        do i = 1, nstates
          read (iu_d,*, iostat=ios) idum, jdum, ctemp(1:2)
          if (dbg>2) write (out,*) i, j, ctemp(1), ctemp(2)
          if (idir.ne.2) then
            angmom(i,j,idir) = cmplx (ctemp(1), ctemp(2), kind(KREAL))
          else
            angmom(i,j,idir) = cmplx (-ctemp(2), ctemp(1), kind(KREAL))
             !angmom(i,j,idir) = cmplx (ctemp(1), ctemp(2), kind(KREAL)) ! new
          end if
          if (ios /= 0) then
            write (err,*) 'idir, i, j = ', idir, i, j
            write (err,*) 'error reading spin value from '//trim(cs)
            stop 'error termination'
          end if
        end do ! i
      end do ! j
      
      close (iu_d)
    
    end do ! idir = spin components
    
    write (out,'(1x,a/)') 'successfully read spin matrices from files'

    magdip = ge * angmom ! spin contribution to magnetic moment

  end if ! havespin
  
  ! -----------------------------------------------------
  ! read the angular momentum data in a loop over x, y, z
  ! and add to magdip, if applicable
  ! -----------------------------------------------------
  
  if (haveang) then

    angmom = 0
  
    do idir = 1, 3
      
      write (cs,'(a,i1,a)') 'angmom-',idir,'.txt'
      
      open(unit=iu_d, file=trim(cs), status='old', iostat=ios)
      if (ios /= 0) then
        write (err,*) 'error: file '//trim(cs)//' does not exist'
        stop 'error termination'
      end if
      
      read(iu_d,*) cstemp
      
      ! n.b. inner loop must be the row index
      do j = 1, nstates
        do i = 1, nstates
          read (iu_d,*, iostat=ios) idum, jdum, ctemp(1:2)
          if (dbg>2) write (out,*) i, j, ctemp(1), ctemp(2)
          ! angmom is -r x del in Molcas, missing a factor of i
          ! relative to r x p. 
          angmom(i,j,idir) = cmplx (-ctemp(2), ctemp(1), kind(KREAL))
          if (ios /= 0) then
            write (err,*) 'idir, i, j = ', idir, i, j
            write (err,*) 'error reading angmom value from '//trim(cs)
            stop 'error termination'
          end if
        end do ! i
      end do ! j
      
      close (iu_d)
      
    end do ! idir = angular momentum components
    
    write (out,'(1x,a/)') 'successfully read angular momentum data files'

    magdip = magdip + angmom ! add orbital angular momentum to magnetic moment
    
  end if ! haveang
  
  if (havespin .or. haveang) deallocate (angmom)

  ! -----------------------------------------
  ! read velocity data in a loop over x, y, z
  ! -----------------------------------------

  if (havevel) then

    if (.not. allocated(velocity)) &
      stop 'read_data_files: velocity array not allocated, though it should be'

    velocity = 0
    
    do idir = 1, 3
      
      write (cs,'(a,i1,a)') 'velocity_dipole-',idir,'.txt'
      
      open(unit=iu_d, file=trim(cs), status='old', iostat=ios)
      if (ios /= 0) then
        write (err,*) 'error: file '//trim(cs)//' does not exist'
        stop 'error termination'
      end if
      
      read(iu_d,*) cstemp
      
      ! n.b. inner loop must be the row index
      do j = 1, nstates
        do i = 1, nstates
          read (iu_d,*, iostat=ios) idum, jdum, ctemp(1:2)
          !read (iu_d,'(I4,I4,1x,E25.16,1x,E25.16)', iostat=ios) & 
          !  idum, jdum, ctemp(1), ctemp(2)
          if (dbg>2) write (out,*) i, j, ctemp(1), ctemp(2)
          velocity(i,j,idir) = cmplx (-ctemp(2), ctemp(1), kind(KREAL))
          if (ios /= 0) then
            write (err,*) 'idir, i, j = ', idir, i, j
            write (err,*) 'error reading velocity value from '//trim(cs)
            stop 'error termination'
          end if
        end do ! i
      end do ! j
     
      close (iu_d)
      
    end do ! idir = velocity (del) components
    
    write (out,'(1x,a/)') 'successfully read velocity data files'

  end if ! havevel 

  ! --------------------------------------
  ! now process electric multipole moments
  ! --------------------------------------

  ! ---------------------------------------------------------
  ! read the electric dipole data in a loop over x, y, z. The
  ! electric dipole elements include a factor -e = -1 au
  ! ---------------------------------------------------------

  if (havedip) then

    if (.not. allocated(eldip)) &
      stop 'read_data_files: eldip not allocated, though it should be'

    eldip = 0
    
    do idir = 1, 3
      
      write (cs,'(a,i1,a)') 'dipole-',idir,'.txt'
      
      open(unit=iu_d, file=trim(cs), status='old', iostat=ios)
      if (ios /= 0) then
        write (err,*) 'error: file '//trim(cs)//' does not exist'
        stop 'error termination'
      end if
      
      read(iu_d,*) cstemp
      
      ! n.b. inner loop must be the row index
      do j = 1, nstates
        do i = 1, nstates
          read (iu_d,*, iostat=ios) idum, jdum, ctemp(1:2)
          !read (iu_d,'(I4,I4,1x,E25.16,1x,E25.16)', iostat=ios) & 
          !  idum, jdum, ctemp(1), ctemp(2)
          if (dbg>2) write (out,*) i, j, ctemp(1), ctemp(2)
          eldip(i,j,idir) = cmplx (ctemp(1), ctemp(2), kind(KREAL))
          if (ios /= 0) then
            write (err,*) 'idir, i, j = ', idir, i, j
            write (err,*) 'error reading dipole value from '//trim(cs)
            stop 'error termination'
          end if
        end do ! i
      end do ! j
     
      close (iu_d)
      
    end do ! idir = electric dipole components
    
    write (out,'(1x,a/)') 'successfully read electric dipole data files'

  end if ! havedip

  ! -------------------------------------------------------------
  ! read the electric quadrupole data in a loop over x, y, z. The
  ! matrix elements include a factor -e = -1 au
  ! -------------------------------------------------------------
  
  
  if (havequad) then

    if (.not. allocated(elquad)) &
      stop 'read_data_files: elquad not allocated, though it should be'

    elquad = 0
    
    do idir = 1, 6
      
      write (cs,'(a,i1,a)') 'quadrupole-',idir,'.txt'
      
      open(unit=iu_d, file=trim(cs), status='old', iostat=ios)
      if (ios /= 0) then
        write (err,*) 'error: file '//trim(cs)//' does not exist'
        stop 'error termination'
      end if
      
      read(iu_d,*) cstemp
      
      ! n.b. inner loop must be the row index
      do j = 1, nstates
        do i = 1, nstates
          read (iu_d,*, iostat=ios) idum, jdum, ctemp(1:2)
          !read (iu_d,'(I4,I4,1x,E25.16,1x,E25.16)', iostat=ios) & 
          !  idum, jdum, ctemp(1), ctemp(2)
          if (dbg>2) write (out,*) i, j, ctemp(1), ctemp(2)
          elquad(i,j,idir) = cmplx (ctemp(1), ctemp(2), kind(KREAL))
          if (ios /= 0) then
            write (err,*) 'idir, i, j = ', idir, i, j
            write (err,*) 'error reading quad. value from '//trim(cs)
            stop 'error termination'
          end if
        end do ! i
      end do ! j
      
      close (iu_d)
      
    end do ! idir = quadrupole tensor components
    
    write (out,'(1x,a/)') 'successfully read electric quadrupole data files'

    ! ------------------------------------------------------
    ! form traceless version of the quadrupole operator, use
    ! quadtmp for temp. storage. The matrix does NOT
    ! include the commonly used factor of 3/2, so this may
    ! need to be taken care of in the calling routine.
    ! ------------------------------------------------------

    allocate (quadtmp(nstates,nstates,6))
    
    quadtmp = 0
    
    quadtmp(:,:,:) = elquad(:,:,:)
    
    ! accumulate (1/3)trace in the elquad ZZ components:
    do idir = 1,2
      elquad(:,:,qindex(3,3)) = elquad(:,:,qindex(3,3)) &
        + elquad(:,:,qindex(idir,idir))
    end do
    elquad(:,:,qindex(3,3)) = elquad(:,:,qindex(3,3)) * third
    
    !subtract trace to get  [u v - (1/3) delta(u,v) r^2] with u,v, = x,y,z
    do idir = 1,3
      quadtmp(:,:,qindex(idir,idir)) = &
        quadtmp(:,:,qindex(idir,idir)) - elquad(:,:,qindex(3,3))
    end do
    
    elquad(:,:,:) = quadtmp(:,:,:)
    
    ! test for violations of traceless condition

    quadtmp = 0
    do idir = 1,3
      quadtmp(:,:,1) = quadtmp(:,:,1) + &
        elquad(:,:,qindex(idir,idir))
    end do

    if (maxval(abs(quadtmp(:,:,1))).gt.tiny) then
      write (cstemp,'(1x,a,1x,f15.13)') &
        'WARNING: quadrupole array not traceless within ',tiny
      write (out,*) trim(cstemp)
      write (out,*) maxval(abs(quadtmp(:,:,1)))
      write (err,*) trim(cstemp)
      write (err,*) maxval(abs(quadtmp(:,:,1)))
    end if

    deallocate(quadtmp)

    write (out,'(1x,a/)') 'quadrupoles transformed to traceless form'
    
  end if ! havequad
  
  ! --------------------------------------
  ! print some of the assembled data, if requested

  if (print_m) then
    write (out,'(/1x,a/)') 'magnetic dipole moment matrix elements in au (x,y,z)'
    do i = 1,nstates
      do j = 1,nstates
        write (out,'(i5,1x,i5,1x,3("("F15.10,SP,F15.10,"i)"))') i,j, &
          -half*magdip(i,j,1:3)
      end do
    end do
  end if ! print_m

  if (print_d) then
    write (out,'(/1x,a/)') 'electric dipole moment matrix elements in au (x,y,z)'
    do i = 1,nstates
      do j = 1,nstates
        write (out,'(i5,1x,i5,1x,3("("F15.10,SP,F15.10,"i)"))') i,j, &
          eldip(i,j,1:3)
      end do
    end do
  end if ! print_d

  ! ============================================================================

  return
  
end subroutine read_data_files
