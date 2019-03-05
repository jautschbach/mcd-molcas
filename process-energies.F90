subroutine process_energies (iu_e, out, nstates, nlevels, do_group, &
  ddelta, skip, degen, energy, elevel, &
  deglist, levels, accl)

  use definitions

  implicit none

  ! subroutine arguments:

  logical, intent(in) :: do_group
  integer(KINT), intent(in) :: iu_e, out, nstates
  real(KREAL), intent(in) :: ddelta
  integer(KINT), intent(inout) :: skip, degen
  
  integer(KINT), intent(out) :: nlevels
  integer(KINT), intent(out) :: levels(nstates), accl(nstates), deglist(nstates)
  real(KREAL), intent(out) :: energy(nstates), elevel(nstates)

  ! numerical constants:
  
  real(KREAL), parameter :: zero=0.0_KREAL
  real(KREAL), parameter :: one=1.0_KREAL
  real(KREAL), parameter :: two=2.0_KREAL
  real(KREAL), parameter :: three=3.0_KREAL
  real(KREAL), parameter :: half=0.5_KREAL
  real(KREAL), parameter :: small=1.0E-5_KREAL
  real(KREAL), parameter :: tiny=1E-10_KREAL

  ! local variables

  integer(KINT) :: i, j, ios
  logical :: found
  character(len=1) :: junk

  ! ===========================================================================

  ! --------------------------------------
  ! read the energies from the energy file
  ! --------------------------------------

  open(unit=iu_e, file='energies.txt', status='old', iostat=ios)
  if (ios /= 0) then
    stop 'error: file energies.txt does not exist'
  end if

  read (iu_e,*) junk
  read (iu_e,*, iostat=ios) energy(1:nstates)
  if (ios /= 0) then
    stop 'error: failure reading energies from energies.txt'
  end if

  close(iu_e)

  write (out,'(1x,a,1x,i4,1x,a/)') &
    'read energies for',nstates,'states from energies.txt'

  if (do_group) write (out,'(1x,a,f10.7/)') &
    'Criterion to detect degeneracies (au):',ddelta

  ! re-set energies relative to E(1) if not already done
  
  do i = 2, nstates
    energy(i) = energy(i) - energy(1)
  end do
  energy(1) = zero

  ! determine degeneracy of levels if grouping is requested

  nlevels = nstates

  if (do_group) then
    deglist(1) = 1
    do i = 2, nstates
      deglist(i) = 0 ! hope this triggers a crash if we don't assign properly
    end do
    
    j = 1
    do i = 2,nstates
      if (abs(energy(i)-energy(i-1)).gt.ddelta) then
        j = j+1
      end if
      deglist(i) = j
    end do

    if (j < nstates) nlevels = j

  else
    do i = 1, nstates
      deglist(i) = i
    end do
  end if

  ! next, we want to assign an array with level information, to
  ! simplify the electronic state loops, and average the energies
  ! within a given level.
  !
  ! accl(i) is the accumulated number of states below level i such
  ! that if we have a loop
  !
  ! j = 1, nlevels
  !   k = 1, levels(nlevels)
  !  ...
  !
  ! then accl(j) + k refers to the level component in the
  ! original state list

    ! write (out,'(1x,a/(1x,i5))') 'deglist:',deglist
  
  j = 1
  levels(1) = 1
  accl(:) = 1
  
  do i = 2, nstates
    
    if (deglist(i).eq.deglist(i-1)) then
      levels(j) = levels(j) + 1
      accl(j) = accl(j) + 1
    else if (deglist(i).gt.deglist(i-1)) then
      j = j + 1
      if (j .gt. nlevels) then
        write (out,*) 'i, nstates, j,nlevels', i, nstates, j,nlevels
        write (out,'(1x,a/(1x,i5))') 'levels:',levels
        write (out,'(1x,a/(1x,i5))') 'accl:',accl
        stop 'count error 1 for j in levels loop'
      end if
      levels(j) = 1
      accl(j) = accl(j-1) + 1
    else
      write (out,*) 'deglist: ',deglist
      stop 'degeneracies list is broken'
    end if
    
  end do
  
  if (j.ne.nlevels) stop 'count error 2 for j in levels loop'

  !write (out,'(1x,a/6(1x,i5))') 'accl:',accl


  ! clean up from the grouping. We must still deal with a degeneracy
  ! of the GS if it was set by input, no matter if state grouping was
  ! requested. Also, if we let the code determine degen we have to take
  ! care of that.

  ! if degen wasn't set by input we have to set it here based
  ! on the calculated state grouping. This is the GS degeneracy:
  
  if (do_group .and. degen.lt.1) then
    degen = levels(1)
  end if

  if (do_group) then
    
    ! accl(i) is currently the number of states up to an including
    ! level no. i. We want to change that so that levels(i) is
    ! not included:
    
    if (degen.gt.0 .and. levels(1).ne.degen) stop &
      'count error levels(1) vs degen'
    if (accl(1).ne.degen) stop 'count error acc(1) vs degen'
    do i = 1,nlevels
      accl(i) = accl(i) - levels(i)
    end do

  else if (.not.do_group .and. degen > 0) then

    ! no grouping, and degen defined manually: set the levels
    ! variables accordingly:
    
    levels(1) = degen
    nlevels = nstates - degen + 1
    accl(1) = 0
    do j = 1,degen
      deglist(j) = 1
    end do
    do j = degen+1,nstates
      deglist(j) = j - degen + 1
    end do
    do j = 2,nlevels
      accl(j) = degen + j -2
    end do
    
  end if ! do_group or not

  !write (out,*) 'nlevels = ',nlevels
  !write (out,'(1x,a/6(1x,i5))') 'accl:',accl(1:nlevels)
  !write (out,'(1x,a/6(1x,i5))') 'levels:',levels(1:nlevels)

  if (accl(2).ne.degen  .and. ddelta >0) stop 'count error acc(2) vs degen'

  write(out,'(/1x,a,i5/)') 'Ground state degeneracy: ', degen

  ! now assign the energies per grouped level. While we're at it,
  ! the energies within a grouped level are averaged:

  if (do_group) then
    do i = 1, nlevels
      elevel(i) = zero
      do j = 1, levels(i)
        elevel(i) = elevel(i) + energy(accl(i)+j)
      end do
      elevel(i) = elevel(i) / real(levels(i))
    end do
  else
    elevel(1) = zero
    do i = 1, degen
      elevel(1) = elevel(1) + energy(i)
    end do
    elevel(1) = elevel(1) / real(degen)
    if (elevel(1).lt.zero) stop 'elevel(1) error'
    do i = 2,nlevels
      elevel(i) = energy(accl(i)+1)
    end do
  end if

  ! make sure the number of skipped states is compatible with
  ! the state grouping

  if (skip .gt. 0) then
    found = .false.
    do i =2,nlevels
      if (skip+degen.eq.accl(i)) then
        found = .true.
        exit
      end if
    end do
    if (found) then
      write (out,'(/1x,a,1x,i5,1x,a/1x,a,1x,i5,1x,a/)') 'Will skip',&
        skip,'excited states,','corresponding to', &
        i-2,'levels above the ground level'
      skip = i-2
    else
      stop 'skip parameter not compatible with state grouping'
    end if
  end if ! skip states?

  
  ! all done

  ! ===========================================================================
  
  return

end subroutine process_energies
