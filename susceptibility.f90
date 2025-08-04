program susceptibility

  ! calculate van Vleck susceptibility from Molcas data

  ! (c) 2025, Jochen Autschbach, SUNY Buffalo

  ! this programs's options are controlled by a Fortran namelist
  ! input.  the namelist is called 'options' and read from a file
  ! called 'options.dat' defined in the constants_parameters module.
  ! Set trange(1:2) in options.dat for the start and end temperature
  ! to be used in this calculation.
  !
  ! use
  !
  ! grep 'T,chi' <outputfile> | awk '{print $2,$3}'
  !
  ! to make a list of T, chi data for plotting.
  ! For further details see the code & comments below.

  use definitions

  use namelist_module

  use constants_parameters

  use shared_variables

  implicit none

  ! conversion parameters
  ! Avogadro number:
  real(KREAL), parameter :: avogadro = 6.02214076E+23_KREAL
  ! Bohr magneton in erg/G = esu cm
  real(KREAL), parameter :: muB = 9.2740100657E-21_KREAL
  ! Hartree energy in erg:
  real(KREAL), parameter :: au2erg = 4.3597447222060E-11_KREAL
  ! Conversion factor for susceptibility
  real(KREAL), parameter :: convchi = (muB**2)*avogadro/au2erg

  complex(KREAL) :: chi, chitmp

  integer(KINT) :: idir, jdir, i, j

  integer(KINT) ::  ilevel, jlevel, is, js

  real(KREAL) :: kTm1, q0, ei, ej, bfac, chitens(3,3)

  ! in-line functions

  !real(KREAL) :: waveno, evolt
  !waveno(rtemp) = rtemp * au2cm
  !evolt(rtemp) = rtemp * au2ev

  ! ============================================================================

  write (out,'(/1x,a/)') 'Van Vleck Susceptibility Program'

  if (dbg>0) then
    write (out,'(/1x,a/)') 'namelist default values'
    write (out,options)
  end if

  ! N.B. the namelist variables set below are defined in namelist-module.F90

  ! define if we want to use spin, angmom, quadrupole

  nospin    = .false.
  noangmom  = .false.
  nodip     = .true.
  noquad    = .true.

  print_d = .false. ! options for printing detailed data
  print_m = .false. ! for analysis purposes
  print_q = .false.

  magdiag = .false. !  require a diagonalization of the magnetic
  ! moment operator in the basis of the degenerate ground state
  ! components.

  ! default delta criterion for degeneracy, in au
  ddelta = 1E-5_KREAL

  ! options below will lead to a crash, this is to make sure
  ! that actual data for the electronic states are given in the input
  ! file:

  degen = 0 ! this can stay 0 if ddelta is not zero
  nstates = 0
  temp = zero
  skip = 0
  theta = .false.

  ! ----------------------------------------------
  ! read namelist input for molecule-specific data
  ! ----------------------------------------------

  call read_options

  if (dbg>0) then
    write (out,'(/1x,a/)') 'namelist values after prog. startup and input'
    write (out,options)
  end if

  ! sanity checks:
  ! the 'temp' valiable will be ignored. Instead, we use the 'trange' array:
  ! trange(1:2) = start-T, end-T (we use 1K increments)

  if (nstates.le.degen) then
    close (iu_op)
    stop &
      'nstates < degen in namelist input'
  end if

  if (nospin .and. noangmom) then
    stop 'no spin and no angmom. what am I supposed to do ???'
  end if

  do_group = (ddelta.gt.zero) ! if we should group
  if (degen.lt.1 .and. .not.do_group) stop 'degen < 1 requires ddelta > 0'

  ! define temperatures at which to calculate the susceptibility
  if (trange(1).le.zero) stop 'T(start) <= 0. Aborting'
  if (trange(2).lt.trange(1)) stop 'T(end) < T(start. Aborting'

  ! we're defining these for convenience:
  havespin = .not.nospin    ! we will use spin matrices
  haveang  = .not.noangmom  ! we will use angular momentum matrices

  ! we only need magnetic moments in this code, so let's over-ride any
  ! options for electric moments that may be in options.dat

  havedip = .false.
  havequad = .false.
  if (.not.nodip) write(out,*) 'NODIP=.F. option IGNORED in input'
  if (.not.noquad) write(out,*) 'NOQUAD=.F. option IGNORED in input'

  ! print header:

  write (out, '(/1x,a/1x,a,1x,f7.2,a,f7.2/)') &
    'Susceptibility calculation','T =',trange(1),' --',trange(2)

  if (nospin) write (out,'(1x,a/)') '*** not including S contributions'
  if (noangmom) write (out,'(1x,a/)') '*** not including L contributions'

  ! -----------------------------------------
  ! allocate memory other than scratch arrays
  ! -----------------------------------------

  allocate (energy(nstates))

  allocate (deglist(nstates))
  allocate(levels(nstates), accl(nstates), elevel(nstates))

  ! --------------------------------------------
  ! read the energies from file and process them
  ! --------------------------------------------

  call process_energies (nstates, nlevels, &
    do_group, ddelta, skip, degen, energy, elevel, &
    deglist, levels, accl)

  ! at this point, we have possibly re-defined the variables degen and
  ! skip to match the state grouping. However, if the GS is degenerate
  ! but not grouped properly the calculation will probably fail.

  ! -----------------------------------------------------------------
  ! read transition moment data from the data files, and assemble the
  ! magnetic moment operator matrix elements from (with u = x,y,z)
  ! L_u + 2 S_u. We also attach the factor
  ! -e\hbar / (2 m_e) = -1/2 au.
  ! -----------------------------------------------------------------

  allocate (magdip(nstates,nstates,3))
  magdip = 0

  call read_data_files

  ! magdip = -half * magdip ! not used so we can use the L+2S formula

  ! ----------------------------
  ! Main part of the computation
  ! ----------------------------

  write (out,'(//1x,40(''-'')/1x,a/)') 'running Susceptibility calculation ...'

  ! ---------------------------------------------
  ! define initial temperature T and then execute
  ! a loop over T in increments of 1 kelvin
  ! ---------------------------------------------

  temp = trange(1)

  do while (temp < trange(2))

    ! Boltzmann factor
    kT = temp * boltzau
    if (kT.le.zero) stop 'kT < 0. Aborting'
    kTm1 = one/kT

    chi = c0 ! initialize isotropic suscep. with complex zeros
    chitens = zero

    ! -----------------------------------------------------
    ! start loop over the components of the B-field components
    ! -----------------------------------------------------

    q0 = zero ! initialize partition function

    ! -------------------------------
    ! loop over all electronic levels
    ! -------------------------------

    do ilevel = 1,nlevels

      ei = elevel(ilevel)

      ! calculate Boltzmann factor for level no. i
      rtemp = kTm1 * ei
      if (rtemp.lt.zero) stop 'ei/kT < 0! aborting'
      if (rtemp.lt.explim) then
        bfac = exp(-rtemp)
      else
        bfac = zero
      end if

      ! cycle this loop if the Boltzmann factor is effectively zero
      if (bfac.lt.tiny) then
        if (dbg>1) write(out,*) ilevel,'ilevel loop cycled (small B)'
        cycle
      end if

      ! add Boltzmann factor times degeneracy to the partition function
      q0 = q0 + bfac * float(levels(ilevel))

      ! -------------------------------------------
      ! Curie term:
      ! double loop over components of level ilevel
      ! -------------------------------------------

      if (dbg>0) then
        write(out,*) 'Temperature (K): ',temp
        write(out,*) 'Level ',ilevel,' has ',levels(ilevel),' states'
        write (out,*) '  Boltzmann factor:', bfac
      end if

      do i = 1, levels(ilevel)

        is = accl(ilevel) + i ! refers to un-grouped set of states

        if (is.gt.nstates .or. is.lt.1) &
          stop 'is out of bounds'

        do j = 1,levels(ilevel)

          js = accl(ilevel) + j ! refers to un-grouped set of states

          if (js.gt.nstates .or. js.lt.1) &
            stop 'js out of bounds'

          do idir = 1,3

            chitmp = magdip(is,js,idir) * magdip(js,is,idir)

            ! output debug code if the previous product is not purely
            ! real, as it is supposed to be
            if (abs(aimag(chitmp)).gt.small) then
              write(out,*) 'Curie is,js,idir,chitmp',is,js,idir,chitmp
              write(out,*) magdip(is,js,idir), magdip(js,is,idir)
            end if

            chitmp = chitmp * kTm1
            chi = chi + bfac * chitmp

            do jdir=1,3
              chitmp = magdip(is,js,idir) * magdip(js,is,jdir) * kTm1
              chitens(idir,jdir) = chitens(idir,jdir)  + bfac * real(chitmp)
            end do ! jdir

          end do ! idir

        end do ! j
      end do ! i

      ! --------------------------
      ! SOS term:
      ! loop over all other levels
      ! --------------------------

      do jlevel = 1,nlevels

        if (jlevel.eq.ilevel) cycle

        ej = elevel(jlevel)

        deltae = ej - ei

        if (abs(deltae).lt.small) then
          write (out,*) &
            'ilevel,jlevel,E(j-i):',ilevel,jlevel,deltae
          stop 'Delta E too small in TIP/SOS term'
        end if

        ! nested loop over components of level i
        ! and components of level j

        do i = 1,levels(ilevel)

          is = accl(ilevel) + i ! refers to un-grouped set of states

          if (is.gt.nstates .or. is.lt.1) &
            stop 'is out of bounds'

          do j = 1,levels(jlevel)

            js = accl(jlevel) + j ! refers to un-grouped set of states

            if (js.gt.nstates .or. js.lt.1) &
              stop 'js out of bounds'

            do idir = 1,3

              chitmp = magdip(is,js,idir) * magdip(js,is,idir)

              ! output debug code if the previous product is not purely
              ! real, as it is supposed to be
              if (abs(aimag(chitmp)).gt.small) then
                write(out,*) 'SOS is,js,idir,chitmp',is,js,idir,chitmp
                write(out,*) magdip(is,js,idir), magdip(js,is,idir)
              end if

              chitmp = chitmp / deltae
              chi = chi + two * bfac * chitmp

              do jdir=1,3
                chitmp = magdip(is,js,idir) * magdip(js,is,jdir) / deltae
                chitens(idir,jdir) = chitens(idir,jdir) + two* bfac * real(chitmp)
              end do ! jdir

            end do ! idir

          end do ! j

        end do ! i

      end do ! jlevel

    end do ! ilevel

    ! --------------------------------
    ! done loop over electronic levels
    ! --------------------------------

    write(out,*) 'Temperature (K): ', temp
    write(out,*) 'Partition Function: ', q0
    if (abs(q0).gt.small) then
      chi = third * chi / q0
      chitens(:,:) = chitens(:,:) / q0
    else
      write(out,*) real(chi), aimag(chi)
      stop 'partition function too small! aborting.'
    end if

    write(out,*) 'Susceptibility results (Re/Im):'
    write(out,'(1x,a,f6.2,1x,f15.8)') 'T,chi(T): ', temp, real(chi)*convchi
    write(out,*)

    write(out,*) 'Susceptibility tensor:'
    write(out,*) chitens(1:3,1)*convchi
    write(out,*) chitens(1:3,2)*convchi
    write(out,*) chitens(1:3,3)*convchi
    rtemp = zero
    do idir = 1,3
      rtemp = rtemp + chitens(idir,idir)
    end do
    write(out,*) 'Average: ',rtemp*third*convchi

    write(out,*)

    ! increment temperature by 1 K:
    temp = temp + one

  end do ! while temp < trange(2)

  ! ---------------------------
  ! done loop over temperatures
  ! ---------------------------

  call print_constants

  ! -----------------
  ! deallocate arrays
  ! -----------------

  deallocate(energy, magdip, deglist, levels, elevel, accl)

  ! --------
  ! all done
  ! --------

  stop 'normal termination of Susceptibility code'

  ! ============================================================================


end program susceptibility
