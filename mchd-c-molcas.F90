program mchd_c_molcas

  ! calculate MChD C-term intensities from Molcas data
  ! (magnetochiral dichroism)

  ! (c) Jochen Autschbach, SUNY Buffalo

  ! Implementation according to equations 6.4.2e and 6.4.2h on page on
  ! page 328 of Barron's Molecular Light Scattering anf Optical
  ! Activity, 2nd edition

  ! from Barron's description it is not clear if the ground state
  ! components have to diagonalize the Zeeman operator, but it very
  ! much looks like it. We keep the option, therefore

  ! code based on mcd-c-molcas.F90 where further information can be
  ! found

  ! this programs's options are controlled by a Fortran namelist
  ! input.  the namelist is called 'options' and read from a file
  ! called 'mcd-options'. For details see the code & comments below.

  use definitions

  use namelist_module
  
  implicit none

  ! numerical constants:
  
  real(KREAL), parameter :: zero=0.0_KREAL
  real(KREAL), parameter :: one=1.0_KREAL
  real(KREAL), parameter :: two=2.0_KREAL
  real(KREAL), parameter :: three=3.0_KREAL
  real(KREAL), parameter :: half=0.5_KREAL
  real(KREAL), parameter :: fourth=0.25_KREAL
  real(KREAL), parameter :: oneby15=(one/15.0_KREAL)
  real(KREAL), parameter :: small=1.0E-5_KREAL
  real(KREAL), parameter :: tiny=1E-10_KREAL

  ! physical constants:
  
  real(KREAL), parameter :: ge=2.00231930436182_KREAL   ! g_e factor
  real(KREAL), parameter :: debye = 2.5417463E0_KREAL   ! au -> Debye
  real(KREAL), parameter :: aucm=2.194746313705E+5_KREAL ! au -> cm**(-1)
  real(KREAL), parameter :: boltzcm=0.69503457_KREAL    ! cm**(-1) / K
  real(KREAL), parameter :: thunit=5.02883E5_KREAL      ! C -> [theta]
  real(KREAL), parameter :: cspeed=137.035999_KREAL     ! speed of light in a.u.
  real(KREAL), parameter :: csq = cspeed*cspeed         ! c**2

  ! allocatable arrays:
  
  real(KREAL), dimension(:), allocatable :: energy, elevel

  integer(KINT), dimension(:), allocatable :: deglist, levels, accl
  
  complex(KREAL), dimension(:,:,:), allocatable :: magdip, magdip_mag, eldip, &
    eldip_mag, angmom, elquad, elquad_mag
  
  complex(KREAL), dimension(:,:), allocatable :: eigv, tmpmat1, tmpmat2, d_tmp

  complex(KREAL), dimension(:), allocatable :: cgav, cglist, caav, calist

  ! other variables:
  
  integer(KINT), parameter :: maxrecs = 1000000, nmax=10000
  
  integer(KINT) :: idir, jdir, kdir, i, j, ios, dbg, idum, jdum, &
    ilevel, jlevel, is1, is2, js1, i1, i2, j1

  integer(KINT) :: qindex(3,3)
  
  character(len=LCHARS) :: cs, cstemp, mcdfile(0:3)
  
  logical :: do_group
  
  integer(KINT), parameter :: out = 6, err = 0, iu_e = 7, iu_d = 8, iu_m = 9, &
    iu_op=5

  integer(KINT) :: iu_mcd(0:3)
  data iu_mcd(0), iu_mcd(1), iu_mcd(2), iu_mcd(3) /10, 11, 12, 13/

  real(KREAL) :: deltae, rtemp, ctemp(2), kT
  integer(KINT) :: ntemp, skip1
  complex(KREAL) :: cg, ctmp, ca, preg, prea, lc(3,3,3)
  logical havespin, haveang, havequad
  

  ! in-line functions

  real(KREAL) :: waveno
  waveno(rtemp) = rtemp * aucm

  ! ============================================================================
  
  ! debug level:

  dbg = 1

  ! N.B. the namelist variables set below are defined in namelist-module.F90

  ! define if we want to use spin, angmom, quadrupole

  nospin    = .false. 
  noangmom  = .false.
  noquad    = .false.
  
  print_d = .false. ! options for printing detailed data
  print_m = .false. ! for analysis purposes
  print_q = .false.

  magdiag = .false. !  require a diagonalization of the magnetic
  ! moment operator in the basis of the degenerate ground state
  ! components.

  ! default delta criterion for degeneracy, in au
  ddelta = 1E-5_KREAL

  ! Not used. Only here so the same mcd-option used by mcd-b-molcas can
  ! be used here without altering:
  states_sos = 0 

  ! default options below will lead to a crash, this is to make sure
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

  open (iu_op, file='mcd-options', status='old', iostat=ios)
    if (ios /= 0) then
    stop 'error: file mcd-options does not exist'
  end if


  read(iu_op,options, iostat=ios)
  if (ios /= 0) then
    close (iu_op)
    stop 'error reading namelist input from file mcd-options'
  end if

  ! sanity checks: 

  if (temp.le.zero .or. nstates.le.degen .or. &
    (nstates - skip).le.degen .or. skip.lt.0) then
    close (iu_op)
    stop &
      'one or more namelist parameters in mcd-options are unreasonable'
  end if

  if (nospin .and. noangmom) then
    close (iu_op)
    stop 'no spin and no angmom. what am I supposed to do ???'
  end if

  ! looks OK so far
  
  close (iu_op)

  ! we can let the code figure out the GS degeneracy (if degen = 0) but then
  ! ddelta must be set to a finite value

  do_group = (ddelta.gt.zero) ! if we should group 
  if (degen.lt.1 .and. .not.do_group) stop 'degen < 1 requires ddelta > 0'

  ! we're defining these for convenience:
  havespin = .not.nospin    ! we will use spin matrices
  haveang  = .not.noangmom  ! we will use angular momentum matrices
  havequad = .not.noquad    ! we will use the quadrupole terms

  ! print header:

  write (out, '(/1x,a/1x,a,1x,f7.2/1x,a,1x,i4/)') &
    'MChD C-term calculation','T =',temp

  if (nospin) write (out,'(1x,a/)') '*** not including S contributions'
  if (noangmom) write (out,'(1x,a/)') '*** not including L contributions'
  if (noquad) write (out,'(1x,a/)') '*** not including quad. contributions'

  ! Boltzmann factor

  kT = temp * boltzcm
  if (kT.le.zero) stop 'kT < 0. Aborting'
  write(out,'(/1x,a,1x,f12.6/)') 'kT in cm**(-1)/K:',kT

  ! Levi-Civita symbol, define here as complex so we can use it later
  ! easily in multiplications with other complex variables without getting
  ! complaints about type mismatches (depening on compiler)

  lc(:,:,:) =  cmplx( zero, zero, kind(KREAL))
  lc(1,2,3) =  cmplx( one,  zero, kind(KREAL))
  lc(1,3,2) =  cmplx(-one,  zero, kind(KREAL))
  lc(2,3,1) =  cmplx( one,  zero, kind(KREAL))
  lc(2,1,3) =  cmplx(-one,  zero, kind(KREAL))
  lc(3,1,2) =  cmplx( one,  zero, kind(KREAL))
  lc(3,2,1) =  cmplx(-one,  zero, kind(KREAL))

  ! index array for the quadrupole operator elements in
  ! packed storage:
  ! XX=1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5, ZZ=6

  qindex(1,1) = 1
  qindex(1,2) = 2
  qindex(2,1) = 2
  qindex(1,3) = 3
  qindex(3,1) = 3
  qindex(2,2) = 4
  qindex(2,3) = 5
  qindex(3,2) = 5
  qindex(3,3) = 6
  
  ! -----------------------------------------
  ! allocate memory other than scratch arrays
  ! -----------------------------------------

  allocate (energy(nstates))

  allocate (deglist(nstates))
  allocate(levels(nstates), accl(nstates), elevel(nstates))

  ! --------------------------------------------
  ! read the energies from file and process them
  ! --------------------------------------------

  call process_energies (iu_e, out, nstates, nlevels, &
    do_group, ddelta, skip, degen, energy, elevel, &
    deglist, levels, accl)

  ! at this point, we have possibly re-defined the variables degen and
  ! skip to match the state grouping. However, if the GS is degenerate
  ! but not grouped properly the calculation will probably fail.

  ! print energies and the level number they belong to

  write (out,'(/1x,a/)') 'Energies relative to state no. 1 and grouping:'
  do i = 1, nstates
    write (out,'(1x,i5,a,1x,f15.8,3x,f15.8,3x,i5)') i,':', &
      energy(i),waveno(energy(i)),deglist(i)
  end do

  write (out,'(/1x,a,1x,i5,1x,a/)') 'There are',nlevels,'grouped levels'
  write (out,*)

  write (out,'(/1x,a/)') 'Levels relative to level no. 1 and degeneracy:'
  do i = 1, nlevels
    write (out,'(1x,i5,a,1x,f15.8,3x,f15.8,3x,i5)') i,':', &
      elevel(i),waveno(elevel(i)),levels(i)
  end do

  write (out,*)

  ! for C-term spectra, degen must be > 1

  if (degen.lt.2) stop &
    'For C-term spectra GS degeneracy must be >= 2. Not detected or set'

  ! allocate some more arrays
  
  allocate (eldip(nstates,nstates,3))
  allocate (magdip(nstates,nstates,3))
  allocate (elquad(nstates,nstates,6))
  ! allocate (d_tmp(3,degen))
  
  
  ! ------------------------------------------------------------
  ! assume that we have spin or angular momentum matrices,
  ! assemble the magnetic moment operator matrix elements as
  ! (L_u + 2 S_u), with dimensionless angular momentum operators
  ! and we exclude the factor -e\hbar / (2 m_e) = -1/2 au.
  ! with that, the magnetic moment is defined the same as in
  ! Piepho & Schatz. We'l keep that in mind further below.
  ! ------------------------------------------------------------

  ! sanity check for options:
  
  if (.not.( &
    (havespin .and. haveang) .or. &
    (havespin .and. .not.haveang) .or. &
    (.not.havespin .and. haveang))) then
    stop 'option combination in magnetic moment assembly not allowed'
  end if

  magdip = 0
  
  allocate (angmom(nstates,nstates,3))
  
  ! ----------------------------------------------------
  ! read the electron spin matrices in a loop over x, y, z
  ! ----------------------------------------------------
  
  angmom = 0

  if (havespin) then
    
    do idir = 1, 3
      
      write (cs,'(a,i1,a)') 'spin-',idir,'.txt'
      
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
          if (idir.ne.2) then
            angmom(i,j,idir) = cmplx (ctemp(1), ctemp(2), kind(KREAL))
          else
            angmom(i,j,idir) = cmplx (-ctemp(2), ctemp(1), kind(KREAL))
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
  
  ! ----------------------------------------------------
  ! read the angular momentum data in a loop over x, y, z
  ! ----------------------------------------------------
  
  angmom = 0

  if (haveang) then
  
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
          ! angmom is -Im, Re in Molcas, missing a factor of -i
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

  ! memory deallocations:
  
  deallocate (angmom)

  ! --------------------------------------
  ! now process electric multipole moments
  ! --------------------------------------

  ! -----------------------------------------------------------
  ! read the electric dipole data in a loop over x, y, z. The
  ! electric dipole elements will exclude the factor -e = -1 au
  ! -----------------------------------------------------------

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

  ! -------------------------------------------------------------
  ! read the electric quadrupole data in a loop over x, y, z. The
  ! matrix elements will exclude the factor -e = -1 au
  ! -------------------------------------------------------------
  
  elquad = 0
  
  if (havequad) then 
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
    
  end if ! havequad
  
  ! --------------------------------------
  ! print the assembled data, if requested

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
          -eldip(i,j,1:3)
      end do
    end do
  end if ! print_d
  
  ! memory allocations:

  allocate (magdip_mag(nstates,nstates,3))
  allocate (eldip_mag(nstates,nstates,3))
  allocate (elquad_mag(nstates,nstates,6))
  allocate (d_tmp(3,degen))
  allocate (eigv(degen,degen))
  allocate (cgav(nlevels), cglist(nlevels),caav(nlevels), calist(nlevels))

  ! ------------------------------------------------------
  ! form traceless version of the quadrupole operator, use
  ! elquad_mag for temp. storage
  ! ------------------------------------------------------

  elquad_mag = 0

  if (havequad) then

    elquad_mag(:,:,:) = three * half * elquad(:,:,:)
    
    ! accumulate (1/2)trace in the elquad ZZ components:
    do idir = 1,2
      elquad(:,:,qindex(3,3)) = elquad(:,:,qindex(3,3)) &
        + elquad(:,:,qindex(idir,idir))
    end do
    elquad(:,:,qindex(3,3)) = elquad(:,:,qindex(3,3)) * half

    !subtract trace to get (1/2) [3 u v - delta(u,v) r^2] with u,v, = x,y,z
    do idir = 1,3
      elquad_mag(:,:,qindex(idir,idir)) = &
        elquad_mag(:,:,qindex(idir,idir)) - elquad(:,:,qindex(3,3))
    end do

    elquad(:,:,:) = elquad_mag(:,:,:)

    ! test for violations of traceless condition

    elquad_mag = 0
    do idir = 1,3
      elquad_mag(:,:,1) = elquad_mag(:,:,1) + &
        elquad(:,:,qindex(idir,idir))
    end do

    if (maxval(abs(elquad_mag(:,:,1))).gt.tiny) then
      write (err,*) 'WARNING: quadrupole array not traceless within 1e-10'
      write (err,*) maxval(abs(elquad_mag(:,:,1)))
    end if

    elquad_mag = 0

  end if ! havequad

  ! ------------------
  ! C-term computation
  ! ------------------
  
  write (out,'(//1x,40(''-'')/1x,a/)') &
    'C-term calculation'

  ! the prefactor for the C(G)-term is 1/(degeneracy of ground state GS)
  ! the prefactor for the C(A)-term is omega/(15 times degeneracy of GS)
  ! We do NOT include the omega here in C(A).
  ! Also, THE C TERMA ARE MULTIPLIED BY c**2 TO AVOID LOSS OF
  ! PRECISION IN THE PRINTED RESULTS. Consider this as absorbing a
  ! mu0/(4pi) from the intensity in the terms
  preg = cmplx(one, zero, kind(KREAL)) / (degen)
  rtemp = oneby15
  prea = cmplx(rtemp, zero, kind(KREAL)) / (degen)
  preg = preg * csq ! here is the c**2 factor
  prea = prea * csq
  if (dbg>0) write (out,*) 'preg, prea = ', preg, prea
  
  cgav(:) = cmplx(zero, zero, kind(KREAL)) ! isotropic C(G)-term
  caav(:) = cmplx(zero, zero, kind(KREAL)) ! isotropic C(A)-term

  ! create output files, open, and write namelist
  ! input for plot program

  do idir = 0,3
    
    write(mcdfile(idir),'(a,i1)') 'mchdspectrum-',idir
    open (iu_mcd(idir), file=trim(mcdfile(idir)), status='unknown', iostat=ios)
    if (ios /= 0) then
      write (err,*) 'problem opening file '//trim(mcdfile(idir))
      stop 'error with one of the mchd output files'
    end if
    
    ntemp = nlevels - skip - 1
    !write (out,*) 'ntemp, nlevels, skip', ntemp, nlevels, skip
    if (ntemp.lt.1) stop 'attempting to print data for less than 1 level'
        write(iu_mcd(idir),'(a,i7,a/,a,i7,a,f7.2,a,l,a/a/a/a)') &
      '# &plot nsyme(1)=',ntemp, &
      ', ndegen(1)=1, sigma=1000, sharpen=1, npoints=300,', &
      '# nexcit=',ntemp,', invert=F, waveno=T, term=''C'', temp=',temp,&
      ', theta=',theta,' /', &
      '# the last 2 columns should be zero and are printed for debug purposes', &
      '# all C terms contain a factor c**2 to avoid loss of printed precision', &
      '#  E(cm**-1), Re-C(G), Im-C(A'')/omega, Im-C(G)/omega, Re-C(A'')'

  end do ! idir
  
  ! -----------------------------------------------------
  ! start loop over the components of the B-field components
  ! -----------------------------------------------------

  do idir = 1,3

    write (out,'(1x,a,1x,i1/)') 'B-field direction',idir
    

    write (out,'(/1x,a,1x,a/1x,a/)') 'MChD C-terms for 0 -> f',&
      & 'in au' ,&
      & 'The data will be written to file '//trim(mcdfile(idir))

    ! we're wasting a bit of memory here by just copying
    ! <array> to <array>_mag.
    ! We may diagonalize Zeeman Hamiltonian in the ground
    ! state (GS), leading to a unitary transformation of the
    ! GS components which then also affected the other dipole and quad.
    ! matrix elements between the GS and the excited states (ESs).
    ! As we may have to do this 3 times, for the different components
    ! of the Zeeman operator, it is convenient to keep the original data
    ! un-transformed in the respective arrays
    
    eldip_mag = eldip
    magdip_mag = magdip
    elquad_mag = elquad
    
    ! use old code if magdiag is set:
    
    if (magdiag) then
      call diagonalize_magdip_gs
      write (out,'(1x,a/)') 'Ground state now diagonalizes Zeeman operator'
      write (out,'(1x,a)') 'Complex magnetic moment matrix elements for GS:'
      do i = 1,degen
        write (out,'(/1x,a,1x,i2)') 'GS component',i
        do jdir = 1,3
          write (out,'(1x,F14.8,SP,F14.8,a)') magdip(i,i,jdir),'i'
        end do
      end do ! i
      write (out,*)
    end if ! magdiag
    
    cglist(:) = cmplx(zero, zero, kind(KREAL)) ! C(G) terms
    calist(:) = cmplx(zero, zero, kind(KREAL)) ! C(A) terms
    
    if (skip > 0) then
      write (out,'(1x,a,1x,i7,1x,a/)') 'The lowest',skip,'excited levels&
        & will be skipped'
    end if

    d_tmp = cmplx(zero, zero, kind(KREAL))

    ! ------------------------
    ! loop over excited levels 
    ! ------------------------

    skip1 = skip
    
    do jlevel = 2,nlevels
      
      if (skip1.gt.0) then
        skip1 = skip1 - 1
        cycle
      end if
      
      cg = cmplx(zero, zero, kind(KREAL))
      ca = cmplx(zero, zero, kind(KREAL))
      
      ! ------------------------------
      ! loop over EL components
      ! ------------------------------
      
      do j1 = 1,levels(jlevel)
        
        js1 = accl(jlevel) + j1 ! refers to un-grouped set of states
        
        if (js1.gt.nstates .or. js1.lt.1) &
          stop 'js1 out of bounds'                                   
        
        ! ----------------------------------------
        ! single or double loop over GL components
        ! depending on whether magdiag=.T. or not
        ! ----------------------------------------
        
        ilevel = 1
        
        do i1 = 1, levels(ilevel)
          
          is1 = accl(ilevel) + i1
          
          if (is1.gt.nstates .or. is1.lt.1) &
            stop 'is1 out of bounds'

          ! Note: We have a combo of
          ! three electron transition moments, and none of them contains the -e
          ! = -1 au factor that is present in all of the operators. Moreover,
          ! the magnetic moment operator is missing the 1/2 from the Bohr
          ! magneton. I.e. the correct result in atomic units requires an
          ! additional factor -1/4 for C(G) and -1/2 for C(A).
          ! We apply these factors below when calculating cg and ca
          
          if (magdiag) then
            
            ! assume the Zeeman operator is diagonal in the GS:

            ! (a) calculate C(G)

            do jdir = 1,3
              do kdir = 1,3
                ctmp = preg * lc(jdir,kdir,idir) * &
                  magdip_mag(is1,is1,idir) * &
                  eldip_mag(is1,js1,jdir) * magdip_mag(js1,is1,kdir)
                cg = cg - fourth * ctmp
              end do
            end do

            ! (b) calculate C(A)

            do jdir = 1,3
              ctmp = prea * magdip_mag(is1,is1,idir) * &
                ( three * eldip_mag(is1,js1,jdir) * &
                elquad_mag(js1,is1,qindex(jdir,idir)) - &
                eldip_mag(is1,js1,idir) * &
                elquad_mag(js1,is1,qindex(jdir,jdir)) )
              ca = ca - half *ctmp
            end do ! jdir
            
          else
            
            ! NEED TO TEST THIS LATER !!!
            
            do i2 = 1, levels(ilevel)
              
              is2 = accl(ilevel) + i2
              
              if (is2.gt.nstates .or. is2.lt.1) &
                stop 'is2 out of bounds'
              
              ! (a) calculate C(A)
              
              do jdir = 1,3
                do kdir = 1,3
                  ctmp = preg * lc(jdir,kdir,idir) * &
                    magdip_mag(is2,is1,idir) * &
                    eldip_mag(is1,js1,jdir) * magdip_mag(js1,is2,kdir)
                  cg = cg - fourth * ctmp
                end do
              end do
              
              ! (b) calculate C(A)
              
              do jdir = 1,3
                ctmp = prea * magdip_mag(is2,is1,idir) * &
                  ( three * eldip_mag(is1,js1,jdir) * &
                  elquad_mag(js1,is2,qindex(jdir,idir)) - &
                  eldip_mag(is1,js1,idir) * &
                  elquad_mag(js1,is2,qindex(jdir,jdir)) )
                ca = ca - half * ctmp
              end do ! jdir
              
            end do ! i2
            
          end if ! magdiag
          
        end do ! i1
      end do ! j1
      
      cglist(jlevel) = cglist(jlevel) + cg
      !write(out,*) 'cglist for jlevel now',jlevel,cglist(jlevel)      
      cgav(jlevel) = cgav(jlevel) + cg

      calist(jlevel) = calist(jlevel) + ca
      !write(out,*) 'calist for jlevel now',jlevel,calist(jlevel)      
      caav(jlevel) = caav(jlevel) + ca
      
    end do ! jlevel
    
    ! -----------------------------
    ! done loop over excited levels
    ! -----------------------------
    
    ! write idir spectrum to file
    
    do ilevel = 2+skip,nlevels
      deltae = elevel(ilevel) - elevel(1)
      write (iu_mcd(idir),'(1x,f14.2,3x,4(f20.8,2x))') &
        waveno(deltae), &
        real(cglist(ilevel)),  aimag(calist(ilevel)), &
        aimag(cglist(ilevel)), real(calist(ilevel))
    end do
    

  end do ! idir = magnetic field directions

  ! ----------------------------------------
  ! done loop over magnetic field components
  ! ----------------------------------------
  
  ! write isotropic spectrum to file
  
  do ilevel = 2+skip,nlevels
    deltae = elevel(ilevel) - elevel(1)
    write (iu_mcd(0),'(1x,f14.2,3x,4(f20.8,2x))') &
              waveno(deltae), &
        real(cgav(ilevel)),  aimag(caav(ilevel)), &
        aimag(cgav(ilevel)), real(caav(ilevel))
  end do
  
  
  ! --------------------
  ! close mcd data files
  ! --------------------
  
  do idir = 0,3
    close (iu_mcd(idir))
    write (out,'(1x,a,1x,i7,1x,a,1x,a)') 'wrote C-term data for', &
      nstates-degen-skip,'transitions to file',trim(mcdfile(idir))  
  end do
  
  ! --------------------------------------------------
  ! deallocate arrays, clean up if necessary, and exit
  ! --------------------------------------------------
  
  deallocate(energy, eldip, eldip_mag, magdip, magdip_mag, &
    elquad, elquad_mag, cgav, caav, &
    deglist, cglist, calist, levels, elevel, accl, d_tmp, eigv)
  
  stop 'normal termination of mchd'
  
  ! ============================================================================


contains

  subroutine diagonalize_magdip_gs

    ! -------------------------------------------------
    ! diagonalize the Zeeman Hamiltonian component idir
    ! within the degenerate GS
    ! -------------------------------------------------

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
    ! note: here and below for the electric moments, we use the copies
    ! <array>_mag and transform all components of the operators, leaving the
    ! original arrays un-touched so we can transform them again with respect to
    ! the diagonalization of another component idir of the Zeeman operator
    
    do jdir = 1,3
      allocate(tmpmat2(degen, nstates))    
      tmpmat2 = matmul(tmpmat1, magdip_mag(1:degen,:,jdir))
      magdip_mag(1:degen,:, jdir) = tmpmat2
      deallocate (tmpmat2)
      allocate (tmpmat2(nstates, degen))
      tmpmat2 = matmul(magdip_mag(:, 1:degen, jdir), eigv)
      magdip_mag(:, 1:degen, jdir) = tmpmat2(:, 1:degen)
      deallocate (tmpmat2)
    end do

    call print_rec_matrix(out, degen, real(magdip_mag(1:degen,1:degen,idir)),&
      & 'Transformed Zeeman Hamiltonian REAL part')

    call print_rec_matrix(out, degen, aimag(magdip_mag(1:degen,1:degen,idir)),&
      & 'Transformed Zeeman Hamiltonian IMAG part')

    ! -------------------------------------------------------
    ! transform the dipole and quad matrix elements between the GS
    ! components and the ESs accordingly, store in *_mag
    ! -------------------------------------------------------

    do jdir = 1,3
      allocate(tmpmat2(degen, nstates))         
      tmpmat2 = matmul(tmpmat1, eldip_mag(1:degen,:,jdir))
      eldip_mag(1:degen,:, jdir) = tmpmat2
      deallocate (tmpmat2)      
      allocate (tmpmat2(nstates, degen))      
      tmpmat2 = matmul(eldip_mag(:, 1:degen, jdir), eigv)
      eldip_mag(:, 1:degen, jdir) = tmpmat2(:, 1:degen)
      deallocate (tmpmat2)
    end do ! jdir

    do jdir = 1,6
      allocate(tmpmat2(degen, nstates))          
      tmpmat2 = matmul(tmpmat1, elquad_mag(1:degen,:,jdir))
      elquad_mag(1:degen,:, jdir) = tmpmat2
      deallocate (tmpmat2)      
      allocate (tmpmat2(nstates, degen))      
      tmpmat2 = matmul(elquad_mag(:, 1:degen, jdir), eigv)
      elquad_mag(:, 1:degen, jdir) = tmpmat2(:, 1:degen)
      deallocate (tmpmat2)
    end do ! jdir

    deallocate(tmpmat1)

  end subroutine diagonalize_magdip_gs
  
end program mchd_c_molcas
