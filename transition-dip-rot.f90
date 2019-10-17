program transition_dip_rot

  ! calculate oscillator and rotatory strengths from Molcas data

  ! (c) 2019, Jochen Autschbach, SUNY Buffalo

  ! NOT YET COMPLETE, so some of the description below is aspirational.
  
  ! the dipole (D), magnetic (M), and quadrupole (Q) intensities
  ! are those from S. DeBeer et al. Inorganica Chimica Acta 361 (2008) 965,
  ! and Bernadotte et al., JCP 137 (2012), 204106.
  ! For now, we lack the magnetic quadrupole and electric octupole terms
  ! needed to render the intensities beyond the electric dipole approximation
  ! origin-invariant, and we don't yet have the velocity variants of the
  ! higher order electric multipoles to overcome the origin dependence due
  ! to the basis set incompleteness. I.e. the higher order terms are
  ! calculated as in DeBeer er al., assuming atom-centered transitions
  ! with the atom in question being at the coordinate origin (typically
  ! a metal center in a complex)

  ! note: for the traceless electric quadrupole, we use the same
  ! definition as DeBeer et al., without the commonly used factor of 3/2

  ! note: the electric multipole integrals from Molcas include a factor
  ! of -1 for the electron charge. We therefore generate the magnetic moment
  ! matrix elements also with a -1 factor, unlike in the MCD codes.

  ! this programs's options are controlled by a Fortran namelist
  ! input.  the namelist is called 'options' and read from a file
  ! called 'options.dat' defined in the constants_parameters module.
  ! For further details see the code & comments below.

  use definitions

  use namelist_module

  use constants_parameters

  use shared_variables
  
  implicit none

  complex(KREAL), dimension(:), allocatable :: cdav, crav
  
  complex(KREAL), dimension(:,:), allocatable :: cdlist, crlist
  
  complex(KREAL) :: cd, ctmp, cr, pred, prer

  integer(KINT) :: idir, jdir, kdir, i, j

  integer(KINT) ::  ilevel, jlevel, is, js

  ! in-line functions

  real(KREAL) :: waveno, evolt
  waveno(rtemp) = rtemp * au2cm
  evolt(rtemp) = rtemp * au2ev

  ! ============================================================================

  write (out,'(/1x,a/)') 'Transition Dipole and Rotatory Strengths'
  
  ! debug level:

  dbg = 0

  ! N.B. the namelist variables set below are defined in namelist-module.F90

  ! by default, we want to use spin, angmom, dipole, quadrupole

  nospin    = .false. 
  noangmom  = .false.
  nodip     = .false.
  noquad    = .false.
  
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

  ! Boltzmann factor
  kT = temp * boltzcm
  if (kT.le.zero) stop 'kT < 0. Aborting'

  ! sanity checks: 

  if (temp.le.zero .or. nstates.le.degen .or. &
    (nstates - skip).le.degen .or. skip.lt.0) then
    close (iu_op)
    stop &
      'one or more namelist parameters in options file are unreasonable'
  end if

  if (nospin .and. noangmom .and. nodip .and. noquad) then
    stop 'all NO-operator options set. There is nothing to do ...'
  end if

  if (nospin .and. noangmom) then
    write (out,'(/1x,a,//)') &
      '*** no spin and no angmom requested => no rotatory strengths ***'
  end if

  do_group = (ddelta.gt.zero) ! if we should group 
  if (degen.lt.1 .and. .not.do_group) stop 'degen < 1 requires ddelta > 0'

  if (magdiag) write (out,'(/1x,a,//)') &
    '*** NOTE: magdiag is set to T in options file, but ignored here ***'

  ! we're defining these for convenience:
  havespin = .not.nospin    ! we will use spin matrices
  haveang  = .not.noangmom  ! we will use angular momentum matrices
  havedip  = .not.nodip     ! we will use the el. dipole terms
  havequad = .not.noquad    ! we will use the quadrupole terms

  ! print header:

  write (out, '(/1x,a/1x,a,1x,f7.2/1x,a,1x,f12.6/)') &
    'Dipole and Rotatory Strength calculation','T =',temp,&
    'kT in cm**(-1)/K:',kT

  if (nospin) write (out,'(1x,a/)') '*** not including S contributions'
  if (noangmom) write (out,'(1x,a/)') '*** not including L contributions'
  if (nodip) write (out,'(1x,a/)') '*** not including dipole contributions'
  if (noquad) write (out,'(1x,a/)') '*** not including quad. contributions'

  if (polnotprop) then
    write (out,'(/1x,66(''*'')/4(1x,a/),1x,66(''*''))') &
      'POLNOTPROP option set. The printed oscillator strengths correspond',&
      'to the polarization directions 1,2,3 = x,y,z and the average (0).',&
      'The rotatory strengths in files 1,2,3 are the products x-x, y-y,',&
      ' and z-z of the components of the electric and magnetic TDM.'
  end if


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
  ! We operate under the assumption that the electric dipole and
  ! quadrupoles include -e = -1 au factors.  As we are dealing with
  ! transition dipoles, there is no need to remove any nuclear
  ! contributions to the moments.  Upon return from read_data, the
  ! quadrupole is traceless.
  ! -----------------------------------------------------------------
  
  if (havedip) then
    allocate (eldip(nstates,nstates,3))
    eldip = 0
  end if
  if (havespin .or. haveang) then
    allocate (magdip(nstates,nstates,3))
    magdip = 0
  end if
  if (havequad) then
    allocate (elquad(nstates,nstates,6))
    elquad = 0
  end if

  call read_data_files

  if (havespin .or. haveang) then
    magdip = -half * magdip
  end if


  ! ----------------------------
  ! Main part of the computation
  ! ----------------------------
  
  write (out,'(//1x,40(''-'')/1x,a/)') 'D and R calculation'

  ! memory allocations:
  
  allocate (cdav(nlevels), crav(nlevels))
  allocate (cdlist(nlevels,3), crlist(nlevels,3))


  ! common pre-factors for the oscillator and rotatory strength.
  ! Note that for R we pick a factor 3, because in the isotropic
  ! spectrum it gets averaged again, and there we want to get the full
  ! scalar product of mu and m, without a factor 1/3.
  
  pred = cp1 
  prer = cp1 * three
 
  if (dbg>0) write (out,*) 'pred, prer = ', pred, prer

  ! initialize result data arrays with complex-zero:
  
  cdav(:) = c0 ! isotropic D
  crav(:) = c0 ! isotropic R

  cdlist(:,:) = c0 ! directional D
  crlist(:,:) = c0 ! directional R

  ! create output files, open, and write namelist
  ! input for plot program

  do idir = 0,3
    
    write(outfile(idir),'(a,i1)') 'spectrum-',idir
    open (iu_out(idir), file=trim(outfile(idir)), status='unknown', iostat=ios)
    if (ios /= 0) then
      write (err,*) 'problem opening file '//trim(outfile(idir))
      stop 'error with one of the spectrum output files'
    end if
    
    ntemp = nlevels - skip - 1
    !write (out,*) 'ntemp, nlevels, skip', ntemp, nlevels, skip
    if (ntemp.lt.1) stop 'attempting to print data for less than 1 level'
        write(iu_out(idir),'(a,i7,a/,a,i7,a,f7.2,a,l,a/a/a/a)') &
      '&plot nsyme(1)=',ntemp, &
      ', ndegen(1)=1, sigma=1000, sharpen=1, npoints=300,', &
      ' nexcit=',ntemp,', invert=F, waveno=T, temp=',temp,&
      ', theta=',theta,' /', &
      '# oscillator strength f, and rotatory strength R in 1E-40 esu**2 cm**2', &
      '# the last 2 columns should be zero and are printed for debug purposes', &
      '#  E(cm**-1), f, R, Im[f], Im[R]'

  end do ! idir
  
  ! --------------------------------------------------
  ! start loop over the components of the polarization
  ! or propagation directions
  ! --------------------------------------------------

  do idir = 1,3

    write (out,'(1x,a,1x,i1/)') 'Direction',idir
    

    write (out,'(/1x,a,1x,a/1x,a/)') 'D and R for 0 -> f',&
      & 'in 1E-40 cgs units' ,&
      & 'The data will be written to file '//trim(outfile(idir))
    
    if (skip > 0) then
      write (out,'(1x,a,1x,i7,1x,a/)') 'The lowest',skip,'excited levels&
        & will be skipped'
    end if

    ! ---------------------------------------------
    ! loop over excited levels and their components
    ! ---------------------------------------------
    
    skip1 = skip
    
    do jlevel = 2,nlevels
      
      if (skip1.gt.0) then
        skip1 = skip1 - 1
        cycle
      end if
      
      cd = c0 ! set temp variables to complex zero
      cr = c0
      
      do j = 1,levels(jlevel)
        
        js = accl(jlevel) + j ! refers to un-grouped set of states
        
        if (js.gt.nstates .or. js.lt.1) &
          stop 'js1 out of bounds'                                   

        ! ---------------------------------
        ! loop over ground level components
        ! ---------------------------------
        
        ilevel = 1
        
        do i = 1, levels(ilevel)

          deltae = elevel(jlevel) - elevel(ilevel)

          if (dbg>0) write (out,*) 'jlevel,ilevel,delta-E',jlevel,ilevel,deltae
          
          is = accl(ilevel) + i
          
          if (is.gt.nstates .or. is.lt.1) &
            stop 'is1 out of bounds'
          
          ! (a) calculate oscillator strength for polarization idir

          if (havedip) then
            ctmp = pred * two * deltae *  &
              eldip(is,js,idir) * eldip(js,is,idir)
            cd = cd + ctmp
          end if
          
          
          ! (b) calculate rotatory strength matching the
          ! Buckingham-Dunn convention of the OR tensor, diagonal
          ! elements.  We take care of taking the imaginary part of
          ! the electric-magnetic dipole contribution by multiplying
          ! it with -i = -sqm1

          if ( (havespin.or.haveang) .and. havedip) then
            ctmp = prer *  &
              eldip(is,js,idir) * magdip(js,is,idir)
          
            cr = cr - sqm1 * ctmp
          end if

          ! now add the quadrupole part. We take Equation (9) from my
          ! paper ChemPhysChem 12 (2011), 3224, for alpha = beta = idir,
          ! gamma = jdir, delta = kdir,
          ! and adjust for the lack of a factor 3/2 in the quadrupole terms
          ! with an overall prefactor of (1/2)(1/3)(3/2) = 1/4.
          ! For alpha=beta the two quadrupole terms are the same, but
          ! I leave the equation below in case we want to get the whole
          ! rotatory strength tensor.

          if (havequad .and. havedip) then
            do jdir = 1,3
              do kdir = 1,3
                ctmp = prer * fourth * deltae * ( &
                  lc(idir,jdir,kdir) * eldip(is,js,jdir) * &
                  elquad(js,is,qindex(kdir,idir)) + &
                  lc(idir,jdir,kdir) * eldip(is,js,jdir) * &
                  elquad(js,is,qindex(kdir,idir)) )
                cr = cr + ctmp
              end do
            end do
          end if
          
          
        end do ! i1
      end do ! j1
      
      cdlist(jlevel,idir) = cdlist(jlevel,idir) + cd
      !write(out,*) 'cdlist for jlevel now',jlevel,cdlist(jlevel,idir)      

      crlist(jlevel,idir) = crlist(jlevel,idir) + cr
      !write(out,*) 'crlist for jlevel now',jlevel,crlist(jlevel,idir)      
      
    end do ! jlevel
    
    ! -----------------------------
    ! done loop over excited levels
    ! -----------------------------

  end do ! idir = polarization or propagation directions

  ! -------------------------
  ! done loop over directions
  ! -------------------------

  crlist = crlist * rotconv

  ! -----------------------------------------------------------
  ! write spectral data to files. idir = 1,2,3 is for the light
  ! beam propagating in x, y, or z direction.
  ! we use cdav and crav for temp storage.
  ! -----------------------------------------------------------

  ! x-direction

  idir = 1

  if (polnotprop) then
    cdav(:) = cdlist(:,1) 
    crav(:) = crlist(:,1)
  else
    cdav(:) = half * (cdlist(:,2) + cdlist(:,3))
    crav(:) = half * (crlist(:,2) + crlist(:,3))
  end if
    
  do ilevel = 2+skip,nlevels
    deltae = elevel(ilevel) - elevel(1)
    write (iu_out(idir),'(1x,f14.2,3x,4(f20.8,2x))') &
      waveno(deltae), &
      real(cdav(ilevel)),  real(crav(ilevel)), &
      aimag(cdav(ilevel)), aimag(cdav(ilevel))
  end do
  
  close (iu_out(idir))
  write (out,'(1x,a,1x,i7,1x,a,1x,a)') 'wrote spectral data for', &
    nstates-degen-skip,'transitions to file',trim(outfile(idir))  

  ! y-direction

  idir = 2

  if (polnotprop) then
    cdav(:) = cdlist(:,2) 
    crav(:) = crlist(:,2)
  else
    cdav(:) = half * (cdlist(:,1) + cdlist(:,3))
    crav(:) = half * (crlist(:,1) + crlist(:,3))
  end if
  
    
  do ilevel = 2+skip,nlevels
    deltae = elevel(ilevel) - elevel(1)
    write (iu_out(idir),'(1x,f14.2,3x,4(f20.8,2x))') &
      waveno(deltae), &
      real(cdav(ilevel)),  real(crav(ilevel)), &
      aimag(cdav(ilevel)), aimag(cdav(ilevel))
  end do

  close (iu_out(idir))
  write (out,'(1x,a,1x,i7,1x,a,1x,a)') 'wrote spectral data for', &
    nstates-degen-skip,'transitions to file',trim(outfile(idir))  

  ! z-direction

  idir = 3

  if (polnotprop) then
    cdav(:) = cdlist(:,3) 
    crav(:) = crlist(:,3)
  else
    cdav(:) = half * (cdlist(:,1) + cdlist(:,2))
    crav(:) = half * (crlist(:,1) + crlist(:,2))
  end if
    
  do ilevel = 2+skip,nlevels
    deltae = elevel(ilevel) - elevel(1)
    write (iu_out(idir),'(1x,f14.2,3x,4(f20.8,2x))') &
      waveno(deltae), &
      real(cdav(ilevel)),  real(crav(ilevel)), &
      aimag(cdav(ilevel)), aimag(cdav(ilevel))
  end do

  close (iu_out(idir))
  write (out,'(1x,a,1x,i7,1x,a,1x,a)') 'wrote spectral data for', &
    nstates-degen-skip,'transitions to file',trim(outfile(idir))  

  ! isotropic spectrum

  idir = 0
  
  cdav(:) = third * (cdlist(:,1) + cdlist(:,2) + cdlist(:,3))
  crav(:) = third * (crlist(:,1) + crlist(:,2) + crlist(:,3))
    
  do ilevel = 2+skip,nlevels
    deltae = elevel(ilevel) - elevel(1)
    write (iu_out(idir),'(1x,f14.2,3x,4(f20.8,2x))') &
      waveno(deltae), &
      real(cdav(ilevel)),  real(crav(ilevel)), &
      aimag(cdav(ilevel)), aimag(cdav(ilevel))
  end do

  close (iu_out(idir))
  write (out,'(1x,a,1x,i7,1x,a,1x,a)') 'wrote spectral data for', &
    nstates-degen-skip,'transitions to file',trim(outfile(idir))  

  call print_constants

  
  ! --------------------------------------------------
  ! deallocate arrays, clean up if necessary, and exit
  ! --------------------------------------------------
  
  deallocate(energy, cdav, crav, &
    deglist, cdlist, crlist, levels, elevel, accl)

  if (havespin .or. haveang) deallocate(magdip)
  if (havedip) deallocate(eldip)
  if (havequad) deallocate(elquad)
  
  stop 'normal termination of transitions-dip-rot'
  
  ! ============================================================================

  
end program transition_dip_rot
