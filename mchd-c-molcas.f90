program mchd_c_molcas

  ! calculate MChD C-term intensities from Molcas data
  ! (magnetochiral dichroism)

  ! (c) 2019, Jochen Autschbach, SUNY Buffalo

  ! Implementation according to equations 6.4.2e and 6.4.2h on page on
  ! page 328 of Barron's Molecular Light Scattering anf Optical
  ! Activity, 2nd edition

  ! this programs's options are controlled by a Fortran namelist
  ! input.  the namelist is called 'options' and read from a file
  ! called 'options.dat' defined in the constants_parameters module.
  ! For further details see the code & comments below.

  use definitions

  use namelist_module

  use constants_parameters

  use shared_variables
  
  implicit none

  complex(KREAL), dimension(:,:,:), allocatable :: eldip_orig, magdip_orig, &
    elquad_orig

  complex(KREAL), dimension(:), allocatable :: cgav, cglist, caav, calist

  complex(KREAL) :: cg, ctmp, ca, preg, prea

  integer(KINT) :: idir, jdir, kdir, i, i1, i2, j
  
  integer(KINT) ::  ilevel, jlevel, is1, is2, js

  real(KREAL) :: kTm1
  
  ! in-line functions

  real(KREAL) :: waveno, evolt
  waveno(rtemp) = rtemp * au2cm
  evolt(rtemp) = rtemp * au2ev

  ! ============================================================================

  write (out,'(/1x,a/)') 'MChD PROGRAM'
  
  ! debug level:

  dbg = 0

  if (dbg>0) then
    write (out,'(/1x,a/)') 'namelist default values'
    write (out,options)
  end if

  ! N.B. the namelist variables set below are defined in namelist-module.F90

  ! define if we want to use spin, angmom, quadrupole

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
  
  if (dbg>0) then
    write (out,'(/1x,a/)') 'namelist values after prog. startup and input'
    write (out,options)
  end if

  ! Boltzmann factor
  kT = temp * boltzau
  if (kT.le.zero) stop 'kT < 0. Aborting'
  kTm1 = one/kT

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
    stop 'no spin and no angmom. what am I supposed to do ???'
  end if

  if (nodip) then
    stop 'no electric dipole. what am I supposed to do ???'
  end if

  do_group = (ddelta.gt.zero) ! if we should group 
  if (degen.lt.1 .and. .not.do_group) stop 'degen < 1 requires ddelta > 0'

  ! we're defining these for convenience:
  havespin = .not.nospin    ! we will use spin matrices
  haveang  = .not.noangmom  ! we will use angular momentum matrices
  havedip  = .not.nodip     ! we will use the el. dipole terms
  havequad = .not.noquad    ! we will use the quadrupole terms

  ! note: below we will not further check for havedip and haveang.or.havespin,
  ! as there is nothing to do in either case, and we have already an error
  ! exit in that situation. 

  ! print header:

  write (out, '(/1x,a/1x,a,1x,f7.2/1x,a,1x,f20.15/)') &
    'MChD C-term calculation','T =',temp,'kT in au/K:',kT

  if (nospin) write (out,'(1x,a/)') '*** not including S contributions'
  if (noangmom) write (out,'(1x,a/)') '*** not including L contributions'
  if (nodip) write (out,'(1x,a/)') '*** not including dipole contributions'
  if (noquad) write (out,'(1x,a/)') '*** not including quad. contributions'


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

  ! for C-term spectra, degen must be > 1

  if (degen.lt.2) stop &
    'For C-term spectra GS degeneracy must be >= 2. Not detected or set'

  
  ! -----------------------------------------------------------------
  ! read transition moment data from the data files, and assemble the
  ! magnetic moment operator matrix elements from (with u = x,y,z)
  ! L_u + 2 S_u. We also attach the factor
  ! -e\hbar / (2 m_e) = -1/2 au. 
  ! We operate under the assumption that the electric dipole and
  ! quadrupoles include -e = -1 au factors.  As we are dealing with
  ! transition dipoles, there is no need to remove any nuclear
  ! contributions to the moments.  Upon return from read_data, the
  ! quadrupole is traceless, and it does NOT contain the factor of
  ! 3/2 that it has in Barron's book. We therefore fix that, too.
  ! -----------------------------------------------------------------
  
  allocate (eldip(nstates,nstates,3))
  eldip = 0
  
  allocate (magdip(nstates,nstates,3))
  magdip = 0
  
  if (havequad) then
    allocate (elquad(nstates,nstates,6))
    elquad = 0
  end if

  call read_data_files

  magdip = -half * magdip

  if (havequad) elquad = elquad * three * half

  
  ! ----------------------------
  ! Main part of the computation
  ! ----------------------------

  write (out,'(//1x,40(''-'')/1x,a/)') 'C-term calculation'
  
  ! if magdiag is set, then save the transition dipole matrices in the
  ! original basis of states in <array>_orig so we can restore the
  ! arrays for each run of idir = 1,3. Not absolutely needed, but convenient.

  if (magdiag) then
    allocate (eldip_orig(nstates,nstates,3))
    allocate (magdip_orig(nstates,nstates,3))
    
    eldip_orig  = eldip
    magdip_orig = magdip
    
    if (havequad) then
      allocate (elquad_orig(nstates,nstates,6))
      elquad_orig  = elquad
    end if
    
  end if
  
  allocate (cgav(nlevels), cglist(nlevels), caav(nlevels), calist(nlevels))

  ! the prefactor for the C(G)-term is 1/(degeneracy of ground state GS)
  ! the prefactor for the C(A)-term is omega/(15 times degeneracy of GS)
  ! We do NOT include the omega here in C(A).
  
  preg = cp1 / degen
  rtemp = oneby15
  prea = rtemp * cp1 / degen
  preg = preg * kTm1 ! division by kT
  prea = prea * kTm1
  if (dbg>0) write (out,*) 'preg, prea = ', preg, prea
  
  cgav(:) = c0 ! initialize isotropic C(G)-terms with complex zeros
  caav(:) = c0 ! initialize isotropic C(A)-term

  ! create output files, open, and write namelist
  ! input for plot program

  do idir = 0,3
    
    write(outfile(idir),'(a,i1)') 'mchd-c-spectrum-',idir
    open (iu_out(idir), file=trim(outfile(idir)), status='unknown', iostat=ios)
    if (ios /= 0) then
      write (err,*) 'problem opening file '//trim(outfile(idir))
      stop 'error with one of the mchd output files'
    end if
    
    ntemp = nlevels - skip - 1
    !write (out,*) 'ntemp, nlevels, skip', ntemp, nlevels, skip
    if (ntemp.lt.1) stop 'attempting to print data for less than 1 level'
        write(iu_out(idir),'(a,i7,a/,a,i7,a,f7.2,a,l,a/a/a/a)') &
      '# &plot nsyme(1)=',ntemp, &
      ', ndegen(1)=1, sigma=1000, sharpen=1, npoints=300,', &
      '# nexcit=',ntemp,', invert=F, waveno=T, term=''C'', temp=',temp,&
      ', theta=',theta,' /', &
      '# the last 2 columns should be effectively zero', &
      '# all C terms contain a factor 1/kT in atomic units', &
      '#  E(cm**-1), Re-C(G), Im-C(A'')/omega, Im-C(G)/omega, Re-C(A'')'

  end do ! idir
  
  ! -----------------------------------------------------
  ! start loop over the components of the B-field components
  ! -----------------------------------------------------

  do idir = 1,3

    write (out,'(1x,a,1x,i1/)') 'B-field direction',idir
    
    write (out,'(/1x,a,1x,a/1x,a/)') 'MChD C-terms for 0 -> f',&
      & 'in au' ,&
      & 'The data will be written to file '//trim(outfile(idir))


    ! select GS components to diagonalize B(idir) Zeeman Hamiltonian if
    ! magdiag is set
    
    if (magdiag) then
      
      eldip = eldip_orig
      magdip = magdip_orig
      if (havequad) elquad = elquad_orig
      
      call diagonalize_magdip_gs(idir)
      
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
    
    cglist(:) = c0 ! initialize C(G) terms with complex zeros
    calist(:) = c0 ! initialize C(A) terms
    
    if (skip > 0) then
      write (out,'(1x,a,1x,i7,1x,a/)') 'The lowest',skip,'excited levels&
        & will be skipped'
    end if

    ! ------------------------
    ! loop over excited levels 
    ! ------------------------

    skip1 = skip
    
    do jlevel = 2,nlevels
      
      if (skip1.gt.0) then
        skip1 = skip1 - 1
        cycle
      end if
     
      cg = c0 ! set to zero
      ca = c0
      
      ! ------------------------------
      ! loop over EL components
      ! ------------------------------
      
      do j = 1,levels(jlevel)
        
        js = accl(jlevel) + j ! refers to un-grouped set of states
        
        if (js.gt.nstates .or. js.lt.1) &
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
          
          if (magdiag) then
            
            ! assume the Zeeman operator is diagonal in the GS:

            ! (a) calculate C(G)

            do jdir = 1,3
              do kdir = 1,3
                ctmp = preg * lc(jdir,kdir,idir) * &
                  magdip(is1,is1,idir) * &
                  eldip(is1,js,jdir) * magdip(js,is1,kdir)
                cg = cg + ctmp
              end do
            end do

            ! (b) calculate C(A)

            do jdir = 1,3
              ctmp = prea * magdip(is1,is1,idir) * &
                ( three * eldip(is1,js,jdir) * &
                elquad(js,is1,qindex(jdir,idir)) - &
                eldip(is1,js,idir) * &
                elquad(js,is1,qindex(jdir,jdir)) )
              ca = ca + ctmp
            end do ! jdir
            
          else
            
            do i2 = 1, levels(ilevel)
              
              is2 = accl(ilevel) + i2
              
              if (is2.gt.nstates .or. is2.lt.1) &
                stop 'is2 out of bounds'
              
              ! (a) calculate C(A)
              
              do jdir = 1,3
                do kdir = 1,3
                  ctmp = preg * lc(jdir,kdir,idir) * &
                    magdip(is2,is1,idir) * &
                    eldip(is1,js,jdir) * magdip(js,is2,kdir)
                  cg = cg + ctmp
                end do
              end do
              
              ! (b) calculate C(A)
              
              do jdir = 1,3
                ctmp = prea * magdip(is2,is1,idir) * &
                  ( three * eldip(is1,js,jdir) * &
                  elquad(js,is2,qindex(jdir,idir)) - &
                  eldip(is1,js,idir) * &
                  elquad(js,is2,qindex(jdir,jdir)) )
                ca = ca + ctmp
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
      write (iu_out(idir),'(1x,f14.4,3x,4(f25.15,2x))') &
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
    write (iu_out(0),'(1x,f14.4,3x,4(f25.15,2x))') &
              waveno(deltae), &
        real(cgav(ilevel)),  aimag(caav(ilevel)), &
        aimag(cgav(ilevel)), real(caav(ilevel))
  end do
  
  ! --------------------
  ! close mcd data files
  ! --------------------
  
  do idir = 0,3
    close (iu_out(idir))
    write (out,'(1x,a,1x,i7,1x,a,1x,a)') 'wrote C-term data for', &
      nstates-degen-skip,'transitions to file',trim(outfile(idir))  
  end do

  call print_constants
  
  ! -----------------
  ! deallocate arrays
  ! -----------------
  
  deallocate(energy, eldip, magdip, cgav, caav, &
    deglist, cglist, calist, levels, elevel, accl)
  if (magdiag) deallocate (eldip_orig, magdip_orig)
  if (havequad) deallocate (elquad)
  if (havequad .and. magdiag) deallocate (elquad_orig)

  ! -----------------
  ! some debug output
  ! -----------------

  if (dbg>0) then
    write (out,'(/1x,a)') 'Levi Civita Tensor nonzero elements'
    do idir = 1,3
      do jdir = 1,3
        do kdir = 1,3
          if (lc(idir,jdir,kdir).ne.c0) &
            write (out,'(1x,3(1x,i2),a,1x,F4.1)') &
            idir, jdir, kdir, ': ', real(lc(idir,jdir,kdir))
        end do
      end do
    end do
    write (out,'(/1x,a)') 'quadrupole indices'
    do idir = 1,3
      do jdir = 1,3
        write (out,'(1x,2(1x,i2),a,1x,i2)') &
          idir, jdir,': ', qindex(idir,jdir)
      end do
    end do
  end if ! dbg

  ! --------
  ! all done
  ! --------
  
  stop 'normal termination of mchd-c'
  
  ! ============================================================================

  
end program mchd_c_molcas
