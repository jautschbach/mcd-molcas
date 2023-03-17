program mcd_b_molcas

  ! calculate MCD B-term intensities from Molcas data

  ! (c) 2018, 2019 Jochen Autschbach, SUNY Buffalo

  ! Implementation according to Piepho & Schatz (1983), pages 79 - 88.

  ! the code, as set up, requires a modified version of Molcas that writes
  ! dipole, spin, and angular momentum matrix elements to files named
  ! dipole-X.txt, spin-X.txt and angmom-X.txt, rescpectively, with
  ! X = 1,2,3 for the Cartesian components x, y, z. We did this for two
  ! reasons: 1. convenience 2. to get machine precision data instead of
  ! parsing the Molcas output.

  ! this programs's options are controlled by a Fortran namelist input.
  ! the namelist is called 'options' and read from a file called
  ! 'options.dat'. For details see the code & comments below.

  ! Note that we use a finite sum over states for the B-terms, and
  ! therefore the B-terms are only approximate. 

  use definitions

  use namelist_module

  use constants_parameters
  
  use shared_variables
  
  implicit none
  
  complex(KREAL), dimension(:), allocatable :: ctav, clist

  complex(KREAL) :: ct, cttmp, prefac, vec(3,3)

  integer(KINT) :: idir, jdir, i, j, k

  integer(KINT) ::  ilevel, jlevel, klevel, is, js, ks

  ! in-line functions

  real(KREAL) :: waveno
  waveno(rtemp) = rtemp * au2cm

  ! ============================================================================

  write (out,'(/1x,a/)') 'MCD B-Term PROGRAM'
  
  if (dbg>0) then
    write (out,*) 'namelist default values'
    write (out,options)
  end if

  ! N.B. the namelist variables set below are defined in namelist-module.F90
  
  ! define if we have spin matrices and angmom and quadrupoles available

  nospin    = .false. 
  noangmom  = .false.
  nodip     = .false.
  noquad    = .true.  

  print_d = .false. ! options for printing detailed data
  print_m = .false. ! for analysis purposes

  magdiag = .false. 
  ! We use the equations from Piepho & Schatz, pp.
  ! 84 - 86, which do not require that degenerate states diagonalize
  ! the Zeeman operator. The option is ignored in this code even if T


  usemag = .false. ! determines the B-term by replacing the electronic dipole with the magnetic dipole.
  !The contribution of the magetic dipole is usually negectable. This has been implemented for all terms
  !by only debugged for the C-term

  ! default delta criterion for degeneracy, in au
  ddelta = 1E-5_KREAL
  
  ! default options below will lead to a crash, this is to make sure
  ! that actual data for the electronic states are given in the input file:
  
  degen = 0 ! this can stay 0 if ddelta is not zero
  nstates = 0
  states_sos = 0
  temp = zero
  skip = 0
  theta=.false.
  
  ! ----------------------------------------------
  ! read namelist input for molecule-specific data
  ! ----------------------------------------------

  call read_options
  
  if (dbg>0) then
    write (out,*) 'namelist values after prog. startup and input'
    write (out,options)
  end if

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

  if (magdiag) write (out,'(/1x,a/)') &
    '*** magdiag set .T. in options file, but not available here. Disabled.'

  if (nospin .and. noangmom .and. nodip .and. noquad) then
    stop 'all NO-operator options set. There is nothing to do ...'
  end if

  if (nospin .and. noangmom) then
    stop 'no spin and no angmom. what am I supposed to do ???'
  end if

  if (nodip) then
    stop 'no electric dipole. what am I supposed to do ???'
  end if

  if (states_sos .eq.0) then
    ! assume the SOS limit of states was not set in the input, or set to 0
    states_sos = nstates
  else if (states_sos .lt. 0) then
    stop 'states_sos variable smaller than 0. aborting'
  else if (states_sos.gt.nstates) then
    stop 'states_sos variable > nstates. aborting'
  end if

  write (out,'(/1x,a,i5/)') 'No. of states used in SOS: ',states_sos

  ! we can let the code figure out the GS degeneracy (if degen = 0) but then
  ! ddelta must be set to a finite value

  do_group = (ddelta.gt.zero) ! if we should group 
  if (degen.lt.1 .and. .not.do_group) stop 'degen < 1 requires ddelta > 0'

  ! we're defining these for convenience:
  havespin = .not.nospin    ! we will use spin matrices
  haveang  = .not.noangmom  ! we will use angular momentum matrices
  havedip  = .not.nodip     ! we will use the el. dipole terms
  havequad = .not.noquad    ! not relevant here

  ! print header:

  write (out, '(/1x,a/1x,a,1x,f7.2/1x,a,1x,f12.6/)') &
    'MCD B-term calculation','T =',temp,'kT in cm**(-1)/K:',kT

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
  ! L_u + 2 S_u. We DO NOT attach the factor
  ! -e\hbar / (2 m_e) = -1/2 au, to correspond to Piepho & Schatz.
  ! We operate under the assumption that the electric dipole elements
  ! include  -e = -1 au factors.  As we are dealing with
  ! transition dipoles, there is no need to remove any nuclear
  ! contributions to the moments. 
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
    write (out,*) &
      'WARNING: noquad=.F. in options file, but quadrupole terms not used.'
    havequad = .false.
  end if

  call read_data_files

  ! note: below we will not further check for havedip and haveang.or.havespin,
  ! as there is nothing to do otherwise, and we have already an error
  ! exit in that situation.

  ! ----------------------------
  ! Main part of the computation
  ! ----------------------------

  write (out,'(//1x,40(''-'')/1x,a/)') 'B-term calculation'
  
  ! memory allocations:

  allocate (ctav(nlevels), clist(nlevels))  
    

  ! in Piepho and Schatz, p. 88, 
  ! the prefactor for the B-term is +2/(degeneracy of ground state)
  ! and then we have to take the imaginary part of the result
  
  prefac = two * cp1 / (degen)
  if (dbg>0) write (out,*) 'prefac = ', prefac
  
  ctav(:) = c0 ! isotropic B-term

  ! create output files, open, and write namelist
  ! input for plot program

  do idir = 0,3
    
    write(outfile(idir),'(a,i1)') 'mcd-b-spectrum-',idir
    open (iu_out(idir), file=trim(outfile(idir)), status='unknown', iostat=ios)
    if (ios /= 0) then
      write (err,*) 'problem opening file '//trim(outfile(idir))
      stop 'error with one of the mcd output files'
    end if
    
    ntemp = nlevels - skip - 1
    if (ntemp.lt.1) stop 'attempting to print data for less than 1 level'
        write(iu_out(idir),'(a,i7,a/1x,a,i7,a,f7.2,a,l,a///a)') &
      '&plot nsyme(1)=',ntemp, &
      ', ndegen(1)=1, sigma=1000, sharpen=1, npoints=300,', &
      'nexcit=',ntemp,', invert=F, waveno=T, term=''B'', temp=',temp,&
      ', theta=',theta,' /', &
      '#  E(cm**-1), Re-B (au), Re-B (D**2), Im-B (au), Im-B (D**2)'

  end do ! idir
  
  ! -----------------------------------------------------
  ! start loop over the components of the B-field components
  ! -----------------------------------------------------
  
  do idir = 1,3
    
    write (out,'(1x,a,1x,i1/)') 'B-field direction',idir
    
    write (out,'(/1x,a,1x,a/1x,a/1x,a/)') 'B-terms for 0 -> f',&
      & 'in au and  Debye**2' , 'The&
      & Im part should be zero and is printed for debug purposes',&
      & 'The data will be written to file '//trim(outfile(idir))
    
    clist(:) = c0 ! initialize B-terms with complex zero
    
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
      
      ct = cmplx(zero, zero, kind(KREAL))
      
      ! ------------------------------
      ! loop over EL components
      ! ------------------------------
      
      do j = 1,levels(jlevel)
        
        js = accl(jlevel) + j ! refers to un-grouped set of states
        
        if (js.gt.nstates .or. js.lt.1) &
          stop 'js out of bounds'
        
        ! -----------------------
        ! loop over GL components
        ! -----------------------
        
        ilevel = 1
        
        do i = 1, levels(ilevel)
          
          is = accl(ilevel) + i
          
          if (is.gt.nstates .or. is.lt.1) &
            stop 'is out of bounds'
          
          ! ---------------------------------------------------
          ! loop over intermediate level k components, incl. GL
          ! this is the SOS in the B-term
          ! ---------------------------------------------------
          
          do klevel = 1, states_sos ! the sum may be limited to < nstates
            
            do k = 1, levels(klevel)
              
              ks = accl(klevel) + k
              
              if (ks.gt.nstates .or. ks.lt.1) &
                stop 'ks out of bounds'
              
              
              ! -----------------
              ! part 1: k .ne. GL
              ! -----------------
              
              if (klevel.ne.ilevel) then
                
                do jdir = 1,3
                  vec(jdir,1) = eldip(is,js,jdir)
                  vec(jdir,2) = eldip(js,ks,jdir)
                  vec(jdir,3) = c0
                end do
                call vector_product_cmplx(vec(:,1),vec(:,2),vec(:,3))
                cttmp = prefac * magdip(ks,is,idir) * vec(idir,3)
                deltae = elevel(klevel) - elevel(ilevel)
                
                if (deltae.lt.small) then
                  write (out,*) &                         
                    'ilevel,jlevel,klevel,E(k-i):',ilevel,jlevel,klevel,deltae
                  stop 'Delta E too small in B-term SOS 1'
                  cttmp = 0
                  deltae = 1
                end if
                ct = ct + cttmp/deltae
                
              end if
              
              ! -----------------
              ! part k: k .ne. EL
              ! -----------------
              
              if (klevel.ne.jlevel) then
                
                do jdir = 1,3
                  vec(jdir,1) = eldip(is,js,jdir)
                  vec(jdir,2) = eldip(ks,is,jdir)
                  vec(jdir,3) = c0
                end do
                call vector_product_cmplx(vec(:,1),vec(:,2),vec(:,3))
                cttmp = prefac * magdip(js,ks,idir) * vec(idir,3)
                deltae = elevel(klevel) - elevel(jlevel)
                
                if (abs(deltae).lt.small) then
                  write (out,*) &                         
                    'ilevel,jlevel,klevel,E(k-i):',ilevel,jlevel,klevel,deltae
                  stop 'Delta E too small in B-term SOS 2'
                end if
                ct = ct + cttmp/deltae 
              end if
            end do ! k
          end do ! klevels
        end do ! i
      end do ! j
      
      
      clist(jlevel) = clist(jlevel) + ct
      
      ctav(jlevel) = ctav(jlevel) + ct/three        
      
    end do ! jlevel
    
    ! -----------------------------
    ! done loop over excited levels
    ! -----------------------------
    
    ! write idir spectrum to file
    
    do ilevel = 2+skip,nlevels
      deltae = elevel(ilevel) - elevel(1)
      write (iu_out(idir),'(1x,f14.2,3x,4(f20.8,2x))') &
        waveno(deltae), aimag(clist(ilevel)), aimag(clist(ilevel)*(debye**2)), &
        real(clist(ilevel)), real(clist(ilevel)*(debye**2))
    end do
    
    
  end do ! idir = magnetic field directions
  
  ! ----------------------------------------
  ! done loop over magnetic field components
  ! ----------------------------------------
  
  ! write isotropic spectrum to file
  
  do ilevel = 2+skip,nlevels
    deltae = elevel(ilevel) - elevel(1)
    write (iu_out(0),'(1x,f14.2,3x,4(f20.8,2x))') &
      waveno(deltae), aimag(ctav(ilevel)), aimag(ctav(ilevel)*(debye**2)), &
      real(ctav(ilevel)), real(ctav(ilevel)*(debye**2))
  end do
  
  
  ! --------------------
  ! close mcd data files
  ! --------------------
  
  do idir = 0,3
    close (iu_out(idir))
    write (out,'(1x,a,1x,i7,1x,a,1x,a)') 'wrote B-term data for', &
      nstates-degen-skip,'transitions to file',trim(outfile(idir))  
  end do

  call print_constants
  
  ! --------------------------------------------------
  ! deallocate arrays, clean up if necessary, and exit
  ! --------------------------------------------------
  
  deallocate(energy, eldip, magdip, ctav, &
    deglist, clist, levels, elevel, accl)
  
  stop 'normal termination of mcd-b'
  
  ! ============================================================================
  
end program mcd_b_molcas
