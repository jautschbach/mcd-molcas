program mcd_a_molcas

  ! calculate MCD A-term intensities from Molcas data

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


  use definitions

  use namelist_module

  use constants_parameters
  
  use shared_variables
  
  implicit none
  
  complex(KREAL), dimension(:), allocatable :: ctav, clist

  complex(KREAL) :: ct, cttmp, prefac, vec(3,3)

  integer(KINT) :: idir, jdir, i1, i2, j1, j2
  
  integer(KINT) ::  ilevel, jlevel, is1, is2, js1, js2
 
  ! in-line functions

  real(KREAL) :: waveno
  waveno(rtemp) = rtemp * au2cm

  ! ============================================================================

  write (out,'(/1x,a/)') 'MCD A-Term PROGRAM'
  
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

  magdiag = .false. ! 
  ! if true, require that degenerate states diagonalize
  ! the Zeeman operator. The option is ignored in this code even if .true.

  ! default delta criterion for degeneracy, in au
  ddelta = 1E-5_KREAL

  usemag = .false. ! determines the A-term by replacing the electronic dipole with the magnetic dipole.
  !The contribution of the magetic dipole is usually negectable. This has been implemented for all terms
  !by only debugged for the C-term

  ! default options below will lead to a crash, this is to make sure
  ! that actual data for the electronic states are given in the input file:
  
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
    'MCD A-term calculation','T =',temp,'kT in cm**(-1)/K:',kT
  
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

  write (out,'(//1x,40(''-'')/1x,a/)') 'A-term calculation'
  
  ! memory allocations:

  
  allocate (ctav(nlevels), clist(nlevels))

  ! in Piepho and Schatz, p. 88, 
  ! the prefactor for the A-term is +i/(degeneracy of ground state)
  ! However, in my notes I changed that around so that we can keep
  ! the same pre-factor as in the C term code
  
  prefac = -sqm1 / (degen)
  if (dbg>0) write (out,*) 'prefac = ', prefac
  
  ctav(:) = c0 ! isotropic A-term

  ! create output files, open, and write namelist
  ! input for plot program

  do idir = 0,3
    
    write(outfile(idir),'(a,i1)') 'mcd-a-spectrum-',idir
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
      'nexcit=',ntemp,', invert=F, waveno=T, term=''A'', temp=',temp,&
      ', theta=',theta,' /', &
      '#  E(cm**-1), Re-A (au), Re-A (D**2), Im-A (au), Im-A (D**2), f'

  end do ! idir
  
  ! -----------------------------------------------------
  ! start loop over the components of the B-field components
  ! -----------------------------------------------------

  do idir = 1,3

    write (out,'(1x,a,1x,i1/)') 'B-field direction',idir
    
    
    write (out,'(/1x,a,1x,a/1x,a/1x,a/)') 'A-terms for 0 -> f',&
      & 'in au and  Debye**2' , 'The&
      & Im part should be zero and is printed for debug purposes',&
      & 'The data will be written to file '//trim(outfile(idir))
    
    clist(:) = c0 ! initialize A-terms with complex zero
    
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

      ct = c0

      ! ------------------------------
      ! double loop over EL components
      ! ------------------------------
      
      do j1 = 1,levels(jlevel)
        
        js1 = accl(jlevel) + j1 ! refers to un-grouped set of states
        
        do j2 = 1, levels(jlevel)
          
          js2 = accl(jlevel) + j2
          
          ! ------------------------------
          ! double loop over GL components
          ! ------------------------------
          
          ilevel = 1
          
          do i1 = 1, levels(ilevel)
            
            is1 = accl(ilevel) + i1
            
            do i2 = 1, levels(ilevel)
              
              is2 = accl(ilevel) + i2
              
              ! make sure we have OK values for the state indices
              
              if (is1.gt.nstates .or. is1.lt.1) &
                stop 'is1 out of bounds'             
              if (is2.gt.nstates .or. is2.lt.1) &
                stop 'is2 out of bounds'              
              if (js1.gt.nstates .or. js1.lt.1) &
                stop 'js1 out of bounds'              
              if (js2.gt.nstates .or. js2.lt.1) &
                stop 'js2 out of bounds'
              
              if (j1.eq.j2) then
                do jdir = 1,3
                  vec(jdir,1) = eldip(is1,js1,jdir)
                  vec(jdir,2) = eldip(js2,is2,jdir)
                  vec(jdir,3) = c0
                end do
                call vector_product_cmplx(vec(:,1),vec(:,2),vec(:,3))
                cttmp = prefac * magdip(is2,is1,idir) * vec(idir,3)          
                ct = ct + cttmp
              end if
              if (i1.eq.i2) then
                do jdir = 1,3
                  vec(jdir,1) = eldip(is1,js1,jdir)
                  vec(jdir,2) = eldip(js2,is2,jdir)
                  vec(jdir,3) = c0
                end do
                call vector_product_cmplx(vec(:,1),vec(:,2),vec(:,3))
                cttmp = prefac * magdip(js1,js2,idir) * vec(idir,3)          
                ct = ct - cttmp
             end if
            end do ! i2
          end do ! i1
        end do ! j2
      end do ! j1
     
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
        waveno(deltae), real(clist(ilevel)), real(clist(ilevel)*(debye**2)), &
        aimag(clist(ilevel)), aimag(clist(ilevel)*(debye**2))
    end do

    
  end do ! idir = magnetic field directions

  ! ----------------------------------------
  ! done loop over magnetic field components
  ! ----------------------------------------

  ! write isotropic spectrum to file

  do ilevel = 2+skip,nlevels
    deltae = elevel(ilevel) - elevel(1)
    write (iu_out(0),'(1x,f14.2,3x,4(f20.8,2x))') &
      waveno(deltae), real(ctav(ilevel)), real(ctav(ilevel)*(debye**2)), &
      aimag(ctav(ilevel)), aimag(ctav(ilevel)*(debye**2))
  end do


  ! --------------------
  ! close mcd data files
  ! --------------------
  
  do idir = 0,3
    close (iu_out(idir))
    write (out,'(1x,a,1x,i7,1x,a,1x,a)') 'wrote A-term data for', &
      nstates-degen-skip,'transitions to file',trim(outfile(idir))  
  end do

  call print_constants
    
  
  ! --------------------------------------------------
  ! deallocate arrays, clean up if necessary, and exit
  ! --------------------------------------------------

  deallocate(energy, eldip, magdip, ctav, &
    deglist, clist, levels, elevel, accl)

  stop 'normal termination of mcd-a'
  
  ! ============================================================================

end program mcd_a_molcas
