program mcd_b_molcas

  ! calculate MCD B-term intensities from Molcas data

  ! (c) Jochen Autschbach, SUNY Buffalo

  ! Implementation according to Piepho & Schatz (1983), pages 79 - 88.

  ! the code, as set up, requires a modified version of Molcas that writes
  ! dipole, spin, and angular momentum matrix elements to files named
  ! dipole-X.txt, spin-X.txt and angmom-X.txt, rescpectively, with
  ! X = 1,2,3 for the Cartesian components x, y, z. We did this for two
  ! reasons: 1. convenience 2. to get machine precision data instead of
  ! parsing the Molcas output.

  ! this programs's options are controlled by a Fortran namelist input.
  ! the namelist is called 'options' and read from a file called
  ! 'mcd-options'. For details see the code & comments below.

  ! Note that we use a finite sum over states for the B-terms, and
  ! therefore the B-terms are only approximate. 

  use definitions

  use namelist_module
  
  implicit none

  ! numerical constants:
  
  real(KREAL), parameter :: zero=0.0_KREAL
  real(KREAL), parameter :: one=1.0_KREAL
  real(KREAL), parameter :: two=2.0_KREAL
  real(KREAL), parameter :: three=3.0_KREAL
  real(KREAL), parameter :: half=0.5_KREAL
  real(KREAL), parameter :: small=1.0E-5_KREAL
  real(KREAL), parameter :: tiny=1E-10_KREAL

  ! physical constants:
  
  real(KREAL), parameter :: ge=2.00231930436182_KREAL       ! g_e factor
  real(KREAL), parameter :: debye = 2.5417463E0_KREAL       ! au -> Debye
  real(KREAL), parameter :: aucm=2.194746313705E+5_KREAL    ! au -> cm**(-1)
  real(KREAL), parameter :: boltzcm=0.69503457_KREAL        ! cm**(-1) / K
  real(KREAL), parameter :: thunit=5.02883E5_KREAL          ! C -> [theta]
  real(KREAL), parameter :: cspeed=137.035999_KREAL         ! speed of light in a.u. 

  ! allocatable arrays:
  
  real(KREAL), dimension(:), allocatable :: energy, elevel

  integer(KINT), dimension(:), allocatable :: deglist, levels, accl
  
  complex(KREAL), dimension(:,:,:), allocatable :: magdip, eldip, spin, &
    angmom(:,:,:)
  
  complex(KREAL), dimension(:), allocatable :: ctav, clist

  ! other variables:
  
  integer(KINT), parameter :: maxrecs = 1000000, nmax=10000
  
  integer(KINT) :: idir, jdir, ios, dbg, idum, jdum, &
    ilevel, jlevel, i,j, is, js, klevel, ks, k
  
  character(len=LCHARS) :: cs, cstemp, mcdfile(0:3)
  
  logical :: do_group
  
  integer(KINT), parameter :: out = 6, err = 0, iu_e = 7, iu_d = 8, iu_m = 9, &
    iu_op=5

  integer(KINT) :: iu_mcd(0:3)
  data iu_mcd(0), iu_mcd(1), iu_mcd(2), iu_mcd(3) /10, 11, 12, 13/

  real(KREAL) :: deltae, rtemp, ctemp(2), kT
  integer(KINT) :: ntemp, skip1
  complex(KREAL) :: ct, cttmp, prefac, vec(3,3)
  logical havespin, haveang

  ! in-line functions

  real(KREAL) :: waveno
  waveno(rtemp) = rtemp * aucm

  ! ============================================================================
  
  ! debug level:

  dbg = 1

  ! N.B. the namelist variables set below are defined in namelist-module.F90
  
  ! define if we have spin matrices and angmom and quadrupoles available

  nospin    = .false. 
  noangmom  = .false.
  noquad    = .true.  

  print_d = .false. ! options for printing detailed data
  print_m = .false. ! for analysis purposes

  magdiag = .false. ! the original version of this code required a
  ! diagonalization of the magnetic moment operator in the basis of 
  ! the degenerate ground state components.
  ! The default is now to use the equations from Piepho & Schatz, pp.
  ! 84 - 86, which bypass that diagonalization.

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

  ! print header:

  write (out, '(/1x,a/)') &
    'MCD B-term calculation'

  if (nospin) write (out,'(1x,a/)') '*** not including S contributions'
  if (noangmom) write (out,'(1x,a/)') '*** not including L contributions'

  ! Boltzmann factor

  kT = temp * boltzcm
  if (kT.le.zero) stop 'kT < 0. Aborting'
  write(out,'(/1x,a,1x,f12.6/)') 'kT in cm**(-1)/K:',kT
  
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

  ! allocate some more arrays
  
  allocate (eldip(nstates,nstates,3))
  allocate (magdip(nstates,nstates,3))
  ! allocate (d_tmp(3,degen))

  ! -------------------------------------------------------
  ! assume that we have spin and angular momentum matrices.
  ! assemble the magnetic moment operator from those.
  ! -------------------------------------------------------
  
  
  allocate (spin(nstates,nstates,3), angmom(nstates,nstates,3))
  
  ! ----------------------------------------------------
  ! read the electron spin matrices in a loop over x, y, z
  ! ----------------------------------------------------
  
  spin = 0
  
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
          spin(i,j,idir) = cmplx (ctemp(1), ctemp(2), kind(KREAL))
        else
          spin(i,j,idir) = cmplx (-ctemp(2), ctemp(1), kind(KREAL))
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
  
  ! ----------------------------------------------------
  ! read the angular momentum data in a loop over x, y, z
  ! ----------------------------------------------------
  
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

  ! -----------------------------------------------------------
  ! We use the Piepho & Schatz unit conventions, so the magnetic
  ! matrix elements we need are (L_u + 2 S_u), or use g_e S_u
  ! -----------------------------------------------------------
  
  do idir = 1,3
    do i = 1,nstates
      do j = 1,nstates
        if (havespin .and. haveang) then
          magdip(i,j,idir) = one* ( &
            angmom(i,j,idir) + ge * spin(i,j,idir) )
        else if (havespin .and. .not.haveang) then
          magdip(i,j,idir) = one* ( &
            ge * spin(i,j,idir) )
        else if (.not.havespin .and. haveang) then
          magdip(i,j,idir) = one* ( &
            angmom(i,j,idir)  )
        else
          stop 'option combination in magnetic moment assembly not allowed'
        end if
      end do
    end do
  end do

  
  ! ----------------------------------------------------
  ! read the electric dipole data in a loop over x, y, z
  ! ----------------------------------------------------

  eldip = 0

  if (usemag) then
     do idir = 1, 3
        ! n.b. inner loop must be the row index
        do j = 1, nstates
           do i = 1, nstates
              eldip(i,j,idir) =  -(half/cspeed)* ( &
                   magdip(i,j,idir) )
           end do ! j
        end do ! i
     end do ! idir
     
  else 
     
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
  end if
  write (out,'(1x,a/)') 'successfully read electric dipole data files'
  

  
  
  ! memory de/allocations:
  
  deallocate (angmom, spin)
  allocate (ctav(nlevels), clist(nlevels))  
    

  ! ------------------
  ! B-term computation
  ! ------------------
  
  write (out,'(//1x,40(''-'')/1x,a/)') &
    'B-term calculation'

  ! in Piepho and Schatz, p. 88, 
  ! the prefactor for the B-term is +2/(degeneracy of ground state)
  ! and then we have to take the imaginary part of the result
  
  prefac = cmplx(two, zero, kind(KREAL)) / (degen)
  if (dbg>0) write (out,*) 'prefac = ', prefac
  
  ctav(:) = cmplx(zero, zero, kind(KREAL)) ! isotropic B-term

  if (states_sos .eq.0) then
    ! assume the SOS limit of states was not set in the input
    states_sos = nstates
  else if (states_sos .lt. 0) then
    stop 'states_sos variable smaller than 0. aborting'
  else if (states_sos.gt.nstates) then
    stop 'states_sos variable > nstates. aborting'
  end if

  write (out,'(/1x,a,i5/)') 'No. of states used in SOS: ',states_sos

  ! create output files, open, and write namelist
  ! input for plot program

  do idir = 0,3
    
    write(mcdfile(idir),'(a,i1)') 'mcdspectrum-',idir
    open (iu_mcd(idir), file=trim(mcdfile(idir)), status='unknown', iostat=ios)
    if (ios /= 0) then
      write (err,*) 'problem opening file '//trim(mcdfile(idir))
      stop 'error with one of the mcd output files'
    end if
    
    ntemp = nlevels - skip - 1
    if (ntemp.lt.1) stop 'attempting to print data for less than 1 level'
        write(iu_mcd(idir),'(a,i7,a/1x,a,i7,a,f7.2,a,l,a///a)') &
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
          & 'The data will be written to file '//trim(mcdfile(idir))
    
     clist(:) = cmplx(zero, zero, kind(KREAL)) ! B-term
     
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
                  vec(jdir,3) = cmplx(zero, zero, kind(KREAL))
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
                     vec(jdir,3) = cmplx(zero, zero, kind(KREAL))
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
      write (iu_mcd(idir),'(1x,f14.2,3x,4(f20.8,2x))') &
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
    write (iu_mcd(0),'(1x,f14.2,3x,4(f20.8,2x))') &
      waveno(deltae), aimag(ctav(ilevel)), aimag(ctav(ilevel)*(debye**2)), &
      real(ctav(ilevel)), real(ctav(ilevel)*(debye**2))
  end do


  ! --------------------
  ! close mcd data files
  ! --------------------
  
  do idir = 0,3
    close (iu_mcd(idir))
    write (out,'(1x,a,1x,i7,1x,a,1x,a)') 'wrote B-term data for', &
      nstates-degen-skip,'transitions to file',trim(mcdfile(idir))  
  end do 
  
  ! --------------------------------------------------
  ! deallocate arrays, clean up if necessary, and exit
  ! --------------------------------------------------

  deallocate(energy, eldip, magdip, ctav, &
    deglist, clist, levels, elevel, accl)

  stop 'normal termination of mcd'
  
  ! ============================================================================
  
end program mcd_b_molcas
