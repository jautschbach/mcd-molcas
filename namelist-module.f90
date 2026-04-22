
module namelist_module

  ! this module is part of J. Autschbach's set of programs to process
  ! Molcas data for the generation of various types of spectral
  ! intensities and other properties

  ! (c) 2019-2026 Jochen Autschbach (JA), SUNY Buffalo

  use definitions

  implicit none

  real(KREAL) :: temp, ddelta, trange(2)

  integer(KINT) :: degen, nstates, nlevels, skip, &
    states_sos, debuglevel

  logical :: magdiag, nospin, noangmom, nodip, noquad, novel, &
    print_d, print_m, print_q, theta, polnotprop, &
    oldspiny, specgen

  namelist /options/ &
    degen, &      ! [int] degeneracy of the ground state. determined automaticcaly if not set
    temp, &       ! [real] absolute temperature in kelvin. used for C-term spectra
    nstates, &    ! [int] number of states used for the input data
    skip, &       ! [int] number of excited states to skip before collecting data
    magdiag, &    ! [log] explicit diagonalization of Zeeman operator (x,y,z) for degenerate states
    ddelta, &     ! [real] criterion (Hartree) to determine state degeneracy
    debuglevel, & ! [int] 0: no debug outpout. 1: some output. 2 or higher: LOTS of debug output
    nospin, &     ! [log] if .T., do not process matrix elemenst of the spin operator
    noangmom, &   ! [log] if .T., do not process matrix elements of the orbital angular momentum
    nodip, &      ! [log] if .T., do not process electric dipole matrix elements
    noquad, &     ! [log] if .T., do not process electric quadrupole matrix elements
    novel, &      ! [log] if .T., do not process velocity-gauge matrix elements
    print_d, &    ! [log] print electric dipole matrix elements (only some codes)
    print_m, &    ! [log] print magnetic dipole matrix elements (only some codes)
    print_q, &    ! [log] print electric quadrupole matrix elements (only some codes)
    theta, &      ! [log] output MCD as ellipticity (option passed on to plotspec inputs)
    states_sos, & ! [int] number of states to include in the SOS for B-terms
    polnotprop, & ! [log] output absorption spectra for x,y,z polarization instead of propagation
    oldspiny, &   ! [log] treat Spin-y matrices as in Molcas from 2029 or prior (see read-data-files.f90)
    trange, &     ! [real(2)] temperature range (K) for susceptibility code
    specgen       ! [log] create output data files for use with spec-gen (from JA GitHub)

  ! populate namelist with default options at compile time

  data &
    degen, temp, nstates, states_sos, skip, magdiag, ddelta, &
    debuglevel, &
    nospin, noangmom, nodip, noquad, novel, &
    print_d, print_m, print_q, &
    theta, polnotprop, oldspiny, specgen / &
    0, 0.0_KREAL , 0, 0, 0, .false., 1E-5_KREAL, &
    0, &
    .true., .true.,  .true., .true., .true., &
    .false., .false., .false., &
    .false., .false., .false., .false. /

  data trange /1.0_KREAL,301.0_KREAL/

  save degen, temp, nstates, states_sos, skip, magdiag, ddelta, &
    debuglevel, &
    nospin, noangmom, nodip, noquad, novel, print_d, print_m, print_q, &
    theta, polnotprop, oldspiny, trange, specgen

end module namelist_module
