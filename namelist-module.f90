
module namelist_module

  ! this module is part of J. Autschbach's set of programs to process
  ! Molcas data for the generation of various types of spectral
  ! intensities
  
  ! (c) 2019-2022 Jochen Autschbach, SUNY Buffalo

  use definitions

  implicit none
  
  real(KREAL) :: temp, ddelta
  
  integer(KINT) :: degen, nstates, nlevels, skip, &
    states_sos, debuglevel
  
  logical :: magdiag, nospin, noangmom, nodip, noquad, novel, &
    print_d, print_m, print_q, theta, usemag, polnotprop, &
    oldspiny
  
  namelist /options/ degen, temp, nstates, skip, magdiag, ddelta, &
    debuglevel, &
    nospin, noangmom, nodip, noquad, novel, print_d, print_m, print_q, &
    theta, states_sos, usemag, polnotprop, oldspiny

  ! populate namelist with default options at compile time

  data &
    degen, temp, nstates, states_sos, skip, magdiag, ddelta, &
    debuglevel, &
    nospin, noangmom, nodip, noquad, novel, &
    print_d, print_m, print_q, &
    theta, usemag, polnotprop, oldspiny / &
    0, 0.0_KREAL , 0, 0, 0, .false., 1E-5_KREAL, &
    0, & 
    .true., .true.,  .true., .true., .true., &
    .false., .false., .false., &
    .false., .false., .false. , .false. /

  save degen, temp, nstates, states_sos, skip, magdiag, ddelta, &
    debuglevel, & 
    nospin, noangmom, nodip, noquad, novel, print_d, print_m, print_q, &
    theta, usemag, polnotprop, oldspiny
  
end module namelist_module
  
