
module namelist_module

  use definitions

  implicit none
  
  real(KREAL) :: temp, ddelta
  
  integer(KINT) :: degen, nstates, nlevels, skip, &
    states_sos
  
  logical :: magdiag, nospin, noangmom, nodip, noquad, &
    print_d, print_m, print_q, theta, usemag, polnotprop
  
  namelist /options/ degen, temp, nstates, skip, magdiag, ddelta, &
    nospin, noangmom, nodip, noquad, print_d, print_m, print_q, &
    theta, states_sos, usemag, polnotprop

  ! populate namelist with default options at compile time

  data degen, temp, nstates, skip, magdiag, ddelta, &
    nospin, noangmom, nodip, noquad, print_d, print_m, print_q, &
    theta, states_sos, usemag, polnotprop / &
    0, 0.0_KREAL , 0, 0, .false., 1E-5_KREAL, .true., .true., &
    .true., .true., .false., .false., .false., .false., 0, .false., .false. /

  save degen, temp, nstates, skip, magdiag, ddelta, &
    nospin, noangmom, nodip, noquad, print_d, print_m, print_q, &
    theta, states_sos, usemag, polnotprop
  
end module namelist_module
  
