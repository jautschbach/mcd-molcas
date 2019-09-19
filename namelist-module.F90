
module namelist_module

  use definitions

  implicit none
  
  real(KREAL) :: temp, ddelta
  integer(KINT) :: degen, nstates, nlevels, skip, &
    states_sos
  logical :: magdiag, nospin, noangmom, noquad, print_d, print_m, print_q, &
    theta, usemag
  namelist /options/ degen, temp, nstates, skip, magdiag, ddelta, &
    nospin, noangmom, noquad, print_d, print_m, print_q, theta, states_sos, &
    usemag
  
end module namelist_module
  
