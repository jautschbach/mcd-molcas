
subroutine print_constants

  use definitions

  use constants_parameters

  use shared_variables

  implicit none

  ! ============================================================================
  
  write (out, '(/1x,a/1x,a//, 9(1x,a,1x,f24.15/))') &
    'Some of the constants that this code may have used internally are:', &
    '(there may be some trailing numerical noise)', &
    'c = speed of light in au:      ', cspeed, &
    'c**2 :                         ', csq, &
    'Hartree to cm**(-1) conversion:', au2cm, &
    'Hartree to eV conversion:      ', au2ev, &
    'Boltzmann constant in au/K     ', boltzau, &
    'Boltzmann constant in cm^-1/K  ', boltzcm, &
    'Dipole au to Debye conversion  ', debye, &
    'free electron g factor:        ', ge, &
    'au to 1E-40 esu**2 cm**2:      ', d2au2cgs
    
end subroutine print_constants
    
