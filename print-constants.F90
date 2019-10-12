
subroutine print_constants

  use definitions

  use constants_parameters

  use shared_variables

  implicit none

  ! ============================================================================
  
  write (out, '(/1x,a//, 8(1x,a,1x,f22.13/))') &
    'Some of the constants that this code may have used internally are:', &
    'c = speed of light in au:      ', cspeed, &
    'c**2 :                         ', csq, &
    'Hartree to cm**(-1) conversion:', au2cm, &
    'Hartree to eV conversion:      ', au2ev, &
    'Dipole au to Debye conversion  ', debye, &
    'free electron g factor:        ', ge, &
    'au to 1E-40 esu**2 cm**2:      ', d2au2cgs
    
end subroutine print_constants
    
