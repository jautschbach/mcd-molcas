
module definitions

  implicit none

  ! define the byte length of integer and real variables:
  integer, parameter       :: KINT = kind(1)
  integer(KINT), parameter :: KREAL = kind(1.0d0)
  integer(KINT), parameter :: LCHARS = 160

end module definitions
  
