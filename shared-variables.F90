
module shared_variables

  use definitions

  implicit none

  ! allocatable arrays:
  
  real(KREAL), dimension(:), allocatable :: energy, elevel
  
  integer(KINT), dimension(:), allocatable :: deglist, levels, accl
  
  complex(KREAL), dimension(:,:,:), allocatable :: magdip, eldip, &
    elquad, angmom

  ! other variables:
  
  integer(KINT) :: ios, dbg
  
  character(len=LCHARS) :: cs, cstemp, outfile(0:3)

  real(KREAL) :: deltae, rtemp, ctemp(2), kT
  integer(KINT) :: ntemp, skip1

  logical havespin, haveang, havedip, havequad, do_group
  data havespin, haveang, havedip, havequad, do_group / &
    .false., .false., .false., .false., .false. /

  save havespin, haveang, havedip, havequad, do_group

  save energy, elevel

  save deglist, levels, accl

  save magdip, eldip, elquad

end module shared_variables
  
