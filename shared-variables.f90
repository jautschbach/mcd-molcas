
module shared_variables

  ! this module is part of J. Autschbach's set of programs to process
  ! Molcas data for the generation of various types of spectral
  ! intensities
  
  ! (c) 2019-2022 Jochen Autschbach, SUNY Buffalo  
  
  use definitions

  implicit none

  ! allocatable arrays:
  
  real(KREAL), dimension(:), allocatable :: energy, elevel
  
  integer(KINT), dimension(:), allocatable :: deglist, levels, accl
  
  complex(KREAL), dimension(:,:,:), allocatable :: magdip, eldip, &
    elquad, angmom, velocity, veloquad

  ! other variables:
  
  integer(KINT) :: ios, dbg
  
  character(len=LCHARS) :: cs, cstemp, outfile(0:3)

  real(KREAL) :: deltae, rtemp, ctemp(2), kT
  integer(KINT) :: ntemp, skip1

  logical havespin, haveang, havedip, havevel, havequad, do_group
  data havespin, haveang, havedip, havequad, do_group / &
    .false., .false., .false., .false., .false. /

  save havespin, haveang, havedip, havequad, do_group

  save energy, elevel

  save deglist, levels, accl

  save magdip, eldip, elquad, velocity, veloquad

  save dbg

end module shared_variables
  
