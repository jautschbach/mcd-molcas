
subroutine read_options

  ! this routine is part of J. Autschbach's set of programs to process
  ! Molcas data for the generation of various types of spectral
  ! intensities
  
  ! (c) 2019-2022 Jochen Autschbach, SUNY Buffalo

  use definitions

  use constants_parameters

  use namelist_module

  use shared_variables

  implicit none


  ! ===========================================================================

  open (iu_op, file=optfile, status='old', iostat=ios)
  if (ios /= 0) then
    stop 'error: options file does not exist'
  end if
  
  
  read (iu_op,options, iostat=ios)
  if (ios /= 0) then
    close (iu_op)
    stop 'error reading namelist input from options file'
  end if

  close (iu_op)

  dbg = debuglevel

  if (dbg>0) write (out,'(/1x,a,i2)') 'Debug level = ',dbg 
  
  return
  

end subroutine read_options
