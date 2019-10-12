
module constants_parameters

  use definitions

  implicit none

  ! numerical constants:
  
  real(KREAL), parameter :: zero=0.0_KREAL
  real(KREAL), parameter :: one=1.0_KREAL
  real(KREAL), parameter :: two=2.0_KREAL
  real(KREAL), parameter :: three=3.0_KREAL
  real(KREAL), parameter :: half=0.5_KREAL
  real(KREAL), parameter :: third=one/three
  real(KREAL), parameter :: fourth=0.25_KREAL
  real(KREAL), parameter :: oneby15=(one/15.0_KREAL)
  real(KREAL), parameter :: small=1.0E-5_KREAL
  real(KREAL), parameter :: tiny=1E-10_KREAL

  complex(KREAL), parameter :: c0=cmplx( zero, zero, kind(KREAL))
  complex(KREAL), parameter :: cp1=cmplx( one,  zero, kind(KREAL))
  complex(KREAL), parameter :: cm1=cmplx(-one,  zero, kind(KREAL))
  complex(KREAL), parameter :: sqm1=cmplx( zero, one, kind(KREAL))

  ! physical constants, most from 2014 CODATA:
  
  real(KREAL), parameter :: ge=2.00231930436182_KREAL   ! g_e factor
  real(KREAL), parameter :: debye = 2.5417463E0_KREAL   ! au -> Debye
  real(KREAL), parameter :: au2cm=2.194746313702E+5_KREAL ! au -> cm**(-1)
  real(KREAL), parameter :: au2ev=27.21138602E0_KREAL   ! au -> eV
  real(KREAL), parameter :: boltzcm=0.69503457_KREAL    ! cm**(-1) / K
  real(KREAL), parameter :: thunit=5.02883E5_KREAL      ! C -> [theta]
  real(KREAL), parameter :: cspeed=137.035999139_KREAL  ! speed of light in a.u.
  real(KREAL), parameter :: csq = cspeed*cspeed         ! c**2

  ! conversion factors for dipole and rotstrength from
  ! atomic units to 1E-40 CGS, accounting for (1/c) not being included
  ! in the magnetic moments:

  real(KREAL), parameter :: d2au2cgs=64604.75024E0_KREAL
  real(KREAL), parameter :: rotconv=d2au2cgs/cspeed

  ! levi-civita symbol as complex variable:

  complex(KREAL), dimension(3,3,3), parameter :: &
    lc = RESHAPE((/ &
      c0, c0, c0, c0, c0,cm1, c0,cp1, c0,  &
      c0, c0,cp1, c0, c0, c0,cm1, c0, c0,  &
      c0,cm1, c0,cp1, c0, c0, c0, c0, c0   &
      /) ,shape(lc))

  ! index array for the quadrupole operator elements in
  ! packed storage:
  ! XX=1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5, ZZ=6

  integer(KINT), dimension(3,3), parameter :: &
    qindex = RESHAPE((/ 1, 2, 3, 2, 4, 5, 3, 5, 6/) ,shape(qindex))  

  ! other parameters, incl. file I/O units

  integer(KINT), parameter :: maxrecs = 1000000, nmax=10000
  
  integer(KINT), parameter :: out = 6, err = 0, iu_e = 7, iu_d = 8, iu_m = 9, &
    iu_op=5
  
  integer(KINT), dimension(0:3), parameter :: & 
    iu_out = (/10, 11, 12, 13/)

  ! name of the file from where we read the namelist with program options
  character(len=11), parameter :: optfile='options.dat'

  
end module constants_parameters
  
