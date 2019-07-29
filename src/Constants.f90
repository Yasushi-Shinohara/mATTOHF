!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Module Constants
  implicit none
!Mathematicl constants
  real(8), parameter :: pi = 4.d0*atan(1.0d0)
  real(8), parameter :: tpi = 2.d0*pi
  real(8), parameter :: fpi = 4.d0*pi
  complex(8), parameter :: zI = (0.d0,1.d0)
!Physical constants
  real(8), parameter :: sol = 137.0d0 !speed of light
  real(8), parameter :: aB = 0.0529177d0 !nanometer
  real(8), parameter :: Hartree = 27.2116d0 !eV
  real(8), parameter :: Atomtime = 0.02419d0 !fs
  real(8), parameter :: Atomfield = 514.2d0 !V/nm
  real(8), parameter :: ch = 1240.0d0 !eV * nm
  real(8), parameter :: chbar = 197.3d0 ! eV * nm
  real(8), parameter :: halfepsc = 3.509d16 ! W/cm^2 \frac{1}{2}*\epsilon_0 * c
  real(8), parameter :: Atomfluence = halfepsc*Atomtime*1.0d-15 ! J/cm^2 ,W/cm^2 * fs = femto J/cm^2
!LDA Constants
  real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
  real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
  real(8),parameter :: CU=0.002d0,DU=-0.0116d0

End Module Constants
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
