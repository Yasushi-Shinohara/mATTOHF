!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Following values never changes after initialization
Module Global_variables
  implicit none
  character(256) :: input_data
  character(256) :: sys_name = ''
!Default values for options
  character(256) :: dir_name = './'
  character(32) :: DFT_option = 'hand'
  logical :: Band_option = .false.
  character(32) :: TD_option = 'RK4'
  character(256) :: functional = 'HF' !'HF', 'PZ', 'PBE'
  integer :: Nscf=100 !number of SCF 
  integer :: NCG=5 !number of CG sequnence within single SCF cycle
  integer :: NT = 2000 !number of time step
  real(8) :: dt = 0.5 ![a.u.]
  integer :: Nomega = 4000 !number of frequency for FT
  real(8) :: omegaM = 4.0 ![a.u.]
  character(2) :: field_option = 'LP' !'LR' or 'LP', Linear-Response or Laser-Pulse
  real(8) :: dAc = 0.0001     !Distortion for LR calculation
! Ac(t) = Acamp*sin(omegac*(t-T/2)+phiCEP)*sin(pi*t/T)**nenvelope
  real(8) :: omegac = 1.55d0 ![eV]
  real(8) :: Tpulse = 16.0d0 ![fs]
  real(8) :: phiCEP = 0.5d0 ![2 pi] !The name is modified from SBE.py
  integer :: nenvelope = 4
  real(8) :: E0 = 1.0d0 ![V/nm]

!Spatial information
  integer :: NG, NK
  real(8),allocatable :: rx(:)
  real(8),allocatable :: Gx(:)
  complex(8),allocatable :: eiGrx(:,:), eiGrxc(:,:)
  integer :: nGzero
  real(8),allocatable :: kx(:)
  real(8),allocatable :: vCGkk(:,:)
  real(8),allocatable :: vCG(:)
  integer :: NB, NBocc, NKK, NKKs, NKKe
  integer,allocatable :: NKKsarray(:), NKKearray(:)
  integer,allocatable :: ikktable(:,:),ikrtable(:),ikctable(:)
  integer,allocatable :: ibbocctable(:,:),iboccrtable(:),iboccctable(:)
!Arrays
  real(8) :: vcell,dvcell
  real(8) :: acell, ax, bx
  real(8),allocatable :: occbk(:,:)
!Pseudopotential
  integer :: Nspecies !Number of species in a cell
!  integer,allocatable :: Zatom(:) !Charges
  real(8),allocatable :: Zeffatom(:) !Charges
  real(8),allocatable :: alphaatom(:) !Softening parameters
  integer :: Nion !Number of total ion
  real(8),allocatable :: Rion(:)
  integer,allocatable :: Kion(:)
  complex(8),allocatable :: rhoionG(:)
  real(8),allocatable :: vion(:)

!FFTW staff
  integer(8) :: planf,planb
!MPI
  include 'mpif.h'
  integer :: Myrank, Nprocs, ierr, irequest, istatus(MPI_STATUS_SIZE)


End Module Global_Variables
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
