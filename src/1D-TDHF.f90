!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Program Main
!  use Global_variables, only : input_data,NT
  use Global_variables 
  use Functions
!$ use Omp_lib
  implicit none
  integer :: iter, ixyz
  integer :: ikk
  real(8) :: ts,te, tFocks,tFocke, tpreps,tprepe, tGSs, tGSe, tRTs, tRTe
  integer :: irank, numk=0, itag
!Following arrays change dynamically in a sequence, GS converenge and time-propagation
!Orbital function, ubk(1:NG,1:NG2,1:NG3,NB,NK)
  complex(8),allocatable :: ubk(:,:,:),tempbk(:,:,:),ubkGS(:,:,:)
  real(8),allocatable :: dns(:)
!Blocked matrix for Fock-term operation, VFbbkk(1:NG,1:NG2,1:NG3,NB,NB,NK,NK)
  complex(8),allocatable :: VFbbkk(:,:,:,:,:,:,:)
  complex(8),allocatable :: VFbb(:,:,:)
!
  real(8) :: Etot, A0

  call MPI_init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Myrank,ierr)
  Allocate(NKKsarray(0:Nprocs-1), NKKearray(0:Nprocs-1))

  If (Myrank==0) then 
    Call Print_date('s')
    Write(*,'("System information +++++++++")')
!$  If(Nprocs*0 == 0) then
!$    Write(*,'("parallel = Hybrid")')
!$  Else
      Write(*,'("parallel = Flat MPI")')
!$  End if
    Write(*,'("Number of MPI process =",i8)') Nprocs
!$  Write(*,'("Number of OpenMP thread =",i8)')  omp_get_max_threads()
    Write(*,'("++++++++++++++++++++++++++++")')
  End If
  ts = MPI_WTIME()

  tpreps = MPI_WTIME()
!A MPI policy: read and preconditinings for epsbk and pmat are only performed at Myrank=0, then they are broadcasted into whole ranks.
  If (Myrank==0) then
    Write(*,*) 'You can specify a file for input data just below(If you do not have, just input [""]'
    Read(*,*) input_data
    Write(*,*) trim(input_data)
    If (input_data == '') then
      input_data = 'Default value'
      write(*,*) 'Defalut value will be used.'
      write(*,*) '++++++++'
    Else
      write(*,*) 'input_data = ',trim(input_data)
      write(*,*) '++++++++'
    End if

    If (input_data/='Default value') Call Read_input
  End If

  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  Call Initialization
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  Allocate(ubk(1:NG,1:NB,1:NK))
  Allocate(tempbk(1:NG,1:NB,1:NK))
  Allocate(dns(1:NG))
  Allocate(VFbb(1:NG,1:NBocc,1:NBocc))

  If (Myrank == 0) then
    Write(*,'("Memory consumption+++++++++++++++++++++++++++++++++++++++")') 
    Write(*,'("ubk(1:NG,1:NB,1:NK): ",es12.3," GByte in each process")') 16.d0*dble(NG)*dble(NB)*dble(NK)/1.0e9
    Write(*,'("(ubk*NCG (only in GS calc.): ",es12.3," GByte in each process)")') 16.d0*dble(NG)*dble(NB)*dble(NK)*dble(NCG)/1.0e9
!$  Write(*,'("(ubk in OpenMP copy: ",es12.3," GByte in each process)")') &
      16.d0*dble(NG)*dble(NB)*dble(NK)*dble(omp_get_max_threads())/1.0e9
!$  Write(*,'("(ubk*NCG in OpenMP copy: ",es12.3," GByte in each process)")') &
      16.d0*dble(NG)*dble(NB)*dble(NK)*dble(omp_get_max_threads())*dble(NCG)/1.0e9 
    If (trim(functional)=='HF') then
      Write(*,'("vCGkk(1:NG,NKKs:NKKe): ",es12.3," GByte in each process (maximum value)")') &
           16.d0*dble(NG)*dble(maxval(NKKearray(:)-NKKsarray(:)+1))/1.0e9
      Write(*,'("vCGkk : ",es12.3," GByte in total")') &
           16.d0*dble(NG)*dble(NKK)/1.0e9
      Write(*,'("ikktable(NK,NK): ",es12.3," GByte in each process")') 4.d0*dble(NK*NK)/1.0e9
      Write(*,'("ikrtable(NKK): ",es12.3" GByte in each process")') 4.d0*dble(NKK)/1.0e9
      Write(*,'("ikctable(NKK): ",es12.3," GByte in each process")') 4.d0*dble(NKK)/1.0e9
      Write(*,*)
      Write(*,'("All above: ",es12.3," GByte in each process")') &
           (16.d0*dble(NG)*dble(NB)*dble(NK) + 16.d0*dble(NG)*dble(maxval(NKKearray(:)-NKKsarray(:)+1)) + 3.d0*4.d0*dble(NKK))/1.0e9
    End If
    Write(*,'("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")') 
  End If
  Call Init_wf(ubk,tempbk)
  Call Gram_Schmidt(ubk)
  tprepe = MPI_WTIME()

  tGSs = MPI_WTIME()
  dns = ubk_dns(ubk)
  If (Myrank==0) Write(*,*) 'Number of particle per cell',sum(dns)*vcell/dble(NG)
  Call Prep_periodic_ion_potential !Constructing vion
  A0 = 0.d0
  Etot = ubk_Etot(ubk,A0)
  if(Myrank == 0) write(*,*) 'Total energy at initial =', Etot
  if(Myrank == 0) write(*,*) 'This is the end of preparation for ground state calculation'

  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
  Call GS_calculation(ubk)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

!  If(Band_option) Call Band_map(ubk)
  tGSe = MPI_WTIME()

  tRTs = MPI_WTIME()
  Allocate(ubkGS(1:NG,1:NB,1:NK))
  ubkGS = ubk
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
  Call RT_calculation(ubk,ubkGS)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
  tRTe = MPI_WTIME()

  te = MPI_WTIME()
  If (Myrank == 0) then
    Write(*,'("The execution is finished successfully!")') 
    Write(*,'("Time for preparation =",f12.8," [sec]")') tprepe-tpreps
    Write(*,'("Time for GS calculatoin =",e12.4," [sec] =",f12.8," [min]")') tGSe-tGSs, (tGSe-tGSs)/60.d0
    Write(*,'("Time for RT calculatoin =",e12.4," [sec] =",f12.8," [hour]")') tRTe-tRTs, (tRTe-tRTs)/3600.d0
    Write(*,'("Time for an interaction in RT =",e12.4," [sec]")') (tRTe-tRTs)/dble(NT)
    Write(*,'("Time for total  =",e12.4," [sec] =", f12.8," [hour]")') te-ts, (te-ts)/3600.d0
    Call Print_date('e')
  End If

  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
  call MPI_FINALIZE(ierr) 
  stop
End Program Main
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Gram_Schmidt(ubk)
  use Constants
  use Global_variables, only : NG,NG,vcell,NB,NK, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr
  implicit none
  complex(8), intent(inout) :: ubk(1:NG,1:NB,1:NK)
  integer :: ib,ibt,ik
  real(8) :: s
  complex(8) :: zov

!!!$omp parallel do private(ib,ibt,zov,s)
  Do ik = 1, NK
    Do ib = 1, NB
      Do ibt=1,ib-1
        zov=sum(conjg(ubk(:,ibt,ik))*ubk(:,ib,ik))*vcell/dble(NG)
        ubk(:,ib,ik)=ubk(:,ib,ik)-ubk(:,ibt,ik)*zov
      End Do
      s=sum(abs(ubk(:,ib,ik))**2)*vcell/dble(NG)
      ubk(:,ib,ik)=ubk(:,ib,ik)/sqrt(s)
    End Do
  End Do
  Call MPI_BCAST(ubk,NG*NB*NK,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

  return
End Subroutine Gram_Schmidt
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Init_wf(ubk,tempbk)
  use Constants
  use Functions
  use Global_variables
  implicit none
  complex(8), intent(out) :: ubk(1:NG,1:NB,1:NK),tempbk(1:NG,1:NB,1:NK)
  integer :: iseed
  real(8) :: rnd
  real(8) :: rxc,ryc,rzc,r2(1:NG),rnd1,rnd2,rnd3
  integer :: ib,ik

!Initial guess for the orbital
  iseed=123
  Do ik=1,NK
    Do ib=1,NB
      Call quickrnd(iseed,rnd)
      rnd1 = rnd
      rxc = rnd1*acell
      r2(:) = (rx(:)-rxc)**2
      ubk(:,ib,ik) = exp(-0.5d0*r2(:))
    End Do
  End Do
  tempbk(:,:,:) = ubk(:,:,:)

  return
End Subroutine Init_wf
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine quickrnd(iseed,rnd)
  implicit none
  integer,parameter :: im=6075,ia=106,ic=1283
  integer, intent(inout) :: iseed
  real(8), intent(inout) :: rnd

  iseed=mod(iseed*ia+ic,im)
  rnd=dfloat(iseed)/dfloat(im)

  return
End Subroutine quickrnd
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Read_input
  use Constants
  use Global_variables
  implicit none
  character(256) :: str
  integer :: ispe,ii

  Open(10,file=input_data,status='old',action='read') 
  Do while(.true.)
    read(10,*,end=1000) str
    Select case(trim(str))
      case('sys_name')
        read(10,*) sys_name
      case('dir_name')
        read(10,*) dir_name
      case('DFT_option')
        read(10,*) DFT_option
      case('Band_option')
        read(10,*) Band_option
      case('TD_option')
        read(10,*) TD_option
      case('acell')
        read(10,*) acell
      case('NB')
        read(10,*) NB
      case('NBocc')
        read(10,*) NBocc
      case('NG')
        read(10,*) NG
      case('NK')
        read(10,*) NK
      case('atoms')
        read(10,*) Nspecies
!        Allocate(Zatom(Nspecies),alphaatom(Nspecies))
!        read(10,*) Zatom(1:Nspecies)
        Allocate(Zeffatom(Nspecies),alphaatom(Nspecies))
        read(10,*) Zeffatom(1:Nspecies)
        read(10,*) alphaatom(1:Nspecies)
        read(10,*) Nion
        Allocate(Rion(Nion),Kion(Nion))
        Do ii = 1, Nion
          read(10,*) Rion(ii), Kion(ii)
        End Do
      case('Nscf')
        read(10,*) Nscf
      case('NCG')
        read(10,*) NCG
      case('NT')
        read(10,*) NT
      case('dt')
        read(10,*) dt
      case('Nomega')
        read(10,*) Nomega
      case('omegaM')
        read(10,*) omegaM
      case('field_option')
        read(10,*) field_option
      case('dAc')
        read(10,*) dAc
      case('omegac')
        read(10,*) omegac
      case('Tpulse')
        read(10,*) Tpulse
      case('phiCEP')
        read(10,*) phiCEP
      case('nenvelope')
        read(10,*) nenvelope
      case('E0')
        read(10,*) E0
    End select
  End do
1000 continue
  close(10)
  return
End Subroutine Read_input
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Print_date(se)
  implicit none
  character(1), intent(in) :: se
  character(8) :: date
  character(10) :: time
  character(5) :: zone
  integer :: ia(8)

  Call Date_and_time(date, time, zone, ia)
  Write(*,*) '======================================================'
  If (se=='s')  Write(*,*) 'This is calculation is startred at YYYY/MM/DD HH:MM:SS'
  If (se=='e')  Write(*,*) 'This is calculation is ended at YYYY/MM/DD HH:MM:SS'
  Write(*,'(2x,i4,"/",i2,"/",i2,"    ",i2,":",i2,":"i2)') ia(1),ia(2),ia(3),ia(5),ia(6),ia(7)
  Write(*,*) 'Time zone is ',zone,' from UTC'
  Write(*,*) '======================================================'

  return
End Subroutine Print_date
