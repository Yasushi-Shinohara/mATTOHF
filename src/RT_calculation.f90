!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine RT_calculation(ubk,ubkGS)
  use Constants
  use Global_variables, only : sys_name, NG, NB, NK, dvcell,&
       NT, dt, &
       Myrank, MPI_COMM_WORLD,ierr, MPI_WTIME
  use Functions
!  use Functions_MPI
  implicit none
  complex(8),intent(inout) :: ubk(NG,NB,NK)
  complex(8),intent(in) :: ubkGS(NG,NB,NK)
  character(256) :: RTinfo
  real(8) :: dns(NG)
  real(8) :: time(1:NT), Etot(1:NT), Nele(1:NT), Jtot(1:NT), E(1:NT), E2(1:NT), A(1:NT), Aave(1:NT)
  real(8) :: Etotm1, Nelem1, Am1, Jtotm1
  integer :: it, NTpulse
  real(8) :: tRTs, tRTe

  tRTs = MPI_WTIME()
  If (Myrank == 0) Write(*,'("########Real-time propagation is started.########")')
  If (Myrank == 0) then
    RTinfo = trim(sys_name)//'_EtotNJEAt.out'
    Open(10,file=trim(RTinfo),action='write')
    Write(10,'("#time(it), Etot(it)-Etot0, N(it)-N0, J(it), E(it), A(it)")')
  End If

  call Define_Afield(time,A,Aave,E,E2)

  Am1 = 0.d0
  dns(:) = ubk_dns(ubk)
!  dns(:,:,:) = ubk_dns_MPI(ubk)
  Nelem1 = sum(dns)*dvcell
  Etotm1 = ubk_Etot(ubk,Am1)
  Jtotm1 = ubk_Je(ubk,Am1)
!  Etot0 = ubk_Etot_MPI(ubk,A0)
!  Jtot0 = ubk_Je_MPI(ubk,A0)
  If (Myrank == 0) Write(*,'("#time(it), Etot(it)-Etot0, N(it)-N0, J(it), E(it), A(it)")')
  Do it = 1, NT
!    Call RT_TE_ubk(ubk,Aave(it))
!    Call RT_TEPC_ubk(ubk,Aave(it))
!    Call RT_TEPC2_ubk(ubk,Aave(it))
    Call RT_RK4_ubk(ubk,Aave(it),A(it))
    dns(:) = ubk_dns(ubk)
!    dns(:,:,:) = ubk_dns_MPI(ubk)
    Nele(it) = sum(dns)*dvcell
    Etot(it) = ubk_Etot(ubk,A(it))
    Jtot(it) = ubk_Je(ubk,A(it))
!    Etot(it) = ubk_Etot_MPI(ubk,A(1:3,it))
!    Jtot(:,it) = ubk_Je_MPI(ubk,A(1:3,it))
    If (Myrank == 0) then
      tRTe = MPI_WTIME()
      Write(10,'(20e20.12)') time(it), Etot(it), Etot(it)-Etotm1, Nele(it), Nele(it)-Nelem1, &
           Jtot(it), E(it), A(it)
      If (mod(it,100) == 0) Write(*,'(i6,20e14.6)') it, time(it), Etot(it), Etot(it)-Etotm1, &
           Nele(it), Nele(it)-Nelem1, Jtot(it), E(it), A(it)
      If (mod(it,1000) == 0) &
           Write(*,'("## calc-time for ",i4," iterations :",e12.4," [sec] =",f12.8," [hour]##")') &
           it,tRTe-tRTs, (tRTe-tRTs)/3600.d0
    End If
  End Do
  If (Myrank == 0) Write(*,'("########Real-time propagation is finished.########")')
  If (Myrank == 0) close(10)

  If (Myrank == 0) Call FourierTr(time,Jtot,Jtotm1,E,A)
  If (Myrank == 0) Call WindowFourierTr(time,Jtot,Jtotm1,E,A)
!  Ecomps = ubk_Ecomps(ubk,A0)

  return
End Subroutine RT_calculation
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Real-time evolution with Taylor-expansion (TE)
Subroutine RT_TE_ubk(ubk,Aave)
  use Constants
  use Global_variables, only : NG, NB, NK, dt
  use Functions
!  use Functions_MPI
  implicit none
  complex(8),intent(inout) :: ubk(NG,NB,NK)
  real(8),intent(in) :: Aave
  integer :: iter, iorder, Norder = 4
  complex(8) :: workbk(NG,NB,NK),ubk0(NG,NB,NK)

  ubk0 = ubk
  workbk = ubk
  Do iorder = 1,Norder
    workbk = -(zI*dt)*ubkwbk_hwbk(ubk0,workbk,Aave)/dble(iorder)
!    workbk = -(zI*dt)*ubkwbk_hwbk_MPI(ubk,workbk,Aave)/dble(iorder)
    ubk = ubk + workbk
  End Do
      
  return
End Subroutine RT_TE_ubk
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Real-time evolution with Taylor-expansion (TE) withint predictor-corrector (PC)
Subroutine RT_TEPC_ubk(ubk,Aave)
  use Constants
  use Global_variables, only : NG, NB, NK, dt
  use Functions
!  use Functions_MPI
  implicit none
  complex(8),intent(inout) :: ubk(NG,NB,NK)
  real(8),intent(in) :: Aave
  integer :: iter, iorder, Norder = 4
  complex(8) :: workbk(NG,NB,NK),ubk0(NG,NB,NK),ubk1(NG,NB,NK)

  ubk0 = ubk
  ubk1 = ubk
  workbk = ubk
  Do iorder = 1,Norder
    workbk = -(zI*dt)*ubkwbk_hwbk(ubk0,workbk,Aave)/dble(iorder)
    ubk1 = ubk1 + workbk
  End Do

  workbk = ubk
  Do iorder = 1,Norder
    workbk = -(zI*dt*0.5d0)*(ubkwbk_hwbk(ubk0,workbk,Aave)+ubkwbk_hwbk(ubk1,workbk,Aave))/dble(iorder)
    ubk = ubk + workbk
  End Do
  
  return
End Subroutine RT_TEPC_ubk
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Real-time evolution with Taylor-expansion (TE) withint predictor-corrector (PC)
Subroutine RT_TEPC2_ubk(ubk,Aave)
  use Constants
  use Global_variables, only : NG, NB, NK, dt
  use Functions
!  use Functions_MPI
  implicit none
  complex(8),intent(inout) :: ubk(NG,NB,NK)
  real(8),intent(in) :: Aave
  integer :: iter, iorder, Norder = 4
  complex(8) :: workbk(NG,NB,NK),ubk0(NG,NB,NK),ubk1(NG,NB,NK)

  ubk0 = ubk
  ubk1 = ubk
  workbk = ubk
  Do iorder = 1,Norder
    workbk = -(zI*dt)*ubkwbk_hwbk(ubk0,workbk,Aave)/dble(iorder)
    ubk1 = ubk1 + workbk
  End Do

  workbk = ubk
  ubk0 = 0.5d0*(ubk0+ubk1)
  Do iorder = 1,Norder
    workbk = -(zI*dt)*ubkwbk_hwbk(ubk0,workbk,Aave)/dble(iorder)
    ubk = ubk + workbk
  End Do
  
  return
End Subroutine RT_TEPC2_ubk
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Real-time evolution with Taylor-expansion (TE)
Subroutine RT_RK4_ubk(ubk,Aave,A)
  use Constants
  use Global_variables, only : NG, NB, NK, dt
  use Functions
!  use Functions_MPI
  implicit none
  complex(8),intent(inout) :: ubk(NG,NB,NK)
  real(8),intent(in) :: Aave, A
  real(8) :: Am
  complex(8) :: k1(NG,NB,NK),k2(NG,NB,NK),k3(NG,NB,NK),k4(NG,NB,NK)

  Am = 2.0d0*Aave - A
  k1 = -zI*ubkwbk_hwbk(ubk,ubk,Am)
  k2 = -zI*ubkwbk_hwbk(ubk+(0.5*dt)*k1,ubk+(0.5*dt)*k1,Aave)
  k3 = -zI*ubkwbk_hwbk(ubk+(0.5*dt)*k2,ubk+(0.5*dt)*k2,Aave)
  k4 = -zI*ubkwbk_hwbk(ubk+dt*k3,ubk+dt*k3,A)
  ubk = ubk + (k1+2.0d0*k2+2.0d0*k3+k4)*dt/6.d0
      
  return
End Subroutine RT_RK4_ubk
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Define_Afield(time,A,Aave,E,E2)
  use Constants
  use Global_variables, only : NT, dt, omegac, Tpulse, E0, dAc, field_option, phiCEP, nenvelope, &
       Myrank
  implicit none
  real(8), intent(out) :: time(1:NT), A(1:NT), Aave(1:NT), E(1:NT), E2(1:NT)
  integer :: it, ixyz, NTpulse
  real(8) ::  Am1, Aamp

!Time and field construction
!time(1) = 0.d0; time-origin but not "initial condition", initial condition is the GS
  Do it = 1,NT
    time(it) = dt*(it-1)
  End Do
  Am1 = 0.d0
  A(:) = 0.d0
  If (Myrank == 0) then
    Write(*,'(" Time and field infromatio+++++++++++++++++++")')
    Write(*,'(" Time-step: ",f12.4," [a.u.] =",f12.4," [fs]")') dt, dt*Atomtime
    Write(*,'(" The corresponding energy: ",f12.4," [a.u.] =",f12.4," [eV]")') &
         tpi/dt, tpi/dt*Hartree
    Write(*,'(" Total time: ",f12.4," [a.u.] =",f12.4," [fs]")') dble(NT)*dt, dble(NT)*dt*Atomtime
    Write(*,'(" The corresponding energy: ",f12.4," [a.u.] =",f12.4," [meV]")') &
         tpi/(dble(NT)*dt), tpi/(dble(NT)*dt)*Hartree*1000.d0
  End If

  If (field_option == 'LR') then
    If (Myrank == 0) Write(*,*) 'This calculation is the linear-response  with impulsive kick, whose amplitude', dAc, '[a.u.]'
    A(:) = dAc
!    write(*,*) 'The initial current is approximated as',J0-N0*Ac(:,1)/vcell
!    write(*,*) 'The initial energy is approximated as',Etot0+ 0.5*N0*(Ac(1,1)**2+Ac(2,1)**2+Ac(3,1)**3)

  Else if (field_option == 'LP') then
!Unit conversion 
    omegac = omegac/Hartree
    Tpulse = Tpulse/Atomtime
    E0 = E0/Atomfield
    If (Myrank == 0) Write(*,*) 'This calculation is under a laser pulse.'
    Aamp = E0/omegac
    If (Myrank == 0) Write(*,*) 'Amplitude of the A-field:',Aamp, '[a.u.]'
    NTpulse = int(min(Tpulse/dt,dble(NT)))
! A(t) = Aamp*sin(omegac*(t-T/2)+phi)*sin(pi*t/T)**nenvelope
    Do it = 1,NTpulse
      A(it) = Aamp*sin(omegac*(time(it)-0.5*Tpulse)+2.0*pi*phiCEP)*(sin(pi*time(it)/Tpulse))**nenvelope
    End do

  Else
    stop ('# Error stop: substituted field_option  is not acceptable.')
  End if

!Construct an averaging A-field: Aave(it) = (A(it)+A(it-1))/2
!Convert E-field from A-field: E = - d(Ac)/dt
  Aave(1) = 0.5d0*(A(1)+Am1)
  Do it = 2,NT
    Aave(it) = 0.5d0*(A(it)+A(it-1))
  End do
  E(1) = -(A(1)-Am1)/dt  !It may be replaced more accurate derivative formula
  Do it = 2,NT
    E(it) = -(A(it)-A(it-1))/dt  !It may be replaced more accurate derivative formula
  End do
  Do it = 1,NT
    E2(it) = E(it)**2
  End do
  If (Myrank == 0) then
    Write(*,*) 'Fluence of the field:', sum(E2)*dt*Atomfluence,' J/cm^2'
    Write(*,'("++++++++++++++++++++++++++++++++++++++++++++")')
  End If

  return
End Subroutine Define_Afield
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine FourierTr(time,Jtot,Jtotm1,E,A)
  use Constants
  use Global_variables, only : sys_name, NT,dt, omegaM,Nomega, dAc, field_option
  implicit none
  real(8), intent(in) :: time(1:NT), Jtot(1:NT), Jtotm1, E(1:NT), A(1:NT)
  character(256) :: FTinfo, LRinfo
  integer :: it,iomega
  real(8) ::  filter(1:NT), omega(1:Nomega), Jtemp(1:NT) 
  complex(8) :: JF(1:Nomega), EF(1:Nomega), AF(1:Nomega), eiwt(1:NT)
  complex(8) :: sigma(1:Nomega), eps(1:Nomega)

!Filter construction
  filter(:) = 1.0d0-3.0d0*(time(:)/time(NT))**2+2.0d0*(time(:)/time(NT))**3

!Frequency construction
  Do iomega = 1,Nomega
    omega(iomega) = dble(iomega)/dble(Nomega)*omegaM
  End Do

  Write(*,*) 'Jtotm1 =', Jtotm1
  Jtemp(:) = filter(:)*(Jtot(:) - Jtotm1)

  Do iomega = 1,Nomega
!     eiwt(:) = exp(+zI*omega(iomega)*(time(:)+0.5d0*dt))
     eiwt(:) = exp(+zI*omega(iomega)*time(:))
     JF(iomega) = sum(eiwt(:)*Jtemp(:))*dt
     EF(iomega) = sum(eiwt(:)*filter(:)*E(:))*dt
     AF(iomega) = sum(eiwt(:)*filter(:)*A(:))*dt
  End Do

  FTinfo = trim(sys_name)//'_JEAF.out'
  Open(10,file=trim(FTinfo),action='write')
  Write(10,'("#omega(iomega), JF(iomega), EF(iomega), AF(iomega)")')
  Do iomega = 1, Nomega
    Write(10,'(20e20.12)') omega(iomega), &
         real(JF(iomega)),aimag(JF(iomega)), &
         real(EF(iomega)),aimag(EF(iomega)), &
         real(AF(iomega)),aimag(AF(iomega))
  End Do
  Close(10)

  If (field_option == 'LR') then
    sigma(:) = -JF(:)/dAc
    eps(:) = fpi*zI*sigma(:)/omega(:)

    LRinfo = trim(sys_name)//'_sigmaeps.out'
    Open(10,file=trim(LRinfo),action='write')
    Write(10,'("#omega(iomega), sigma(iomega), eps(iomega)")')
    Do iomega = 1, Nomega
      Write(10,'(20e20.12)') omega(iomega), &
           real(sigma(iomega)),aimag(sigma(iomega)), &
           real(eps(iomega)),aimag(eps(iomega))
    End Do
  Close(10)
  End If

  return
End Subroutine FourierTr
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine WindowFourierTr(time,Jtot,Jtotm1,E,A)
  use Constants
  use Global_variables, only : sys_name, NT,dt, omegaM,Nomega, field_option, Tpulse, nenvelope
  implicit none
  real(8), intent(in) :: time(1:NT), Jtot(1:NT), Jtotm1, E(1:NT), A(1:NT)
  character(256) :: FTinfo, LRinfo
  integer :: it,iomega, NTpulse
  real(8) ::  filter(1:NT), omega(1:Nomega), Jtemp(1:NT) 
  complex(8) :: JF(1:Nomega), EF(1:Nomega), AF(1:Nomega), eiwt(1:NT)

  If (field_option == 'LR') return

!Filter construction
  NTpulse = int(min(Tpulse/dt,dble(NT)))
  filter(:) = 0.0d0
  Do it = 1,NTpulse
    filter(it) = (sin(pi*time(it)/Tpulse))**nenvelope
  End do

!Frequency construction
  Do iomega = 1,Nomega
    omega(iomega) = dble(iomega)/dble(Nomega)*omegaM
  End Do

  Jtemp(:) = filter(:)*(Jtot(:) - Jtotm1)

  Do iomega = 1,Nomega
!     eiwt(:) = exp(+zI*omega(iomega)*(time(:)+0.5d0*dt))
     eiwt(:) = exp(+zI*omega(iomega)*time(:))
     JF(iomega) = sum(eiwt(:)*Jtemp(:))*dt
     EF(iomega) = sum(eiwt(:)*filter(:)*E(:))*dt
     AF(iomega) = sum(eiwt(:)*filter(:)*A(:))*dt
  End Do

  FTinfo = trim(sys_name)//'_JEAF_window.out'
  Open(10,file=trim(FTinfo),action='write')
  Write(10,'("#omega(iomega), JF(iomega), EF(iomega), AF(iomega)")')
  Do iomega = 1, Nomega
    Write(10,'(20e20.12)') omega(iomega), &
         real(JF(iomega)),aimag(JF(iomega)), &
         real(EF(iomega)),aimag(EF(iomega)), &
         real(AF(iomega)),aimag(AF(iomega))
  End Do
  Close(10)

  return
End Subroutine WindowFourierTr
