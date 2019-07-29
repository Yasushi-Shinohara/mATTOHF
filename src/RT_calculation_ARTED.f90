!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine RT_calculation(ubk,ubkGS)
  use Constants
  use Global_variables, only : sys_name, NG1,NG2,NG3,NG, NB, NK, dvcell,&
       NT, dt, &
       Myrank, MPI_COMM_WORLD,ierr, MPI_WTIME
  use Functions
  use Functions_MPI
  implicit none
  complex(8),intent(inout) :: ubk(NG1,NG2,NG3,NB,NK)
  complex(8),intent(in) :: ubkGS(NG1,NG2,NG3,NB,NK)
  character(256) :: RTinfo
  real(8) :: dns(NG1,NG2,NG3)
!A should have negative time coordinate to take the average to a step before without if branching
  real(8) :: time(1:NT), Etot(1:NT), Nele(1:NT), Jtot(3,1:NT), E(3,1:NT), E2(1:NT), A(3,1:NT), Aave(3,1:NT)
  real(8) :: Etot0, Nele0, A0(3), Jtot0(3)
  integer :: it, NTpulse, ixyz
  real(8) :: tRTs, tRTe

  tRTs = MPI_WTIME()
  If (Myrank == 0) Write(*,'("########Real-time propagation is started.########")')
  If (Myrank == 0) then
    RTinfo = trim(sys_name)//'_EtotNJEAt.out'
    Open(10,file=trim(RTinfo),action='write')
    Write(10,'("#time(it), Etot(it)-Etot0, N(it)-N0, J(it), E(3,it), A(3,it)")')
  End If

  call Define_Afield(time,A,Aave,E,E2)

  A0 = 0.d0
!  dns(:,:,:) = ubk_dns(ubk)
  dns(:,:,:) = ubk_dns_MPI(ubk)
  Nele0 = sum(dns)*dvcell
!  Etot0 = ubk_Etot(ubk,A0)
!  Jtot0 = ubk_Je(ubk,A0)
  Etot0 = ubk_Etot_MPI(ubk,A0)
  Jtot0 = ubk_Je_MPI(ubk,A0)
  If (Myrank == 0) Write(*,'("#time(it), Etot(it)-Etot0, N(it)-N0, J(it), E(3,it), A(3,it)")')
  Do it = 1, NT
    Call RT_ubk(ubk,Aave(:,it))
!    dns(:,:,:) = ubk_dns(ubk)
    dns(:,:,:) = ubk_dns_MPI(ubk)
    Nele(it) = sum(dns)*dvcell
!    Etot(it) = ubk_Etot(ubk,A(1:3,it))
!    Jtot(:,it) = ubk_Je(ubk,A(1:3,it))
    Etot(it) = ubk_Etot_MPI(ubk,A(1:3,it))
    Jtot(:,it) = ubk_Je_MPI(ubk,A(1:3,it))
    If (Myrank == 0) then
      tRTe = MPI_WTIME()
      Write(10,'(20e20.12)') time(it), Etot(it), Etot(it)-Etot0, Nele(it), Nele(it)-Nele0, Jtot(:,it)
      If (mod(it,100) == 0) Write(*,'(i6,20e14.6)') it, time(it), Etot(it), Etot(it)-Etot0, Nele(it), Nele(it)-Nele0, Jtot(:,it)
      If (mod(it,1000) == 0) Write(*,'("## calc-time for ",i4," iterations :",e12.4," [sec] =",f12.8," [hour]##")') it,tRTe-tRTs, (tRTe-tRTs)/3600.d0
    End If
  End Do
  If (Myrank == 0) Write(*,'("########Real-time propagation is finished.########")')
  If (Myrank == 0) close(10)

  If (Myrank == 0) Call FourierTr(time,Jtot,Jtot0,E,A)
!  Ecomps = ubk_Ecomps(ubk,A0)

  return
End Subroutine RT_calculation
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Imaginary time evolution with given dns
Subroutine RT_ubk(ubk,Aave)
  use Constants
  use Global_variables, only : NG1,NG2,NG3, NB, NK, dt
  use Functions
  use Functions_MPI
  implicit none
  complex(8),intent(inout) :: ubk(NG1,NG2,NG3,NB,NK)
  real(8),intent(in) :: Aave(3)
  integer :: iter, iorder, Norder = 4
  complex(8) :: workbk(NG1,NG2,NG3,NB,NK)

  workbk = ubk
  Do iorder = 1,Norder
!    workbk = -(zI*dt)*ubkwbk_hwbk(ubk,workbk,Aave)/dble(iorder)
    workbk = -(zI*dt)*ubkwbk_hwbk_MPI(ubk,workbk,Aave)/dble(iorder)
    ubk = ubk + workbk
  End Do
      
  return
End Subroutine RT_ubk
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Define_Afield(time,A,Aave,E,E2)
  use Constants
  use Global_variables, only : NT, dt, omegac, Tpulse, E0, theta, phi, dAc, field_option, phiCEP, nenvelope, &
       Myrank
  implicit none
  real(8), intent(out) :: time(1:NT), A(3,1:NT), Aave(3,1:NT), E(3,1:NT), E2(1:NT)
  integer :: it, ixyz, NTpulse
  real(8) ::  evec(3), A0(3), Aamp

!Time and field construction
  Do it = 1,NT
    time(it) = dt*(it-1)
  End Do
  A0(:) = 0.d0
  A(:,:) = 0.d0
!Field direction ==This should be taken in more suitable place==
  evec(1) = cos(phi*tpi)*sin(theta*tpi)
  evec(2) = sin(phi*tpi)*sin(theta*tpi)
  evec(3) = cos(theta*tpi)
  If (Myrank == 0) then
    Write(*,'(" Time and field infromatio+++++++++++++++++++")')
    Write(*,'(" Time-step: ",f12.4," [a.u.] =",f12.4," [fs]")') dt, dt*Atomtime
    Write(*,'(" The corresponding energy: ",f12.4," [a.u.] =",f12.4," [eV]")') &
         tpi/dt, tpi/dt*Hartree
    Write(*,'(" Total time: ",f12.4," [a.u.] =",f12.4," [fs]")') dble(NT)*dt, dble(NT)*dt*Atomtime
    Write(*,'(" The corresponding energy: ",f12.4," [a.u.] =",f12.4," [meV]")') &
         tpi/(dble(NT)*dt), tpi/(dble(NT)*dt)*Hartree*1000.d0
    Write(*,'(" Field direction in spherical coordinate:theta =",f12.4," [deg], phi",f12.4," [deg]")') theta, phi
    Write(*,'(" Field direction in Cartesian coordinate",3f12.4)') evec(:)
  End If

  If (field_option == 'LR') then
    If (Myrank == 0) Write(*,*) 'This calculation is the linear-response  with impulsive kick, whose amplitude', dAc, '[a.u.]'
    Do ixyz = 1,3
      A(ixyz,:) = dAc*evec(ixyz)
    End do
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
    Do ixyz = 1,3
      Do it = 1,NTpulse
        A(ixyz,it) = Aamp*sin(omegac*(time(it)-0.5*Tpulse)+2.0*pi*phiCEP)*(sin(pi*time(it)/Tpulse))**nenvelope*evec(ixyz)
      End do
    End do

  Else
    stop ('# Error stop: substituted field_option  is not acceptable.')
  End if

!Construct an averaging A-field: Aave = 
!Convert E-field from A-field: E = - d(Ac)/dt
  Do ixyz = 1,3
    Aave(ixyz,1) = 0.5d0*(A(ixyz,1)+A0(ixyz))
    Do it = 2,NT
      Aave(ixyz,it) = 0.5d0*(A(ixyz,it)+A(ixyz,it-1))
    End do
  End do
  Do ixyz = 1,3
    E(ixyz,1) = -(A(ixyz,1)-A0(ixyz))/dt  !It may be replaced more accurate derivative formula
    Do it = 2,NT
      E(ixyz,it) = -(A(ixyz,it)-A(ixyz,it-1))/dt  !It may be replaced more accurate derivative formula
    End do
  End do
  Do it = 1,NT
    E2(it) = E(1,it)**2 + E(2,it)**2 + E(3,it)**2
  End do
  If (Myrank == 0) then
    Write(*,*) 'Fluence of the field:', sum(E2)*dt*Atomfluence,' J/cm^2'
    Write(*,'("++++++++++++++++++++++++++++++++++++++++++++")')
  End If

  return
End Subroutine Define_Afield
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine FourierTr(time,Jtot,Jtot0,E,A)
  use Constants
  use Global_variables, only : sys_name, NT,dt, omegaM,Nomega, theta, phi, dAc, field_option
  implicit none
  real(8), intent(in) :: time(1:NT), Jtot(3,1:NT), Jtot0(3), E(3,1:NT), A(3,1:NT)
  character(256) :: FTinfo, LRinfo
  integer :: it,ixyz,iomega
  real(8) ::  filter(1:NT), omega(1:Nomega), Jtemp(3,1:NT) , evec(3)
  complex(8) :: JF(1:3,1:Nomega), EF(1:3,1:Nomega), AF(1:3,1:Nomega), eiwt(1:NT)
  complex(8) :: sigma(1:3,1:Nomega), eps(1:3,1:Nomega)

!Filter construction
  filter(:) = 1.0d0-3.0d0*(time(:)/time(NT))**2+2.0d0*(time(:)/time(NT))**3

!Frequency construction
  Do iomega = 1,Nomega
    omega(iomega) = dble(iomega)/dble(Nomega)*omegaM
  End Do

  Write(*,*) 'Jtot0 =', Jtot0
  Do ixyz = 1,3
    Jtemp(ixyz,:) = filter(:)*(Jtot(ixyz,:) - Jtot0(ixyz))
  End Do

  Do iomega = 1,Nomega
!     eiwt(:) = exp(+zI*omega(iomega)*(time(:)+0.5d0*dt))
     eiwt(:) = exp(+zI*omega(iomega)*time(:))
     Do ixyz = 1,3
       JF(ixyz,iomega) = sum(eiwt(:)*Jtemp(ixyz,:))*dt
       EF(ixyz,iomega) = sum(eiwt(:)*filter(:)*E(ixyz,:))*dt
       AF(ixyz,iomega) = sum(eiwt(:)*filter(:)*A(ixyz,:))*dt
     End Do
  End Do

  FTinfo = trim(sys_name)//'_JEAF.out'
  Open(10,file=trim(FTinfo),action='write')
  Write(10,'("#omega(iomega), JF(3,iomega), EF(3,iomega), AF(3,iomega)")')
  Do iomega = 1, Nomega
    Write(10,'(20e20.12)') omega(iomega), &
         (real(JF(ixyz,iomega)),aimag(JF(ixyz,iomega)),ixyz=1,3), &
         (real(EF(ixyz,iomega)),aimag(EF(ixyz,iomega)),ixyz=1,3), &
         (real(AF(ixyz,iomega)),aimag(AF(ixyz,iomega)),ixyz=1,3)
  End Do
  Close(10)

  If (field_option == 'LR') then
    evec(1) = cos(phi*tpi)*sin(theta*tpi)
    evec(2) = sin(phi*tpi)*sin(theta*tpi)
    evec(3) = cos(theta*tpi)
    sigma(:,:) = -JF(:,:)/dAc
    Do ixyz = 1,3
      eps(ixyz,:) = evec(ixyz)+fpi*zI*sigma(ixyz,:)/omega(:)
    End Do

    LRinfo = trim(sys_name)//'_sigmaeps.out'
    Open(10,file=trim(LRinfo),action='write')
    Write(10,'("#omega(iomega), sigma(3,iomega), eps(3,iomega)")')
    Do iomega = 1, Nomega
      Write(10,'(20e20.12)') omega(iomega), &
           (real(sigma(ixyz,iomega)),aimag(sigma(ixyz,iomega)),ixyz=1,3), &
           (real(eps(ixyz,iomega)),aimag(eps(ixyz,iomega)),ixyz=1,3)
    End Do
  Close(10)
  End If

  return
End Subroutine FourierTr
