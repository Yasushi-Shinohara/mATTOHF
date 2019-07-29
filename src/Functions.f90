!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Module Functions
  implicit none
contains
!==
  Function cross(a,b) result(c)
    implicit none
    real(8) :: c(1:3)
    real(8), intent(in) :: a(1:3), b(1:3)
    
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
    return
  End Function cross
!==
!==
  Function bandtrace(array) result(trace)
    use Global_Variables, only : NB
    implicit none
    complex(8), intent(in) :: array(1:NB,1:NB)
    complex(8) :: trace 
    integer :: ib
    
    trace = (0.d0,0.d0)
    Do ib = 1,NB
      trace = trace + array(ib,ib)
    End do
    return
  End Function bandtrace
!==
!==
  Function k_vCk(k) result(vCk)
    use Constants
    use Global_Variables, only : NK,acell,NG
    implicit none
    Real(8), intent(in) :: k
    Integer, parameter :: Nx = 100000
    Real(8) :: vCk, kx, dx
    integer ix

    dx = acell/dble(NG)/20.0d0
    vCk = 0.d0
    ix = 0
    Do While ((ix*dx) < (NK*acell))
      vCk = vCk + 0.5d0*(cos(k*(ix*dx))/sqrt((ix*dx)**2+1.0d0**2) + cos(k*((ix+1)*dx))/sqrt(((ix+1)*dx)**2+1.0d0**2))
      ix = ix + 1
    End Do
    vCk = 2.0d0*vCk*dx

    return
  End Function k_vCk
!==
!==
  Function ubk_dns(ubk) result(dns) !########Several 3D array is used as 1D array here for OMP.########
    use Global_Variables, only : NG,NB,NK,occbk
    implicit none
    complex(8), intent(in) :: ubk(1:NG,1:NB,1:NK)
    real(8) :: dns1D(1:NG), dns(1:NG)
    integer :: ib,ik,ig
    
    dns = 0.d0
    Do ig= 1,NG
      Do ik= 1,NK
        Do ib =1,NB
          dns(ig) = dns(ig) + occbk(ib,ik)*conjg(ubk(ig,ib,ik))*ubk(ig,ib,ik)
        End do
      End do
    End do
    return
  End Function ubk_dns
!==
!==
  Function ubk_Etot(ubk,A) result(Etot)
    use Global_Variables, only : NG,NB,NK,functional
    implicit none
    complex(8), intent(in) :: ubk(NG,NB,NK)
    real(8), intent(in) :: A
    real(8) :: Etot, Ekin, Eion, EH, Exc, dns(NG)
    
    dns = ubk_dns(ubk)
    Ekin = ubk_Ekin(ubk,A)
    Eion = dns_Eion(dns)
    EH = dns_EH(dns)
    Select case (trim(functional))
      Case('HF')
        Exc = ubk_Ex(ubk)
    End Select
    Etot = Ekin + Eion + EH + Exc
    return
  End Function ubk_Etot
!==
!==
  Function ubk_Ekin(ubk,A) result(Ekin)
    use Global_Variables, only : NG,NB,NK,dvcell,Gx,kx,occbk, planf, planb
    implicit none
    complex(8), intent(in) :: ubk(NG,NB,NK)
    real(8), intent(in) :: A
    real(8) :: Ekin
    complex(8) :: work(NG),work2(NG)
    real(8) :: kAx,kA2
    integer :: ib,ik
!FFTW staff
    include 'fftw3.f'

    Ekin = 0.d0
    Do ik= 1,NK
      kAx = kx(ik)+A
      kA2 = kAx**2
      Do ib =1,NB
        Call dfftw_execute_dft(planf,ubk(:,ib,ik),work)
        work(:) = 0.5d0*(Gx(:)**2 + 2.0d0*Gx(:)*kAx + kA2)*work(:)/dble(NG)
        Call dfftw_execute_dft(planb,work,work2(:))
        Ekin = Ekin + occbk(ib,ik)*sum(conjg(ubk(:,ib,ik))*work2(:))*dvcell
      End do
    End do
    
    return
  End Function ubk_Ekin
!==
!==
  Function dns_Eion(dns) result(Eion)
    use Global_Variables, only : NG,vion,dvcell
    implicit none
    real(8), intent(in) :: dns(NG)
    real(8) :: Eion

    Eion = sum(vion(:)*dns(:))*dvcell
    
    return
  End Function dns_Eion
!==
!==
  Function dns_EH(dns) result(EH)
    use Global_Variables, only : NG,dvcell
    implicit none
    real(8), intent(in) :: dns(NG)
    real(8) :: vH(NG),EH

    vH = dns_vH(dns)
    EH = 0.5d0*sum(vH(:)*dns(:))*dvcell

    return
  End Function dns_EH
!==
!==
  Function dns_vH(dns) result(vH)
    use Global_Variables, only : NG,vCG, planf,planb
    implicit none
    real(8), intent(in) :: dns(1:NG)
    real(8) :: vH(1:NG)
    complex(8) :: work1(1:NG),work2(1:NG)
!FFTW staff
    include 'fftw3.f'

    work1 = dns
    Call dfftw_execute_dft(planf,work1,work2)
    work2 = vCG*work2
    Call dfftw_execute_dft(planb,work2,work1)
    vH = dble(work1)/dble(NG)
    return
  End Function dns_vH
!==
!==
  Function ubk_Ex(ubk) result(Ex)
    use Global_Variables, only : NG,NB,NK, NBocc, occbk, dvcell
    implicit none
    complex(8), intent(in) :: ubk(NG,NB,NK)
    complex(8) :: vFubk(NG,NB,NK)
    real(8) :: Ex
    integer :: ib,ik

    vFubk = ubkwbk_vFwbk(ubk,ubk)
    Ex = 0.d0
    Do ik = 1,NK
      Do ib = 1,NBocc
        Ex = Ex + occbk(ib,ik)*sum(conjg(ubk(:,ib,ik))*vFubk(:,ib,ik))*dvcell
      End Do
    End Do
    Ex = 0.5d0*Ex
    
    return
  End Function ubk_Ex
!==
!==
!Hamiltonian operation to wbk with given ubk
  Function ubkwbk_hwbk(ubk,wbk,A) result(hwbk)
    use Global_Variables, only : NG,NB,NK,functional,vion
    implicit none
    complex(8), intent(in) :: ubk(NG,NB,NK), wbk(NG,NB,NK)
    real(8), intent(in) :: A
    complex(8) :: hwbk(NG,NB,NK)
    real(8) :: vH(NG),vxc(NG),dns(NG)
    integer :: ib,ik

    hwbk = ubk_tubk(wbk,A)
    dns = ubk_dns(ubk)
    vH = dns_vH(dns)
    Do ik = 1,NK
      Do ib = 1,NB
        hwbk(:,ib,ik) = hwbk(:,ib,ik) + (vH(:)+vion(:))*wbk(:,ib,ik)
      End Do
    End Do
    Select case (trim(functional))
      Case('HF')
        hwbk = hwbk + ubkwbk_vFwbk(ubk,wbk)
    End Select
    return
  End Function ubkwbk_hwbk
!==
!==
  Function ubk_tubk(ubk,A) result(tubk)
    use Global_Variables, only : NG, NB, NK, Gx, kx, planf,planb
    implicit none
    complex(8), intent(in) :: ubk(1:NG,1:NB,1:NK)
    real(8), intent(in) :: A
    complex(8) :: tubk(1:NG,1:NB,1:NK)
    complex(8) :: work(1:NG)
    real(8) :: kAx,kA2
    integer :: ib,ik
!FFTW staff
    include 'fftw3.f'

    Do ik= 1,NK
      kAx = kx(ik)+A
      kA2 = kAx**2
      Do ib =1,NB
        Call dfftw_execute_dft(planf,ubk(:,ib,ik),work)
        work(:) = 0.5d0*(Gx(:)**2 + 2.d0*Gx(:)*kAx + kA2)*work(:)/dble(NG)
        Call dfftw_execute_dft(planb,work,tubk(:,ib,ik))
      End do
    End do
    return
  End Function ubk_tubk
!==
!==
  Function ubkwbk_vFwbk(ubk,wbk) result(vFwbk)
    use Global_Variables, only : NG, NB,NBocc, NK,occbk, NKKs,NKKe,ikrtable,ikctable,vCGkk, planf,planb &
         , MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD,ierr, vcell
    implicit none
    complex(8), intent(in) :: ubk(1:NG,1:NB,1:NK), wbk(1:NG,1:NB,1:NK)
    complex(8) :: vFwbk(1:NG,1:NB,1:NK),vFwbk_l(1:NG,1:NB,1:NK)
    complex(8) :: pairrho (NG),pairrhoG(NG), VF(1:NG)
    integer :: ibr,ibc,ikk

    vFwbk_l(:,:,:) = 0.d0
    Do ikk = NKKs,NKKe
      Do ibr = 1,NB
        Do ibc = 1,NBocc
          pairrho(:) = wbk(:,ibr,ikrtable(ikk))*conjg(ubk(:,ibc,ikctable(ikk)))/dble(NK)
          Call dfftw_execute_dft(planf, pairrho, pairrhoG)
          pairrhoG(:) = -vCGkk(:,ikk)*pairrhoG(:)/dble(NG)
          Call dfftw_execute_dft(planb, pairrhoG, VF)
          vFwbk_l(:,ibr,ikrtable(ikk)) = vFwbk_l(:,ibr,ikrtable(ikk)) + VF(:)*ubk(:,ibc,ikctable(ikk))
        End Do
      End Do
    End Do
    Call MPI_ALLREDUCE(vFwbk_l,vFwbk,NG*NB*NK,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

    return
  End Function ubkwbk_vFwbk
!==
!==
!Energy computation w.r.t. each components
  Function ubk_Ecomps(ubk,A) result(Ecomps)
    use Global_Variables, only : NG,NB,NK,functional
    implicit none
    complex(8), intent(in) :: ubk(NG,NB,NK)
    real(8), intent(in) :: A
    real(8) :: Ecomps(4), Ekin, Epsl, Epsnl, EH, Exc, Eion, Etot, dns(NG)
    
    dns = ubk_dns(ubk)
    Ekin = ubk_Ekin(ubk,A)
    Eion = dns_Eion(dns)
    EH = dns_EH(dns)
    Select case (trim(functional))
      Case('HF')
        Exc = ubk_Ex(ubk)
    End Select
    Etot = Ekin + Eion + EH + Exc
    Ecomps(1) = Ekin
    Ecomps(2) = Eion
    Ecomps(3) = EH
    Ecomps(4) = Exc
    return
  End Function ubk_Ecomps
!==
!==
  Function ubk_Je(ubk,A) result(Je)
    use Constants
    use Global_Variables, only : NG, NB,NK, dvcell,vcell, Gx,kx, occbk, planf, planb
    implicit none
    complex(8), intent(in) :: ubk(NG,NB,NK)
    real(8), intent(in) :: A
    complex(8) :: JpkA
    complex(8) :: JpkAx
    real(8) :: Je
    complex(8) :: work1(NG),work2(NG)
    real(8) :: GkAx(NG)
    integer :: ib,ik
!FFTW staff
    include 'fftw3.f'

    JpkAx = 0.d0
    Do ik = 1,NK
      GkAx(:) = Gx(:) + kx(ik) + A
      Do ib = 1,NB
        Call dfftw_execute_dft(planf,ubk(:,ib,ik),work1)
        work1 = work1/dble(NG)
        Call dfftw_execute_dft(planb,GkAx(:)*work1(:),work2)
        JpkAx = JpkAx + occbk(ib,ik)*sum(conjg(ubk(:,ib,ik))*work2(:))
      End do
    End do
    JpkA = JpkAx*dvcell/vcell

    Je = -real(JpkA)
    return
  End Function ubk_Je
!==
!==
  Function ubk_epsbk(ubk) result(epsbk)
    use Global_Variables, only : NG,NB,NK,dvcell
    implicit none
    complex(8), intent(in) :: ubk(NG,NB,NK)
    complex(8) :: hubk(NG,NB,NK)
    real(8) :: epsbk(NB,NK), A0=0.d0
    integer :: ib,ik

    hubk = ubkwbk_hwbk(ubk,ubk,A0)
    Do ik = 1,NK
      Do ib = 1,NB
        epsbk(ib,ik) = real(sum(conjg(ubk(:,ib,ik))*hubk(:,ib,ik)))*dvcell
      End Do
    End Do
    
    return
  End Function ubk_epsbk
!==
!==
  Function ubk_epserrorbk(ubk) result(epserrorbk)
    use Global_Variables, only : NG,NB,NK,dvcell
    implicit none
    complex(8), intent(in) :: ubk(NG,NB,NK)
    complex(8) :: hubk(NG,NB,NK)
    real(8) :: epsbk(NB,NK),epserrorbk(NB,NK), A0=0.d0
    integer :: ib,ik

    hubk = ubkwbk_hwbk(ubk,ubk,A0)
    Do ik = 1,NK
      Do ib = 1,NB
        epsbk(ib,ik) = real(sum(conjg(ubk(:,ib,ik))*hubk(:,ib,ik)))*dvcell
        hubk(:,ib,ik) = hubk(:,ib,ik) - epsbk(ib,ik)*ubk(:,ib,ik)
        epserrorbk(ib,ik) = sqrt(sum(abs(hubk(:,ib,ik)))*dvcell)
      End Do
    End Do
    
    return
  End Function ubk_epserrorbk
!==
End Module Functions
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
