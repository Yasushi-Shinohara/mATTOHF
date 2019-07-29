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
  Function ubk_dns(ubk) result(dns) !########Several 3D array is used as 1D array here for OMP.########
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,occbk
    implicit none
!    complex(8), intent(in) :: ubk(1:NG1,NG2,NG3,1:NB,1:NK)
!    real(8) :: dns(1:NG1,1:NG2,1:NG3)
    complex(8), intent(in) :: ubk(1:NG1*NG2*NG3,1:NB,1:NK)
    real(8) :: dns1D(1:NG1*NG2*NG3), dns(1:NG1,1:NG2,1:NG3)
    integer :: ib,ik,ig
    integer :: ig1,ig2,ig3
    
    dns1D = 0.d0
!$omp parallel
!$omp do private(ik,ib,ig)
    Do ig= 1,NG1*NG2*NG3
      Do ik= 1,NK
        Do ib =1,NB
          dns1D(ig) = dns1D(ig) + occbk(ib,ik)*conjg(ubk(ig,ib,ik))*ubk(ig,ib,ik)
!        dns(:,:,:) = dns(:,:,:) + occbk(ib,ik)*conjg(ubk(:,:,:,ib,ik))*ubk(:,:,:,ib,ik)
!        dns1D(:) = dns1D(:) + occbk(ib,ik)*conjg(ubk(:,ib,ik))*ubk(:,ib,ik)
        End do
      End do
    End do
!$omp end do
!$omp end parallel
    dns = reshape(dns1D,(/NG1,NG2,NG3/))
    return
  End Function ubk_dns
!==
!==
  Function ubk_tubk(ubk,A) result(tubk)
    use Global_Variables, only : NG1,NG2,NG3,NG, NB, NK, Gx,Gy,Gz, kx,ky,kz, planf,planb
    implicit none
    complex(8), intent(in) :: ubk(1:NG1,NG2,NG3,1:NB,1:NK)
    real(8), intent(in) :: A(1:3)
    complex(8) :: tubk(1:NG1,NG2,NG3,1:NB,1:NK)
    complex(8) :: work(1:NG1,NG2,NG3)
    real(8) :: kAx,kAy,kAz,kA2
    integer :: ib,ik
!FFTW staff
!    integer(8) :: planf,planb
    include 'fftw3.f'

!Make the plan
!    Call dfftw_plan_dft_3d(planf, NG1, NG2, NG3, ubk(:,:,:,1,1) , work,FFTW_FORWARD ,FFTW_ESTIMATE)
!    Call dfftw_plan_dft_3d(planb, NG1, NG2, NG3, work, tubk(:,:,:,1,1) ,FFTW_BACKWARD,FFTW_ESTIMATE)
    
!$omp parallel
!$omp do private(ik,kAx,kAy,kAz,kA2,ib,work)
    Do ik= 1,NK
      kAx = kx(ik)+A(1)
      kAy = ky(ik)+A(2)
      kAz = kz(ik)+A(3)
      kA2 = kAx**2+kAy**2+kAz**2
      Do ib =1,NB
        Call dfftw_execute_dft(planf,ubk(:,:,:,ib,ik),work)
        work(:,:,:) = 0.5d0*(Gx(:,:,:)**2 + Gy(:,:,:)**2 + Gz(:,:,:)**2 + &
             2.d0*(Gx(:,:,:)*kAx + Gy(:,:,:)*kAy + Gz(:,:,:)*kAz) + kA2)*work(:,:,:)/dble(NG)
!        work(:,:,:) = (Gx(:,:,:)**2+Gy(:,:,:)**2+Gz(:,:,:)**2+2.d0*(Gx(:,:,:)*kAx+Gy(:,:,:)*kAy+Gz(:,:,:)*kAz)+kA2)*work(:,:,:)/dble(NG)/2.d0
        Call dfftw_execute_dft(planb,work,tubk(:,:,:,ib,ik))
      End do
    End do
!$omp end do
!$omp end parallel
    return
  End Function ubk_tubk
!==
!==
  Function dns_vH(dns) result(vH)
    use Global_Variables, only : NG1,NG2,NG3,NG,fpiG2inv
    implicit none
    real(8), intent(in) :: dns(1:NG1,NG2,NG3)
    real(8) :: vH(1:NG1,NG2,NG3)
    complex(8) :: work1(1:NG1,NG2,NG3),work2(1:NG1,NG2,NG3)
!FFTW staff
    integer(8) :: planf,planb
    include 'fftw3.f'

!Make the plan
!    Call dfftw_plan_dft_r2c_3d(planf, NG1, NG2, NG3, dns , work,FFTW_FORWARD ,FFTW_ESTIMATE)
!    Call dfftw_plan_dft_c2r_3d(planb, NG1, NG2, NG3, work, vH ,FFTW_BACKWARD,FFTW_ESTIMATE)
    Call dfftw_plan_dft_3d(planf, NG1, NG2, NG3, work1, work2, FFTW_FORWARD ,FFTW_ESTIMATE)
    Call dfftw_plan_dft_3d(planb, NG1, NG2, NG3, work2, work1 ,FFTW_BACKWARD,FFTW_ESTIMATE)
    
    work1 = dns
    Call dfftw_execute_dft(planf,work1,work2)
    work2 = fpiG2inv*work2
    Call dfftw_execute_dft(planb,work2,work1)
    vH = dble(work1)/dble(NG)
    return
  End Function dns_vH
!==
!==
  Function ubk_vPSubk(ubk,A) result(vPSubk)
    use Constants
    use Global_Variables, only : NG1,NG2,NG3,NG,NBocc,NB,NK,a1,a2,a3,dvcell&
         ,a_tbl,Mps,Nps,Nion,uV,iuV,Nlma,Ips,Ipsn1,Ipsn2,Ipsn3&
         ,kx,ky,kz,rx,ry,rz
    implicit none
    complex(8), intent(in) :: ubk(1:NG1,NG2,NG3,1:NB,1:NK)
    real(8), intent(in) :: A(1:3)
    complex(8) :: vPSubk(1:NG1,NG2,NG3,1:NB,1:NK),uVpsi,ekr(Nps,Nion)
    real(8) :: kA(3),kr,x,y,z
    integer :: ib,ik
    integer :: ilma,ia,j,i1,i2,i3

    vPSubk=0.d0
!$omp parallel
!$omp do private(ik,kA,ia,j,i1,i2,i3,x,y,z,kr,ekr,ilma,uVpsi,ib)
    Do ik = 1,NK
      kA(1) = kx(ik)+A(1); kA(2) = ky(ik)+A(2); kA(3) = kz(ik)+A(3)
      Do ia=1,Nion
        Do j=1,Mps(ia)
          i1 = Ips(1,j,ia); i2 = Ips(2,j,ia); i3 = Ips(3,j,ia)
          x = rx(i1,i2,i3) - Ipsn1(j,ia)*a1(1) - Ipsn2(j,ia)*a2(1) - Ipsn3(j,ia)*a3(1)
          y = ry(i1,i2,i3) - Ipsn1(j,ia)*a1(2) - Ipsn2(j,ia)*a2(2) - Ipsn3(j,ia)*a3(2)
          z = rz(i1,i2,i3) - Ipsn1(j,ia)*a1(3) - Ipsn2(j,ia)*a2(3) - Ipsn3(j,ia)*a3(3)
          kr = kA(1)*x + kA(2)*y + kA(3)*z
          ekr(j,ia)=exp(zI*kr)
        End Do
      End Do
      Do ib = 1,NB
        Do ilma=1,Nlma
          ia=a_tbl(ilma)
          uVpsi=0.d0
          Do j=1,Mps(ia)
            i1 = Ips(1,j,ia); i2 = Ips(2,j,ia); i3 = Ips(3,j,ia)
            uVpsi=uVpsi+uV(j,ilma)*ekr(j,ia)*ubk(i1,i2,i3,ib,ik)
          End Do
          uVpsi=uVpsi*dvcell*iuV(ilma)
          Do j=1,Mps(ia)
            i1 = Ips(1,j,ia); i2 = Ips(2,j,ia); i3 = Ips(3,j,ia)
            vPSubk(i1,i2,i3,ib,ik) = vPSubk(i1,i2,i3,ib,ik) + conjg(ekr(j,ia))*uVpsi*uV(j,ilma)
          End Do
        End Do
      End Do
    End Do
!$omp end do
!$omp end parallel

    return
  End Function ubk_vPSubk
!==
!==
  Function ubk_vFubk(ubk) result(vfubk)
    use Global_Variables, only : NG1,NG2,NG3,NG,NB,NK,occbk,Gx,Gy,Gz,kx,ky,kz
    implicit none
    complex(8), intent(in) :: ubk(1:NG1,NG2,NG3,1:NB,1:NK)
    complex(8) :: vFubk(1:NG1,NG2,NG3,1:NB,1:NK)
    integer :: ib,ik

    vFubk = ubk

    return
  End Function ubk_vFubk
!==
!==
  Function ubkwbk_vFwbk(ubk,wbk) result(vFwbk)
    use Global_Variables, only : NG1,NG2,NG3,NG, NB,NBocc, NK,occbk, NKKs,NKKe,ikrtable,ikctable,fpiGkk2inv, planf,planb &
         , MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD,ierr, vcell
    implicit none
    complex(8), intent(in) :: ubk(1:NG1,NG2,NG3,1:NB,1:NK), wbk(1:NG1,NG2,NG3,1:NB,1:NK)
    complex(8) :: vFwbk(1:NG1,NG2,NG3,1:NB,1:NK),vFwbk_l(1:NG1,NG2,NG3,1:NB,1:NK)
    complex(8) :: pairrho (NG1,NG2,NG3),pairrhoG(NG1,NG2,NG3), VF(1:NG1,NG2,NG3)
    integer :: ibr,ibc,ikk

    vFwbk_l(:,:,:,:,:) = 0.d0
!$omp parallel
!!!$omp do private(ikk,ibr,ibc,pairrho,pairrhoG,VF) reduction(+:vFwbk_l) !Seems right but leading NaN
!$omp do private(ikk,ibr,ibc,pairrho,pairrhoG,VF)
    Do ikk = NKKs,NKKe
!      Do ibr = 1,NBocc
      Do ibr = 1,NB
        Do ibc = 1,NBocc
          pairrho(:,:,:) = wbk(:,:,:,ibr,ikrtable(ikk))*conjg(ubk(:,:,:,ibc,ikctable(ikk)))/dble(NK)
          Call dfftw_execute_dft(planf, pairrho, pairrhoG)
!          pairrhoG(:,:,:) = -fpiGkk2inv(:,:,:,ikrtable(ikk),ikctable(ikk))*pairrhoG(:,:,:)/dble(NG)
          pairrhoG(:,:,:) = -fpiGkk2inv(:,:,:,ikk)*pairrhoG(:,:,:)/dble(NG)
          Call dfftw_execute_dft(planb, pairrhoG, VF)
          vFwbk_l(:,:,:,ibr,ikrtable(ikk)) = vFwbk_l(:,:,:,ibr,ikrtable(ikk)) + VF(:,:,:)*ubk(:,:,:,ibc,ikctable(ikk))
        End Do
      End Do
    End Do
!$omp end do
!$omp end parallel
    Call MPI_ALLREDUCE(vFwbk_l,vFwbk,NG*NB*NK,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!    vFwbk = vFwbk/vcell !debug

    return
  End Function ubkwbk_vFwbk
!==
!==
  Function dns_vxcPZ(dns) result(vxc)
    use Constants
    use Global_Variables, only : NG1,NG2,NG3
    implicit none
    real(8), intent(in) :: dns(1:NG1,NG2,NG3)
    real(8) :: vxc(1:NG1,1:NG2,1:NG3)
    real(8) :: exc(1:NG1,1:NG2,1:NG3)    !XC-energy density
    real(8) :: dexc(1:NG1,1:NG2,1:NG3)   !d(exc)/d(dns)
    real(8) :: rs(1:NG1,1:NG2,1:NG3),rssq(1:NG1,1:NG2,1:NG3),rsln(1:NG1,1:NG2,1:NG3)
    integer :: i1,i2,i3

    rs(:,:,:) = (3.d0/dns(:,:,:)/fpi)**(1.d0/3.d0)
    exc(:,:,:) = -0.4582d0/rs(:,:,:)
    dexc(:,:,:) = exc(:,:,:)/dns(:,:,:)/3.d0
    rssq=sqrt(rs)
    rsln=log(rs)
    Do i3=1,NG3
    Do i2=1,NG2
    Do i1=1,NG1
      If (rs(i1,i2,i3)>1.d0) then
        exc(i1,i2,i3) = exc(i1,i2,i3) + gammaU/(1.d0 + beta1U*rssq(i1,i2,i3) + beta2U*rs(i1,i2,i3))
        dexc(i1,i2,i3) = dexc(i1,i2,i3) + gammaU*(0.5d0*beta1U*rssq(i1,i2,i3) + beta2U*rs(i1,i2,i3))&
             /(3.d0*dns(i1,i2,i3))/(1.d0 + beta1U*rssq(i1,i2,i3) + beta2U*rs(i1,i2,i3))**2
      Else
        exc(i1,i2,i3) = exc(i1,i2,i3) + AU*rsln(i1,i2,i3) + BU + CU*rs(i1,i2,i3)*rsln(i1,i2,i3) + DU*rs(i1,i2,i3)
        dexc(i1,i2,i3) = dexc(i1,i2,i3) - rs(i1,i2,i3)/(3.d0*dns(i1,i2,i3))&
             *(AU/rs(i1,i2,i3) + CU*(rsln(i1,i2,i3)+1.d0) + DU)
      End If
    End Do
    End Do
    End Do

    vxc(:,:,:) = exc(:,:,:) + dns(:,:,:)*dexc(:,:,:)
    return
  End Function dns_vxcPZ
!==
!==
  Function dns_vxcPBE(dns) result(vxc)
    use Global_Variables, only : NG1,NG2,NG3
    implicit none
    real(8), intent(in) :: dns(1:NG1,NG2,NG3)
    real(8) :: vxc(1:NG1,NG2,NG3)

    vxc = dns !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    return
  End Function dns_vxcPBE
!==
!==
  Function ubk_hubk(ubk,A) result(hubk)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,functional,vpsl
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: A(3)
    complex(8) :: hubk(NG1,NG2,NG3,NB,NK)
    real(8) :: dns(NG1,NG2,NG3),vH(NG1,NG2,NG3),vxc(NG1,NG2,NG3)
    integer :: ib,ik

    hubk = ubk_tubk(ubk,A)
    dns = ubk_dns(ubk)
    vH = dns_vH(dns)
    Do ik = 1,NK
      Do ib = 1,NB
        hubk(:,:,:,ib,ik) = hubk(:,:,:,ib,ik) + (vH(:,:,:)+vpsl(:,:,:))*ubk(:,:,:,ib,ik)
      End Do
    End Do
    hubk = hubk + ubk_vPSubk(ubk,A)
    Select case (trim(functional))
      Case('PZ')
        vxc = dns_vxcPZ(dns)
        Do ik = 1,NK
          Do ib = 1,NB
            hubk(:,:,:,ib,ik) = hubk(:,:,:,ib,ik) + vxc(:,:,:)*ubk(:,:,:,ib,ik)
          End Do
        End Do
      Case('PBE')
        vxc = dns_vxcPBE(dns)
        Do ik = 1,NK
          Do ib = 1,NB
            hubk(:,:,:,ib,ik) = hubk(:,:,:,ib,ik) + vxc(:,:,:)*ubk(:,:,:,ib,ik)
          End Do
        End Do
      Case('HF')
        hubk = hubk + ubk_vFubk(ubk)
    End Select
    return
  End Function ubk_hubk
!==
!==
!Hamiltonian operation to ubk with given dns
  Function dnsubk_hubk(dns,ubk,A) result(hubk)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,functional,vpsl
    implicit none
    real(8), intent(in) :: dns(NG1,NG2,NG3)
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: A(3)
    complex(8) :: hubk(NG1,NG2,NG3,NB,NK)
    real(8) :: vH(NG1,NG2,NG3),vxc(NG1,NG2,NG3)
    integer :: ib,ik

    hubk = ubk_tubk(ubk,A)
    vH = dns_vH(dns)
    Do ik = 1,NK
      Do ib = 1,NB
        hubk(:,:,:,ib,ik) = hubk(:,:,:,ib,ik) + (vH(:,:,:)+vpsl(:,:,:))*ubk(:,:,:,ib,ik)
      End Do
    End Do
    hubk = hubk + ubk_vPSubk(ubk,A)
    Select case (trim(functional))
      Case('PZ')
        vxc = dns_vxcPZ(dns)
        Do ik = 1,NK
          Do ib = 1,NB
            hubk(:,:,:,ib,ik) = hubk(:,:,:,ib,ik) + vxc(:,:,:)*ubk(:,:,:,ib,ik)
          End Do
        End Do
      Case('PBE')
        vxc = dns_vxcPBE(dns)
        Do ik = 1,NK
          Do ib = 1,NB
            hubk(:,:,:,ib,ik) = hubk(:,:,:,ib,ik) + vxc(:,:,:)*ubk(:,:,:,ib,ik)
          End Do
        End Do
      Case('HF')
        hubk = hubk + ubk_vFubk(ubk)
    End Select
    return
  End Function dnsubk_hubk
!==
!==
  Function ubk_Etot(ubk,A) result(Etot)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,functional
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: A(3)
    real(8) :: Etot, Ekin, Epsl, Epsnl, EH, Exc, dns(NG1,NG2,NG3)
    
    dns = ubk_dns(ubk)
    Ekin = ubk_Ekin(ubk,A)
    Epsl = dns_Epsl(dns)
    Epsnl = ubk_Epsnl(ubk,A)
    EH = dns_EH(dns)
    Select case (trim(functional))
      Case('PZ')
        Exc = dns_ExcPZ(dns)
      Case('PBE')
        Exc = dns_ExcPBE(dns)
      Case('HF')
        Exc = ubk_Ex(ubk)
    End Select
    Etot = Ekin + Epsl + Epsnl + EH + Exc
    return
  End Function ubk_Etot
!==
!==
!Energy computation w.r.t. each components
  Function ubk_Ecomps(ubk,A) result(Ecomps)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,functional
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: A(3)
    real(8) :: Ecomps(7), Ekin, Epsl, Epsnl, EH, Exc, Eion, Etot, dns(NG1,NG2,NG3)
    
    dns = ubk_dns(ubk)
    Ekin = ubk_Ekin(ubk,A)
    Epsl = dns_Epsl(dns)
    Epsnl = ubk_Epsnl(ubk,A)
    EH = dns_EH(dns)
    Select case (trim(functional))
      Case('PZ')
        Exc = dns_ExcPZ(dns)
      Case('PBE')
        Exc = dns_ExcPBE(dns)
      Case('HF')
        Exc = ubk_Ex(ubk)
    End Select
    Eion = 0.d0
    Etot = Ekin + Epsl + Epsnl + EH + Exc + Eion
    Ecomps(1) = Ekin
    Ecomps(2) = Epsl
    Ecomps(3) = Epsnl
    Ecomps(4) = EH
    Ecomps(5) = Exc
    Ecomps(6) = Eion
    Ecomps(7) = Etot
    return
  End Function ubk_Ecomps
!==
!==
  Function ubk_Ekin(ubk,A) result(Ekin)
    use Global_Variables, only : NG1,NG2,NG3,NG,NB,NK,dvcell,Gx,Gy,Gz,kx,ky,kz,occbk, planf, planb
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: A(3)
    real(8) :: Ekin
    complex(8) :: work(NG1,NG2,NG3),work2(NG1,NG2,NG3)
    real(8) :: kAx,kAy,kAz,kA2
    integer :: ib,ik
!FFTW staff
!    integer(8) :: planf,planb
    include 'fftw3.f'

!Make the plan
!    Call dfftw_plan_dft_3d(planf, NG1, NG2, NG3, ubk(:,:,:,1,1) , work,FFTW_FORWARD ,FFTW_ESTIMATE)
!    Call dfftw_plan_dft_3d(planb, NG1, NG2, NG3, work, work2(:,:,:) ,FFTW_BACKWARD,FFTW_ESTIMATE)

    Ekin = 0.d0
!$omp parallel
!$omp do private(ik,kAx,kAy,kAz,kA2,ib,work,work2) reduction(+:Ekin)
    Do ik= 1,NK
      kAx = kx(ik)+A(1)
      kAy = ky(ik)+A(2)
      kAz = kz(ik)+A(3)
      kA2 = kAx**2+kAy**2+kAz**2
      Do ib =1,NB
        Call dfftw_execute_dft(planf,ubk(:,:,:,ib,ik),work)
        work(:,:,:) = 0.5d0*(Gx(:,:,:)**2+Gy(:,:,:)**2+Gz(:,:,:)**2 + &
             2.d0*(Gx(:,:,:)*kAx + Gy(:,:,:)*kAy + Gz(:,:,:)*kAz) + kA2)*work(:,:,:)/dble(NG)
        Call dfftw_execute_dft(planb,work,work2(:,:,:))
        Ekin = Ekin + occbk(ib,ik)*sum(conjg(ubk(:,:,:,ib,ik))*work2(:,:,:))*dvcell
      End do
    End do
!$omp end do
!$omp end parallel
    
    return
  End Function ubk_Ekin
!==
!==
  Function dns_Epsl(dns) result(Epsl)
    use Global_Variables, only : NG1,NG2,NG3,vpsl,dvcell
    implicit none
    real(8), intent(in) :: dns(NG1,NG2,NG3)
    real(8) :: Epsl

    Epsl = sum(vpsl(:,:,:)*dns(:,:,:))*dvcell
    
    return
  End Function dns_Epsl
!==
!==
  Function ubk_Epsnl(ubk,A) result(Epsnl)
    use Constants
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,rx,ry,rz,a1,a2,a3,kx,ky,kz&
         ,Ips, Ipsn1, Ipsn2, Ipsn3, Nion, Nps, Nlma, Mps, a_tbl, uV, iuV&
         ,NBocc, occbk, dvcell
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: A(3)
    real(8) :: Epsnl,kA(3),kr,x,y,z
    complex(8) :: uVpsi, ekr(Nps,Nion)
    integer :: ib,ik,ia,j,ilma,i1,i2,i3

    Epsnl = 0.d0
!$omp parallel
!$omp do private(ik,kA,ia,j,i1,i2,i3,x,y,z,kr,ekr,ib,ilma,uVpsi) reduction(+:Epsnl)
    Do ik = 1,NK
      kA(1) = kx(ik)+A(1); kA(2) = ky(ik)+A(2); kA(3) = kz(ik)+A(3)
      Do ia=1,Nion
        Do j=1,Mps(ia)
          i1 = Ips(1,j,ia); i2 = Ips(2,j,ia); i3 = Ips(3,j,ia)
          x = rx(i1,i2,i3) - Ipsn1(j,ia)*a1(1) - Ipsn2(j,ia)*a2(1) - Ipsn3(j,ia)*a1(1)
          y = ry(i1,i2,i3) - Ipsn1(j,ia)*a1(2) - Ipsn2(j,ia)*a2(2) - Ipsn3(j,ia)*a1(2)
          z = rz(i1,i2,i3) - Ipsn1(j,ia)*a1(3) - Ipsn2(j,ia)*a2(3) - Ipsn3(j,ia)*a1(3)
          kr = kA(1)*x + kA(2)*y + kA(3)*z
          ekr(j,ia)=exp(zI*kr)
        End Do
      End Do
      Do ib = 1,NB
        Do ilma=1,Nlma
          ia=a_tbl(ilma)
          uVpsi=0.d0
          Do j=1,Mps(ia)
            i1 = Ips(1,j,ia); i2 = Ips(2,j,ia); i3 = Ips(3,j,ia)
            uVpsi=uVpsi+uV(j,ilma)*ekr(j,ia)*ubk(i1,i2,i3,ib,ik)
          End Do
          uVpsi=uVpsi*dvcell
          Epsnl = Epsnl + occbk(ib,ik)*abs(uVpsi)**2*iuV(ilma)
        End Do
      End Do
    End Do
!$omp end do
!$omp end parallel
    
    return
  End Function ubk_Epsnl
!==
!==
  Function dns_EH(dns) result(EH)
    use Global_Variables, only : NG1,NG2,NG3,dvcell
    implicit none
    real(8), intent(in) :: dns(NG1,NG2,NG3)
    real(8) :: vH(NG1,NG2,NG3),EH

    vH = dns_vH(dns)
    EH = 0.5d0*sum(vH(:,:,:)*dns(:,:,:))*dvcell
    
    return
  End Function dns_EH
!==
!==
  Function dns_ExcPZ(dns) result(Exc_l)
    use Constants
    use Global_Variables, only : NG1,NG2,NG3,dvcell
    implicit none
    real(8), intent(in) :: dns(NG1,NG2,NG3)
    real(8) :: exc(NG1,NG2,NG3),Exc_l !XC energy density and its spatial average
    real(8) :: rs(1:NG1,1:NG2,1:NG3),rssq(1:NG1,1:NG2,1:NG3),rsln(1:NG1,1:NG2,1:NG3)
    integer :: i1,i2,i3

    rs(:,:,:) = (3.d0/dns(:,:,:)/fpi)**(1.d0/3.d0)
    exc(:,:,:) = -0.4582d0/rs(:,:,:)
    rssq=sqrt(rs)
    rsln=log(rs)
    Do i3=1,NG3
    Do i2=1,NG2
    Do i1=1,NG1
      If (rs(i1,i2,i3)>1.d0) then
        exc(i1,i2,i3) = exc(i1,i2,i3) + gammaU/(1.d0+beta1U*rssq(i1,i2,i3) + beta2U*rs(i1,i2,i3))
      Else
        exc(i1,i2,i3) = exc(i1,i2,i3) + AU*rsln(i1,i2,i3) + BU  + CU*rs(i1,i2,i3)*rsln(i1,i2,i3) + DU*rs(i1,i2,i3)
      End If
    End Do
    End Do
    End Do

    Exc_l = sum(dns(:,:,:)*exc(:,:,:))*dvcell
    
    return
  End Function dns_ExcPZ
!==
!==
  Function dns_ExcPBE(dns) result(Exc_l)
    use Global_Variables, only : NG1,NG2,NG3
    implicit none
    real(8), intent(in) :: dns(NG1,NG2,NG3)
    real(8) :: exc(NG1,NG2,NG3),Exc_l

    Exc_l = 0.d0
    
    return
  End Function dns_ExcPBE
!==
!==
  Function ubk_Ex(ubk) result(Ex)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK, NBocc, occbk, dvcell
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    complex(8) :: vFubk(NG1,NG2,NG3,NB,NK)
    real(8) :: Ex
    integer :: ib,ik

    vFubk = ubkwbk_vFwbk(ubk,ubk)
    Ex = 0.d0
!$omp parallel
!$omp do private(ik,ib) reduction(+:Ex)
    Do ik = 1,NK
      Do ib = 1,NBocc
        Ex = Ex + occbk(ib,ik)*sum(conjg(ubk(:,:,:,ib,ik))*vFubk(:,:,:,ib,ik))*dvcell
      End Do
    End Do
!$omp end do
!$omp end parallel
    Ex = 0.5d0*Ex
    
    return
  End Function ubk_Ex
!==
!==
!Total energy from ubk with given dns
  Function dnsubk_Etot(dns,ubk,A) result(Etot)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,functional
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: dns(NG1,NG2,NG3)
    real(8), intent(in) :: A(3)
    real(8) :: Etot, Ekin, Epsl, Epsnl, EH, Exc
    
    Ekin = ubk_Ekin(ubk,A)
    Epsl = dns_Epsl(dns)
    Epsnl = ubk_Epsnl(ubk,A)
    EH = dns_EH(dns)
    Select case (trim(functional))
      Case('PZ')
        Exc = dns_ExcPZ(dns)
      Case('PBE')
        Exc = dns_ExcPBE(dns)
      Case('HF')
        Exc = ubk_Ex(ubk)
    End Select
    Etot = Ekin + Epsl + Epsnl + EH + Exc
    return
  End Function dnsubk_Etot
!==
!==
  Function dnsubk_epsbk(dns,ubk,A0) result(epsbk)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,dvcell
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: A0(3)
    real(8), intent(in) :: dns(NG1,NG2,NG3)
    complex(8) :: hubk(NG1,NG2,NG3,NB,NK)
    real(8) :: epsbk(NB,NK)
    integer :: ib,ik

    hubk = dnsubk_hubk(dns,ubk,A0)
    Do ik = 1,NK
      Do ib = 1,NB
        epsbk(ib,ik) = real(sum(conjg(ubk(:,:,:,ib,ik))*hubk(:,:,:,ib,ik)))*dvcell
      End Do
    End Do
    
    return
  End Function dnsubk_epsbk
!==
!==
  Function dnsubk_epserrorbk(dns,ubk,A0) result(epserrorbk)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,dvcell
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: A0(3)
    real(8), intent(in) :: dns(NG1,NG2,NG3)
    complex(8) :: hubk(NG1,NG2,NG3,NB,NK)
    real(8) :: epsbk(NB,NK),epserrorbk(NB,NK)
    integer :: ib,ik

    hubk = dnsubk_hubk(dns,ubk,A0)
    Do ik = 1,NK
      Do ib = 1,NB
        epsbk(ib,ik) = real(sum(conjg(ubk(:,:,:,ib,ik))*hubk(:,:,:,ib,ik)))*dvcell
        hubk(:,:,:,ib,ik) = hubk(:,:,:,ib,ik) - epsbk(ib,ik)*ubk(:,:,:,ib,ik)
        epserrorbk(ib,ik) = sqrt(sum(abs(hubk(:,:,:,ib,ik)))*dvcell)
      End Do
    End Do
    
    return
  End Function dnsubk_epserrorbk
!==
!==
!Hamiltonian operation to ub with given dns for each k-index
  Function dnsub_hub(dns,ub,ik,A) result(hub)
    use Global_Variables, only : NG1,NG2,NG3,NB,functional,vpsl
    implicit none
    real(8), intent(in) :: dns(NG1,NG2,NG3)
    complex(8), intent(in) :: ub(NG1,NG2,NG3,NB)
    integer, intent(in) :: ik
    real(8), intent(in) :: A(3)
    complex(8) :: hub(NG1,NG2,NG3,NB)
    real(8) :: vH(NG1,NG2,NG3),vxc(NG1,NG2,NG3)
    integer :: ib

    hub = ub_tub(ub,ik,A)
    vH = dns_vH(dns)
    Do ib = 1,NB
      hub(:,:,:,ib) = hub(:,:,:,ib) + (vH(:,:,:)+vpsl(:,:,:))*ub(:,:,:,ib)
    End Do
    hub = hub + ub_vPSub(ub,ik,A)
    Select case (trim(functional))
      Case('PZ')
        vxc = dns_vxcPZ(dns)
        Do ib = 1,NB
          hub(:,:,:,ib) = hub(:,:,:,ib) + vxc(:,:,:)*ub(:,:,:,ib)
      End Do
      Case('PBE')
        vxc = dns_vxcPBE(dns)
        Do ib = 1,NB
          hub(:,:,:,ib) = hub(:,:,:,ib) + vxc(:,:,:)*ub(:,:,:,ib)
        End Do
      Case('HF')
        stop ('HF can not be applied with dnsub_hub function.')
    End Select
    return
  End Function dnsub_hub
!==
!==
  Function ub_tub(ub,ik,A) result(tub)
    use Global_Variables, only : NG1,NG2,NG3,NG,NB,Gx,Gy,Gz,kx,ky,kz
    implicit none
    complex(8), intent(in) :: ub(1:NG1,NG2,NG3,1:NB)
    integer, intent(in) :: ik
    real(8), intent(in) :: A(1:3)
    complex(8) :: tub(1:NG1,NG2,NG3,1:NB)
    complex(8) :: work(1:NG1,NG2,NG3)
    real(8) :: kAx,kAy,kAz,kA2
    integer :: ib
!FFTW staff
    integer(8) :: planf,planb
    include 'fftw3.f'

!Make the plan
    Call dfftw_plan_dft_3d(planf, NG1, NG2, NG3, ub(:,:,:,1) , work,FFTW_FORWARD ,FFTW_ESTIMATE)
    Call dfftw_plan_dft_3d(planb, NG1, NG2, NG3, work, tub(:,:,:,1) ,FFTW_BACKWARD,FFTW_ESTIMATE)
    
    kAx = kx(ik)+A(1)
    kAy = ky(ik)+A(2)
    kAz = kz(ik)+A(3)
    kA2 = kAx**2+kAy**2+kAz**2
    Do ib =1,NB
      Call dfftw_execute_dft(planf,ub(:,:,:,ib),work)
      work(:,:,:) = 0.5d0*(Gx(:,:,:)**2+Gy(:,:,:)**2+Gz(:,:,:)**2 + &
           2.d0*(Gx(:,:,:)*kAx + Gy(:,:,:)*kAy + Gz(:,:,:)*kAz) + kA2)*work(:,:,:)/dble(NG)
      Call dfftw_execute_dft(planb,work,tub(:,:,:,ib))
    End do
    return
  End Function ub_tub
!==
!==
  Function ub_vPSub(ub,ik,A) result(vPSub)
    use Constants
    use Global_Variables, only : NG1,NG2,NG3,NG,NB,a1,a2,a3,dvcell&
         ,a_tbl,Mps,Nps,Nion,uV,iuV,Nlma,Ips,Ipsn1,Ipsn2,Ipsn3&
         ,kx,ky,kz,rx,ry,rz
    implicit none
    complex(8), intent(in) :: ub(1:NG1,NG2,NG3,1:NB)
    integer, intent(in) :: ik
    real(8), intent(in) :: A(1:3)
    complex(8) :: vPSub(1:NG1,NG2,NG3,1:NB),uVpsi,ekr(Nps,Nion)
    real(8) :: kA(3),kr,x,y,z
    integer :: ib
    integer :: ilma,ia,j,i1,i2,i3

    vPSub=0.d0
    kA(1) = kx(ik)+A(1); kA(2) = ky(ik)+A(2); kA(3) = kz(ik)+A(3)
    Do ia=1,Nion
      Do j=1,Mps(ia)
        i1 = Ips(1,j,ia); i2 = Ips(2,j,ia); i3 = Ips(3,j,ia)
        x = rx(i1,i2,i3) - Ipsn1(j,ia)*a1(1) - Ipsn2(j,ia)*a2(1) - Ipsn3(j,ia)*a3(1)
        y = ry(i1,i2,i3) - Ipsn1(j,ia)*a1(2) - Ipsn2(j,ia)*a2(2) - Ipsn3(j,ia)*a3(2)
        z = rz(i1,i2,i3) - Ipsn1(j,ia)*a1(3) - Ipsn2(j,ia)*a2(3) - Ipsn3(j,ia)*a3(3)
        kr = kA(1)*x + kA(2)*y + kA(3)*z
        ekr(j,ia)=exp(zI*kr)
      End Do
    End Do
    Do ib = 1,NB
      Do ilma=1,Nlma
        ia=a_tbl(ilma)
        uVpsi=0.d0
        Do j=1,Mps(ia)
          i1 = Ips(1,j,ia); i2 = Ips(2,j,ia); i3 = Ips(3,j,ia)
          uVpsi=uVpsi+uV(j,ilma)*ekr(j,ia)*ub(i1,i2,i3,ib)
        End Do
        uVpsi=uVpsi*dvcell*iuV(ilma)
        Do j=1,Mps(ia)
          i1 = Ips(1,j,ia); i2 = Ips(2,j,ia); i3 = Ips(3,j,ia)
          vPSub(i1,i2,i3,ib) = vPSub(i1,i2,i3,ib) + conjg(ekr(j,ia))*uVpsi*uV(j,ilma)
        End Do
      End Do
    End Do

    return
  End Function ub_vPSub
!==
!Hamiltonian operation to ub with given dns for each k-index: Thread Safe version
  Function dnsub_hub_TS(dns,ub,ik,A) result(hub)
    use Global_Variables, only : NG1,NG2,NG3,NB,functional,vpsl, planf,planb
    implicit none
    real(8), intent(in) :: dns(NG1,NG2,NG3)
    complex(8), intent(in) :: ub(NG1,NG2,NG3,NB)
    integer, intent(in) :: ik
    real(8), intent(in) :: A(3)
    complex(8) :: hub(NG1,NG2,NG3,NB)
    real(8) :: vH(NG1,NG2,NG3),vxc(NG1,NG2,NG3)
    integer :: ib

    hub = ub_tub_TS(ub,ik,A)
    vH = dns_vH_TS(dns)
    Do ib = 1,NB
      hub(:,:,:,ib) = hub(:,:,:,ib) + (vH(:,:,:)+vpsl(:,:,:))*ub(:,:,:,ib)
    End Do
    hub = hub + ub_vPSub(ub,ik,A)
    Select case (trim(functional))
      Case('PZ')
        vxc = dns_vxcPZ(dns)
        Do ib = 1,NB
          hub(:,:,:,ib) = hub(:,:,:,ib) + vxc(:,:,:)*ub(:,:,:,ib)
      End Do
      Case('PBE')
        vxc = dns_vxcPBE(dns)
        Do ib = 1,NB
          hub(:,:,:,ib) = hub(:,:,:,ib) + vxc(:,:,:)*ub(:,:,:,ib)
        End Do
      Case('HF')
        stop ('HF can not be applied with dnsub_hub function.')
    End Select
    return
  End Function dnsub_hub_TS
!==
!==
!Thread Safe version
  Function ub_tub_TS(ub,ik,A) result(tub)
    use Global_Variables, only : NG1,NG2,NG3,NG,NB,NK,Gx,Gy,Gz,kx,ky,kz, planf,planb
    implicit none
    complex(8), intent(in) :: ub(1:NG1,NG2,NG3,1:NB)
    integer, intent(in) :: ik
    real(8), intent(in) :: A(1:3)
    complex(8) :: tub(1:NG1,NG2,NG3,1:NB)
    complex(8) :: work(1:NG1,NG2,NG3)
    real(8) :: kAx,kAy,kAz,kA2
    integer :: ib 
!FFTW staff
    include 'fftw3.f'

    kAx = kx(ik)+A(1)
    kAy = ky(ik)+A(2)
    kAz = kz(ik)+A(3)
    kA2 = kAx**2+kAy**2+kAz**2
    Do ib =1,NB
      Call dfftw_execute_dft(planf,ub(:,:,:,ib),work)
      work(:,:,:) = 0.5d0*(Gx(:,:,:)**2+Gy(:,:,:)**2+Gz(:,:,:)**2 + &
           2.d0*(Gx(:,:,:)*kAx + Gy(:,:,:)*kAy + Gz(:,:,:)*kAz) + kA2)*work(:,:,:)/dble(NG)
      Call dfftw_execute_dft(planb,work,tub(:,:,:,ib))
    End do
    return
  End Function ub_tub_TS
!==
!==
  Function dns_vH_TS(dns) result(vH)
    use Global_Variables, only : NG1,NG2,NG3,NG,fpiG2inv, planf,planb
    implicit none
    real(8), intent(in) :: dns(1:NG1,NG2,NG3)
    real(8) :: vH(1:NG1,NG2,NG3)
    complex(8) :: work1(1:NG1,NG2,NG3),work2(1:NG1,NG2,NG3)
!FFTW staff
    include 'fftw3.f'

    work1 = dns
    Call dfftw_execute_dft(planf,work1,work2)
    work2 = fpiG2inv*work2
    Call dfftw_execute_dft(planb,work2,work1)
    vH = dble(work1)/dble(NG)
    return
  End Function dns_vH_TS
!==
!==
!Hamiltonian operation to wbk with given ubk
  Function ubkwbk_hwbk(ubk,wbk,A) result(hwbk)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,functional,vpsl
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK), wbk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: A(3)
    complex(8) :: hwbk(NG1,NG2,NG3,NB,NK)
    real(8) :: vH(NG1,NG2,NG3),vxc(NG1,NG2,NG3),dns(NG1,NG2,NG3)
    integer :: ib,ik

    hwbk = ubk_tubk(wbk,A)
    dns = ubk_dns(ubk)
    vH = dns_vH(dns)
!$omp parallel
!$omp do private(ik,ib)
    Do ik = 1,NK
      Do ib = 1,NB
        hwbk(:,:,:,ib,ik) = hwbk(:,:,:,ib,ik) + (vH(:,:,:)+vpsl(:,:,:))*wbk(:,:,:,ib,ik)
      End Do
    End Do
!$omp end do
!$omp end parallel
    hwbk = hwbk + ubk_vPSubk(wbk,A)
    Select case (trim(functional))
      Case('PZ')
        vxc = dns_vxcPZ(dns)
!$omp parallel
!$omp do private(ik,ib)
        Do ik = 1,NK
          Do ib = 1,NB
            hwbk(:,:,:,ib,ik) = hwbk(:,:,:,ib,ik) + vxc(:,:,:)*wbk(:,:,:,ib,ik)
          End Do
        End Do
!$omp end do
!$omp end parallel
      Case('PBE')
        vxc = dns_vxcPBE(dns)
        Do ik = 1,NK
          Do ib = 1,NB
            hwbk(:,:,:,ib,ik) = hwbk(:,:,:,ib,ik) + vxc(:,:,:)*wbk(:,:,:,ib,ik)
          End Do
        End Do
      Case('HF')
        hwbk = hwbk + ubkwbk_vFwbk(ubk,wbk)
    End Select
    return
  End Function ubkwbk_hwbk
!==
!==
  Function ubk_epsbk(ubk) result(epsbk)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,dvcell
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), parameter :: A0(3)=0.d0
    complex(8) :: hubk(NG1,NG2,NG3,NB,NK)
    real(8) :: epsbk(NB,NK)
    integer :: ib,ik

    hubk = ubkwbk_hwbk(ubk,ubk,A0)
    Do ik = 1,NK
      Do ib = 1,NB
        epsbk(ib,ik) = real(sum(conjg(ubk(:,:,:,ib,ik))*hubk(:,:,:,ib,ik)))*dvcell
      End Do
    End Do
    
    return
  End Function ubk_epsbk
!==
!==
  Function ubk_epserrorbk(ubk) result(epserrorbk)
    use Global_Variables, only : NG1,NG2,NG3,NB,NK,dvcell
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), parameter :: A0(3)=0.d0
    complex(8) :: hubk(NG1,NG2,NG3,NB,NK)
    real(8) :: epsbk(NB,NK),epserrorbk(NB,NK)
    integer :: ib,ik

    hubk = ubkwbk_hwbk(ubk,ubk,A0)
    Do ik = 1,NK
      Do ib = 1,NB
        epsbk(ib,ik) = real(sum(conjg(ubk(:,:,:,ib,ik))*hubk(:,:,:,ib,ik)))*dvcell
        hubk(:,:,:,ib,ik) = hubk(:,:,:,ib,ik) - epsbk(ib,ik)*ubk(:,:,:,ib,ik)
        epserrorbk(ib,ik) = sqrt(sum(abs(hubk(:,:,:,ib,ik)))*dvcell)
      End Do
    End Do
    
    return
  End Function ubk_epserrorbk
!==
!==
  Function ubk_Je(ubk,A) result(Je)
    use Constants
    use Global_Variables, only : NG1,NG2,NG3,NG, NB,NK, dvcell,vcell, Gx,Gy,Gz,GxfJ,GyfJ,GzfJ,kx,ky,kz, occbk, planf, planb &
         ,a_tbl,Mps,Nps,Nion,uV,iuV,Nlma,Ips,Ipsn1,Ipsn2,Ipsn3 &
         ,a1,a2,a3, rx,ry,rz
    implicit none
    complex(8), intent(in) :: ubk(NG1,NG2,NG3,NB,NK)
    real(8), intent(in) :: A(3)
    complex(8) :: JpkA(3), Jps(3)
    complex(8) :: JpkAx, JpkAy, JpkAz, Jpsx, Jpsy, Jpsz
    real(8) :: Je(3)
    complex(8) :: work1(NG1,NG2,NG3),work2(NG1,NG2,NG3)
    real(8) :: GkAx(NG1,NG2,NG3),GkAy(NG1,NG2,NG3),GkAz(NG1,NG2,NG3)
    integer :: ib,ik
    complex(8) :: ekr(Nps,Nion), uVpsi,uVpsix,uVpsiy,uVpsiz
    real(8) :: kA(3),kr,x,y,z
    integer :: ilma,ia,j,i1,i2,i3
!FFTW staff
    include 'fftw3.f'

    JpkAx = 0.d0;    JpkAy = 0.d0;    JpkAz = 0.d0
!$omp parallel
!$omp do private(ik,GkAx,GkAy,GkAz,ib,work1,work2) reduction(+:JpkAx,JpkAy,JpkAz)
    Do ik = 1,NK
!      GkAx(:,:,:) = GxfJ(:,:,:) + kx(ik) + A(1)
!      GkAy(:,:,:) = GyfJ(:,:,:) + ky(ik) + A(2)
!      GkAz(:,:,:) = GzfJ(:,:,:) + kz(ik) + A(3)
      GkAx(:,:,:) = Gx(:,:,:) + kx(ik) + A(1)
      GkAy(:,:,:) = Gy(:,:,:) + ky(ik) + A(2)
      GkAz(:,:,:) = Gz(:,:,:) + kz(ik) + A(3)
      Do ib = 1,NB
        Call dfftw_execute_dft(planf,ubk(:,:,:,ib,ik),work1)
        work1 = work1/dble(NG)
        Call dfftw_execute_dft(planb,GkAx(:,:,:)*work1(:,:,:),work2)
        JpkAx = JpkAx + occbk(ib,ik)*sum(conjg(ubk(:,:,:,ib,ik))*work2(:,:,:))
        Call dfftw_execute_dft(planb,GkAy(:,:,:)*work1(:,:,:),work2)
        JpkAy = JpkAy + occbk(ib,ik)*sum(conjg(ubk(:,:,:,ib,ik))*work2(:,:,:))
        Call dfftw_execute_dft(planb,GkAz(:,:,:)*work1(:,:,:),work2)
        JpkAz = JpkAz + occbk(ib,ik)*sum(conjg(ubk(:,:,:,ib,ik))*work2(:,:,:))
      End do
    End do
!$omp end do
!$omp end parallel
    JpkA(1) = JpkAx*dvcell/vcell
    JpkA(2) = JpkAy*dvcell/vcell
    JpkA(3) = JpkAz*dvcell/vcell

    Jpsx = 0.d0;    Jpsy = 0.d0;    Jpsz = 0.d0
!$omp parallel
!$omp do private(ik,kA,ia,j,i1,i2,i3,x,y,z,kr,ekr,ib,ilma,uVpsi,uVpsix,uVpsiy,uVpsiz) reduction(+:Jpsx,Jpsy,Jpsz)
    do ik = 1,NK
      kA(1) = kx(ik)+A(1); kA(2) = ky(ik)+A(2); kA(3) = kz(ik)+A(3)
      do ia = 1,Nion
        do j = 1,Mps(ia)
          i1 = Ips(1,j,ia); i2 = Ips(2,j,ia); i3 = Ips(3,j,ia)
          x = rx(i1,i2,i3) - Ipsn1(j,ia)*a1(1) - Ipsn2(j,ia)*a2(1) - Ipsn3(j,ia)*a3(1)
          y = ry(i1,i2,i3) - Ipsn1(j,ia)*a1(2) - Ipsn2(j,ia)*a2(2) - Ipsn3(j,ia)*a3(2)
          z = rz(i1,i2,i3) - Ipsn1(j,ia)*a1(3) - Ipsn2(j,ia)*a2(3) - Ipsn3(j,ia)*a3(3)
          kr = kA(1)*x + kA(2)*y + kA(3)*z
          ekr(j,ia) = exp(zI*kr)
!          i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
!          kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
!          ekr(j,ia)=exp(zI*kr)
       enddo
      end do
      do ib = 1,NB
        do ilma = 1,Nlma
          ia = a_tbl(ilma)
          uVpsi = 0.d0; uVpsix = 0.d0; uVpsiy = 0.d0; uVpsiz = 0.d0
          do j = 1,Mps(ia)
            i1 = Ips(1,j,ia); i2 = Ips(2,j,ia); i3 = Ips(3,j,ia)
            x = rx(i1,i2,i3) - Ipsn1(j,ia)*a1(1) - Ipsn2(j,ia)*a2(1) - Ipsn3(j,ia)*a3(1)
            y = ry(i1,i2,i3) - Ipsn1(j,ia)*a1(2) - Ipsn2(j,ia)*a2(2) - Ipsn3(j,ia)*a3(2)
            z = rz(i1,i2,i3) - Ipsn1(j,ia)*a1(3) - Ipsn2(j,ia)*a2(3) - Ipsn3(j,ia)*a3(3)
            uVpsi  = uVpsi  + uV(j,ilma)*ekr(j,ia)  *ubk(i1,i2,i3,ib,ik)
            uVpsix = uVpsix + uV(j,ilma)*ekr(j,ia)*x*ubk(i1,i2,i3,ib,ik)
            uVpsiy = uVpsiy + uV(j,ilma)*ekr(j,ia)*y*ubk(i1,i2,i3,ib,ik)
            uVpsiz = uVpsiz + uV(j,ilma)*ekr(j,ia)*z*ubk(i1,i2,i3,ib,ik)
!            i = Jxyz(j,ia); ix = Jxx(j,ia); iy = Jyy(j,ia); iz = Jzz(j,ia)
!            x = Lx(i)*Hx - ix*aLx
!            y = Ly(i)*Hy - iy*aLy
!            z = Lz(i)*Hz - iz*aLz
!            uVpsi  = uVpsi  + uV(j,ilma)*ekr(j,ia)  *zu(i,ib,ik)
!            uVpsix = uVpsix + uV(j,ilma)*ekr(j,ia)*x*zu(i,ib,ik)
!            uVpsiy = uVpsiy + uV(j,ilma)*ekr(j,ia)*y*zu(i,ib,ik)
!            uVpsiz = uVpsiz + uV(j,ilma)*ekr(j,ia)*z*zu(i,ib,ik)
          end do
          uVpsi  = uVpsi *dvcell*iuV(ilma)
          uVpsix = uVpsix*dvcell
          uVpsiy = uVpsiy*dvcell
          uVpsiz = uVpsiz*dvcell
          Jpsx = Jpsx + occbk(ib,ik)*conjg(uVpsix)*uVpsi
          Jpsy = Jpsy + occbk(ib,ik)*conjg(uVpsiy)*uVpsi
          Jpsz = Jpsz + occbk(ib,ik)*conjg(uVpsiz)*uVpsi
!          jxt = jxt + occ(ib,ik)/aLxyz*2*imag(conjg(uVpsix)*uVpsi)
!          jyt = jyt + occ(ib,ik)/aLxyz*2*imag(conjg(uVpsiy)*uVpsi)
!          jzt = jzt + occ(ib,ik)/aLxyz*2*imag(conjg(uVpsiz)*uVpsi)
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel
    Jps(1) = 2.0d0*Jpsx/vcell
    Jps(2) = 2.0d0*Jpsy/vcell
    Jps(3) = 2.0d0*Jpsz/vcell

    Je = -(real(JpkA)+imag(Jps))
    return
  End Function ubk_Je
!==
End Module Functions
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
