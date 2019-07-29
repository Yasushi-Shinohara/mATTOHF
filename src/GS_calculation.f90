!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine GS_calculation(ubk)
  use Constants
  use Global_variables, only : sys_name, NG,Gx,NB, NBocc, NK, kx, Nscf, Myrank,dvcell, MPI_COMM_WORLD,ierr, MPI_WTIME
  use Functions
!  use Functions_MPI
  implicit none
  complex(8),intent(inout) :: ubk(NG,NB,NK)
  integer, parameter :: NDoS=10000
  character(256) :: GSinfo, DoSfile, bandfile
  real(8) :: beta !Inverse temperature
  real(8) :: dns(NG), dnsscf(NG,Nscf+1), dnsmix = 0.5d0
  real(8) :: epsbk(NB,NK), epserrorbk(NB,NK)
  real(8) :: Etot, A0=0.d0, Etotscf(0:Nscf), Ecomps(4), JeGS
  real(8) :: epsDoS(NDoS), DoS(NDoS), DoSocc(NDoS), epsmin, epsmax
  integer :: iter, iDoS, ib,ik
  real(8) :: tGSs, tGSe

  tGSs = MPI_WTIME()
  If (Myrank == 0) then
    GSinfo = trim(sys_name)//'_GSinfo.out'
    Open(10,file=trim(GSinfo),action='write')
    Write(10,'("#iter, Etotscf(iter), Etoscf(iter)-Etotscf(iter-1), ,sum(abs(dnsscf(:,iter+1)-dnsscf(:,iter)))/dble(NG) &
         sum(epserrorbk)/dble(NB*NK),maxval(epserrorbk),Fundamental gap, Optical Gap")')
  End If
  dnsscf(:,1) = ubk_dns(ubk)
  Etotscf(0) = ubk_Etot(ubk,A0)
  Do iter = 1,Nscf
    Call Gram_Schmidt(ubk)
!    Call KS_ubk2_SVD_MPI(ubk,A0)
    Call KS_ubk2_SVD(ubk,A0)
    Call Gram_Schmidt(ubk)
    Call Subspace_diag(ubk,A0)

!    Etotscf(iter) = ubk_Etot_MPI(ubk,A0)
    Etotscf(iter) = ubk_Etot(ubk,A0)
!    JeGS = ubk_Je_MPI(ubk,A0)
    JeGS = ubk_Je(ubk,A0)

!    epsbk = ubk_epsbk_MPI(ubk)
!    epserrorbk = ubk_epserrorbk_MPI(ubk)
    epsbk = ubk_epsbk(ubk)
    epserrorbk = ubk_epserrorbk(ubk)
    dnsscf(:,iter+1) = dnsmix*ubk_dns(ubk)+(1.d0-dnsmix)*dnsscf(:,iter)
    If (Myrank == 0) then
      tGSe = MPI_WTIME()
      Write(*,'("iter, Etotscf(iter), JeGS =",i4,f14.8,e12.4,2x,3e12.4)') iter, Etotscf(iter), Etotscf(iter)-Etotscf(iter-1), JeGS
      Write(*,'(4(i3,f12.6,2x))') (ib,epsbk(ib,1),ib=1,NB)
      Write(10,'(2x,i5,100e20.12)') iter, Etotscf(iter), Etotscf(iter)-Etotscf(iter-1)  &
           ,sum(abs(dnsscf(:,iter+1)-dnsscf(:,iter)))/dble(NG) &
           ,sum(epserrorbk)/dble(NB*NK),maxval(epserrorbk) &
           ,minval(epsbk(NBocc+1:,:)) - maxval(epsbk(:NBocc,:)),minval(epsbk(NBocc+1,:)-epsbk(NBocc,:))
      If (mod(iter,20) == 0) &
           Write(*,'("##GS calc-time for ",i4," iterations :",e12.4," [sec] =",f12.8," [min]##")') &
           iter,tGSe-tGSs, (tGSe-tGSs)/60.d0
    End If
  End Do
  If (Myrank == 0) Write(*,'("########SCF procedure is finished.########")')

  Ecomps = ubk_Ecomps(ubk,A0)

  If (Myrank == 0) then
    Write(10,'("#Ekin  at SCF:",f12.8)') Ecomps(1)
    Write(10,'("#Eion  at SCF:",f12.8)') Ecomps(2)
    Write(10,'("#EH    at SCF:",f12.8)') Ecomps(3)
    Write(10,'("#Exc   at SCF:",f12.8)') Ecomps(4)
    Write(*,'(" Min. value of the spectra: minval(epsbk)    =",f12.8," =",f12.8," [eV]")') &
         minval(epsbk), minval(epsbk)*Hartree
    Write(*,'(" Max. value of the spectra: maxval(epsbk)    =",f12.8," =",f12.8," [eV]")') &
         maxval(epsbk), maxval(epsbk)*Hartree
    Write(*,'(" Valence top: maxval(epsbk(:NBocc,:))        =",f12.8," =",f12.8," [eV]")') &
         maxval(epsbk(:NBocc,:)), maxval(epsbk(:NBocc,:))*Hartree
    Write(*,'(" Conduction bottom: minval(epsbk(NBocc+1:,:))=",f12.8," =",f12.8," [eV]")') &
         maxval(epsbk(NBocc+1:,:)), maxval(epsbk(NBocc+1:,:))*Hartree
    Write(*,'(" Fundamental gap: ",f12.8," =",f12.8," [eV]")') &
         minval(epsbk(NBocc+1:,:)) - maxval(epsbk(:NBocc,:)), (minval(epsbk(NBocc+1:,:)) - maxval(epsbk(:NBocc,:)))*Hartree
    Write(*,'(" Min. gap among same k-point: ",f12.8," =",f12.8," [eV]")') &
         minval(epsbk(NBocc+1,:)-epsbk(NBocc,:)), (minval(epsbk(NBocc+1,:)-epsbk(NBocc,:)))*Hartree
    Write(10,'("#Min. value of the spectra: minval(epsbk)    =",f12.8," =",f12.8," [eV]")') &
         minval(epsbk), minval(epsbk)*Hartree
    Write(10,'("#Max. value of the spectra: maxval(epsbk)    =",f12.8," =",f12.8," [eV]")') &
         maxval(epsbk), maxval(epsbk)*Hartree
    Write(10,'("#Valence top: maxval(epsbk(:NBocc,:))        =",f12.8," =",f12.8," [eV]")') &
         maxval(epsbk(:NBocc,:)), maxval(epsbk(:NBocc,:))*Hartree
    Write(10,'("#Conduction bottom: minval(epsbk(NBocc+1:,:))=",f12.8," =",f12.8," [eV]")') &
         maxval(epsbk(NBocc+1:,:)), maxval(epsbk(NBocc+1:,:))*Hartree
    Write(10,'("#Fundamental gap: ",f12.8," =",f12.8," [eV]")') &
         minval(epsbk(NBocc+1:,:)) - maxval(epsbk(:NBocc,:)), (minval(epsbk(NBocc+1:,:)) - maxval(epsbk(:NBocc,:)))*Hartree
    Write(10,'("#Min. gap among same k-point: ",f12.8," =",f12.8," [eV]")') &
         minval(epsbk(NBocc+1,:)-epsbk(NBocc,:)), (minval(epsbk(NBocc+1,:)-epsbk(NBocc,:)))*Hartree
    Close(10)
  End If
  epsmin = minval(epsbk) - 0.1d0*(maxval(epsbk) - minval(epsbk)); epsmax = maxval(epsbk) + 0.1d0*(maxval(epsbk) - minval(epsbk))
  Do iDoS = 1, NDoS
     epsDoS(iDoS) = epsmin + (epsmax-epsmin)*dble(iDoS-1)/dble(NDoS)
  End Do
  Call Make_DoS(epsbk,NDoS,epsDoS,DoS,DoSocc)

  DoSfile = trim(sys_name)//'_DoS.out'
  bandfile = trim(sys_name)//'_band.out'
  If (Myrank == 0) then
    Open(10,file=trim(DoSfile),action='write')
    Write(10,*) '# ik,kx(ik),  sum(epserrorbk(:,ik))/dble(NB),  epsbk(:,ik)'
    Do ik = 1,NK
      Write(10,'("# ",i5,3e20.12,4x,1e20.12,4x,1000e20.12)') ik,kx(ik),sum(epserrorbk(:,ik))/dble(NB),epsbk(:,ik)
    End Do
    Do iDoS = 1, NDoS
       write(10,'(3e20.12e3)') epsDoS(iDoS), DoS(iDoS), DoSocc(iDoS)
    End Do
    Close(10)
    Write(*,*) 'DoS check =', sum(DoS)*(epsDoS(2)-epsDoS(1)), sum(DoSocc)*(epsDoS(2)-epsDoS(1))
    Open(10,file=trim(bandfile),action='write')
    Write(10,*) '# ik,kx(ik),  sum(epserrorbk(:,ik))/dble(NB),  epsbk(:,ik)'
    Do ik = 1,NK
      Write(10,'(i5,3e20.12,4x,1e20.12,4x,1000e20.12)') ik,kx(ik),sum(epserrorbk(:,ik))/dble(NB),epsbk(:,ik)
    End Do
    Close(10)
  End If

  return
End Subroutine GS_calculation
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Diagonarization of orbitals based on Krylov-subspace
Subroutine KS_ubk2_SVD(ubk,A0)
  use Constants
  use Global_variables, only : NG, NB, NK, NBocc, NCG, dvcell, Gx &
         , MPI_DOUBLE_COMPLEX, MPI_SUM, Myrank, MPI_COMM_WORLD,ierr
  use Functions
!$  use omp_lib
  implicit none
  complex(8),intent(inout) :: ubk(NG,NB,NK)
  real(8),intent(in) :: A0
  complex(8) :: work(NG,NB*NCG,NK), coef(NB*NCG,NB*NCG) ,hwork(NG,NB*NCG,NK), work2(NG,NB*NCG)
  real(8) :: eigs(NB*NCG)
  complex(8) :: whw(NB*NCG,NB*NCG)
  integer :: ik, iCG, ibr, ibc, iter
  real(8) :: s
  complex(8) :: zov
  real(8) :: eshift
!For ZHEEV
  integer :: LWORK_EV
  complex(8) :: WORK_EV(2*(2*NB*NCG-1))
  real(8) :: RWORK_EV(3*NB*NCG-2)
  integer :: INFO_EV
!For ZGESVD
  integer :: INFO_SVD, LWORK_SVD
  Real(8) :: S_SVD(NB*NCG), RWORK_SVD(5*Min(NG,NB*NCG))
  Complex(8) :: WORK_SVD(10*(2*Min(NG,NB*NCG)+Max(NG,NB*NCG))),U_SVD(1),VT_SVD(1)
  LWORK_EV = 2*(2*NB*NCG-1)
  LWORK_SVD =10*(2*Min(NG,NB*NCG)+Max(NG,NB*NCG))

  eshift = 0.d0*eshift

  work(:,1:NB,:) = ubk(:,1:NB,:)
  Do iCG = 2,NCG
    work(:,NB*(iCG-1)+1:NB*iCG,:) = ubkwbk_hwbk(ubk,work(:,NB*(iCG-2)+1:NB*(iCG-1),:),A0) &
         -eshift*work(:,NB*(iCG-2)+1:NB*(iCG-1),:)
  End Do

!Multiple Gram-Schmidt procedures to get accurate orthogonality
  Do ik = 1,NK
    Call ZGESVD('O','N',NG,NB*NCG,work(:,:,ik),NG,S_SVD,U_SVD,1,VT_SVD,1,&
         WORK_SVD,LWORK_SVD,RWORK_SVD,INFO_SVD)
    s = sum(abs(work(:,1,ik))**2)*dvcell
    work(:,:,ik) = work(:,:,ik)/sqrt(s)
  End Do

!Operation of the Hamiltonian for making matrix elmenet
  Do iCG = 1,NCG
    hwork(:,NB*(iCG-1)+1:NB*iCG,:) = ubkwbk_hwbk(ubk,work(:,NB*(iCG-1)+1:NB*iCG,:),A0)
  End Do

  Do ik = 1,NK
!Constructing matrix within the subspace     
    Do ibr = 1,NB*NCG
      Do ibc = 1,ibr-1
        whw(ibr,ibc) = sum(conjg(work(:,ibr,ik))*hwork(:,ibc,ik))*dvcell
        whw(ibc,ibr) = conjg(whw(ibr,ibc))
      End Do
      whw(ibr,ibr) = real(sum(conjg(work(:,ibr,ik))*hwork(:,ibr,ik)))*dvcell
    End Do
    coef = whw
    Call ZHEEV('V','U',NB*NCG,coef,NB*NCG,eigs,WORK_EV,LWORK_EV,RWORK_EV,INFO_EV)
    work2 = 0.d0
    Do ibr = 1,NB*NCG
      Do ibc = 1,NB*NCG
        work2(:,ibr) = work2(:,ibr) + work(:,ibc,ik)*coef(ibc,ibr)
      End Do
    End Do
    ubk(:,1:NB,ik) = work2(:,1:NB)
  End Do

  return
End Subroutine KS_ubk2_SVD
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Imaginary time evolution with given dns
Subroutine Subspace_diag(ubk,A0)
  use Constants
  use Global_variables, only : NG, NB, NK, dvcell
  use Functions
  implicit none
  complex(8),intent(inout) :: ubk(NG,NB,NK)
  real(8),intent(in) :: A0
  complex(8) :: workbk(NG,NB,NK), coef(NB,NB), workbk2(NG,NB,NK)
  real(8) :: eigs(NB)
  complex(8) :: ubhub(NB,NB)
  integer :: ibr,ibc,ik
!For ZHEEV
  integer :: LWORK_EV
  complex(8) :: WORK_EV(2*(2*NB-1))
  real(8) :: RWORK_EV(3*NB-2)
  integer :: INFO_EV
  LWORK_EV = 2*(2*NB-1)

  workbk = ubkwbk_hwbk(ubk,ubk,A0)
  workbk2 = 0.d0
  Do ik = 1,NK
    Do ibr = 1,NB
      Do ibc = 1,ibr-1
        ubhub(ibr,ibc) = sum(conjg(ubk(:,ibr,ik))*workbk(:,ibc,ik))*dvcell
        ubhub(ibc,ibr) = conjg(ubhub(ibr,ibc))
      End Do
      ubhub(ibr,ibr) = real(sum(conjg(ubk(:,ibr,ik))*workbk(:,ibr,ik)))*dvcell
    End Do
    coef = ubhub
    Call ZHEEV('V','U',NB,coef,NB,eigs,WORK_EV,LWORK_EV,RWORK_EV,INFO_EV)
    Do ibr = 1,NB
      Do ibc = 1,NB
        workbk2(:,ibr,ik) = workbk2(:,ibr,ik) + ubk(:,ibc,ik)*coef(ibc,ibr)
      End Do
    End Do
  End Do
  ubk = workbk2

      
  return
End Subroutine Subspace_diag
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Make_DoS(epsbk,NDoS,epsDoS,DoS,DoSocc)
  use Constants
  use Global_variables, only : NB, NK, occbk
  use Functions
  implicit none
  integer,intent(in) :: NDoS
  real(8),intent(in) :: epsbk(NB,NK), epsDoS(NDoS)
  real(8),intent(out) :: DoS(NDoS), DoSocc(NDoS)
  real(8) :: ewidth, val(NDoS), occmax
  integer :: ib,ik,iDoS

  ewidth = 5.d0*(epsDoS(2) - epsDoS(1))
  DoS = 0.d0; DoSocc = 0.d0
  Do ik = 1,NK
    occmax = maxval(occbk(:,ik))
    Do ib = 1,NB
      val(:) = exp(-(epsbk(ib,ik)-epsDoS(:))**2/ewidth**2)
      DoS(:) = DoS(:) + occmax*val(:)
      DoSocc(:) = DoSocc(:) + occbk(ib,ik)*val(:)
    End Do
  End Do
!Normalization: delta(x) = exp(-(x-x0)/e**2)/(e*sqrt(pi))
   DoS(:) = DoS(:)/(ewidth*sqrt(pi))
   DoSocc(:) = DoSocc(:)/(ewidth*sqrt(pi))
      
  return
End Subroutine Make_DoS
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
