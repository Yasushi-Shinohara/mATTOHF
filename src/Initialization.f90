!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Initialization
  use Constants
  use Functions
  use Global_variables
  implicit none
  integer i1,ir,ig,ncount,ib,ik,ikr,ikc,ikk,ibr,ibc,irank,numkk,ii
  real(8) :: x,y,z
  real(8), allocatable :: Gkk2(:)
  complex(8), allocatable :: work1(:), work2(:)
!Fock-term renorm
  real(8) :: Rc
!FFTW staff
  include 'fftw3.f'

!System information and options
  input_data = trim(input_data)
  sys_name = trim(sys_name)
  dir_name = trim(dir_name)
  DFT_option = trim(DFT_option)
  TD_option = trim(TD_option)
  functional = trim(functional)
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(input_data,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(sys_name,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(dir_name,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(DFT_option,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(Band_option,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(TD_option,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(functional,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  If (Myrank==0) then
    Write(*,'("System infromation and options+++++++++++++++++++++++++++")') 
    Write(*,*) 'input_data: ', trim(input_data)
    Write(*,*) 'sys_name: ', trim(sys_name)
    Write(*,*) 'dir_name: ', trim(dir_name)
    Write(*,*) 'DFT_option: ', trim(DFT_option)
    Write(*,*) 'Band_option: ', Band_option
    Write(*,*) 'TD_option: ', trim(TD_option)
    Write(*,*) 'functional: ', trim(functional)
    Write(*,'("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")') 
  End If

!====
!Crystal structure
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(acell,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  ax = acell
  vcell = abs(acell)
  bx = tpi/ax
  If (Myrank==0) then
    Write(*,'("Lattice infromation +++++++++++++++++++++++++++++++++++++")') 
    write(*,*) 'acell = ', acell
    write(*,*) 'volume of cell in real-space, vcell = ',vcell
    write(*,*) 'bx = ', bx
  End If

!====
!Grid info.
  Call MPI_BCAST(NG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(NK,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  NKK = NK**2
  dvcell = vcell/dble(NG)

  numkk=0
  Do irank = 0,Nprocs-1
    NKKsarray(irank) = numkk + 1
    If (irank <= mod(NKK,Nprocs)-1) then
      numkk =  numkk + (NKK/Nprocs + 1)
    Else if (irank > mod(NKK,Nprocs)-1) then
      numkk =  numkk + NKK/Nprocs
    End if
    NKKearray(irank) = numkk
    If (Myrank == 0 ) Write(*,'("rank ID, NKKs, NKKe =",3i12)') irank, NKKsarray(irank), NKKearray(irank)
  End do
  NKKs = NKKsarray(Myrank)
  NKKe = NKKearray(Myrank)

  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  If (Myrank == 0 ) then
    Write(*,'("NG =",3i4)')  NG
    Write(*,'("Volume elment: vcell/NG =",f12.8)') dvcell
    Write(*,'("NK =",3i4)') NK
  End If
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  Allocate(rx(1:NG))
  Allocate(Gx(1:NG))
  Allocate(eiGrx(1:NG,1:NG),eiGrxc(1:NG,1:NG))
  Allocate(kx(NK))
  Allocate(vCG(1:NG))

!Grid construction
  do i1=1,NG
    rx(i1)=dble(i1-1)/dble(NG)*acell
  enddo

  Do i1=1,NG/2
    Gx(i1)=(i1-1)*bx
  End Do
  Do i1=NG/2+1,NG
    Gx(i1)=(i1-1-NG)*bx
  End Do
  If (Myrank == 0) then
    write(*,*) 'maxval(Gx**2)', maxval(Gx**2)
    Write(*,'("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")') 
  End If

  ncount=0
  do i1=1,NK
    ncount = ncount + 1
    kx(ncount)=(-0.5d0+(2.d0*i1-1.0d0)/dble(2*NK))*bx
  enddo


!Rlevant allocatoin for Hartree-Fock
  If (trim(functional)=='HF') then
    Rc = acell*dble(NK)
    Allocate(vCGkk(1:NG,NKKs:NKKe))
    Allocate(ikktable(NK,NK),ikrtable(NKK),ikctable(NKK))
    ncount=0
    Do ikr = 1,NK
      Do ikc = 1,NK
        ncount = ncount + 1
        ikktable(ikr,ikc)=ncount
        ikrtable(ncount)=ikr
        ikctable(ncount)=ikc
      End Do
    End Do
    Allocate(Gkk2(1:NG))
    Do ikk=NKKs,NKKe
      ikr = ikrtable(ikk)
      ikc = ikctable(ikk)
      Gkk2(:)=(Gx(:)+kx(ikr)-kx(ikc))**2
      Do ig = 1,NG
        vCGkk(ig,ikk) = k_vCk(Gx(ig)+kx(ikr)-kx(ikc))
      End Do
    End Do
    Deallocate(Gkk2)
    Write(*,*) 'For debug, maxval(vCGkk)=', maxval(vCGkk) !Debug
    Write(*,*) 'For debug, minval(vCGkk)=', minval(vCGkk) !Debug
!
  End If
  Do ig = 1,NG
    vCG(ig) = k_vCk(Gx(ig))
  End Do
  vCG(1) = 0.d0

!This is a check for Fourier transform of 1D-Coulomb potential
  Open(10,file='vk.out')
  Do ig = 1,NG
    Do ik = 1,NK
      Write(10,*) Gx(ig)+kx(ik), k_vCk(Gx(ig)+kx(ik))
    End Do
  End Do
  Close(10)

  Do ir=1,NG
    Do ig=1,NG/2
      eiGrx(ig,ir)=exp(zI*tpi*dble((ir-1)*(ig-1))/dble(NG))
    End Do
    Do ig=NG/2-1,NG
      eiGrx(ig,ir)=exp(zI*tpi*dble((ir-1)*(ig-1-NG))/dble(NG))
    End Do
  End Do
  eiGrxc=conjg(eiGrx)

!====
!Band occupation
  Call MPI_BCAST(NB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(NBocc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  If (Myrank == 0) then
    Write(*,'("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")') 
    Write(*,*) 'NBocc, NB =',NBocc, NB
    Write(*,'("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")') 
  End If
  Allocate(occbk(NB,NK))
  Allocate(ibbocctable(NBocc,NBocc),iboccrtable(NBocc**2),iboccctable(NBocc**2))

  Do ik = 1,NK
    Do ib = 1,NBocc
      occbk(ib,ik)=2.d0/dble(NK)
    End Do
    Do ib = NBocc+1,NB
      occbk(ib,ik)=0.d0
    End Do
  End Do

  ncount=0
  Do ibr = 1,NBocc
    Do ibc = 1,NBocc
      ncount = ncount + 1
      ibbocctable(ibr,ibc)=ncount
      iboccrtable(ncount)=ibr
      iboccctable(ncount)=ibc
    End Do
  End Do
  
!====
!Pseudo potential part
!Conversion from relative to Cartesian coordinate
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(Nspecies,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(Nion,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  If (Myrank /= 0)  Allocate(Zatom(Nspecies),alphaatom(Nspecies))
  If (Myrank /= 0)  Allocate(Zeffatom(Nspecies),alphaatom(Nspecies))
  If (Myrank /= 0)  Allocate(Rion(Nion),Kion(Nion))
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  Call MPI_BCAST(Zatom,Nspecies,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(Zeffatom,Nspecies,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(alphaatom,Nspecies,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(Rion,Nion,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(Kion,Nion,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!Zatom, Lref, Rion, Kion are already substituted by the input file.
  If (Myrank == 0) then
    Write(*,'("Ion info+++++++++++++++++++++++++++++++++++++++++++++++++")') 
    write(*,*) 'Nspecies =', Nspecies
!    write(*,*) 'Zatom(:), alphaatom(:) =',Zatom(:), alphaatom(:)
    write(*,*) 'Zeffatom(:), alphaatom(:) =',Zeffatom(:), alphaatom(:)
    write(*,*) 'Nion =', Nion
    Write(*,'("ii ,    Rion(ii),     Kion(ii)")') 
    Do ii = 1,Nion
      Write(*,'(i5,2x,f12.8,2x,i5)') ii, Rion(ii), Kion(ii)
    End Do
    Write(*,'("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")') 
  End If

  Allocate(rhoionG(1:NG))
  Allocate(vion(1:NG))
 
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(Nscf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(NCG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  If (Myrank == 0) then
    Write(*,'("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")') 
    write(*,*) 'Number of SCF iteration: ', Nscf
    write(*,*) 'Number of CG iteration within a SCF cycle: ', NCG
    Write(*,'("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")') 
  End If
  Call MPI_BCAST(NT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(Nomega,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(omegaM,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(field_option,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(dAc,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(omegac,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(Tpulse,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(phiCEP,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(nenvelope,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(E0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  If (Myrank == 0) then
    Write(*,'("Parameter for RT propagation+++++++++++++++++++++++++++++")') 
    write(*,*) 'NT =' ,NT
    write(*,*) 'dt =', dt, '[a.u.]'
    write(*,*) 'Nomega =', Nomega
    write(*,*) 'omegaM =', omegaM, '[a.u.]'
    write(*,*) 'field_option =',field_option
    write(*,*) 'dAc =', dAc, '[a.u.]'
    write(*,*) 'omegac =',omegac, '[eV]'
    write(*,*) 'Tpulse =', Tpulse, '[fs]'
    write(*,*) 'phiCEP =',phiCEP
    write(*,*) 'nenvelope =', nenvelope
    write(*,*) 'E0 = ',E0, '[V/nm]'
    Write(*,'("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")') 
  End If
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!plan preparation for FFTW
  Allocate(work1(NG),work2(NG))
  Call dfftw_plan_dft_1d(planf, NG, work1(:), work2(:), FFTW_FORWARD , FFTW_ESTIMATE)
  Call dfftw_plan_dft_1d(planb, NG, work2(:), work1(:), FFTW_BACKWARD, FFTW_ESTIMATE)
  Deallocate(work1,work2)

  return
End Subroutine Initialization
