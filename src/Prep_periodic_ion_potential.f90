!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Prep_periodic_ion_potential
  use Constants
  use Global_Variables, only : vion, NK, NG, Nion, Zeffatom, alphaatom, Rion, Kion, rx, ax
  implicit none
  integer :: icell,ii
  real(8) :: Za, dist(1:NG), pottop
  integer :: ig

  vion = 0.d0
  Do ii = 1, Nion
!    Za = dble(Zatom(Kion(ii)))
    Za = dble(Zeffatom(Kion(ii)))
    Do icell = -4*NK, 4*NK
      dist(:) = abs(rx(:)-ax*Rion(ii)-ax*icell)
      vion(:) = vion(:) - Za/sqrt(dist(:)**2+alphaatom(Kion(ii))**2)
    End Do
  End Do
  pottop = maxval(vion)
  vion(:) = vion(:) - pottop
  
!Debug
  Open(10,file='potential.out')
  Write(*,*) NG
  Do ig = 1,NG
    Write(10,*) rx(ig), vion(ig)
  End Do
!Debug

  return
End Subroutine Prep_periodic_ion_potential
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
