program test_decoherence
  use types
  use constants
  use read
  use crystal
  use interactions
  use deco
  implicit none

  double precision :: t

  double precision, allocatable :: pseudo_angle(:),eigen_ener(:)
  type (rot), allocatable :: matrot(:)
  type (rot), allocatable :: matrottrans(:)
  

  ! test rotmat
  nb_pairs = 1

  allocate (pseudo_angle(nb_pairs))
  allocate (matrot(nb_pairs))
  allocate (matrottrans(nb_pairs))

  pseudo_angle(1) = pi
  
  call rotmat (pseudo_angle, matrot, matrottrans)

  write(*,*)'matrot',matrot(1)%elements
  write(*,*)'matrottrans',matrottrans(1)%elements

  write(*,*)matmul(matrot(1)%elements, matrottrans(1)%elements)
 
  ! test Z_gate
  allocate (eigen_ener(nb_pairs),Zgate(nb_pairs))
  eigen_ener(1) = 1.d0
  t = 2.d0
  call Z_gate (eigen_ener, t)

  write(*,*)'Zgate',Zgate(1)%elements

  ! Test pseudo
  nb_pairs = 1
  deallocate (Zgate)
  allocate (C12(nb_pairs),DJ(nb_pairs))

  C12(1) = 1.d0
  DJ(1) = 1.d0
  polar = 1.d0

  Tmax = 10.d0
  nb_pts_t = 10

  call pseudo
  

end program test_decoherence

