program test_decoherence
  use types
  implicit none
 
  ! test eigen_energies
  call test_eigen_energies

  ! test pseudo_angles
  call test_pseudo_angles

  ! test rotmat
  call test_rotmat

  ! test Z_gate
  call test_Z_gate

  ! test FID
  call test_FID

  ! test Hahn
  call test_Hahn

end program test_decoherence

subroutine test_Hahn
  use types
  use deco
  implicit none

  integer :: i,k
  double precision :: t, t2, dt, L
  double precision, allocatable :: eigen_ener_up(:)
  double precision, allocatable :: eigen_ener_down(:)
  double precision, allocatable :: pseudo_angle_up(:)
  double precision, allocatable :: pseudo_angle_down(:)
  double precision, allocatable :: L_pairs(:)
  double precision, allocatable :: A_0(:), Ax(:), Ay(:), Az(:), A_02(:)
  double precision, allocatable :: Axu2(:), Axl2(:), Azu2(:), Azl2(:)
  double precision, allocatable :: theta_uf(:), theta_lf(:), eps(:)
  character(len=100) :: filename

  nb_pairs = 5

  allocate (eigen_ener_up(nb_pairs), eigen_ener_down(nb_pairs))
  allocate (pseudo_angle_up(nb_pairs), pseudo_angle_down(nb_pairs))
  allocate (A_0(nb_pairs), Ax(nb_pairs), A_02(nb_pairs))
  allocate (Ay(nb_pairs), Az(nb_pairs))
  allocate (L_pairs(nb_pairs))
  allocate (Axu2(nb_pairs), Axl2(nb_pairs), Azu2(nb_pairs), Azl2(nb_pairs))
  allocate (theta_uf(nb_pairs), theta_lf(nb_pairs), eps(nb_pairs))

  eigen_ener_up = (/(i, i=1,nb_pairs)/)
  eigen_ener_down = (/(i/2.d0, i=1,nb_pairs)/)
  pseudo_angle_up = (/(i**2, i=1,nb_pairs)/)
  pseudo_angle_down = (/(i**3, i=1,nb_pairs)/)

  Qubittype = "electron"
  Dynadeco  = "Hahn"
  CP_seq    = 1
  Tmax = 10.d0
  nb_pts_t = 1000

  call Hahn (eigen_ener_up, eigen_ener_down, pseudo_angle_up, &
                   pseudo_angle_down, CP_seq)

  write(filename, '(a)') "Averaged_decay_electron_IFF_Hahn_test.dat"
  open(16, file=filename)

  dt = Tmax / dble(nb_pts_t)
  t  = - dt

  do k=1,nb_pts_t + 1
     t = t + dt
     t2 = t / 2.d0
     A_0 = cos(eigen_ener_up * t2) * cos(eigen_ener_down * t2)
     A_0 = A_0 - sin(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          cos(pseudo_angle_up - pseudo_angle_down)
     Ax = cos(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          sin(pseudo_angle_down) 
     Ax = Ax + cos(eigen_ener_down * t2) * sin(eigen_ener_up * t2) * &
          sin(pseudo_angle_up)
     Ay = - sin(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          sin(pseudo_angle_up - pseudo_angle_down)
     Az = cos(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          cos(pseudo_angle_down) 
     Az = Az + cos(eigen_ener_down * t2) * sin(eigen_ener_up * t2) * &
          cos(pseudo_angle_up)

     L_pairs = abs(1.d0 - 2.d0 * Ay**2 + 2.d0 * dcmplx(0.d0, Ax * Ay))
     L_pairs = 0.5d0 + 0.5d0 * L_pairs
     
     L = product(L_pairs)
     write(16, "(es20.10e3, es20.10e3)")t, L
  end do
  close(16)
     
  write(*,*) "Both Averaged_decay_electron_IFF_Hahn.dat and Averaged_decay_electron_IFF_Hahn_test.dat are ready to be compared..."


  Qubittype = "electron"
  Dynadeco  = "CP"
  CP_seq    = 2

  call Hahn (eigen_ener_up, eigen_ener_down, pseudo_angle_up, &
                   pseudo_angle_down, CP_seq)

  write(filename, '(a)') "Averaged_decay_electron_IFF_CP2_test.dat"
  open(16, file=filename)

  t  = - dt
  do k=1,nb_pts_t + 1
     t = t + dt
     t2 = t / 4.d0
     A_0 = cos(eigen_ener_up * t2) * cos(eigen_ener_down * t2)
     A_0 = A_0 - sin(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          cos(pseudo_angle_up - pseudo_angle_down)
     Ax = cos(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          sin(pseudo_angle_down) 
     Ax = Ax + cos(eigen_ener_down * t2) * sin(eigen_ener_up * t2) * &
          sin(pseudo_angle_up)
     Ay = - sin(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          sin(pseudo_angle_up - pseudo_angle_down)
     Az = cos(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          cos(pseudo_angle_down) 
     Az = Az + cos(eigen_ener_down * t2) * sin(eigen_ener_up * t2) * &
          cos(pseudo_angle_up)

     L_pairs = abs(1.d0 - 8.d0 * Ay**2 * (Ax**2 + Az**2) &
          + 4.d0 * dcmplx(0.d0, 1.d0) * (1.d0 - 2.d0 * (Ax**2 + Az**2)) * &
          Ax * Ay)
     L_pairs = 0.5d0 + 0.5d0 * L_pairs
     
     L = product(L_pairs)
     write(16, "(es20.10e3, es20.10e3)")t, L
  end do
  close(16)
     
  write(*,*) "Both Averaged_decay_electron_IFF_CP2.dat and Averaged_decay_electron_IFF_CP2_test.dat are ready to be compared..."

  Qubittype = "electron"
  Dynadeco  = "CP"
  CP_seq    = 6

  call Hahn (eigen_ener_up, eigen_ener_down, pseudo_angle_up, &
                   pseudo_angle_down, CP_seq)

  write(filename, '(a)') "Averaged_decay_electron_IFF_CP6_test.dat"
  open(16, file=filename)

  t  = - dt
  do k=1,nb_pts_t + 1
     t = t + dt
     t2 = t / 12.d0
     A_0 = cos(eigen_ener_up * t2) * cos(eigen_ener_down * t2)
     A_0 = A_0 - sin(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          cos(pseudo_angle_up - pseudo_angle_down)
     Ax = cos(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          sin(pseudo_angle_down) 
     Ax = Ax + cos(eigen_ener_down * t2) * sin(eigen_ener_up * t2) * &
          sin(pseudo_angle_up)
     Ay = - sin(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          sin(pseudo_angle_up - pseudo_angle_down)
     Az = cos(eigen_ener_up * t2) * sin(eigen_ener_down * t2) * &
          cos(pseudo_angle_down) 
     Az = Az + cos(eigen_ener_down * t2) * sin(eigen_ener_up * t2) * &
          cos(pseudo_angle_up)

     A_02 = 1.d0 - 2.d0 * (Ax**2 + Az**2)
     Axu2 = 2.d0 * (Ax * A_0 - Ay * Az)
     Axl2 = 2.d0 * (Ax * A_0 + Ay * Az)
     Azu2 = 2.d0 * (Az * A_0 + Ax * Ay)
     Azl2 = 2.d0 * (Az * A_0 - Ax * Ay)

     eps = acos(A_02)

     theta_uf = atan2(Axu2 , Azu2)
     theta_lf = atan2(Axl2 , Azl2)
     

     L_pairs = 1.d0 - sin((theta_uf - theta_lf)/2.d0)**2 * &
          (1.d0 - cos(6.d0 * eps))
     L_pairs = abs(L_pairs + dcmplx(0.d0,1.d0) * sin(6.d0 * eps) * &
          sin((theta_uf - theta_lf)/2.d0) * sin((theta_uf + theta_lf)/2.d0))

     L_pairs = 0.5d0 + 0.5d0 * L_pairs
     
     L = product(L_pairs)
     write(16, "(es20.10e3, es20.10e3)")t, L
  end do
  close(16)
     
  write(*,*) "Both Averaged_decay_electron_IFF_CP6.dat and Averaged_decay_electron_IFF_CP6_test.dat are ready to be compared..."
  
  deallocate (eigen_ener_up, eigen_ener_down)
  deallocate (pseudo_angle_up, pseudo_angle_down)
  deallocate (A_0, Ax)
  deallocate (Ay, Az)
  deallocate (L_pairs)
  deallocate (Axu2, Axl2, Azu2, Azl2)
  deallocate (theta_uf, theta_lf, eps)

end subroutine test_Hahn

subroutine test_FID
  use types
  use deco
  implicit none
  integer :: i,k
  double precision :: t, dt, L
  double precision, allocatable :: eigen_ener_up(:)
  double precision, allocatable :: eigen_ener_down(:)
  double precision, allocatable :: pseudo_angle_up(:)
  double precision, allocatable :: pseudo_angle_down(:)
  double precision, allocatable :: L_pairs(:)
  type (rot), allocatable :: matrot_u(:), matrottrans_u(:)
  type (rot), allocatable :: matrot_d(:), matrottrans_d(:)
  type (rot), allocatable :: Zgate_u(:), Zgate_d(:)
  type (rot), allocatable :: Tu(:), Td(:)
  character(len=100) :: filename

  nb_pairs = 5

  allocate (eigen_ener_up(nb_pairs), eigen_ener_down(nb_pairs))
  allocate (pseudo_angle_up(nb_pairs), pseudo_angle_down(nb_pairs))
  allocate (matrot_u(nb_pairs), matrottrans_u(nb_pairs))
  allocate (matrot_d(nb_pairs), matrottrans_d(nb_pairs))
  allocate (Zgate_u(nb_pairs), Zgate_d(nb_pairs))
  allocate (Tu(nb_pairs), Td(nb_pairs), L_pairs(nb_pairs))

  eigen_ener_up = (/(i, i=1,nb_pairs)/)
  eigen_ener_down = (/(i/2.d0, i=1,nb_pairs)/)
  pseudo_angle_up = (/(i**2, i=1,nb_pairs)/)
  pseudo_angle_down = (/(i**3, i=1,nb_pairs)/)
  
  matrot_u%elements(1,1)=cos(pseudo_angle_up/2.d0)
  matrot_u%elements(1,2)=sin(pseudo_angle_up/2.d0)  
  matrot_u%elements(2,1)=-sin(pseudo_angle_up/2.d0)
  matrot_u%elements(2,2)=cos(pseudo_angle_up/2.d0)

  matrottrans_u%elements(1,1)=cos(pseudo_angle_up/2.d0)
  matrottrans_u%elements(1,2)=-sin(pseudo_angle_up/2.d0)  
  matrottrans_u%elements(2,1)=sin(pseudo_angle_up/2.d0)
  matrottrans_u%elements(2,2)=cos(pseudo_angle_up/2.d0)

  matrot_d%elements(1,1)=cos(pseudo_angle_down/2.d0)
  matrot_d%elements(1,2)=sin(pseudo_angle_down/2.d0)  
  matrot_d%elements(2,1)=-sin(pseudo_angle_down/2.d0)
  matrot_d%elements(2,2)=cos(pseudo_angle_down/2.d0)

  matrottrans_d%elements(1,1)=cos(pseudo_angle_down/2.d0)
  matrottrans_d%elements(1,2)=-sin(pseudo_angle_down/2.d0)  
  matrottrans_d%elements(2,1)=sin(pseudo_angle_down/2.d0)
  matrottrans_d%elements(2,2)=cos(pseudo_angle_down/2.d0)

  Qubittype = "electron"

  Tmax = 10.d0
  nb_pts_t = 1000

  call FID (eigen_ener_up, eigen_ener_down, pseudo_angle_up, &
                  pseudo_angle_down)

  write(filename, '(a)') "Averaged_decay_electron_IFF_FID_test.dat"
  open(16, file=filename)

  dt = Tmax / dble(nb_pts_t)
  t  = - dt

  do k=1,nb_pts_t + 1
     t = t + dt

     do i=1,nb_pairs

        Zgate_u%elements(1,1)=exp(-dcmplx(0.d0,eigen_ener_up)*t)
        Zgate_u%elements(1,2)=0.d0
        Zgate_u%elements(2,1)=0.d0
        Zgate_u%elements(2,2)=exp(dcmplx(0.d0,eigen_ener_up)*t)
        
        Zgate_d%elements(1,1)=exp(-dcmplx(0.d0,eigen_ener_down)*t)
        Zgate_d%elements(1,2)=0.d0
        Zgate_d%elements(2,1)=0.d0
        Zgate_d%elements(2,2)=exp(dcmplx(0.d0,eigen_ener_down)*t)

        Tu(i)%elements = matmul(Zgate_u(i)%elements, matrot_u(i)%elements)
        Tu(i)%elements = matmul(matrottrans_u(i)%elements, Tu(i)%elements)

        Td(i)%elements = matmul(Zgate_d(i)%elements, matrot_d(i)%elements)
        Td(i)%elements = matmul(matrottrans_d(i)%elements, Td(i)%elements)

        L_pairs(i) = abs(conjg(Td(i)%elements(1, 1))*Tu(i)%elements(1, 1) &
               - Td(i)%elements(1, 2) * Tu(i)%elements(2, 1))

        
        L_pairs(i) = 0.5d0 + 0.5d0*L_pairs(i)

     end do
     L = product(L_pairs) 
     write(16, "(es20.10e3, es20.10e3)")t, L

  end do
  close(16)

  write(*,*) "Both Averaged_decay_electron_IFF_FID.dat and Averaged_decay_electron_IFF_FID_test.dat are ready to be compared..."

  deallocate (eigen_ener_up, eigen_ener_down)
  deallocate (pseudo_angle_up, pseudo_angle_down)
  deallocate (matrot_u, matrottrans_u)
  deallocate (matrot_d, matrottrans_d)
  deallocate (Zgate_u, Zgate_d)
  deallocate (Tu, Td, L_pairs)

end subroutine test_FID

subroutine test_Z_gate
  use types
  use deco
  implicit none 
  integer :: i 
  double precision :: t
  double precision, allocatable :: eigen_ener_up(:)
  double precision, allocatable :: eigen_ener_down(:)
  type (rot), allocatable :: Zgate_u(:), Zgate_d(:)
  logical :: test

  nb_pairs = 5

  allocate (eigen_ener_up(nb_pairs), eigen_ener_down(nb_pairs))
  allocate (Zgate_u(nb_pairs), Zgate_d(nb_pairs))

  eigen_ener_up = (/(i, i=1,nb_pairs)/)
  eigen_ener_down = (/(i**2, i=1,nb_pairs)/)

  test = .true.
  do i=1,nb_pairs
     t = 10.d0
     call Z_gate (eigen_ener_up, t, Zgate_u)
     call Z_gate (eigen_ener_down, t, Zgate_d)

     if (Zgate_u(i)%elements(1,1) .ne. exp(-dcmplx(0.d0,dble(i)*t))) &
          test = .false.
     if (Zgate_u(i)%elements(1,2) .ne. dcmplx(0.d0,0.d0)) test = .false.
     if (Zgate_u(i)%elements(2,1) .ne. dcmplx(0.d0,0.d0)) test = .false.
     if (Zgate_u(i)%elements(2,2) .ne. exp(dcmplx(0.d0,dble(i)*t))) &
          test = .false.

     if (Zgate_d(i)%elements(1,1) .ne. exp(-dcmplx(0.d0,dble(i**2)*t))) &
          test = .false.
     if (Zgate_d(i)%elements(1,2) .ne. dcmplx(0.d0,0.d0)) test = .false.
     if (Zgate_d(i)%elements(2,1) .ne. dcmplx(0.d0,0.d0)) test = .false.
     if (Zgate_d(i)%elements(2,2) .ne. exp(dcmplx(0.d0,dble(i**2)*t))) &
          test = .false.

  end do

  if (test) write(*,*) "test Z_gate ok..."
  if (test == .false.) write(*,*) "Pb in Z_gate"

  deallocate (eigen_ener_up, eigen_ener_down)
  deallocate (Zgate_u, Zgate_d)
 
end subroutine test_Z_gate

subroutine test_rotmat
  use types
  use deco
  implicit none
  integer :: i 
  double precision, allocatable :: pseudo_angle_up(:)
  double precision, allocatable :: pseudo_angle_down(:)
  type (rot), allocatable :: matrot_u(:), matrottrans_u(:)
  type (rot), allocatable :: matrot_d(:), matrottrans_d(:)
  logical :: test

  nb_pairs = 5

  allocate (pseudo_angle_up(nb_pairs), pseudo_angle_down(nb_pairs))
  allocate (matrot_u(nb_pairs), matrottrans_u(nb_pairs))
  allocate (matrot_d(nb_pairs), matrottrans_d(nb_pairs))

  pseudo_angle_up = (/(i, i=1,nb_pairs)/)
  pseudo_angle_down  = (/(i**2, i=1,nb_pairs)/)

  call rotmat (pseudo_angle_up, matrot_u, matrottrans_u)
  call rotmat (pseudo_angle_down, matrot_d, matrottrans_d)

  test = .true.
  do i=1,nb_pairs
     if (matrot_u(i)%elements(1,1) .ne. cos(i * 0.5d0)) test = .false.
     if (matrot_u(i)%elements(1,2) .ne. sin(i * 0.5d0)) test = .false.
     if (matrot_u(i)%elements(2,1) .ne. -sin(i * 0.5d0)) test = .false.
     if (matrot_u(i)%elements(2,2) .ne. cos(i * 0.5d0)) test = .false.

     if (matrottrans_u(i)%elements(1,1) .ne. cos(i * 0.5d0)) test = .false.
     if (matrottrans_u(i)%elements(1,2) .ne. -sin(i * 0.5d0)) test = .false.
     if (matrottrans_u(i)%elements(2,1) .ne. sin(i * 0.5d0)) test = .false.
     if (matrottrans_u(i)%elements(2,2) .ne. cos(i * 0.5d0)) test = .false.

     if (matrot_d(i)%elements(1,1) .ne. cos(i**2 * 0.5d0)) test = .false.
     if (matrot_d(i)%elements(1,2) .ne. sin(i**2 * 0.5d0)) test = .false.
     if (matrot_d(i)%elements(2,1) .ne. -sin(i**2 * 0.5d0)) test = .false.
     if (matrot_d(i)%elements(2,2) .ne. cos(i**2 * 0.5d0)) test = .false.

     if (matrottrans_d(i)%elements(1,1) .ne. cos(i**2 * 0.5d0)) test = .false.
     if (matrottrans_d(i)%elements(1,2) .ne. -sin(i**2 * 0.5d0)) test = .false.
     if (matrottrans_d(i)%elements(2,1) .ne. sin(i**2 * 0.5d0)) test = .false.
     if (matrottrans_d(i)%elements(2,2) .ne. cos(i**2 * 0.5d0)) test = .false.
  end do

  if (test) write(*,*) "test rotmat ok..."
  if (test == .false.) write(*,*) "Pb in rotmat"

  deallocate (pseudo_angle_up, pseudo_angle_down)
  deallocate (matrot_u, matrottrans_u)
  deallocate (matrot_d, matrottrans_d)
 
end subroutine test_rotmat

subroutine test_pseudo_angles
  use types
  use deco
  implicit none
  integer :: i 
  double precision :: tan_up, tan_down
  double precision, allocatable :: pseudo_angle_up(:)
  double precision, allocatable :: pseudo_angle_down(:)
  logical :: test
  
  nb_pairs = 5 

  allocate (pseudo_angle_up(nb_pairs), pseudo_angle_down(nb_pairs))
  allocate (C12(nb_pairs), DJ(nb_pairs))

  polar_up   = 0.5d0
  polar_down = -0.5d0

  C12 = (/(i, i=1,nb_pairs)/)
  DJ  = (/(i**2, i=1,nb_pairs)/)

  call pseudo_angles (pseudo_angle_up, pseudo_angle_down)

  test = .true.
  do i=1,nb_pairs
     tan_up = datan(abs(i / (polar_up * i**2)))
     tan_down = datan(abs(i / (polar_down * i**2)))
     if (i / (polar_up * i**2) .lt. 0.d0) tan_up = pi - tan_up
     if (i / (polar_down * i**2) .lt. 0.d0) tan_down = pi - tan_down
      if (pseudo_angle_up(i) .ne. tan_up) test = .false.
     if (pseudo_angle_down(i) .ne. tan_down) test = .false.

  end do
 
  if (test) write(*,*) "test pseudo_angles ok..."
  if (test == .false.) write(*,*) "Pb in pseudo_angles"

  deallocate(pseudo_angle_up, pseudo_angle_down)
  deallocate (C12, DJ)
  
end subroutine test_pseudo_angles

subroutine test_eigen_energies
  use types
  use deco
  implicit none
  integer :: i 
  double precision, allocatable :: eigen_ener_up(:)
  double precision, allocatable :: eigen_ener_down(:)
  logical :: test
  
  nb_pairs = 5 

  allocate (eigen_ener_up(nb_pairs), eigen_ener_down(nb_pairs))
  allocate (C12(nb_pairs), DJ(nb_pairs))

  polar_up   = 0.5d0
  polar_down = -0.5d0

  C12 = (/(i, i=1,nb_pairs)/)
  DJ  = (/(i**2, i=1,nb_pairs)/)

  call eigen_energies (eigen_ener_up, eigen_ener_down)

  test = .true.
  do i=1,nb_pairs
     if (eigen_ener_up(i) .ne. dsqrt(i**2 + (0.5d0 * i**2)**2)/4.d0) &
          test = .false.
     if (eigen_ener_down(i) .ne. dsqrt(i**2 + (-0.5d0 * i**2)**2)/4.d0) &
          test = .false.
  end do

  if (test) write(*,*) "test eigen_energies ok..."
  if (test == .false.) write(*,*) "Pb in eigen_energies"

  deallocate(eigen_ener_up, eigen_ener_down)
  deallocate(C12, DJ)

end subroutine test_eigen_energies
