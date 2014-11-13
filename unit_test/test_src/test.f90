 use crystal
  use interactions

  ! Generate the crystal
  !call generate_crystal

  ! test function F for hyperfine couplings
  !call test_F

  ! test subroutine hyperfine
  !call test_hyperfine

  ! test subroutine dipolar
  call test_dipolar
subroutine test_dipolar
  use types
  use interactions
  implicit none
  integer :: i, x0, y0, z0, proj
  double precision :: Orient_norm, dist, cos12

  allocate (lattice_imp(3))

  nb_imp = 3

  lattice_imp%x = (/1, 2, 3/)
  lattice_imp%y = (/4, 5, 6/)
  lattice_imp%z = (/7, 8, 9/)
  
  x0 = 0
  y0 = 0
  z0 = 0

  Orientation%x = 1
  Orientation%y = 0
  Orientation%z = 0

  Orient_norm = dsqrt(dble(Orientation%x**2) + &
       dble(Orientation%y**2) + &
       dble(Orientation%z**2))

  call dipolar (x0, y0, z0, Orient_norm)

  print*,"Testing..."

  do i=1,10

  dist = dsqrt(dble(i-x0)**2 + dble(0-y0)**2 + dble(0-y0)**2)
  print*,dist
  proj = Orientation%x*(i-x0) + Orientation%y*(0-y0) + Orientation%z*(0-y0)
  print*,proj
  cos12  = proj / (Orient_norm * dist)
  print*, cos12
  C12 = pref * (1.d0 - 3.d0 * cos12**2) / ((scaling * dist)**3)
  print*,C12
  write(10,*)scaling * dist, C12

  end do

  

end subroutine test_dipolar

subroutine test_hyperfine
  use types
  use interactions
  implicit none
  double precision :: x1, y1, z1, arg1, F1
  double precision :: x2, y2, z2, arg2, F2
  double precision :: x3, y3, z3, arg3, F3
  double precision :: t1, t2, t3
  
  allocate (lattice_imp(3))

  nb_imp = 3

  lattice_imp%x = (/1, 2, 3/)
  lattice_imp%y = (/4, 5, 6/)
  lattice_imp%z = (/7, 8, 9/)

  call hyperfine

  print*,"Testing.."

  x1 = 1.d0*scaling
  y1 = 4.d0*scaling
  z1 = 7.d0*scaling

  print*,x1
  print*,y1
  print*,z1

  arg1 = dsqrt( (x1 / (n*b))**2 + (y1**2 + z1**2) / ((n*a)**2) )
  F1   = dexp(-arg1) / dsqrt( pi * (n*a)**2 * (n*b) )
  t1 = F1 * cos(k0 * x1)
  print*, t1

  arg2 = dsqrt( (y1 / (n*b))**2 + (z1**2 + x1**2) / ((n*a)**2) )
  F2   = dexp(-arg2) / dsqrt( pi * (n*a)**2 * (n*b) )
  t2 = F2 * cos(k0 * y1)
  print*, t2

  arg3 = dsqrt( (z1 / (n*b))**2 + (x1**2 + y1**2) / ((n*a)**2) )
  F3   = dexp(-arg3) / dsqrt( pi * (n*a)**2 * (n*b) )
  t3 = F3 * cos(k0 * z1)
  print*, t3
  print*, t1 + t2 + t3
  print*, p*(t1 + t2 + t3)**2


  !x2 = 2.d0*scaling
  !y2 = 5.d0*scaling
  !z2 = 8.d0*scaling

  !print*,x2
  !print*,y2
  !print*,z2

  !arg2 = dsqrt( (x2 / (n*b))**2 + (y2**2 + z2**2) / ((n*a)**2) )
  !F2   = dexp(-arg2) / dsqrt( pi * (n*a)**2 * (n*b) )

  !t2 = F2 * cos(k0 * x1)
  !print*, t1

end subroutine test_hyperfine

subroutine test_F
  use types
  use interactions
  implicit none
  double precision :: X(3), Y(3), Z(3)
  double precision :: x1, y1, z1, arg1, F1
  double precision :: x2, y2, z2, arg2, F2
  double precision :: x3, y3, z3, arg3, F3

  nb_imp = 3

  X = (/1.d0, 2.d0, 3.d0/)
  Y = (/4.d0, 5.d0, 6.d0/)
  Z = (/7.d0, 8.d0, 9.d0/)
  
  print*,a,b,n,pi
  print*,F(X, Y, Z)

  x1 = 1.d0
  y1 = 4.d0
  z1 = 7.d0

  arg1 = dsqrt( (x1 / (n*b))**2 + (y1**2 + z1**2) / ((n*a)**2) )
  F1   = dexp(-arg1) / dsqrt( pi * (n*a)**2 * (n*b) )
  print*,f1

  x2 = 2.d0
  y2 = 5.d0
  z2 = 8.d0

  arg2 = dsqrt( (x2 / (n*b))**2 + (y2**2 + z2**2) / ((n*a)**2) )
  F2   = dexp(-arg2) / dsqrt( pi * (n*a)**2 * (n*b) )
  print*,f2

  x3 = 3.d0
  y3 = 6.d0
  z3 = 9.d0

  arg3 = dsqrt( (x3 / (n*b))**2 + (y3**2 + z3**2) / ((n*a)**2) )
  F3   = dexp(-arg3) / dsqrt( pi * (n*a)**2 * (n*b) )
  print*,f3

end subroutine test_F

