program test_interactions
  implicit none

  ! test F
  call test_F

  ! Test subroutine hyperfine
  call test_hyperfine

  ! Test subroutine dipolar
  call test_dipolar
  
end program test_interactions

subroutine test_dipolar
  use types
  use constants
  use read
  use crystal
  use interactions
  implicit none
  integer :: i,l,k
  double precision :: R12

  Orientation%x = 1
  Orientation%y = 1
  Orientation%z = 1
  
  nb_imp = 10
  nb_pairs = nb_imp * (nb_imp - 1) / 2
  allocate(lattice_imp(nb_imp))
  do i=1,nb_imp
     lattice_imp(i)%x = i
     lattice_imp(i)%y = i
     lattice_imp(i)%z = i
  end do
  call dipolar
  open(10,file='C12.dat')
  open(20,file='C12_analytical.dat')
  l = 0
  do i=1,nb_imp - 1
     do k=i + 1,nb_imp
        l = l + 1
        write(10,*)scaling*sqrt( dble(lattice_imp(i)%x - lattice_imp(k)%x)**2 &
             & + dble(lattice_imp(i)%y - lattice_imp(k)%y)**2 &
             & + dble(lattice_imp(i)%z - lattice_imp(k)%z)**2),C12(l)
     end do
  end do

  do i=1,500
     R12 = dble(i) / 10.d0
     write(20,*)R12, pref * (-2.d0) / R12**3
  end do

  close(10)
  close(20)
  write(*,*)'C12.dat C12_analytical.dat generated'

end subroutine test_dipolar

subroutine test_hyperfine
  use types
  use constants
  use read
  use crystal
  use interactions
  implicit none
  integer :: i,k,l
  logical :: test

  double precision, allocatable :: F_test1(:), F_test2(:), F_test3(:)
  double precision, allocatable :: J_test(:)
  double precision, allocatable :: x_test(:), y_test(:), z_test(:)
  double precision :: arg1, arg2, arg3

  nb_imp = 3
  nb_pairs = nb_imp * (nb_imp - 1) / 2
  allocate (lattice_imp(nb_imp))
  allocate (F_test1(nb_imp), F_test2(nb_imp), F_test3(nb_imp))
  allocate (J_test(nb_imp))
  allocate (x_test(nb_imp), y_test(nb_imp), z_test(nb_imp))
  do i=1,nb_imp
     lattice_imp(i)%x = i
     lattice_imp(i)%y = 2*i
     lattice_imp(i)%z = i**2
  end do

  call hyperfine

  do i=1,nb_imp

     x_test(i) = scaling*lattice_imp(i)%x
     y_test(i) = scaling*lattice_imp(i)%y
     z_test(i) = scaling*lattice_imp(i)%z
 
     arg1 = x_test(i)**2 / ((n * b)**2) + &
          (y_test(i)**2 + z_test(i)**2) /((n * a)**2)
     arg2 = y_test(i)**2 / ((n * b)**2) + &
          (z_test(i)**2 + x_test(i)**2) /((n * a)**2)
     arg3 = z_test(i)**2 / ((n * b)**2) + &
          (x_test(i)**2 + y_test(i)**2) /((n * a)**2)

     F_test1(i) = dexp(-dsqrt(arg1)) / dsqrt(pi * (n * a)**2 * n * b)
     F_test2(i) = dexp(-dsqrt(arg2)) / dsqrt(pi * (n * a)**2 * n * b)
     F_test3(i) = dexp(-dsqrt(arg3)) / dsqrt(pi * (n * a)**2 * n * b)

     J_test(i)  = F_test1(i)*dcos(k0*x_test(i)) + &
                  F_test2(i)*dcos(k0*y_test(i)) + &
                  F_test3(i)*dcos(k0*z_test(i))

     J_test(i)  = p*J_test(i)**2

  end do

  deallocate (lattice_imp)

  test = .true.
  do i=1,nb_imp
     if (abs(J(i) - J_test(i)) .ne. 0.d0) test = .false.
     !write(*,*)J(i),J_test(i)
  end do

  if (test) then
     write(*,*)'J ok...'
  else
     write(*,*)'Pb in J'
  end if

  test = .true.

  l = 0
  do i=1,nb_imp - 1
     do k=i + 1,nb_imp
        l = l + 1
        if (abs(DJ(l) - (J_test(i) - J_test(k)) ) .ne. 0.d0) test = .false.
     end do
  end do

  if (test) then
     write(*,*)'DJ ok...'
  else
     write(*,*)'Pb in DJ'
  end if

end subroutine test_hyperfine

subroutine test_F
  use types
  use constants
  use read
  use crystal
  use interactions
  implicit none
  integer :: i
  logical :: test

  double precision, allocatable :: F_test1(:), F_test2(:), F_test3(:)
  double precision, allocatable :: J_test(:)
  double precision, allocatable :: x_test(:), y_test(:), z_test(:)
  double precision :: arg1

  ! Test function F

  nb_imp = 3
  allocate (F_test1(nb_imp), F_test2(nb_imp), F_test3(nb_imp))
  allocate (J_test(nb_imp))
  allocate (x_test(nb_imp), y_test(nb_imp), z_test(nb_imp))

  do i=1,nb_imp
     x_test(i)=dble(i)
     y_test(i)=dble(i)+1.d0
     z_test(i)=dble(i)**2
     !x_test(i)=0.d0
     !y_test(i)=0.d0
     !z_test(i)=0.d0

     arg1 = x_test(i)**2 / ((n * b)**2) + &
          (y_test(i)**2 + z_test(i)**2) /((n * a)**2)

     F_test1(i) = dexp(-dsqrt(arg1)) / dsqrt(pi * (n * a)**2 * n * b)

  end do
  
  F_test2 = F(x_test, y_test, z_test)
  
  test = .true.
  do i=1,nb_imp
     if (abs(F_test1(i) - F_test2(i)) .ne. 0.d0) test = .false.
     !write(*,*)F_test1(i),F_test2(i)
  end do

  if (test) then
     write(*,*)'F function ok...'
  else
     write(*,*)'Pb in F function'
  end if

end subroutine test_F

