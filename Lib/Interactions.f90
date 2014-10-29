module interactions
  use types
  use constants
  implicit none

contains

  subroutine couplings
    implicit none
    double precision :: pref
 
    ! define the number of impurity pairs
    nb_pairs = nb_imp * (nb_imp - 1) / 2

    ! calculate hyperfine coupling values
    call hyperfine

    ! set the prefactor for the dipolar interaction
    select case (Qubittype)
    case ("electron")
       pref = pref_e
     case ("nuclear")
       pref = pref_n
    case ("mixed")
       write(*,*)"To be continued..." 
       stop
    case default
       write(*,*)"No valid Qubit type in Dynamics.inp" 
       stop
    end select
       
    ! calculate dipolar coupling values between all impurity pairs
    call dipolar (pref)
    
  end subroutine couplings

  subroutine dipolar (pref)
    implicit none    
    double precision, intent(in) :: pref
    ! Local variables and arrays
    integer :: i, j, k
    double precision :: Orient_norm
    integer :: proj(nb_pairs)
    integer :: X_dist(nb_pairs)
    integer :: Y_dist(nb_pairs)
    integer :: Z_dist(nb_pairs)
    double precision :: dist(nb_pairs)
    double precision :: cos_t12(nb_pairs)

    ! allocations
    allocate (C12(nb_pairs))

    ! Orientation norm
    Orient_norm = dsqrt(dble(Orientation%x**2) + &
                        dble(Orientation%y**2) + &
                        dble(Orientation%z**2))
    
    ! coordinate difference between every impurity spin pairs
    ! avoiding double counting
    k = 0
    do i=1,nb_imp - 1
       do j=i + 1,nb_imp
          k = k + 1
          X_dist(k) = lattice_imp(j)%x - lattice_imp(i)%x
          Y_dist(k) = lattice_imp(j)%y - lattice_imp(i)%y
          Z_dist(k) = lattice_imp(j)%z - lattice_imp(i)%z
       end do
    end do
  
    ! distance between every spin partners of impurity pairs
    dist    = dsqrt(dble(X_dist**2) + dble(Y_dist**2) + dble(Z_dist**2))
    ! Projection
    proj    = Orientation%x*X_dist + Orientation%y*Y_dist + Orientation%z*Z_dist
    ! cos(theta_12)
    cos_t12 = proj / (Orient_norm * dist)
    
    ! array of C12 values for each impurity pair
    C12 = pref * (1.d0 - 3.d0 * cos_t12**2) / ((scaling * dist)**3)
  
  end subroutine dipolar

  subroutine hyperfine
    implicit none
    ! Local variables and arrays
    integer :: k, l, m
    double precision :: X_imp(nb_imp), Y_imp(nb_imp), Z_imp(nb_imp)
    double precision :: term1(nb_imp), term2(nb_imp), term3(nb_imp)
 
    allocate (J(nb_imp))
    allocate (DJ(nb_pairs))

    ! scale impurities coordinates in angstroms
    X_imp = scaling * lattice_imp%x
    Y_imp = scaling * lattice_imp%y
    Z_imp = scaling * lattice_imp%z
   
    ! compute the terms of the KLW expression
    term1 = F(X_imp, Y_imp, Z_imp) * dcos(k0 * X_imp)
    term2 = F(Y_imp, Z_imp, X_imp) * dcos(k0 * Y_imp)
    term3 = F(Z_imp, X_imp, Y_imp) * dcos(k0 * Z_imp)

    ! array of hyperfine values for each impurity
    J = term1 + term2 + term3
    J = p * J**2
    
    ! array of hyperfine differences for each impurity pair
    m = 0
    do k=1,nb_imp - 1
       do l=k + 1,nb_imp
          m = m + 1
          DJ(m) = J(k) - J(l)
       end do
    end do
           
  end subroutine hyperfine
 
  function F(x, y, z)
    implicit none
    double precision :: F(nb_imp)
    double precision, intent(in) :: x(nb_imp), y(nb_imp), z(nb_imp)
    ! Local array
    double precision :: arg(nb_imp)

    arg = dsqrt((x / (n * b))**2 + (y**2 + z**2) / ((n * a)**2))
    F   = dexp(-arg) / dsqrt(pi * (n * a)**2 * n * b)

  end function F
  
end module interactions
