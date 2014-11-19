module write
  use types
  implicit none
  ! writing formats
  character (len=*), parameter :: fmt01 = "(a1, t10, a, t25, a, t40, a)"
  character (len=*), parameter :: fmt02 = "(a1, t5, a, t10, a, t15, a,&
                                            t20, a, t25, a, t30, a, t40, a)"
  character (len=*), parameter :: fmt03 = "(a1, t10, a, t25, a, t41, a,&
                                            t61, a, t81, a)"
  character (len=*), parameter :: fmt04 = "(a1, t5, a, t10, a, t15, a,&
                                            t20, a, t25, a, t30, a, &
                                            t41, a, t61, a, t81, a)"
  character (len=*), parameter :: fmt05 = "(a1, t11, a, t31, a, t51, a, t65, a)"
  character (len=*), parameter :: fmt1 = "(i5, i5, i5)"
  character (len=*), parameter :: fmt2 = "(i5, i5, i5, i5, i5, i5, es20.10e3)"
  character (len=*), parameter :: fmt3 = "(i5, i5, i5, i5, i5, i5, es20.10e3,&
       es20.10e3, es20.10e3)"
  character (len=*), parameter :: fmt4 = "(es20.10e3, es20.10e3, es20.10e3)"
  
contains

  subroutine couplings_nuclear_outputs (CA1)
    implicit none
    ! local variables
    integer :: k,l,m  
    double precision, intent(in) :: CA1(nb_imp - 1)
    character(len=100) :: filename

    write(filename, '(a, es10.3e3, a)') "couplings_nuclear_J=",Jnuc,"_rad.dat"
    open(14,file=filename)
    
    ! Header
    write(14, fmt05)"#","CA1","CA2","C12","[rad/s]"

    m = 0
    do k=1,nb_imp - 2
       do l=k + 1,nb_imp - 1
          m = m + 1
          write(14, fmt4)CA1(k), CA1(l), C12(m)
       end do
    end do

  end subroutine couplings_nuclear_outputs

  subroutine couplings_outputs
    implicit none
    ! local variables
    integer :: k, l, m
    
    open(14,file='C12.dat')
    open(15,file='J.dat')
    open(16,file='couplings.dat')
    
    ! Headers
    write(14, fmt01)"#","spin 1","spin 2","C12"
    write(14, fmt02)"#","nx,","ny,","nz","nx,","ny,","nz","[rad/s]"
    write(15, fmt03)"#","spin 1","spin 2","J1","J2","DJ"
    write(15, fmt04)"#","nx,","ny,","nz","nx,","ny,","nz",&
                    "[rad/s]","[rad/s]","[rad/s]"
    write(16, fmt05)"#","J1","J2","C12","[rad/s]"

    m = 0
    do k=1,nb_imp - 1
       do l=k + 1,nb_imp
          m = m + 1 
          write(14, fmt2)lattice_imp(k)%x, lattice_imp(k)%y, lattice_imp(k)%z,&
                        lattice_imp(l)%x, lattice_imp(l)%y, lattice_imp(l)%z,&
                        C12(m)
          write(15, fmt3)lattice_imp(k)%x, lattice_imp(k)%y, lattice_imp(k)%z,&
                        lattice_imp(l)%x, lattice_imp(l)%y, lattice_imp(l)%z,&
                        J(k), J(l), DJ(m)       
          write(16, fmt4)J(k), J(l), C12(m)
       end do
    end do
    close(14)
    close(15)
    close(16)
  end subroutine couplings_outputs

  subroutine crystal_ouputs
    implicit none
    ! local variables
    integer :: i
    
    open(11, file="Crystal_sites_coord.dat")
    do i=1,nb_sites
       write(11, fmt1)lattice_sites(i)%x, lattice_sites(i)%y, lattice_sites(i)%z
    end do
    close(11)
  end subroutine crystal_ouputs

  subroutine impurities_ouputs
    implicit none
    integer :: i
    
    open(12, file="Impurities_coord.dat")
    do i=1,nb_imp
       write(12, fmt1)lattice_imp(i)%x, lattice_imp(i)%y, lattice_imp(i)%z
    end do
    close(12)
  end subroutine impurities_ouputs

end module write
