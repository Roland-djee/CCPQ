module write
  use types
  implicit none
  ! writing formats
  character (len=*), parameter :: fmt1 = "(i5, i5, i5)"
  character (len=*), parameter :: fmt2 = "(i5, i5, i5, i5, i5, i5, es20.10e3)"
  character (len=*), parameter :: fmt3 = "(i5, i5, i5, i5, i5, i5, es20.10e3,&
       es20.10e3, es20.10e3)"
  
contains

  subroutine couplings_outputs
    implicit none
    ! local variables
    integer :: k, l, m
    
    open(14,file='C12.dat')
    open(15,file='J.dat')
    open(16,file='couplings.dat')

    m = 0
    do k=1,nb_imp - 1
       do l=k + 1,nb_imp
          m = m + 1
          write(14,fmt2)lattice_imp(k)%x, lattice_imp(k)%y, lattice_imp(k)%z,&
                        lattice_imp(l)%x, lattice_imp(l)%y, lattice_imp(l)%z,&
                        C12(m)
          write(15,fmt3)lattice_imp(k)%x, lattice_imp(k)%y, lattice_imp(k)%z,&
                        lattice_imp(l)%x, lattice_imp(l)%y, lattice_imp(l)%z,&
                        J(k), J(l), DJ(m)       
          write(16, "(es20.10e3, es20.10e3, es20.10e3)")C12(m), J(k), J(l)
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
