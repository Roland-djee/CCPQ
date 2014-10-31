module read
  use types
  implicit none
  character (len=*), parameter :: fmt1 = "(t50, i, i, i)"
  character (len=*), parameter :: fmt2 = "(t50, i)"
  character (len=*), parameter :: fmt3 = "(t50, es)"
  character (len=*), parameter :: fmt4 = "(t50, a)"
  character (len=*), parameter :: fmt5 = "(t50, l)"
    
contains

  subroutine read_outputs
    implicit none
    integer :: Output=11
    
    open(unit=Output, file="../input/Output.inp")
    read(Output, fmt5) crys_sites
    read(Output, fmt5) crys_imp
    read(Output, fmt5) coup_val
    read(Output, fmt5) deco_val
    close(Output)
    
  end subroutine read_outputs

  subroutine read_inputs_generate
    implicit none
    ! local variables
    integer :: i, System=10
    
    open(unit=System, file="../input/System.inp")
    read(System, fmt1) Orientation%x, Orientation%y, Orientation%z
    read(System, fmt2) nb_basisvec
    allocate(basis_vectors(nb_basisvec))
    do i=1,nb_basisvec
       read(System, fmt1)basis_vectors(i)%x, basis_vectors(i)%y, &
                         basis_vectors(i)%z
    end do
    read(System, fmt2) lattice%modulo
    read(System, fmt3) lattice%imp
    read(System, fmt4) lattice%impdist
    read(System, fmt2) lattice%range
    read(System, fmt3) lattice%param
    close(System)

  end subroutine read_inputs_generate

  subroutine decoherence_inputs
    implicit none
    ! local variables
    integer :: Dynamics=13
 
    open(unit=Dynamics, file="../input/Dynamics.inp")
    read(Dynamics, fmt4) Systype
    read(Dynamics, fmt4) Qubittype
    read(Dynamics, fmt4) Theotype
    read(Dynamics, fmt4) Decotype
    read(Dynamics, fmt3) Tmax
    read(Dynamics, fmt2) nb_pts_t
    read(Dynamics, fmt4) Dynadeco
    if (Dynadeco == "Hahn" .or. Dynadeco == "DD") read(Dynamics, fmt3) Tau
    if (Qubittype == "nuclear") read(Dynamics, fmt3) Jnuc
    close(Dynamics)

  end subroutine decoherence_inputs

end module read
