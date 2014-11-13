program test_read

  implicit none

  call test_read_inputs_generate

  call test_decoherence_inputs

  call test_read_outputs

end program test_read

subroutine test_read_outputs
  use types
  use read
  implicit none

  logical :: test

  ! test read_ouputs
  test = .true.  
  call read_outputs

  if (crys_sites) then
     write(*,*)'crys_sites ok...'
  else
     write(*,*)'Problem in crys_sites'
     test = .false.
  end if
  if (crys_imp) then
     write(*,*)'crys_imp ok...'
  else
     write(*,*)'Problem in crys_imp'
     test = .false.
  end if
  if (coup_val) then
     write(*,*)'coup_val ok...'
  else
     write(*,*)'Problem in coup_val'
     test = .false.
  end if
  if (deco_val) then
     write(*,*)'deco_va ok...'
  else
     write(*,*)'Problem in deco_val'
     test = .false.
  end if

  if (test) then
     write(*,*)'test read_ouputs ok...'
  else
     write(*,*)'test read_ouputs failed...'
  end if

end subroutine test_read_outputs

subroutine test_decoherence_inputs
  use types
  use read
  implicit none

  logical :: test

  ! test decoherence_inputs with a Dynamics.inp test file
  ! in ./features
  call decoherence_inputs

  test = .true.
  if (Systype == 'Donor_Silicon') then
     write(*,*)'Donor_Silicon ok...'
  else
     write(*,*)'Problem in Donor_Silicon'
     test = .false.
  end if
  if (Qubittype == 'nuclear') then
     write(*,*)'Qubittype ok...'
  else
     write(*,*)'Problem in Qubittype'
     test = .false.
  end if
  if (Decotype == 'IFF') then
     write(*,*)'Decotype ok...'
  else
     write(*,*)'Problem in Decotype'
     test = .false.
  end if
  if (Theotype == 'Pseudospin') then
     write(*,*)'Theotype ok...'
  else
     write(*,*)'Problem in Theotype'
     test = .false.
  end if
  if (Tmax == 1.D0) then
     write(*,*)'Tmax ok...'
  else
     write(*,*)'Problem in Tmax'
     test = .false.
  end if
  if (nb_pts_t == 1000) then
     write(*,*)'nb_pts_t ok...'
  else
     write(*,*)'Problem in nb_pts_t'
     test = .false.
  end if
  if (Dynadeco == 'CP') then
     write(*,*)'Dynadeco ok...'
  else
     write(*,*)'Problem in Dynadeco'
     test = .false.
  end if
  if (CP_seq == 2) then
     write(*,*)'CP_seq ok...'
  else
     write(*,*)'Problem in CP_seq'
     write(*,*)CP_seq
     test = .false.
  end if
  if (Jnuc == 1.D6) then
     write(*,*)'Jnuc ok...'
  else
     write(*,*)'Problem in Jnuc'
     test = .false.
  end if

  if (test) then
     write(*,*)'test decoherence_inputs ok...'
  else
     write(*,*)'test decoherence_inputs failed...'
  end if

end subroutine test_decoherence_inputs

subroutine test_read_inputs_generate
  use types
  use read
  implicit none

  logical :: test
  
  ! test read_inputs_generate with a System.inp test file
  ! in ./features
  call read_inputs_generate

  test = .true.
  if (Orientation%x == 3 .and. Orientation%y == 2 .and. &
       &Orientation%z == 1) then
     write(*,*)'Orientation coordinates ok...'
  else
     write(*,*)'Problem in Orientation coordinates'
     test = .false.
  end if
  if (nb_basisvec == 8) then
     write(*,*)'nb_basisvec ok...'
  else
     write(*,*)'Problem in nb_basisvec'
     test = .false.
  end if
  if (basis_vectors(1)%x == 0 .and. basis_vectors(1)%y == 0 &
       &.and. basis_vectors(1)%z == 0) then
     write(*,*)'Vector 1 coordinates ok...'
  else
     write(*,*)'Problem in Vector 1 coordinates'
     test = .false.
  end if
  if (basis_vectors(2)%x == 0 .and. basis_vectors(2)%y == 2 &
       &.and. basis_vectors(2)%z == 2) then
     write(*,*)'Vector 2 coordinates ok...'
  else
     write(*,*)'Problem in Vector 2 coordinates'
     test = .false.
  end if
  if (basis_vectors(3)%x == 2 .and. basis_vectors(3)%y == 0 &
       &.and. basis_vectors(3)%z == 2) then
     write(*,*)'Vector 3 coordinates ok...'
  else
     write(*,*)'Problem in Vector 3 coordinates'
     test = .false.
  end if
  if (basis_vectors(4)%x == 2 .and. basis_vectors(4)%y == 2 &
       &.and. basis_vectors(4)%z == 0) then
     write(*,*)'Vector 4 coordinates ok...'
  else
     write(*,*)'Problem in Vector 4 coordinates'
     test = .false.
  end if
  if (basis_vectors(5)%x == 3 .and. basis_vectors(5)%y == 3 &
       &.and. basis_vectors(5)%y == 3) then
     write(*,*)'Vector 5 coordinates ok...'
  else
     write(*,*)'Problem in Vector 5 coordinates'
     test = .false.
  end if
  if (basis_vectors(6)%x == 3 .and. basis_vectors(6)%y == 1 &
       &.and. basis_vectors(6)%z == 1) then
     write(*,*)'Vector 6 coordinates ok...'
  else
     write(*,*)'Problem in Vector 6 coordinates'
     test = .false.
  end if
  if (basis_vectors(7)%x == 1 .and. basis_vectors(7)%y == 3 &
       &.and. basis_vectors(7)%z == 1) then
     write(*,*)'Vector 7 coordinates ok...'
  else
     write(*,*)'Problem in Vector 7 coordinates'
     test = .false.
  end if
  if (basis_vectors(8)%x == 1 .and. basis_vectors(8)%y == 1 &
       &.and. basis_vectors(8)%z == 3) then
     write(*,*)'Vector 8 coordinates ok...'
  else
     write(*,*)'Problem in Vector 8 coordinates'
     test = .false.
  end if
  if (lattice%modulo == 4) then
     write(*,*)'lattice%modulo ok...'
  else
     write(*,*)'Problem in lattice%modulo'
     test = .false.
  end if
  if (lattice%imp == 4.67D0) then
     write(*,*)'lattice%imp ok...'
  else
     write(*,*)'Problem in lattice%imp'
     test = .false.
  end if
  if (lattice%impdist == 'Uniform') then
     write(*,*)'lattice%impdist ok...'
  else
     write(*,*)'Problem in lattice%impdist'
     test = .false.
  end if
  if (lattice%range == 2) then
     write(*,*)'lattice%range ok...'
  else
     write(*,*)'Problem in lattice%range'
     test = .false.
  end if
  if (lattice%param == 5.43D0) then
     write(*,*)'lattice%range ok...'
  else
     write(*,*)'Problem in lattice%param'
     test = .false.
  end if

  if (test) then
     write(*,*)'test read_inputs_generate ok...'
  else
     write(*,*)'test read_inputs_generate failed...'
  end if

end subroutine test_read_inputs_generate

