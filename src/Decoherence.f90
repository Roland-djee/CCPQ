program decoherence
  use types
  use read
  use crystal
  use write
  use interactions
  use deco
  implicit none
  
  !======= GENERATING THE CRYSTAL =======

  ! Inputs for generating the crystal
  call read_inputs_generate

  write(*,*)'Read inputs for generating the crystal...'

  ! Outputs parameters
  call read_outputs

  write(*,*)'Read outputs parameters...'

  ! Generate the crystal site coordinates
  call generate_crystal
  ! Generate the impurities coordinate
  call generate_impurities

  write(*,*)'Crystal site and impurities generated...'

  ! Write crystal ouputs
  if (crys_sites) call crystal_ouputs
  ! Write impurities ouputs
  if (crys_imp) call impurities_ouputs

  if (crys_sites .or. crys_imp) write(*,*)'Outputs written...'

  !======= DECOHERENCE DYNAMICS ========

  ! Inputs for decoherence
  call decoherence_inputs

  write(*,*)'Read inputs for decoherence...'

  ! Generate couplings
  call couplings

  write(*,*)'Couplings computed...'

  ! Write coupling values
  if (coup_val) call couplings_outputs
  if (coup_val) write(*,*)'Couplings written...'
  
  ! Calculate decoherence
  call decohere

  write(*,*)'Decoherence computed and written...Enjoy !'

end program decoherence
