program compare_decoherence
  use types
  use constants
  use write
  use read
  use crystal
  use interactions
  use deco
  implicit none

  !======= GENERATING THE CRYSTAL =======

  ! Inputs for generating the crystal
  call read_inputs_generate

  write(*,*)'Read inputs for generating the crystal...'

  ! Generate the crystal site coordinates
  call generate_crystal
  ! Generate the impurities coordinate
  call generate_impurities

  write(*,*)'Crystal site and impurities generated...'

  ! Write crystal ouputs
  call crystal_ouputs
  ! Write impurities ouputs
  call impurities_ouputs

  write(*,*)'Outputs written...'

  !======= DECOHERENCE DYNAMICS ========

  ! Inputs for decoherence
  call decoherence_inputs

  write(*,*)'Read inputs for decoherence...'

  ! Generate couplings
  call couplings

  write(*,*)'Couplings computed...'

  ! Write coupling values
  call couplings_outputs

  write(*,*)'Couplings written...'
  
  ! Calculate decoherence
  call decohere

  write(*,*)'Decoherence computed and written...Enjoy !'

end program compare_decoherence

