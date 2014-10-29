program test_crystal
  use types
  use read
  use crystal
 
  implicit none
  integer :: k,l,m,n
  logical :: test

  ! test generate_crystal

  lattice%range  = 2    ! For the sake of testing lattice%range >= 2
  lattice%modulo = 4
  nb_basisvec    = 1
  allocate (basis_vectors(nb_basisvec))
  basis_vectors(1)%x = 1
  basis_vectors(1)%y = 1
  basis_vectors(1)%z = 1
  
  call generate_crystal

  test = .true.
  n = 0
  do k=-lattice%range,lattice%range-1
     do l=-lattice%range,lattice%range-1
        do m=-lattice%range,lattice%range-1
           n = n + 1 
           if (lattice_sites(n)%x == basis_vectors(1)%x + k*lattice%modulo) then
              write(*,*)'lattice_site%x ok...'
           else
              write(*,*)'Pb in lattice_site%x...'
              test = .false.
           end if
           if (lattice_sites(n)%y == basis_vectors(1)%y + l*lattice%modulo) then
              write(*,*)'lattice_site%y ok...'
           else
              write(*,*)'Pb in lattice_site%y...'
              test = .false.
           end if
           if (lattice_sites(n)%z == basis_vectors(1)%z + m*lattice%modulo) then
              write(*,*)'lattice_site%z ok...'
           else
              write(*,*)'Pb in lattice_site%z...'
              test = .false.
           end if           
        end do
     end do
  end do

  if (test) then
     write(*,*)'lattice_sites ok..'
  else
     write(*,*)'Pb in lattice_site'
  end if

  ! test generate_impurities
  lattice%imp = 50.d0
  call generate_impurities

  test = .true.
  do k=1,nb_imp
     if(lattice_imp(k)%x == 0 .and. lattice_imp(k)%y == 0 &
          .and. lattice_imp(k)%z == 0) then
        write(*,*)'Impurity on the origin'
        test = .false.
     end if
     do l=k+1,nb_imp
        if(lattice_imp(k)%x == lattice_imp(l)%x .and. &
           lattice_imp(k)%y == lattice_imp(l)%y .and. &
           lattice_imp(k)%z == lattice_imp(l)%z) then
           write(*,*)'Identical impurities'
           test = .false. 
           !write(*,*)'k,l',k,l
        end if
     end do
  end do
     
  if (test) then
     write(*,*)'lattice_imp ok...'
  else
     write(*,*)'Pb in lattice_imp...'    
  end if

end program test_crystal

