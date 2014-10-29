module crystal
  use types
  implicit none

contains 

  subroutine generate_crystal
    implicit none
    ! local variables
    integer :: i, j, k, l, m, N0, N1
    integer :: volume

    ! define the lattice volume and the number of sites
    volume   = (2*lattice%range)**3
    if (lattice%range == 0) volume = 1
    nb_sites = nb_basisvec*volume
    
    ! allocate the sites coordinates array
    allocate (lattice_sites(nb_sites))

    ! generate the crystal sites coordinates
    m = 0
    N0 = lattice%range
    N1 = lattice%range
    if (lattice%range == 0) then
       N0 = 0
       N1 = 1
    end if
    do i=-N0,N1-1
       do j=-N0,N1-1
          do k=-N0,N1-1
             do l=1,nb_basisvec
                m = m + 1
                lattice_sites(m)%x = basis_vectors(l)%x + i*lattice%modulo
                lattice_sites(m)%y = basis_vectors(l)%y + j*lattice%modulo
                lattice_sites(m)%z = basis_vectors(l)%z + k*lattice%modulo
                ! Select the origin indice
                if (lattice_sites(m)%x == 0 .and. &
                    lattice_sites(m)%y == 0 .and. &
                    lattice_sites(m)%z == 0) origin = m
             end do
          end do
       end do
    end do

  end subroutine generate_crystal

  subroutine generate_impurities
    implicit none
    ! local variables
    integer :: i
    integer, allocatable :: imp_temp(:), index(:)
    double precision :: rand2(1)
    double precision, allocatable :: rand1(:)
    logical, allocatable :: mask(:)

    ! define the number of impurities
    nb_imp = nint(lattice%imp*nb_sites/1.D2)

    ! allocate the impurities , random 
    allocate (imp(nb_imp))
    allocate (rand1(nb_imp))
    allocate (mask(nb_imp))

    ! define randomly nb_imp sites as impurities among the nb_sites
    ! excluding 0
    call random_seed()
    call random_number(rand1)
    imp = nint( rand1 * (nb_sites - 1) + 1)

    ! if so, remove the origin and replace by another random site
    call random_number(rand2)
    where (imp == origin) 
       imp = nint ( rand2 * (nb_sites - 1) + 1)
    end where

    ! check and remove duplicates if any
    mask = .true.
    do i=nb_imp,2,-1
       mask(i) = .not.(any(imp(:i-1) == imp(i)))
    end do
    allocate (index, source=pack([(i,i=1,nb_imp)], mask))
    allocate (imp_temp, source=imp(index))

    ! reset parameters
    nb_imp = size(imp_temp)
    deallocate(imp)
    allocate(imp(nb_imp))
    imp = imp_temp
    
    ! allocate coordinates arrays
    allocate (lattice_imp(nb_imp))

    ! impurities coordinates
    lattice_imp%x = lattice_sites(imp)%x
    lattice_imp%y = lattice_sites(imp)%y
    lattice_imp%z = lattice_sites(imp)%z

    deallocate (rand1)
    deallocate (mask)
    deallocate (imp_temp)
    deallocate (index)

  end subroutine generate_impurities

end module crystal
