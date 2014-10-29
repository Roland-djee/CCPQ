module types
  implicit none
  
  ! IO units, inputs 
  character(len=20) :: Systype, Qubittype, Decotype, Theotype
  character(len=20) :: Dynadeco
  integer :: nb_basisvec, nb_sites, nb_imp, nb_pairs
  integer :: origin

  integer :: nb_pts_t
  double precision :: Tmax, polar, Jnuc

  logical :: crys_sites
  logical :: crys_imp
  logical :: coup_val
  logical :: deco_val

  ! Arrays
  integer, allocatable :: imp(:)
  double precision, allocatable :: J(:), DJ(:), C12(:)

  ! Vector type
  type vectors
     integer :: x, y, z
  end type vectors
  type (vectors) :: Orientation
  type (vectors), allocatable :: basis_vectors(:)
  type (vectors), allocatable :: lattice_sites(:)
  type (vectors), allocatable :: lattice_imp(:)

  ! Lattice type
  type lat
     integer           :: modulo, range   
     double precision  :: imp, param
     character(len=20) :: impdist
  end type lat
  type (lat) :: lattice

  ! 2x2 Rotation matrices
  type rot
     double complex    :: elements(2, 2) 
  end type rot
  type (rot), allocatable :: matrot_u(:), matrottrans_u(:)
  type (rot), allocatable :: matrot_d(:), matrottrans_d(:)
  type (rot), allocatable :: Zgate_u(:), Zgate_d(:)
  type (rot), allocatable :: Tu(:), Td(:)

end module types
