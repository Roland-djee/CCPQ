module types
  implicit none
  
  ! IO units, inputs 
  character(len=20) :: Systype, Qubittype, Decotype, Theotype
  character(len=20) :: Dynadeco
  integer :: nb_basisvec, nb_sites, nb_imp, nb_pairs
  integer :: origin

  integer :: nb_pts_t, CP_seq
  double precision :: Tmax
  double precision :: expval_up, expval_down
  double precision :: polar_up, polar_down
  double precision :: Jnuc

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

end module types
