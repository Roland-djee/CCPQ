module deco
  use types
  implicit none

contains

  subroutine decohere
    implicit none

    ! set the polarisation of the qubit under study
    select case (Qubittype)
    case ("electron")
       polar = 0.5d0
    case ("nuclear")
       polar = 0.5d0
    case ("mixed")
       write(*,*)"To be continued..." 
       stop
    case default
       write(*,*)"No valid Qubit type in Dynamics.inp" 
       stop
    end select

    ! set the dynamics type
    select case (Theotype)
    case ("Pseudospin")
       call pseudo
    case ("CCE2")
       write(*,*)"To be continued..." 
       stop
    case default
       write(*,*)"No valid Theory type in Dynamics.inp" 
       stop
    end select

  end subroutine decohere

  subroutine pseudo
    implicit none
    double precision :: CA1(nb_imp - 1)

    ! set the dynamics type
    select case (Decotype)
    case ("DFF")
       if (Qubittype == "electron" .or. Qubittype == "mixed") then
          write(*,*)"DFF is not supported for electron/mixed decoherence..."
          stop
       else
          call select_nuc (CA1)
          call DFF (CA1) 
       end if
    case ("IFF")
       if (Qubittype == "electron" .or. Qubittype == "mixed") then
          call IFF
       else
          call select_nuc (CA1)
          call IFF
       end if
    case ("IDFF")
       if (Qubittype == "electron" .or. Qubittype == "mixed") then
          write(*,*)"IDFF is only supported for nuclear decoherence..."
          stop
       else
          call select_nuc (CA1) 
          call DFF (CA1) 
          call IFF
       end if
     case default
        write(*,*)"No valid Decoherence type in Dynamics.inp" 
        stop
    end select
  
  end subroutine pseudo

  subroutine select_nuc (CA1)
    implicit none
    double precision, intent(out) :: CA1(nb_imp - 1)
    ! Local variables
    integer :: loc, k, l, m, n
    integer :: Np
    double precision, allocatable :: C12_tmp(:)
    double precision, allocatable :: DC_tmp(:)

    Np = (nb_imp - 1) * (nb_imp - 2) / 2
    allocate(C12_tmp(Np), DC_tmp(Np))
  
    ! locate the nearest impurity to the input J value
    loc  = minloc(abs(J - Jnuc), 1)
    Jnuc = J(loc)

    ! select corresponding C12 values between qubit A and
    ! all other impurities for DFF
    m = 0
    n = 0
    do k=1,nb_imp - 1
       do l=k + 1,nb_imp
          m = m + 1
          if ((k .ne. loc) .and. (l == loc)) then
             n = n + 1
             CA1(n) = C12(m)
             C12(m) = 0.d0
          else if (k == loc) then
             n = n +1
             CA1(n) = C12(m)
             C12(m) = 0.d0
          end if
       end do
    end do

    ! creating a temporary array of C12 differences between qubit A and
    ! all impurity pair
    m = 0
    do k=1,nb_imp - 2
       do l=k + 1,nb_imp - 1
          m = m + 1
          DC_tmp(m) = CA1(k) - CA1(l)
       end do
    end do

    ! creating a temporary C12 array excluding previous CA1 values
    C12_tmp = pack (C12, C12 .ne. 0.d0)

    deallocate (C12, DJ)
    allocate (C12(Np), DJ(Np))

    ! reset C12, J arrays and nb_pairs for IFF
    C12      = C12_tmp
    DJ       = DC_tmp
    nb_pairs = Np 

    deallocate(C12_tmp, DC_tmp)
  
  end subroutine select_nuc
  
  subroutine DFF (CA1)
    implicit none
    double precision, intent(in) :: CA1(nb_imp - 1)
    ! Local variables
    integer :: i
    double precision :: dt, t, L_av
    double precision :: L_pairs(nb_imp - 1)
    character(len=*), parameter :: fmt="(es20.10e3, es20.10e3)"
    character(len=100) :: filename

    ! Loop over time window
    dt = Tmax / dble(nb_pts_t)
    t  = - dt
    
    write(filename, '(a, es10.3e3, a)') "Averaged_decay_nuclear_DFF_J=",&
                                         Jnuc,"_rad.dat"
    open(16, file=filename)

    do i=1,nb_pts_t + 1
       t = t + dt       
       L_pairs = dcos(0.5d0 * CA1 * t)
       L_av    = product (L_pairs)
       write(16,fmt)t, abs(L_av)
    end do
    close(16)
    
  end subroutine DFF

  subroutine IFF
    implicit none
    ! Local variables and arrays
    integer :: i, j
    double precision :: L

    double precision :: eigen_ener(nb_pairs)
    double precision :: pseudo_angle(nb_pairs)
    double precision :: L_pairs(nb_pairs)
    character(len=*), parameter :: fmt="(es20.10e3, es20.10e3)"
   
    allocate (matrot_u(nb_pairs))
    allocate (matrot_d(nb_pairs))
    allocate (matrottrans_u(nb_pairs))
    allocate (matrottrans_d(nb_pairs))
    allocate (Zgate_u(nb_pairs), Zgate_d(nb_pairs))
    allocate (Tu(nb_pairs), Td(nb_pairs))
    
    ! Compute the eigenenergies for all pairs
    call eigen_energies (eigen_ener)
    
    ! Compute the pseudospin angles for all pairs
    call pseudo_angles (pseudo_angle)

    !open(17, file='eigen_energies.dat')
    !open(18, file='pseudo_angles.dat')
    !m = 0
    !do k=1,nb_imp - 1
    !   do l1=k + 1,nb_imp
    !      m = m + 1
    !      write(17,fmt)eigen_ener(m)
    !      write(18,fmt)pseudo_angle(m)
    !   end do
    !end do

    select case (Dynadeco)
    case ("FID")
       call FID (eigen_ener, pseudo_angle)
    case ("Hahn")
       call Hahn
    case ("DD")
       write(*,*)"To be finished..."
       stop
    case default
       write(*,*)"No valid Dynamical decoupling type in Dynamics.inp" 
       stop
    end select

    ! Compute rotation matrices for all pairs in up/down states
    !call rotmat (pseudo_angle, matrot_u, matrottrans_u)
    !call rotmat (-pseudo_angle, matrot_d, matrottrans_d)

    !open(16, file='FID_Check.dat')

    ! Loop over time window
    !dt = Tmax / dble(nb_pts_t)
    !t  = - dt
    !do j=1,nb_pts_t + 1
    !   t = t + dt

       ! Compute Zgates for all pairs in up/down states
    !   call Z_gate (eigen_ener, t, Zgate_u)
    !   call Z_gate (eigen_ener, t, Zgate_d)

       ! work out transition matrices in up/down states
       ! and compute the decoherence
    !   do i=1,nb_pairs
          ! rotate to eigenbasis and propagate
    !      Tu(i)%elements = matmul(Zgate_u(i)%elements, matrot_u(i)%elements)
          ! rotate back to bath basis
    !      Tu(i)%elements = matmul(matrottrans_u(i)%elements, Tu(i)%elements)
          ! Idem down state
    !      Td(i)%elements = matmul(Zgate_d(i)%elements, matrot_d(i)%elements)
    !      Td(i)%elements = matmul(matrottrans_d(i)%elements, Td(i)%elements)
          ! Decoherence from initial |down-up> bath state
    !      L_pairs(i) = abs(conjg(Td(i)%elements(1, 1)) * Tu(i)%elements(1, 1) &
    !           + conjg(Td(i)%elements(2, 1)) * Tu(i)%elements(2, 1))
          
          ! average over the bath states
    !      L_pairs(i) = 0.5d0 + 0.5d0 * L_pairs(i)
    !   end do

       ! Final decay as the product over all pair decays
    !   L = product(L_pairs) 
       
       !write the output
    !   write(16, fmt)t, L

    !end do
    !close(16)

    deallocate (matrot_u)
    deallocate (matrot_d)
    deallocate (matrottrans_u)
    deallocate (matrottrans_d)
    deallocate (Zgate_u, Zgate_d)
    deallocate (Tu, Td)
    
  end subroutine IFF
  
  subroutine Hahn
    implicit none

  end subroutine Hahn

  subroutine FID (eigen_ener, pseudo_angle)
    implicit none
    double precision, intent(in) :: eigen_ener(nb_pairs)
    double precision, intent(in) :: pseudo_angle(nb_pairs)
    ! Local variables and arrays
    integer :: i
    double precision :: dt, t, L_av
    double precision :: A(nb_pairs), B(nb_pairs), C(nb_pairs)
    double precision :: cos_t(nb_pairs), sin_t(nb_pairs)
    double precision :: L(nb_pairs)
    character(len=*), parameter :: fmt="(es20.10e3, es20.10e3)"
    character(len=100) :: filename

    if (Qubittype == "electron") then
       write(filename, '(a)') "Averaged_decay_electron_IFF_FID.dat"
    else if (Qubittype == "nuclear") then
       write(filename, '(a, es10.3e3, a)') "Averaged_decay_nuclear_IFF_J=",&
                                         Jnuc,"_rad_FID.dat"
    end if
    open(16, file=filename)

    ! Loop over time window
    dt = Tmax / dble(nb_pts_t)
    t  = - dt

    do i=1,nb_pts_t + 1
       t = t + dt
       A = dcos(eigen_ener * t)
       C = dsin(eigen_ener * t)
       cos_t = dcos(pseudo_angle/2.d0)
       sin_t = dsin(pseudo_angle/2.d0)
       B = C * (sin_t**2 - cos_t**2)

       ! Decoherence
       L = abs(A**2 + B**2 - (2.d0 * cos_t * sin_t * C)**2)
       ! Average decoherence over bath states
       L = 0.5d0 + 0.5d0 * L
       ! Product over all pairs
       L_av = product(L)
       write(16, fmt)t, L_av    
    end do   
    close(16)

  end subroutine FID

  subroutine Z_gate (eigen_ener, t, Zgate)
    implicit none
    double precision, intent(in) :: eigen_ener(nb_pairs)
    double precision, intent(in) :: t
    type (rot), intent(out) :: Zgate(nb_pairs)
    ! Local array 
    double complex :: phase(nb_pairs)
    
    phase = dcmplx(0.d0, eigen_ener)
    Zgate%elements(1, 1) = exp(-phase*t)
    Zgate%elements(1, 2) = dcmplx(0.d0, 0.d0)
    Zgate%elements(2, 1) = dcmplx(0.d0, 0.d0)
    Zgate%elements(2, 2) = exp(phase*t)
        
  end subroutine Z_gate

  subroutine rotmat (pseudo_angle, matrot, matrottrans)
    implicit none
    double precision, intent(in) :: pseudo_angle(nb_pairs)
    type (rot), intent(out) :: matrot(nb_pairs)
    type (rot), intent(out) :: matrottrans(nb_pairs)
    
    ! set all 2x2 elements matrix components
    matrot%elements(1, 1) = dcmplx(dcos(pseudo_angle/2.d0), 0.d0) 
    matrot%elements(1, 2) = dcmplx(dsin(pseudo_angle/2.d0), 0.d0) 
    matrot%elements(2, 1) = -matrot%elements(1, 2)
    matrot%elements(2, 2) = matrot%elements(1, 1)

    ! set transpose matrix
    matrottrans%elements(1, 1) = matrot%elements(1, 1)
    matrottrans%elements(1, 2) = matrot%elements(2, 1)
    matrottrans%elements(2, 1) = matrot%elements(1, 2)
    matrottrans%elements(2, 2) = matrot%elements(1, 1)

  end subroutine rotmat

  subroutine pseudo_angles (pseudo_angle)
    implicit none
    ! Local arrays
    double precision, intent(out) :: pseudo_angle(nb_pairs)

    ! Array of pseudospin angles
    pseudo_angle = datan (C12 / (polar*DJ))
    !print*, 'pseudo_angle',(pseudo_angle(i),i=1,nb_pairs)
    
  end subroutine pseudo_angles

  subroutine eigen_energies (eigen_ener)
    implicit none
    ! Local arrays
    double precision, intent(out) :: eigen_ener(nb_pairs)
    
    ! Array of eigenenergies
    eigen_ener = 0.25d0 * dsqrt(C12**2 + (polar*DJ)**2)
    !print*, 'eigen_ener',eigen_ener  

  end subroutine eigen_energies

end module deco


