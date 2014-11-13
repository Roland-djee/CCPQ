module deco
  use types
  use constants
  implicit none

contains

  subroutine decohere
    implicit none

    ! set the polarisation of the qubit under study
    select case (Qubittype)
    case ("electron")
       polar_up   = 0.5d0
       polar_down = -0.5d0
    case ("nuclear")
       polar_up   = 0.5d0
       polar_down = -0.5d0
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
 
    double precision :: eigen_ener_up(nb_pairs)
    double precision :: eigen_ener_down(nb_pairs)
    double precision :: pseudo_angle_up(nb_pairs)
    double precision :: pseudo_angle_down(nb_pairs)
   
    ! Compute the eigenenergies for all pairs
    call eigen_energies (eigen_ener_up, eigen_ener_down)
    
    ! Compute the pseudospin angles for all pairs in up/down qubit states
    call pseudo_angles (pseudo_angle_up, pseudo_angle_down)

    select case (Dynadeco)
    case ("FID")
       call FID (eigen_ener_up, eigen_ener_down, pseudo_angle_up, &
                 pseudo_angle_down)
    case ("Hahn")
       CP_seq = 1
       call Hahn (eigen_ener_up, eigen_ener_down, pseudo_angle_up, &
                 pseudo_angle_down, CP_seq)
    case ("CP")
       call Hahn (eigen_ener_up, eigen_ener_down, pseudo_angle_up, &
                 pseudo_angle_down, CP_seq)
    case default
       write(*,*)"No valid Dynamical decoupling type in Dynamics.inp" 
       stop
    end select
   
  end subroutine IFF
  
  subroutine Hahn (eigen_ener_up, eigen_ener_down, pseudo_angle_up, &
                   pseudo_angle_down, CP_seq)
    implicit none
    integer, intent(in) :: CP_seq
    double precision, intent(in) :: eigen_ener_up(nb_pairs)
    double precision, intent(in) :: eigen_ener_down(nb_pairs)
    double precision, intent(in) :: pseudo_angle_up(nb_pairs)
    double precision, intent(in) :: pseudo_angle_down(nb_pairs)
    ! Local variables and arrays
    integer :: i, j, k
    double precision :: t, dt, L
    character(len=*), parameter :: fmt="(es20.10e3, es20.10e3)"
    character(len=100) :: filename

    type (rot) :: matrot_u(nb_pairs), matrottrans_u(nb_pairs)
    type (rot) :: matrot_d(nb_pairs), matrottrans_d(nb_pairs)
    type (rot) :: Zgate_u(nb_pairs), Zgate_d(nb_pairs)
    type (rot) :: Tu(nb_pairs), Td(nb_pairs)
    type (rot) :: Tud(nb_pairs), Tdu(nb_pairs)
    type (rot) :: Tfu(nb_pairs), Tfd(nb_pairs)

    double precision :: L_pairs(nb_pairs)
    
    ! Initialize rotation matrices for all pairs in up/down states
    call rotmat (pseudo_angle_up, matrot_u, matrottrans_u)
    call rotmat (pseudo_angle_down, matrot_d, matrottrans_d)

    if (Qubittype == "electron") then
       if (Dynadeco == "Hahn") then
          write(filename, '(a)') "Averaged_decay_electron_IFF_Hahn.dat"
       else if (Dynadeco == "CP") then
          write(filename, '(a, i1, a)') "Averaged_decay_electron_IFF_CP",&
               CP_seq,".dat"
       end if
    else if (Qubittype == "nuclear") then
       if (Dynadeco == "Hahn") then
          write(filename, '(a, es10.3e3, a)') &
               "Averaged_decay_nuclear_IFF_J=",Jnuc,"_rad_Hahn.dat"
       else if (Dynadeco == "CP") then
          write(filename, '(a, es10.3e3, a, i1, a)') &
               "Averaged_decay_nuclear_IFF_J=",Jnuc,"_rad_CP",CP_seq,".dat"
       end if
    end if
    open(16, file=filename)

    ! Loop over time window, t correponds now to 2 * CP_seq
    ! CP_seq being the number of pi-pulses in the Carr-Purcell sequence

    dt = Tmax / dble(nb_pts_t)
    t  = - dt
    
    do j=1,nb_pts_t + 1
       t = t + dt
       ! Compute Zgates for all pairs in up/down states
       call Z_gate (eigen_ener_up, t / dble(2 * CP_seq), Zgate_u)
       call Z_gate (eigen_ener_down, t / dble(2 * CP_seq), Zgate_d)

       ! Initialize identity matrices
       Tfu%elements(1,1) = dcmplx(1.d0, 0.d0)
       Tfu%elements(1,2) = dcmplx(0.d0, 0.d0)
       Tfu%elements(2,1) = dcmplx(0.d0, 0.d0)
       Tfu%elements(2,2) = dcmplx(1.d0, 0.d0)
       Tfd%elements(1,1) = dcmplx(1.d0, 0.d0)
       Tfd%elements(1,2) = dcmplx(0.d0, 0.d0)
       Tfd%elements(2,1) = dcmplx(0.d0, 0.d0)
       Tfd%elements(2,2) = dcmplx(1.d0, 0.d0)

       ! work out transition matrices in up/down states
       ! and compute the decoherence
       do i=1,nb_pairs
          ! rotate to eigenbasis and propagate
          Tu(i)%elements = matmul(Zgate_u(i)%elements, matrot_u(i)%elements)
          ! rotate back to bath basis
          Tu(i)%elements = matmul(matrottrans_u(i)%elements, Tu(i)%elements)
          ! Idem down state
          Td(i)%elements = matmul(Zgate_d(i)%elements, matrot_d(i)%elements)
          Td(i)%elements = matmul(matrottrans_d(i)%elements, Td(i)%elements)

          ! Initialize Tud and Tdu
          Tud(i)%elements = matmul(Tu(i)%elements, Td(i)%elements)
          Tdu(i)%elements = matmul(Td(i)%elements, Tu(i)%elements)

         do k=1,CP_seq
             if(modulo(k, 2) .ne. 0) then
                Tfu(i)%elements = matmul(Tfu(i)%elements, Tud(i)%elements)
                Tfd(i)%elements = matmul(Tfd(i)%elements, Tdu(i)%elements)
             else if(modulo(k, 2) == 0) then
                Tfu(i)%elements = matmul(Tfu(i)%elements, Tdu(i)%elements)
                Tfd(i)%elements = matmul(Tfd(i)%elements, Tud(i)%elements)
             end if
          end do

          ! Decoherence from initial |down-up> bath state
          L_pairs(i) = abs(conjg(Tfd(i)%elements(1, 1))*Tfu(i)%elements(1, 1) &
               - Tfd(i)%elements(1, 2) * Tfu(i)%elements(2, 1))
          
          ! average over the bath states
          L_pairs(i) = 0.5d0 + 0.5d0 * L_pairs(i)
       end do

       ! Final decay as the product over all pair decays
       L = product(L_pairs) 
       
       ! write the output
       write(16, fmt)t, L

    end do
    close(16)

  end subroutine Hahn

  subroutine FID (eigen_ener_up, eigen_ener_down, pseudo_angle_up, &
                  pseudo_angle_down)
    implicit none
    double precision, intent(in) :: eigen_ener_up(nb_pairs)
    double precision, intent(in) :: eigen_ener_down(nb_pairs)
    double precision, intent(in) :: pseudo_angle_up(nb_pairs)
    double precision, intent(in) :: pseudo_angle_down(nb_pairs)
    ! Local variables and arrays
    integer :: i
    double precision :: dt, t, L_av
    double precision :: a_u(nb_pairs), b_u(nb_pairs)
    double precision :: c_u(nb_pairs), d_u(nb_pairs)
    double precision :: a_d(nb_pairs), b_d(nb_pairs)
    double precision :: c_d(nb_pairs), d_d(nb_pairs)
    double complex   :: z_u(nb_pairs), z_d(nb_pairs)
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

    ! Initialize matrices
    a_u = dcos(pseudo_angle_up * 0.5d0)
    b_u = dsin(pseudo_angle_up * 0.5d0)
    a_d = dcos(pseudo_angle_down * 0.5d0)
    b_d = dsin(pseudo_angle_down * 0.5d0)

    do i=1,nb_pts_t + 1
       t = t + dt
       c_u = dcos(eigen_ener_up * t)
       d_u = dsin(eigen_ener_up * t)
       c_d = dcos(eigen_ener_down * t)
       d_d = dsin(eigen_ener_down * t)

       z_u = dcmplx(c_u, d_u * (b_u**2 - a_u**2))
       z_d = dcmplx(c_d, - d_d * (b_d**2 - a_d**2))

       ! Decoherence
       L = abs(z_d * z_u  + 4.d0 * a_u * a_d * b_u * b_d * d_u * d_d)
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
    Zgate%elements(1, 1) = exp(- phase * t)
    Zgate%elements(1, 2) = dcmplx(0.d0, 0.d0)
    Zgate%elements(2, 1) = dcmplx(0.d0, 0.d0)
    Zgate%elements(2, 2) = exp(phase * t)
        
  end subroutine Z_gate

  subroutine rotmat (pseudo_angle, matrot, matrottrans)
    implicit none
    double precision, intent(in) :: pseudo_angle(nb_pairs)
    type (rot), intent(out) :: matrot(nb_pairs)
    type (rot), intent(out) :: matrottrans(nb_pairs)
    
    ! set all 2x2 elements matrix components
    matrot%elements(1, 1) = dcmplx(dcos(pseudo_angle * 0.5d0), 0.d0) 
    matrot%elements(1, 2) = dcmplx(dsin(pseudo_angle * 0.5d0), 0.d0) 
    matrot%elements(2, 1) = -matrot%elements(1, 2)
    matrot%elements(2, 2) = matrot%elements(1, 1)

    ! set transpose matrix
    matrottrans%elements(1, 1) = matrot%elements(1, 1)
    matrottrans%elements(1, 2) = matrot%elements(2, 1)
    matrottrans%elements(2, 1) = matrot%elements(1, 2)
    matrottrans%elements(2, 2) = matrot%elements(1, 1)

  end subroutine rotmat

  subroutine pseudo_angles (pseudo_angle_up, pseudo_angle_down)
    implicit none
    ! Local arrays
    double precision, intent(out) :: pseudo_angle_up(nb_pairs)
    double precision, intent(out) :: pseudo_angle_down(nb_pairs)

    ! Array of pseudospin angles
    pseudo_angle_up   = datan (abs(C12 / (polar_up * DJ)))
    pseudo_angle_down = datan (abs(C12 / (polar_down * DJ)))
    where (C12 / (polar_up * DJ) .lt. 0.d0) &
         pseudo_angle_up = pi - pseudo_angle_up
    where (C12 / (polar_down * DJ) .lt. 0.d0) &
         pseudo_angle_down = pi - pseudo_angle_down

  end subroutine pseudo_angles

  subroutine eigen_energies (eigen_ener_up, eigen_ener_down)
    implicit none
    ! Local arrays
    double precision, intent(out) :: eigen_ener_up(nb_pairs)
    double precision, intent(out) :: eigen_ener_down(nb_pairs)
    
    ! Array of eigenenergies
    eigen_ener_up   = 0.25d0 * dsqrt(C12**2 + (polar_up * DJ)**2)
    eigen_ener_down = 0.25d0 * dsqrt(C12**2 + (polar_down * DJ)**2)

  end subroutine eigen_energies

end module deco


