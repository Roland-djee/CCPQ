!  Hamiltonian in units of 2*pi* megahertz (H=>H/hbar):
! H= B_0(mu_e S_z + mu_n I I_z) + 2pi ACO I.S
! if also have a B-field along x then:
! add H--> H + Bx (mu_eS_x+mu_n I_x)
! if have oscillating microwave field:
! add + om_mw S_z cos omega t

! external field along Z in units of Tesla
! B0 in Tesla
! mu_e=2*pi*28d3 Mega-Hertz/Tesla
! Cmus=mu_e/2*pi
! ACS=2pi*1.4754GHz,AC=1.4754GHz
! Bx =field along x in Tesla
! omega=frequency of oscillating field
! all the frequencies below must be multiplied by 2 pi 
! to turn them into angular frequencies omega=2 pi f
!OWP=B0=0.187975d0
  PARAMETER(B0=0.35d0,Bx=0.d0,AC=117.53d0,AB0=0.d0,&
	    Cmus=28.0249526619551d3,Cmun=-17.1d0,&
! moments for Si 29
            Bmus=0.d0,&
            Bmun=-8.4654991695408d0,&
 ! spin magnitudes :electronic & nuclear 
           xs=0.5d0,xi=0.5d0,&
           xsb=0.0d0,xib=0.5d0,&
! nloc= (2*xs+1)*(2*xi+1)
! here the second atom Si 29 has S=0 
! NRED=4 if we only use eg 4 states from the bismuth 20-state space
       NCS=(2*xs+1)*(2*xi+1),NB=(2*xsb+1)*(2*xib+1),NRED=4,&
!  is bath a set of simple spin 1/2? if so ibfree=1,or ibfree=0
	    nsite=2,NDIM=NRED*NB**nsite,ibfree=1,&
! npulse=number of pi pulses
! NMAX matrix storage dimension
!iEPR1 or 2= the central spin levels involved in EPR transition
	    npulse=400,NMAX=50000)

