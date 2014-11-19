
! DOES DD with CPMG for Hahn and EVEN pulses only
	IMPLICIT NONE  
        integer::NPULSE    
        integer  ::ll,kk,it,ii,jj,icluster,npairs,itime,ib
        integer::ibath,ilower,iup
        integer::iclusterN,NBM,NTOT
        double precision::time,CPMG,BMAG,temp,temp1,T2
! analytical pseudospin
      double precision:: XUP,ZUP,XD,ZD,OMUP,OMD,PHIUP,PHID
     double precision::XJ1,XJ2,C12,XMD,XMUP,POLD,POLUP
! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::ntau,NMAX
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun

        include 'DDPAR.h'
        double precision:: PI,PI2
! time propagation
        double precision::Del_tau

!  cluster coefficients
   integer, dimension(nsite):: NUMREAL
  double precision:: C12_STORE(nsite,20000)
   double precision:: JCOEFFS(nsite,20000)
! Decoherence
           double complex::THAHN(ntau),CPMG2(ntau),EVEN_N(ntau)
     double precision, dimension(NMAX,ntau):: TILDH,TILD2,TILDE
 
     double precision, dimension(ntau+1):: DHAHN,DCPMG2,DCPMGN
 

!*****************************************
! write the pulse number for DD
! NPULSE=0 => FID; NPULSE=1 = hahn; NPULSE=6 -> 6 CPMG 
    write(6,*)'How many CPMG pulses? 1=Hahn'
      read(5,*)NPULSE
     if(NPULSE.LT.0)write(6,*)'negative number of pulses!!'
!*****************************            
 	       pi=dacos(-1.d0)
               pi2=2.d0*pi 
               ACS=AC*pi2
! Now input timestep in milliseconds then below 
   write(6,*)'time step in tau, eg about 0.01 for 1ms T2'
   read(5,*) Del_tau
! turn into microseconds if using MRads units for couplings.
              Del_tau=Del_tau*1.d3 
!************************
      OPEN(8,file='HAHN.dat',status='unknown')
     OPEN(14,file='CPMG2.dat',status='unknown')
      OPEN(16,file='CPMGN.dat',status='unknown')
      OPEN(10,file='T2EQ9.dat',status='unknown')
      OPEN(12,file='T2.dat',status='unknown')
!******************************************************************
      NTOT=NRED*NB**nsite
     write(6,*)'NTOT=',NTOT
! Read in data on the clusters
   
  Call ReadClusters(NUMREAL,JCOEFFS,C12_STORE)

! loop over magnetic field if you want
    do ib=1,1
     BMAG=B0+(ib-1)*0.005

! Zero trace arrays
                TILDH=0.d0
                TILD2=0.d0
                TILDE=0.d0    

!**********************************************************
!DONORS IN SILICON
!***********************************************************          
! work out  POLUP,POLD=2<S_z> in upper or lower state for T_e of donors
! you need to provide m=m_s+m_I values for m --> m-1 transition and
! index numbers of transition; eg bismuth has 20 states; phosphorus has 4 states.
! Gary's OWP line was 14--> 7 in bismuth
      iup=14
      ilower=7
! m =ms+m_ifor upper and lower state
      XMUP=-1.d0
      XMD=-2.d0
     Call  POLAR_donor(iup,ilower,NTOT,BMAG,XMD,XMUP,POLD,POLUP)

! or set POLUP, POLD for T_e of  NV centres
! POLUP=\pm 2.d0, POLD=0.0 since m_s=\pm1 or m_s=0
! watch out for factor of 2 in our convention

! or POLUP =+ 1, POLD=-1  for 2<I_z> for T_n of  Si29 nuclei
POLUP = 1.d0
POLD  =-1.d0
!!****************************************

! Loop over clusters
       npairs=NUMREAL(2)
! size of cluster=2 for pseudospins
        iclusterN=2 
       DO icluster=1,NPAIRS
! Note that below we assume couplings in Mrads. Convert (eg *2pi if necessary
     C12=C12_STORE(1,icluster) * 1.d-6
     XJ1=JCOEFFS(1,icluster) * 1.d-6
     XJ2=JCOEFFS(2,icluster) * 1.d-6
     
100 format(2I5,D14.6)

! work out pseudospin tilt angles and precession frequencies
! 
! Donors in silicon case T_e
!********************************************** 
         ZD =(XJ1-XJ2)*POLD
         ZUP=(XJ1-XJ2)*POLUP
! if wanted correct for electron mediated coupling
!       omb=BMAG*Cmus*2*pi
!     C12=C12+XJ1*XJ2/omb
         XUP=C12
         XD=C12
!****************************************
! can add options for NV centres T_e
! POLUP=\pm 2.d0, POLD=0.0 since m_s=\pm1 or m_s=0
! watch out for factor of 2 in our convention 
 
!   ZD=...
!   ZUP=...
!   XD= 0.do for m=0 state
!   XUP=...

! Si29 nuclei T_n
!************************************
! should add hyperfine mediated correction of Yao sham Liu PRB 74 2006
! if JA,JB1,JB2 is large
! = J1J2/4 sigma_e  where sigma e is Zeeman energy of electron 9 GHz*2pi
! for symmetric sites XJB1=XJB2=XJB
!     Sigma4=9.7d9*2*pi
!     C12=C12+XJB1*XJB2/sigma4
!     CB1=CB1+XJB1*XJA/sigma4
!     CB2=CB2+XJB2*XJA/sigma4

! Scale C12 and XJ2 by 1.e-6 because Roland uses rads UNITS 
! code is written in Mrads units
!     C12=1.d-6*C12
!     CB1=1.d-6*CB1
!     CB2=1.d-6*CB2
!     POLD=-1.d0
!     POLUP=1.d0
!         ZD =(CB1-CB2)*POLD
!         ZUP=(CB1-CB2)*POLUP
!         XUP=C12
!         XD=C12
!*********************************************************
! from here on, the code is system independent
!  work out pseudospin angles and frequencies 
  call PSEUDOSPIN(XUP,ZUP,XD,ZD,OMUP,OMD,PHIUP,PHID)

 Call Propagator(NPULSE,OMUP,OMD,PHIUP,PHID,Del_tau,THAHN,CPMG2,EVEN_N)

!**********************************************************
! propagate each realisation in time: 
!**************************************************************
! average over the initial ibath states
! and store the results in tildeL

! HAHN   
   do itime=1,ntau
! average over 4 bath states...for the up up or d d, then L=1
       TILDH(icluster,itime)=0.5+0.5d0*abs(THAHN(itime))
       TILD2(icluster,itime)=0.5+0.5d0*abs(CPMG2(itime))
       TILDE(icluster,itime)=0.5+0.5d0*abs(EVEN_N(itime))
! end loop over time
   enddo
! finish loop over clusters
    enddo

!*************************************************************
! finally combine all clusters
!*************************************************************
               call DECOHERE(iclusterN,NUMREAL,TILDH,DHAHN)
               call DECOHERE(iclusterN,NUMREAL,TILD2,DCPMG2)
              call DECOHERE(iclusterN,NUMREAL,TILDE,DCPMGN)
! Write out decays

!  Hahn
       TEMP=1.d0
       time=0.d0
        write(8,102)time,TEMP,BMAG

  do itime=1,ntau
      time=2*itime*Del_tau*1.d-3
  write(8,102)time,DHAHN(itime),BMAG     
 enddo
!  CPMG2
       TEMP=1.d0
       time=0.d0
        write(14,102)time,TEMP,BMAG
  do itime=1,ntau
      time=4*itime*Del_tau*1.d-3
  write(14,102)time,DCPMG2(itime),BMAG     
 enddo

! write EVEN N
       TEMP=1.d0
       time=0.d0
        write(16,102)time,TEMP,BMAG
  do itime=1,ntau
! scale up time depending on number of Tmatrices used
! eg N=6 [T_ul T_lu]^NPULSE/2 and each T propagates for 2tau
     CPMG=2*NPULSE
      time=CPMG*itime*Del_tau*1.d-3
  write(16,102)time,DCPMGN(itime),BMAG     
 enddo



! END LOOP OVER MAGNETIC FIELDS
        enddo
102 format(6E14.6)

        stop
        end           
!************************************************************
Subroutine Propagator(NPULSE,OMUP,OMD,PHIUP,PHID,Del_tau,THAHN,CPMG2,EVEN_N)
	IMPLICIT NONE
           integer::itime,NBM 
! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::ntau,NMAX
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun
        double precision:: Del_tau
        include 'DDPAR.h'
     integer::ii,jj,kk,N,NPULSE,i2N
       DOUBLE precision::omut,omdt,a0,ax,ay,az,csu,csd,t1,t2,t3
     DOUBLE precision::A20,A2XU,A2XD,A2ZU,A2ZD,Fphase,thetu,thetd,thetdif
       DOUBLE precision:: TIME,OMUP,OMD,phiup,phid
         double complex ::term1,term2
          DOUBLE COMPLEX, dimension(ntau):: THAHN,CPMG2,EVEN_N
! zero array
       THAHN=CMPLX(0.d0,0.d0)
       CPMG2=CMPLX(0.d0,0.d0)
       EVEN_N=CMPLX(0.d0,0.d0)
       i2N=NPULSE
! begin time evolution for different pulse delays (tau=time)         
        do jj=1,ntau
        TIME=jj*Del_tau
         omut=omup*time
         omdt=omd*time
         a0=cos(omut)*cos(omdt)-sin(omut)*sin(omdt)*cos(phiup-phid)
         csu=cos(omut)*sin(omdt)
         csd=cos(omdt)*sin(omut)
         ax=csd*sin(phiup)+csu*sin(phid)
         ay=-(sin(omut))*sin(omdt)*sin(phiup-phid)
         az=csd*cos(phiup)+csu*cos(phid)
        
! hahn        
      term2=CMPLX(0.d0,2*ax*ay)
       term1= CMPLX((1.d0-2*ay*ay),0.d0)
      THAHN(jj)=term1+term2 
      THAHN(jj)=abs(THAHN(jj))
! CPMG2   
     t1=1.d0-8.*ay*ay*(ax**2+az**2)
     t2=4.d0*ax*ay*(1.d0-2.d0*(ax**2+az**2))
     term1=CMPLX(0.d0,t2)
     term2=CMPLX(t1,0.d0)
      CPMG2(jj)=term1+term2
      CPMG2(jj)=abs(CPMG2(jj))

! even pulse CPMG
       A20=a0**2+ay**2-ax**2-az**2
       A2XU=2.d0*(ax*a0-ay*az)
       A2XD=2.d0*(ax*a0+ay*az)      
       A2ZU=2.d0*(a0*az+ax*ay)
       A2ZD=2.d0*(a0*az-ax*ay)
      t1=sqrt(1.d0-A20**2) 
      Fphase=acos(A20)
      t2=A2ZU**2+A2XU**2
      t1=A2ZU/sqrt(t2)
      thetu=ACOS(t1)

      t2=A2ZD**2+A2XD**2
      t1=A2ZD/sqrt(t2)
      thetd=ACOS(t1)
     thetdif=(thetu-thetd)*0.5d0
     t1=(cos(thetdif))**2 
     t2=(sin(thetdif))**2*cos(i2N*Fphase)
      t3=sin(i2N*Fphase)*sin(thetdif)*sin(0.5*(thetu+thetd))

      term1=CMPLX(t1+t2,0.d0)
     term2=CMPLX(0.d0,t3)
      EVEN_N(jj)=term1+term2   
!       EVEN_N(jj)=1.d0-2.d0*(sin(thetdif))**2*sin(i2N*Fphase/2.d0)**2
    EVEN_N(jj)=abs(EVEN_N(jj))
        enddo

102  format(I3,4D16.8)
100  format(6D14.6)
      RETURN
      END
!*********************************************************
 
!*************************************************************
! work out pseudospin tilt angles and precession frequencies
!
  Subroutine PSEUDOSPIN(XUP,ZUP,XD,ZD,OMUP,OMD,PHIUP,PHID)
	IMPLICIT NONE
           integer::itime,NBM
 
! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::ntau,NMAX,iEPR1,iEPR2,iup,id
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun
        double precision:: Del_tau
        include 'DDPAR.h'
       DOUBLE precision, intent(in)::XUP,ZUP,XD,ZD
       DOUBLE precision, intent(out)::OMUP,OMD,PHIUP,PHID

       DOUBLE PRECISION::DELJD,DELJUP,EPSU,EPSL,PI
       PI=dacos(-1.d0)      
        OMUP=sqrt(ZUP**2+XUP**2)
        OMD=sqrt(ZD**2+XD**2)
! note that real eigenvalues of the flip flop states are 
! eigenu= (deljup \pm epsu)/4
!eigenl = ( deljd \pm epsl)/4
       OMUP=OMUP/4.d0
       OMD=OMD/4.d0
!! NB!!! that above we have 1/4 as per old convention while in the
! papers we have 1/2. This is because both C12 and POLup/POLD are 
! rescaled by 2 as requested by John/Gary for neatness.
! John preferred POLUP =<S_z>  and C12--> C12/2
! ANGLES
         phiup=atan(abs(XUP)/abs(ZUP))
         phid=atan(abs(XD)/abs(ZD))
         if(ZUP.lt.0.d0)phiup=pi-phiup
         if(ZD.lt.0.d0)phid=pi-phid

! if you generate rotation matrices halve the angle!!!
!         phiup=phiup/2
!         phid=phid/2
       RETURN
       END

 
!**************************************************************
!***************************************************
 Subroutine Readclusters(NUMREAL,JCOEFFS,C12_STORE)
   IMPLICIT NONE
   integer:: ii,ll,jj,nn,icluster,isite,jsite,iclusterN
   integer::ispin1,ispin2,ncount,ncouple,index,npairs
! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::ntau,NMAX,iEPR1,iEPR2
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun
        include 'DDPAR.h'

   integer, dimension(nsite):: NUMREAL
   double precision:: C12_STORE(nsite,20000)
   double precision:: JCOEFFS(nsite,20000)

     NUMREAL=0
     iclusterN=2
 
     do ii=1,nsite
     do jj=1,NMax 
    JCOEFFS(ii,jj)=0.d0
    C12_STORE(ii,jj)=0.d0
     enddo
     enddo

  Open(7,file='coupling.dat',status='unknown')
   read(7,*) npairs
!  IF(nn /= iclusterN) THEN
!  write(6,*) 'Problem with data file'
!  STOP
!  ENDIF
  NUMREAL(iclusterN)=npairs
    write(6,*)iclusterN,npairs,npairs
        DO icluster=1,NPAIRS
      read(7,*)JCOEFFS(1,icluster),JCOEFFS(2,icluster),C12_STORE(1,icluster)

        ENDDO
        RETURN
        END

!***********************************************
! combine results to work out decay in time   
   SUBROUTINE DECOHERE(iclusterN,NUMREAL,TILDL,DECAY)
   IMPLICIT NONE

! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::ntau,NMAX,iEPR1,iEPR2
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun
        include 'DDPAR.h'
      integer:: NT,NTOT,NBM,ii,jj,icn,itime
       integer::icluster,npairs,iclusterN
        integer, dimension(nsite):: NUMREAL 
          double precision::TEMP
          double precision, dimension(NMAX,ntau):: TILDL
         double precision, dimension(ntau):: DECAY
       DECAY=0.d0       
       npairs=NUMREAL(iclusterN)
       do itime=1,ntau
       TEMP=1.d0
       do icluster=1,npairs 
       TEMP=TEMP*TILDL(icluster,itime)
       enddo
       DECAY(itime)=TEMP
       enddo
       RETURN
       END
 
 
   SUBROUTINE POLAR_donor(iEPR1,iEPR2,NTOT,BMAG,XMD,XMUP,POLD,POLUP)
   IMPLICIT NONE
        double precision::pi2,pi,XMA,XM,anorm,POLD,POLUP,XMD,XMUP
       double precision::COEF1,COEF2,AP,AM,RM,RM2,BMAG
! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::ntau,NMAX,iEPR1,iEPR2
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun
        include 'DDPAR.h'
      integer:: NT,NTOT,NBM,ii,jj,ll,itime
       double precision::Cms(NCS),Cmi(NCS)
       double precision:: POL(NCS),POL1(NCS),EVALS(NCS)

       pi=dacos(-1.d0)
       pi2=2.d0*pi
       ACS=AC*pi2
       XMA=BMAG/AC*(CMUS-CMUN)
        POLUP=0.d0
        POLD=0.d0
        
!        XM=-3.d0
        XM=XMUP
        RM2=(XM+XMA)**2+(xi+0.5d0)**2-XM**2
        RM=SQRT(RM2)
        coef1=(RM+ (XM+XMA))/2./RM
        coef1=sqrt(coef1)
        anorm=2*RM*(RM+XM+XMA)
        COEF2=sqrt((xi+0.5d0)**2-XM**2)/sqrt(anorm)
        POLUP=coef1**2-coef2**2
! if dominant state is m=-1/2 then reverse polarity
       jj=ncs/2+1
         IF(iEPR1 <jj) POLUP=-POLUP
!        XM=-4.d0
         XM=XMD
        RM2=(XM+XMA)**2+(xi+0.5d0)**2-XM**2
        RM=SQRT(RM2)
        coef1=(RM+ (XM+XMA))/2./RM
        coef1=sqrt(coef1)
        anorm=2*RM*(RM+XM+XMA)
        COEF2=sqrt((xi+0.5d0)**2-XM**2)/sqrt(anorm)
        POLD=coef1**2-coef2**2
! if dominant state is m=-1/2 then reverse polarity
       jj=ncs/2+1
         IF(iEPR2 <jj) POLD=-POLD
        
!        write(6,100)BMAG,POLD,POLUP
     
100   format(3D14.6,I3)
       return
       end
!
