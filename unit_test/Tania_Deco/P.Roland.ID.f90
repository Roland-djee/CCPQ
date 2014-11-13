! Code to calculate analytical T matriics for hahn sequence
! using only pair correlation
!  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Inputs: 1) adjust file 'P.Si.h'
!         2) choose time size Deltat (in microseconds) in line 78 ie we evolve for npulse* Deltat
!     so if your decay lasts 1 second and npulse=100 you want Delta=10000.
!     maybe relocate Deltat to P.Si.h if you prefer.

!    3) in line 135 input the value of J_1 XJ1=1.0013832327 for 1MHz example
! ie we use MHz units for J and C12 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	IMPLICIT NONE      
        integer  ::ll,kk,it,ii,jj,ireal,NPAIR,itime,ib,NREAL,ipair
        integer:: isite,jsite,icount,ncouple,i1,i2
        integer::ibath,index1,index2,id,iup
        integer::ncluster,NBM,NTOT,NZ,NT
        double precision::time,COJ1,COJ2,OT2,OTOT,BMAG,temp,temp1,T2
! analytics
     double precision::XJ1,XJ2,C12,XMD,XMUP,POLD,POLUP,Sigma4
! assume NCS> or =NB (dim of CS bigger or eq to Bath)

! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::npulse,NMAX,iEPR1,iEPR2
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun

        include 'P.SI.h'
        double precision:: PI,PI2
! time propagation
        double precision::DELTAT,test,TT2,T2TOT

        
! configurations for central spin and 1 bath spin
         INTEGER, DIMENSION(NRED):: ICS

      double precision, dimension(nsite,nsite)::CHYP
      double precision, dimension(nsite)::CSHF

! Cluster data and spin labels
 

   integer, dimension(10000):: NUMREAL
  double precision:: HYPHOLD(3000,8000)
   double precision:: SHFHOLD(3000,8000)
! Decoherence
         double complex::TRACEBATH
           double complex::TRACEB(npulse)
     double precision, dimension(NMAX,npulse):: TILDL
     double precision:: XL(npulse,NB**nsite)
     double precision, dimension(npulse+1):: DECAY,DECAYav
     double precision, dimension(10000,npulse+1):: DROP
 
!&&&&&&&&&&&&&&&&&&& April 2014
! this is set up for phosphorus now
! state 4 has m=+1; state 3 m=0; state 2 m=-1; state 1 m=0
! | 1/2 1/2>=4; | 1/2 -1/2>=3; | -1/2 -1/2>=2; |-1/2 1/2>=1;
! But check the numbering above!
!&&&&&&&&&&&&&&&&&&
           DATA ICS/1,2,3,4/
     NT=4
! m values for m --> m-1 transition 
! set the quantum numbers of the upper and lower state
      iEPR1=4
      iEPR2=1
      iup=iEPR1
      id=iEPR2
       write(6,*)iup,id
! Set quantum numbers of upper and lower states to get 
! the polarizations
    DECAYav=0.d0
     XMUP=1.d0
      XMD=0.d0

      write(6,*)'iup,id,m_up,m_d',iup,id,XMUP,XMD
      !stop
! set time step in microseconds. About 2000 at OWP
! since the decays are very long. About 4.d0 far 
! from OWP. Note that Hahn fails at OWP- ifyou want OWP use the
! FID code.
!              DELTAT=10000.d0
!               DELTAT=4.d0
      write(6,*)'time step deltat in mu s? use ~10000 if 1 s decay; 4. if 0.4ms'
        read(5,*) deltat
      write(6,*) 'number of realisations <20 '
      read(5,*) NPAIR
      NREAL=1
      !deltat = 10000.D0
      !NREAL = 576
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! sort out index numbers of EPR transitions 
    do ll=1,NRED
    if(ICS(ll)==iEPR1)Index1=ll
    if(ICS(ll)==iEPR2)Index2=ll
    enddo
! Zero arrays
                TILDL=0.d0
                XL=0.d0
            
 	       pi=dacos(-1.d0)
               pi2=2.d0*pi 
               ACS=AC*pi2

!************************
      OPEN(11,file='couplings.dat',status='unknown')
      OPEN(8,file='DECAY.dat',status='unknown')
      OPEN(10,file='DECAYav.dat',status='unknown')
      OPEN(12,file='T2.dat',status='unknown')
!******************************************************************
      NTOT=NRED*NB**nsite
     write(6,*)NTOT,NTOT
! Read in data on the clusters
   
!  Call ReadClusters(NREAL,NUMREAL,SHFHOLD,HYPHOLD)
 
! loop over realisations
    do ireal=1,NREAL
 
! Zero trace arrays
                TILDL=0.d0
             
 ! polarisation of the states
            Call  POLAR(iEPR1,iEPR2,NT,NTOT,BMAG,XMD,XMUP,POLD,POLUP)
! Loop over realisations
!       NPAIR=NUMREAL(ireal)
      !NPAIR=1
      ncluster=2 

       DO ipair=1,NPAIR
!       DO ipair=1,1
!      write(6,*)'realization=',ipair
! Load up the cluster coupling arrays CSHF(ncluster),CHYP(ncluster,ncluster)
          read(11,*)C12,XJ1,XJ2
          !if(i1 == i2)then
          !   C12 = 0.D0
          !end if
!          print*,i1,i2,i3,C12,XJ1
!          stop
!          XJ2=XJ1
          Sigma4=9.7d9*2*pi
 !         C12=C12+XJ1*XJ2/sigma4
 
!     C12=HYPHOLD(ipair,ireal)
!     XJ2=SHFHOLD(ipair,ireal)
!     XJ1=1.0543049487D6
!     XJ1=1.0013832327d+00
!    XJ1=5.4612417111d+00
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Scale C12 and XJ2 by 1.e-6 because Roland uses Herz units
! code is written in Mrads
     C12=1.d-6*C12
     XJ2=1.d-6*XJ2
     XJ1=1.d-6*XJ1
    write(6,*)C12,XJ2,XJ2-XJ1
    !stop
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
100 format(2I5,D14.6)

    open(15,file='pseudo_angles.dat')
    open(16,file='eigen_energies.dat')

  Call DIAGONALISER(id,iup,NTOT,NT,NBM,ICS,XJ1,XJ2,C12,POLD,POLUP,DELTAT,TRACEB,OT2)
     
!**********************************************************
! propagate each realisation in time: 
!**************************************************************
! average the XL(ibath,itime) over the initial ibath states
! and store the results in tildeL
   
   do itime=1,npulse

! average over bath states...up up/dd L=1
       TILDL(ipair,itime)=0.5+0.5d0*(TRACEB(itime))
       !TILDL(ipair,itime)=(TRACEB(itime))
       !print*,TILDL(ipair,itime)
   enddo
   !stop

! finish loop over pairs
    enddo



       test=16./22.*0.5*(abs(pold)+abs(polup))
     write(6,102)BMAG,test/abs(POLD-POLUP),POLD,POLUP
!     write(10,102)BMAG,test/abs(POLD-POLUP),POLD,POLUP
102 format(6D14.6)

      
 

!*************************************************************
! finally combine all clusters- multiply all |L(t)| together
!*************************************************************
         call DECOHERE(NPAIR,NUMREAL,TILDL,DECAY)
       time=0.d0
!        write(8,*)time,TEMP,BMAG
!        DROP(ib,1)=TEMP
       do itime=1,npulse
       time=(itime)*DeltaT*1.d-3
       DROP(ib,itime)=DECAY(itime)
       DECAYav(itime)=DECAYav(itime)+DECAY(itime)/NReal
! print out total decay
!       write(8,*)2*time,DECAY(itime)
!
       enddo


! END LOOP OVER REALIZATIONS
        enddo
! write averaged decay
           do itime=1,npulse
       time=(itime)*DeltaT*1.d-3
! double the time because Hahn is for 2*Tau
              write(10,*)time,abs(DECAYav(itime))
        enddo

        stop
        end






   Subroutine DIAGONALISER(id,iup,NTOT,NT,NBM,ICS,COJ1,COJ2,C12,POLD,POLUP,DELTAT,TRACE,OT2)
	IMPLICIT NONE
           integer::itime,NBM
 
! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::npulse,NMAX,iEPR1,iEPR2,iup,id
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun
        double precision:: DELTAT,PI,DIF
        include 'P.SI.h'
     integer::NTOT,NT,ii,jj,icn,index

       DOUBLE precision::PH,COJ1,COJ2,C12,HOLDM,HOLDP,ANORM,test
       DOUBLE PRECISION::DP1,DM1,RP1,RM1,DP2,DM2,RP2,RM2,DS1,DS2,RS1,RS2
       DOUBLE PRECISION::gcmud,gcpud,gsmud,gspud,gcm,gcp,gsm,gsp,phiup,phid
       DOUBLE precision:: COEF1,COEF2,TIME,POLD,POLUP,OT2,TEM,TEM1,TEM2,TEM3
         INTEGER, DIMENSION(NRED):: ICS
        double precision::POL(NCS)
        double precision::ENECP(NT),ENECM(NT),ENE(NT),EPS(NT)
        double precision::EIGM(NT),EIGP(NT)
         double complex ::PHASE,DPP,DMM,RPP,RMM,DPP2,DMM2,RPP2,RMM2,RSM,RSP
          DOUBLE COMPLEX:: TRACE(npulse)
       DOUBLE COMPLEX, DIMENSION(2,2)::HOLD,TU,TL,ROTU,ROTL,ROTRU,ROTRL
       DOUBLE COMPLEX, DIMENSION(2,2)::ZGU,ZGL,HOLD1,ZGUL,ZGLU
 !        analytical
       PI=dacos(-1.d0)
       POL=0.d0
       POL(id)=POLD
       POL(iup)=POLUP

       TRACE=CMPLX(0.d0,0.d0)

       ROTU=CMPLX(0.d0,0.d0)
       ROTL=CMPLX(0.d0,0.d0)
       ROTRU=CMPLX(0.d0,0.d0)
       ROTRL=CMPLX(0.d0,0.d0)
       HOLD=CMPLX(0.d0,0.d0)
       HOLD1=CMPLX(0.d0,0.d0)
       ZGU=CMPLX(0.d0,0.d0)
       ZGL=CMPLX(0.d0,0.d0)
! zero arrays
       do ii=1,NT
       ENE(ii)=0.d0
       ENECP(ii)=0.d0
       ENECM(ii)=0.d0
       EPS(ii)=0.d0

       enddo
      PHASE=CMPLX(0.d0,0.d0)

       !do ii=1,NRED
       !index=ICS(ii)
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! truncate J
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        DIF=(COJ1-COJ2)
        !if(abs(dif).lt. 1.d-2) DIF=0.d0
        !POL(index) = 0.5d0
        !HOLDM=DIF*POL(index)
        !HOLDM=DIF
        !ENE(index)=HOLDM
        !EPS(index)=sqrt(HOLDM**2+C12**2)
        !ENECP(index)=HOLDM+EPS(index)
        !ENECM(index)=HOLDM-EPS(index)
!       ! write(6,*)index,EPS(index)
        ! print*,'POL(index)',POL(index)
        ! print*,'index',index
        ! print*,'iup,id',iup,id        
      !enddo
     
      !print*,'ENE(iup),ENE(id)',ENE(iup),ENE(id)
      !stop
      !print*,'POLU',POL(iup)
      !print*,'POLD',POL(id)
      


! work out the pseudo spin angles
         phiup=atan(c12/ (0.5d0*DIF))
         phid=-phiup
         !if(ENE(iup).lt.0.d0)phiup=pi-phiup
         !if(ENE(id).lt.0.d0)phid=pi-phid

         !write(15,*)c12,DIF,phiup,phid
         !write(16,*)c12,DIF,EPS(iup),EPS(id)


         phiup=phiup/2
         phid=phid/2

! upper state : rotation matrices
         ROTU(1,1)=CMPLX(cos(phiup),0.d0)
         ROTU(1,2)=CMPLX(sin(phiup),0.d0)
         ROTU(2,1)=CMPLX(-sin(phiup),0.d0)
         ROTU(2,2)=ROTU(1,1)
         ROTRU(1,1)=ROTU(1,1)
         ROTRU(2,2)=ROTU(2,2)
         ROTRU(1,2)=-ROTU(1,2)
         ROTRU(2,1)=-ROTU(2,1)
!    lower state: rotation matrices
         ROTL(1,1)=CMPLX(cos(phid),0.d0)
         ROTL(1,2)=CMPLX(sin(phid),0.d0)
         ROTL(2,1)=CMPLX(-sin(phid),0.d0)
         ROTL(2,2)=ROTL(1,1)
         ROTRL(1,1)=ROTL(1,1)
         ROTRL(2,2)=ROTL(2,2)
         ROTRL(1,2)=-ROTL(1,2)
         ROTRL(2,1)=-ROTL(2,1)

        do jj=1,npulse

        ZGU=CMPLX(0.d0,0.d0)
        ZGL=CMPLX(0.d0,0.d0)
        TU=CMPLX(0.d0,0.d0)
        TL=CMPLX(0.d0,0.d0)

        TIME=jj*deltat
!        print*,'TIME,jj,deltat',TIME,jj,deltat

        EPS(iup) = 0.25d0 * dsqrt(C12**2 + (0.5d0*DIF)**2)
        EPS(id)  = 0.25d0 * dsqrt(C12**2 + (0.5d0*DIF)**2)

! Work out time evolution matrices
   
! ZGU and ZGL are the Z gates for upper and lower states

         phase=CMPLX(0.d0,-(EPS(iup))*TIME/4.)
         ZGU(1,1)=EXP(PHASE)
         phase=CMPLX(0.d0,(EPS(iup))*TIME/4.)
         ZGU(2,2)=EXP(PHASE)
         phase=CMPLX(0.d0,-(EPS(id))*TIME/4.)
         ZGL(1,1)=EXP(PHASE)
         phase=CMPLX(0.d0,(EPS(id))*TIME/4.)
         ZGL(2,2)=EXP(PHASE)
! work out time evolutin matrices for FID
         HOLD=CMPLX(0.d0,0.d0)

         HOLD=MATMUL(ZGU,ROTU)
         TU=MATMUL(ROTRU,HOLD)
         HOLD=CMPLX(0.d0,0.d0)
         HOLD=MATMUL(ZGL,ROTL)
         TL=MATMUL(ROTRL,HOLD)
! decoherence of Si29 by direct flip flop
!   FID
!      TRACE(jj)=abs(TU(1,1))
!  Hahn
!        TRACE(jj)=(abs(TU(1,1)))**2

! Hahn +ID 
!        TRACE(jj)=abs(TU(1,1))**2+TU(1,2)**2
!       TRACE(jj)=abs(TRACE(jj))
TRACE(jj)=abs(conjg(TL(1,1))*TU(1,1) + conjg(TL(2,1))*TU(2,1))

        enddo

102  format(I3,4D16.8)
100  format(6D14.6)
!**************************************************

      RETURN
      END

 
!**************************************************************
!***************************************************
 Subroutine Readclusters(NREAL,NUMREAL,SHFHOLD,HYPHOLD)
   IMPLICIT NONE
   integer:: ii,ll,jj,nn,ipair,isite,jsite,ncluster
   integer::ispin1,ispin2,ncount,ncouple,index,NPAIR,ireal,NREAL
! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::npulse,NMAX,iEPR1,iEPR2,Noused
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB,XX
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun
        include 'P.SI.h'

   integer, dimension(20):: NUMREAL
   double precision:: HYPHOLD(3000,8000)
   double precision:: SHFHOLD(3000,8000)

     NUMREAL=0
     ncluster=2
  
                  do ireal=1,NREAL

   
     do jj=1,NMax 
    SHFHOLD(jj,ireal)=0.d0
    HYPHOLD(jj,ireal)=0.d0
     enddo
    

!  Open(7,file='nuclear.coupling.dat',status='unknown')
  Open(7,file='Roland.coupling.dat',status='unknown')
   read(7,*) NPAIR
!  IF(nn /= ncluster) THEN
!  write(6,*) 'Problem with data file'
!  STOP
!  ENDIF
  NUMREAL(ireal)=NPAIR
        write(6,*)NPAIR,NREAL
        DO ipair=1,NPAIR
!           print*,ipair
      read(7,*)Noused,Noused,Noused,HYPHOLD(ipair,ireal),SHFHOLD(ipair,ireal)
        ENDDO

                       enddo
        RETURN
        END

!***********************************************
! combine results to work out decay in time   
   SUBROUTINE DECOHERE(NPAIR,NUMREAL,TILDL,DECAY)
   IMPLICIT NONE

! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::npulse,NMAX,iEPR1,iEPR2
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun
        include 'P.SI.h'
      integer:: NT,NTOT,NBM,ii,jj,icn,itime
       integer::ipair,NPAIR,ncluster
        integer, dimension(nsite):: NUMREAL 
          double precision::TEMP
          double precision, dimension(NMAX,npulse):: TILDL
         double precision, dimension(npulse):: DECAY
       DECAY=0.d0       
       
       do itime=1,npulse
       TEMP=1.d0
       do ipair=1,NPAIR
          TEMP=TEMP*TILDL(ipair,itime)
          !TEMP=TEMP + TILDL(ipair,itime)
       enddo
       DECAY(itime)=TEMP
       enddo
       RETURN
       END
 
 
   SUBROUTINE POLAR(iEPR1,iEPR2,NT,NTOT,BMAG,XMD,XMUP,POLD,POLUP)
   IMPLICIT NONE
        double precision::pi2,pi,XMA,XM,anorm,POLD,POLUP,XMD,XMUP
       double precision::COEF1,COEF2,AP,AM,RM,RM2,BMAG
! parameter file 
        integer:: NCS,NB,NRED,NDIM,ibfree,nsite
        integer::npulse,NMAX,iEPR1,iEPR2
! hf constants of central and bth spins
        double precision::ACS,AB,AC,AB0
! fields
        double precision:: B0,BX,OMX,FB
        double precision::xs,xi,xsb,xib
        double precision::Cmus,Cmun,Bmus,Bmun
        include 'P.SI.h'
      integer:: NT,NTOT,NBM,ii,jj,ll,itime
       double precision::Cms(NCS),Cmi(NCS)
       double precision:: POL(NCS),POL1(NCS),EVALS(NCS)

       pi=dacos(-1.d0)
       pi2=2.d0*pi
       ACS=AC*pi2
       XMA=BMAG/AC*(CMUS-CMUN)
        POLUP=0.d0
        POLD=0.d0

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
!        write(6,*)'magnetic field, Polarisation (down,up)'
!        write(6,100)BMAG,POLD,POLUP
     
100   format(3D14.6,I3)
       return
       end
!
