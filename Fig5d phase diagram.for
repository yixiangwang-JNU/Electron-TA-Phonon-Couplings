*******************************************************************
*     Solve the mean-field gap with the gradient descent method
*******************************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.1416 
      real(8), parameter :: xm=2.5               ! Dirac mass
      real(8), parameter :: xi=120               ! in-plane band inversion parameter
      real(8), parameter :: xiz=-200             ! z-direction band inversion parameter
      real(8), parameter :: v=295.313            ! in-plane Fermi velocity 
      real(8), parameter :: vz=0                 ! z-direction Fermi velocity      
      real(8), parameter :: xmuB=0.05788         ! Bohr magneton 
      real(8), parameter :: g1=-6                ! Zeeman parameter1
      real(8), parameter :: g2=10                ! Zeeman parameter2
      real(8), parameter :: dkz=1E-5             ! precision of kz
      real(8), parameter :: xLz=1/dkz            ! length in the z direction 
      real(8), parameter :: gammL=200            ! step length for LA coupling
      real(8), parameter :: gammT=100            ! step length for TA coupling 
cc      
      common xB,xlB,xmu,xkF1,xkF2,QL1,QL2,QT1,QT1b,   
     &       g,gg,xmm,deg,
     &       DL1old,DL2old,DT1old
      end module value
cc
**********************
*     main program
**********************
      program main
	use value
      implicit double precision (a-h,o-z)
      integer, parameter :: Num=241              ! number of magnetic field
      integer, parameter :: Num1=118
      character(100) filename1,filename2,filename3,
     &               filename4,filename5,filename6,
     &               filename7,filename8,filename9,filename10    
cc      
      write(filename1,'(a20)') 'mu1.dat' 
      write(filename2,'(a20)') 'mu2.dat'
      write(filename3,'(a20)') 'results.dat'     
      write(filename4,'(a20)') 'DL2.dat'
      write(filename5,'(a20)') 'DL2results.dat'      
      write(filename6,'(a20)') 'DLco.dat'
      write(filename7,'(a20)') 'DLcoresults1.dat'
      write(filename8,'(a20)') 'DLcoresults2.dat'
      write(filename9,'(a20)') 'DT1.dat' 
      write(filename10,'(a20)') 'DT1results.dat'       
      open(11,file=filename1)
      open(12,file=filename2)
      open(13,file=filename3) 
      open(21,file=filename4)
      open(22,file=filename5)
      open(31,file=filename6)
      open(32,file=filename7)
      open(33,file=filename8)
      open(41,file=filename9)
      open(42,file=filename10)      
cc
      do 100 i=1,Num1
      write(*,*) i
      read(11,*) xB,xmu                          ! magnetic field & chemical potential 
      xlB=25.6/sqrt(xB)                          ! magnetic length, in unit of nm     
      deg=2*Pi*xlB**2*xLz                        ! degeneracy factor
cc
cc    solve some variables
      xmm=xi/(xlB*xlB)-xm-g1*xmuB*xB/2
      aa=(xmu-xmm-g2*xmuB*xB/2)/xiz
      bb=(-xmu-xmm+g2*xmuB*xB/2)/xiz
      xkF1=0                                     ! Fermi wavevector for the valence band          
      if(aa>=0) xkF1=sqrt(aa)         
      xkF2=0                                     ! Fermi wavevector for the conduction band    
      if(bb>=0) xkF2=sqrt(bb) 
      QL1=2*xkF1                                 ! LA1 phonon wavevector 
      QL2=2*xkF2                                 ! LA2 phonon wavevector                  
      QT1=0                                      ! TA phonon wavevector
      if(QL1/=0 .and. QL2/=0) QT1=xkF2-xkF1   
      write(13,'(I8,100f16.8)') i,xB,xkF1,xkF2 
cc    
cc    solve the e-ph couplings      
      do g=1.5,15,0.02                           ! e-ph couling strength, in unit of eV
         gg=g*g*118.081 
         call iterate_LA2(DL2new,tenL2,iterL2) 
         write(21,'(4f16.8,I8)') xB,g,DL2new,tenL2,iterL2 
      end do
cc
100   continue      
cc      
cc      do 200 i=Num1+1,Num 
cc      write(*,*) i          
cc      read(12,*) xB,xmu                          ! magnetic field & chemical potential                               
cc      xlB=25.6/sqrt(xB)                          ! magnetic length, in unit of nm     
cc      deg=2*Pi*xlB**2*xLz                        ! degeneracy factor
cc
cc    solve some variables
cc      xmm=xi/(xlB*xlB)-xm-g1*xmuB*xB/2
cc      aa=(xmu-xmm-g2*xmuB*xB/2)/xiz     
cc      bb=(-xmu-xmm+g2*xmuB*xB/2)/xiz
cc      xkF1=0                                     ! Fermi wavevector for the valence band          
cc      if(aa>=0) xkF1=sqrt(aa)         
cc      xkF2=0                                     ! Fermi wavevector for the conduction band    
cc      if(bb>=0) xkF2=sqrt(bb) 
cc      QL1=2*xkF1                                 ! LA1 phonon wavevector 
cc      QL2=2*xkF2                                 ! LA2 phonon wavevector                  
cc      QT1=0                                      ! TA phonon wavevector
cc      if(QL1/=0 .and. QL2/=0) QT1=xkF2-xkF1   
cc      write(13,'(I8,100f16.8)') i,xB,xkF1,xkF2 
cc
cc      do 200 g=1.3,15,0.02                       ! e-ph couling strength, in unit of eV
cc         gg=g*g*118.081   
cc         call iterate_LA(DL1new,DL2new,tenL,iterL1,iterL2)  
cc         write(31,'(5f16.8,2I8)') xB,g,DL1new,DL2new,tenL,iterL1,iterL2
cc         
cc         if(QT1>0) then
cc           call iterate_TA1a(DT1new,tenT1,iterT1)
cc         else if(QT1<0) then
cc           QT1b=-QT1
cc           call iterate_TA1b(DT1new,tenT1,iterT1)
cc        else
cc           DT1new=0
cc           tenT1=0
cc         end if  
cc         write(41,'(4f16.8,I8)') xB,g,DT1new,tenT1,iterT1 
cc
cc200   continue       
cc
      end program main
cc
***************************************      
*     subroutine: iterate solve DL2
***************************************
      subroutine iterate_LA2(DL2new,tenL2,iterL2)
      use value
      implicit double precision (a-h,o-z)      
cc
      DL2old=0.5                                 ! initial value 
      do i=1,100000                  
         call solve_LA2(tenL2,tdenL2)
         DL2new=DL2old-gammL*tdenL2
         error=abs(DL2new-DL2old)
         write(22,'(4f16.8,I8)') xB,g,DL2new,tdenL2,i 
         if(error<1E-8) exit
         DL2old=DL2new         
      end do   
      iterL2=i 
cc
      end subroutine iterate_LA2
cc      
*********************************************************
*     subroutine: solve the gap for e-LA2 ph coupling
*********************************************************
      subroutine solve_LA2(tenL2,tdenL2)
      use value
      implicit double precision (a-h,o-z)
cc
      tdenL2=0                                   ! derivative of total energy
cc       
      ten1=0                                  
      do 100 zk=-xkF2,0,dkz
      en=-xmm+g2*xmuB*xB/2-xiz*zk*zk-xiz*QL2**2/2-xiz*zk*QL2        ! varepsilon_{-kF2} 
     &    -0.5*sqrt(xiz**2*(QL2**2+2*zk*QL2)**2+4*DL2old**2) 
      den=-2*DL2old/sqrt(xiz**2*(QL2**2+2*zk*QL2)**2+4*DL2old**2) 
cc
      ten1=ten1+(en-xmu)/deg
      tdenL2=tdenL2+den/deg 
100   continue
      ten1=ten1*2                                ! inversion symmetry
      tdenL2=tdenL2*2+2*DL2old/gg                ! inversion symmetry
cc
      ten2=0
      do 200 zk=-0.5,-xkF1,dkz                   ! varepsilon_{-kF1}
         en=xiz*zk*zk+xmm+g2*xmuB*xB/2
         ten2=ten2+(en-xmu)/deg 
200   continue
      ten2=ten2*2                                ! inversion symmetry
      tenL2=ten1+ten2+DL2old**2/gg  
cc      
      end subroutine solve_LA2
cc
**************************************************    
*     subroutine: iterate solve both DL1 & DL2
**************************************************
      subroutine iterate_LA(DL1new,DL2new,tenL,iter1,iter2)
      use value
      implicit double precision (a-h,o-z)
cc
      DL1old=0.1                                 ! initial value 
      do i1=1,100000
         call solve_LA1b(tenL1,tdenL1)    
         DL1new=DL1old-gammL*tdenL1  
         error=abs(DL1new-DL1old)
         write(32,'(4f16.8,I8)') xB,g,DL1new,tdenL1,i1          
         if(error<1E-8) exit
         DL1old=DL1new         
      end do
      iter1=i1
cc
      DL2old=0.1                                 ! initial value 
      do i2=1,100000                  
         call solve_LA2b(tenL2,tdenL2)
         DL2new=DL2old-gammL*tdenL2
         error=abs(DL2new-DL2old)
         write(33,'(4f16.8,I8)') xB,g,DL2new,tdenL2,i2 
         if(error<1E-8) exit
         DL2old=DL2new         
      end do  
      iter2=i2
cc      
      tenL=tenL1+tenL2     
cc      
      end subroutine iterate_LA
cc
**********************************************************************
*     subroutine: solve the gap for both e-LA1 & e-LA2 ph coupling
**********************************************************************
      subroutine solve_LA1b(tenL1,tdenL1)
      use value
      implicit double precision (a-h,o-z)
cc
      tenL1=0                               
      tdenL1=0                                   ! derivative of total energy
cc    
      do 100 zk=-0.5,-xkF1,dkz
      en=xmm+g2*xmuB*xB/2+xiz*zk*zk+xiz*QL1**2/2+xiz*zk*QL1      ! varepsilon_{-kF1}
     &   -0.5*sqrt(xiz**2*(QL1**2+2*zk*QL1)**2+4*DL1old**2) 
      den=-2*DL1old/sqrt(xiz**2*(QL1**2+2*zk*QL1)**2+4*DL1old**2) 
cc
      tenL1=tenL1+(en-xmu)/deg
      tdenL1=tdenL1+den/deg
100   continue
cc
      tenL1=tenL1*2+DL1old**2/gg                 ! inversion symmetry
      tdenL1=tdenL1*2+2*DL1old/gg                ! inversion symmetry
cc
      end subroutine solve_LA1b
cc
**********************************************************************
*     subroutine: solve the gap for both e-LA1 & e-LA2 ph coupling
**********************************************************************
      subroutine solve_LA2b(tenL2,tdenL2)
      use value
      implicit double precision (a-h,o-z)
cc
      tenL2=0                                  
      tdenL2=0                                   ! derivative of total energy
cc       
      do 100 zk=-xkF2,0,dkz
      en=-xmm+g2*xmuB*xB/2-xiz*zk*zk-xiz*QL2**2/2-xiz*zk*QL2        ! varepsilon_{-kF2} 
     &    -0.5*sqrt(xiz**2*(QL2**2+2*zk*QL2)**2+4*DL2old**2) 
      den=-2*DL2old/sqrt(xiz**2*(QL2**2+2*zk*QL2)**2+4*DL2old**2) 
cc
      tenL2=tenL2+(en-xmu)/deg
      tdenL2=tdenL2+den/deg 
100   continue
      tenL2=tenL2*2+DL2old**2/gg                 ! inversion symmetry
      tdenL2=tdenL2*2+2*DL2old/gg                ! inversion symmetry
cc
      end subroutine solve_LA2b
cc
********************************************************************
*     subroutine: iterate solve DT1 with mu above the Weyl nodes
********************************************************************
      subroutine iterate_TA1a(DT1new,tenT1,iterT1)
      use value
      implicit double precision (a-h,o-z)   
cc           
      DT1old=0.5                                 ! initial value 
      do i=1,100000                                      
         call solve_TA1a(tenT1,tdenT1)    
         DT1new=DT1old-gammT*tdenT1  
         error=abs(DT1new-DT1old)
         if(error<1E-8) exit 
         DT1old=DT1new         
      end do
      iterT1=i
cc
      end subroutine iterate_TA1a
cc
*********************************************************
*     subroutine: solve the gap for e-TA1 ph coupling
*********************************************************
      subroutine solve_TA1a(tenT1,tdenT1)
      use value
      implicit double precision (a-h,o-z)
cc
      tenT1=0                                    ! total energy
      tdenT1=0                                   ! derivative of total energy
cc    
      do 100 zk=-xkF2,0,dkz                      ! varepsilon_{-kF2}
         aa=xiz*(2*zk**2+2*zk*QT1+QT1**2)+2*xmm
         en=xiz*QT1*(QT1+2*zk)/2-sqrt(aa**2+8*DT1old**2)/2+g2*xmuB*xB/2
         den=-4*DT1old/sqrt(aa**2+8*DT1old**2)
cc      
         tenT1=tenT1+(en-xmu)/deg  
         tdenT1=tdenT1+den/deg         
100   continue   
cc
      do 200 zk=-0.5,-xkF1,dkz                   ! varepsilon_{-kF1}
         bb=xiz*(2*zk**2-2*zk*QT1+QT1**2)+2*xmm 
         en=-xiz*QT1*(QT1-2*zk)/2-sqrt(bb**2+8*DT1old**2)/2+g2*xmuB*xB/2
         den=-4*DT1old/sqrt(bb**2+8*DT1old**2) 
cc
         tenT1=tenT1+(en-xmu)/deg 
         tdenT1=tdenT1+den/deg 
200   continue
cc
      tenT1=tenT1*2+2*DT1old**2/gg               ! inversion symmetry
      tdenT1=tdenT1*2+4*DT1old/gg                ! inversion symmetry 
cc      
      end subroutine solve_TA1a
cc
********************************************************************
*     subroutine: iterate solve DT1 with mu below the Weyl nodes
********************************************************************
      subroutine iterate_TA1b(DT1new,tenT1,iterT1)
      use value
      implicit double precision (a-h,o-z)   
cc           
      DT1old=0.5                                 ! initial value 
      do i=1,100000                                      
         call solve_TA1b(tenT1,tdenT1)    
         DT1new=DT1old-gammT*tdenT1  
         error=abs(DT1new-DT1old)
         if(error<1E-8) exit 
         DT1old=DT1new         
      end do
      iterT1=i 
cc
      end subroutine iterate_TA1b
cc
*********************************************************
*     subroutine: solve the gap for e-TA1 ph coupling
*********************************************************
      subroutine solve_TA1b(tenT1,tdenT1)
      use value
      implicit double precision (a-h,o-z)
cc
      tenT1=0                                    ! total energy
      tdenT1=0                                   ! derivative of total energy
cc    
      do 100 zk=-xkF2,0,dkz                      ! varepsilon_{-kF2}
         aa=xiz*(2*zk**2-2*zk*QT1b+QT1b**2)+2*xmm
         en=-xiz*QT1b*(2*zk-QT1b)/2-sqrt(aa**2+8*DT1old**2)/2
     &      +g2*xmuB*xB/2
         den=-4*DT1old/sqrt(aa**2+8*DT1old**2)
cc      
         tenT1=tenT1+(en-xmu)/deg  
         tdenT1=tdenT1+den/deg
100   continue   
cc
      do 200 zk=-0.5,-xkF1,dkz                   ! varepsilon_{-kF1}
         bb=xiz*(2*zk**2+2*zk*QT1b+QT1b**2)+2*xmm 
         en=-xiz*QT1b*(2*zk+QT1b)/2-sqrt(bb**2+8*DT1old**2)/2
     &      +g2*xmuB*xB/2
         den=-4*DT1old/sqrt(bb**2+8*DT1old**2) 
cc
         tenT1=tenT1+(en-xmu)/deg 
         tdenT1=tdenT1+den/deg 
200   continue
cc
      tenT1=tenT1*2+2*DT1old**2/gg               ! inversion symmetry
      tdenT1=tdenT1*2+4*DT1old/gg                ! inversion symmetry 
cc      
      end subroutine solve_TA1b
cc
