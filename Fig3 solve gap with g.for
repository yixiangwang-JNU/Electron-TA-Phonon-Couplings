*******************************************************************
*     Solve the mean-field gap with the gradient descent method
*******************************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.1416 
      real(8), parameter :: xm=5                 ! Dirac mass
      real(8), parameter :: xi=100               ! in-plane band inversion parameter
      real(8), parameter :: xiz=200              ! z-direction band inversion parameter
      real(8), parameter :: v=393.75             ! in-plane Fermi velocity 
      real(8), parameter :: vz=0                 ! z-direction Fermi velocity      
      real(8), parameter :: xmuB=0.05788         ! Bohr magneton 
      real(8), parameter :: g1=-8                ! Zeeman parameter1
      real(8), parameter :: g2=10                ! Zeeman parameter2
      real(8), parameter :: xB=6.5               ! magnetic field
      real(8), parameter :: xlB=25.6/sqrt(xB)    ! magnetic length, in unit of nm
      real(8), parameter :: xmu=3.73             ! chemical potential 
      real(8), parameter :: dkz=1E-5             ! precision of kz
      real(8), parameter :: xLz=1/dkz            ! length in the z direction 
cc    real(8), parameter :: gg=13200             ! interaction strength 
      real(8), parameter :: gammL=200            ! step length
      real(8), parameter :: gammA=100
cc      
      common xmm,xkF1,xkF2,QL1,QL2,QT1,gg,DL1old,DL2old,DT1old
      end module value
cc
**********************
*     main program
**********************
      program main
	use value
      implicit double precision (a-h,o-z)
      character(100) filename1,filename2,filename3,filename4,
     &               filename5,filename6 
      write(filename1,'(a20)') 'DL1.dat' 
      write(filename2,'(a20)') 'resultsDL1.dat'
      write(filename3,'(a20)') 'DL2.dat'       
      write(filename4,'(a20)') 'resultsDL2.dat'
      write(filename5,'(a20)') 'DT1.dat'       
      write(filename6,'(a20)') 'resultsDT1.dat' 
      open(10,file=filename1)     
      open(20,file=filename2)    
      open(30,file=filename3)     
      open(40,file=filename4)     
      open(50,file=filename5)     
      open(60,file=filename6)  
cc
      xmm=xi/(xlB*xlB)-xm-g1*xmuB*xB/2      
      xkF1=sqrt((-xmu-xmm+g2*xmuB*xB/2)/xiz)    ! Fermi wavevectors    
      xkF2=sqrt((xmu-xmm-g2*xmuB*xB/2)/xiz)      
      QL1=2*xkF1                                ! LA phonon wavevectors       
      QL2=2*xkF2                                     
      QT1=xkF2-xkF1                             ! TA phonon wavevectors    
cc 
      DLnew0=0.5
      DTnew0=0.5
      do 100 g=1.3,20,0.1                  ! the e-ph coupling strength, in unit of eV 
      write(*,*) g 
      gg=g*g*132.231 
cc 
      DL1old=DLnew0
      do 200 i=1,100000                    ! LA1                   
         call solve_LA1(ten,tden)    
         DL1new=DL1old-gammL*tden  
         error=abs(DL1new-DL1old)
         write(10,'(3f16.8,I8)') g,DL1new,tden,i 
         if(error<1E-8) exit
         DL1old=DL1new         
200   continue   
      write(20,'(3f16.8,I8)') g,DL1new,ten,i
cc
      DL2old=DLnew0      
      do 300 i=1,200000                    ! LA2
         call solve_LA2(ten,tden)
         DL2new=DL2old-gammL*tden
         error=abs(DL2new-DL2old)
         write(30,'(3f16.8,I8)') g,DL2new,tden,i 
         if(error<1E-8) exit
         DL2old=DL2new         
300   continue   
      write(40,'(3f16.8,I8)') g,DL2new,ten,i    
cc      
      DT1old=DTnew0
      do 400 i=1,100000                    ! TA1                     
         call solve_TA1(ten,tden)    
         DT1new=DT1old-gammA*tden  
         error=abs(DT1new-DT1old)
         write(50,'(4f16.8,I8)') g,DT1new,tden,i 
         if(error<1E-8) exit 
         DT1old=DT1new         
400   continue
      write(60,'(3f16.8,I8)') g,DT1new,ten,i 
cc
100   continue 
cc
      end program main
cc
********************************************* 
*     solve the gap for e-LA1 ph coupling
********************************************* 
      subroutine solve_LA1(ten,tden)
      use value
      implicit double precision (a-h,o-z)
cc
      ten1=0                              ! total energy
      tden=0                              ! derivative of total energy
cc    
      do 100 zk=-0.5,-xkF1,dkz
      en=-xmm+g2*xmuB*xB/2-xiz*zk*zk-xiz*QL1**2/2-xiz*zk*QL1    ! varepsilon_{-kF1}
     &    -0.5*sqrt(xiz**2*(QL1**2+2*zk*QL1)**2+4*DL1old**2) 
      den=-2*DL1old/sqrt(xiz**2*(QL1**2+2*zk*QL1)**2+4*DL1old**2) 
cc
      ten1=ten1+(en-xmu)/(2*Pi*xlB**2*xLz) 
      tden=tden+den/(2*Pi*xlB**2*xLz)
100   continue
cc
      ten1=ten1*2                      ! inversion symmetry
      tden=tden*2+2*DL1old/gg          ! inversion symmetry
cc
      ten2=0                     
      do 200 zk=-xkF2,xkF2,dkz          ! varepsilon_{-kF2}--varepsilon_{kF2}
         en=xiz*zk*zk+xmm+g2*xmuB*xB/2
         ten2=ten2+(en-xmu)/(2*Pi*xlB**2*xLz) 
200   continue
      ten=ten1+ten2+DL1old**2/gg 
cc
      end subroutine solve_LA1
cc
********************************************* 
*     solve the gap for e-LA2 ph coupling
********************************************* 
      subroutine solve_LA2(ten,tden)
      use value
      implicit double precision (a-h,o-z)
cc
      ten1=0                              ! total energy
      tden=0                              ! derivative of total energy
cc       
      do 100 zk=-xkF2,0,dkz
      en=xmm+g2*xmuB*xB/2+xiz*zk*zk+xiz*QL2**2/2+xiz*zk*QL2        ! varepsilon_{-kF2}
     &    -0.5*sqrt(xiz**2*(QL2**2+2*zk*QL2)**2+4*DL2old**2) 
      den=-2*DL2old/sqrt(xiz**2*(QL2**2+2*zk*QL2)**2+4*DL2old**2) 
cc
      ten1=ten1+(en-xmu)/(2*Pi*xlB**2*xLz) 
      tden=tden+den/(2*Pi*xlB**2*xLz)
100   continue
cc
      ten1=ten1*2                    ! inversion symmetry
      tden=tden*2+2*DL2old/gg        ! inversion symmetry
cc
      ten2=0
      do 200 zk=-0.5,-xkF1,dkz       ! varepsilon_{-kF1}
         en=-xiz*zk*zk-xmm+g2*xmuB*xB/2
         ten2=ten2+(en-xmu)/(2*Pi*xlB**2*xLz) 
200   continue
      ten2=ten2*2                    ! inversion symmetry
      ten=ten1+ten2+DL2old**2/gg  
cc      
      end subroutine solve_LA2
cc       
********************************************* 
*     solve the gap for e-TA1 ph coupling
********************************************* 
      subroutine solve_TA1(ten,tden)
      use value
      implicit double precision (a-h,o-z)
cc
      ten=0                                ! total energy
      tden=0                               ! derivative of total energy
cc    
      do 100 zk=-xkF2,0,dkz                ! varepsilon_{-kF2}
         aa=xiz*(2*zk**2+2*zk*QT1+QT1**2)+2*xmm
         en=-xiz*QT1*(QT1+2*zk)/2-sqrt(aa**2+8*DT1old**2)/2+g2*xmuB*xB/2
         den=-4*DT1old/sqrt(aa**2+8*DT1old**2)
cc      
         ten=ten+(en-xmu)/(2*Pi*xlB**2*xLz) 
         tden=tden+den/(2*Pi*xlB**2*xLz)         
100   continue   
cc
      do 200 zk=-0.5,-xkF1,dkz             ! varepsilon_{-kF1}
         bb=xiz*(2*zk**2-2*zk*QT1+QT1**2)+2*xmm 
         en=xiz*QT1*(QT1-2*zk)/2-sqrt(bb**2+8*DT1old**2)/2+g2*xmuB*xB/2 
         den=-4*DT1old/sqrt(bb**2+8*DT1old**2) 
cc
         ten=ten+(en-xmu)/(2*Pi*xlB**2*xLz) 
         tden=tden+den/(2*Pi*xlB**2*xLz)
200   continue
cc
      ten=ten*2+2*DT1old**2/gg             ! inversion symmetry
      tden=tden*2+4*DT1old/gg              ! inversion symmetry 
cc      
      end subroutine solve_TA1      
      