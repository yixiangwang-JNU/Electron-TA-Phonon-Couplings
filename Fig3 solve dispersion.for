*********************************************************
*     Solve the renormalized bands with e-ph coupling  
*********************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.14159 
      real(8), parameter :: xm=5                 ! Dirac mass
      real(8), parameter :: xi=100               ! in-plane band inversion parameter
      real(8), parameter :: xiz=200              ! z-direction band inversion parameter
      real(8), parameter :: v=393.75             ! in-plane Fermi velocity 
      real(8), parameter :: vz=0                 ! z-direction Fermi velocity      
      real(8), parameter :: xmuB=0.05788         ! Bohr magneton 
      real(8), parameter :: g1=-8                ! Zeeman parameter1
      real(8), parameter :: g2=10                ! Zeeman parameter2
      real(8), parameter :: xB=9                 ! magnetic field
      real(8), parameter :: xlB=25.6/sqrt(xB)    ! magnetic length 
      real(8), parameter :: xmu=2.11             ! chemical potential 
      real(8), parameter :: DL1=0.3              ! Delta_L1  
      real(8), parameter :: DL2=0.3
      real(8), parameter :: DT1=0.212
cc
      end module value
cc
**********************
*     main program
**********************
      program main
	use value
      implicit double precision (a-h,o-z)
      character(100) filename1,filename2,filename3,filename4,filename5
      write(filename1,'(a20)') 'DL1a.dat'
      write(filename2,'(a20)') 'DL1b.dat'
      write(filename3,'(a20)') 'DL2a.dat'
      write(filename4,'(a20)') 'DL2b.dat'
      write(filename5,'(a20)') 'DT1.dat'   
      open(10,file=filename1)     
      open(20,file=filename2) 
      open(30,file=filename3)     
      open(40,file=filename4) 
      open(50,file=filename5)     
cc        
      call solve_en_LA1() 
      call solve_en_LA2() 
      call solve_en_TA1()       
cc 
      end program main
cc
****************************************************
*     solve the dispersion for LA1 e-ph coupling 
****************************************************
      subroutine solve_en_LA1()
      use value
      implicit double precision (a-h,o-z)
cc    
      xmm=xi/(xlB*xlB)-xm-g1*xmuB*xB/2 
      xkF1=sqrt((-xmu-xmm+g2*xmuB*xB/2)/xiz)        ! Fermi wavevectors   
      QL1=2*xkF1                                    ! LA phonon wavevectors 
      vF1=2*xiz*xkF1                                ! Fermi velocities 
cc
      do zk=-0.3,0.3,1E-4    
         en=xiz*zk*zk+xmm+g2*xmuB*xB/2
         write(10,'(100f16.8)') zk,en
      end do       
cc       
      do 100 zk=-0.3,0,1E-4
      if(zk+xkF1>0) then
      en1=-xmm+g2*xmuB*xB/2-xiz*zk*zk-xiz*QL1**2/2-xiz*zk*QL1 
     &    +0.5*sqrt(xiz**2*(QL1**2+2*zk*QL1)**2+4*DL1**2) 
      write(20,'(100f16.8)') zk,en1
      end if
cc  
      if(zk+xkF1<0) then
      en1=-xmm+g2*xmuB*xB/2-xiz*zk*zk-xiz*QL1**2/2-xiz*zk*QL1 
     &    -0.5*sqrt(xiz**2*(QL1**2+2*zk*QL1)**2+4*DL1**2) 
      write(20,'(100f16.8)') zk,en1
      end if
100   continue
cc      
      do 200 zk=1E-4,0.3,1E-4  
      if(zk-xkF1>0) then
      en2=-xmm+g2*xmuB*xB/2-xiz*zk*zk-xiz*QL1**2/2+xiz*zk*QL1 
     &    -0.5*sqrt(xiz**2*(QL1**2-2*zk*QL1)**2+4*DL1**2) 
      write(20,'(100f16.8)') zk,en2  
      end if
cc
      if(zk-xkF1<0) then
      en2=-xmm+g2*xmuB*xB/2-xiz*zk*zk-xiz*QL1**2/2+xiz*zk*QL1 
     &    +0.5*sqrt(xiz**2*(QL1**2-2*zk*QL1)**2+4*DL1**2)   
      write(20,'(100f16.8)') zk,en2 
      end if 
200   continue       
cc
      end subroutine solve_en_LA1
cc
****************************************************
*     solve the dispersion for LA2 e-ph coupling 
****************************************************
      subroutine solve_en_LA2()
      use value
      implicit double precision (a-h,o-z)
cc    
      xmm=xi/(xlB*xlB)-xm-g1*xmuB*xB/2        
      xkF2=sqrt((xmu-xmm-g2*xmuB*xB/2)/xiz)      ! Fermi wavevectors                  
      QL2=2*xkF2                                 ! LA phonon wavevectors                          
cc      
      do zk=-0.3,0.3,1E-4    
         en=-(xiz*zk*zk+xmm)+g2*xmuB*xB/2
         write(30,'(100f16.8)') zk,en
      end do 
cc          
      do 100 zk=-0.3,0,1E-4
      if(zk+xkF2>0) then
      en1=xmm+g2*xmuB*xB/2+xiz*zk*zk+xiz*QL2**2/2+xiz*zk*QL2 
     &    -0.5*sqrt(xiz**2*(QL2**2+2*zk*QL2)**2+4*DL2**2) 
      write(40,'(100f16.8)') zk,en1 
      end if
cc  
      if(zk+xkF2<0) then
      en2=xmm+g2*xmuB*xB/2+xiz*zk*zk+xiz*QL2**2/2+xiz*zk*QL2
     &    +0.5*sqrt(xiz**2*(QL2**2+2*zk*QL2)**2+4*DL2**2) 
      write(40,'(100f16.8)') zk,en2
      end if
100   continue
cc      
      do 200 zk=1E-4,0.3,1E-4 
      if(zk-xkF2>=0) then
      en1=xmm+g2*xmuB*xB/2+xiz*zk*zk+xiz*QL2**2/2-xiz*zk*QL2 
     &    +0.5*sqrt(xiz**2*(QL2**2-2*zk*QL2)**2+4*DL2**2) 
      write(40,'(100f16.8)') zk,en1  
      end if
cc
      if(zk-xkF2<0) then
      en2=xmm+g2*xmuB*xB/2+xiz*zk*zk+xiz*QL2**2/2-xiz*zk*QL2 
     &    -0.5*sqrt(xiz**2*(QL2**2-2*zk*QL2)**2+4*DL2**2) 
      write(40,'(100f16.8)') zk,en2 
      end if 
200   continue       
cc
      end subroutine solve_en_LA2
cc      
****************************************************
*     solve the dispersion for TA1 e-ph coupling 
****************************************************
      subroutine solve_en_TA1()
      use value
      implicit double precision (a-h,o-z)
cc    
      xmm=xi/(xlB*xlB)-xm-g1*xmuB*xB/2 
      xkF1=sqrt((-xmu-xmm+g2*xmuB*xB/2)/xiz)           ! Fermi wavevectors    
      xkF2=sqrt((xmu-xmm-g2*xmuB*xB/2)/xiz)      
      QT1=xkF1-xkF2                                    ! TA phonon wavevectors  
cc       
      do 100 zk=-0.3,0,1E-4  
      aa=xiz*(2*zk**2-2*zk*QT1+QT1**2)+2*xmm
      if(zk+xkF2>0) then                               ! varepsilon_{-kF2} 
         en1=xiz*QT1*(2*zk-QT1)/2-sqrt(aa**2+8*DT1**2)/2+g2*xmuB*xB/2
      else
         en1=xiz*QT1*(2*zk-QT1)/2+sqrt(aa**2+8*DT1**2)/2+g2*xmuB*xB/2 
      end if
cc  
      bb=xiz*(2*zk**2+2*zk*QT1+QT1**2)+2*xmm
      if(zk+xkF1>0) then                               ! varepsilon_{-kF1}   
         en2=xiz*QT1*(2*zk+QT1)/2+sqrt(bb**2+8*DT1**2)/2+g2*xmuB*xB/2  
      else
         en2=xiz*QT1*(2*zk+QT1)/2-sqrt(bb**2+8*DT1**2)/2+g2*xmuB*xB/2 
      end if
      write(50,'(100f16.8)') zk,en1,en2
cc
100   continue
cc      
      do 200 zk=1E-4,0.3,1E-4 
      aa=xiz*(2*zk**2-2*zk*QT1+QT1**2)+2*xmm
      if(zk-xkF1>0) then                               ! varepsilon_{kF1}        
         en1=xiz*QT1*(QT1-2*zk)/2-sqrt(aa**2+8*DT1**2)/2+g2*xmuB*xB/2
      else
         en1=xiz*QT1*(QT1-2*zk)/2+sqrt(aa**2+8*DT1**2)/2+g2*xmuB*xB/2
      end if
cc
      bb=xiz*(2*zk**2+2*zk*QT1+QT1**2)+2*xmm 
      if(zk-xkF2>0) then                               ! varepsilon_{kF2}   
         en2=xiz*QT1*(-QT1-2*zk)/2+sqrt(bb**2+8*DT1**2)/2+g2*xmuB*xB/2
      else
         en2=xiz*QT1*(-QT1-2*zk)/2-sqrt(bb**2+8*DT1**2)/2+g2*xmuB*xB/2
      end if 
cc
      write(50,'(100f16.8)') zk,en2,en1        
200   continue       
cc
      end subroutine solve_en_TA1
