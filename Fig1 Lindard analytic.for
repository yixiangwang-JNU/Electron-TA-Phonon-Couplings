*************************************************
*     Calculate the Lindard response function 
*************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.14159 
      real(8), parameter :: xmu=3.33             ! chemical potential 
      real(8), parameter :: xm=5                 ! Dirac mass
      real(8), parameter :: xi=100               ! in-plane band inversion parameter
      real(8), parameter :: xiz=200              ! z-direction band inversion parameter
      real(8), parameter :: v=393.75             ! in-plane Fermi velocity 
      real(8), parameter :: vz=0                 ! z-direction Fermi velocity      
      real(8), parameter :: xmuB=0.05788         ! Bohr magneton 
      real(8), parameter :: xB=6.5               ! magnetic field
      real(8), parameter :: xlB=25.6/sqrt(xB)    ! magnetic length 
      real(8), parameter :: deg=2*Pi*xlB*xlB     ! degeneracy factor   
      real(8), parameter :: g1=-8                ! Zeeman parameter1
      real(8), parameter :: g2=10                ! Zeeman parameter2
cc      
      end module value
cc      
**********************
*     main program
**********************
      program main
	use value
      implicit double precision (a-h,o-z)
      character(100) filename 
      write(filename,'(a100)') 'results.dat'
      open(10,file=filename)        
cc
      xmm=xi/(xlB*xlB)-xm-g1*xmuB*xB/2
      aa=(-xmu-xmm+g2*xmuB*xB/2)/xiz
      bb=(xmu-xmm-g2*xmuB*xB/2)/xiz
      xkF1=0                                     ! Fermi wavevector for the valence band          
      if(aa>=0) xkF1=sqrt(aa)         
      xkF2=0                                     ! Fermi wavevector for the conduction band    
      if(bb>=0) xkF2=sqrt(bb) 
      QL1=2*xkF1                                 ! LA1 phonon wavevector 
      QL2=2*xkF2                                 ! LA2 phonon wavevector                  
      write(*,*) QL1,QL2 
cc    
      do 100 qz=0.001,0.4,0.001     
      aa=(qz+QL1)/(qz-QL1) 
      chi1=log(abs(aa))/(deg*xiz*qz)
cc    
      bb=(qz+QL2)/(qz-QL2)
      chi2=log(abs(bb))/(deg*xiz*qz)
      write(10,'(100f16.8)') qz,chi1,chi2 
100   continue          
cc      
      end program main
cc      