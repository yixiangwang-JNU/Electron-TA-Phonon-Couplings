*****************************************************
*     Calculate the renormalized phonon frequency
*****************************************************
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
      real(8), parameter :: eta=0.01             ! linewidth parameter 
      real(8), parameter :: xkb=0.08614          ! Boltzmann constant    
      real(8), parameter :: Tem=1E-6             ! temperature 
cc      
      end module value
cc      
**********************
*     main program
**********************
      program main
	use value
      implicit double precision (a-h,o-z)
      character(100) filename1,filename2 
      real(8) gLArray(5),gTArray(5),romegaL(5),romegaT(5)
      write(filename1,'(a2,f12.6,a10)') 'T=',tem,',LA.dat'
      write(filename2,'(a2,f12.6,a10)') 'T=',tem,',TA.dat' 
      open(10,file=filename1)        
      open(20,file=filename2)        
cc
      gLArray(1)=1
      gLArray(2)=2
      gLArray(3)=3
      gLArray(4)=4
      gLArray(5)=5      
cc      
      do 100 qz=0.001,0.4,0.001     
      do 101 i=1,5
      gL=gLarray(i)
cc
cc    Case1: alpha=1,beta=1 
      chi1=0 
      do zk1=-0.5,0.5001,0.001
cc    do zk1=-0.5,0,0.001          
         zk2=zk1+qz
         if(zk2>0.5) zk2=zk2-1  
         call solve_em(em1,zk1)
         call solve_em(em2,zk2)
cc       Nf1=0
cc       Nf2=0
cc       if(em1<xmu) Nf1=1
cc       if(em2<xmu) Nf2=1        
         call Ferdis(em1,Fermi1)
         call Ferdis(em2,Fermi2)
         aa=em1-em2
         chi1=chi1+aa*(Fermi1-Fermi2)/(aa*aa+eta*eta)
      end do
      chi1=chi1/(1000*deg)            
cc             
cc    Case2: alpha=2,beta=2 
      chi2=0
      do zk1=-0.5,0.5001,0.001
cc    do zk1=-0.5,0,0.001
         zk2=zk1+qz
         if(zk2>0.5) zk2=zk2-1 
         call solve_ep(ep1,zk1)
         call solve_ep(ep2,zk2)
cc       Nf1=0
cc       Nf2=0
cc       if(ep1<xmu) Nf1=1
cc       if(ep2<xmu) Nf2=1
         call Ferdis(ep1,Fermi1)
         call Ferdis(ep2,Fermi2)
         aa=ep1-ep2
         chi2=chi2+aa*(Fermi1-Fermi2)/(aa*aa+eta*eta)
      end do  
      chi2=chi2/(1000*deg)    
cc      
      aa=1+132.397*gL*gL*(chi1+chi2)  
      romegaL(i)=0
      if(aa>0) romegaL(i)=sqrt(aa)*qz 
101   continue
      write(10,'(100f16.8)') qz,(romegaL(i),i=1,5) 
100   continue 
cc      
      gTArray(1)=1
      gTArray(2)=2
      gTArray(3)=3
      gTArray(4)=4
      gTArray(5)=5   
      do 200 qz=0.001,0.4,0.001   
      do 201 i=1,5
      gT=gLarray(i)    
cc    Case3: alpha=1,beta=2
      chi3=0 
      do zk1=-0.5,0.5001,0.001
         zk2=zk1+qz
         if(zk2>0.5) zk2=zk2-1 
         call solve_em(em,zk1)
         call solve_ep(ep,zk2)
cc       Nf1=0
cc       Nf2=0
cc       if(em<xmu) Nf1=1
cc       if(ep<xmu) Nf2=1
         call Ferdis(em,Fermi1)
         call Ferdis(ep,Fermi2)
         aa=em-ep
         chi3=chi3+aa*(Fermi1-Fermi2)/(aa*aa+eta*eta)
      end do  
      chi3=chi3/(1000*deg)        
cc
cc    Case4: alpha=2,beta=1
      chi4=0
      do zk1=-0.5,0.5001,0.001
         zk2=zk1+qz
         if(zk2>0.5) zk2=zk2-1          
         call solve_ep(ep,zk1)
         call solve_em(em,zk2)
cc       Nf1=0
cc       Nf2=0
cc       if(ep<xmu) Nf1=1
cc       if(em<xmu) Nf2=1
         call Ferdis(ep,Fermi1)
         call Ferdis(em,Fermi2)         
         aa=ep-em
         chi4=chi4+aa*(Fermi1-Fermi2)/(aa*aa+eta*eta)
      end do
      chi4=chi4/(1000*deg) 
cc
      bb=1+132.397*gT*gT*(chi3+chi4)  
      romegaT(i)=0
      if(bb>0) romegaT(i)=sqrt(bb)*qz 
201   continue      
      write(20,'(100f16.8)') qz,romegaT
cc
200   continue
cc      
      end program main
cc      
***************************
*     n=0, lambda=-1 LL
***************************
      subroutine solve_em(em,zk)
      use value
      implicit double precision (a-h,o-z)
cc
      aa=xm-xiz*zk*zk-xi/(xlB*xlB)+xmuB*xB*g1/2
      em=aa+xmuB*xB*g2/2
cc
      end subroutine solve_em
cc      
**************************
*     n=0, lambda=1 LL
**************************
      subroutine solve_ep(ep,zk)
      use value
      implicit double precision (a-h,o-z)
cc   
      aa=xm-xiz*zk*zk-xi/(xlB*xlB)+xmuB*xB*g1/2
      ep=-aa+xmuB*xB*g2/2 
cc
      end subroutine solve_ep
cc
*************************************
*     Fermi distribution function    
*************************************
      subroutine Ferdis(en,Fermi)
      use value
      implicit double precision (a-h,o-z)       
cc
      if(T<1E-5.and.en<xmu) Fermi=1         ! zero temperature case 
      if(T<1E-5.and.en>xmu) Fermi=0       
      aa=exp((en-xmu)/(xkb*Tem))
      Fermi=1/(aa+1)
cc      
      end subroutine Ferdis 
