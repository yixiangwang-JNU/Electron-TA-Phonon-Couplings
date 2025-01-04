**************************************************
*     solve the chemical potential at fixed n0
**************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.14159 
      real(8), parameter :: xm=2.5               ! Dirac mass
      real(8), parameter :: xi=120               ! in-plane band inversion parameter
      real(8), parameter :: xiz=-200             ! z-direction band inversion parameter
      real(8), parameter :: v=295.313            ! in-plane Fermi velocity 
      real(8), parameter :: vz=0                 ! z-direction Fermi velocity      
      real(8), parameter :: xmuB=0.05788         ! Bohr magneton 
      real(8), parameter :: g1=-6                ! Zeeman parameter1
      real(8), parameter :: g2=10                ! Zeeman parameter2
cc    real(8), parameter :: xB=1.2               ! magnetic field
cc    real(8), parameter :: xlB=25.6/sqrt(xB)    ! magnetic length
      real(8), parameter :: xn0=14               ! carrier density 
      real(8), parameter :: gamm=0.01            ! linewidth 
cc    
      common xB,xlB,n  
cc
      end module value
cc      
**********************
*     main program
**********************
      program main
      include 'link_fnl_shared.h'
      use TWODQ_INT
	use value
      implicit double precision (a-h,o-z)
      external f0p,f0m,fnpp,fnpm,G,H
      open(10,file='10.dat')  
      open(20,file='20.dat')        
cc
      xmu1=2.8 
cc    do 100 xB=8.35,15,0.05                     ! magnetic field
      xB=15
      xlB=25.6/sqrt(xB)                          ! magnetic length
cc      
      Res2=100
      do 200 xmu=xmu1-0.2,80,0.01                ! chemical potential 
         write(*,*) xmu 
         Res1=Res2
         A=0                                     ! lower limit of integration
         B=xmu                                   ! upper limit of integration
cc         
         call TWODQ(f0p,A,B,G,H,Res0p)           ! n=0, lambda=1 LL
         call TWODQ(f0m,A,B,G,H,Res0m)           ! n=0, lambda=-1 LL
         Res0p=Res0p*gamm*1D5/(4*Pi*Pi*Pi*xlB*xlB)
         Res0m=Res0m*gamm*1D5/(4*Pi*Pi*Pi*xlB*xlB)
         Res2=Res0p+Res0m
cc         
         n=1
         call TWODQ(fnpp,A,B,G,H,Resnpp)         ! n=1, s=1, lambda=1 LL
         call TWODQ(fnpm,A,B,G,H,Resnpm)         ! n=1, s=1, lambda=-1 LL
         Resnpp=Resnpp*gamm*1D5/(4*Pi*Pi*Pi*xlB*xlB)
         Resnpm=Resnpm*gamm*1D5/(4*Pi*Pi*Pi*xlB*xlB)
         Res2=Res2+Resnpp+Resnpm 
cc
         write(10,'(100f16.8)') xB,xmu,Res2
         if(Res2>xn0) then                       ! judge which one should we choose
            error1=xn0-Res1                  
            error2=Res2-xn0
            if(error1>error2) then 
               xmu1=xmu 
               write(20,'(100f16.8)') xB,xmu1               
            else 
               xmu1=xmu-0.01
               write(20,'(100f16.8)') xB,xmu1               
            end if
cc
            exit                                 ! exit from the present cycle 
         end if
200   continue
      write(10,*)
cc
100   continue
cc      
      end program main
cc      
**************************
*     n=0, lambda=1 LL
**************************
      real function f0p(epsilon,zk)
      use value
      implicit double precision (a-h,o-z)
cc   
      aa=xm-xiz*zk*zk-xi/(xlB*xlB)+xmuB*xB*g1/2
      e0p=-aa+xmuB*xB*g2/2
cc      
      if(e0p>0) then
         temp=(epsilon-e0p)**2+gamm**2
         f0p=1/temp
      else 
         f0p=0 
      end if
cc
      end function f0p
cc
***************************
*     n=0, lambda=-1 LL
***************************
      real function f0m(epsilon,zk)
      use value
      implicit double precision (a-h,o-z)
cc
      aa=xm-xiz*zk*zk-xi/(xlB*xlB)+xmuB*xB*g1/2
      e0m=aa+xmuB*xB*g2/2
cc      
      if(e0m>0) then
        temp=(epsilon-e0m)**2+gamm**2
        f0m=1/temp
      else
        f0m=0  
      end if 
cc
      end function f0m
cc      
*******************************
*     n>0, s=1, lambda=1 LL 
*******************************
      real function fnpp(epsilon,zk)
      use value
      implicit double precision (a-h,o-z)
cc   
      aa=xm-2.*xi*n/(xlB*xlB)-xiz*zk*zk-xmuB*xB*g2/2
      enp=sqrt(aa**2+2*n*v*v/(xlB*xlB))+xi/(xlB*xlB)-xmuB*xB*g1/2
cc
      if(enp>0) then
        temp=(epsilon-enp)**2+gamm**2  
        fnpp=1/temp
      else 
        fnpp=0
      end if 
cc      
      end function fnpp
cc
********************************
*     n>0, s=1, lambda=-1 LL 
********************************
      real function fnpm(epsilon,zk)
      use value
      implicit double precision (a-h,o-z)
cc   
      aa=xm-2.*xi*n/(xlB*xlB)-xiz*zk*zk+xmuB*xB*g2/2
      enm=sqrt(aa**2+2*n*v*v/(xlB*xlB))-xi/(xlB*xlB)+xmuB*xB*g1/2
cc
      if(enm>0) then
        temp=(epsilon-enm)**2+gamm**2
        fnpm=1/temp
      else 
        fnpm=0
      end if
cc      
      end function fnpm
cc
******************************************
*     the lower limit of the variable zk
******************************************
      real function G(zk)
      use value
      implicit double precision (a-h,o-z)
      G=-Pi      
      end function G
cc
******************************************
*     the upper limit of the variable z
******************************************
      real function H(zk)
      use value      
      implicit double precision (a-h,o-z)
      H=Pi     
      end function H
cc 