*******************************************************************
*     determine the chemical potential at fixed carrier density
*******************************************************************
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
      real(8), parameter :: xn0=8                ! carrier density 
      real(8), parameter :: eta=0.01             ! linewidth 
cc      
      common xmu,xB,xlB,n  
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
      character(100) filename1,filename2 
      external f0p,f0m,fnpp,fnpm,G,H
      write(filename1,'(a10)') '10.dat'
      write(filename2,'(a20)') 'results.dat'
      open(10,file=filename1)     
      open(20,file=filename2)     
cc
      xmu1=3.38   
      do 100 xB=6.41,15.01,0.05 
      xlB=25.6/sqrt(xB)
cc
      Res2=100                                     ! initial value  
      do 200 xmu=xmu1,50,0.01                      ! chemical potential 
         Res1=Res2 
         A=0                                       ! lower limit of integration
         B=xmu                                     ! upper limit of integration
cc         
         call TWODQ(f0p,A,B,G,H,Res0p)             ! n=0, s=1 LL
         call TWODQ(f0m,A,B,G,H,Res0m)             ! n=0, s=-1 LL
         Res0p=Res0p*eta*1D5/(4*Pi*Pi*Pi*xlB*xlB)
         Res0m=Res0m*eta*1D5/(4*Pi*Pi*Pi*xlB*xlB)            
cc
         n=1                
         call TWODQ(fnpp,A,B,G,H,Res1pp)           ! n=1, s=1, lambda=1 LL
         call TWODQ(fnpm,A,B,G,H,Res1pm)           ! n=1, s=1, lambda=-1 LL
         Res1pp=Res1pp*eta*1D5/(4*Pi*Pi*Pi*xlB*xlB)
         Res1pm=Res1pm*eta*1D5/(4*Pi*Pi*Pi*xlB*xlB)
cc
cc       n=2
cc       call TWODQ(fnpp,A,B,G,H,Res2pp)           ! n=2, s=1, lambda=1 LL
cc       call TWODQ(fnpm,A,B,G,H,Res2pm)           ! n=2, s=1, lambda=-1 LL
cc       Res2pp=Res2pp*eta*1D5/(4*Pi*Pi*Pi*xlB*xlB)
cc       Res2pm=Res2pm*eta*1D5/(4*Pi*Pi*Pi*xlB*xlB)
cc         
cc       n=3
cc       call TWODQ(fnpp,A,B,G,H,Res3pp)           ! n=3, s=1, lambda=1 LL
cc       call TWODQ(fnpm,A,B,G,H,Res3pm)           ! n=3, s=1, lambda=-1 LL
cc       Res3pp=Res3pp*eta*1D5/(4*Pi*Pi*Pi*xlB*xlB)
cc       Res3pm=Res3pm*eta*1D5/(4*Pi*Pi*Pi*xlB*xlB) 
cc         
cc       n=4
cc       call TWODQ(fnpp,A,B,G,H,Res4pp)           ! n=4, s=1, lambda=1 LL
cc       call TWODQ(fnpm,A,B,G,H,Res4pm)           ! n=4, s=1, lambda=-1 LL
cc       Res4pp=Res4pp*eta*1D5/(4*Pi*Pi*Pi*xlB*xlB)
cc       Res4pm=Res4pm*eta*1D5/(4*Pi*Pi*Pi*xlB*xlB)    
cc
         Res2=Res0p+Res0m+Res1pp+Res1pm
cc    &       +Res2pp+Res2pm+Res3pp+Res3pm+Res4pp+Res4pm
         write(*,*) xB,xmu,Res2         
         write(10,'(100f16.8)') xB,xmu,Res2,Res0p,Res0m,Res1pp,Res1pm
cc    &                       Res2pp,Res2pm,Res3pp,Res3pm,Res4pp,Res4pm  
cc
         if(Res2>xn0) then                         ! judge which one should we choose
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
            exit                                   ! exit from the present cycle 
         end if
200   continue  
      write(10,*)      
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
         temp=(epsilon-e0p)**2+eta**2
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
        temp=(epsilon-e0m)**2+eta**2
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
      enpp=sqrt(aa**2+2*n*v*v/(xlB*xlB))+xi/(xlB*xlB)-xmuB*xB*g1/2  
cc      
      if(enpp>0) then
        temp=(epsilon-enpp)**2+eta**2  
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
      enpm=sqrt(aa**2+2*n*v*v/(xlB*xlB))-xi/(xlB*xlB)+xmuB*xB*g1/2
cc      
      if(enpm>0) then
        temp=(epsilon-enpm)**2+eta**2
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