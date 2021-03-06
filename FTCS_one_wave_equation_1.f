ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                         cc 
cc     Program for solving one-dimentional wave equation                   cc
cc     with finite differencd method (CTCS scheme)                         cc
cc     Version 1.0                                                         cc
cc                                                                         cc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

c------(using parameter)-----------------------------------------------------
c
c      i : 
c      j : 
c     jm : number of discrete points
c     dt : time-step interval
c     dx : distance between two grid points
c      H : height
c     mu : speed of this system
c     UB : Fluid speed values for the past time-step
c      U : Fluid speed values for the present time-step
c     UA : Fluid speed values for the future time-step
c  Eta_A : See level  values for the past time-step
c  Eta   : See level  values for the present time-step
c  Eta_B : See level  values for the future time-step
c          (j = 0, jm+1 : dummy for the cyclic boundary condition)
c           周期境界条件のためのダミー       
c   imax : maximum time step for integration
c   iint : timestep interval between data output
c
c============================================================================

c      program main
c      implicit none 
c      integer i,j
c      real a,b
      parameter (g=9.80655)
      parameter (jm=101)
      parameter (dt=0.1, dx=2.0, mu=1.0, H=1.0)
      parameter (imax= 2000, iint= 20)
c      dimension v(0:jm)
c      dimention 
      dimension UB(0:jm+1,1:imax),U(0:jm+1,1:imax),UA(0:jm+1,1:imax)
      dimention Eta_A(0:jm+1,1:imax),Eta(0:jm+1,1:imax),Eta_B(0:jm+1,1:imax)
      open(5000,file="result_wave_U_sta_1.dat",status="new")
      open(6000,file="result_wave_Eta_sta_1.dat",status="new")
c
c------(set initial condition)------------------------------------------------
c       Gaussian distribution centered at i=50
c       with the amplitude of 1.0
c-----------------------------------------------------------------------------
c
      do 100 j=0,jm+1
        U(j,2) = exp(-(float(j-20)/10.)**2)
        if (U(j,2).lt.0.01) U(j,2) = 0.0
        UB(j,2) = U(j,2)
        UA(j,2) = 0.0
  100 continue

      do 200 j=0,jm+1
        Eta(j,2) = exp(-(float(j-20)/10.)**2)
        if (Eta(j,2).lt.0.01) Eta(j,2) = 0.0
        Eta_B(j,2) = Eta(j,2)
        Eta_A(j,2) = 0.0
  200 continue
c
c------(initial condition detaset)--------------------------------------------
c     Write initial condition to dataset 'fort.1'
c-----------------------------------------------------------------------------
c
      U(0,1)=0
      do 110 j=1,jm+1
      U(j,1)=U(j-1,1)+1
 110  continue

      Eta(0,1)=0
      do 210 j=1,jm+1
      Eta(j,1)=Eta(j-1,1)+1
 210  continue


c      do 110 j=1,jm
c        write(iw,'(2f10.4)') C(j,1),C(j,2)
c  110 continue
c
c-------------------------------------------------------
c     Calculate parameter
c-------------------------------------------------------
c
      a=(mu*dt)/(dx) 
      b=(g*dt)/(dx)
      c=(2*g*dt)/(dx)
      d=(2*mu*dt)/(dx) 
      write(*,*) a,b,c,d
c
c==============================================
c     Calculate numerical integration
c==============================================
c
      do 1000 its=3,imax

!      write(*,*) a,b
c
c---------------------------------------------
c     Calculate next time-step
c---------------------------------------------
c
      do 1001 j=1,jm
         UA(j,its)=U(j,its-2) - a*(U(j+1,its-1)-U(j-1,its-1))
     *   -b*(U(j+1,its-1)-U(j-1,its-1))
         Eta_A(j,its)=Eta(j,its-2) - a*(Eta(j+1,its-1)-Eta(j-1,its-1))
     *   -b*(Eta(j+1,its-1)-Eta(j-1,its-1))   
 1001 continue   

!      write(*,*) CA
      
c     
c---------------------------------------------
c     Set boundary condition
c    (cyclic boundary condition)
c---------------------------------------------
c
      CA(0,its) = CA(jm,its)
      CA(jm+1,its) = CA(1,its)
      Eta_A(0,its) = Eta_A(jm,its)
      Eta_A(jm+1,its) = Eta_A(1,its)

c
c---------------------------------------------
c     Shift one time-step for the next loop
c---------------------------------------------
c
      do 130 j=0,jm+1
        CB(j,its) = C(j,its)  !間違いではない
        C(j,its) = CA(j,its)
 130  continue
      do 230 j=0,jm+1
        Eta_B(j,its) = Eta(j,its)  !間違いではない
        Eta(j,its) = Eta_A(j,its)
 230  continue
c
c---------------------------------------------
c     Write data to dataset 'fort.(iw)'
c---------------------------------------------
 1000  continue

c      if (mod(its,iint).eq.0) then
c        iw=iw+1

        do 500 j=1,jm
          write(5000,50) (U(j,m),m=1,imax)
 500   continue
        do 600 j=1,jm
          write(6000,50) (Eta(j,m),m=1,imax)
        
c      end if

c 1000 continue
c      write(5000,'(f20.4)') C(1:its,1:jm)
      write(*,*)"(^-^)/~~ CAL END"
 50   format('',1x,f5.0,1999(f10.5,1x))
      close(5000)
      close(6000) 
      stop
      end
