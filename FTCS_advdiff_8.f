ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                         cc 
cc     Program for solving one-dimentional advective diffusion equation    cc
cc     with finite differencd method (FTCS scheme)                         cc
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
c      D : diffusion coefficient
c      U : velocity of the basic flow
c     CB : values for the past time-step
c      C : values for the present time-step
c     CA : values for the future time-step
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
      parameter (jm=101)
      parameter (dt=600, dx=300, U=0.2, D=20)
      parameter (imax= 2000, iint= 20)
      dimension v(0:jm)
      dimension CB(0:jm+1,1:imax),C(0:jm+1,1:imax),CA(0:jm+1,1:imax)
      open(5000,file="result6.dat",status="replace")
c
c------(set initial condition)------------------------------------------------
c       Gaussian distribution centered at i=50
c       with the amplitude of 1.0
c-----------------------------------------------------------------------------
c
      do 100 j=0,jm+1
        C(j,2) = exp(-(float(j-20)/10.)**2)
c        if (C(j,2).lt.0.01) C(j,2) = 0.0
        CB(j,2) = C(j,2)
        CA(j,2) = 0.0
  100 continue
c
c------(initial condition detaset)--------------------------------------------
c     Write initial condition to dataset 'fort.1'
c-----------------------------------------------------------------------------
c
      C(0,1)=0
      do 130 j=1,jm+1
      C(j,1)=C(j-1,1)+1
 130  continue


c      do 110 j=1,jm
c        write(iw,'(2f10.4)') C(j,1),C(j,2)
c  110 continue
c
c-------------------------------------------------------
c     Calculate parameter
c-------------------------------------------------------
c
      a=(U*dt)/(2*dx) 
      b=(D*dt)/(dx**2)
      write(*,*) a,b,C(0,2)
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
      do 200 j=1,jm
         CA(j,its)=C(j,its-1) - a*(C(j+1,its-1)-C(j-1,its-1))
     *   +b*(C(j+1,its-1)-2*C(j,its-1)+C(j-1,its-1))
 200  continue   

!      write(*,*) CA
      
c     
c---------------------------------------------
c     Set boundary condition
c    (cyclic boundary condition)
c---------------------------------------------
c
      CA(0,its) = CA(jm,its)
      CA(jm+1,its) = CA(1,its)
c
c---------------------------------------------
c     Shift one time-step for the next loop
c---------------------------------------------
c
      do 210 j=0,jm+1
        CB(j,its) = C(j,its)  !間違いではない
        C(j,its) = CA(j,its)
  210 continue
c
c---------------------------------------------
c     Write data to dataset 'fort.(iw)'
c---------------------------------------------
 1000  continue

c      if (mod(its,iint).eq.0) then
c        iw=iw+1
        do 220 j=1,jm
          write(5000,50) (C(j,m),m=1,imax)
  220   continue
c      end if

c 1000 continue
c      write(5000,'(f20.4)') C(1:its,1:jm)
      write(*,*)"(^-^)/~~ CAL END"
 50   format('',1x,f5.0,1999(f10.5,1x))
      close(5000)
      stop
      end
