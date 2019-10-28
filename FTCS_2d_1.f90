!===============================================================
!
!     Program for solving one-dimensional diffusion equation 
!     (dc/dt-a*d2c/dx2=0)
!     with finite difference method (FTCS scheme).
!
!---------------------------------------------------------------
!
!     jm : number of discrete points
!     dt : time-step interval
!     dx : distance between two grid points
!      a : diffusion coefficient
!      c : values for the present time-step
!     ca : values for the future time-step
!          (j = 0, jm+1 : dummy for the cyclic boundary condition)
!   imax : maximum time step for integration
!   iint : timestep interval between data output
!
!===============================================================
      parameter (jm=101,lm=101)                            !number of separate points
      parameter (dt=1.0, dx=2.0, dy=2.0, a=1.0)     !Indicades time-step interval , distance between two grid points , diffusion coefficient 
      parameter (imax= 200, iint= 20)               !imax is maximum time step for integration ,iint is timestep interval between deta output
!
      dimension c(jm+1,lm+1),ca(jm+1,lm+1)                !c and ca has line 
!
!-------------------------------------------------------
!     Set initial condition
!     Gaussian distribution centered at i = 50
!     with the amplitude of 1.0
!-------------------------------------------------------
!
       do 100 l=0,lm+1
          select case (l)
             case (50:101)
                  do 101 j=0,jm+1
                     select case (j) 
                        case (50:101)  
                             c(j,l)=10
                        case default
                             c(j,l)=0
                     end select
 101              continue
             case default 
                  do 102 j=0,jm+1
                     c(j,l)=0
 102              continue   
          end select
 100   continue     


!      do 100 j=0,jm+1
!         select case (j)
!            case (45:55)
!               c(j)=10
!            case (95:101)
!               c(j)=100
!            case default 
!               c(j)=0
!         end select
!        c(j) = exp(-(float(j-50)/10.)**2)
!        if (c(j).lt.0.01) c(j) = 0.0
!        ca(j) = 0.0
!  100 continue

!
!--------------------------------------------------------
!     Write initial condition to dataset 'fort.10'
!--------------------------------------------------------
!
      iw=29
!c      do 110 l=1,lm
         do 111 j=1,jm
            write(iw,'(1x,f10.4)') (c(j,l),l=1,lm)  !Without deciding filename , made filename is decidid [fort.○] automatically
 111     continue
!c 110  continue   
!c   
!c-------------------------------------------------------
!c     Calculate parameter related to the CFL condition
!c-------------------------------------------------------
!c
      cfl = a*dt/(((dx+dy)/2)**2)          !CFL condition for one-dimensional diffusion equation
!c
!c==============================================
!c     Calculate numerical integration
!c==============================================
!c
      do 1000 its=1,imax          !
!c
!c---------------------------------------------
!c     Calculate next time-step with FTCS scheme
!c---------------------------------------------
!c

      do 200 l=1,lm
         do 201 j=1,jm
            ca(j,l)=c(j,l)+cfl*((c(j+1,l)+c(j-1,l)-2*c(j,l))&
                    +(c(j,l+1)+c(j,l-1)-2*c(j,l)))
 201     continue   
 200  continue

!c      do 200 j=1,jm
!c        ca(j) = c(j) + cfl*(c(j+1)+c(j-1)-2.*c(j))           !FTCS scheme
!c  200 continue

!c
!c---------------------------------------------
!c     Set boundary condition
!c    (fixed boundary condition)
!c---------------------------------------------
!c
      ca(lm,jm) = 0.0                 !boundary condition
      ca(lm,1)  = 0.0                 !boundary condition
      ca(1,jm)  = 0.0                 !boundary condition
      ca(1,1)   = 0.0                 !boundary condition
!c
!c---------------------------------------------
!c     Shift one time-step for the next loop
!c---------------------------------------------
!c
      do 210 l=0,lm+1
         do 211 j=0,jm+1             
            c(j,l) = ca(j,l)
 211     continue
 210  continue

!c
!c---------------------------------------------
!c     Write data to dataset 'fort.(iw)'
!c---------------------------------------------
!c
      if (mod(its,iint).eq.0) then  !mod(its,iint)=itsをiintで割った余り
         iw=iw+1
 
      do 220  j=1,jm
         write(iw,10) (c(j,l),l=1,lm)
 220  continue   

!c        do 220 l=1,lm
!c            do 221 j=1,jm
!c              write(iw,'(f10.4)') c(l,j)  !file output
!c 221       continue
!c 220     continue 

      write(*,*)iw


      end if

!c
 1000 continue
 10   format('',f10.4)
!c
      stop
      end
