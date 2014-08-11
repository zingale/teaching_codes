! solve the inviscid Burger's equation on a finite-volume grid using
! second-order (linear) reconstruction
!
! u_t + (0.5 u^2)_x = 0
!

program upwind

  implicit none

  integer, parameter :: nx = 256
  integer, parameter :: ng = 2

  double precision, dimension(2*ng + nx) :: x, u, ul, ur, uinit, f

  double precision, parameter :: xmin = 0.d0
  double precision, parameter :: xmax = 1.d0
  double precision :: dx
  
  double precision :: dt, time
  double precision, parameter :: tmax = 0.3d0


  ! specify the CFL number
  double precision, parameter :: cfl = 0.8d0

  double precision, parameter :: pi = 3.14159d0

  integer :: imin, imax

  integer :: i, n
  
  ! initial condition type
  integer, parameter :: init = 2
  

  ! slope type (1=godunov, 2=centered, 3=minmod limiting, 4=MC limiting, 5=superbee
  integer, parameter :: islopetype = 4


  ! create the grid
  dx = (xmax - xmin)/dble(nx)

  imin = ng+1
  imax = ng+nx

  do i = 1, 2*ng+nx
     x(i) = (i-ng-0.5d0)*dx + xmin
  enddo


  ! initialize
  do i = imin, imax
     
     if (init == 1) then
        ! sin wave
        u(i) = sin(2.d0*pi*x(i))

     else if (init == 2) then
        ! square wave
        if (x(i) > 0.333d0 .and. x(i) < 0.666d0) then
           u(i) = 1.d0
        else
           u(i) = -1.d0
        endif
     endif

  enddo

  call fillBC(nx, ng, u)
  uinit(:) = u(:)

  print *, '# init = ', init
  print *, '# islopetype = ', islopetype
  print *, '# cfl = ', cfl

  time = 0.d0

  do while (time < tmax)

     ! fill the boundary conditions -- just the innermost GC
     call fillBC(nx, ng, u)


     ! get the timestep
     call timestep(nx, ng, dx, cfl, u, dt)

     ! compute the interface states
     call states(nx, ng, dx, dt, islopetype, u, ul, ur)

     
     ! solve the Riemann problem
     call riemann(nx, ng, ul, ur, f)


     ! do the conservative update
     do i = imin, imax
        u(i) = u(i) + (dt/dx)*(f(i) - f(i+1))
     enddo

     time = time + dt

  enddo
     
  ! print out
  do i = imin, imax
     print *, x(i), uinit(i), u(i)
  enddo


end program upwind
  


subroutine fillBC(nx, ng, u)

  implicit none

  integer :: nx, ng
  double precision, dimension(2*ng+nx) :: u

  integer :: i, imin, imax

  imin = ng+1
  imax = ng+nx

  ! left boundary
  do i = 1, imin-1
     u(i) = u(imin)
  enddo

  
  ! right boundary
  do i = imax+1, 2*ng+nx
     u(i) = u(imax)
  enddo

  return
end subroutine fillBC
  
  

subroutine timestep(nx, ng, dx, cfl, u, dt)
  
  implicit none

  integer :: nx, ng
  double precision :: dx
  double precision :: a, cfl
  double precision, dimension(2*ng+nx) :: u
  double precision :: dt
  
  integer :: i, imin, imax
  double precision, parameter :: SMALL = 1.d-12

  imin = ng+1
  imax = ng+nx
  
  
  ! in the linear advection equation, the timestep is trivial
  dt = 1.e33
  do i = imin, imax
     dt = min(dt, dx/(abs(u(i)) + SMALL))
  enddo
  dt = cfl*dt

    
  return
end subroutine timestep
  


subroutine states(nx, ng, dx, dt, islopetype, u, ul, ur)

  implicit none

  integer :: nx, ng
  integer :: islopetype

  double precision, dimension(2*ng+nx) :: u, ul, ur
  double precision, dimension(2*ng+nx) :: slope

  double precision :: slope1, slope2

  double precision :: dx, dt

  double precision :: minmod, maxmod

  integer :: i, imin, imax

  imin = ng+1
  imax = ng+nx


  ! compute the centered difference for linear slopes
  do i = imin-1, imax+1

     if (islopetype == 1) then

        ! Godunov's method (piecewise constant)
        slope(i) = 0.d0
        
     else if (islopetype == 2) then

        ! centered difference (equivalent to Fromm's method)
        slope(i) = 0.5*(u(i+1) - u(i-1))/dx

     else if (islopetype == 3) then

        ! minmod limited slope
        slope(i) = minmod((u(i) - u(i-1))/dx, (u(i+1) - u(i))/dx)

     else if (islopetype == 4) then

        ! MC limiter
        slope(i) = minmod(minmod(2.d0*(u(i) - u(i-1))/dx, &
                                 2.d0*(u(i+1) - u(i))/dx), &
                          0.5d0*(u(i+1) - u(i-1))/dx)

     else if (islopetype == 5) then

        ! superbee
        slope1 = minmod((u(i+1) - u(i))/dx, &
                        2.d0*(u(i) - u(i-1))/dx)

        slope2 = minmod(2.d0*(u(i+1) - u(i))/dx, &
                        (u(i) - u(i-1))/dx)

        slope(i) = maxmod(slope1, slope2)

     endif
  enddo


  ! for each interface, we want to construct the left and right
  ! states.  Here, interface i refers to the left edge of zone i
  
  ! interfaces imin to imax+1 affect the data in zones [imin,imax]
  do i = imin, imax+1

     ! the left state on the current interface comes from
     ! zone i-1.  
     ul(i) = u(i-1) + 0.5*dx*(1.d0 - u(i-1)*(dt/dx))*slope(i-1) 

     
     ! the right state on the current interface comes from
     ! zone i
     ur(i) = u(i) - 0.5*dx*(1.d0 + u(i)*(dt/dx))*slope(i)

  enddo

  return
end subroutine states
  


subroutine riemann(nx, ng, ul, ur, f)

  implicit none

  integer :: nx, ng
  double precision, dimension(2*ng+nx) :: ul, ur, f
  double precision :: S, us

  integer :: i, imin, imax

  imin = ng+1
  imax = ng+nx


  ! loop over all the interfaces and solve the Riemann problem.  
  do i = imin, imax+1
     if (ul(i) > ur(i)) then

        ! shock
        S = 0.5*(ul(i) + ur(i))
        
        if (S >= 0) then
           us = ul(i)
        else
           us = ur(i)
        endif

     else

        if (ul(i) >= 0.0) then
           us = ul(i)
        else if (ur(i) <= 0.0) then
           us = ur(i)
        else
           us = 0.d0
        endif

     endif
        
     f(i) = 0.5*us*us

  enddo

  return
end subroutine riemann




function minmod(a,b)

  implicit none

  double precision :: a, b
  double precision :: minmod

  if (abs(a) < abs(b) .and. a*b > 0.d0) then
     minmod = a
  else if (abs(b) < abs(a) .and. a*b > 0) then
     minmod = b
  else
     minmod = 0.d0
  endif

  return 
end function minmod



function maxmod(a,b)

  implicit none

  double precision :: a, b
  double precision :: maxmod

  if (abs(a) > abs(b) .and. a*b > 0.d0) then
     maxmod = a
  else if (abs(b) > abs(a) .and. a*b > 0) then
     maxmod = b
  else
     maxmod = 0.d0
  endif

  return 
end function maxmod
