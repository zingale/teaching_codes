! solve the linear advection equation on a finite-volume grid using
! using Godunov's method
!
! u_t + a u_x = 0
!

program advect_godunov

  implicit none

  ! the number of zones (nx) and number of guardcells (ng)
  integer, parameter :: nx = 64
  integer, parameter :: ng = 2

  ! the domain size
  double precision, parameter :: xmin = 0.d0
  double precision, parameter :: xmax = 1.d0

  ! the advection velocity
  double precision, parameter :: a = 1.d0

  ! the CFL number
  double precision, parameter :: cfl = 0.8d0

  ! initial condition type
  integer, parameter :: inittype = 2
  
  ! maximum simulation time
  double precision, parameter :: tmax = 5.d0


  double precision, dimension(2*ng + nx) :: x
  double precision, dimension(2*ng + nx) :: u, ul, ur, f

  double precision :: dx
  double precision :: time, dt


  ! setup the grid and set the initial conditions
  call grid(nx, ng, xmin, xmax, dx, x)

  call init(nx, ng, inittype, x, u)


  time = 0.d0

  call output(nx, ng, inittype, a, time, x, u)


  ! evolution loop -- construct the interface states, solve the Riemann 
  ! problem, and update the solution to the new time level
  do while (time < tmax)

     call fillBC(nx, ng, u)

     call timestep(nx, ng, dx, a, cfl, u, dt)

     call states(nx, ng, u, ul, ur)

     call riemann(nx, ng, a, ul, ur, f)

     call update(nx, ng, dx, dt, u, f)

     time = time + dt

  enddo
     
  call output(nx, ng, inittype, a, time, x, u)


end program advect_godunov
  


!============================================================================
! grid: create the grid
!============================================================================
subroutine grid(nx, ng, xmin, xmax, dx, x)

  implicit none

  integer :: nx, ng
  double precision :: xmin, xmax

  double precision :: dx
  double precision, dimension(2*ng+nx) :: x

  integer :: i

  ! create the grid
  dx = (xmax - xmin)/dble(nx)

  do i = 1, 2*ng+nx
     x(i) = (i-ng-0.5d0)*dx + xmin
  enddo

  return
end subroutine grid



!============================================================================
! init: set the initial conditions
!============================================================================
subroutine init(nx, ng, inittype, x, u)

  implicit none

  integer :: nx, ng
  integer :: inittype

  double precision, dimension(2*ng+nx) :: x, u

  integer :: i, imin, imax

  double precision, parameter :: pi = 3.14159d0

  imin = ng+1
  imax = ng+nx

  ! loop over all the zones and set the initial conditions.  To be
  ! consistent with the finite-volume discretization, we should store
  ! the zone averages here, but, to second-order accuracy, it is
  ! sufficient to evaluate the initial conditions at the zone center.

  do i = imin, imax
     
     if (inittype == 1) then
        ! sin wave
        u(i) = sin(2.d0*pi*x(i))

     else if (inittype == 2) then
        ! square wave
        if (x(i) > 0.333d0 .and. x(i) < 0.666d0) then
           u(i) = 1.d0
        else
           u(i) = 0.d0
        endif
     endif

  enddo

  return
end subroutine init



!============================================================================
! output: write out the solution
!============================================================================
subroutine output(nx, ng, inittype, a, time, x, u)

  implicit none

  integer :: nx, ng

  integer :: inittype

  double precision :: a

  double precision :: time

  double precision, dimension(2*ng+nx) :: x, u

  character (len=4) :: time_string

  integer :: i, imin, imax

  imin = ng+1
  imax = ng+nx

  
  ! open the output file
  write (time_string, '(f4.2)') time

  open(unit=10, file="advect_t="//time_string, status="unknown")

  write (10,*) "# advection problem: u_t + a u_x = 0"
  write (10,*) "# a = ", a
  write (10,*) "# init = ", inittype
  write (10,*) "# time = ", time


  do i = imin, imax
     write (10,*) x(i), u(i)
  enddo

  return
end subroutine output



!============================================================================
! fillBC: fill the boundary conditions
!============================================================================
subroutine fillBC(nx, ng, u)

  implicit none

  integer :: nx, ng
  double precision, dimension(2*ng+nx) :: u

  integer :: i, imin, imax

  imin = ng+1
  imax = ng+nx

  ! left boundary
  do i = 1, imin-1
     u(i) = u(imax-ng+i)
  enddo

  
  ! right boundary
  do i = imax+1, 2*ng+nx
     u(i) = u(i-imax+ng)
  enddo

  return
end subroutine fillBC
  
  

!============================================================================
! timestep: compute the new timestep
!============================================================================
subroutine timestep(nx, ng, dx, a, cfl, u, dt)
  
  implicit none

  integer :: nx, ng
  double precision :: dx
  double precision :: a, cfl
  double precision, dimension(2*ng+nx) :: u
  double precision :: dt
  
  integer :: i, imin, imax

  imin = ng+1
  imax = ng+nx
  
  
  ! in the linear advection equation, the timestep is trivial
  dt = cfl*dx/abs(a)
    
  return
end subroutine timestep
  


!============================================================================
! states: compute the interface states used in solving the Riemann problem
!============================================================================
subroutine states(nx, ng, u, ul, ur)

  implicit none

  integer :: nx, ng

  double precision, dimension(2*ng+nx) :: u, ul, ur

  integer :: i, imin, imax

  imin = ng+1
  imax = ng+nx


  ! for each interface, we want to construct the left and right
  ! states.  Here, interface i refers to the left edge of zone i
  
  ! interfaces imin to imax+1 affect the data in zones [imin,imax]
  do i = imin, imax+1

     ! the left state on the current interface comes from zone i-1.  
     ul(i) = u(i-1)
    
     ! the right state on the current interface comes from zone i
     ur(i) = u(i)

  enddo

  return
end subroutine states
  


!============================================================================
! riemann: solve the Riemann problem
!============================================================================
subroutine riemann(nx, ng, a, ul, ur, f)

  implicit none

  integer :: nx, ng
  double precision, dimension(2*ng+nx) :: ul, ur, f

  double precision :: a

  integer :: i, imin, imax

  imin = ng+1
  imax = ng+nx


  ! loop over all the interfaces and solve the Riemann problem.  Here,
  ! since we are doing the linear advection eq, we just use the advection
  ! velocity to tell us which direction is upwind, and use that state
  
  if (a >= 0.d0) then
     
     do i = imin, imax+1
        f(i) = a*ul(i)
     enddo
     
  else

     do i = imin, imax+1
        f(i) = a*ur(i)
     enddo
     
  endif

  return
end subroutine riemann



!============================================================================
! update: conservatively update the solution to the new time level
!============================================================================
subroutine update(nx, ng, dx, dt, u, f)

  implicit none

  integer :: nx, ng

  double precision :: dx, dt

  double precision, dimension(2*ng+nx) :: u, f

  integer :: i, imin, imax

  imin = ng+1
  imax = ng+nx

  do i = imin, imax
     u(i) = u(i) + (dt/dx)*(f(i) - f(i+1))
  enddo
  
  return
end subroutine update



