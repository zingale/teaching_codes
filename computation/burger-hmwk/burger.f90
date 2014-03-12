! solve Burger's equation on a finite-difference grid
!
! set conserved = .false. to solve in non-conservative form:
! u_t + u u_x = 0
!
! set conserved = .true. to solve in conservative form:
! u_t + [0.5 u^2]_x = 0
!
!
! here our grid is
!
!      X                                            X
! +----+----+----+-- // -+----+----+----+-- // +----+----+
! -1   0    1    2      i-1   i   i+1  i+2    nx-2 nx-1  nx
!
! X marks the physical boundary -- we have one guardcell on the 
! EACH end.
!
! outflow boundaries mean that u
!
! We will update from i=0 to i=nx-1, and apply the boundary conditions
! u(-1) = u(0)
! u(nx) = u(nx-1)
!
! these allow us to use the same stencil throughout the entire
! physical domain.
!
program advection

  implicit none

  integer, parameter :: nx = 64

  double precision, dimension(-1:nx) :: x, u, unew, uinit

  double precision, parameter :: xmin = 0.0
  double precision, parameter :: xmax = 1.0
  double precision :: dx

  ! in the linear advection equation, the courant number
  ! is a dt/dx, so we don't need to specify a or dt
  double precision, parameter :: cfl = 0.8

  logical, parameter :: conserved = .true.

  integer :: i, imin, imax
  
  integer, parameter :: init=2

  double precision :: time, dt
  double precision, parameter :: tmax = 0.2
  double precision, parameter :: SMALL = 1.d-12


  ! create the grid
  dx = (xmax - xmin)/dble(nx-1)

  ! the range of points that we update is 0 to nx-1
  imin = 0
  imax = nx-1


  do i = -1, nx
     x(i) = i*dx + xmin
  enddo



  ! initialize
  do i = imin, imax
     
     if (init == 1) then
        ! first set of boundary condtions
        if (x(i) < 0.5d0) then
           u(i) = 2.d0
        else
           u(i) = 1.d0
        endif

     else if (init == 2) then
        ! second set of boundary conditions
        if (x(i) < 0.5d0) then
           u(i) = 1.d0
        else
           u(i) = 2.d0
        endif
     endif

  enddo

  ! store the initial conditions for output later
  uinit(:) = u(:)
  
  time = 0.0

  do while (time < tmax)

     ! fill the boundary conditions -- just the innermost GC
     u(-1) = u(0)
     u(nx) = u(nx-1)

     ! compute the timestep
     dt = 1.e33
     do i = imin, imax
        dt = min(dt, dx/(abs(u(i)) + SMALL))
     enddo
     dt = CFL*dt

     ! evolve
     do i = imin, imax

        if (conserved) then
           unew(i) = u(i) - 0.5*(dt/dx)*(u(i)**2 - u(i-1)**2)
        else
           unew(i) = u(i) - (dt/dx)*u(i)*(u(i) - u(i-1))
        endif

     enddo

     u(:) = unew(:)
     
     time = time + dt
  enddo
     
  ! print out
  do i = imin, imax
     print *, i, x(i), uinit(i), u(i)
  enddo


end program advection
  
  
  
