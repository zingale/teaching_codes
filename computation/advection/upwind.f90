! solve the linear advection equation via. first order upwinding
!
! u_t + a u_x = 0

program upwind

  implicit none

  integer, parameter :: nx = 64
  integer, parameter :: ng = 2

  double precision, dimension(2*ng + nx) :: x, u, unew, uinit

  double precision, parameter :: xmin = 0.0
  double precision, parameter :: xmax = 1.0
  double precision :: dx

  ! in the linear advection equation, the courant number
  ! is a dt/dx, so we don't need to specify a or dt
  double precision, parameter :: cfl = 0.8

  double precision, parameter :: pi = 3.14159

  integer :: nsteps

  integer :: imin, imax

  integer :: i, n
  
  integer, parameter :: init=2

  ! create the grid
  dx = (xmax - xmin)/dble(nx)

  imin = ng+1
  imax = ng+nx

  do i = 1, 2*ng+nx
     x(i) = (i-ng-0.5d0)*dx + xmin
  enddo

  ! compute the number of steps needed to advect once
  nsteps = floor((xmax-xmin)/(cfl*dx))
  print *, '# nsteps = ', nsteps, ' nx = ', nx

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
           u(i) = 0.d0
        endif
     endif

  enddo

  uinit(:) = u(:)
  
  do n = 1, nsteps

     ! fill the boundary conditions -- just the innermost GC
     u(imin-2) = u(imax-1)
     u(imin-1) = u(imax)

     u(imax+1) = u(imin)
     u(imax+2) = u(imin+1)

     ! evolve
     do i = imin, imax

        ! FTCS
        ! unew(i) = u(i) - 0.5*cfl*(u(i+1) - u(i-1))

        ! FTBS (upwind)
        ! unew(i) = u(i) - cfl*(u(i) - u(i-1))
        
        ! FTFS (downwind)
        ! unew(i) = u(i) - cfl*(u(i+1) - u(i))

        ! Lax-Friedrichs
        ! unew(i) = 0.5*(1.d0 + cfl)*u(i-1) + 0.5*(1.d0 - cfl)*u(i+1)

        ! Lax-Wendroff
        ! unew(i) = u(i) - 0.5*cfl*(u(i+1) - u(i-1)) + 0.5*cfl*cfl*(u(i+1) - 2.d0*u(i) + u(i-1))

        ! Beam-Warming
        ! unew(i) = 0.5*cfl*(cfl-1.d0)*u(i-2) + cfl*(2.d0-cfl)*u(i-1) + 0.5*(cfl-1.d0)*(cfl-2.d0)*u(i)

        ! Fromm's
        unew(i) = -0.25*(1.d0-cfl)*cfl*u(i-2) + 0.25*(5.d0-cfl)*cfl*u(i-1) + &
             0.25*(1.d0-cfl)*(4.d0+cfl)*u(i) - 0.25*(1.d0-cfl)*cfl*u(i+1)
     enddo

     u(:) = unew(:)
     
  enddo
     
  ! print out
  do i = imin, imax
     print *, x(i), uinit(i), u(i)
  enddo


end program upwind
  
  
  
