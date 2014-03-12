! solve the linear advection equation via. FTCS, on a finite-difference grid
!
! u_t + a u_x = 0

program advection

  implicit none

  integer, parameter :: nx = 64
  integer, parameter :: ng = 1

  double precision, dimension(0:ng+nx-1) :: x, u, unew, uinit

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
  
  integer, parameter :: init=1

  ! create the grid
  dx = (xmax - xmin)/dble(nx-1)

  ! x = 0 is the same as x = nx-1, so we don't loop over it
  imin = 1
  imax = nx-1

  do i = 0, ng+nx-1
     x(i) = i*dx + xmin
  enddo

  ! compute the number of steps needed to advect once
  nsteps = (xmax-xmin)/(cfl*dx)+ 1

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
     u(0) = u(nx-1)
     u(nx) = u(1)

     ! evolve
     do i = imin, imax

        ! FTCS
        !unew(i) = u(i) - 0.5*cfl*(u(i+1) - u(i-1))
        
        ! FTBS (upwind)
        unew(i) = u(i) - cfl*(u(i) - u(i-1))
        
        ! FTFS (downwind)
        ! unew(i) = u(i) - cfl*(u(i+1) - u(i))

     enddo

     u(:) = unew(:)
     
  enddo
     
  ! print out
  do i = 0, imax
     print *, x(i), uinit(i), u(i)
  enddo


end program advection
  
  
  
