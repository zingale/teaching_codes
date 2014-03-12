! solve the linear advection equation via L-W, computing the timestep
! on the fly in a f-d method.  This is a precursor for doing the method
! for a non-linear equation.
!
! u_t + a u_x = 0
!
! here our grid is
!
! X                                            X
! +----+----+-- // -+----+----+----+-- // +----+----+
! 0    1    2      i-1   i   i+1  i+2    nx-2 nx-1  nx
!
! X marks the physical boundary -- we have one guardcell on the 
! right end.
!
! periodic boundary conditions mean u_0 = u_{nx-1}, so we only
! need to update one of these points -- we choose the latter.
!
! We will update from i=1 to i=nx-1, and apply the boundary conditions
! u(0) = u(nx-1)
! u(nx) = u(1)
!
! these allow us to use the same stencil throughout the entire
! physical domain.
!
program advection

  implicit none

  integer, parameter :: nx = 64
  integer, parameter :: ng = 1

  double precision, dimension(0:ng+nx-1) :: x, u, unew, uinit

  double precision, parameter :: xmin = 0.0
  double precision, parameter :: xmax = 1.0
  double precision :: dx

  double precision, parameter :: cfl = 0.8

  double precision, parameter :: a = 1.d0

  double precision :: c

  double precision, parameter :: pi = 3.141592653

  double precision, parameter :: small = 1.d-12

  double precision, parameter :: tmax = 1.d0

  integer :: imin, imax

  integer :: i, n

  double precision :: time, dt

  integer, parameter :: init=2

  ! create the grid
  dx = (xmax - xmin)/dble(nx-1)

  ! x = 0 is the same as x = nx-1, so we don't loop over it
  imin = 1
  imax = nx-1

  print *,'limits: ', imin, imax

  do i = 0, ng+nx-1
     x(i) = i*dx + xmin
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
           u(i) = 0.d0
        endif
     endif

  enddo

  ! fill the boundary conditions -- just the innermost GC
  u(0) = u(nx-1)
  u(nx) = u(1)

  uinit(:) = u(:)

  time = 0.0

  ! compute the initial timestep
  dt = 1.e33
  do i = imin, imax
     dt = min(dt, dx/abs(a + small))
  enddo
  dt = cfl*dt
  print *, 'initial dt = ', dt
  
  do while (time < tmax)

     c = a*dt/dx

     ! evolve
     do i = imin, imax

        ! L-W
        unew(i) = 0.5*c*(1.0+c)*u(i-1) + (1.0-c**2)*u(i) - 0.5*c*(1.0-c)*u(i+1)

     enddo

     u(:) = unew(:)

     ! fill the boundary conditions -- just the innermost GC
     u(0) = u(nx-1)
     u(nx) = u(1)

     ! update the time
     time = time + dt

     ! compute the new timestep
     dt = 1.e33
     do i = imin, imax
        dt = min(dt, dx/abs(a + small))
     enddo
     dt = cfl*dt
     print *, time, dt

  enddo
     
  print *, '# time = ', time
  ! print out
  do i = 0, nx-1
     print *, x(i), uinit(i), u(i)
  enddo


end program advection
  
  
  
