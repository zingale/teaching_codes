! solve the linear advection equation via. FTCS, on a finite-difference grid
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

  ! in the linear advection equation, the courant number
  ! is a dt/dx, so we don't need to specify a or dt
  double precision, parameter :: cfl = 0.8

  double precision, parameter :: pi = 3.141592653

  integer :: nsteps

  integer :: imin, imax

  integer :: i, n
  
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

  ! fill the boundary conditions -- just the innermost GC
  u(0) = u(nx-1)
  u(nx) = u(1)

  uinit(:) = u(:)
  
  do n = 1, nsteps

     ! evolve
     do i = imin, imax

        ! FTCS
!        unew(i) = u(i) - 0.5*cfl*(u(i+1) - u(i-1))
        unew(i) = u(i) - cfl*(u(i+1) - u(i))
        
     enddo

     u(:) = unew(:)

     ! fill the boundary conditions -- just the innermost GC
     u(0) = u(nx-1)
     u(nx) = u(1)
     
  enddo
     
  ! print out
  do i = 0, nx+ng-1
     print *, i, x(i), uinit(i), u(i)
  enddo


end program advection
  
  
  
