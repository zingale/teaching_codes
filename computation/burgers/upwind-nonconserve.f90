! solve the inviscid Burger's equation using conservative upwinding.
!
! here we assume that u > 0 always
!
! u_t + u u_x = 0

program upwind

  implicit none

  integer, parameter :: nx = 128
  integer, parameter :: ng = 1

  double precision, dimension(2*ng + nx) :: x, u, unew, uinit

  double precision, parameter :: xmin = 0.0
  double precision, parameter :: xmax = 1.0
  double precision :: dx

  double precision :: time

  ! here, cfl is dt/dx
  double precision, parameter :: cfl = 0.2

  double precision, parameter :: pi = 3.14159

  integer :: nsteps

  integer :: imin, imax

  integer :: i, n
  
  integer, parameter :: init=2

  character (len=4) :: filenum

  ! create the grid
  dx = (xmax - xmin)/dble(nx)

  imin = ng+1
  imax = ng+nx

  do i = 1, 2*ng+nx
     x(i) = (i-ng-0.5d0)*dx + xmin
  enddo

  ! compute the number of steps needed to advect once
  nsteps = (xmax-xmin)/(cfl*dx)

  ! initialize
  do i = imin, imax
     
     if (init == 1) then
        ! sin wave
        u(i) = sin(2.d0*pi*x(i)) + 1.5

     else if (init == 2) then
        ! square wave
        if (x(i) > 0.5d0) then
           u(i) = 1.d0
        else
           u(i) = 2.d0
        endif
     endif

  enddo

  uinit(:) = u(:)
  
  time = 0.d0


  open (unit = 10, file="burger-nc.0000", status="unknown")
  write (10,*) "#time = ", time
  do i = imin, imax
     write (10,*) x(i), u(i)
  enddo

  close (10)


  do n = 1, nsteps

     ! fill the boundary conditions -- outflow
     u(imin-1) = u(imin)
     u(imax+1) = u(imax)

     ! evolve
     do i = imin, imax

        ! FTBS (upwind)
        unew(i) = u(i) - cfl*u(i)*(u(i) - u(i-1))


     enddo

     time = time + cfl*dx
     
     if (mod(n,50) == 0) then
        write (filenum, '(i4.4)') n
        open (unit = 10, file="burger-nc."//filenum, status="unknown")
        write (10,*) "# time = ", time
         do i = imin, imax
           write (10,*) x(i), u(i)
        enddo

        close (10)
     endif

     u(:) = unew(:)
     
  enddo
     
  ! print out
  do i = imin, imax
     print *, i, uinit(i), u(i)
  enddo


end program upwind
  
  
  
