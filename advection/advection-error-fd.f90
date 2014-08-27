! solve the linear advection equation via. first order upwinding
!
! u_t + a u_x = 0

program upwind

  implicit none

  integer, parameter :: nx_max = 1024
  integer, parameter :: ng = 1

  integer :: nx
  
  character (len=4) :: ifile

  double precision, dimension(ng + nx_max) :: x, u, unew, uinit, e
  double precision :: e1, e2, einfty

  double precision, parameter :: xmin = 0.0
  double precision, parameter :: xmax = 1.0
  double precision :: dx

  ! in the linear advection equation, the courant number
  ! is a dt/dx, so we don't need to specify a or dt
  double precision, parameter :: cfl = 0.8

  double precision, parameter :: pi = 3.14159

  integer :: nsteps, ngrids

  integer :: imin, imax

  integer :: i, n, m
  
  integer, parameter :: init=1


  ngrids = int(log(dble(nx_max))/log(2.d0)) - 1
  print *, 'ngrids = ', ngrids

  nx = 4

  do m = 1, ngrids


     ! create the grid
     dx = (xmax - xmin)/dble(nx-2)

     imin = ng+1
     imax = ng+nx-1

     do i = 1, ng+nx
        x(i) = (i-ng-1.d0)*dx + xmin
     enddo

     ! compute the number of steps needed to advect once
     nsteps = (xmax-xmin)/(cfl*dx)

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
        u(imin-1) = u(imax)
        u(imax+1) = u(imin)

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
           unew(i) = u(i) - 0.5*cfl*(u(i+1) - u(i-1)) + 0.5*cfl*cfl*(u(i+1) - 2.d0*u(i) + u(i-1))

        enddo

        u(:) = unew(:)
     
     enddo
     
     ! compute the error
     do i = imin, imax
        e(i) = u(i) - uinit(i)
     enddo

     ! compute the norms of the error
     e1 = 0.d0
     e2 = 0.d0
     einfty = -1.e33

     do i = imin, imax
        e1 = e1 + abs(e(i))
        e2 = e2 + e(i)*e(i)
        einfty = max(einfty, e(i))
     enddo
     
     e1 = e1/nx
     e2 = sqrt(e2/nx)
     einfty = einfty/nx

     write(ifile,'(i4.4)'), nx
     open (unit=10, file="advect."//ifile, status="unknown")
     do i = imin, imax
        write (10,*), x(i), u(i)
     enddo
     close(unit=10)

     print *, nx, e1, e2, einfty

     nx = 2*nx
  enddo

end program upwind
  
  
  
