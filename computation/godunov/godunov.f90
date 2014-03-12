! so the Euler equations using Godunov's method


program godunov

  implicit none

  integer :: i

  integer, parameter :: nx = 64
  integer, parameter :: ng = 1

  integer :: nstep

  double precision, parameter :: gamma = 1.4d0
  double precision, parameter :: cfl = 0.8d0

  double precision, parameter :: xmin = 0.d0
  double precision, parameter :: xmax = 1.d0

  ! which of the test problems do we run -- this corresponds to 
  ! Toro Ch. 4
  integer, parameter :: itest = 4

  double precision, parameter :: tmax = 1.0


  double precision, dimension(2*ng+nx) :: xl, x, xr
  double precision :: dx

  double precision :: time, dt, dtdx

  double precision, dimension(2*ng+nx) :: dens, xmom, ener
  double precision, dimension(2*ng+nx) :: dens_l, xmom_l, ener_l
  double precision, dimension(2*ng+nx) :: dens_r, xmom_r, ener_r
  double precision, dimension(2*ng+nx) :: flux_dens, flux_xmom, flux_ener

  double precision :: timestep

  ! create the grid
  call grid(nx, ng, xmin, xmax, dx, xl, x, xr)

  ! initialize the grid
  call init(nx, ng, dx, xmin, xmax, x, itest, gamma, dens, xmom, ener)

  ! fill the boundary conditions
  call boundary(nx, ng, dens, xmom, ener)


  ! compute the initial timestep
  dt = timestep(nx, ng, dx, cfl, gamma, dens, xmom, ener)

  nstep = 0
  call output(nstep, nx, ng, x, gamma, dens, xmom, ener)

  ! evolution loop
  time = 0.0
  do while (time < tmax)

     ! compute the left and right states -- here we just do piecewise
     ! constant

     do i = ng+1, ng+nx+1
        dens_l(i) = dens(i-1)
        xmom_l(i) = xmom(i-1)
        ener_l(i) = ener(i-1)

        dens_r(i) = dens(i)
        xmom_r(i) = xmom(i)
        ener_r(i) = ener(i)
     enddo


     ! compute the fluxes through the cell interfaces
     do i = ng+1, ng+nx+1
        call riemann(gamma, &
             dens_l(i), xmom_l(i), ener_l(i), &
             dens_r(i), xmom_r(i), ener_r(i), &
             flux_dens(i), flux_xmom(i), flux_ener(i))
     enddo

     dtdx = dt/dx

     ! do the conservative updating
     do i = ng+1, ng+nx
        dens(i) = dens(i) + dtdx*(flux_dens(i) - flux_dens(i+1))
        xmom(i) = xmom(i) + dtdx*(flux_xmom(i) - flux_xmom(i+1))
        ener(i) = ener(i) + dtdx*(flux_ener(i) - flux_ener(i+1))
     enddo


     ! fill the boundary conditions
     call boundary(nx, ng, dens, xmom, ener)

     
     time = time + dt

     ! compute the new timestep
     dt = timestep(nx, ng, dx, cfl, gamma, dens, xmom, ener)

     if (time + dt > tmax) dt = tmax - time

     nstep = nstep + 1

  enddo

  call output(nstep, nx, ng, x, gamma, dens, xmom, ener)


end program godunov




subroutine output(nstep, nx, ng, x, gamma, dens, xmom, ener)

  implicit none

  integer :: nstep
  integer :: nx, ng
  double precision :: gamma
  double precision, dimension(2*ng+nx) :: x, dens, xmom, ener
  
  character (len = 4) :: istep

  integer :: i


  write (istep, '(i4.4)') nstep
  open (unit=10, file="godunov."//istep, status="unknown")

  do i = ng+1, ng+nx+1
     write (10,*) x(i), dens(i), xmom(i)/dens(i), &
          (ener(i)-0.5*xmom(i)*xmom(i)/dens(i))*(gamma-1.d0), &
          (ener(i)-0.5*xmom(i)*xmom(i)/dens(i))/dens(i)
  enddo

  return
end subroutine output



subroutine init(nx, ng, dx, xmin, xmax, x, itest, gamma, dens, xmom, ener)

  implicit none

  integer :: nx, ng

  integer :: itest


  double precision :: x_interface
  
  integer :: i

  double precision :: gamma
  double precision :: dx
  double precision, dimension(2*ng+nx) :: x
  double precision :: xmin, xmax

  double precision, dimension(2*ng+nx) :: dens, xmom, ener

  ! specify the initial conditions in terms of primitive variables
  double precision :: dens_l, dens_r, u_l, u_r, p_l, p_r

  if (itest == 1) then
     dens_l = 1.d0
     dens_r = 0.125d0

     u_l = 0.d0
     u_r = 0.d0

     p_l = 1.d0
     p_r = 0.1d0

  else if (itest == 2) then
     dens_l = 1.d0
     dens_r = 1.d0

     u_l = -2.d0
     u_r = 2.d0

     p_l = 0.4d0
     p_r = 0.4d0

  else if (itest == 3) then
     dens_l = 1.d0
     dens_r = 1.d0

     u_l = 0.d0
     u_r = 0.d0

     p_l = 1000.d0
     p_r = 0.1d0

  else if (itest == 4) then
     dens_l = 5.6698d0
     dens_r = 1.d0

     u_l = -1.4701d0
     u_r = -10.5d0

     p_l = 100.d0
     p_r = 1.d0

  endif


  x_interface = 0.5*(xmin + xmax)

  do i = ng+1, ng+nx

     if (x(i) < x_interface) then

        ! left state
        dens(i) = dens_l
        xmom(i) = dens_l*u_l
        ener(i) = 0.5*dens_l*u_l*u_l + p_l/(gamma - 1.0)

     else

        ! right state
        dens(i) = dens_r
        xmom(i) = dens_r*u_r
        ener(i) = 0.5*dens_r*u_r*u_r + p_r/(gamma - 1.0)

     endif

  enddo

  return

end subroutine init




subroutine grid(nx, ng, xmin, xmax, dx, xl, x, xr)

  ! create the grid.  Here, the number of zones, n, and guardcells, ng, as
  ! well as the domain extrema, xmin and xmax are input.  
  !
  ! dx, xl, x, and xr are output

  implicit none

  integer :: nx, ng
  double precision :: xmin, xmax
  
  double precision :: dx
  double precision, dimension(2*ng+nx) :: xl, x, xr

  integer :: i

  dx = (xmax - xmin)/nx

  do i = 1, 2*ng+nx
     xl(i) = (i-ng-1)*dx + xmin
     xr(i) = (i-ng)*dx + xmin
     x(i) = 0.5*(xl(i) + xr(i))
  enddo

  return
end subroutine grid




subroutine boundary(nx, ng, dens, xmom, ener)

  ! apply the boundary conditions -- here we just do outflow (zero-gradient)

  implicit none

  integer :: nx, ng
  double precision, dimension(2*ng+nx) :: dens, xmom, ener

  integer :: i

  ! left boundary
  do i = 1, ng
     dens(i) = dens(ng+1)
     xmom(i) = xmom(ng+1)
     ener(i) = ener(ng+1)
  enddo

  
  ! right boundary
  do i = ng+nx+1, 2*ng+nx
     dens(i) = dens(ng+nx)
     xmom(i) = xmom(ng+nx)
     ener(i) = ener(ng+nx)
  enddo

  return
end subroutine boundary

  


function timestep(nx, ng, dx, cfl, gamma, dens, xmom, ener)

  implicit none

  integer :: nx, ng
  double precision :: dx
  double precision :: cfl

  double precision :: gamma

  double precision, dimension(2*nx+ng) :: dens, xmom, ener

  double precision :: eint, p, cs, u
  double precision :: dt

  double precision :: timestep
  integer :: i

  dt = 1.e33

  do i = ng+1, ng+nx

     ! compute the sound speed for this zone
     eint = ener(i) - 0.5*xmom(i)*xmom(i)/dens(i)
     p = eint*(gamma - 1.0)

     cs = sqrt(gamma*p/dens(i))
     u = xmom(i)/dens(i)

     dt = min(dt, dx/(abs(u) + cs))

  enddo

  timestep = dt
  return 

end function timestep
