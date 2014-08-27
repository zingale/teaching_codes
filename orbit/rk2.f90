! compute the orbit of earth around the Sun using first-order accurate
! Euler's method

! here, we work in units of AU, yr, and M_sun, for convienence.

program orbit

  implicit none

  double precision, parameter :: G = 39.47757d0   ! AU^3 M_sun^{-1} yr^{-2}
  double precision, parameter :: M_sun = 1.d0

  double precision :: x_0, y_0, x, y, x_half, y_half
  double precision :: vx_0, vy_0, vx, vy, vx_half, vy_half

  double precision :: ax_0, ay_0, ax_half, ay_half
  double precision :: r
  
  double precision :: t, dt
  
  integer :: n, nsteps

  ! set the initial conditions
  x_0 = 1.d0
  y_0 = 0.d0

  vx_0 = 0.d0
  vy_0 = 6.283185307d0  ! 2 pi

  
  ! set the initial timestep
  print *, 'Enter number of steps per year: '
  read *, nsteps
  dt = 1.0d0/nsteps

  t = 0.d0
  
  ! compute the number of steps needed for 2 years of integration
  nsteps = 2.d0*nsteps

  print *, "# dt = ", dt
  print *, "# t, x, y, vx, vy"

  print *, real(t), real(x_0), real(y_0), real(vx_0), real(vy_0)
  
  do n = 1, nsteps

     r = sqrt(x_0**2 + y_0**2)
     ax_0 = -G*M_sun*x_0/r**3
     ay_0 = -G*M_sun*y_0/r**3

     x_half = x_0 + vx_0*dt/2.d0
     y_half = y_0 + vy_0*dt/2.d0

     vx_half = vx_0 + ax_0*dt/2.d0
     vy_half = vy_0 + ay_0*dt/2.d0

     r = sqrt(x_half**2 + y_half**2)
     ax_half = -G*M_sun*x_half/r**3
     ay_half = -G*M_sun*y_half/r**3

     x = x_0 + vx_half*dt
     y = y_0 + vy_half*dt

     vx = vx_0 + ax_half*dt
     vy = vy_0 + ay_half*dt


     x_0 = x
     y_0 = y
     vx_0 = vx
     vy_0 = vy
     
     t = t + dt

     print *, real(t), real(x), real(y), real(vx), real(vy), real(sqrt(x**2+y**2))

  enddo

  print *, 'dt, error r, error dx = ', dt, abs(sqrt(x**2 + y**2) - 1.d0), sqrt((x-1.d0)**2 + (y-0.d0)**2)


end program orbit
