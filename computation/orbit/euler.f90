! compute the orbit of earth around the Sun using first-order accurate
! Euler's method

! here, we work in units of AU, yr, and M_sun, for convienence.

program orbit

  implicit none

  double precision, parameter :: G = 39.47757d0   ! AU^3 M_sun^{-1} yr^{-2}
  double precision, parameter :: M_sun = 1.d0

  double precision :: x_0, y_0, x, y
  double precision :: vx_0, vy_0, vx, vy
  double precision :: ax_0, ay_0
  double precision :: r
  
  double precision :: t, dt
  
  integer :: n, nsteps

  ! set the initial conditions
  x_0 = 1.d0
  y_0 = 0.d0

  vx_0 = 0.d0
  vy_0 = 6.283185307d0   ! 2 pi

  
  ! set the initial timestep
  print *, 'Enter number of timesteps per year: '
  read *, nsteps
  dt = 1.d0/nsteps

  t = 0.d0
  
  ! compute the number of steps needed for 2 years of integration
  nsteps = 2*nsteps

  print *, "# dt = ", dt
  print *, "# t, x, y, vx, vy"

  print *, real(t), real(x_0), real(y_0), real(vx_0), real(vy_0)
  
  do n = 1, nsteps

     r = sqrt(x_0**2 + y_0**2)
     ax_0 = -G*M_sun*x_0/r**3
     ay_0 = -G*M_sun*y_0/r**3

     x = x_0 + vx_0*dt
     y = y_0 + vy_0*dt

     vx = vx_0 + ax_0*dt
     vy = vy_0 + ay_0*dt

     x_0 = x
     y_0 = y
     vx_0 = vx
     vy_0 = vy
     
     t = t + dt

     print *, real(t), real(x), real(y), real(vx), real(vy)

  enddo

  print *, 'dt, error r, error dx = ', dt, abs(sqrt(x**2 + y**2) - 1.d0), sqrt((x-1.d0)**2 + (y-0.d0)**2)

end program orbit
