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
  
  integer :: i, n, nsteps, nsteps0

  ! set the initial conditions

  
  ! set the initial timestep
  nsteps0 = 1

  do i = 1, 12

     x_0 = 1.d0
     y_0 = 0.d0
     
     vx_0 = 0.d0
     vy_0 = 6.283185307d0   ! 2 pi
     
     dt = 1.d0/nsteps0

     t = 0.d0
  
     ! compute the number of steps needed for 2 years of integration
     nsteps = 2*nsteps0

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

     enddo

     print *, real(dt), real(abs(sqrt(x**2 + y**2) - 1.d0)), real(sqrt((x-1.d0)**2 + (y-0.d0)**2))

     nsteps0 = 2*nsteps0
     
  enddo

end program orbit
