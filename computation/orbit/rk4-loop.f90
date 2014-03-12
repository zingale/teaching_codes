! compute the orbit of earth around the Sun using first-order accurate
! Euler's method

! here, we work in units of AU, yr, and M_sun, for convienence.

program orbit

  implicit none

  double precision :: k1(4), k2(4), k3(4), k4(4)
  double precision :: y(4), f(4)
  double precision :: t, dt


  integer :: i, n, nsteps, nsteps0


  
  ! set the initial timestep
  nsteps0 = 1

  do i = 1, 20

     ! set the initial conditions
     y(1) = 1.d0   ! x
     y(2) = 0.d0   ! y
     
     y(3) = 0.d0      ! v_x
     y(4) = 6.283185307d0 ! v_y

     dt = 1.d0/nsteps0

     t = 0.d0
  
     ! compute the number of steps needed for 2 years of integration
     nsteps = 2*nsteps0
     
     do n = 1, nsteps
        
        call rhs(t, y, f)
        k1(:) = dt*f(:)
        
        call rhs(t+0.5*dt, y(:)+0.5*k1(:), f)
        k2(:) = dt*f(:)
        
        call rhs(t+0.5*dt, y(:)+0.5*k2(:), f)
        k3(:) = dt*f(:)
        
        call rhs(t+dt, y(:)+k3(:), f)
        k4(:) = dt*f(:)

        y(:) = y(:) + (1.0d0/6.d0)*(k1(:) + 2.d0*k2(:) + 2.d0*k3(:) + k4(:))
        
        t = t + dt

     enddo

     print *, real(dt), real(abs(sqrt(y(1)**2 + y(2)**2) - 1.d0)), real(sqrt((y(1)-1.d0)**2 + (y(2)-0.d0)**2))

     nsteps0 = 2*nsteps0

  enddo

end program orbit



subroutine rhs(t,y, f)

  implicit none

  ! in these units, G = v^2 = (2 pi)^2
  double precision, parameter :: G = 39.478417602

  double precision, parameter :: M_sun = 1.d0

  double precision :: t, y(4), f(4)

  double precision :: r

  ! y(1) = x, y(2) = y, y(3) = v_x, y(4) = v_y

  f(1) = y(3)
  f(2) = y(4)
  r = sqrt(y(1)**2 + y(2)**2)
  f(3) = -G*M_sun*y(1)/r**3
  f(4) = -G*M_sun*y(2)/r**3

  return

end subroutine rhs
