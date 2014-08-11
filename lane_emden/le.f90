! integrate the Lane-Emden equation:
!
!  (1/x**2) (x**2 phi')' = -phi**n
!
! here, we define two variables,
!
!  f1 = phi
!  f2 = phi' (= d(phi)/dx)
!
! We'll use 2nd order RK
!
program le

  implicit none

  double precision :: f1, f1_new, f1_half, f2, f2_new, f2_half
  double precision :: xi, xi_half

  ! step size
  double precision :: h_0
  double precision :: h

  double precision :: tol = 1.e-12
  double precision :: R_est

  double precision :: n 

!  print *, 'enter step size: '
!  read *, h_0

!  print *, 'enter polytropic index, n: '
!  read *, n
  h_0 = 1.d-2
  n = 1.5


  ! set the boundary conditions
  f1 = 1.d0
  f2 = 0.d0

  xi = 0.d0

  ! set the initial stepsize.  We will use a constant stepsize until we
  ! approach the boundary of the star, where we will decrease it until
  ! the step size is smaller than the tolerance, tol.
  h = h_0

  ! integrate outward until the density (f1) goes negative
  do while (h > tol)

     f1_half = f1 + 0.5*h*f2

     if (xi == 0.d0) then

        ! use a Taylor expansion of phi about the origin (see Clayton
        ! p. 158)
        f2_half = f2 + 0.5*h*((2.d0/3.d0) - f1**n)
     else
        f2_half = f2 + 0.5*h*(-f1**n - 2.d0*f2/xi)
     endif

     xi_half = xi + 0.5*h

     f1_new = f1 + h*f2_half
     f2_new = f2 + h*(-f1_half**n - 2.d0*f2_half/xi_half)

     xi = xi + h
     f1 = f1_new
     f2 = f2_new

     ! set the new stepsize.  Since phi'' < 0, we are always convex,
     ! and the intersection of phi' with the x-axis will always be
     ! a conservative estimate of the radius of the star.  Make sure
     ! that the stepsize does not take us past that
     R_est = xi - f1/f2

     if (xi + h > R_est) then
        h = -f1/f2
     endif


     print *, xi, f1, f2, R_est, h


  enddo

end program le
