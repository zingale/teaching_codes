! integrate the Lane-Emden equation:
!
!  (1/x**2) (x**2 phi')' = -phi**n
!
! here, we define two variables,
!
!  y1 = phi
!  y2 = phi' (= d(phi)/dx)
!
! We'll use 4th order RK
!
program le

  implicit none

  double precision :: k1(2), k2(2), k3(2), k4(2)
  double precision :: y(2), f(2)
  double precision :: xi

  ! step size
  double precision :: h_0
  double precision :: h

  double precision :: tol = 1.e-12
  double precision :: R_est

  integer :: i

  double precision, parameter :: n = 3.0


  h_0 = 1.0

  do i = 1, 15

     ! set the boundary conditions
     y(1) = 1.d0
     y(2) = 0.d0

     xi = 0.d0

     ! set the initial stepsize.  We will use a constant stepsize until we
     ! approach the boundary of the star, where we will decrease it until
     ! the step size is smaller than the tolerance, tol.
     h = h_0

     ! integrate outward until the density (f1) goes negative
     do while (h > tol)
        
        call rhs(n, xi, y, f)
        k1(:) = h*f(:)
        
        call rhs(n, xi+0.5*h, y(:)+0.5*k1(:), f)
        k2(:) = h*f(:)
        
        call rhs(n, xi+0.5*h, y(:)+0.5*k2(:), f)
        k3(:) = h*f(:)
        
        call rhs(n, xi+h, y(:)+k3(:), f)
        k4(:) = h*f(:)
        
        y(:) = y(:) + (1.0d0/6.d0)*(k1(:) + 2.d0*k2(:) + 2.d0*k3(:) + k4(:))
        
        xi = xi + h
        

        ! set the new stepsize.  Since phi'' < 0, we are always convex,
        ! and the intersection of phi' with the x-axis will always be
        ! a conservative estimate of the radius of the star.  Make sure
        ! that the stepsize does not take us past that
        R_est = xi - y(1)/y(2)
        
        if (xi + h > R_est) then
           h = -y(1)/y(2)
        endif

     enddo

     print *, h_0, xi

     h_0 = h_0/2.0

  enddo
end program le



subroutine rhs(n, xi ,y, f)

  implicit none

  double precision :: n
  double precision :: xi, y(4), f(4)

  f(1) = y(2)

  if (xi == 0.0) then
     f(2) = (2.d0/3.d0) - y(1)**n
  else
     f(2) = -2.d0*y(2)/xi - y(1)**n
  endif

  return
end subroutine rhs
