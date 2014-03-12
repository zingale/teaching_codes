! compute the derivative of sin(x) for several different difference
! methods.

program derivative

  implicit none

  double precision :: f, fprime
  double precision :: x_0, dx

  integer, parameter :: nsub = 5
  integer :: i

  double precision :: true

  double precision :: dplus, dminus, d0, d3

  x_0 = 1.d0
  dx = .1d0
  
  true = fprime(x_0)

  print *, '#   h             D+             D-             D0             D3'
  do i = 1, nsub

     print *, real(dx), real(dplus(x_0, dx)-true), real(dminus(x_0, dx)-true), real(d0(x_0, dx)-true), real(d3(x_0, dx)-true)
     
     dx = 0.5d0*dx
  
  enddo

end program derivative



! f is the function whose derivative we wish to approximate

! fprime is the analytic derivative of f, which will be used in computing
! the error of the approximation

function f(x)
  implicit none

  double precision f, x

  f = sin(X)
  return
end function f


function fprime(x)
  implicit none

  double precision fprime, x

  fprime = cos(x)
  return
end function fprime


! dplus is a first-order accurate forward difference derivative
function dplus(x, dx)

  implicit none

  double precision dplus, x, dx
  double precision f

  dplus = (f(x+dx) - f(x))/dx

  return
end function dplus


! dminus is a first-order accurate, backward difference derivative
function dminus(x, dx)

  implicit none

  double precision dminus, x, dx
  double precision f

  dminus = (f(x) - f(x-dx))/dx

  return
end function dminus


! d2 is a second-order accurate derivative
function d0(x, dx)

  implicit none

  double precision d0, x, dx
  double precision f

  d0 = (f(x+dx) - f(x-dx))/(2.d0*dx)

  return
end function d0


! d3 is a third-order accurate derivative
function d3(x, dx)

  implicit none

  double precision d3, x, dx
  double precision f

  d3 = (2.d0*f(x+dx) + 3.d0*f(x) - 6.d0*f(x-dx) + f(x-2.0*dx))/(6.d0*dx)

  return
end function d3


