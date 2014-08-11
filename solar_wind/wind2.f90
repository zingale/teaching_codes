program wind

  implicit none

  integer, parameter :: npts = 20
  
  double precision, parameter :: xi_min = 0.4
  double precision, parameter :: xi_max = 2.0

  double precision, parameter :: C = -3.0

  double precision :: xi, dxi
  
 
  double precision :: f
  double precision :: M1, M2, Mnew
  double precision :: f1, f2, fnew

  integer, parameter :: MAX_ITER = 25
  double precision, parameter :: tol = 1.d-4

  logical :: converged

  integer :: i, iter

  dxi = (xi_max - xi_min)/(npts - 1)
  
  xi = xi_min

  do i = 1, npts


     M1 = 1.01
     M2 = 10.0

     f1 = f(M1, xi, C)
     f2 = f(M2, xi, C)

     do iter = 1, MAX_ITER
        
        Mnew = 0.5*(M1 + M2)
        fnew = f(Mnew, xi, C)

        if (f1 == 0.0 .or. f2 == 0.0 .or. fnew == 0.0) exit

        if ((fnew < 0.0 .and. f1 > 0.0) .or. (fnew > 0.0 .and. f1 < 0.0)) then 
           ! solution is between M1 and Mnew
           M2 = Mnew
           
        else if ((fnew < 0.0 .and. f2 > 0.0) .or. (fnew > 0.0 .and. f2 < 0.0)) then
           ! solution is between Mnew and M2
           M1 = Mnew

        else
           print *, 'no root', xi
           exit
        endif
        

     enddo

     print *, xi, Mnew, (M2-M1)/M1
     xi = xi + dxi

  enddo

end program wind


function f(M, xi, C)

  double precision f, M, xi, C

  f = M**2 - log(M**2) - 4.d0*log(xi) - 4.d0/xi - C

  return
end function f


  
