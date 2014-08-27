program wind

  implicit none

  integer, parameter :: npts = 20
  
  double precision, parameter :: xi_min = 0.4
  double precision, parameter :: xi_max = 2.0

  double precision, parameter :: C = -3.0

  double precision :: xi, dxi
  
  double precision :: M, dM
  double precision :: f, fprime

  integer, parameter :: MAX_ITER = 100
  double precision, parameter :: tol = 1.d-4

  logical :: converged

  integer :: i, iter

  dxi = (xi_max - xi_min)/(npts - 1)
  
  xi = xi_min

  do i = 1, npts


     ! initial guess at M
     M = 0.2d0

     converged = .false.

     do iter = 1, MAX_ITER
        
        f = M**2 - log(M**2) - 4.d0*log(xi) - 4.d0/xi - C
        fprime = 2.d0*M - 1.d0/M**2

        dM = -f/fprime

        M = min(1.1*M, max(M + dM, 0.9*M))


        if (abs(dM) < tol*M) then
           print *, xi, M
           converged = .true.
           exit
        endif

     enddo

     if (.not. converged) then
        print *, 'ERROR: did not converge for xi = ', xi, dM
     endif


     xi = xi + dxi

  enddo

end program wind



  
