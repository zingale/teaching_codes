program diffuse

  ! implement 1-d implicit (backward-Euler) diffusion of a Gaussian
  ! with Neumann BCs on a cell-centered grid.
  !
  ! For the cell-centered grid variables, the overall grid looks like this:
  !
  !    
  !    |     |      |     X     |     |      |     |     X     |      |     |
  !    +--*--+- // -+--*--X--*--+--*--+- // -+--*--+--*--X--*--+- // -+--*--+
  !      -ng          -1     0     1           ...  nx-1    nx        nx+ng-1
  !                        (lo)                     (hi)    
  !    |                  |                              |                  |
  !    |<- ng ghostcells->|<---- nx interior zones ----->|<- ng ghostcells->|
  !    |                  |                              |                  |
  !
  ! Our discretization is:
  !
  !     n+1       n           n+1           n+1      n+1
  !  phi_i   - phi_i       phi_{i+1} - 2 phi_i  + phi_{i-1}
  !  --------------- =  k ---------------------------------
  !        dt                            dx**2
  !
  ! defining a = k dt/dx**2, this is:
  !
  !       n+1                 n+1        n+1         n
  ! -a phi_{i+1} + (1 + 2a) phi_i - a phi_{i-1} = phi_i
  !
  ! Neumann BCs are done as
  !
  !  phi_{lo-1} = phi_{lo}
  !  phi_{hi+1} = phi_{hi}
  !
  ! so at the left boundary, our equation is:
  !
  !            n+1          n+1          n
  ! (1 + a) phi_{lo} - a phi_{lo+1} = phi_{lo}
  !
  ! and at the right boundary, it is:
  !
  !        n+1                  n+1
  ! - a phi_{hi-1} + (1 + a) phi_{hi} = phi_{hi}
  !
  ! This can be expressed as an nx x nx matrix A that is tridiagonal.
  !

  implicit none

  integer :: nx

  double precision :: cfl
  double precision :: t, dt, tmax

  double precision :: xmin, xmax, xcenter, dx
  integer :: lo, hi

  double precision :: k, alpha

  double precision :: phi_0, phi_max, t_0
  double precision :: phi_true, phi_exact

  double precision, allocatable :: a(:), b(:), c(:), f(:)
  double precision, allocatable :: phi_new(:), phi_old(:)
  double precision, allocatable :: x(:)

  
  integer :: i

  ! set parameters
  nx = 128

  cfl = 1.0d0
  tmax = 0.01d0

  k = 1.0d0

  xmin = 0.0d0
  xmax = 1.0d0


  ! problem initial condition parameters
  t_0 = 0.001
  phi_0   = 1.0
  phi_max = 2.0


  ! initialize the grid properties
  dx = (xmax - xmin)/nx
  xcenter = 0.5d0*(xmin + xmax)

  lo = 0
  hi = nx-1

  allocate(x(lo:hi))
  do i = lo, hi
     x(i) = xmin + dble(i + 0.5d0)*dx
  enddo


  ! initialize a Gaussian
  allocate(phi_old(lo:hi))
  allocate(phi_new(lo:hi))

  do i = lo, hi
     phi_old(i) = phi_exact(x(i), xcenter, 0.0d0, t_0, k, phi_0, phi_max)
  enddo


  ! allocate array and RHS
  allocate(a(lo:hi), b(lo:hi), c(lo:hi), f(lo:hi))


  ! evolution loop
  t = 0.0d0
  do while (t < tmax)

     ! get timestep
     dt = cfl * 0.5d0*dx**2/k
     if (t + dt > tmax) dt = tmax - t

     alpha = k*dt/dx**2

     ! setup the matrix and RHS
     a(:) = 0.0d0
     b(:) = 0.0d0
     c(:) = 0.0d0

     ! left boundary
     b(lo) = 1.0d0 + alpha
     c(lo) = -alpha
     
     ! interior
     do i = lo+1, hi-1
        a(i) = -alpha
        b(i) = 1.0d0 + 2.0d0*alpha
        c(i) = -alpha
     enddo

     ! right boundary
     a(hi) = -alpha
     b(hi) = 1.0d0 + alpha

     
     ! RHS
     f(:) = phi_old(:)


     ! solve A phi = b
     call tridiag(a, b, c, f, phi_new, nx)

     phi_old(:) = phi_new(:)

     t = t + dt

  enddo

  ! output -- compare with the analytic solution
  do i = lo, hi
     phi_true = phi_exact(x(i), xcenter, t, t_0, k, phi_0, phi_max)
     print *, x(i), phi_true, phi_new(i)
  enddo

end program diffuse


function phi_exact(x, xcenter, t, t_0, k, phi_0, phi_max)

  implicit none

  double precision :: x, xcenter, t, t_0, k, phi_0, phi_max
  double precision :: phi_exact

  phi_exact = (phi_max - phi_0)*sqrt(t_0/(t + t_0))* &
       exp(-0.25d0*(x-xcenter)**2/(k*(t + t_0))) + phi_0

  return
end function phi_exact
  
  
subroutine tridiag(a,b,c,r,u,n)

  ! tridiagonal solver (from BoxLib and/or NR)
  !
  ! a_i u_{i-1} + b_i u_i + c_i x_{i+1} = r_i
  !

  implicit none

  integer         , intent(in   ) ::  n
  double precision, intent(in   ) :: a(n), b(n), c(n), r(n)
  double precision, intent(inout) :: u(n)

  integer j
  double precision, allocatable :: gam(:)
  double precision :: bet

  allocate(gam(n))

  if (b(1) == 0.0d0) then
     print *, 'tridiag: CANT HAVE B(1) = ZERO'
     stop
  endif

  bet = b(1)

  print *, 'in tridiag: ', bet

  u(1) = r(1)/bet

  do j = 2,n
     gam(j) = c(j-1)/bet
     bet = b(j) - a(j)*gam(j)

     if (bet == 0.0) then
        print *, 'tridiag: TRIDIAG FAILED'
        stop
     endif

     u(j) = (r(j)-a(j)*u(j-1))/bet
  end do

  do j = n-1,1,-1
     u(j) = u(j) - gam(j+1)*u(j+1)
  end do

end subroutine tridiag

