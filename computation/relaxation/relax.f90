! Do brute force relaxation on u'' = g(x), u(0) = 0, u(1) = 1, with 
! g(x) = cos(x)

! M. Zingale (2005-03-29)

program fmgfv

  implicit none

  integer, parameter :: nx = 256          ! number of interior zones
  integer, parameter :: ng = 1           ! number of guardcells

  integer, parameter :: nsmooth = 100000   ! number of smoothing blocks

  ! set the left and right boundary conditions
  double precision, parameter :: lbc = 0.d0
  double precision, parameter :: rbc = 1.d0

  ! imin and imax will always point to the starting and ending index of 
  ! the interior zones on the current level
  integer :: imin, imax                   

  double precision :: xmin, xmax, dx
  double precision :: source_norm

  integer :: i, j, n
  
  double precision, dimension(:), allocatable :: v, f
  double precision, dimension(:), allocatable :: w

  ! the function g holds the RHS, the function true holds the analytic 
  ! solution.
  double precision :: g, true, error
  double precision :: temp, temp2

  
  ! initialize the solution and rhs arrays
  allocate(f(nx+2*ng))
  allocate(v(nx+2*ng))
  allocate(w(nx+2*ng))

  f(:) = 0.d0
  v(:) = 0.d0

  
  ! setup the grid
  xmin = 0.d0
  xmax = 1.d0
  dx = (xmax - xmin)/dble(nx)

  imin = ng + 1           ! we are using 1 based indexing
  imax = ng + nx 


  ! fill the RHS of the finest level the true RHS
  do i = imin, imax
     f(i) = g(dble(i - ng - 1 + 0.5d0)*dx + xmin)
  enddo


  ! compute the source norm -- we will use this for error estimating
  source_norm = error(nx, ng, dx, f(:))
  print *, 'Source norm = ', source_norm

  open (unit =10, file = "relax.out", status= "replace")
  write (10,*) '# cycle,   L2 error,   L-infinity error'


  ! relax
  do j = 1, nsmooth
     call smooth(nx, ng, dx, lbc, rbc, v, f, 1)


     ! compare to the true solution
     w(:) = 0.d0
     do i = imin, imax
        w(i) = true(dble(i - ng - 1 + 0.5d0)*dx + xmin) - v(i)
     enddo

     ! print out every 10th iteration
     if (mod(j, 10) == 0) then
        temp = error(nx, ng, dx, w)
        write (10,*) j, temp, maxval(w)
     endif

  enddo

  close (unit=10)

100 format(1x, 1g13.6, 1x, 1g13.6, 1x, 1g13.6, 1x, 1g13.6)

  temp = error(nx, ng, dx, w)
  write (10,*) j, temp

  print *, 'after ', nsmooth, ' iterations, L2 norm of the error is ', temp

 
  
end program fmgfv
 


!=============================================================================
! g
!=============================================================================

function g(x)

  ! the RHS of the Poisson equation we are solving

  implicit none

  double precision g, x

  g = cos(x)

  return
end function g



!=============================================================================
! true
!=============================================================================

function true(x)

  ! the analytic solution to our equation

  implicit none

  double precision true, x
  double precision g

  true = -g(x) + x*g(1.d0) + 1.d0
  
  return
end function true



!=============================================================================
! error
!=============================================================================

function error(nx, ng, dx, v)

  ! compute the L2 norm 
  
  implicit none

  integer :: nx, ng
  double precision :: dx
  double precision, dimension(nx + 2*ng) :: v

  integer :: i, imin, imax
  double precision error

  imin = ng + 1
  imax = ng + nx

  error = 0.d0

  do i = imin, imax
     error = error + v(i)**2
  enddo
  
  error = dx*error    ! make it grid invariant
  error = sqrt(error)

  return
end function error



!=============================================================================
! smooth
!=============================================================================

subroutine smooth(nx, ng, dx, lbc, rbc, v, f, nsmooth)

  ! given a solution vector, v, and a RHS vector, f,
  ! smooth v to better satisfy the equation.  This is
  ! done in place, using Red-Black Gauss-Seidel

  ! lbc and rbc are the left and right Dirichlet boundary conditions.
  ! Because we are finite-volume, and therefore, cell-centered, we 
  ! need to extrapolate to match the desired Dirichlet BC.

  implicit none

  integer :: nx, ng, nsmooth
  double precision :: dx
  double precision :: lbc, rbc
  double precision, dimension(nx + 2*ng) :: v, f

  integer :: i, m, ioff, color
  integer :: imin, imax

  imin = ng + 1
  imax = ng + nx

  ! do some smoothing -- Red-Black Gauss-Seidel
  do m = 1, nsmooth

     ! set the guardcells to give the proper boundary condition, using
     ! extrapolation
     v(ng) = 2*lbc - v(imin)
     v(ng+nx+1) = 2*rbc - v(imax)

    ioff = 0
    do color = 0, 1
       do i = imin+ioff, imax, 2
          v(i) = 0.5d0*(v(i-1) + v(i+1) - dx*dx*f(i))
       enddo
       
       ioff = 1 - ioff
    enddo
       
 enddo

  return
end subroutine smooth



