! Multigrid V-cycle Example
!
! Do multigrid V-cycles for the problem  u'' = g(x), u(0) = 0, u(1) = 1
! with g(x) = cos(x)
!
! Any other RHS can be done by changing the function g(x) below.  Also,
! for diagnostics, the analytic solution, computed in true(x) will need
! to be updated.
!
! We use a finite-volume discretiation.  The grid is required to have a 
! power-of-two number of interior zones + atleast 1 guardcell on each 
! end.
!
!  |     |      |     X     |     |      |     |     X     |      |     |
!  +-----+- // -+-----X-----+-----+- // -+-----+-----X-----+- // -+-----+
!     1           ng    ng+1  ng+2         ...  ng+nx ng+nx+1      2ng+nx
!  
!                       
!  |<- ng guardcells->|<---- nx interior zones ----->|<- ng guardcells->|
!
!
! multigrid requires a lot of operations that go between meshes of different
! sizes.  Here, we always coarsen by a factor of two.
!
! The coarse mesh corresponding to above is
!
!
!   ... -+------------X-----------+- // -+-----------X------------+- ...
!             ng           ng+1             ng+nx       ng+nx+1
!
! Given a zone i on the fine mesh, the corresponding coarse mesh zone 
! index is (i-ng-1)/2 + ng + 1 
!
! Given a zone i on the coarse mesh, the two corresponding fine mesh 
! zone indices are 2*(i-ng-1) + ng + 1 and 2*(i-ng-1) + ng + 2
!
! The restriction process is a simple averaging of the fine grid zones
! into the coarse zone that contains them.  The prolongation uses a 
! linear reconstruction of the coarse data and integrates this over
! the corresponding child zones.  Both procedures are conservative.
!
! We use a standard 3-point stencil:
!
!  u_{i+1} - 2 u_{i} + u_{i-1}
! ----------------------------  = f_{i}
!            dx**2
!
! This is second order accurate.   
!
! Special care needs to be taken enforcing the Dirichlet boundary
! conditions, since we are using a finite volume discretization.
! The boundary conditions need to be enforced on the boundary of the
! control volume, which means that we need to extrapolate to fill
! the guardcell.
!
!              X
!      |       X       |       |
!      +-------X-------+-------+- ...
!         ng      ng+1   ng+2
!
!
! The boundary condition needs to be satisfied at ng + 1/2.
! Through extrapolation, we have
!
!   u_{ng+1/2} - u_{ng+1} = (u_{ng+1} - u_{ng})/dx  * (xmin - (xmin + dx/2))
!
! Here we have taken the function averages at the cell centers.
! Then we can solve for the guardcell value:
!
!   u_ng = 2 u_{ng+1/2} - u_{ng+1}
!
! Similarly for the right boundary.
!
!
!
! Outline:
!
! This is just simple V-cycles.  The problem is smoothed on the fine
! mesh and then the residual is computed and restricted to a coarser
! mesh where the error equation is smoothed and then its residual is
! restricted, etc., until we get to a single zone, where we solve
! analytically.  The error is then prolonged up the 'V' and the solution
! on that level is corrected, etc.
!
! grid
!
! 5 x               x               x
!    \             / .             .
! 4   o           o   o           o
!      \         /     .         .
! 3     o       o       o       o
!        \     /         .     .
! 2       o   o           o   o
!          \ /             . .
! 1         x               x
!
!
!   |<--initial V-->|<- subsequent V-cycles ->
!
! In reality, we do not know the true solution to the problem at hand,
! so we compute the ratio of the residual error to the source norm
! and look for this to be smaller than the tolerance.
!
! The norm of the true error, the residual error, and the relative error
! of the solution is output in multigrid.out for each V-cycle.
!
! M. Zingale (2006-05-01)

program vcycle_test

  implicit none

  integer, parameter :: nx = 256            ! number of interior zones
  integer, parameter :: ng = 1              ! number of guardcells

  integer, parameter :: nsmooth = 3         ! number of smoothing iterations
  integer, parameter :: num_v_cycles = 100   ! maximum number of  V-cycles

  double precision, parameter :: tol = 1.d-10  ! tolerance for residual error

  double precision, parameter :: SMALL = 1.e-20 ! used to prevent / by 0

  ! set the left and right boundary conditions
  double precision, parameter :: lbc = 0.d0
  double precision, parameter :: rbc = 1.d0

  ! imin and imax will always point to the starting and ending index of 
  ! the interior zones on the current level
  integer :: imin, imax                   

  double precision :: xmin, xmax, dx
  double precision :: source_norm

  integer :: nlevels
  integer :: i, j,  n
  
  double precision, dimension(:,:), allocatable :: v, f
  double precision, dimension(:), allocatable :: w, vold

  ! the function g holds the RHS, the function true holds the analytic 
  ! solution.
  double precision :: g, true, error

  double precision :: true_error, res_error, solution_error

  
  ! make sure that the number of interior zones is a power of 2 and
  ! compute the number of levels
  do n = 1, 100
     if (2**n == nx) then
        nlevels = n
        exit
     endif
  enddo

  ! we want level 1 to be nx = 1
  nlevels = nlevels + 1

  ! initialize the solution and rhs arrays
  allocate(f(nx+2*ng, nlevels))
  allocate(v(nx+2*ng, nlevels))
  allocate(w(nx+2*ng))
  allocate(vold(nx+2*ng))

  f(:,:) = 0.d0
  v(:,:) = 0.d0
  vold(:) = 0.d0
  
  ! setup the finest grid
  xmin = 0.d0
  xmax = 1.d0
  dx = (xmax - xmin)/dble(nx)

  imin = ng + 1           ! we are using 1 based indexing
  imax = ng + nx 


  ! fill the RHS of the finest level the true RHS
  do i = imin, imax
     f(i, nlevels) = g(dble(i - ng - 1 + 0.5d0)*dx + xmin)
  enddo


  ! compute the source norm -- we will use this for error estimating
  source_norm = error(nx, ng, dx, f(:,nlevels))
  print *, 'Source norm = ', source_norm

  open (unit=10, file="multigrid.out", status="replace")
  write (10,*) '# multigrid results, nx = ', nx
  write (10,*) '# cycle,   true error,   residual error,    solution error'

  solution_error = 1.d33
  res_error = 1.d33

  j = 1
  do while (j < num_v_cycles .and. res_error > tol)

     call vcycle(nx, ng, xmin, xmax, lbc, rbc, nlevels, nsmooth, &
          v(:,nlevels), f(:,nlevels))

     
     ! compare to the true solution
     w(:) = 0.d0
     do i = imin, imax
        w(i) = true(dble(i - ng - 1 + 0.5d0)*dx + xmin) - v(i,nlevels)
     enddo

     print *, ' '
     true_error = error(nx, ng, dx, w)
     print *, 'maximum true error = ', true_error


     ! compare the residual to the source norm
     w(:) = 0.d0
     call residual(nx, ng, dx, v(:,nlevels), f(:,nlevels), w(:))
  
     res_error = error(nx, ng, dx, w)/source_norm
     print *, 'ratio of residual error to source norm = ', res_error


     ! compare the current solution to the previous solution
     w(:) = (v(:,nlevels) - vold(:))/(v(:,nlevels) + SMALL)

     solution_error = error(nx, ng, dx, w)
     print *, 'solution error = ', solution_error

     vold(:) = v(:,nlevels)

     write (10,*) j, true_error, res_error, solution_error
     
     j = j + 1
  enddo

  close (unit=10)
  
end program vcycle_test
 


!=============================================================================
! vcycle
!=============================================================================

subroutine vcycle(nx, ng, xmin, xmax, lbc, rbc, nlevels, nsmooth, &
     v_pass, f_pass)

  implicit none

  integer :: nx, ng, nlevels, nsmooth
  double precision :: xmin, xmax

  double precision, dimension(nx + 2*ng) :: v_pass, f_pass
  
  double precision, dimension(nx+2*ng, nlevels) :: v, f
  double precision, dimension(nx+2*ng) :: w
  
  double precision :: dx
  double precision :: lbc, rbc, lbc_level, rbc_level


  integer :: n, i
  integer :: nx_level

  integer :: imin, imax

  double precision :: error

  ! initialize 
  f(:,:) = 0.d0
  f(:,nlevels) = f_pass(1:nx+2*ng)

  v(:,:) = 0.d0
  v(:,nlevels) = v_pass(1:nx+2*ng)

  dx = (xmax - xmin)/dble(nx)

  imin = ng + 1
  imax = ng + nx

  print *, ' '
  print *, '<<< beginning V-cycle >>>'


  !------------------------------------------------------------------------
  ! Descending part of V-cycle
  !------------------------------------------------------------------------
  nx_level = nx
     
  do n = nlevels, 2, -1
     
     print *, ' '
     print *, 'level = ', n, ' nx = ', nx_level, ' dx = ', dx

     ! compute the residual -- this is for debugging only!
     w(:) = 0.d0
     call residual(nx_level, ng, dx, v(:,n), f(:,n), w(:))

     print *, 'before G-S, residual L2 norm = ', error(nx_level, ng, dx, w)


     ! do some smoothing
     if (n == nlevels) then
        lbc_level = lbc
        rbc_level = rbc
     else
        lbc_level = 0.d0
        rbc_level = 0.d0
     endif

     call smooth(nx_level, ng, dx, lbc_level, rbc_level, &
          v(:,n), f(:,n), nsmooth)

     
     ! compute the residual
     w(:) = 0.d0
     call residual(nx_level, ng, dx, v(:,n), f(:,n), w(:))

     print *, 'after G-S, residual L2 norm = ', error(nx_level, ng, dx, w)
     

     ! restrict the residual onto the coarser level
     nx_level = nx_level/2
     imax = imin + nx_level - 1
     dx = dx*2.d0

     call restrict(nx_level, ng, w(:), f(:,n-1))
          
  enddo

  
  !------------------------------------------------------------------------
  ! solve the coarse problem exactly
  !------------------------------------------------------------------------

  print *, ' '
  print *, '<<< bottom solve >>>'

  ! we know the numerical solution exactly when nx = 1
  !v(imin,1) = -0.5d0*(f(imin,1)*dx*dx - v(imin-1,1) - v(imin+1,1))
  v(imin,1) = -0.25*f(imin,1)*dx*dx
  
  !------------------------------------------------------------------------
  ! Ascending part of the V-cycle
  !------------------------------------------------------------------------
  do n = 2, nlevels, 1
        
     nx_level = 2*nx_level
     imax = imin + nx_level - 1
     dx = 0.5d0*dx
        
     print *, ' '
     print *, 'level = ', n, ' nx = ', nx_level, ' dx = ', dx
        
     
     ! prolong the error up from the coarser grid -- just do linear
     ! reconstruction for now
     w(:) = 0.d0
     call prolong(nx_level, ng, dx, v(:,n-1), w(:))

     
     ! add the correction to the previous solution
     do i = imin, imax
        v(i,n) = v(i,n) + w(i)
     enddo
     
     
     ! compute the residual -- this is for debugging only
     w(:) = 0.d0
     call residual(nx_level, ng, dx, v(:,n), f(:,n), w(:))
     
     print *, 'before G-S, residual L2 norm = ', error(nx_level, ng, dx, w)


     ! do some smoothing -- here we need to be careful with the boundary
     ! conditions -- if we are solving the original problem, then
     ! we use the user defined boundary conditions.  Otherwise, we are 
     ! solving the equation for the error, in which case, the boundary
     ! conditions on the error are 0.
     if (n == nlevels) then
        lbc_level = lbc
        rbc_level = rbc
     else
        lbc_level = 0.d0
        rbc_level = 0.d0
     endif

     call smooth(nx_level, ng, dx, lbc_level, rbc_level, &
          v(:,n), f(:,n), nsmooth)

     
     ! compute the residual -- this is for debugging only
     w(:) = 0.d0
     call residual(nx_level, ng, dx, v(:,n), f(:,n), w(:))

     print *, 'after G-S, residual L2 norm = ', error(nx_level, ng, dx, w)
        
  enddo

  v_pass(1:nx+2*ng) = v(:,nlevels)

  return
end subroutine vcycle



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
! residual
!=============================================================================

subroutine residual(nx, ng, dx, v, f, r)
  
  ! given a guess at the solution, v, and the RHS vector, f,
  ! compute the residual, r

  implicit none

  integer :: nx, ng
  double precision :: dx
  double precision, dimension(nx + 2*ng) :: v, f, r

  integer :: i, imin, imax

  imin = ng + 1
  imax = ng + nx

  do i = imin, imax
     r(i) = f(i) - (v(i-1) - 2.d0*v(i) + v(i+1))/(dx*dx)
  enddo

  return
end subroutine residual



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
  double precision, dimension(nx + 2*ng) :: v, f
  double precision :: lbc, rbc

  integer :: i, m, ioff, color
  integer :: imin, imax

  imin = ng + 1
  imax = ng + nx

  ! do some smoothing -- Red-Black Gauss-Seidel
  do m = 1, nsmooth

     v(imin-1) = 2*lbc - v(imin)
     v(imax+1) = 2*rbc - v(imax)

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



!=============================================================================
! restrict
!=============================================================================

subroutine restrict(nx, ng, U_fine, U_coarse)

  ! restrict the data in U_fine by a factor of two, and store
  ! in the corresponding locations in U_coarse
  !
  ! nx and ng correspond to the coarse mesh.
  !
  ! the fine mesh has nx_fine = 2*nx

  implicit none

  integer :: nx, ng

  double precision, dimension(nx+2*ng) :: U_coarse
  double precision, dimension(2*nx+2*ng) :: U_fine

  integer :: i, imin, imax

  ! set the bounds of interior zones for the coarse grid
  imin = ng + 1
  imax = ng + nx

  do i = imin, imax
     U_coarse(i) = 0.5d0*(U_fine(2*(i-ng-1)+ng+1) + U_fine(2*(i-ng-1)+ng+2))
  enddo

  return
end subroutine restrict



!=============================================================================
! prolong
!=============================================================================

subroutine prolong(nx, ng, dx, U_coarse, U_fine)

  ! prolong the data U_coarse by a factor of two and store
  ! in the corresponding locations in U_fine.
  !
  ! nx, ng, and dx correspond to the fine mesh
  !
  ! the coarse mesh haas nx_coarse = nx/2
  !
  ! We do a linear reconstruction of the data in the coarse zones,
  ! using f(x) = mx + b, where <f> = b = coarse zone value.
  ! The slope is computed using the neighboring coarse zones.
  ! The fine data is the result of integrating over [-dx/2,0] and
  ! [0,dx/2], normalized by dx/2.

  implicit none

  integer :: nx, ng
  double precision :: dx

  double precision, dimension(nx/2 + 2*ng) :: U_coarse
  double precision, dimension(nx + 2*ng) :: U_fine

  integer :: i, imin, imax

  double precision :: slope, dx_coarse
  integer :: iparent

  ! set the bounds of the interior zones for the fine grid
  imin = ng + 1
  imax = ng + nx

  dx_coarse = 2.d0*dx

  do i = imin, imax, 2

     ! find the location of the coarse zone corresponding to fine
     ! zones i and i+1
     iparent = (i-ng-1)/2 + ng + 1

     ! compute the slope
     slope = (U_coarse(iparent+1) - U_coarse(iparent-1))/(2*dx_coarse)

     ! fill the two fine cells
     U_fine(i) =  -0.25d0*slope*dx_coarse + U_coarse(iparent)
     U_fine(i+1) = 0.25d0*slope*dx_coarse + U_coarse(iparent)

     ! alternately, do just direct injection
     ! U_fine(i) = U_coarse(iparent)
     ! U_fine(i+1) = U_coarse(iparent)

  enddo
  
  return
end subroutine prolong
