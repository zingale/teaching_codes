! integrate a white dwarf, assuming an equation of state 
! P = K rho^n
!
! where n = 5/3 for a non-relativistic degenerate gas
!
! We'll use 2nd order discretization of the HSE equation
!
! p_{i+1} = p_i + (1/2) (rho_i + rho_{i+1}) g dr
!
! information is stored at zone centers and the gravitational
! acceleration is computed on the interface between zones.


program le

  implicit none

  double precision :: dr

  double precision :: p, p_old
  double precision :: rho, rho_old

  double precision :: dpdrho, drho, p_want

  double precision :: rho_c

  double precision :: mass_enclosed, g_interface
  double precision :: r_left, r_right

  ! Newton's constant in CGS
  double precision, parameter :: Gconst = 6.67428d-8

  double precision, parameter :: pi = 3.14159265358979323846

  double precision, parameter :: M_solar = 2.d33

  logical :: converged, star_edge

  integer :: i

  integer :: iter
  integer, parameter :: MAX_ITERS = 100
  double precision, parameter :: TOL = 1.d-8


  ! set the WD model parameters -- all in CGS units

  ! central density
  rho_c = 1.d6

  ! resolution
  dr = 4.d6


  ! set the first zone (i = 1)
  i = 1
  rho_old = rho_c
  call eos(rho_old, p_old, dpdrho)


  ! output the first zone's information
  print *, (dble(i-1) + 0.5d0)*dr, rho_old, p_old

  mass_enclosed = (4.d0/3.d0)*pi*dr**3*rho_old

  ! integrate outward until the density goes negative
  star_edge = .false.
  i = 2
  do while (.not. star_edge)    

     ! find the physical coordinate of the interface between the
     ! current and previous zone
     r_left = dble(i-1)*dr

     ! the coordinate of the right edge of the zone
     r_right = dble(i)*dr

     ! find the gravitational acceleration at the interface
     g_interface = -Gconst*mass_enclosed/r_left**2


     ! guess at the density in the new cell
     rho = rho_old

     ! iterate to find the pressure and density for the new zone
     converged = .false.
     do iter = 1, MAX_ITERS

        ! find the pressure from HSE
        p_want = p_old + 0.5d0*(rho_old + rho)*g_interface*dr
        
        ! find the corresponding pressure from the EOS
        call eos(rho, p, dpdrho)

        ! Newton method to find the zero of (p_want - p)
        drho = - (p_want - p)/(0.5d0*g_interface*dr - dpdrho)

        ! restrict the change to be 10% -- this makes the Newton
        ! iterations better behaved and allows us to more reliably
        ! detect if we are at the edge of the star
        if (abs(drho) > 0.1d0*rho) then
           rho = rho + 0.1*drho
        else
           rho = rho + drho
        endif

        ! check for convergence
        if (abs(drho) < TOL*rho) then
           converged = .true.
           exit
        endif

        ! check to see if we are at the edge of the star
        if (rho < 0.d0) then
           star_edge = .true.
           rho = rho + drho
           exit
        endif
        
     enddo

     if (.not. converged .and. .not. star_edge) then
        print *, "ERROR: failure to converge"
        stop
     endif

     ! output
     if (rho > 0.d0) then 
        print *, (dble(i-1) + 0.5d0)*dr, rho, p

        ! store the information for the next zone
        rho_old = rho
        p_old = p

        ! compute the new mass enclosed
        mass_enclosed = mass_enclosed + (4.d0/3.d0)*pi*dr* &
             (r_left**2 + r_left*r_right + r_right**2)*rho_old
     
        i = i+1
     endif

  enddo


  print *, "WD mass = ", mass_enclosed/M_solar


end program le

subroutine eos(rho, p, dpdrho) 

  implicit none

  double precision, intent(in) :: rho
  double precision, intent(out) :: p
  double precision, intent(out) :: dpdrho

  ! parameters for a non-relativistic degenerate gas, in CGS
  double precision, parameter :: K = 1.d13
  double precision :: n = 5./3.
  double precision :: Z_over_A = 0.5d0

  p = K*(Z_over_A*rho)**n
  dpdrho = n*K*(Z_over_A)**n*rho**(n-1.0d0)

end subroutine eos
