!  Solve riemann shock tube problem for a general equation of state using 
!  the method of Colella and Glaz.  Use a two shock approximation, and
!  linearly interpolation between the head and tail of a rarefaction to
!  treat rarefactions.
!
!  The Riemann problem for the Euler's equation produces 4 states, 
!  separated by the three characteristics (u - cs, u, u + cs):
!
!
!        l_1      t    l_2       l_3
!         \       ^   .       /
!          \  *L  |   . *R   /
!           \     |  .     /
!            \    |  .    /
!        L    \   | .   /    R
!              \  | .  /
!               \ |. / 
!                \|./
!       ----------+----------------> x
!
!       l_1 = u - cs   eigenvalue
!       l_2 = u        eigenvalue (contact)
!       l_3 = u + cs   eigenvalue
!
!       only density jumps across l_2
!
!  References:
!
!   CG:   Colella & Glaz 1985, JCP, 59, 264.
!
!   CW:   Colella & Woodward 1984, JCP, 54, 174.
!
!   Fry:  Fryxell et al. 2000, ApJS, 131, 273.
!
!   Toro: Toro 1999, ``Riemann Solvers and Numerical Methods for Fluid
!         Dynamcs: A Practical Introduction, 2nd Ed.'', Springer-Verlag
!

subroutine riemann(gamma, &
                   dens_l, xmom_l, ener_l, &
                   dens_r, xmom_r, ener_r, &
                   F_dens, F_xmom, F_ener)

  implicit none

  integer :: n

  double precision :: gamma
  double precision :: dens_l, xmom_l, ener_l
  double precision :: dens_r, xmom_r, ener_r
  double precision :: F_dens, F_xmom, F_ener

  double precision :: riemann_tol;
  integer :: nriem;

  double precision :: scratch, scratch2;
  double precision :: ustar_sgn;

  double precision :: vlft, ulft, plft, clft
  double precision :: vrght, urght, prght, crght

  double precision :: pstar, pstar1, pstar2
  double precision :: wlft, wlft1
  double precision :: wrght, wrght1

  double precision :: wes, westar;
  double precision :: rhos, rhostr, us, ustar, ps;

  double precision :: smlrho, smallp, smallu;

  double precision :: ustrl1, ustrr1, ustrl2, ustrr2, delu1, delu2;

  double precision :: rhoav, uav, pav;

  double precision :: eps;

  double precision :: pres_err;
  logical :: has_converged;

  double precision :: vs, vstar, ces, cestar, ws;

  ! some parameters 
  riemann_tol = 1.d-5;
  nriem = 15;

  eps = 1.d-8;

  smlrho = 1.d-10;
  smallp = 1.d-10;
  smallu = 1.d-10;

  
  ! switch to primitive variables
  vlft = 1.d0/max(dens_l, smlrho)
  vrght = 1.d0/max(dens_r, smlrho)

  ulft = xmom_l/dens_l
  urght = xmom_r/dens_r

  plft = (ener_l - 0.5*xmom_l*xmom_l/dens_l)*(gamma - 1.0)
  plft = max(smallp, plft);

  prght = (ener_r - 0.5*xmom_r*xmom_r/dens_r)*(gamma - 1.0)
  prght = max(smallp, prght)

  clft = sqrt(gamma*plft*dens_l);
  crght = sqrt(gamma*prght*dens_r);


  ! construct first guess for secant iteration by assuming that the
  ! nonlinear wave speed is equal to the sound speed -- the resulting
  ! expression is the same as Toro, Eq. 9.28 in the Primitive Variable
  ! Riemann Solver (PVRS).  See also Fry Eq. 72.
  pstar1 = prght - plft - crght*(urght - ulft)
  pstar1 = plft + pstar1*(clft/(clft + crght))
  pstar1 = max(smallp, pstar1)


  ! calculate nonlinear wave speeds for the left and right moving
  ! waves based on the first guess for the pressure jump.  Again,
  ! there is a left and a right wave speed.  Compute this using CG
  ! Eq. 34.

  ! note -- we simplify a lot here, assuming constant gamma 
  wlft1 = pstar1 + 0.5*(gamma - 1.0)*(pstar1 + plft)
  wlft1 = sqrt(dens_l*abs(wlft1))
       
  wrght1 = pstar1 + 0.5*(gamma - 1.0)*(pstar1 + prght)
  wrght1 = sqrt(dens_r*abs(wrght1))


  ! construct second guess for the pressure using the nonlinear wave
  ! speeds from the first guess.  This is basically the same thing we
  ! did to get pstar1, except now we are using the better wave speeds
  ! instead of the sound speed.
  pstar2 = prght - plft - wrght1*(urght - ulft)
  pstar2 = plft + pstar2*wlft1/(wlft1 + wrght1)
  pstar2 = max(smallp, pstar2)
  

  ! begin the secant iteration -- see CG Eq. 17 for details.  We will
  ! continue to interate for convergence until the error falls below
  ! tol (in which case, things are good), or we hit nriem iterations
  ! (in which case we have a problem, and we spit out an error).
  has_converged = .false.

  do n = 1, nriem
          
     ! new nonlinear wave speeds, using CG Eq. 34 
     wlft = pstar2 + 0.5*(gamma  - 1.d0)*(pstar2 + plft)
     wlft = sqrt(dens_l*abs(wlft))

     wrght = pstar2 + 0.5*(gamma - 1.d0)*(pstar2 + prght)
     wrght = sqrt(dens_r*abs(wrght))
          

     ! compute the velocities in the "star" state -- using CG
     ! Eq. 18 -- ustrl2 and ustrr2 are the velocities they define
     ! there.  ustrl1 and ustrl2 seem to be the velocities at the
     ! last time, since pstar1 is the old 'star' pressure, and
     ! wlft1 is the old wave speed.
     ustrl1 = ulft - (pstar1 - plft)/wlft1
     ustrr1 = urght + (pstar1 - prght)/wrght1
     ustrl2 = ulft - (pstar2 - plft)/wlft
     ustrr2 = urght + (pstar2 - prght)/wrght
          
     delu1 = ustrl1 - ustrr1
     delu2 = ustrl2 - ustrr2

     scratch = delu2  - delu1
          
     if (abs(pstar2 - pstar1) <= smallp) scratch = 0.d0;

     if (abs(scratch) < smallu) then
        delu2 = 0.d0
        scratch = 1.d0
     endif

     ! pressure at the "star" state -- using CG Eq. 18
     pstar = pstar2 - delu2*(pstar2 - pstar1)/scratch
     pstar = max(smallp, pstar)

     ! check for convergence of iteration
     pres_err = abs(pstar - pstar2)/pstar
     if (pres_err < riemann_tol) then
        has_converged = .true.
        exit
     endif
          
     ! reset variables for next iteration 
     pstar1 = pstar2
     pstar2 = pstar

     wlft1 = wlft
     wrght1 = wrght
          
  enddo

  if (.not. has_converged) then
       
     print *, " "
     print *, "Nonconvergence in subroutine rieman!"
     print *, "Pressure error = ", pres_err
     print *, "pL = ", plft, " pR = ", prght
     print *, "uL = ", ulft, " uR = ", urght
     print *, "cL = ", clft, " crght = ", crght
     print *, "Terminating execution"
     stop

  endif

  ! end of secant iteration 

  ! calculate fluid velocity for the "star" state -- this comes from
  ! the shock jump equations, Fry Eq. 68 and 69.  The ustar velocity
  ! can be computed using either the jump eq. for a left moving or
  ! right moving shock -- we use the average of the two.

  scratch = ulft - (pstar - plft)/wlft
  scratch2 = urght + (pstar - prght)/wrght
  ustar = 0.5*(scratch + scratch2)

  ustar_sgn = sign(1.d0, ustar)

  ! decide which state is located at the zone iterface based on
  ! the values of the wave speeds.  This is just saying that if
  ! ustar > 0, then the state is U_L.  if ustar < 0, then the
  ! state on the axis is U_R.
  scratch = 0.50*(1.0 + ustar_sgn)
  scratch2 = 0.5e0*(1.0 - ustar_sgn)
       
  ps = plft*scratch + prght*scratch2
  us = ulft*scratch + urght*scratch2
  vs = vlft*scratch + vrght*scratch2

  rhos = 1.0/vs
  rhos = max(smlrho, rhos)

  vs = 1.0/rhos
  ws = wlft*scratch + wrght*scratch2
  ces = sqrt(gamma*ps*vs)

  ! compute rhostar, using the shock jump condition (Fry Eq. 80)
  vstar = vs - (pstar - ps)/(ws*ws)
  rhostr = 1.0/ vstar
  cestar = sqrt(gamma*pstar*vstar)

  ! compute some factors, Fry Eq. 81 and 82 
  wes = ces - ustar_sgn*us
  westar = cestar - ustar_sgn*ustar

  scratch = ws*vs - ustar_sgn*us
       
  if (pstar - ps >= 0.d0) then
     wes = scratch
     westar = scratch
  endif

      
  ! compute correct state for rarefaction fan by linear interpolation 
  scratch = max(wes - westar, wes + westar)
  scratch = max(scratch, smallu)

  scratch = (wes + westar)/scratch
         
  scratch = 0.5*(1.0 + scratch)
  scratch2 = 1.0 - scratch
         
  rhoav = scratch*rhostr + scratch2*rhos
  uav = scratch*ustar + scratch2*us
  pav = scratch*pstar + scratch2*ps

  if (westar >= 0.0) then
     rhoav  = rhostr
     uav = ustar
     pav = pstar
  endif
       
  if (wes < 0.0) then
     rhoav = rhos
     uav = us
     pav = ps
  endif

      
  ! now compute the fluxes 
  F_dens = rhoav*uav
  F_xmom = rhoav*uav*uav + pav
  F_ener = uav*(pav/(gamma - 1.0) + 0.5*rhoav*uav*uav + pav);

  return
end subroutine riemann
