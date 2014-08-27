program epsilon_test

  implicit none

  real :: eps


  eps = 1.0

  do while (1.0 + eps /= 1.0)
     eps = eps/2.0
  enddo

  print *, eps, epsilon(eps)

end program epsilon_test
  
