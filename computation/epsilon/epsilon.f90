program epsilon_test

  implicit none

  double precision :: eps


  eps = 1.d0

  do while (1.d0 + eps /= 1.d0)
     eps = eps/2.d0
  enddo

  print *, eps, epsilon(eps)

end program epsilon_test
  
