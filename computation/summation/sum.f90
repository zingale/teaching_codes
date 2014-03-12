! demonstrate roundoff error by adding smaller and smaller numbers to
! one until we get 1

program summation

  real :: epsilon

  epsilon = 1.0e0

  do while (1.0e0 + epsilon /= 1.0e0)  
     epsilon = epsilon/2.0e0
  enddo

  print *, 'epsilon = ', epsilon

end program summation
