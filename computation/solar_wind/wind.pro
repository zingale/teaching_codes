pro wind

npts = 20

xi_min = 0.4
xi_max = 2

dxi = (xi_max - xi_min)/(npts -1)

xi = findgen(npts)*dxi + xi_min

print, xi



end
