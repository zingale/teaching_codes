
loadct, 0
!p.charsize = 1.75
!p.background = 255
!p.color = 0
!x.thick=2
!y.thick=2
!p.thick=2

data = readtab('convergence.out', /double)

n = data[0,*]
dx = 1.0/n
e2 = data[1,*]

plot, dx, e2, /xlog, /ylog, xtitle="!4D!3x", ytitle="error", thick=3, xstyle=1, psym=4
oplot, dx, e2[0]*(dx/dx[0])^(2.0), linestyle=2, thick=2


end
