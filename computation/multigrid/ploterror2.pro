
loadct, 0
!p.charsize = 1.75
!p.background = 255
!p.color = 0
!x.thick=2
!y.thick=2
!p.thick=2

data256 = readtab('multigrid.out.256', /double)

n256 = data256[0,*]
e256 = data256[1,*]
r256 = data256[2,*]

plot, n256, e256, /ylog, xtitle="# of V-cycles", ytitle="error", thick=2, xstyle=1, linestyle=0, yrange=[1.e-12,1000], ystyle=1
oplot, n256, r256, thick=2, linestyle=2

xyouts, 5, 0.1, '|r|!d2'
xyouts, 5, 1.e-5, '|e|!d2'

end
