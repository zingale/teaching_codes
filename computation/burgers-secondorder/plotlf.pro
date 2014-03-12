pro plotlf

loadct, 0
!p.charsize = 1.25
!p.background = 255
!p.color = 0

data = readtab('burgers.out', /double)

x = data[0,*]
u0 = data[1,*]
un = data[2,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-1.25,1.25], ystyle=1
oplot, x, u0, thick=3, linestyle=1
end
