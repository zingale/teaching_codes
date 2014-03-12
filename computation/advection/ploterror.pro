pro ploterror

loadct, 0
!p.charsize = 1.25
!p.background = 255
!p.color = 0

data = readtab('bw-0.8.out', /double)

;data = readtab('upwind-0.8.out', /double)

x = data[0,*]
u0 = data[1,*]
un = data[2,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.5,1.5], ystyle=1
oplot, x, u0, thick=3, linestyle=1

end
