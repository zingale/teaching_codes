pro plotlf

loadct, 0
!p.charsize = 1.25
!p.background = 255
!p.color = 0

d1 = 'conserved-shock-256.out'
d2 = 'conserved-shock-128.out'
d3 = 'conserved-shock-64.out'

data1 = readtab(d1, /double)
data2 = readtab(d2, /double)
data3 = readtab(d3, /double)


plot, data1[1,*], data[2, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.5,1.5], ystyle=1
oplot, x, u0, thick=3, linestyle=1
oplot, x, uold, thick=3, linestyle=0, color=128
end
