pro plotfromm

loadct, 0
!p.charsize = 1.25
!p.background = 255
!p.color = 0

data = readtab('fromm.out', /double)
data2 = readtab('lw.out', /double)
data3 = readtab('bw.out', /double)

x = data[0,*]
u0 = data[1,*]
un = data[2,*]

ulw = data2[2,*]
ubw = data3[2,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.5,1.5], ystyle=1
oplot, x, u0, thick=3, linestyle=1
oplot, x, ulw, thick=3, linestyle=0, color=128
oplot, x, ubw, thick=3, linestyle=2, color=128
oplot, x, un, thick=3, linestyle=0

xyouts, 0.84, 1.3, 'Fromm', align=1
oplot, [0.85,0.95], [1.3,1.3], linestyle=0, thick=3

xyouts, 0.84, 1.2, 'Lax-Wendroff', align=1
oplot, [0.85,0.95], [1.2,1.2], linestyle=0, thick=3, color=128

xyouts, 0.84, 1.1, 'Beam-Warming', align=1
oplot, [0.85,0.95], [1.1,1.1], linestyle=2, thick=3, color=128




end
