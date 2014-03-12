pro plothat

!p.multi=[0,3,2]

loadct, 0
!p.charsize = 2.0
!p.background = 255
!p.color = 0

window, xsize=1200, ysize=800

; upwind
data = readtab('advect-upwind-tophat-t=0.00', /double)
data2 = readtab('advect-upwind-tophat-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="unsplit"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0

; Fromm
data = readtab('advect-centered-tophat-t=0.00', /double)
data2 = readtab('advect-centered-tophat-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="Fromm"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; L-W
data = readtab('advect-LW-tophat-t=0.00', /double)
data2 = readtab('advect-LW-tophat-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="Lax-Wendroff"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; minmod
data = readtab('advect-minmod-tophat-t=0.00', /double)
data2 = readtab('advect-minmod-tophat-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="minmod limiter"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; MC
data = readtab('advect-MC-tophat-t=0.00', /double)
data2 = readtab('advect-MC-tophat-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="MC limiter"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; superbee
data = readtab('advect-superbee-tophat-t=0.00', /double)
data2 = readtab('advect-superbee-tophat-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="superbee limiter"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


end
