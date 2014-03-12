pro plotburger

!p.multi=[0,3,2]

loadct, 0
!p.charsize = 2.0
!p.background = 255
!p.color = 0

window, xsize=1200, ysize=800

; upwind
data = readtab('burger-upwind-tophat-nx=128-t=0.00', /double)
data2 = readtab('burger-upwind-tophat-nx=128-t=0.30', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="unsplit"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0

; Fromm
data = readtab('burger-centered-tophat-nx=128-t=0.00', /double)
data2 = readtab('burger-centered-tophat-nx=128-t=0.30', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="Fromm"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; L-W
data = readtab('burger-LW-tophat-nx=128-t=0.00', /double)
data2 = readtab('burger-LW-tophat-nx=128-t=0.30', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="Lax-Wendroff"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; minmod
data = readtab('burger-minmod-tophat-nx=128-t=0.00', /double)
data2 = readtab('burger-minmod-tophat-nx=128-t=0.30', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="minmod limiter"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; MC
data = readtab('burger-MC-tophat-nx=128-t=0.00', /double)
data2 = readtab('burger-MC-tophat-nx=128-t=0.30', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="MC limiter"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; superbee
data = readtab('burger-superbee-tophat-nx=128-t=0.00', /double)
data2 = readtab('burger-superbee-tophat-nx=128-t=0.30', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-0.25,1.25], ystyle=1, title="superbee limiter"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


end
