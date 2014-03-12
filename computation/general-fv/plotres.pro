pro plotres

!p.multi=[0,3,2]

loadct, 0
!p.charsize = 2.0
!p.background = 255
!p.color = 0

window, xsize=1200, ysize=800

; upwind
data = readtab('advect-MC-packet-nx=16-t=0.00', /double)
data2 = readtab('advect-MC-packet-nx=16-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-1.25,1.25], ystyle=1, title="nx = 16"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0

; Fromm
data = readtab('advect-MC-packet-nx=32-t=0.00', /double)
data2 = readtab('advect-MC-packet-nx=32-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-1.25,1.25], ystyle=1, title="nx = 32"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; L-W
data = readtab('advect-MC-packet-nx=64-t=0.00', /double)
data2 = readtab('advect-MC-packet-nx=64-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-1.25,1.25], ystyle=1, title="nx = 64"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; minmod
data = readtab('advect-MC-packet-nx=128-t=0.00', /double)
data2 = readtab('advect-MC-packet-nx=128-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-1.25,1.25], ystyle=1, title="nx = 128"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; MC
data = readtab('advect-MC-packet-nx=256-t=0.00', /double)
data2 = readtab('advect-MC-packet-nx=256-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-1.25,1.25], ystyle=1, title="nx = 256"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


; superbee
data = readtab('advect-MC-packet-nx=512-t=0.00', /double)
data2 = readtab('advect-MC-packet-nx=512-t=5.00', /double)

x = data[0,*]
u0 = data[1,*]
un = data2[1,*]

plot, x, un, thick=3, xtitle="x", ytitle="u", xrange=[0,1.0], xstyle=1, yrange=[-1.25,1.25], ystyle=1, title="nx = 512"
oplot, x, u0, thick=3, linestyle=0, color=168
oplot, x, un, thick=3, linestyle=0, color=0


end
