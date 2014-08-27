pro ploterror

!p.charsize = 1.25

data = readtab('derivatives.out', /double)

h = data[0,*]
dp = data[1,*]
dm = data[2,*]
d0 = data[3,*]
d3 = data[4,*]

loadct, 12
plot, h, abs(dp), /xlog, /ylog, thick=3, xtitle="h", ytitle="error", yrange=[1.e-8,0.1],xrange=[0.005,0.2], xstyle=1
oplot, h, abs(dm), thick=3, color=128
oplot, h, abs(d0), thick=3, color=192
oplot, h, abs(d3), thick=3, color=16

x = findgen(10)*(1-1.e-4)/(10) + 1.e-4
print, x

oplot, x, 4.29e-2*(10.*x), linestyle=1
oplot, x, 9.e-4*(10.*x)^2, linestyle=1
oplot, x, 6.82e-5*(10.*x)^3, linestyle=1

end
