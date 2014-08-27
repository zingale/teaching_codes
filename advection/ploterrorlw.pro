pro ploterrorlw

loadct, 0
!p.charsize = 1.25
!p.background = 255
!p.color = 0

data = readtab('error-square-lw.out', /double)

x = data[1,*]
e1 = data[2,*]
e2 = data[3,*]
einfty = data[4,*]

data = readtab('error-square-upwind.out', /double)

x = data[1,*]
ue1 = data[2,*]
ue2 = data[3,*]
ueinfty = data[4,*]

plot, x, e1, /xlog, /ylog, xtitle="!4D!3x", ytitle="error", thick=3, xrange=[9.e-4,0.3], xstyle=1

oplot, x, e2, thick=3, linestyle=1
oplot, x, einfty, thick=3, linestyle=2

oplot, x, ue1, thick=3, linestyle=0, color=128
oplot, x, ue2, thick=3, linestyle=1, color=128
oplot, x, ueinfty, thick=3, linestyle=2, color=128

xyouts, 0.09, 4.5e-2, '1-norm', align=1
oplot, [0.1,0.2], [4.5e-2,4.5e-2], linestyle=0, thick=3
xyouts, 0.09, 3.e-2, '2-norm', align=1
oplot, [0.1,0.2], [3.e-2,3.e-2], linestyle=1, thick=3
xyouts, 0.09, 2.e-2, 'infinity-norm', align=1
oplot, [0.1,0.2], [2.e-2,2.e-2], linestyle=2, thick=3

end
