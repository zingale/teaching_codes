pro plotriemann

ipost = 1

if (ipost) then begin
    set_plot, 'PS'
    outfile = 'test4-64.ps'
    xsize = 10.0
    ysize = 7.5

    device, FILE = outfile, XSIZE = xsize, YSIZE = ysize, $
      XOFF = 0.5, YOFF = 10.5, /INCH, /COLOR, $
      BITS_PER_PIXEL = 8, LANDSCAPE=1, language_level=2
endif


loadct, 12
!p.charsize = 1.0

!p.multi=[0,2,2]


data = readtab('test4-exact.out', /double)
datag = readtab('godunov.test4.64', /double)
title = "test4, nx = 64"

x   = data[0,*]
rho = data[1,*]
u   = data[2,*]
p   = data[3,*]
e   = data[4,*]

x_god   = datag[0,*]
rho_god = datag[1,*]
u_god   = datag[2,*]
p_god   = datag[3,*]
e_god   = datag[4,*]

print, max(rho)
print, max(rho_god)

print, max(rho) > max(rho_god)

plot, x_god, rho_god, xtitle='x', ytitle='!4q!3', thick=3, ystyle=1, yrange=[min(rho < rho_god), max(rho) > max(rho_god)], xstyle=1
oplot, x, rho, thick=2, linestyle=0, color=192
oplot, x_god, rho_god, thick=3, linestyle=0

plot, x_god, u_god, xtitle='x', ytitle='u', thick=3, ystyle=1,yrange=[min(u < u_god), max(u) > max(u_god)], xstyle=1
oplot, x, u, thick=2, linestyle=0, color=192
oplot, x_god, u_god, thick=3, linestyle=0

plot, x_god, p_god, xtitle='x', ytitle='p', thick=3, ystyle=1, yrange=[min(p < p_god), max(p) > max(p_god)], xstyle=1
oplot, x, p, thick=2, linestyle=0, color=192
oplot, x_god, p_god, thick=3, linestyle=0

plot, x_god, e_god, xtitle='x', ytitle='e', thick=3, ystyle=1, yrange=[min(e < e_god), max(e) > max(e_god)], xstyle=1
oplot, x, e, thick=2, linestyle=0, color=192
oplot, x_god, e_god, thick=3, linestyle=0

xyouts, 0.5, .97, title, alignment=0.5, /normal



if (ipost) then begin
    device, /close
    set_plot, 'X'
endif


end
