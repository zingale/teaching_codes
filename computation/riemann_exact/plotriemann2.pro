pro plotriemann2

loadct, 0
!p.charsize = 2.5
!p.background = 255
!p.color = 0

!p.multi=[0,1,3]
!x.margin=[12,3]
!y.margin=[3,2]

window, xsize=600, ysize=675

data = readtab('test2-exact.out', /double)

x   = data[0,*]
rho = data[1,*]
u   = data[2,*]
p   = data[3,*]
e   = data[4,*]

; create the characteristic variables
gamma = 1.4
c = sqrt(gamma*p/rho)

w1 = u - 2.0*c/(gamma-1.0)
w2 = 8.26e7*alog(p/rho^gamma)  ; entropy k_B/m_p = 8.26e7
w3 = u + 2.0*c/(gamma-1.0)

print, p/rho^gamma

plot, x, w1, xtitle='x', ytitle='w1 = u - 2c/(!4c!3-1)', thick=3, ystyle=1
plot, x, w2, xtitle='x', ytitle='w2 = s', thick=3, ystyle=1
plot, x, w3, xtitle='x', ytitle='w3 = u + 2c/(!4c!3-1)', thick=3, ystyle=1


end
