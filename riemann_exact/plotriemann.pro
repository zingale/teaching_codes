pro plotriemann

loadct, 0
!p.charsize = 1.25
!p.background = 255
!p.color = 0

!p.multi=[0,2,2]
!x.margin=[7,3]
!y.margin=[3,2]

data = readtab('test3-exact.out', /double)

x   = data[0,*]
rho = data[1,*]
u   = data[2,*]
p   = data[3,*]
e   = data[4,*]


plot, x, rho, xtitle='x', ytitle='!4q!3', thick=3, ystyle=1
plot, x, u, xtitle='x', ytitle='u', thick=3, ystyle=1
plot, x, p, xtitle='x', ytitle='p', thick=3, ystyle=1
plot, x, e, xtitle='x', ytitle='e', thick=3, ystyle=1

end
