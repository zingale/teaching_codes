


function phi_l, p, rho_l, u_l, p_l

gamma = 1.4

c_l = sqrt(gamma*p_l/rho_l)

if (p LE p_l) then begin
    ; rarefaction
    u = u_l + (2.0*c_l/(gamma-1.0))*(1.0 - (p/p_l)^((gamma-1.0)/(2.0*gamma)))
endif else begin
    ; shock
    beta = (gamma+1)/(gamma-1)
    u = u_l + (2.0*c_l/sqrt(2.0*gamma*(gamma-1.0)))*(1.0 - p/p_l)/sqrt(1.0 + beta*p/p_l)
endelse

return, u

end



function phi_r, p, rho_r, u_r, p_r

gamma = 1.4

c_r = sqrt(gamma*p_r/rho_r)

if (p LE p_r) then begin
    ; rarefaction
    u = u_r - (2.0*c_r/(gamma-1.0))*(1.0 - (p/p_r)^((gamma-1.0)/(2.0*gamma)))
endif else begin
    beta = (gamma+1.0)/(gamma-1.0)
    u = u_r - (2.0*c_r/sqrt(2.0*gamma*(gamma-1.0)))*(1.0- p/p_r)/sqrt(1.0 + beta*p/p_r)
endelse

return, u

end


loadct, 12
!p.charsize = 1.25
!p.background = 255
!p.color = 0
!p.thick=2
!x.thick=2
!y.thick=2

window, xsize=600, ysize=500

; plot the Hugoniot loci -- see LeVeque 

;rho_l = 1.0
;u_l = 0.0
;p_l = 1.0

;rho_r = 0.125
;u_r = 0.0
;p_r = 0.1

;pmax = 1.5
;umin = -4.0
;umax = 4.0
;pskip = 0.05
;uskip = 0.1

; -------------------------------------------------
; double rarefaction (Toro test 2)
rho_l = 1.0
u_l = -2.0
p_l = 0.4

rho_r = 1.0
u_r = 2.0
p_r = 0.4

; plotting parameters
pmax = 0.5
umin = -2.5
umax = 2.5
pskip = 0.01
uskip=0.1
;---------------------------------------------------


;rho_l = 1.0
;u_l = 0.0
;p_l = 1000.0

;rho_r = 1.0
;u_r = 2.0
;p_r = 0.01

; plotting parameters
;pmax = 1200.0
;umin = -20
;umax = 80
;pskip = 20
;uskip=0.1

npts = 1000

p = pmax*findgen(npts)/npts

s1 = dblarr(npts)
s2 = dblarr(npts)

for i = 0, npts-1 do begin
    s1[i] = phi_l(p[i], rho_l, u_l, p_l)
    s2[i] = phi_r(p[i], rho_r, u_r, p_r)
endfor

; find out which states in the left are due to a rarefaction (p <= p_l)
index_l = (where(p GT p_l))[0]
index_r = (where(p GT p_r))[0]


plot, [0], [umin], xrange=[0,pmax], yrange=[umin,umax], xstyle=1, ystyle=1, xtitle='p', ytitle='u'

; plot the left state
oplot, p[0:index_l], s1[0:index_l], linestyle=2, color=96, thick=3
oplot, p[index_l:*], s1[index_l:*], linestyle=0, color=96, thick=3

; plot the right state
oplot, p[0:index_r], s2[0:index_r], linestyle=2, color=192, thick=3
oplot, p[index_r:*], s2[index_r:*], linestyle=0, color=192, thick=3

; plot the labels
oplot, [p_l], [u_l], psym=4
xyouts, p_l+pskip, u_l+uskip, 'left', charsize=1.25
oplot, [p_r], [u_r], psym=4
xyouts, p_r+pskip, u_r-uskip, 'right', charsize=1.25


end








