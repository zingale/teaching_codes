
; make some descriptive PPM plots.  Use a small number of ODD points,
; with the center value = 0.0, and the zone width = 1, for simplicity.
; the plots will be labeled as i-1/2, etc., but the -1/2 in the
; subscript is always the real value in this code.

!p.charthick=1
!p.thick=2

; start by defining the grid, number of interior points + ng
npts = 5  ; for these plots, we require npts to be odd
ng = 4    ; ppm always has 4 gc

imin = ng
imax = ng+npts-1

; for simplicity in the graphs, dx = 1
dx = 1

; x[npts/2] will be 0.0
xl = (findgen(npts+2*ng) - ng - npts/2 - 0.5)*dx + xmin
xr = (findgen(npts+2*ng) - ng - npts/2 + 0.5)*dx + xmin
x = 0.5*(xl + xr)

print, x[imin:imax]

; define the initial data
a = fltarr(npts+2*ng)  
;a = 8*x*x + 3.0
a[imin+npts/2] = 0.6
a[imin:imin+npts/2-2] = 1.0
a[imin+npts/2-1] = 0.95

a[imin+npts/2+1] = 0.1
a[imin+npts/2+2:*] = 0.05

; fill the GC
a[0:imin-1] = a[imin]
a[imax+1:*] = a[imax]


xticknames = strarr(npts+1)
xtickvalue = fltarr(npts+1)

for i = 0, npts do begin
    num = i*2 - npts
    xticknames[i] = '!4n!3!Ii' + strcompress(string(num,format="(I+)"), /remove_all) + '/2!N'
    xtickvalue[i] = float(num)/2
endfor

; plot the zone averaged initial data
plot, x[imin:imax], a[imin:imax], psym=3, xtickn=xticknames, xtickv=xtickvalue, xticks=npts, $
  charsize=2, yrange=[0.0, 1.5*max(a[imin:imax])], xrange=[xtickvalue[0],xtickvalue[npts]], $
  xticklen=0.001, xstyle=9, ystyle=9, yticks=1, ytickformat='nolabel', $
  xtitle = '!4n!3', ytitle = 'a(!4n!3)', yticklen=0.001, xthick=2, ythick=2


; plot the zone averages
for i = imin, imax do begin
    oplot, [xl[i],xr[i]], [a[i],a[i]], linestyle=2
endfor

; plot the vertical lines
for i = imin, imax do begin
    oplot, [xr[i], xr[i]], [0, 1.5*max(a[imin:imax])], linestyle=1
endfor

;-----------------------------------------------------------------------------
; now get the interface values
;-----------------------------------------------------------------------------

al = fltarr(npts+2*ng)
ar = fltarr(npts+2*ng)

; al[i] is the i-1/2 interface
for i = imin, imax+1 do begin

    ; compute the interface values using 1.6 - 1.8, 
    ; assuming constant grid spacing
    da_i = 0.5*(a[i+1] - a[i-1])
    if ((a[i+1] - a[i])*(a[i] - a[i-1]) GT 0) then begin
        da_i = sgn(da_i)*(abs(da_i) < 2.0*(abs(a[i] - a[i-1])) < 2.0*(abs(a[i] - a[i+1])))
    endif else begin
        da_i = 0.0
    endelse

    da_im1 = 0.5*(a[i] - a[i-2])
    if ((a[i] - a[i-1])*(a[i-1] - a[i-2]) GT 0) then begin
        da_im1 = sgn(da_im1)*(abs(da_im1) < 2.0*(abs(a[i-1] - a[i-2])) < 2.0*(abs(a[i-1] - a[i])))
    endif else begin
        da_im1 = 0.0
    endelse

    al[i] = a[i-1] + 0.5*(a[i] - a[i-1]) + (1./6.)*(-da_i + da_im1)
    ar[i] = al[i]

;    al[i] = (7./12.)*(a[i-1] + a[i]) - (1./12.)*(a[i+1] + a[i-2])  ; CW 1.9
;    ar[i] = al[i]
    oplot, [xl[i]], [al[i]], psym=4
endfor


;-----------------------------------------------------------------------------






end





