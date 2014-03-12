; quick hack to get color gif off the screen 
; get the plot on screen, and then type 
; IDL> color_gif, 'filename.gif'

pro color_bitmap, filename

; save the various plot options -- so we don't mess anything up . . .
;old_plot = !p

;set_plot, 'Z'
;tvlct, red, green, blue, /get

;set_plot, 'X'
;!p = old_plot

a = tvrd(TRUE=1)

image = color_quan(a, 1, red, green, blue)

if (!VERSION.RELEASE LE 5.3) then begin
    write_gif, filename, image, red, green, blue
endif else begin
    write_png, filename, image, red, green, blue
endelse

end


