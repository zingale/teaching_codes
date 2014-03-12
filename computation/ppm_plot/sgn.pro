

function sgn, x

if x LT 0.0 then begin
    return, -1.0
endif else if x EQ 0.0 then begin
    return, 0.0
endif else begin
    return, 1.0
endelse

end


