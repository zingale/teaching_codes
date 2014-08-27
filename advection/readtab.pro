; mz -- modified readtab.pro from  Helge Rebhan (Helge.Rebhan@gmx.net),
;       now assume the maximum array size instead of doing an
;       expensive copy -- this makes it much, much faster (and hey, 
;       memory is cheap).


;liese naechste Datenzeile, Kommentarzeilen mit ';' oder '!' skippen
FUNCTION readnextline,lun,line
  ok=0
  while not eof(lun) and ok eq 0 do begin
      readf,lun,line
      bl= byte(strtrim(line,2))
      ok=1                      ;suche Kommentarzeilen
      if bl(0) eq 59 then ok=0  ; ;
      if bl(0) eq 33 then ok=0  ; !
      IF bl(0) eq 35 then ok=0  ; #
      if strlen(line) eq 0 then ok=0
  end
  return,ok
end
;
FUNCTION readtab, filename, $
                  cols=cols, $
                  MONOTONIC=monotonic, $
                  SKIPLINES=skiplines, $
                  DOUBLE=double
;+
; NAME:
;
;       READTAB
;
; PURPOSE:
;
;       Read any kind of ASCII-tables into floating point array
;       Empty Lines or lines starting with ';','#' or '!' are ignored
;
; CATEGORY:
;
;       Input/Output
;
; CALLING SEQUENCE:
;
;       var= READTAB(filename, [cols])
;
; INPUTS:
;
;       Filename  -  Input Datafile
;
;
; KEYWORD PARAMETERS:
;
;       COLS - Vektor with colum numbers to pick, otherwise all
;              colums are returned
;
;       MONOTONIC -- if defined, assume that the first coordinate is
;                    the independent variable, and enforce that x[i+1]
;                    > x[i] by dropping any values that fail this
;                    test.  This allows you to plot a flash.dat file
;                    with restarts without getting any hiccups. (mz)
;
;       DOUBLE - use double precision
;
; OUTPUTS:
;
;       Floating point array with dimension specified by Colums and
;       Rows in the input file
;
;-

;  on_error, 2
  if n_params() eq 0 then message,'*** Aufruf var=READTAB(filename)'
  openr,lun,filename,/get_lun
  line=''

  if n_elements(monotonic) EQ 0 then monotonic = 0
  if n_elements(skiplines) EQ 0 then skiplines = 1
  if n_elements(double) EQ 0 then double = 0

;erste Zeile lesen
  if readnextline(lun,line) eq 0 then message,'*** Keine Daten !'
  bl=byte(line)
  s= size(bl) & n= s(1) -1
; Anzahl der Spalten ermitteln
  ncol=0 & ws=0 & wsold=1
  for i=0,n do begin
      if bl(i) eq 32 or bl(i) eq 09 then ws=1 else ws=0
      if ws ne wsold then begin
          if ws eq 0 then ncol=ncol+1
          wsold=ws
      end
  end
  print,'File hat ',ncol,' Spalten'
  if (NOT double) then begin
      datline= fltarr(ncol)
  endif else begin
      datline= dblarr(ncol)
  endelse

; Zeilen lesen
  point_lun,lun,0
  i=0l
  iread = 0l

; mz -- assume the maximum number of rows
  max_rows = 200000
  
  if (NOT double) then begin
      data = fltarr(ncol,max_rows)
  endif else begin
      data = dblarr(ncol,max_rows)
  endelse

; for efficiency, don't but the monotonic test in the loop
  if (NOT monotonic) then begin
      while not eof(lun) do begin
          if readnextline(lun,line) eq 1 then begin
              reads,line,datline

              if (iread mod skiplines EQ 0) then begin
                  data[*,i] = datline
                  i = i + 1
              endif

              iread = iread + 1
          end
      end
  endif else begin
      while not eof(lun) do begin
          if readnextline(lun,line) eq 1 then begin
              reads,line,datline

; make sure we aren't going backwards in the dependent variable
              if (datline[0] GE data[0,i-1 > 0] AND iread mod skiplines EQ 0) then begin
                  data[*,i] = datline
                  i = i + 1
              endif

              iread = iread + 1
          end
      end
  endelse

  data = temporary(data[*,0:i-1])

  print,'      ',i,' Zeilen gelesen'
  close,lun & free_lun,lun

  if keyword_set(cols) then begin
      cols=cols-1               ;Spaltennr ab 1 !
      data=data([cols],*)
  end
  return,data
end








