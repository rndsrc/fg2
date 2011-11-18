;; Copyright (C) 2011 Chi-kwan Chan
;; Copyright (C) 2011 NORDITA

;; This file is part of fg2.

;; Fg2 is free software: you can redistribute it and/or modify it
;; under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; Fg2 is distributed in the hope that it will be useful, but WITHOUT
;; ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
;; or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
;; License for more details.

;; You should have received a copy of the GNU General Public License
;; along with fg2.  If not, see <http://www.gnu.org/licenses/>.

pro vis, i, all=all, png=png

  if not keyword_set(all) then all = 0
  if not keyword_set(png) then png = 0

  openr, lun, string(i, format='(i04)') + '.raw', /get_lun
  size = lonarr(4)
  readu, lun, size
  if size[3] eq 8 then data = dblarr(size[2], size[1], size[0]) $
  else                 data = fltarr(size[2], size[1], size[0])
  readu, lun, data
  close, lun & free_lun, lun

  x = (dindgen(size[0]) + 0.5) / size[0] - 0.5
  y = (dindgen(size[1]) + 0.5) / size[1] - 0.5
  data[0,*,*] = exp(data[0,*,*])
  data[3,*,*] = exp(data[3,*,*])

  if png then begin
    dsaved = !d.name
    set_plot, 'z'
    device, set_resolution=[size[0],size[1]], decompose=0
  endif else begin
    window, xSize=size[0], ySize=size[1]
  endelse

  if all then begin
    psaved=!p.multi
    !p.multi=[0,2,2]
    shade_surf, transpose(reform(data[0,*,*])), x, y, charsize=2
    shade_surf, transpose(reform(data[1,*,*])), x, y, charsize=2
    shade_surf, transpose(reform(data[2,*,*])), x, y, charsize=2
    shade_surf, transpose(reform(data[3,*,*])), x, y, charsize=2
    !p.multi=psaved
  endif else begin
    tvscl, transpose(reform(data[0,*,*]))
  endelse

  if png then begin
    write_png, string(i, format='(i04)') + '.png', tvrd()
    set_plot, dsaved
  endif

end
