;; Copyright (C) 2011 Chi-kwan Chan
;; Copyright (C) 2011 NORDITA
;;
;; This file is part of fg2.
;;
;; Fg2 is free software: you can redistribute it and/or modify it
;; under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.
;;
;; Fg2 is distributed in the hope that it will be useful, but WITHOUT
;; ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
;; or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
;; License for more details.
;;
;; You should have received a copy of the GNU General Public License
;; along with fg2.  If not, see <http://www.gnu.org/licenses/>.

pro vis, i, rtheta=rtheta, all=all, png=png

  if not keyword_set(rtheta) then rtheta = 0
  if not keyword_set(all)    then all    = 0
  if not keyword_set(png)    then png    = 0

  if rtheta and all then begin
    print, '/rtheta and /all are incompatible; choose only one of them'
    return
  endif

  name = string(i, format='(i04)') + '.raw'
  data = load(name, rtheta=rtheta)
  print, 'Loaded "' + name + '"'

  if png then begin
    dsaved = !d.name
    set_plot, 'z'
    device, decompose=0, set_resolution=[1024, 1024], set_pixel_depth=24
  endif else if !d.window eq -1 then begin
    window, xSize=1024, ySize=1024, retain=2
    device, decompose=0
  endif

  if all then begin
    n = size(data.d) & n = n[3]
    psaved=!p.multi
    if n gt 1 then !p.multi=[0,2,(n-1)/2+1]
    x = data.x
    y = data.y
    for i = 0, n-1 do shade_surf, data.d[*,*,i], x, y, charsize=2
    !p.multi=psaved
  endif else begin
    dx = data.x[1] - data.x[0]
    if rtheta then begin
      d = sp2c(data.d[*,*,0], size=1024)
      tvscl, reverse(d)
      tvscl, d, 1
    endif else begin
      tvscl, congrid(data.d[*,*,0], 1024, 1024)
    endelse
  endelse

  if png then begin
    write_png, string(i, format='(i04)') + '.png', tvrd(/true)
    set_plot, dsaved
  endif

end
