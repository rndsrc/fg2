;; Copyright (C) 2011,2012 Chi-kwan Chan
;; Copyright (C) 2011,2012 NORDITA
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

function load, name

  openr, lun, name, /get_lun
  size = lonarr(6)
  readu, lun, size
  if size[3] eq 8 then begin
    time = dcomplex(0.0, 0.0) ; sizeof(long double) = 16 on my machines
    info = dblarr(3)
    data = dblarr(size[2], size[1], size[0])
  endif else begin
    time = double(0.0)
    info = fltarr(3)
    data = fltarr(size[2], size[1], size[0])
  endelse
  readu, lun, time
  readu, lun, info
  readu, lun, data
  close, lun & free_lun, lun

  x = info[0] * (dindgen(size[0]) + 0.5) / size[0]
  y = info[1] * (dindgen(size[1]) + 0.5) / size[1]
  d = transpose(data, [2,1,0])

  return, {x:x, y:y, d:d}

end
