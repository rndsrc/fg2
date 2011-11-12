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

pro vis, i

  openr, lun, string(i, format='(i04)') + '.raw', /get_lun
  size = lonarr(4)
  readu, lun, size
  if size[3] eq 8 then data = dblarr(size[2], size[1], size[0]) $
  else                 data = fltarr(size[2], size[1], size[0])
  readu, lun, data
  close, lun & free_lun, lun

  h = size[0] / 2
  x = (findgen(size[0]) + 0.5) / size[0]
  y = 0.5 * (data[0,*,h-1] + data[0,*,h])

  !p.multi=[0,2,2]
  shade_surf, data[0,*,*], charsize=2
  plot, x, y
  shade_surf, data[1,*,*], charsize=2
  shade_surf, data[2,*,*], charsize=2

end