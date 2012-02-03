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

function sp2c, in, size=size, zoom=zoom

  rmin = 3.0d ;; hardwired, see C code disk.cuh

  n  = size(in, /dim)
  n1 = n[0]
  n2 = n[1]

  hd2 = 0.5d * !pi / n2
  hd1 = alog(hd2 + sqrt(hd2 * hd2 + 1.0d))
  dx  = 2.0d * hd1

  if not keyword_set(size) then size = n1
  if not keyword_set(zoom) then zoom = (size / 2) / (rmin * exp(n1 * dx))

  out = dblarr(size / 2, size) + 1.1 * min(in)
  x   = dindgen(size / 2) + 0.5d

  for j = 0, size-1 do begin
     z = size / 2 - 0.5d - j
     r = sqrt(x * x + z * z)
     t = atan(x, z)

     r = fix(alog(r / rmin / zoom) / dx)
     t = fix(n2 * t / !pi)

     i = where(0 le r and r lt n1, c)
     if c ne 0 then out[i, j] = in[r[i], t[i]]
  endfor

  return, out

end
