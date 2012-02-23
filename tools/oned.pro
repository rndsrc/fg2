;; Copyright (C) 2012 Chi-kwan Chan
;; Copyright (C) 2012 NORDITA
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

function oned, k, E

  kG = min(size(k, /dimensions)) / 2 + 0.99 ; Galerkin cutoff wavenumber
  n  = ceil(kG)
  kb = findgen(n + 1) + 0.5
  kc = findgen(n) + 1
  nb = lonarr(n)
  Ec = fltarr(n)

  for i = 0, n-1 do begin
    case i of
      0    : j = where(0.0   lt k and k lt kb[i+1], count)
      n-1  : j = where(kb[i] le k and k le kb[i+1], count)
      else : j = where(kb[i] le k and k lt kb[i+1], count)
    endcase
    if count ne 0 then begin
      Ec[i] = total(E[j])
      nb[i] = count
    endif
  endfor

  return, {E:Ec, k:kc, n:nb, b:kb, c:kG}

end
