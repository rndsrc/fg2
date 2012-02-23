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

function getk, f, zeronyquist=zeronyquist

  n  = size(f, /dimensions)
  n1 = n[0]
  n2 = n[1]

  k1 = [findgen(n1-n1/2), -reverse(findgen(n1/2)+1)]
  k2 = [findgen(n2-n2/2), -reverse(findgen(n2/2)+1)]
  if keyword_set(zeroNyquist) then begin
    if n1 mod 2 eq 0 then k1[n1/2] = 0
    if n2 mod 2 eq 0 then k2[n2/2] = 0
  endif

  return, {k1:rebin(k1,n1,n2), k2:rebin(transpose(k2),n1,n2)}

end
