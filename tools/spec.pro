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

pro spec, i, vorticity=vorticity, eps=eps, png=png

  name = string(i, format='(i04)')
  print, 'Loading "' + name + '.raw"'
  data = load(name + '.raw')

  u1 = fft(data.d[*,*,1])
  u2 = fft(data.d[*,*,2])
  E  = (abs(u1)^2 + abs(u2)^2) / 2

  j  = getk(E)
  k1 = j.k1
  k2 = j.k2
  k  = sqrt(k1^2 + k2^2)
  s  = oned(k, E)

  dsaved = !d.name
  if keyword_set(eps) then begin
    set_plot, 'ps'
    device, filename=name + '.eps', /encap, /color, /inch, xSize=4, ySize=4
    csz = 1
  endif else if keyword_set(png) then begin
    set_plot, 'z'
    device, decompose=0, set_resolution=[640,640], set_pixel_depth=24
    csz = 2
  endif else if !d.window eq -1 then begin
    window, xSize=640, ySize=640, retain=2
    device, decompose=0
    csz = 2
  endif
  loadct, 39

  plot, [1,max(s.b)], [1e-15,1], /nodata, /xStyle, /yStyle, /xLog, /yLog, $
        xTitle='Wavenumber k', title='Frame = ' + name, $
        yTitle='Shell-integrated energy spectrum E(k)', charSize=csz

  oplot, [1,s.c], [1,s.c]^(-5./3), lineStyle=1
  oplot, [1,s.c], [1,s.c]^(-2   ), lineStyle=2
  oplot, [1,s.c], [1,s.c]^(-3   ), lineStyle=3

  oplot, s.k, s.E, thick=2

  if keyword_set(vorticity) then begin
    j = getk(u1, /zeronyquist)
    w = complex(0, j.k1) * u2 - complex(0, j.k2) * u1
    Z = abs(w)^2 / 2
    s = oned(k, Z)
    oplot, s.k, s.E / s.k^2, thick=2, color=254
  endif

  if keyword_set(eps) then begin
    device, /close
  endif else if keyword_set(png) then begin
    write_png, name + '.png', tvrd(/true)
  endif
  set_plot, dsaved

end
