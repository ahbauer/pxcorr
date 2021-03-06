; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2010  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.jpl.nasa.gov
;
; -----------------------------------------------------------------------------
function alm_i2t, index, alm_vec, complex=complex, help=help ;, lmax=lmax, mmax=mmax
;+
; NAME:
;       alm_i2t
;
;
; PURPOSE:
;       turn an indexed list of alm (as generated by fits2alm) into a(l,m) array
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;      alm_tabular = alm_i2t(index, alm_vector, complex=)
;
;
; INPUTS:
;     index: Integer vector of size nl containing the index the 
;            of alm coefficients, related to {l,m} by the relation
;             i = l^2 + l + m + 1
;
;     alm_vector: array of alm coefficients, with dimension (nl, nalm [,nsig])
;            -- corresponding to
;                  nl   = number of {l,m} indices
;                  nalm = 2 for real and imaginary parts of alm coefficients
;                         4 for above plus corresponding error values
;                  nsig = number of signals (usually 1 for any of T E B
;                         or 3 for T,E,B together)
;
; KEYWORD PARAMETERS:
;    /complex: if set, the output array is complex with dimensions
;          (lmax+1, mmax+1, [nalm/2 , nsig]),
;      otherwise, the array is real with dimensions
;          (lmax+1, mmax+1, nalm [, nsig])
;
;    /help: if set, prints out this documentation
;
; OUTPUTS:
;     alm_tabular: real or complex array, containing all the alm for l
;         in [0,lmax] and m in [0,mmax].
;    It has dimension (lmax+1, mmax+1, nalm [, nsig]), or
;       (lmax+1, mmax+1, [nalm/2 , nsig]) if COMPLEX is set.
;
;     The single/double precision alm_tabular matches that of alm_vector
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;     2007-10-04: created
;-

routine = 'alm_i2t'
if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

syntax = 'alm_tabular = '+routine+'(index, alm_vector, COMPLEX=, HELP=)'
if n_params() ne 2 then begin
    print,syntax
    return,-1
endif

index2lm, index, l, m

; find out sizes
lmax = (defined(lmax))? lmax*1L : max(l)
mmax = (defined(mmax))? mmax*1L : max(m)
sz = size(alm_vec)
ndim = sz[0]
n1 = sz[1]
n2 = (ndim gt 1) ? sz[2] : 1L
n3 = (ndim gt 2) ? sz[3] : 1L
ni = n_elements(index)
mp1 = mmax + 1L


if (keyword_set(complex)) then begin
    ;---------------------------------------
    ; complex array
    ;---------------------------------------
    n2 = n2/2
    double = (size(/type, alm_vec) eq size(/type, 1.d0))
    sample = complex(1.0, 0.0, double=double)
    ; create array with convenient shape for IDL l,m indexing
    alm_tab = make_array(lmax+1, mp1*n2*n3, type=size(/type,sample))
    for j=0,n3-1 do begin
        for i=0,n2-1 do begin
            kshift =  mp1 * (i + n2 * j)
            alm_tab[l, m + kshift] = complex(alm_vec[*,2*i,j], alm_vec[*,2*i+1,j], double=double)
        endfor
    endfor
    alm_tab = reform(alm_tab, lmax+1, mp1, n2, n3, /Overwrite)


endif else begin
    ;---------------------------------------
    ; real array
    ;---------------------------------------
    ; create array with convenient shape for IDL l,m indexing
    alm_tab = make_array(lmax+1, mp1*n2*n3, type=size(/type,alm_vec))
    for j=0,n3-1 do begin
        for i=0,n2-1 do begin
            kshift =  mp1 * (i + n2 * j)
            alm_tab[l, m + kshift] = alm_vec[*,i,j]
        endfor
    endfor
endelse


; redimension array in place
alm_tab = reform(alm_tab, lmax+1, mp1, n2, n3, /Overwrite)
alm_tab = reform(alm_tab, /Overwrite) ; remove trailling dimensions if they are 1
   
return, alm_tab
end
