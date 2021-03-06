/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "partpix_map.h"
using namespace std;


template<typename T> void Partpix_Map<T>::Import_degrade
  (const Partpix_Map<T> &orig, const Healpix_Map<double>& resolutionMask, bool pessimistic, double cutoff)
  {
  planck_assert(nside_<orig.nside_,"Import_degrade: this is no degrade");
  int fact = orig.nside_/nside_;
  int factr = nside_/resolutionMask.Nside();
  planck_assert (orig.nside_==nside_*fact,
    "the larger Nside must be a multiple of the smaller one");
  pix2xyf to_xyf = (scheme_==RING) ?
    &Healpix_Base::ring2xyf : &Healpix_Base::nest2xyf;
  xyf2pix from_xyf = (orig.scheme_==RING) ?
    &Healpix_Base::xyf2ring : &Healpix_Base::xyf2nest;
  pix2xyf res_to_xyf = (resolutionMask.Scheme()==RING) ?
    &Healpix_Base::ring2xyf : &Healpix_Base::nest2xyf;

  int minhits = pessimistic ? fact : 1;

  int arrayindex = 0;
  for( int mr=0; mr<resolutionMask.Npix(); ++mr ){
    if( resolutionMask[mr] < 0.5 )
        continue;
    int xr,yr,fr;
    (resolutionMask.*res_to_xyf)(mr,xr,yr,fr);
    for (int jr=factr*yr; jr<factr*(yr+1); ++jr){
      for (int ir=factr*xr; ir<factr*(xr+1); ++ir){
        int thispix = (this->*from_xyf)(ir,jr,fr);
        int x,y,f;
        (this->*to_xyf)(thispix,x,y,f);
        int hits = 0;
        double sum = 0;
        for (int j=fact*y; j<fact*(y+1); ++j){
          for (int i=fact*x; i<fact*(x+1); ++i){
            int origpix = (orig.*from_xyf)(i,j,f);
            if (!approx<double>(orig.partmap_at_highresindex(origpix),Healpix_undef))
              {
              ++hits;
              sum += orig.partmap_at_highresindex(origpix);
              }
          }
        }
        pixel_mapping_arraytohigh[arrayindex] = thispix;
        T result = T((hits<minhits) ? Healpix_undef : sum/hits);
        if( cutoff != Healpix_undef ){
            if( result > cutoff )
                result = T(1);
        }
        partmap[arrayindex] = result;
        ++arrayindex;
    }
  }
  }
  sort_partmap( pixel_mapping_arraytohigh, partmap );
  }
  
  
template<typename T> void Partpix_Map<T>::Import_degrade
  (const Partpix_Map<T> &orig, const Healpix_Map<int>& resolutionMask, bool pessimistic, double cutoff)
  {
  planck_assert(nside_<orig.nside_,"Import_degrade: this is no degrade");
  int fact = orig.nside_/nside_;
  int factr = nside_/resolutionMask.Nside();
  planck_assert (orig.nside_==nside_*fact,
    "the larger Nside must be a multiple of the smaller one");
  pix2xyf to_xyf = (scheme_==RING) ?
    &Healpix_Base::ring2xyf : &Healpix_Base::nest2xyf;
  xyf2pix from_xyf = (orig.scheme_==RING) ?
    &Healpix_Base::xyf2ring : &Healpix_Base::xyf2nest;
  pix2xyf res_to_xyf = (resolutionMask.Scheme()==RING) ?
    &Healpix_Base::ring2xyf : &Healpix_Base::nest2xyf;

  int minhits = pessimistic ? fact : 1;

  int arrayindex = 0;
  for( int mr=0; mr<resolutionMask.Npix(); ++mr ){
    if( resolutionMask[mr] == 0 )
        continue;
    int xr,yr,fr;
    (resolutionMask.*res_to_xyf)(mr,xr,yr,fr);
    for (int jr=factr*yr; jr<factr*(yr+1); ++jr){
      for (int ir=factr*xr; ir<factr*(xr+1); ++ir){
        int thispix = (this->*from_xyf)(ir,jr,fr);
        int x,y,f;
        (this->*to_xyf)(thispix,x,y,f);
        int hits = 0;
        double sum = 0;
        for (int j=fact*y; j<fact*(y+1); ++j){
          for (int i=fact*x; i<fact*(x+1); ++i){
            int origpix = (orig.*from_xyf)(i,j,f);
            if (!approx<double>(orig.partmap_at_highresindex(origpix),Healpix_undef))
              {
              ++hits;
              sum += orig.partmap_at_highresindex(origpix);
              }
          }
        }
        pixel_mapping_arraytohigh[arrayindex] = thispix;
        T result = T((hits<minhits) ? Healpix_undef : sum/hits);
        if( cutoff != Healpix_undef ){
            if( result > cutoff )
                result = T(1);
        }
        partmap[arrayindex] = result;
        ++arrayindex;
    }
  }
  }
  sort_partmap( pixel_mapping_arraytohigh, partmap );
  }

template void Partpix_Map<double>::Import_degrade
  (const Partpix_Map<double> &orig, const Healpix_Map<double>& resolutionMask, bool pessimistic, double cutoff);
template void Partpix_Map<double>::Import_degrade
  (const Partpix_Map<double> &orig, const Healpix_Map<int>& resolutionMask, bool pessimistic, double cutoff);
template void Partpix_Map<int>::Import_degrade
  (const Partpix_Map<int> &orig, const Healpix_Map<int>& resolutionMask, bool pessimistic, double cutoff);
template void Partpix_Map<float>::Import_degrade
  (const Partpix_Map<float> &orig, const Healpix_Map<double>& resolutionMask, bool pessimistic, double cutoff);
template void Partpix_Map<float>::Import_degrade
  (const Partpix_Map<float> &orig, const Healpix_Map<int>& resolutionMask, bool pessimistic, double cutoff);
  
template<typename T> void Partpix_Map<T>::minmax (T &Min, T &Max) const
  {
  Min = T(1e30); Max = T(-1e30);
  for (int m=0; m<npartpix_; ++m)
    {
    T val = partmap[m];
    if (!approx<double>(val,Healpix_undef))
      {
      if (val>Max) Max=val;
      if (val<Min) Min=val;
      }
    }
  }

template void Partpix_Map<float>::minmax (float &Min, float &Max) const;
template void Partpix_Map<double>::minmax (double &Min, double &Max) const;
