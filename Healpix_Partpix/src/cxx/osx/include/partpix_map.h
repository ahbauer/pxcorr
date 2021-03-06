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

/*! \file healpix_map.h
 *  Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PARTPIX_MAP_H
#define PARTPIX_MAP_H

#include <vector>
using std::vector;
#include "arr.h"
#include "healpix_base.h"
#include "healpix_map.h"


/*! A HEALPix map of a given datatype */
template<typename T> class Partpix_Map: public Healpix_Base
  {
  private:
    arr<T> partmap;
    vector<int> pixel_mapping_arraytohigh;
    int npartpix_;
    
  public:
    /*! Constructs an unallocated map. */
    Partpix_Map () {}
    /*! Constructs a map with a given \a order and the ordering
        scheme \a scheme. */
    Partpix_Map (int order, Healpix_Ordering_Scheme scheme)
       : Healpix_Base (order, scheme) {}
    /*! Constructs a map with a given \a nside and the ordering
        scheme \a scheme. */
    // Partpix_Map (int nside, Healpix_Ordering_Scheme scheme, const nside_dummy)
    //   : Healpix_Base (nside, scheme, SET_NSIDE), partmap(npix_) {}
    /*! Constructs a map from the contents of \a data and sets the ordering
        scheme to \a Scheme. The size of \a data must be a valid HEALPix
        map size. */
    // Partpix_Map (const arr<T> &data, Healpix_Ordering_Scheme scheme)
    //   : Healpix_Base (npix2nside(data.size()), scheme, SET_NSIDE), partmap(data) {}

    Partpix_Map( int highResOrder, const Healpix_Map<int>& resolutionMask )
      : Healpix_Base (highResOrder, resolutionMask.Scheme()) {

          planck_assert( highResOrder >= resolutionMask.Order(), "Footprint order is higher than the map order!" );

        // make the mapping between high res pixels and array indices

        // what high-resolution pixel numbers fall inside the mask?
        // add those to the pixel_mapping hash: pixel_mapping[highrespix] = array_index
        // where array_index is the index in arr<T> map

        // how many high-res pixels will we need?
        int n_lowres_pixels = 0;
        for( int i=0; i<resolutionMask.Npix(); ++i ){
            if( resolutionMask[i] > 0.5 )
                ++n_lowres_pixels;
        }
        npartpix_ = int(n_lowres_pixels*double(npix_)/double(resolutionMask.Npix()));

        // initialize the arrays
        partmap = arr<T>(npartpix_);
        //pixel_mapping_arraytohigh = arr<int>(npartpix_);
        pixel_mapping_arraytohigh = vector<int>(npartpix_);

        // assign the mapping
        int fact = nside_/resolutionMask.Nside();
        xyf2pix from_xyf = (scheme_==RING) ?
          &Healpix_Base::xyf2ring : &Healpix_Base::xyf2nest;
      pix2xyf res_to_xyf = (resolutionMask.Scheme()==RING) ?
          &Healpix_Base::ring2xyf : &Healpix_Base::nest2xyf;
        int64 array_index=0;
        int nrpix = resolutionMask.Npix();
        for(int m=0; m<nrpix; ++m){
            if( resolutionMask[m] < 0.5 )
                continue;
          int x,y,f;
          (resolutionMask.*res_to_xyf)(m,x,y,f);
          // if( resolutionMask.Scheme()==RING )
          //     resolutionMask.ring2xyf(m,x,y,f);
          // else
          //     resolutionMask.nest2xyf(m,x,y,f);
          int ymax = fact*(y+1);
          int xmax = fact*(x+1);
          for (int j=fact*y; j<ymax; ++j)
            for (int i=fact*x; i<xmax; ++i){
              pixel_mapping_arraytohigh[array_index] = (this->*from_xyf)(i,j,f);
              ++array_index;
            }
        } // for m
        planck_assert( array_index == npartpix_, "array index != npartpix_" );        
        //sort_partmap( pixel_mapping_arraytohigh, partmap );
        sort( pixel_mapping_arraytohigh.begin(), pixel_mapping_arraytohigh.end() );
    }
    // The same as before, since we don't want to make this a double template for resolutionMask too.
    Partpix_Map( int highResOrder, const Healpix_Map<double>& resolutionMask )
      : Healpix_Base (highResOrder, resolutionMask.Scheme()) {
          
          planck_assert( highResOrder >= resolutionMask.Order(), "Mask order is higher than the map order!" );
          
        // make the mapping between high res pixels and array indices

        // what high-resolution pixel numbers fall inside the mask?
        // add those to the pixel_mapping hash: pixel_mapping[highrespix] = array_index
        // where array_index is the index in arr<T> map

        // how many high-res pixels will we need?
        int n_lowres_pixels = 0;
        int nrpix = resolutionMask.Npix();
        for( int i=0; i<nrpix; ++i ){
            if( resolutionMask[i] > 0.5 )
                ++n_lowres_pixels;
        }
        npartpix_ = int(n_lowres_pixels*double(npix_)/double(nrpix));

        // initialize the arrays
        partmap = arr<T>(npartpix_);
        //pixel_mapping_arraytohigh = arr<int>(npartpix_);
        pixel_mapping_arraytohigh = vector<int>(npartpix_);

        // assign the mapping
        int fact = nside_/resolutionMask.Nside();
        xyf2pix from_xyf = (scheme_==RING) ?
          &Healpix_Base::xyf2ring : &Healpix_Base::xyf2nest;
        pix2xyf res_to_xyf = (resolutionMask.Scheme()==RING) ?
            &Healpix_Base::ring2xyf : &Healpix_Base::nest2xyf;
        int array_index=0;
        for(int m=0; m<nrpix; ++m){
            if( resolutionMask[m] < 0.5 )
                continue;
          int x,y,f;
          (resolutionMask.*res_to_xyf)(m,x,y,f);
          // if( resolutionMask.Scheme()==RING )
          //     resolutionMask.ring2xyf(m,x,y,f);
          // else
          //     resolutionMask.nest2xyf(m,x,y,f);
          int ymax = fact*(y+1);
          int xmax = fact*(x+1);
          for (int j=fact*y; j<ymax; ++j)
            for (int i=fact*x; i<xmax; ++i){
              pixel_mapping_arraytohigh[array_index] = (this->*from_xyf)(i,j,f);
              ++array_index;
            }
        } // for m
        planck_assert( array_index == npartpix_, "array index != npartpix_" );
        //sort_partmap( pixel_mapping_arraytohigh, partmap );
        sort( pixel_mapping_arraytohigh.begin(), pixel_mapping_arraytohigh.end() );
    }
    
    void clear(){
        partmap.dealloc();
        pixel_mapping_arraytohigh.clear();
    }
    
    // how to access the pixel values using the old indices?
    T partmap_at_highresindex( int hrindex ){
        // binary search through the pixel mapping array
        // since we care about memory more than CPU time.
        // vector<int>::iterator pos = lower_bound(pixel_mapping_arraytohigh.begin(), pixel_mapping_arraytohigh.end(), hrindex);
        // planck_assert( *pos == hrindex, "High resolution pixel index not found in the partial map" );
        // return partmap[pos-pixel_mapping_arraytohigh.begin()];
        return partmap_ref_at_highresindex( hrindex );
    }
    const T partmap_at_highresindex( int hrindex ) const{
        return partmap_ref_at_highresindex( hrindex );
    }
    T& partmap_ref_at_highresindex( int hrindex ){
        // binary search through the pixel mapping array
        // since we care about memory more than CPU time.
        vector<int>::iterator pos = lower_bound(pixel_mapping_arraytohigh.begin(), pixel_mapping_arraytohigh.end(), hrindex);
        planck_assert( *pos == hrindex, "High resolution pixel index not found in the partial map" );
        return partmap[pos-pixel_mapping_arraytohigh.begin()];
    }
    const T& partmap_ref_at_highresindex( int hrindex ) const{
        vector<int>::const_iterator pos = lower_bound(pixel_mapping_arraytohigh.begin(), pixel_mapping_arraytohigh.end(), hrindex);
        planck_assert( *pos == hrindex, "High resolution pixel index not found in the partial map" );
        return partmap[pos-pixel_mapping_arraytohigh.begin()];
    }
    
    Healpix_Map<T> to_Healpix( T default_val=Healpix_undef ){
        Healpix_Map<T> outMap( order_, scheme_ );
        outMap.fill(default_val);
        for( int i=0; i<npartpix_; ++i ){
            outMap[pixel_mapping_arraytohigh[i]] = partmap[i];
        }
        return outMap;
    }
    
    /*! Returns the number of EXISTING pixels of the object. */
    int Npartpix() const { return npartpix_; }
    
    /* Returns the pixel number for the array index i
       Useful for iterating through the Partpix map */
    int highResPix(int i) const { return pixel_mapping_arraytohigh[i]; };


    /*! Deletes the old map, creates a map from the contents of \a data and
        sets the ordering scheme to \a scheme. The size of \a data must be a
        valid HEALPix map size.
        \note On exit, \a data is zero-sized! */
    // void Set (arr<T> &data, Healpix_Ordering_Scheme scheme)
    //   {
    //   Healpix_Base::SetNside(npix2nside (data.size()), scheme);
    //   partmap.transfer(data);
    //   }

    /*! Deletes the old map and creates a new map  with a given \a order
        and the ordering scheme \a scheme. */
    // void Set (int order, Healpix_Ordering_Scheme scheme)
    //   {
    //   Healpix_Base::Set(order, scheme);
    //   partmap.alloc(npix_);
    //   }
    /*! Deletes the old map and creates a new map  with a given \a nside
        and the ordering scheme \a scheme. */
    // void SetNside (int nside, Healpix_Ordering_Scheme scheme)
    //   {
    //   Healpix_Base::SetNside(nside, scheme);
    //   partmap.alloc(npix_);
    //   }

    /*! Fills the map with \a val. */
    void fill (const T &val)
      { partmap.fill(val); }

    /*! Imports the map \a orig into the current map, adjusting the
        ordering scheme. \a orig must have the same resolution as the
        current map. */
    void Import_nograde (const Partpix_Map<T> &orig, const Healpix_Map<T>& resolutionMask)
      {
      planck_assert (nside_==orig.nside_,
        "Import_nograde: maps have different nside");
      if (orig.scheme_ == scheme_)
        for (int m=0; m<npartpix_; ++m){
            partmap[m] = orig.partmap[m];
        }
      else
        {
        swapfunc swapper = (scheme_ == NEST) ?
          &Healpix_Base::ring2nest : &Healpix_Base::nest2ring;
//#pragma omp parallel
{
        int m;
        int arrayindex=0;
//#pragma omp for schedule (dynamic,5000)
        for (m=0; m<npix_; ++m){ 
            if( resolutionMask[resolutionMask.ang2pix(Healpix_Base::pix2ang(m))] > 0.5 ){
                pixel_mapping_arraytohigh[arrayindex] = (this->*swapper)(m);
                partmap[arrayindex] = orig.partmap_at_highresindex(m);
                ++arrayindex;
            }
        }
}
        }
      }

    /*! Imports the map \a orig into the current map, adjusting the
        ordering scheme and the map resolution. \a orig must have higher
        resolution than the current map. */
    void Import_upgrade (const Partpix_Map<T> &orig, const Healpix_Map<T>& resolutionMask)
      {
      planck_assert(nside_>orig.nside_,"Import_upgrade: this is no upgrade");
      int fact = nside_/orig.nside_;
      planck_assert (nside_==orig.nside_*fact,
        "the larger Nside must be a multiple of the smaller one");
      pix2xyf to_xyf = (orig.scheme_==RING) ?
        &Healpix_Base::ring2xyf : &Healpix_Base::nest2xyf;
      xyf2pix from_xyf = (scheme_==RING) ?
        &Healpix_Base::xyf2ring : &Healpix_Base::xyf2nest;
      pix2xyf res_to_xyf = (resolutionMask.Scheme()==RING) ?
        &Healpix_Base::ring2xyf : &Healpix_Base::nest2xyf;


    int arrayindex = 0;
    int factr = orig.nside_/resolutionMask.Nside();
    for( int mr=0; mr<resolutionMask.Npix(); ++mr ){
        if( resolutionMask[mr] < 0.5 )
            continue;
        int xr,yr,fr;
        (resolutionMask.*res_to_xyf)(mr,xr,yr,fr);
        int jrmax = factr*(yr+1);
        int irmax = factr*(xr+1);
        for (int jr=factr*yr; jr<jrmax; ++jr){
          for (int ir=factr*xr; ir<irmax; ++ir){
              int origpix = (orig.*from_xyf)(ir,jr,fr);
              // int x,y,f;
              //(orig.*to_xyf)(origpix,x,y,f);
              // planck_assert( x == ir, "x and ir are not the same" );
              // int jmax = fact*(y+1);
              // int imax = fact*(x+1);
              int jmax = fact*(jr+1);
              int imax = fact*(ir+1);
              for (int j=fact*jr; j<jmax; ++j){
                for (int i=fact*ir; i<imax; ++i){
                  int mypix = (this->*from_xyf)(i,j,fr);
                  pixel_mapping_arraytohigh[arrayindex] = mypix;
                  partmap[arrayindex] = orig.partmap_at_highresindex(origpix);
                  ++arrayindex;
                }
              }
          }
        }
      }
      planck_assert( arrayindex == npartpix_, "array index != npartpix_" );  
      sort_partmap( pixel_mapping_arraytohigh, partmap );
    } // Import_upgrade
    
/*            
//#pragma omp parallel
{
      int arrayindex=0;
//#pragma omp for schedule (dynamic,5000)
      for (int m=0; m<orig.npix_; ++m){
        int x,y,f;
        (orig.*to_xyf)(m,x,y,f);
        for (int j=fact*y; j<fact*(y+1); ++j){
          for (int i=fact*x; i<fact*(x+1); ++i){
            int mypix = (this->*from_xyf)(i,j,f);
            if( resolutionMask[resolutionMask.ang2pix(Healpix_Base::pix2ang(mypix))] > 0.5 ){
                pixel_mapping_arraytohigh[arrayindex] = mypix;
                partmap[arrayindex] = orig.partmap_at_highresindex(m);
                ++arrayindex;
            }

          }
        }
      }
}
    }
*/

    /*! Imports the partmap \a orig into the current map, adjusting the
        ordering scheme and the map resolution. \a orig must have higher
        resolution than the current map.
        \a pessimistic determines whether or not
        pixels are set to \a Healpix_undef when not all of the corresponding
        high-resolution pixels are defined.

        This method is instantiated for \a float and \a double only. */
    void Import_degrade (const Partpix_Map<T> &orig, const Healpix_Map<double>& resolutionMask, bool pessimistic=false, double cutoff=Healpix_undef);
    void Import_degrade (const Partpix_Map<T> &orig, const Healpix_Map<int>& resolutionMask, bool pessimistic=false, double cutoff=Healpix_undef);

    /*! Imports the map \a orig into the current map, adjusting the
        ordering scheme and the map resolution if necessary.
        When downgrading, \a pessimistic determines whether or not
        pixels are set to \a Healpix_undef when not all of the corresponding
        high-resolution pixels are defined.

        This method is instantiated for \a float and \a double only. */
    // void Import (const Partpix_Map<T> &orig, bool pessimistic=false)
    //   {
    //   if (orig.nside_ == nside_) // no up/degrading
    //     Import_nograde(orig);
    //   else if (orig.nside_ < nside_) // upgrading
    //     Import_upgrade(orig);
    //   else
    //     Import_degrade(orig,pessimistic);
    //   }

    /*! Returns a constant reference to the pixel with the number \a pix. */
    const T &operator[] (int pix) const { return partmap_ref_at_highresindex(pix); }
    /*! Returns a reference to the pixel with the number \a pix. */
    T &operator[] (int pix) { return partmap_ref_at_highresindex(pix); }

    /*! Swaps the map ordering from RING to NEST and vice versa.
        This is done in-place (i.e. with negligible space overhead). */
    // void swap_scheme()
    //   {
    //   static const int clen[] = { 0,7,5,4,12,10,13,18,14,19,18,17,27,21 };
    //   static const int cycle[][30] = {
    //     { },
    //     { 0,1,8,12,16,21,40 },
    //     { 0,1,2,40,114 },
    //     { 0,4,160,263 },
    //     { 0,4,30,49,51,87,526,1027,1105,1387,1807,2637 },
    //     { 0,8,10,18,39,74,146,307,452,4737 },
    //     { 0,1,2,7,9,17,80,410,1526,1921,32859,33566,38931 },
    //     { 0,5,6,10,12,24,27,95,372,494,924,1409,3492,4248,9137,66043,103369,
    //       156899 },
    //     { 0,1,2,3,4,45,125,351,697,24337,102940,266194,341855,419857 },
    //     { 0,1,2,3,9,16,1705,2082,2126,8177,12753,15410,52642,80493,83235,
    //       88387,99444,1675361,2495125 },
    //     { 0,2,6,8,9,11,20,50,93,152,183,2137,13671,44794,486954,741908,
    //       4803258,5692573 },
    //     { 0,1,5,6,44,53,470,2847,3433,4906,13654,14710,400447,1797382,
    //       2744492,18775974,23541521 },
    //     { 0,4,9,10,16,33,83,117,318,451,5759,10015,128975,171834,211256,
    //       347608,1278690,2154097,2590798,3427694,5581717,21012301,27023976,
    //       72522811,95032729,139166747,171822389 },
    //     { 0,5,10,267,344,363,2968,3159,9083,18437,76602,147614,1246902,
    //       1593138,2035574,6529391,9511830,11340287,29565945,281666026,
    //       677946848 } };
    // 
    //   swapfunc swapper = (scheme_ == NEST) ?
    //     &Healpix_Base::ring2nest : &Healpix_Base::nest2ring;
    // 
    //   planck_assert (order_>=0, "swap_scheme(): need hierarchical map");
    // 
    //   for (int m=0; m<clen[order_]; ++m)
    //     {
    //     int istart = cycle[order_][m];
    // 
    //     T pixbuf = partmap[istart];
    //     int iold = istart, inew = (this->*swapper)(istart);
    //     while (inew != istart)
    //       {
    //       partmap[iold] = partmap[inew];
    //       iold = inew;
    //       inew = (this->*swapper)(inew);
    //       }
    //     partmap[iold] = pixbuf;
    //     }
    //   scheme_ = (scheme_==RING) ? NEST : RING;
    //   }

    /*! performs the actual interpolation using \a pix and \a wgt. */
    T interpolation (const fix_arr<int,4> &pix,
      const fix_arr<double,4> &wgt) const
      {
      return T(partmap_at_highresindex(pix[0])*wgt[0] + partmap_at_highresindex(pix[1])*wgt[1]
             + partmap_at_highresindex(pix[2])*wgt[2] + partmap_at_highresindex(pix[3])*wgt[3]);
      }
    /*! Returns the interpolated map value at \a ptg */
    // T interpolated_value (const pointing &ptg) const
    //   {
    //   fix_arr<int,4> pix;
    //   fix_arr<double,4> wgt;
    //   get_interpol (ptg, pix, wgt);
    //   return interpolation (pix, wgt);
    //   }

    /*! Returns a constant reference to the map data. */
    const arr<T> &Map() const { return partmap; }

    /*! Returns the minimum and maximum value of the map in
        \a Min and \a Max.

        This method is instantiated for \a float and \a double only. */
    void minmax (T &Min, T &Max) const;

    /*! Swaps the contents of two Partpix_Map objects. */
    // void swap (Partpix_Map &other)
    //   {
    //   Healpix_Base::swap(other);
    //   partmap.swap(other.partmap);
    //   }

    /*! Returns the average of all defined map pixels. */
    double average() const
      {
      double avg=0;
      int pix=0;
      for (int m=0; m<npartpix_; ++m)
        if (!approx<double>(partmap[m],Healpix_undef))
          { ++pix; avg+=partmap[m]; }
      return avg/pix;
      }

    /*! Adds \a val to all defined map pixels. */
    void Add (T val)
      {
      for (int m=0; m<npartpix_; ++m)
        if (!approx<double>(partmap[m],Healpix_undef))
          { partmap[m]+=val; }
      }

    /*! Multiplies all defined map pixels by \a val. */
    void Scale (T val)
      {
      for (int m=0; m<npartpix_; ++m)
        if (!approx<double>(partmap[m],Healpix_undef))
          { partmap[m]*=val; }
      }

    /*! Returns the root mean square of the map, not counting undefined
        pixels. */
    double rms() const
      {
      using namespace std;

      double result=0;
      int pix=0;
      for (int m=0; m<npartpix_; ++m)
        if (!approx<double>(partmap[m],Healpix_undef))
          { ++pix; result+=partmap[m]*partmap[m]; }
      return sqrt(result/pix);
      }
    /*! Returns the maximum absolute value in the map, ignoring undefined
        pixels. */
    T absmax() const
      {
      using namespace std;

      T result=0;
      for (int m=0; m<npartpix_; ++m)
        if (!approx<double>(partmap[m],Healpix_undef))
          { result = max(result,abs(partmap[m])); }
      return result;
      }
      


    // sort the first array, but rearrange the second array in the same manner
    // heapsort from numerical recipes
    void sift_down( vector<int>& pixel_mapping, arr<T>& map, const int l, const int r ){
        int a = pixel_mapping[l];
        T a2 = map[l];
        int jold = l;
        int j = 2*l+1;
        while( j <= r ){
            if( j < r && pixel_mapping[j] < pixel_mapping[j+1] )
                ++j;
            if( a >= pixel_mapping[j] )
                break;
            pixel_mapping[jold] = pixel_mapping[j];
            map[jold] = map[j];
            jold = j;
            j = 2*j+1;
        }
        pixel_mapping[jold] = a;
        map[jold] = a2;
    }
    void sort_partmap( vector<int>& pixel_mapping, arr<T>& map ){
        int i;
        int n = pixel_mapping.size();
        for( i=n/2-1; i>=0; --i )
            sift_down( pixel_mapping, map, i, n-1 );
        for( i=n-1; i>0; --i ){
            int tempi = pixel_mapping[0];
            pixel_mapping[0] = pixel_mapping[i];
            pixel_mapping[i] = tempi;
            T temp = map[0];
            map[0] = map[i];
            map[i] = temp;
            sift_down( pixel_mapping, map, 0, i-1 );
        }

        // check the sort
        for( i=0; i<n-1; ++i){
            if( pixel_mapping[i+1] <= pixel_mapping[i] ){
                fprintf( stderr, "pixel %d = %d, %d+1 = %d !\n", i, pixel_mapping[i], i, pixel_mapping[i+1] );
                throw;
            }
        }

    }
    
          
      // // sort the first array, but rearrange the second array in the same manner
      // // from stroustrup
      // void sort_partmap( vector<int>& pixel_mapping, arr<T>& map ){
      //     int n = pixel_mapping.size();
      //     for( int gap=n/2; 0<gap; gap/=2 )
      //         for( int i=gap; i<n; i++ )
      //             for( int j=i-gap; 0<=j; j-=gap ){
      //                 if( pixel_mapping[j+gap] < pixel_mapping[j] ){
      //                     int tempi = pixel_mapping[j];
      //                     pixel_mapping[j] = pixel_mapping[j+gap];
      //                     pixel_mapping[j+gap] = tempi;
      //                     T temp = map[j];
      //                     map[j] = map[j+gap];
      //                     map[j+gap] = temp;
      //                 }
      //             }
      // }
      
  };

#endif
