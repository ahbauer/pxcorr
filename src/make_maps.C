
#include <unistd.h>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include <sstream>
using std::stringstream;
#include <fstream>
using std::ifstream;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <iomanip>
using namespace std;

#include <algorithm>
using std::transform;
#include <ctype.h>
using std::tolower;

#include <healpix_base.h>
#include <healpix_base2.h>
#include <healpix_map.h>
#include <partpix_map.h>
#include <partpix_map2.h>
#include <fitsio.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>
#include <healpix_data_io.h>

#include "H5Cpp.h"
using namespace H5;

using namespace std;

void write_map( Partpix_Map2<float> *dcMap, Partpix_Map2<int> *dcMask, string outfilename, string ftype, float *ang_means_c, float *ang_widths_c, int n_ang_bins );

void read_mask( string mask_filename, vector<int>& mask_pixels, Healpix_Base2& mask_base );

void make_maps( const char *catalog_filename_c, const char *mask_filename_c, float *ang_means_c, int n_ang_bins0, float *ang_widths_c, int n_ang_bins, bool use_counts, bool use_mags, char *suffix_c ){

    //int n_ang_bins = ang_means.size();
    double pixel_multiple = 2.5;

    cerr << "some sizes: " << sizeof(hsize_t) << " " << sizeof(int) << " " << sizeof(int64) << " " << sizeof(float) << " " << sizeof(double) << endl;
    cerr << H5::PredType::NATIVE_ULONG.getSize() << " " << H5::PredType::NATIVE_LLONG.getSize() << " " << H5::PredType::NATIVE_FLOAT.getSize() << endl;

    cerr << setprecision(16);
    
    if( !( use_counts || use_mags) ){
        cerr << "Not using either counts or mags... doing nothing!" << endl;
        throw exception();
    }

    string catalog_filename(catalog_filename_c);
    string mask_filename(mask_filename_c);
    string suffix(suffix_c);
    
    cerr << "make_maps " << catalog_filename << " " << mask_filename << endl;

    int min_footprint_order = 3;
    
    // From the smallest angular bin width, figure out the necessary map resolution.
    float min_width = 1.e9;
    for( int i=0; i<n_ang_bins; ++i ){
        if( ang_widths_c[i] < min_width )
            min_width = ang_widths_c[i];
    }
    int order = 1;
    while(1){
      float pixel_size = sqrt(41253./(12.0*pow(pow(2.0, order), 2.0)));
      if( pixel_multiple*pixel_size < min_width )  // this 0.85 is a bit arbitrary
          break;
      ++order;
    }
    cerr << "Using order " << order << " for minimum angular width " << min_width << endl;
    
    Healpix_Map<int> *footprintMap;
    int mask_order;
    vector<int> mask_pixels;
    Healpix_Base2 mask_base;
    Partpix_Map2<int> *dcMask;
    
    if( catalog_filename.rfind(".fits") != string::npos ){
        
        cerr << "make_maps: Found a Healpix map as input!" << endl;

        if( ! (use_counts == true && use_mags == false) ){
            cerr << "make_maps: Healpix input, but not use_counts=1 and use_mags=0... this doesn't make sense." << endl;
            throw exception();
        }

        Healpix_Map<double> *inputmap = new Healpix_Map<double>(11, RING);
        cerr << "made a new empty healpix map" << endl;
            // read_Healpix_map_from_fits( catalog_filename, *inputmap );
        fitshandle inp;
        cerr << "declared fitshandle" << endl;
        inp.open (catalog_filename);
        cerr << "opened input file" << endl;
        inp.goto_hdu (2);
        cerr << "went to hdu 2" << endl;
        // read_Healpix_map_from_fits( inp, *inputmap, 1 );
        string ordering;
        inp.get_key("ORDERING", ordering);
        cerr << "ordering = " << ordering << endl;
        arr<double> myarr;
        cerr << "reading entire column" << endl;
        inp.read_entire_column(1,myarr);
        cerr << "read in the entire column!" << endl;
        inputmap->Set(myarr, ordering=="RING" ? RING : NEST);
        cerr << "read in the map from file" << endl;
        if( inputmap->Scheme() == NEST )
            inputmap->swap_scheme();
        
        // if there's no extra mask, use all the non-zero healpix pixel as the good area, and to find the footprint.
        if( mask_filename == "None" || mask_filename == "none" || mask_filename == "NONE" || mask_filename == "null" || mask_filename == "NULL" ){
        
            // what's the footprint?
            int footprint_order = min_footprint_order;
            double footprint_area = 50000.;
            Healpix_Map<int> *footprintMap = new Healpix_Map<int>(1, RING);
            while(footprint_order <= inputmap->Order()){
                delete footprintMap;
                footprintMap = new Healpix_Map<int>(footprint_order, RING);
                footprintMap->fill(0);
                for( int i=0; i<inputmap->Npix(); ++i ){
                    if( (*inputmap)[i] != 0. && (*inputmap)[i] != Healpix_undef && ! isnan((*inputmap)[i]) )
                        (*footprintMap)[ footprintMap->ang2pix(inputmap->pix2ang(i)) ] = 1;
                }
                double footprint_area_new = 0.;
                for( int i=0; i<footprintMap->Npix(); ++i ){
                  if( (*footprintMap)[i] == 1 )
                      footprint_area_new += 1.0;
                }
                footprint_area_new *= (41253./footprintMap->Npix());
                if( footprint_area_new <= 0.8*footprint_area ){
                    footprint_area = footprint_area_new;
                    cerr << "Footprint order " << footprint_order << " uses area " << footprint_area << " sq degrees... trying the next order." << endl;
                    ++footprint_order;
                }
                else{
                    break;
                }
            }
            // the previous footprint order was fine.
            // we want this to be small because it's also the minimum jackknife order.
            // so, go back to the last footprint order.
            --footprint_order;
            delete footprintMap;
            footprintMap = new Healpix_Map<int>(footprint_order, RING);
            footprintMap->fill(0);
            for( int i=0; i<inputmap->Npix(); ++i ){
                (*footprintMap)[ footprintMap->ang2pix(inputmap->pix2ang(i)) ] = 1;
            }
            cerr << "Using footprint order " << footprintMap->Order() << " with " << footprint_area << " sq degrees." << endl;
        
            // Now fill in the mask
            dcMask = new Partpix_Map2<int>(inputmap->Order(), *footprintMap);
            dcMask->fill(0);
            for( int i=0; i<inputmap->Npix(); ++i ){
                if( (*inputmap)[i] != 0. && (*inputmap)[i] != Healpix_undef && ! isnan((*inputmap)[i]) ){
                    (*dcMask)[i] = 1;
                }
            }
            cerr << "Finished filling in the mask from the Healpix input itself" << endl;
            
        } // if there's no mask to go with the healpix map
        
        else{
            
            string scheme;
            
            // read in the mask
            read_mask( mask_filename, mask_pixels, mask_base );
            mask_order = mask_base.Order();
            
            // from the masks, determine the best possible healpix footprint
            // requirements: must cover all of the mask area
            // must cover <80% of the area covered by the next lowest resolution
            // must be order <= 6 ?
            int footprint_order = min_footprint_order;
            double footprint_area = 50000.;
            footprintMap = new Healpix_Map<int>(1, RING);
            while(footprint_order <= mask_order){
                delete footprintMap;
                footprintMap = new Healpix_Map<int>(footprint_order, RING);
                footprintMap->fill(0);
                for( unsigned int i=0; i<mask_pixels.size(); ++i ){
                    (*footprintMap)[ footprintMap->ang2pix(mask_base.pix2ang(mask_pixels[i])) ] = 1;
                }
                double footprint_area_new = 0.;
                for( int i=0; i<footprintMap->Npix(); ++i ){
                  if( (*footprintMap)[i] == 1 )
                      footprint_area_new += 1.0;
                }
                footprint_area_new *= (41253./footprintMap->Npix());
                if( footprint_area_new <= 0.8*footprint_area ){
                    footprint_area = footprint_area_new;
                    cerr << "Footprint order " << footprint_order << " uses area " << footprint_area << " sq degrees... trying the next order." << endl;
                    ++footprint_order;
                }
                else{
                    //footprint_area = footprint_area_new;
                    break;
                }
            }
            // the previous footprint order was fine.
            // we want this to be small because it's also the minimum jackknife order.
            // so, go back to the last footprint order.
            --footprint_order;
            delete footprintMap;
            footprintMap = new Healpix_Map<int>(footprint_order, RING);
            footprintMap->fill(0);
            for( unsigned int i=0; i<mask_pixels.size(); ++i ){
                (*footprintMap)[ footprintMap->ang2pix(mask_base.pix2ang(mask_pixels[i])) ] = 1;
            }
            cerr << "Using footprint order " << footprintMap->Order() << " with " << footprint_area << " sq degrees." << endl;
            
            // fill in the Partpix mask(s)
            dcMask = new Partpix_Map2<int>(mask_order, *footprintMap);
            dcMask->fill(0);
            for( unsigned int i=0; i<mask_pixels.size(); ++i ){
              if( ! (*footprintMap)[ footprintMap->ang2pix(mask_base.pix2ang(mask_pixels[i])) ] ){
                  cerr << "Skipping mask pixel " << i << "!" << endl;
                  continue;
              }
              (*dcMask)[mask_pixels[i]] = 1;
            }
            mask_pixels.clear();

            cerr << "Finished filling in the supplied mask" << endl;
            
        } // if there is a mask to go with the healpix map
        
        
        // now make the healpix map into a partpix map
        Partpix_Map2<float> *dcMap = new Partpix_Map2<float>(inputmap->Order(), *footprintMap);
        dcMap->fill(0.);
        for( int i=0; i<inputmap->Npix(); ++i ){
            if( (*inputmap)[i] != 0. && (*inputmap)[i] != Healpix_undef && ! isnan((*inputmap)[i]) ){
                if( (*footprintMap)[ footprintMap->ang2pix(inputmap->pix2ang(i)) ] )
                    (*dcMap)[i] = (*inputmap)[i];
            }
        }
        // free the memory
        delete inputmap;
        
        
        // make sure that the map is the same resolution as the mask
        if( dcMap->Order() > dcMask->Order() ){
            Partpix_Map2<float> *dcMap2 = new Partpix_Map2<float>(mask_order, *footprintMap);
            dcMap2->Import_degrade(*dcMap, *footprintMap);
            delete dcMap;
            dcMap = dcMap2;
        }
        else if( dcMap->Order() < dcMask->Order() ){
            Partpix_Map2<float> *dcMap2 = new Partpix_Map2<float>(mask_order, *footprintMap);
            dcMap2->Import_upgrade(*dcMap, *footprintMap);
            delete dcMap;
            dcMap = dcMap2;
        }
        
    
        // now turn the counts map into a delta map
        double mean_counts = 0.0;
        double n_pix = 0.0;
        for( int64 i1=0; i1<dcMap->Npartpix(); ++i1 ){
            int64 i = dcMap->highResPix(i1);
            int64 maskpix = dcMask->ang2pix(dcMap->pix2ang(i));
            if( (*dcMask)[maskpix] > 0.5 ){
                mean_counts += (*dcMap)[i];
                n_pix += 1.0;
            }
        }
        if( n_pix == 0.0 ){
            cerr << "No unmasked pixels??" << endl;
            throw exception();
        }
        mean_counts /= n_pix;
        if( mean_counts == 0. ){
            cerr << "No counts in the unmasked part of the map??" << endl;
            throw exception();
        }
        for( int64 i1=0; i1<dcMap->Npartpix(); ++i1 ){
            int64 i = dcMap->highResPix(i1);
            (*dcMap)[i] = (*dcMap)[i]/mean_counts - 1.0;
        }
    
        // now make sure the file is the resolution necessary for the smallest angular bin
        // i realize that this is moot if the input res is lower, but it will keep the code from crashing.
        if( order > dcMap->Order() ){
            Partpix_Map2<float> *dcMap2 = new Partpix_Map2<float>(order, *footprintMap);
            Partpix_Map2<int> *dcMask2 = new Partpix_Map2<int>(order, *footprintMap);
            dcMap2->Import_upgrade(*dcMap, *footprintMap);
            dcMask2->Import_upgrade(*dcMask, *footprintMap);
            delete dcMap;
            delete dcMask;
            dcMap = dcMap2;
            dcMask = dcMask2;
        }
        else if( order < dcMap->Order() ){
            Partpix_Map2<float> *dcMap2 = new Partpix_Map2<float>(order, *footprintMap);
            Partpix_Map2<int> *dcMask2 = new Partpix_Map2<int>(order, *footprintMap);
            dcMap2->Import_degrade(*dcMap, *footprintMap);
            dcMask2->Import_degrade(*dcMask, *footprintMap);
            delete dcMap;
            delete dcMask;
            dcMap = dcMap2;
            dcMask = dcMask2;
        }
    
        // and write out the map
        string outfilename = "dc_map" + suffix + ".h5";
        string ftype = "counts";
        write_map( dcMap, dcMask, outfilename, ftype, ang_means_c, ang_widths_c, n_ang_bins );
    
        cerr << "Wrote map to file" << endl;
        delete dcMap;
        delete dcMask;
        
        // we're done!
        return;
        
    } // if healpix map input
    
    
    

    
    // read in the galaxy catalog
    FILE * file = fopen( catalog_filename.c_str(), "rb" );
    if( file == NULL ){
        cerr << "Problem opening " << catalog_filename << endl;
        throw exception();
    }

    vector<pointing> pointings;
    vector<double> mags;
    double invals[3];
    while( fread( invals, sizeof( double ), 3, file ) == 3 ){
        double phi = invals[0]*3.1415926/180.;
        double theta = (90.-invals[1])*3.1415926/180.;
        pointing mypointing(theta, phi);
        pointings.push_back( mypointing );
        if( use_mags )
            mags.push_back( invals[2] );
    }
    fclose( file );    
    cerr << "Read in " << pointings.size() << " objects" << endl;
    
    
    // read in the mask
    string scheme;
    if( mask_filename == "None" || mask_filename == "none" || mask_filename == "NONE" || mask_filename == "null" || mask_filename == "NULL" ){
        cerr << "Making a dummy mask file from the data itself.  This is not the best idea." << endl;
        mask_order = order-2;
        mask_base = Healpix_Base2(mask_order, RING);
        for( unsigned int p=0; p<pointings.size(); ++p ){
            mask_pixels.push_back( mask_base.ang2pix(pointings[p]) );
        }
        unique(mask_pixels.begin(), mask_pixels.end());
    }
    else{
        read_mask(mask_filename, mask_pixels, mask_base);
    }
    mask_order = mask_base.Order();
    cerr << "Finished with the mask: " << mask_pixels.size() << " good pixels with order " << mask_order << endl;
    
    // from the masks, determine the best possible healpix footprint
    // requirements: must cover all of the mask area
    // must cover <80% of the area covered by the next lowest resolution
    // must be order <= 6 ?
    int footprint_order = min_footprint_order;
    double footprint_area = 50000.;
    footprintMap = new Healpix_Map<int>(1, RING);
    while(footprint_order <= mask_order){
        delete footprintMap;
        footprintMap = new Healpix_Map<int>(footprint_order, RING);
        footprintMap->fill(0);
        for( unsigned int i=0; i<mask_pixels.size(); ++i ){
            (*footprintMap)[ footprintMap->ang2pix(mask_base.pix2ang(mask_pixels[i])) ] = 1;
        }
        double footprint_area_new = 0.;
        for( int i=0; i<footprintMap->Npix(); ++i ){
          if( (*footprintMap)[i] == 1 )
              footprint_area_new += 1.0;
        }
        footprint_area_new *= (41253./footprintMap->Npix());
        if( footprint_area_new <= 0.8*footprint_area ){
            footprint_area = footprint_area_new;
            cerr << "Footprint order " << footprint_order << " uses area " << footprint_area << " sq degrees... trying the next order." << endl;
            ++footprint_order;
        }
        else{
            //footprint_area = footprint_area_new;
            break;
        }
    }
    // the previous footprint order was fine.
    // we want this to be small because it's also the minimum jackknife order.
    // so, go back to the last footprint order.
    --footprint_order;
    delete footprintMap;
    footprintMap = new Healpix_Map<int>(footprint_order, RING);
    footprintMap->fill(0);
    for( unsigned int i=0; i<mask_pixels.size(); ++i ){
        (*footprintMap)[ footprintMap->ang2pix(mask_base.pix2ang(mask_pixels[i])) ] = 1;
    }
    cerr << "Using footprint order " << footprintMap->Order() << " with " << footprint_area << " sq degrees." << endl;

    // fill in the Partpix mask(s)
    dcMask = new Partpix_Map2<int>(mask_order, *footprintMap);
    dcMask->fill(0);
    Partpix_Map2<int> *dmMask;
    if( use_mags ){
        dmMask = new Partpix_Map2<int>(mask_order, *footprintMap);
        dmMask->fill(0);
    }
    for( unsigned int i=0; i<mask_pixels.size(); ++i ){
      if( ! (*footprintMap)[ footprintMap->ang2pix(mask_base.pix2ang(mask_pixels[i])) ] ){
          cerr << "Skipping mask pixel " << i << "!" << endl;
          continue;
      }
      (*dcMask)[mask_pixels[i]] = 1;
      if( use_mags )
        (*dmMask)[mask_pixels[i]] = 1;
    }
    mask_pixels.clear();
    
    // convert the mag mask to be the same resolution as the map
    if( use_mags ){
        if( order > dmMask->Order() ){
            cerr << "upgrading mag mask to order " << order << endl;
            Partpix_Map2<int> *dmMask2 = new Partpix_Map2<int>( order, *footprintMap );
            dmMask2->Import_upgrade(*dmMask, *footprintMap);
            delete(dmMask);
            dmMask = dmMask2;
        }
        else if( order < dmMask->Order() ){
            cerr << "degrading mag mask to order " << order << endl;
            Partpix_Map2<int> *dmMask2 = new Partpix_Map2<int>( order, *footprintMap );
            dmMask2->Import_degrade(*dmMask, *footprintMap);
            delete(dmMask);
            dmMask = dmMask2;
        }
    }
    
    cerr << "Finished filling in the mask(s)" << endl;
  
    // system( "rm fp.fits" );
    // fitshandle myfits;
    // myfits.create("fp.fits");
    // write_Healpix_map_to_fits(myfits, *footprintMap, PLANCK_FLOAT64);
    // myfits.close();
    // for( int f=0; f<footprintMap->Npix(); ++f ){
    //     if( (*footprintMap)[f] > 0.5 ){
    //         pointing mypointing = footprintMap->pix2ang(f);
    //         double ra = mypointing.phi*180./3.1415926;
    //         double dec = 90. - mypointing.theta*180./3.1415926;
    //         cout << ra << " " << dec << " " << f << endl;
    //     }
    // }
    
    // make the maps
        
    // always make the delta_counts map since we need that in order to make the delta_mag MASK
    Partpix_Map2<float> *dcMap = new Partpix_Map2<float>(order, *footprintMap);
    dcMap->fill(0.);
    Partpix_Map2<float> *dmMap;
    double mean_mag = 0.;
    int n_mags = 0;
    if( use_mags ){
        dmMap = new Partpix_Map2<float>(order, *footprintMap);
        dmMap->fill(0.);
    }
    for( unsigned int p=0; p<pointings.size(); ++p ){
        // check if objects are outside the footprint map
        // might happen if objects are outside the mask.
        int fppix = footprintMap->ang2pix(pointings[p]);
        if( (*footprintMap)[fppix] < 0.5 )
            continue;
        int64 pixnum = dcMap->ang2pix( pointings[p] );
        (*dcMap)[pixnum] += 1.0;
        if( use_mags ){
            int64 maskpix = dmMask->ang2pix(dmMap->pix2ang(pixnum));
            if( (*dmMask)[maskpix] > 0.5 ){
                (*dmMap)[pixnum] += mags[p];
                mean_mag += mags[p];
                n_mags++;
            }
        }
    }
    
    // now make the delta_mag mask for pixels with counts=0
    // mask out the pixels that don't have any mag info.
    if( use_mags ){
        mean_mag /= n_mags;
        for( int i1=0; i1<dmMap->Npartpix(); ++i1 ){
            int i = dmMap->highResPix(i1);
            if( (*dcMap)[i] > 0. ){
                (*dmMap)[i] = (*dmMap)[i]/(*dcMap)[i] - mean_mag;
            }
            else{
                // possible because the mask and map are the same resolution
                (*dmMask)[i] = 0.0;
            }
        }
    }
    
    // now turn the counts map into a delta map
    double mean_counts = 0.0;
    double n_pix = 0.0;
    for( int64 i1=0; i1<dcMap->Npartpix(); ++i1 ){
        int64 i = dcMap->highResPix(i1);
        int64 maskpix = dcMask->ang2pix(dcMap->pix2ang(i));
        if( (*dcMask)[maskpix] > 0.5 ){
            mean_counts += (*dcMap)[i];
            n_pix += 1.0;
        }
    }
    if( n_pix == 0.0 ){
        cerr << "No unmasked pixels??" << endl;
        throw exception();
    }
    mean_counts /= n_pix;
    if( mean_counts == 0. ){
        cerr << "No counts in the unmasked part of the map??" << endl;
        throw exception();
    }
    for( int64 i1=0; i1<dcMap->Npartpix(); ++i1 ){
        int64 i = dcMap->highResPix(i1);
        (*dcMap)[i] = (*dcMap)[i]/mean_counts - 1.0;
    }

    cerr << "Writing out the maps" << endl;
    if( use_counts ){
        string outfilename = "dc_map" + suffix + ".h5";
        string ftype = "counts";
        write_map( dcMap, dcMask, outfilename, ftype, ang_means_c, ang_widths_c, n_ang_bins );
        
        cerr << "Wrote delta counts map to file" << endl;
        delete dcMap;
        delete dcMask;
    }
    if( use_mags ){
        string outfilename = "dm_map" + suffix + ".h5";
        string ftype = "mag";
        write_map( dmMap, dmMask, outfilename, ftype, ang_means_c, ang_widths_c, n_ang_bins );
        
        cerr << "Wrote delta mag map to file" << endl;
        delete dmMask;
        delete dmMap;
    }
    
    delete footprintMap;

    cerr << "Finished with make_maps!" << endl;

    return;
}


void read_mask( string mask_filename, vector<int>& mask_pixels, Healpix_Base2& mask_base ){
    
    double input_mask_cut = 1.e-5;
        
    FILE * maskfile = fopen( mask_filename.c_str(), "r" );
    if( maskfile == NULL ){
        cerr << "Problem opening " << mask_filename << endl;
        throw exception();
    }
    bool header = true;
    string scheme;
    int mask_order;
    char line_char[256];
    while( fgets( line_char, 256, maskfile ) ){
        if( header ){
            char sch[16];
            int num = sscanf(line_char, "%s %d", sch, &mask_order );
            if( num != 2 ){
                cerr << "Problem reading header from mask file, number of arguments read = " << num << ", not 2" << endl;
                throw exception();
            }
            scheme = string(sch);
            if( (!scheme.compare("RING")) && (!scheme.compare("NEST")) ){
                cerr << "Scheme listed in mask header line is not RING or NEST, but " << scheme << endl;
                throw exception();
            }
            header = false;
            continue;
        }
        int pix;
        double val;
        int num = sscanf(line_char, "%d %lf", &pix, &val );
        if( num != 2 ){
            // allow the implicit assumption of a pixel's existence in the list meaning that the mask is good there.
            val = 1.0;
            num = sscanf(line_char, "%d", &pix );
            if( num != 1 ){
                cerr << "Problem reading from mask file, number of arguments read = " << num << ", not 1" << endl;
                throw exception();
            }
        }
        if( val > input_mask_cut )
            mask_pixels.push_back(pix);
    }
    fclose(maskfile);
    cerr << "Read in mask information, scheme " << scheme << ", order " << mask_order << ", " << mask_pixels.size() << " pixels" << endl;
    mask_base = Healpix_Base2( mask_order, RING );
    vector<int> mask_pixels_nest( mask_pixels );
    if( scheme.compare("NEST") == 0 ){
        cerr << "Changing mask to RING scheme... ";
        for( unsigned int i=0; i<mask_pixels.size(); ++i ){
            mask_pixels[i] = mask_base.nest2ring(mask_pixels_nest[i]);
        }
        scheme = "RING";
        cerr << "Done!" << endl;
    }
    sort( mask_pixels.begin(), mask_pixels.end() );
    mask_pixels_nest.clear();
    
    return;
    
}



void write_map( Partpix_Map2<float> *dcMap, Partpix_Map2<int> *dcMask, string outfilename, string ftype, float *ang_means_c, float *ang_widths_c, int n_ang_bins ){
    
    // write out the maps (hdf5)
    cerr << "writing " << outfilename << endl;

    H5::H5File *file = new H5::H5File( outfilename, H5F_ACC_TRUNC ); // clobber!
    H5::Group* group = new H5::Group( file->createGroup( "/data" ));
    H5::PredType datatype_map( H5::PredType::NATIVE_FLOAT );
    H5::PredType datatype_pix( H5::PredType::NATIVE_ULONG );
    H5::PredType datatype_mask( H5::PredType::NATIVE_ULONG );
    // hsize_t dimsf[2];
    // dimsf[1] = 2;
    // int rank = 2;
    hsize_t dimsf[1];
    int rank = 1;

    // cerr << "sizes " << H5::PredType::NATIVE_FLOAT.getSize() << " " << sizeof(H5::PredType::NATIVE_FLOAT) << " " << sizeof(double) << " " << sizeof(float) << endl;

    // the map dataset
    int64 dsize = dcMap->Npartpix()+1;
    float *data = new float[dsize];
    unsigned long *pix = new unsigned long[dsize];
    pix[0] = dcMap->Order();
    data[0] = RING;
    int64 index = 1;
    for( int64 i1=0; i1<dcMap->Npartpix(); ++i1 ){
        unsigned long i = dcMap->highResPix(i1);
        // only write it out if it's inside the mask
        if( (*dcMask)[dcMask->ang2pix(dcMap->pix2ang(i))] > 0.5 ){
            pix[index] = i;
            data[index] = (*dcMap)[i];
            ++index;
        }
    }
    dimsf[0] = index;
    cerr << "length of array to save = " << index << endl;
    H5::DataSpace *dataspace1 = new H5::DataSpace( rank, dimsf );
    H5::DataSpace *dataspace2 = new H5::DataSpace( rank, dimsf );
    string datasetname1 = "/data/pixels";
    string datasetname2 = "/data/values";
    H5::DataSet *dataset1 = new H5::DataSet( file->createDataSet( datasetname1, datatype_pix, *dataspace1 ) );
    H5::DataSet *dataset2 = new H5::DataSet( file->createDataSet( datasetname2, datatype_map, *dataspace2 ) );
    
    cerr << "Writing the data" << endl;
    cerr << "Values size " << H5::PredType::NATIVE_FLOAT.getSize() << "x" << dimsf[0] << endl;
    dataset2->write( data, H5::PredType::NATIVE_FLOAT );
    
    cerr << "Pixels size " << H5::PredType::NATIVE_ULONG.getSize() << "x" << dimsf[0] << endl;
    dataset1->write( pix, H5::PredType::NATIVE_ULONG );

    cerr << "Deleting info" << endl;

    delete[] data;
    delete[] pix;
    delete dataspace1;
    delete dataspace2;
    delete dataset1;
    delete dataset2;
    
    cerr << "Wrote " << dimsf[0] << "-1 pixels to the map" << endl;

    // the mask dataset
    unsigned long *data_mask = new unsigned long[dcMask->Npartpix()+1];
    unsigned long *pix_mask = new unsigned long[dcMask->Npartpix()+1];
    pix_mask[0] = dcMask->Order();
    data_mask[0] = RING;
    index = 1;
    int nzeros = 0;
    for( int64 i1=0; i1<dcMask->Npartpix(); ++i1 ){
        int64 i = dcMask->highResPix(i1);
        pix_mask[index] = i;
        data_mask[index] = (*dcMask)[i];
        if( data_mask[index] == 0 )
            ++nzeros;
        ++index;
    }
    dimsf[0] = dcMask->Npartpix()+1;
    H5::DataSpace *dataspace_pix = new H5::DataSpace( rank, dimsf );
    H5::DataSpace *dataspace_mask = new H5::DataSpace( rank, dimsf );
    datasetname1 = "/data/mask_pixels";
    datasetname2 = "/data/mask_values";
    H5::DataSet *dataset_mask1 = new H5::DataSet( file->createDataSet( datasetname1, datatype_mask, *dataspace_pix ) );
    H5::DataSet *dataset_mask2 = new H5::DataSet( file->createDataSet( datasetname2, datatype_mask, *dataspace_mask ) );
    dataset_mask1->write( pix_mask, H5::PredType::NATIVE_ULONG );
    dataset_mask2->write( data_mask, H5::PredType::NATIVE_ULONG );
    cerr << "Wrote " << dcMask->Npartpix() << " pixels to the mask, " << dcMask->Npartpix()-nzeros << " non-zero." << endl;
    
    delete dataspace_pix;
    delete dataspace_mask;
    delete dataset_mask1;
    delete dataset_mask2;
    delete[] pix_mask;
    delete[] data_mask;

    // the metadata
    H5::Group *metagroup = new H5::Group( file->createGroup( "/meta" ) );
    H5::PredType metadatatype( H5::PredType::NATIVE_CHAR );
    //H5::DataType metadatatype( H5T_STRING, 1 );
    //rank = 1;
    dimsf[0] = 1;
    H5::DataSpace *metadataspace = new H5::DataSpace( rank, dimsf );
    H5::DataSet *metadataset = new H5::DataSet( file->createDataSet( "/meta/meta", metadatatype, *metadataspace ) );

    string ang_mean_line = "[";
    char buffer[32];
    sprintf( buffer, "%f", ang_means_c[0] );
    ang_mean_line.append(buffer);
    for( int i=1; i<n_ang_bins; ++i ){
        ang_mean_line.append(", ");
        sprintf( buffer,"%f", ang_means_c[i] );
        ang_mean_line.append(buffer);
    }
    ang_mean_line.append("]");

    // Create new dataspace for attribute
    H5::DataSpace attr1_dataspace = H5::DataSpace( H5S_SCALAR );

    // Create new string datatype for attribute
    H5::StrType strdatatype( H5::PredType::C_S1, ang_mean_line.size()+1 ); // of length 256 characters

    // Set up write buffer for attribute
    const H5std_string u_mean_buf( ang_mean_line.c_str() );

    // Create attribute and write to it
    H5::Attribute u_mean_attr = metadataset->createAttribute("u_mean", strdatatype, attr1_dataspace);
    u_mean_attr.write(strdatatype, u_mean_buf);

    string ang_widths_line = "[";
    sprintf( buffer, "%f", ang_widths_c[0] );
    ang_widths_line.append(buffer);
    for( int i=1; i<n_ang_bins; ++i ){
        ang_widths_line.append(", ");
        sprintf( buffer,"%f", ang_widths_c[i] );
        ang_widths_line.append(buffer);
    }
    ang_widths_line.append("]");

    // Create new dataspace for attribute
    H5::DataSpace attr2_dataspace = H5::DataSpace( H5S_SCALAR );

    // Create new string datatype for attribute
    strdatatype = H5::StrType( H5::PredType::C_S1, ang_widths_line.size()+1 ); // of length 256 characters

    // Set up write buffer for attribute
    const H5std_string u_widths_buf( ang_widths_line.c_str() );

    // Create attribute and write to it
    H5::Attribute u_width_attr = metadataset->createAttribute("u_width", strdatatype, attr2_dataspace);
    u_width_attr.write(strdatatype, u_widths_buf);

    // Create new dataspace for attribute
    H5::DataSpace attr3_dataspace = H5::DataSpace( H5S_SCALAR );
    // Create new string datatype for attribute
    strdatatype = H5::StrType( H5::PredType::C_S1, 6 ); // of length 256 characters
    // Set up write buffer for attribute
    const H5std_string fourier_buf( "false" );
    // Create attribute and write to it
    H5::Attribute fourier_attr = metadataset->createAttribute("fourier", strdatatype, attr3_dataspace);
    fourier_attr.write(strdatatype, fourier_buf);


    // Create new dataspace for attribute
    H5::DataSpace attr4_dataspace = H5::DataSpace( H5S_SCALAR );
    // Create new string datatype for attribute
    strdatatype = H5::StrType( H5::PredType::C_S1, 11 ); // of length 256 characters    
    // Set up write buffer for attribute
    string ftype2 = "[\"" + ftype + "\"]";
    const H5std_string ftype_buf( ftype2.c_str() );
    // Create attribute and write to it
    H5::Attribute ftype_attr = metadataset->createAttribute("ftype", strdatatype, attr4_dataspace);
    ftype_attr.write(strdatatype, ftype_buf);

    delete metadataspace;
    delete metadataset;

    delete group;
    delete metagroup;
    delete file;

} // write_map


int main (int argc, char **argv){

    if( argc < 5 ){
        cerr << "Usage: make_maps metadata catalog mask counts mags" << endl;
        cerr << "       for making maps for both counts or mags, or you can just specify one." << endl;
        cerr << "Example input is given in the pxcorr/example_input directory." << endl;
        return 1;
    }
    
    string metadatafilename(argv[1]);
    string catalog_filename(argv[2]);
    string mask_filename(argv[3]);

    // read in the instructions about which maps to make
    int index = 4;
    bool use_counts = false;
    bool use_mags = false;
    while( index < argc ){
        string arg(argv[index]);
        if( arg.find("counts") != string::npos )
            use_counts = true;
        if( arg.find("mags") != string::npos )
            use_mags = true;
        ++index;
    }
    cerr << "use counts: " << use_counts << endl;
    cerr << "use mags: " << use_mags << endl;
    if( !( use_counts || use_mags) ){
        cerr << "Not using either counts or mags... doing nothing!" << endl;
        return 1;
    }
    
    // read in the binning metadata
    ifstream metadatafile( metadatafilename.c_str() );
    vector<float> ang_means;
    vector<float> ang_widths;
    string line;
    getline( metadatafile, line );
    while( !metadatafile.eof()){
        stringstream linestream(line);
        
        // look for ang_mean
        size_t index1 = line.find("ang_mean");
        if(index1 != string::npos){
            ++index1;
            while( index1 < line.length() && line[index1] != '[' ){
                ++index1;
            }
            if( line[index1] != '[' ){
                cerr << "Problem parsing metadata line for ang_mean" << endl << line << endl;
                return 1;
            }
            ++index1;
            while(1){
                size_t index2 = index1+1;
                while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
                    ++index2;
                }
                float val = atof(line.substr(index1, index2-index1).c_str());
                ang_means.push_back(val);
                if( line[index2] == ']' )
                    break;
                ++index2;
                index1 = index2;
            }
        }
        
        // look for ang_width
        index1 = line.find("ang_width");
        if(index1 != string::npos){
            ++index1;
            while( index1 < line.length() && line[index1] != '[' ){
                ++index1;
            }
            if( line[index1] != '[' ){
                cerr << "Problem parsing metadata line for ang_width" << endl << line << endl;
                return 1;
            }
            ++index1;
            while(1){
                size_t index2 = index1+1;
                while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
                    ++index2;
                }
                float val = atof(line.substr(index1, index2-index1).c_str());
                ang_widths.push_back(val);
                if( line[index2] == ']' )
                    break;
                ++index2;
                index1 = index2;
            }
        }
        
        getline( metadatafile, line );
    }
    float ang_means_carray[ang_means.size()];
    cerr << "ang_means: ";
    for( unsigned int i=0; i<ang_means.size(); ++i ){
        ang_means_carray[i] = ang_means[i];
        cerr << ang_means[i] << " ";
    }
    cerr << endl;
    float ang_widths_carray[ang_widths.size()];
    cerr << "ang_widths: ";
    for( unsigned int i=0; i<ang_widths.size(); ++i ){
        ang_widths_carray[i] = ang_means[i];
        cerr << ang_widths[i] << " ";
    }
    cerr << endl;
    if( ang_means.size() == 0 || ang_widths.size() == 0 || ang_means.size() != ang_widths.size() ){
        cerr << "Different length ang_means and ang_widths arrays??" << endl;
        return 1;
    }

    float ang_means_c[ang_means.size()];
    float ang_widths_c[ang_means.size()];
    for( unsigned int i=0; i<ang_means.size(); ++i ){
        ang_means_c[i] = ang_means[i];
        ang_widths_c[i] = ang_widths[i];
    }

    char suffix[4] = "ahb";
    make_maps( catalog_filename.c_str(), mask_filename.c_str(), ang_means_c, ang_means.size(), ang_widths_c, ang_means.size(), use_counts, use_mags, suffix );
    
    return 0;
    

}


