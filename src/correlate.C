#include <unistd.h>
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <sstream>
using std::stringstream;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <iomanip>
using namespace std;

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


#define USE_WEIGHTS 0


int main (int argc, char **argv){
    
    cout << setiosflags(ios::fixed);
    
    int min_footprint_order = 1;
    
    string mapname1( argv[1] );
    string mapname2( argv[2] );

    vector<float> r_lows;
    vector<float> r_mids;
    vector<float> r_highs;
    vector<float> r_wids;
    
    // read in hdf5 maps + masks
    
    // file 1
    H5File file1( mapname1, H5F_ACC_RDONLY );
    H5G_obj_t firstType = file1.getObjTypeByIdx(0);
    if( firstType != H5G_GROUP ){
        cerr << "The map is not standard format;  the first element is not a group." << endl;
        return 1;
    }
    string group_name = file1.getObjnameByIdx(0);
    if( group_name != "data" ){
        cerr << "Group name is " << group_name << ", not data" << endl;
        return 1;
    }
    Group group = file1.openGroup("data");
    string dataset1_name = group.getObjnameByIdx(0);
    if( dataset1_name != "map" ){
        cerr << "Dataset 1 name is " << dataset1_name << ", not map" << endl;
        return 1;
    }
    DataSet dataset1 = group.openDataSet(dataset1_name);
    hsize_t ds1size = dataset1.getStorageSize();
    int npix1 = ds1size/2/sizeof(double);
    double *data1 = (double*) malloc( ds1size );
    dataset1.read(data1, PredType::NATIVE_DOUBLE);
    int map1_order = data1[0];
    int map1_ordering = data1[1];
    cerr << "map 1 info " << map1_order << " " << map1_ordering << endl;
    // for( int i=1; i<npix1; ++i ){
    //     for( int j=0; j<2; ++j ){
    //         cout << data1[2*i+j] << " ";
    //     }
    //     cout << endl;
    // }
    
    string dataset2_name = group.getObjnameByIdx(1);
    if( dataset2_name != "mask" ){
        cerr << "Dataset 2 name is " << dataset2_name << ", not mask" << endl;
        return 1;
    }
    DataSet dataset2 = group.openDataSet(dataset2_name);
    hsize_t ds2size = dataset2.getStorageSize();
    int npix2 = ds2size/2/sizeof(double);
    double *data2 = (double*) malloc( ds2size );
    dataset2.read(data2, PredType::NATIVE_DOUBLE);
    int mask1_order = data2[0];
    int mask1_ordering = data2[1];
    if( mask1_ordering != RING ){
        cerr << "Problem, mask1 ordering is not RING but " << mask1_ordering << endl;
        return 1;
    }
    cerr << "mask 1 info " << mask1_order << " " << mask1_ordering << endl;
    vector<int64> mask1_pixels;
    for( int i=1; i<npix2; ++i ){
        for( int j=0; j<2; ++j ){
            if( data2[2*i+1] > 0.5 ){
                mask1_pixels.push_back(data2[2*i]);
            }
        }
    }
    free(data2);
    Healpix_Base2 mask1_base( mask1_order, RING );
    
    // get the metadata from file 1 (why not)
    group = file1.openGroup("meta");
    DataSet dataSet = group.openDataSet("meta");
    Attribute attr = dataSet.openAttribute("u_mean");
    int attrsize = attr.getStorageSize();
    StrType dataType(PredType::C_S1, attrsize);

    char *u_mean = (char*) malloc(attrsize);
    attr.read( dataType, u_mean );
    cerr << "u_mean: " << u_mean << endl;
    string line(u_mean);
    if( line[0] != '[' ){
        cerr << "Problem parsing metadata line for u_mean" << endl << line << endl;
        return 1;
    }
    size_t index1 = 1;
    while(1){
        size_t index2 = index1+1;
        while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
            ++index2;
        }
        float val = atof(line.substr(index1, index2-index1).c_str());
        r_mids.push_back(val*3.1415926/180.);
        if( line[index2] == ']' )
            break;
        ++index2;
        index1 = index2;
    }

    attr = dataSet.openAttribute("u_width");
    attrsize = attr.getStorageSize();
    dataType = StrType(PredType::C_S1, attrsize);
    char *u_width = (char*) malloc(attrsize);
    attr.read( dataType, u_width );
    cerr << "u_width: " << u_width << endl;
    line = string(u_width);
    if( line[0] != '[' ){
        cerr << "Problem parsing metadata line for u_width" << endl << line << endl;
        return 1;
    }
    index1 = 1;
    while(1){
        size_t index2 = index1+1;
        while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
            ++index2;
        }
        float val = atof(line.substr(index1, index2-index1).c_str());
        r_wids.push_back(val*3.1415926/180.);
        if( line[index2] == ']' )
            break;
        ++index2;
        index1 = index2;
    }
    if( r_mids.size() != r_wids.size() ){
        cerr << "Problem parsing metadata: u_mean and u_width are different sizes:" << endl;
        cerr << u_mean << endl;
        cerr << u_width << endl;
        return 1;
    }

    if( r_mids[0] - r_wids[0]/2.0 < 0.0 )
        r_lows.push_back(0.0001*3.1415926/180.);
    else
        r_lows.push_back(r_mids[0] - r_wids[0]/2.0);
    r_highs.push_back(r_mids[0] + r_wids[0]/2.0);
    for( unsigned int i=1; i<r_mids.size(); ++i ){
        r_lows.push_back(r_mids[i] - r_wids[i]/2.0);
        r_highs.push_back(r_mids[i] + r_wids[i]/2.0);
    }
    cerr << "r_lows: ";
    for( unsigned int i=0; i<r_lows.size(); ++i )
        cerr << r_lows[i] << " ";
    cerr << endl;
    cerr << "r_mids: ";
    for( unsigned int i=0; i<r_mids.size(); ++i )
        cerr << r_mids[i] << " ";
    cerr << endl;
    cerr << "r_highs: ";
    for( unsigned int i=0; i<r_highs.size(); ++i )
        cerr << r_highs[i] << " ";
    cerr << endl;
    
    // file 2
    H5File file2( mapname2, H5F_ACC_RDONLY );
    firstType = file2.getObjTypeByIdx(0);
    if( firstType != H5G_GROUP ){
        cerr << "The map is not standard format;  the first element is not a group." << endl;
        return 1;
    }
    group_name = file2.getObjnameByIdx(0);
    if( group_name != "data" ){
        cerr << "Group name is " << group_name << ", not data" << endl;
        return 1;
    }
    group = file2.openGroup("data");
    dataset1_name = group.getObjnameByIdx(0);
    if( dataset1_name != "map" ){
        cerr << "Dataset 1 name is " << dataset1_name << ", not map" << endl;
        return 1;
    }
    dataset1 = group.openDataSet(dataset1_name);
    ds1size = dataset1.getStorageSize();
    int npix3 = ds1size/2/sizeof(double);
    double *data3 = (double*) malloc( ds1size );
    dataset1.read(data3, PredType::NATIVE_DOUBLE);
    int map2_order = data3[0];
    int map2_ordering = data3[1];
    cerr << "map 2 info " << map2_order << " " << map2_ordering << endl;
    
    dataset2_name = group.getObjnameByIdx(1);
    if( dataset2_name != "mask" ){
        cerr << "Dataset 2 name is " << dataset2_name << ", not mask" << endl;
        return 1;
    }
    dataset2 = group.openDataSet(dataset2_name);
    ds2size = dataset2.getStorageSize();
    int npix4 = ds2size/2/sizeof(double);
    double *data4 = (double*) malloc( ds2size );
    dataset2.read(data4, PredType::NATIVE_DOUBLE);
    int mask2_order = data4[0];
    int mask2_ordering = data4[1];
    if( mask1_ordering != RING ){
        cerr << "Problem, mask2 ordering is not RING but " << mask2_ordering << endl;
        return 1;
    }
    cerr << "mask 2 info " << mask2_order << " " << mask2_ordering << endl;
    vector<int64> mask2_pixels;
    for( int i=1; i<npix4; ++i ){
        for( int j=0; j<2; ++j ){
            if( data2[2*i+1] > 0.5 ){
                mask2_pixels.push_back(data4[2*i]);
            }
        }
    }
    free(data4);
    Healpix_Base2 mask2_base( mask2_order, RING );

    if( map1_order != map2_order ){
        cerr << "Map orders must be the same: are " << map1_order << " and " << map2_order << endl;
        return 1;
    }
    int order = map1_order;

    // from the masks, determine the best possible healpix footprint
    // requirements: must cover all of the mask area
    // must cover <80% of the area covered by the next lowest resolution
    // must be order <= 6 ?
    int footprint_order = min_footprint_order;
    double footprint_area = 50000.;
    Healpix_Map<int> footprintMap;
    while(footprint_order <= 8){
        footprintMap = Healpix_Map<int>(footprint_order, RING);
        footprintMap.fill(0);
        for( unsigned int i=0; i<mask1_pixels.size(); ++i ){
            footprintMap[ footprintMap.ang2pix(mask1_base.pix2ang(mask1_pixels[i])) ] = 1;
        }
        for( unsigned int i=0; i<mask2_pixels.size(); ++i ){
            footprintMap[ footprintMap.ang2pix(mask2_base.pix2ang(mask2_pixels[i])) ] = 1;
        }
        double footprint_area_new = 0.;
        for( int i=0; i<footprintMap.Npix(); ++i ){
          if( footprintMap[i] == 1 )
              footprint_area_new += 1.0;
        }
        footprint_area_new *= (41253./footprintMap.Npix());
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
    footprintMap = Healpix_Map<int>(footprint_order, RING);
    footprintMap.fill(0);
    for( unsigned int i=0; i<mask1_pixels.size(); ++i ){
        footprintMap[ footprintMap.ang2pix(mask1_base.pix2ang(mask1_pixels[i])) ] = 1;
    }
    for( unsigned int i=0; i<mask2_pixels.size(); ++i ){
        footprintMap[ footprintMap.ang2pix(mask2_base.pix2ang(mask2_pixels[i])) ] = 1;
    }
    cerr << "Using footprint order " << footprintMap.Order() << " with " << footprint_area << " sq degrees." << endl;
    
    // system( "rm fp.fits" );
    // fitshandle myfits;
    // myfits.create("fp.fits");
    // write_Healpix_map_to_fits(myfits, footprintMap, PLANCK_FLOAT64);
    // myfits.close();

    // fill in the mask Partpix maps
    Partpix_Map2<int> lowzMask(mask1_order, footprintMap);
    lowzMask.fill(0);
    for( unsigned int i=0; i<mask1_pixels.size(); ++i ){
      if( ! footprintMap[ footprintMap.ang2pix(mask1_base.pix2ang(mask1_pixels[i])) ] ){
          cerr << "Skipping low z mask pixel " << i << "!" << endl;
          continue;
      }
      lowzMask[mask1_pixels[i]] = 1;
    }
    mask1_pixels.clear();
    // make full-resolution versions of the masks so we don't keep having to check for mask indices
    Partpix_Map2<int> lowzMatchedMask = Partpix_Map2<int>(order, footprintMap);
    if( order > lowzMask.Order() ){
        lowzMatchedMask.Import_upgrade( lowzMask, footprintMap );
    }
    else if( order == lowzMask.Order() ){
        lowzMatchedMask.Import_nograde( lowzMask, footprintMap );
    }
    else{
        lowzMatchedMask.Import_degrade( lowzMask, footprintMap );
    }
    lowzMask.clear();

    cerr << "Finished making low z Partpix matched mask" << endl;

    Partpix_Map2<int> highzMask(mask2_order, footprintMap);
    highzMask.fill(0);
    for( unsigned int i=0; i<mask2_pixels.size(); ++i ){
      if( ! footprintMap[ footprintMap.ang2pix(mask2_base.pix2ang(mask2_pixels[i])) ] ){
          cerr << "Skipping high z mask pixel " << i << "!" << endl;
          continue;
      }
      highzMask[mask2_pixels[i]] = 1;
    }
    mask2_pixels.clear();
    Partpix_Map2<int> highzMatchedMask = Partpix_Map2<int>(order, footprintMap);
    if( order > highzMask.Order() ){
        highzMatchedMask.Import_upgrade( highzMask, footprintMap );
    }
    else if( order == highzMask.Order() ){
        highzMatchedMask.Import_nograde( highzMask, footprintMap );
    }
    else{
        highzMatchedMask.Import_degrade( highzMask, footprintMap );
    }
    highzMask.clear();

    cerr << "Finished making high z Partpix matched mask" << endl;


    // make the correlation maps (partpix)
    Partpix_Map2<double> lowzMap(order, footprintMap);
    lowzMap.fill(0.0);
    Partpix_Map2<double> highzMap(order, footprintMap);
    highzMap.fill(0.0);
    cerr << "Created Partpix maps for data" << endl;

#if USE_WEIGHTS
    Partpix_Map2<double> lowzWeightMap(order, footprintMap);
    lowzWeightMap.fill(1.0);  
    Partpix_Map2<double> highzWeightMap(order, footprintMap);
    highzWeightMap.fill(1.0);
    cerr << "Created Partpix maps for weights" << endl;
#endif

    for( int i=1; i<npix1; ++i ){
        for( int j=0; j<2; ++j ){
            lowzMap[data1[2*i]] = data1[2*i+1];
        }
    }
    for( int i=1; i<npix3; ++i ){
        for( int j=0; j<2; ++j ){
            highzMap[data3[2*i]] = data3[2*i+1];
        }
    }
    cerr << "Read in the data maps" << endl;

    // Healpix_Map<double> lowzHMap = lowzMap.to_Healpix( 0. );
    // system( "rm lowzMap.fits" );
    // fitshandle myfits;
    // myfits.create("lowzMap.fits");
    // write_Healpix_map_to_fits(myfits, lowzHMap, PLANCK_FLOAT64);
    // myfits.close();
    // Healpix_Map<double> highzHMap = highzMap.to_Healpix( 0. );
    // system( "rm highzMap.fits" );
    // myfits = fitshandle();
    // myfits.create("highzMap.fits");
    // write_Healpix_map_to_fits(myfits, highzHMap, PLANCK_FLOAT64);
    // myfits.close();
    

    // figure out the jackknife regions using the mask pixels
    // we want >=50 regions, bigger than the max radial bin.
    // insist on 90% completeness within each jackknife pixel.
    int jk_order = footprintMap.Order();
    vector<int> jkpixels;
    Partpix_Map<int> jackknifeMap1(jk_order, footprintMap);
    jackknifeMap1.fill(0);
    double pixel_size = sqrt(41253./jackknifeMap1.Npix());
    while( pixel_size > r_highs[r_highs.size()-1] ){
      if( order < jackknifeMap1.Order() ){ // this is unlikely
          Partpix_Map2<int> mask1(jk_order, footprintMap);
          mask1.Import_upgrade(lowzMatchedMask, footprintMap);
          Partpix_Map2<int> mask2(jk_order, footprintMap);
          mask2.Import_upgrade(highzMatchedMask, footprintMap);
          for( int i1=0; i1<jackknifeMap1.Npartpix(); ++i1 ){
              int i = jackknifeMap1.highResPix(i1);
              if( mask1[i] == 1 && mask2[i] == 1 ){
                  jackknifeMap1[i]++;
                  jkpixels.push_back(i);
              }
          }
          if( jkpixels.size() >= 50 )
              break;
          else{
              cerr << "Jackknife order " << jk_order << " has only " << jkpixels.size() << " good pixels." << endl;
              jkpixels.clear();
              ++jk_order;
              jackknifeMap1 = Partpix_Map<int>(jk_order, footprintMap);
              jackknifeMap1.fill(0.);
              pixel_size = sqrt(41253./jackknifeMap1.Npix());
          }
      }
      else if( order == jackknifeMap1.Order() ){ // also unlikely
          for( int i1=0; i1<jackknifeMap1.Npartpix(); ++i1 ){
              int i = jackknifeMap1.highResPix(i1);
              if( lowzMatchedMask[i] == 1 && highzMatchedMask[i] == 1 ){
                  jackknifeMap1[i]++;
                  jkpixels.push_back(i);
              }
          }
          if( jkpixels.size() >= 50 )
              break;
          else{
              cerr << "Jackknife order " << jk_order << " has only " << jkpixels.size() << " good pixels." << endl;
              jkpixels.clear();
              ++jk_order;
              jackknifeMap1 = Partpix_Map<int>(jk_order, footprintMap);
              jackknifeMap1.fill(0);
              pixel_size = sqrt(41253./jackknifeMap1.Npix());
          }
      }
      else{ // data are higher resolution than the jackknife map
          for( int i1=0; i1<lowzMatchedMask.Npartpix(); ++i1 ){
              int64 i = lowzMatchedMask.highResPix(i1);
              if( lowzMatchedMask[i] == 1 ){
                  int64 j = highzMatchedMask.ang2pix(lowzMatchedMask.pix2ang(i));
                  if( highzMatchedMask[j] == 1 ){
                      jackknifeMap1[jackknifeMap1.ang2pix(lowzMatchedMask.pix2ang(i))]++;
                  }
              }
          }
          double threshold = 0;
          for( int j1=0; j1<jackknifeMap1.Npartpix(); ++j1 ){
              int j = jackknifeMap1.highResPix(j1);
              if( jackknifeMap1[j] > threshold )
                  threshold = jackknifeMap1[j];
          }
          threshold /= 2.0;
          //double threshold = 0.9*(highzMatchedMask.Npix()/jackknifeMap1.Npix());
          for( int i1=0; i1<jackknifeMap1.Npartpix(); ++i1 ){
              int i = jackknifeMap1.highResPix(i1);
              if( jackknifeMap1[i] >= threshold ){
                  jkpixels.push_back(i);
              }
          }
          if( jkpixels.size() >= 50 ){
              system( "rm jackknife.fits" );
              fitshandle myfits = fitshandle();
              myfits.create("jackknife.fits");
              Healpix_Map<int> jkmap = jackknifeMap1.to_Healpix( 0 );
              write_Healpix_map_to_fits(myfits, jkmap, PLANCK_FLOAT64);
              myfits.close();
              break;
          }
          else{
              cerr << "Jackknife order " << jk_order << " has only " << jkpixels.size() << " good pixels." << endl;
              jkpixels.clear();
              ++jk_order;
              jackknifeMap1 = Partpix_Map<int>(jk_order, footprintMap);
              jackknifeMap1.fill(0.);
              pixel_size = sqrt(41253./jackknifeMap1.Npix());
          }
      }
    }
    if( pixel_size < r_highs[r_highs.size()-1] ){
      cerr << "Can not find a jackknife resolution with >50 regions and pixel size <" << 180/3.1415926*r_highs[r_highs.size()-1] << " degrees" << endl;
      return 1;
    }
    Healpix_Base jackknifeMap( jk_order, RING );
    cerr << "Constructed a jackknife map with order " << jk_order << ", "  << jkpixels.size() << " good pixels" << endl;


    int64 n_mask=0;
    for( int64 i1=0; i1<lowzMatchedMask.Npartpix(); ++i1 ){
      int64 i = lowzMatchedMask.highResPix(i1);
      if( lowzMatchedMask[i] == 1 && highzMatchedMask[i] == 1 )
          ++n_mask;
    }
    float area_mask = 41252.962*n_mask/lowzMatchedMask.Npix();
    float area_jk = 41252.962*jkpixels.size()/jackknifeMap.Npix();
    float area_fraction = area_mask/area_jk;
    cerr << "Combined mask (" << 3282.8*area_mask << " sq deg) over jackknife (" << 3282.8*area_jk << " sq deg) area = " << area_fraction << endl;



    vector<double> distances(r_lows.size(), 0.);
    vector<double> correlations(r_lows.size(), 0.);
    vector<double> ns(r_lows.size(), 0.);
    vector< vector<double> > pss( r_lows.size(), vector<double>() );

    cerr << "starting correlation stuff, with map resolutions " << lowzMap.Order() << " and " << highzMap.Order() << endl;

    vector<double> jk_means( r_lows.size(), 0. );
    vector<double> jk_variance( r_lows.size(), 0. );
    vector<int> jk_ns( r_lows.size(), 0. );

    // calculate a jkindex_array to save time in the loops
    vector<int> jkindex_array( jackknifeMap.Npix(), -1 );
    for( unsigned int k=0; k<jkpixels.size(); ++k )
        jkindex_array[jkpixels[k]] = k;


    // for distance
    // # pragma omp parallel for schedule(dynamic, 1)
    ofstream corr_file("correlation");
    for( int rbin=0; rbin<(int)r_lows.size(); ++rbin ){

    float r_low = r_lows[rbin];
    float r_high = r_highs[rbin];

    //vector<double> ni( jkpixels.size()+1, 0. );
    //vector<double> sumi( jkpixels.size()+1, 0. );
    vector<double> sumij( jkpixels.size()+1, 0. );
    //vector<double> sumj( jkpixels.size()+1, 0. );
    //vector<double> nj( jkpixels.size()+1, 0. );
    vector<double> nij( jkpixels.size()+1, 0. );

    double sumdist = 0.0;
    double ndist = 0.0;

    // for pixel in high_res map
    for( int64 i1=0; i1<lowzMap.Npartpix(); ++i1 ){
        int64 i = lowzMap.highResPix(i1);

        pointing pointing_i = lowzMap.pix2ang(i);
        if( ! lowzMatchedMask[i] )
            continue;
        int jkpix_index_i = jkindex_array[jackknifeMap.ang2pix( pointing_i )];

        // get pixel indices within r_low + annulus_width
        vector<int64> listpix_outer;
        pointing center = lowzMap.pix2ang(i);
        lowzMap.query_disc( center, r_high, listpix_outer );

        // get pixel indices within r_low
        vector<int64> listpix_inner;
        lowzMap.query_disc( center, r_low, listpix_inner );

        // sort them both
        sort( listpix_inner.begin(), listpix_inner.end() );
        sort( listpix_outer.begin(), listpix_outer.end() );

        // write a new list of ones in the annulus 
        // (faster than deleting from the big list)
        vector<int64> listpix_annulus;

        if( listpix_inner.size() == 0 )
            listpix_annulus = listpix_outer;
        else{
            unsigned int index_inner = 0;
            unsigned int index_outer = 0;
            while( index_inner < listpix_inner.size() && index_outer < listpix_outer.size() ){
              if( listpix_inner[index_inner] != listpix_outer[index_outer] ){
                listpix_annulus.push_back(listpix_outer[index_outer]);
                index_outer++;
              }
              else{
                index_inner++;
                index_outer++;
              }
            }
            while( index_outer < listpix_outer.size() ){
              listpix_annulus.push_back(listpix_outer[index_outer]);
              index_outer++;
            }
        }

      double sinith = sin(1.5707963-pointing_i.theta);
      double cosith = cos(1.5707963-pointing_i.theta);

      for( unsigned int j=0; j<listpix_annulus.size(); ++j ){

        pointing pointing_j = highzMap.pix2ang(listpix_annulus[j]);

        // make sure this pixel is in the high-resolution area!
        if( footprintMap[footprintMap.ang2pix( pointing_j )] == 0 )
            continue;

        if( ! highzMatchedMask[listpix_annulus[j]] )
            continue;

        int jkpix_index_j= jkindex_array[jackknifeMap.ang2pix( pointing_j )];


        double mult3 = lowzMap[i]*highzMap[listpix_annulus[j]];
        double multw = 1.0;

    #if USE_WEIGHTS
        mult3 *= lowzWeightMap[i]*highzWeightMap[listpix_annulus[j]];
        multw *= lowzWeightMap[i]*highzWeightMap[listpix_annulus[j]];
    #endif
        for( unsigned int k=0; k<jkpixels.size()+1; ++k ){
            if( ((int) k == jkpix_index_j) || ((int) k == jkpix_index_i) )
                continue;
            sumij[k] += mult3;
            nij[k] += multw;
        }

        double sinjth = sin(1.5707963-pointing_j.theta);
        double cosjth = cos(1.5707963-pointing_j.theta);
        double sindphi = sin(pointing_i.phi-pointing_j.phi);
        double cosdphi = cos(pointing_i.phi-pointing_j.phi);
        double disty = sqrt( pow(cosith*sindphi, 2.0) + pow( cosjth*sinith - sinjth*cosith*cosdphi, 2.0 ) );
        double distx = sinjth*sinith + cosjth*cosith*cosdphi;
        double dist = atan2(disty,distx);
        if( isnan(dist) )
            cout << "thetas " << pointing_i.theta  << " " << pointing_j.theta << " phis " << pointing_i.phi << " " << pointing_j.phi << " dist " << dist << endl;
        sumdist += dist*multw;
        ndist += multw;

      } // for j pixels in high-z map

    } // for i pixels in low-z map

    double dist_mean = sumdist/ndist;

    if( nij[jkpixels.size()] > 0 ){
        correlations[rbin] = sumij[jkpixels.size()]/nij[jkpixels.size()];
    }
    else{
        correlations[rbin] = 0.0;
    }
    distances[rbin] = dist_mean*180./3.1415926;

    // now do the jackknife
    for( unsigned int k=0; k<jkpixels.size(); ++k ){
        double corr_wo1 = 0.0;
        if( nij[k] > 0.0 ){
            corr_wo1 = sumij[k]/nij[k];
        }
        else{
            cerr << "Warning, empty jackknife pixel being used..." << endl;
        }
      double ps = jkpixels.size()*correlations[rbin] - (jkpixels.size()-1)*corr_wo1;
      pss[rbin].push_back(ps);
    }    

    // now average up the jackknife results for each radius
    unsigned int r = rbin;
    for( unsigned int p=0; p<pss[r].size(); ++p ){
      jk_means[r] += pss[r][p];
    }
    if( pss[r].size()>2 ){
      jk_means[r] /= pss[r].size();

      for( unsigned int p=0; p<pss[r].size(); ++p ){
    jk_variance[r] += (pss[r][p]-jk_means[r])*(pss[r][p]-jk_means[r]);
      }
      float fraction = 1./(pss[r].size()*(pss[r].size()-1.));
      jk_variance[r] *= fraction;
      jk_variance[r] /= area_fraction;
    }
    else{
      jk_means[r] = 0.;
      jk_variance[r] = 0.;
    }

    corr_file << distances[rbin] << " " << correlations[rbin] << " " << jk_means[r] << " " << sqrt(jk_variance[r]) << endl;
    cerr << distances[rbin] << " " << correlations[rbin] << " " << jk_means[r] << " " << sqrt(jk_variance[r]) << endl;

    } // for annulus radius
    corr_file.close();

    // calculate the covariance matrix
    vector< vector<double> > c( r_lows.size(), vector<double>(r_lows.size(), 0.) );
    cerr << "calculating covariance... " << endl;
    for( unsigned int i=0; i<r_lows.size(); ++i ){
      for( unsigned int j=0; j<r_lows.size(); ++j ){
    if( pss[i].size() != pss[j].size() ){
      cerr << "Problem, jackknife lists " << i << ", " << j << " are not the same size: " << pss[i].size() << ", " << pss[j].size() << endl;
      return 1;
    }
    for( unsigned int k=0; k<pss[i].size(); ++k ){
      c[i][j] += (pss[i][k]-jk_means[i])*(pss[j][k]-jk_means[j]);
    }
    float fraction = 1./(pss[i].size()*(pss[i].size()-1.));
    c[i][j] *= fraction;
    c[i][j] /= area_fraction;

      }
    }

    ofstream cov_file("covariance");
    for( unsigned int k=0; k<r_lows.size(); ++k )
      for( unsigned int l=0; l<r_lows.size(); ++l ){
          cov_file << k << " " << l << " " << c[k][l] << endl;
      }
    cov_file.close();

    
    return 0;
    
}