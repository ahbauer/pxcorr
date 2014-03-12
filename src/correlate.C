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

typedef struct jk_struct{
  int64 pixel;
  float value;
} jk_struct;


template<class T>
string string_ize( const T& val ){
  ostringstream s;
  s << val;
  return s.str();
}

void read_hdf5map( string mapname, unsigned long **map_pix, float **map_data, unsigned long &npix, vector<unsigned long>& mask_pixel_list, int& mask_order, vector<float>& r_mids, vector<float>& r_wids, int& zbin, string& ftype, string& pop, int64& nobj ){
    
    
    H5File file1( mapname, H5F_ACC_RDONLY );
    H5G_obj_t firstType = file1.getObjTypeByIdx(0);
    if( firstType != H5G_GROUP ){
        cerr << "The map is not standard format;  the first element is not a group." << endl;
        throw exception();
    }
    string group_name = file1.getObjnameByIdx(0);
    if( group_name != "data" ){
        cerr << "Group name is " << group_name << ", not data" << endl;
        throw exception();
    }
    Group group = file1.openGroup("data");
    DataSet dataset1 = group.openDataSet("pixels");

    hsize_t ds1size = dataset1.getStorageSize();
    npix = ds1size/sizeof(unsigned long);

    (*map_pix) = new unsigned long[ds1size];

    cerr << "reading map_pix... " << endl;
    dataset1.read(*map_pix, PredType::NATIVE_ULONG);
    cerr << " done!" << endl;

    DataSet dataset1a = group.openDataSet("values");

    ds1size = dataset1a.getStorageSize();
    if( npix != ds1size/sizeof(float) ){
        cerr << "Problem, number of pixels " << npix << " != number of data points " << ds1size/sizeof(float) << endl;
        throw exception();
    }

    (*map_data) = new float[ds1size];
    cerr << "reading map_data...";
    dataset1a.read(*map_data, PredType::NATIVE_FLOAT);
    cerr << " done!" << endl;
    
    cerr << "map has order " << (*map_pix)[0] << " and scheme " << (*map_data)[0] << endl;

    DataSet dataset2 = group.openDataSet("mask_pixels");
    hsize_t ds2size = dataset2.getStorageSize();
    unsigned long npix2 = ds2size/sizeof(unsigned long);
    cerr << "mask 1 has " << npix2-1 << " pixels" << endl;

    unsigned long *mask_pix = new unsigned long[ds2size];
    dataset2.read(mask_pix, PredType::NATIVE_ULONG);

    DataSet dataset2a = group.openDataSet("mask_values");
    ds2size = dataset2.getStorageSize();
    if( npix2 != ds2size/sizeof(unsigned long) ){
        cerr << "Problem, mask has " << npix2 << " pixels but " << ds2size/sizeof(unsigned long) << " values" << endl;
        throw exception();
    }

    unsigned long *mask_data = new unsigned long[ds2size];
    dataset2a.read(mask_data, PredType::NATIVE_ULONG);

    if( mask_data[0] != RING ){
        cerr << "Problem, mask ordering is not RING but " << mask_data[0] << endl;
        throw exception();
    }
    mask_order = mask_pix[0];
    cerr << "mask info " << mask_order << " " << mask_data[0] << endl;
    for( unsigned long i=1; i<npix2; ++i ){
        if( mask_data[i] > 0.5 ){
            mask_pixel_list.push_back(mask_pix[i]);
        }
    }
    
    delete[] mask_pix;
    delete[] mask_data;
    
    // get the metadata from the file
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
        throw exception();
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
        throw exception();
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
        throw exception();
    }
    free(u_mean);
    free(u_width);
    
    attr = dataSet.openAttribute("zbin");
    attrsize = attr.getStorageSize();
    //dataType = StrType(PredType::C_S1, attrsize);
    //char *zbin_char = (char*) malloc(attrsize);
    //attr.read( dataType, zbin_char );
    //zbin = string(zbin_char);
    H5::IntType int_type(H5::PredType::NATIVE_INT);
    attr.read( int_type, &zbin );
    cerr << "zbin: " << zbin << endl;
    
    attr = dataSet.openAttribute("nobj");
    attrsize = attr.getStorageSize();
    H5::IntType int64_type(H5::PredType::NATIVE_ULONG);
    attr.read( int64_type, &nobj );
    cerr << "nobj: " << nobj << endl;
    
    attr = dataSet.openAttribute("ftype");
    attrsize = attr.getStorageSize();
    dataType = StrType(PredType::C_S1, attrsize);
    char *ftype_char = (char*) malloc(attrsize);
    attr.read( dataType, ftype_char );
    ftype = string(ftype_char);
    cerr << "ftype: " << ftype << endl;
    
    attr = dataSet.openAttribute("pop");
    attrsize = attr.getStorageSize();
    dataType = StrType(PredType::C_S1, attrsize);
    char *pop_char = (char*) malloc(attrsize);
    attr.read( dataType, pop_char );
    pop = string(pop_char);
    cerr << "pop: " << pop << endl;
    
    cerr << "Finished reading in the map" << endl;
    
    return;
}

void correlate( char* mapn1, char* mapn2, int r, char* outfilename ){
 
    bool OUTPUT_FITS = false;
    
    cout << setiosflags(ios::fixed) << setprecision(16);
    cerr << setprecision(16);
    
    int min_footprint_order = 3;
    double pixel_multiple = 2.5;
    double jk_threshold_frac = 0.5;
    
    string mapname1( mapn1 );
    string mapname2( mapn2 );

    cerr << "Correlating " << mapname1 << " x " << mapname2 << " rbin " << r << endl;
    
    bool single_rbin = true;
    if( r < 0 ){
        single_rbin = false;
    }

    if( access( mapname1.c_str(), F_OK ) == -1 ){
        cerr << "File " << mapname1 << " is not readable!" << endl;
        throw exception();
    }
    if( access( mapname2.c_str(), F_OK ) == -1 ){
        cerr << "File " << mapname2 << " is not readable!" << endl;
        throw exception();
    }
    
    vector<float> r_mids;
    vector<float> r_wids;
    vector<float> r_lows;
    vector<float> r_highs;
    
    // read in hdf5 maps + masks
    
    cerr << "Reading in file 1" << endl;
    
    vector<unsigned long> mask1_pixels;
    unsigned long *map1_pix;
    float *map1_data;
    unsigned long npix1;
    int mask1_order;
    int zbin1_attr;
    string ftype1_attr;
    string pop1_attr;
    int64 nobj1_attr;
    read_hdf5map( mapname1, &map1_pix, &map1_data, npix1, mask1_pixels, mask1_order, r_mids, r_wids, zbin1_attr, ftype1_attr, pop1_attr, nobj1_attr );
    
    cerr << "returned from reading map 1, zbin1_attr = " << zbin1_attr << endl;
    
    int map1_order = int(map1_pix[0]);

    cerr << "map1 order = " << map1_pix[0] << endl;

    int map1_ordering = int(map1_data[0]);
    
    cerr << "map 1 info " << map1_order << " " << map1_ordering << " " << npix1-1 << " pixels" << endl;

    Healpix_Base2 mask1_base( mask1_order, RING );

    cerr << "Reading in file 2" << endl;
    
    vector<unsigned long> mask2_pixels;
    unsigned long *map2_pix;
    float *map2_data;
    unsigned long npix2;
    int mask2_order;
    
    vector<float> r_mids2;
    vector<float> r_wids2;
    int zbin2_attr;
    string ftype2_attr;
    string pop2_attr;
    int64 nobj2_attr;
    read_hdf5map( mapname2, &map2_pix, &map2_data, npix2, mask2_pixels, mask2_order, r_mids2, r_wids2, zbin2_attr, ftype2_attr, pop2_attr, nobj2_attr );

    cerr << "returned from reading map 1, zbin2_attr = " << zbin2_attr << endl;

    int map2_order = map2_pix[0];
    int map2_ordering = map2_data[0];
    
    cerr << "map 2 info " << map2_order << " " << map2_ordering << " " << npix2-1 << " pixels" << endl;

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
    

    Healpix_Base2 mask2_base( mask2_order, RING );

    if( map1_order != map2_order ){
        cerr << "Map orders must be the same: are " << map1_order << " and " << map2_order << endl;
        throw exception();
    }
    int order = map1_order;

    // from the masks, determine the best possible healpix footprint
    // requirements: must cover all of the mask area
    // must cover <80% of the area covered by the next lowest resolution
    // must be order <= 6 ?
    int footprint_order = min_footprint_order;
    double footprint_area = 50000.;
    Healpix_Map<int> *footprintMap;
    while(footprint_order <= mask1_order){
        footprintMap = new Healpix_Map<int>(footprint_order, RING);
        footprintMap->fill(0);
        for( unsigned int i=0; i<mask1_pixels.size(); ++i ){
            (*footprintMap)[ footprintMap->ang2pix(mask1_base.pix2ang(mask1_pixels[i])) ] = 1;
        }
        for( unsigned int i=0; i<mask2_pixels.size(); ++i ){
            (*footprintMap)[ footprintMap->ang2pix(mask2_base.pix2ang(mask2_pixels[i])) ] = 1;
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
        delete footprintMap;
    }
    // the previous footprint order was fine.
    // we want this to be small because it's also the minimum jackknife order.
    // so, go back to the last footprint order.
    --footprint_order;
    delete footprintMap;
    footprintMap = new Healpix_Map<int>(footprint_order, RING);
    footprintMap->fill(0);
    for( unsigned int i=0; i<mask1_pixels.size(); ++i ){
        (*footprintMap)[ footprintMap->ang2pix(mask1_base.pix2ang(mask1_pixels[i])) ] = 1;
    }
    for( unsigned int i=0; i<mask2_pixels.size(); ++i ){
        (*footprintMap)[ footprintMap->ang2pix(mask2_base.pix2ang(mask2_pixels[i])) ] = 1;
    }
    cerr << "Using footprint order " << footprintMap->Order() << " with " << footprint_area << " sq degrees." << endl;

    if( OUTPUT_FITS && footprintMap->Order()<13 ){
        system( "rm fp.fits" );
        fitshandle myfits;
        myfits.create("fp.fits");
        write_Healpix_map_to_fits(myfits, *footprintMap, PLANCK_FLOAT64);
        myfits.close();
    }

    // fill in the mask Partpix maps
    Partpix_Map2<int> *lowzMask = new Partpix_Map2<int>(mask1_order, *footprintMap);
    lowzMask->fill(0);
    for( unsigned int i=0; i<mask1_pixels.size(); ++i ){
      if( ! (*footprintMap)[ footprintMap->ang2pix(mask1_base.pix2ang(mask1_pixels[i])) ] ){
          cerr << "Skipping low z mask pixel " << i << "!" << endl;
          continue;
      }
      (*lowzMask)[mask1_pixels[i]] = 1;
    }

    mask1_pixels.clear();
    // make full-resolution versions of the masks so we don't keep having to check for mask indices
    Partpix_Map2<int> *lowzMatchedMask = new Partpix_Map2<int>(order, *footprintMap);
    if( order > lowzMask->Order() ){
        lowzMatchedMask->Import_upgrade( *lowzMask, *footprintMap );
    }
    else if( order == lowzMask->Order() ){
        lowzMatchedMask->Import_nograde( *lowzMask, *footprintMap );
    }
    else{
        lowzMatchedMask->Import_degrade( *lowzMask, *footprintMap );
    }

    cerr << "Finished making low z Partpix matched mask" << endl;
    cerr << "Imported from resolution " << mask1_order << " to " << order << endl;

    Partpix_Map2<int> *highzMask = new Partpix_Map2<int>(mask2_order, *footprintMap);
    highzMask->fill(0);
    for( unsigned int i=0; i<mask2_pixels.size(); ++i ){
      if( ! (*footprintMap)[ footprintMap->ang2pix(mask2_base.pix2ang(mask2_pixels[i])) ] ){
          cerr << "Skipping high z mask pixel " << i << "!" << endl;
          continue;
      }
      (*highzMask)[mask2_pixels[i]] = 1;
    }
    mask2_pixels.clear();
    Partpix_Map2<int> *highzMatchedMask = new Partpix_Map2<int>(order, *footprintMap);
    // cerr << "high z mask has order " << order << endl;
    if( order > highzMask->Order() ){
        highzMatchedMask->Import_upgrade( *highzMask, *footprintMap );
        // cerr << "import upgraded" << endl;
    }
    else if( order == highzMask->Order() ){
        highzMatchedMask->Import_nograde( *highzMask, *footprintMap );
        // cerr << "import nograded" << endl;
    }
    else{
        highzMatchedMask->Import_degrade( *highzMask, *footprintMap );
        // cerr << "import degraded" << endl;
    }

    if( OUTPUT_FITS ){
        cerr << "outputing fits!" << endl;        
        if( order < 14 ){ 
            system( "rm lowzMask0.fits" );
            fitshandle myfits;
            myfits.create("lowzMask0.fits");
            Healpix_Map<int> lowzHMask = lowzMask->to_Healpix( 0 );
            write_Healpix_map_to_fits(myfits, lowzHMask, PLANCK_INT32);
            myfits.close();

            system( "rm lowzMask.fits" );
            myfits = fitshandle();
            myfits.create("lowzMask.fits");
            Healpix_Map<int> lowzHMask0 = lowzMatchedMask->to_Healpix( 0 );
            write_Healpix_map_to_fits(myfits, lowzHMask0, PLANCK_INT32);
            myfits.close();

            system( "rm highzMask.fits" );
            myfits = fitshandle();
            myfits.create("highzMask.fits");
            Healpix_Map<int> highzHMask = highzMatchedMask->to_Healpix( 0 );
            write_Healpix_map_to_fits(myfits, highzHMask, PLANCK_INT32);
            myfits.close();
        }
        else{
            system( "rm lowzMask.fits" );
            fitshandle myfits;
            myfits.create("lowzMask.fits");
            Partpix_Map2<int> lowzHMask1( 10, *footprintMap );
            lowzHMask1.Import_degrade( *lowzMatchedMask, *footprintMap );
            Healpix_Map<int> lowzHMask0 = lowzHMask1.to_Healpix( 0 );
            write_Healpix_map_to_fits(myfits, lowzHMask0, PLANCK_INT32);
            myfits.close();

            system( "rm highzMask.fits" );
            myfits = fitshandle();
            myfits.create("highzMask.fits");
            Partpix_Map2<int> highzHMask1( 10, *footprintMap );
            highzHMask1.Import_degrade( *highzMatchedMask, *footprintMap );
            Healpix_Map<int> highzHMask = highzHMask1.to_Healpix( 0 );
            write_Healpix_map_to_fits(myfits, highzHMask, PLANCK_INT32);
            myfits.close();
        }
    }

    delete lowzMask;
    delete highzMask;


    cerr << "Finished making high z Partpix matched mask" << endl;
    cerr << "Imported from resolution " << mask2_order << " to " << order << endl;

    // make the correlation maps (partpix)
    Partpix_Map2<float> *lowzMap = new Partpix_Map2<float>(order, *footprintMap);
    lowzMap->fill(0.0);
    Partpix_Map2<float> *highzMap = new Partpix_Map2<float>(order, *footprintMap);
    highzMap->fill(0.0);
    cerr << "Created Partpix maps for data, order " << order << endl;

#if USE_WEIGHTS
    Partpix_Map2<float> *lowzWeightMap = new Partpix_Map2<float>(order, *footprintMap);
    lowzWeightMap->fill(1.0);  
    Partpix_Map2<float> *highzWeightMap = new Partpix_Map2<float>(order, *footprintMap);
    highzWeightMap->fill(1.0);
    cerr << "Created Partpix maps for weights" << endl;
#endif

    for( unsigned long i=1; i<npix1; ++i ){
        // cerr << "trying low z pixel " << map_data[2*i] << " value " << map_data[2*i+1] << " i = " << i << endl;
        (*lowzMap)[map1_pix[i]] = map1_data[i];
    }
    for( unsigned long i=1; i<npix2; ++i ){
        // cerr << "trying high z pixel " << data3[2*i] << " value " << data3[2*i+1] << " i = " << i << endl;
        (*highzMap)[map2_pix[i]] = map2_data[i];
    }
    cerr << "Read in the data maps" << endl;
    free(map1_pix);
    free(map1_data);
    free(map2_pix);
    free(map2_data);


    if( OUTPUT_FITS && order < 13 ){
        Healpix_Map<float> lowzHMap = lowzMap->to_Healpix( 0. );
        system( "rm lowzMap.fits" );
        fitshandle myfits;
        myfits.create("lowzMap.fits");
        write_Healpix_map_to_fits(myfits, lowzHMap, PLANCK_FLOAT64);
        myfits.close();
        Healpix_Map<float> highzHMap = highzMap->to_Healpix( 0. );
        system( "rm highzMap.fits" );
        myfits = fitshandle();
        myfits.create("highzMap.fits");
        write_Healpix_map_to_fits(myfits, highzHMap, PLANCK_FLOAT64);
        myfits.close();
    }

    // figure out the jackknife regions using the mask pixels
    // we want >=50 regions, bigger than the max radial bin.
    // insist on 90% completeness within each jackknife pixel.
    int jk_order = footprintMap->Order();
    vector<int> jkpixels;
    Partpix_Map2<int> *jackknifeMap1 = new Partpix_Map2<int>(jk_order, *footprintMap);
    jackknifeMap1->fill(0);
    double pixel_size = sqrt(41253./jackknifeMap1->Npix());
    while( pixel_size > r_highs[r_highs.size()-1] ){
      if( order < jackknifeMap1->Order() ){ // this is unlikely
          Partpix_Map2<int> *mask1 = new Partpix_Map2<int>(jk_order, *footprintMap);
          mask1->Import_upgrade(*lowzMatchedMask, *footprintMap);
          Partpix_Map2<int> *mask2 = new Partpix_Map2<int>(jk_order, *footprintMap);
          mask2->Import_upgrade(*highzMatchedMask, *footprintMap);
          for( int i1=0; i1<jackknifeMap1->Npartpix(); ++i1 ){
              int i = jackknifeMap1->highResPix(i1);
              if( (*mask1)[i] == 1 && (*mask2)[i] == 1 ){
                  (*jackknifeMap1)[i]++;
                  jkpixels.push_back(i);
              }
          }
          if( jkpixels.size() >= 50 )
              break;
          else{
              cerr << "Jackknife order " << jk_order << " has only " << jkpixels.size() << " good pixels." << endl;
              jkpixels.clear();
              ++jk_order;
              delete jackknifeMap1;
              jackknifeMap1 = new Partpix_Map2<int>(jk_order, *footprintMap);
              jackknifeMap1->fill(0.);
              pixel_size = sqrt(41253./jackknifeMap1->Npix());
          }
          delete mask1;
          delete mask2;
      }
      else if( order == jackknifeMap1->Order() ){ // also unlikely
          for( int i1=0; i1<jackknifeMap1->Npartpix(); ++i1 ){
              int i = jackknifeMap1->highResPix(i1);
              if( (*lowzMatchedMask)[i] == 1 && (*highzMatchedMask)[i] == 1 ){
                  (*jackknifeMap1)[i]++;
                  jkpixels.push_back(i);
              }
          }
          if( jkpixels.size() >= 50 )
              break;
          else{
              cerr << "Jackknife order " << jk_order << " has only " << jkpixels.size() << " good pixels." << endl;
              jkpixels.clear();
              ++jk_order;
              delete jackknifeMap1;
              jackknifeMap1 = new Partpix_Map2<int>(jk_order, *footprintMap);
              jackknifeMap1->fill(0);
              pixel_size = sqrt(41253./jackknifeMap1->Npix());
          }
      }
      else{ // data are higher resolution than the jackknife map
          for( int i1=0; i1<lowzMatchedMask->Npartpix(); ++i1 ){
              int64 i = lowzMatchedMask->highResPix(i1);
              if( (*lowzMatchedMask)[i] == 1 ){
                  int64 j = highzMatchedMask->ang2pix(lowzMatchedMask->pix2ang(i));
                  if( (*highzMatchedMask)[j] == 1 ){
                      (*jackknifeMap1)[jackknifeMap1->ang2pix(lowzMatchedMask->pix2ang(i))]++;
                  }
              }
          }
          double threshold = 0.;
          for( int j1=0; j1<jackknifeMap1->Npartpix(); ++j1 ){
              int j = jackknifeMap1->highResPix(j1);
              if( (*jackknifeMap1)[j] > threshold )
                  threshold = (*jackknifeMap1)[j];
          }
          threshold *= jk_threshold_frac;
          for( int i1=0; i1<jackknifeMap1->Npartpix(); ++i1 ){
              int i = jackknifeMap1->highResPix(i1);
              if( (*jackknifeMap1)[i] >= threshold ){
                  jkpixels.push_back(i);
              }
          }
          if( jkpixels.size() >= 50 ){
              if( OUTPUT_FITS ){
                  system( "rm jackknife.fits" );
                  fitshandle myfits = fitshandle();
                  myfits.create("jackknife.fits");
                  Healpix_Map<int> jkmap = jackknifeMap1->to_Healpix( 0 );
                  write_Healpix_map_to_fits(myfits, jkmap, PLANCK_FLOAT64);
                  myfits.close();
              }
              break;
          }
          else{
              cerr << "Jackknife order " << jk_order << " has only " << jkpixels.size() << " good pixels." << endl;
              jkpixels.clear();
              ++jk_order;
              delete jackknifeMap1;
              jackknifeMap1 = new Partpix_Map2<int>(jk_order, *footprintMap);
              jackknifeMap1->fill(0);
              pixel_size = sqrt(41253./jackknifeMap1->Npix());
          }
      }
    }
    delete jackknifeMap1;
    
    if( pixel_size < r_highs[r_highs.size()-1] ){
      cerr << "Can not find a jackknife resolution with >50 regions and pixel size <" << 180/3.1415926*r_highs[r_highs.size()-1] << " degrees" << endl;
      throw exception();
    }
    Healpix_Base jackknifeMap( jk_order, RING );
    cerr << "Constructed a jackknife map with order " << jk_order << ", "  << jkpixels.size() << " good pixels" << endl;


    int64 n_mask=0;
    for( int64 i1=0; i1<lowzMatchedMask->Npartpix(); ++i1 ){
      int64 i = lowzMatchedMask->highResPix(i1);
      if( (*lowzMatchedMask)[i] && (*highzMatchedMask)[i] )
          ++n_mask;
    }
    float area_mask = 41252.962*n_mask/lowzMatchedMask->Npix();
    float area_jk = 41252.962*jkpixels.size()/jackknifeMap.Npix();
    float area_fraction = area_mask/area_jk;
    cerr << "Combined mask (" << area_mask << " sq deg) over jackknife (" << area_jk << " sq deg) area = " << area_fraction << " from " << n_mask << " mask pixels" << endl;



    vector<double> distances(r_lows.size(), 0.);
    vector<double> correlations(r_lows.size(), 0.);
    vector<double> ns(r_lows.size(), 0.);
    vector< vector<double> > pss( r_lows.size(), vector<double>() );
    vector< vector<double> > jkcorrs( r_lows.size(), vector<double>() );

    cerr << "starting correlation stuff, with map resolutions " << lowzMap->Order() << " and " << highzMap->Order() << endl;

    vector<double> jk_means( r_lows.size(), 0. );
    vector<double> jk_variance( r_lows.size(), 0. );
    vector<int> jk_ns( r_lows.size(), 0. );

    // calculate a jkindex_array to save time in the loops
    vector<int> jkindex_array( jackknifeMap.Npix(), jkpixels.size() );
    for( unsigned int k=0; k<jkpixels.size(); ++k )
        jkindex_array[jkpixels[k]] = k;


    // make a copy of the original maps so that we can use degraded versions in the correlations
    Partpix_Map2<float> *lowzMap_orig = lowzMap;
    Partpix_Map2<float> *highzMap_orig = highzMap;
    Partpix_Map2<int> *lowzMatchedMask_orig = lowzMatchedMask;
    Partpix_Map2<int> *highzMatchedMask_orig = highzMatchedMask;
    // lowzMap_orig = Partpix_Map2<double>(order, *footprintMap);
    // highzMap_orig = Partpix_Map2<double>(order, *footprintMap);
    // lowzMatchedMask_orig = Partpix_Map2<int>(order, *footprintMap);
    // highzMatchedMask_orig = Partpix_Map2<int>(order, *footprintMap);
    // lowzMap_orig.Import_nograde(lowzMap, *footprintMap);
    // highzMap_orig.Import_nograde(highzMap, *footprintMap);
    // lowzMatchedMask_orig.Import_nograde(lowzMatchedMask, *footprintMap);
    // highzMatchedMask_orig.Import_nograde(highzMatchedMask, *footprintMap);

    vector< vector<double> > sumij_injk;
    vector< vector<double> > nij_injk;
    
    // for distance
    // # pragma omp parallel for schedule(dynamic, 1)
    for( int rbin=0; rbin<(int)r_lows.size(); ++rbin ){

        if( single_rbin ){
            rbin = r;
        }

        float r_low = r_lows[rbin];
        float r_high = r_highs[rbin];

        // make sure that the data resolution is necessary.
        // this is assuming that the bins are never decreasing in width with increasing radius!!
        cerr << "Changing the map orders from " << lowzMap_orig->Order() << " to ";
        order = 1;
        while(1){
          float pixel_size = sqrt(41253./(12.0*pow(pow(2.0, order), 2.0)));
          if( pixel_multiple*pixel_size < 57.3*(r_high-r_low) )  // this 0.85 is a bit arbitrary
              break;
          ++order;
        }
        cerr << "order " << order << "... ";
        if( order < lowzMap_orig->Order() && order > jackknifeMap.Order() ){
            lowzMap = new Partpix_Map2<float>(order, *footprintMap);
            highzMap = new Partpix_Map2<float>(order, *footprintMap);
            lowzMatchedMask = new Partpix_Map2<int>(order, *footprintMap);
            highzMatchedMask = new Partpix_Map2<int>(order, *footprintMap);
            lowzMap->Import_degrade(*lowzMap_orig, *footprintMap);
            highzMap->Import_degrade(*highzMap_orig, *footprintMap);
            lowzMatchedMask->Import_degrade(*lowzMatchedMask_orig, *footprintMap);
            highzMatchedMask->Import_degrade(*highzMatchedMask_orig, *footprintMap);
        }
        else if( order < lowzMap_orig->Order() && order < jackknifeMap.Order()+1 && lowzMap_orig->Order() > jackknifeMap.Order()+1 ){
            lowzMap = new Partpix_Map2<float>(jackknifeMap.Order()+1, *footprintMap);
            highzMap = new Partpix_Map2<float>(jackknifeMap.Order()+1, *footprintMap);
            lowzMatchedMask = new Partpix_Map2<int>(jackknifeMap.Order()+1, *footprintMap);
            highzMatchedMask = new Partpix_Map2<int>(jackknifeMap.Order()+1, *footprintMap);
            lowzMap->Import_degrade(*lowzMap_orig, *footprintMap);
            highzMap->Import_degrade(*highzMap_orig, *footprintMap);
            lowzMatchedMask->Import_degrade(*lowzMatchedMask_orig, *footprintMap);
            highzMatchedMask->Import_degrade(*highzMatchedMask_orig, *footprintMap);
        }
        else if( order == lowzMap_orig->Order() ){
            // lowzMap = Partpix_Map2<double>(order, *footprintMap);
            // highzMap = Partpix_Map2<double>(order, *footprintMap);
            // lowzMatchedMask = Partpix_Map2<int>(order, *footprintMap);
            // highzMatchedMask = Partpix_Map2<int>(order, *footprintMap);
            // lowzMap.Import_nograde(lowzMap_orig, *footprintMap);
            // highzMap.Import_nograde(highzMap_orig, *footprintMap);
            // lowzMatchedMask.Import_nograde(lowzMatchedMask_orig, *footprintMap);
            // highzMatchedMask.Import_nograde(highzMatchedMask_orig, *footprintMap);
            lowzMap = lowzMap_orig;
            highzMap = highzMap_orig;
            lowzMatchedMask = lowzMatchedMask_orig;
            highzMatchedMask = highzMatchedMask_orig;
        }
        else{
            cerr << "This combination of orders is not supported!  lowzMap.Order() = " << lowzMap->Order() << ", order = " << order << ", jackknifeMap.Order() = " << jackknifeMap.Order() << endl;
            throw exception();
        }
        order = lowzMap->Order();
        cerr << "done!" << endl;

        //vector<double> ni( jkpixels.size()+1, 0. );
        //vector<double> sumi( jkpixels.size()+1, 0. );
        vector<double> sumij( jkpixels.size()+1, 0. );
        sumij_injk = vector< vector<double> >( jkpixels.size()+1, vector<double>( jkpixels.size()+1, 0. ) );
        //vector<double> sumj( jkpixels.size()+1, 0. );
        //vector<double> nj( jkpixels.size()+1, 0. );
        vector<double> nij( jkpixels.size()+1, 0. );
        nij_injk = vector< vector<double> >( jkpixels.size()+1, vector<double>( jkpixels.size()+1, 0. ) );

        //double sumdist = 0.0;
        //double ndist = 0.0;

        // for pixel in high_res map
        for( int64 i1=0; i1<lowzMap->Npartpix(); ++i1 ){
            int64 i = lowzMap->highResPix(i1);

            pointing pointing_i = lowzMap->pix2ang(i);
            if( ! (*lowzMatchedMask)[i] )
                continue;
            int jkpix_index_i = jkindex_array[jackknifeMap.ang2pix( pointing_i )];

            // get pixel indices within r_low + annulus_width
            vector<int64> listpix_outer;
            pointing center = lowzMap->pix2ang(i);
            lowzMap->query_disc( center, r_high, listpix_outer );

            // get pixel indices within r_low
            vector<int64> listpix_inner;
            lowzMap->query_disc( center, r_low, listpix_inner );

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

          //double sinith = sin(1.5707963-pointing_i.theta);
          //double cosith = cos(1.5707963-pointing_i.theta);

          for( unsigned int j=0; j<listpix_annulus.size(); ++j ){

            pointing pointing_j = highzMap->pix2ang(listpix_annulus[j]);

            // make sure this pixel is in the high-resolution area!
            if( (*footprintMap)[footprintMap->ang2pix( pointing_j )] == 0 )
                continue;

            if( ! (*highzMatchedMask)[listpix_annulus[j]] )
                continue;

            int jkpix_index_j= jkindex_array[jackknifeMap.ang2pix( pointing_j )];

            if( i == listpix_annulus[j] ){
                continue;
            }

            double mult3 = (*lowzMap)[i]*(*highzMap)[listpix_annulus[j]];
            double multw = 1.0;

    #if USE_WEIGHTS
            mult3 *= (*lowzWeightMap)[i]*(*highzWeightMap)[listpix_annulus[j]];
            multw *= (*lowzWeightMap)[i]*(*highzWeightMap)[listpix_annulus[j]];
    #endif

            sumij_injk[jkpix_index_i][jkpix_index_j] += mult3;
            nij_injk[jkpix_index_i][jkpix_index_j] += multw;

            // double sinjth = sin(1.5707963-pointing_j.theta);
            // double cosjth = cos(1.5707963-pointing_j.theta);
            // double sindphi = sin(pointing_i.phi-pointing_j.phi);
            // double cosdphi = cos(pointing_i.phi-pointing_j.phi);
            // double disty = sqrt( pow(cosith*sindphi, 2.0) + pow( cosjth*sinith - sinjth*cosith*cosdphi, 2.0 ) );
            // double distx = sinjth*sinith + cosjth*cosith*cosdphi;
            // double dist = atan2(disty,distx);
            // if( isnan(dist) )
            //     cout << "thetas " << pointing_i.theta  << " " << pointing_j.theta << " phis " << pointing_i.phi << " " << pointing_j.phi << " dist " << dist << endl;
            // sumdist += dist*multw;
            // ndist += multw;

          } // for j pixels in high-z map

        } // for i pixels in low-z map

        //double dist_mean = sumdist/ndist;
        
        // now rearrange jackknife stuff to exclude the pixel
        for( unsigned int i=0; i<jkpixels.size()+1; ++i ){
            for( unsigned int j=0; j<jkpixels.size()+1; ++j ){
                for( unsigned int k=0; k<jkpixels.size(); ++k ){
                    if( k != i && k != j ){
                        sumij[k] += sumij_injk[i][j];
                        nij[k] += nij_injk[i][j];
                    }
                    sumij[jkpixels.size()] += sumij_injk[i][j];
                    nij[jkpixels.size()] += nij_injk[i][j];
                }
            }
        }

        
        
        if( nij[jkpixels.size()] > 0 ){
            correlations[rbin] = sumij[jkpixels.size()]/nij[jkpixels.size()];
        }
        else{
            correlations[rbin] = 0.0;
        }
        distances[rbin] = r_mids[rbin]*180./3.1415926;  //dist_mean*180./3.1415926;

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
          jkcorrs[rbin].push_back(corr_wo1);
        }    

        // now average up the jackknife results for each radius
        for( unsigned int p=0; p<pss[rbin].size(); ++p ){
          //jk_means[rbin] += pss[rbin][p];
            jk_means[rbin] += jkcorrs[rbin][p];
        }
        if( pss[rbin].size()>2 ){
          //jk_means[rbin] /= pss[rbin].size();
          jk_means[rbin] /= jkcorrs[rbin].size();

          for( unsigned int p=0; p<pss[rbin].size(); ++p ){
            //jk_variance[rbin] += (pss[rbin][p]-jk_means[rbin])*(pss[rbin][p]-jk_means[rbin]);
              jk_variance[rbin] += (jkcorrs[rbin][p]-jk_means[rbin])*(jkcorrs[rbin][p]-jk_means[rbin]);
          }
          //float fraction = 1./(pss[rbin].size()*(pss[rbin].size()-1.));
          float fraction = (jkcorrs[rbin].size()-1.)/jkcorrs[rbin].size();
          jk_variance[rbin] *= fraction;
          jk_variance[rbin] /= area_fraction;
        }
        else{
          jk_means[rbin] = 0.;
          jk_variance[rbin] = 0.;
        }

        cerr << distances[rbin] << " " << correlations[rbin] << " " << jk_means[rbin] << " " << sqrt(jk_variance[rbin]) << endl;

        if( lowzMatchedMask != lowzMatchedMask_orig )
            delete lowzMatchedMask;
        if( highzMatchedMask != highzMatchedMask_orig )
            delete highzMatchedMask;
        if( lowzMap != lowzMap_orig )
            delete lowzMap;
        if( highzMap != highzMap_orig )
            delete highzMap;

        if( single_rbin ){
            rbin = r_lows.size();
        }

    } // for annulus radius
    
    delete lowzMatchedMask_orig;
    delete highzMatchedMask_orig;
    delete lowzMap_orig;
    delete highzMap_orig;
    delete footprintMap;

#if USE_WEIGHTS
    delete lowzWeightMap;
    delete highzWeightMap;
#endif

    if( ! single_rbin ){
        
        stringstream corr_line;
        stringstream jk_line;
        stringstream err_line;
        corr_line << "[" << correlations[0];
        jk_line << "[" << jk_means[0];
        err_line << "[" << sqrt(jk_variance[0]);
        for( unsigned int i=1; i<correlations.size(); ++i ){
            corr_line << ", " << correlations[i];
            jk_line << ", " << jk_means[i];
            err_line << ", " << sqrt(jk_variance[i]);
        }
        corr_line << "]";
        jk_line << "]";
        err_line << "]";
        string corr_filename( "corr" );
        ofstream corr_file(corr_filename.c_str());
        corr_file << "[" << corr_line.str() << ", " << jk_line.str() << ", " << err_line.str() << "]" << endl;
        corr_file.close();

        // calculate the covariance matrix
        vector< vector<double> > c( r_lows.size(), vector<double>(r_lows.size(), 0.) );
        cerr << "calculating covariance... " << endl;
        for( unsigned int i=0; i<r_lows.size(); ++i ){
          for( unsigned int j=0; j<r_lows.size(); ++j ){
            if( jkcorrs[i].size() != jkcorrs[j].size() ){
                cerr << "Problem, jackknife lists " << i << ", " << j << " are not the same size: " << jkcorrs[i].size() << ", " << jkcorrs[j].size() << endl;
                throw exception();
            }
            for( unsigned int k=0; k<pss[i].size(); ++k ){
              //c[i][j] += (pss[i][k]-jk_means[i])*(pss[j][k]-jk_means[j]);
              c[i][j] += (jkcorrs[i][k]-jk_means[i])*(jkcorrs[j][k]-jk_means[j]);
            }
            //float fraction = 1./(pss[i].size()*(pss[i].size()-1.));
            //float fraction = 1./(pss[i].size()-1.);
            float fraction = (jkcorrs[i].size()-1.)/(jkcorrs[i].size());
            c[i][j] *= fraction;
            c[i][j] /= area_fraction;

          }
        }

        string cov_filename = "cov";
        ofstream cov_file(cov_filename.c_str());
        cov_file << "[";
        for( unsigned int k=0; k<r_lows.size(); ++k ){
            if( k == 0 ){
                cov_file << "[" << c[k][0];
            }
            else{
                cov_file << ", [" << c[k][0];
            }
            for( unsigned int l=1; l<r_lows.size(); ++l ){
                cov_file << ", " << c[k][l];
            }
            cov_file << "]";
        }
        cov_file << "]" << endl;
        cov_file.close();
        

    } // if not single_rbin
    
    // write out an hdf file
    else{
        
        H5::H5File *file = new H5::H5File( outfilename, H5F_ACC_TRUNC ); // clobber!
        
        // save the correlation function result
        H5::Group* corr_group = new H5::Group( file->createGroup( "/corr" ));
        H5::PredType datatype_results( H5::PredType::NATIVE_FLOAT );
        hsize_t dimsf[1];
        dimsf[0] = 1;
        int rank = 1;
        float *data = new float[1];
        data[0] = correlations[r];
        H5::DataSpace *dataspace1 = new H5::DataSpace( rank, dimsf );
        string datasetname1 = "/corr/corr0";
        H5::DataSet *dataset1 = new H5::DataSet( file->createDataSet( datasetname1, datatype_results, *dataspace1 ) );
        dataset1->write( data, H5::PredType::NATIVE_FLOAT );
        delete[](data);
        delete dataspace1;
        delete dataset1;
        
        // save the jackknife array in a table
        H5::Group* jk_group = new H5::Group( file->createGroup( "/JK" ));
        jk_struct *jk_data = new jk_struct[jkpixels.size()];
        for( unsigned int k=0; k<jkpixels.size();++k ){
            jk_data[k].pixel = jkpixels[k];
            jk_data[k].value = jkcorrs[r][k];
        }
        
        CompType mtype1( sizeof(jk_struct) );
        const H5std_string MEMBER1( "pixel" );
        const H5std_string MEMBER2( "value" );
        mtype1.insertMember( MEMBER1, HOFFSET(jk_struct, pixel), H5::PredType::NATIVE_ULONG);
        mtype1.insertMember( MEMBER2, HOFFSET(jk_struct, value), H5::PredType::NATIVE_FLOAT);
        
        dimsf[0] = jkpixels.size();
        H5::DataSpace *dataspace2 = new H5::DataSpace( rank, dimsf );
        string datasetname2 = "/JK/JK0";
        H5::DataSet *dataset2 = new H5::DataSet( file->createDataSet( datasetname2, mtype1, *dataspace2 ) );
        dataset2->write( jk_data, mtype1 );
        
        // Create new dataspace for attribute
        H5::DataSpace attr1_dataspace = H5::DataSpace( H5S_SCALAR );
        // Create new datatype for attribute
        // H5::PredType datatype_int( H5::PredType::NATIVE_INT );
        H5::IntType int_type(H5::PredType::NATIVE_INT);
        // Set up write buffer for attribute
        // Create attribute and write to it
        H5::Attribute order_attr = dataset2->createAttribute("order", int_type, attr1_dataspace);
        order_attr.write(int_type, &jk_order);
        
        // Create new dataspace for attribute
        H5::DataSpace attr2_dataspace = H5::DataSpace( H5S_SCALAR );
        // Set up write buffer for attribute
        int scheme = 0;
        // Create attribute and write to it
        H5::Attribute scheme_attr = dataset2->createAttribute("scheme", int_type, attr2_dataspace);
        scheme_attr.write(int_type, &scheme);
        
        delete[](jk_data);
        delete dataspace2;
        delete dataset2;
        
        // save metadata
        H5::Group *metagroup = new H5::Group( file->createGroup( "/meta" ) );
        H5::PredType metadatatype( H5::PredType::NATIVE_CHAR );
        dimsf[0] = 1;
        H5::DataSpace *metadataspace = new H5::DataSpace( rank, dimsf );
        H5::DataSet *metadataset = new H5::DataSet( file->createDataSet( "/meta/meta", metadatatype, *metadataspace ) );
        
        // Create new dataspace for attribute
        H5::DataSpace attr3_dataspace = H5::DataSpace( H5S_SCALAR );
        // Create new string datatype for attribute
        H5::StrType strdatatype = H5::StrType( H5::PredType::C_S1, 6 ); // of length 256 characters
        // Set up write buffer for attribute
        const H5std_string fourier_buf( "false" );
        // Create attribute and write to it
        H5::Attribute fourier_attr = metadataset->createAttribute("fourier", strdatatype, attr3_dataspace);
        fourier_attr.write(strdatatype, fourier_buf);
        
        
        // Create new dataspace for attribute
        H5::DataSpace attr4_dataspace = H5::DataSpace( H5S_SCALAR );
        // Create new string datatype for attribute
        strdatatype = H5::StrType( H5::PredType::C_S1, ftype1_attr.size()+1 ); // of length 256 characters    
        // Set up write buffer for attribute
        // string ftype_string = "[\"" + ftype1 + "\"]";
        const H5std_string ftype0_buf( ftype1_attr );
        // Create attribute and write to it
        H5::Attribute ftype0_attr = metadataset->createAttribute("ftype0", strdatatype, attr4_dataspace);
        ftype0_attr.write(strdatatype, ftype0_buf);
        
        // Create new dataspace for attribute
        attr4_dataspace = H5::DataSpace( H5S_SCALAR );
        // Create new string datatype for attribute
        strdatatype = H5::StrType( H5::PredType::C_S1, ftype2_attr.size()+1 ); // of length 256 characters    
        // Set up write buffer for attribute
        // ftype_string = "[\"" + ftype2 + "\"]";
        const H5std_string ftype1_buf( ftype2_attr );
        // Create attribute and write to it
        H5::Attribute ftype1_attr = metadataset->createAttribute("ftype1", strdatatype, attr4_dataspace);
        ftype1_attr.write(strdatatype, ftype1_buf);

        // Create new dataspace for attribute
        attr4_dataspace = H5::DataSpace( H5S_SCALAR );
        // Create new string datatype for attribute
        strdatatype = H5::StrType( H5::PredType::C_S1, pop1_attr.size()+1 ); // of length 256 characters    
        // Set up write buffer for attribute
        // string pop_string = "[\"" + string_ize(pop1) + "\"]";
        const H5std_string pop0_buf( pop1_attr );
        // Create attribute and write to it
        H5::Attribute pop_attr = metadataset->createAttribute("pop0", strdatatype, attr4_dataspace);
        pop_attr.write(strdatatype, pop0_buf);
        
        // Create new dataspace for attribute
        attr4_dataspace = H5::DataSpace( H5S_SCALAR );
        // Create new string datatype for attribute
        strdatatype = H5::StrType( H5::PredType::C_S1, pop2_attr.size()+1 ); // of length 256 characters    
        // Set up write buffer for attribute
        // pop_string = "[\"" + string_ize(pop2) + "\"]";
        const H5std_string pop1_buf( pop2_attr );
        // Create attribute and write to it
        pop_attr = metadataset->createAttribute("pop1", strdatatype, attr4_dataspace);
        pop_attr.write(strdatatype, pop1_buf);

        // // Create new dataspace for attribute
        // attr4_dataspace = H5::DataSpace( H5S_SCALAR );
        // // Create new string datatype for attribute
        // strdatatype = H5::StrType( H5::PredType::C_S1, 256 ); // of length 256 characters    
        // // Set up write buffer for attribute
        // const H5std_string zbin0_buf( zbin1_attr );
        // // Create attribute and write to it
        // H5::Attribute zbin_attr = metadataset->createAttribute("zbin0", strdatatype, attr4_dataspace);
        // zbin_attr.write(strdatatype, zbin0_buf);
        // 
        // // Create new dataspace for attribute
        // attr4_dataspace = H5::DataSpace( H5S_SCALAR );
        // // Create new string datatype for attribute
        // strdatatype = H5::StrType( H5::PredType::C_S1, 256 ); // of length 256 characters    
        // // Set up write buffer for attribute
        // const H5std_string zbin1_buf( zbin2_attr );
        // // Create attribute and write to it
        // zbin_attr = metadataset->createAttribute("zbin1", strdatatype, attr4_dataspace);
        // zbin_attr.write(strdatatype, zbin1_buf);
        
        
        // Create new dataspace for attribute
        attr4_dataspace = H5::DataSpace( H5S_SCALAR );
        H5::Attribute zbin_attr = metadataset->createAttribute("zbin0", int_type, attr4_dataspace);
        zbin_attr.write(int_type, &zbin1_attr);
        
        // Create new dataspace for attribute
        attr4_dataspace = H5::DataSpace( H5S_SCALAR );
        zbin_attr = metadataset->createAttribute("zbin1", int_type, attr4_dataspace);
        zbin_attr.write(int_type, &zbin2_attr);

        // Create new dataspace for attribute
        H5::IntType int64_type(H5::PredType::NATIVE_ULONG);
        attr4_dataspace = H5::DataSpace( H5S_SCALAR );
        zbin_attr = metadataset->createAttribute("nobj0", int64_type, attr4_dataspace);
        zbin_attr.write(int64_type, &nobj1_attr);
        
        // Create new dataspace for attribute
        attr4_dataspace = H5::DataSpace( H5S_SCALAR );
        zbin_attr = metadataset->createAttribute("nobj1", int64_type, attr4_dataspace);
        zbin_attr.write(int64_type, &nobj2_attr);
        cerr << "done with metadata" << endl;
        
        delete metadataspace;
        delete metadataset;
        
        delete corr_group;
        delete jk_group;
        delete metagroup;
        delete file;
        
    }
    
    return;
    
}


int main (int argc, char **argv){
    
    if( argc != 4 && argc != 5 ){
        cerr << "Usage: correlate map1 map2 (radial bin)" << endl;
        cerr << "       where the maps are hdf5 files produced by make_maps, that each have a map and a mask." << endl;
        return 1;
    }

    int r = -1;
    if( argc == 5 ){
        r = atoi(argv[4]);
    }
    
    char outfilename[128] = "correlation.h5";
    correlate( argv[1], argv[2], r, outfilename );
    
    return 0;
    
}
