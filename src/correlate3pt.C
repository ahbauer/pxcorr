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

template<class T>
string string_ize( const T& val )
{
  ostringstream s;
  s << val;
  return s.str();
}

double distance_between( const double &ra1, const double &dec1, const double &ra2, const double &dec2 ){
    
    double dra2 = 0.0087266461*(ra1-ra2);
    double ddec2 = 0.0087266461*(dec1-dec2);
    return 114.591559*asin( sqrt( sin(ddec2)*sin(ddec2) + cos(0.0174533*dec1)*cos(0.0174533*dec2)*sin(dra2)*sin(dra2) ) ); // degrees
    
}

void correlate3pt( char* mapn1, char* mapn2, char* mapn3, char* sfx ){
    
    bool OUTPUT_FITS = true;
    cout << setiosflags(ios::fixed) << setprecision(16);
    cerr << setprecision(16);
    
    int min_footprint_order = 3;
    
    string mapname1( mapn1 );
    string mapname2( mapn2 );
    string mapname3( mapn3 );
    string suffix( sfx );

    if( access( mapname1.c_str(), F_OK ) == -1 ){
        cerr << "File " << mapname1 << " is not readable!" << endl;
        return;
    }
    if( access( mapname2.c_str(), F_OK ) == -1 ){
        cerr << "File " << mapname2 << " is not readable!" << endl;
        return;
    }
    if( access( mapname3.c_str(), F_OK ) == -1 ){
        cerr << "File " << mapname3 << " is not readable!" << endl;
        return;
    }
    
    vector<float> thetas;
    vector<float> r_12s;
    vector<float> r_23s;
    
    // read in hdf5 maps + masks
    
    // file 1
    H5File file1( mapname1, H5F_ACC_RDONLY );
    H5G_obj_t firstType = file1.getObjTypeByIdx(0);
    if( firstType != H5G_GROUP ){
        cerr << "The map is not standard format;  the first element is not a group." << endl;
        return;
    }
    string group_name = file1.getObjnameByIdx(0);
    if( group_name != "data" ){
        cerr << "Group name is " << group_name << ", not data" << endl;
        return;
    }
    Group group = file1.openGroup("data");
    string dataset1_name = group.getObjnameByIdx(0);
    if( dataset1_name != "map" ){
        cerr << "Dataset 1 name is " << dataset1_name << ", not map" << endl;
        return;
    }
    DataSet dataset1 = group.openDataSet(dataset1_name);
    hsize_t ds1size = dataset1.getStorageSize();
    int npix1 = ds1size/2/sizeof(double);
    double *data1 = (double*) malloc( ds1size );
    dataset1.read(data1, PredType::NATIVE_DOUBLE);
    int map1_order = data1[0];
    int map1_ordering = data1[1];
    cerr << "map 1 info " << map1_order << " " << map1_ordering << " " << npix1-1 << " pixels" << endl;
    // for( int i=1; i<npix1; ++i ){
    //     for( int j=0; j<2; ++j ){
    //         cout << data1[2*i+j] << " ";
    //     }
    //     cout << endl;
    // }
    
    string dataset2_name = group.getObjnameByIdx(1);
    if( dataset2_name != "mask" ){
        cerr << "Dataset 2 name is " << dataset2_name << ", not mask" << endl;
        return;
    }
    DataSet dataset2 = group.openDataSet(dataset2_name);
    hsize_t ds2size = dataset2.getStorageSize();
    long npix2 = ds2size/2/sizeof(long);
    cerr << "mask 1 has " << npix2-1 << " pixels" << endl;
    long *data2 = (long*) malloc( ds2size );
    dataset2.read(data2, PredType::NATIVE_INT64);
    int mask1_order = data2[0];
    int mask1_ordering = data2[1];
    if( mask1_ordering != RING ){
        cerr << "Problem, mask1 ordering is not RING but " << mask1_ordering << endl;
        return;
    }
    cerr << "mask 1 info " << mask1_order << " " << mask1_ordering << endl;
    vector<long> mask1_pixels;
    for( int i=2; i<2*npix2; i+=2 ){
            if( data2[i+1] > 0.5 ){
                mask1_pixels.push_back(data2[i]);
            }
    }
    free(data2);
    Healpix_Base2 mask1_base( mask1_order, RING );
    
    // get the metadata from file 1 (why not)
    group = file1.openGroup("meta");
    DataSet dataSet = group.openDataSet("meta");
    Attribute attr = dataSet.openAttribute("angles");
    int attrsize = attr.getStorageSize();
    StrType dataType(PredType::C_S1, attrsize);

    char *u_mean = (char*) malloc(attrsize);
    attr.read( dataType, u_mean );
    cerr << "angles: " << u_mean << endl;
    string line(u_mean);
    if( line[0] != '[' ){
        cerr << "Problem parsing metadata line for angles" << endl << line << endl;
        return;
    }
    size_t index1 = 1;
    while(1){
        size_t index2 = index1+1;
        while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
            ++index2;
        }
        float val = atof(line.substr(index1, index2-index1).c_str());
        thetas.push_back(val*3.1415926/180.);
        if( line[index2] == ']' )
            break;
        ++index2;
        index1 = index2;
    }

    attr = dataSet.openAttribute("r_12s");
    attrsize = attr.getStorageSize();
    dataType = StrType(PredType::C_S1, attrsize);
    char *u_width = (char*) malloc(attrsize);
    attr.read( dataType, u_width );
    cerr << "r_12s: " << u_width << endl;
    line = string(u_width);
    if( line[0] != '[' ){
        cerr << "Problem parsing metadata line for u_width" << endl << line << endl;
        return;
    }
    index1 = 1;
    while(1){
        size_t index2 = index1+1;
        while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
            ++index2;
        }
        float val = atof(line.substr(index1, index2-index1).c_str());
        r_12s.push_back(val*3.1415926/180.);
        if( line[index2] == ']' )
            break;
        ++index2;
        index1 = index2;
    }
    
    attr = dataSet.openAttribute("r_23s");
    attrsize = attr.getStorageSize();
    dataType = StrType(PredType::C_S1, attrsize);
    char *u_width1 = (char*) malloc(attrsize);
    attr.read( dataType, u_width1 );
    cerr << "r_23s: " << u_width1 << endl;
    line = string(u_width1);
    if( line[0] != '[' ){
        cerr << "Problem parsing metadata line for r_23s" << endl << line << endl;
        return;
    }
    index1 = 1;
    while(1){
        size_t index2 = index1+1;
        while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
            ++index2;
        }
        float val = atof(line.substr(index1, index2-index1).c_str());
        r_23s.push_back(val*3.1415926/180.);
        if( line[index2] == ']' )
            break;
        ++index2;
        index1 = index2;
    }
    
    if( r_12s.size() != r_23s.size() ){
        cerr << "Problem parsing metadata: r_12s and r_23s are different sizes:" << endl;
        cerr << u_width << endl;
        cerr << u_width1 << endl;
        return;
    }
    

    cerr << "thetas: ";
    for( unsigned int i=0; i<thetas.size(); ++i )
        cerr << thetas[i] << " ";
    cerr << endl;
    cerr << "r_12s: ";
    for( unsigned int i=0; i<r_12s.size(); ++i )
        cerr << r_12s[i] << " ";
    cerr << endl;
    cerr << "r_23s: ";
    for( unsigned int i=0; i<r_23s.size(); ++i )
        cerr << r_23s[i] << " ";
    cerr << endl;
    
    
    // file 2
    H5File file2( mapname2, H5F_ACC_RDONLY );
    firstType = file2.getObjTypeByIdx(0);
    if( firstType != H5G_GROUP ){
        cerr << "The map is not standard format;  the first element is not a group." << endl;
        return;
    }
    group_name = file2.getObjnameByIdx(0);
    if( group_name != "data" ){
        cerr << "Group name is " << group_name << ", not data" << endl;
        return;
    }
    group = file2.openGroup("data");
    dataset1_name = group.getObjnameByIdx(0);
    if( dataset1_name != "map" ){
        cerr << "Dataset 1 name is " << dataset1_name << ", not map" << endl;
        return;
    }
    dataset1 = group.openDataSet(dataset1_name);
    ds1size = dataset1.getStorageSize();
    int npix3 = ds1size/2/sizeof(double);
    double *data3 = (double*) malloc( ds1size );
    dataset1.read(data3, PredType::NATIVE_DOUBLE);
    int map2_order = data3[0];
    int map2_ordering = data3[1];
    cerr << "map 2 info " << map2_order << " " << map2_ordering << " " << npix3-1 << " pixels" << endl;
    
    dataset2_name = group.getObjnameByIdx(1);
    if( dataset2_name != "mask" ){
        cerr << "Dataset 2 name is " << dataset2_name << ", not mask" << endl;
        return;
    }
    dataset2 = group.openDataSet(dataset2_name);
    ds2size = dataset2.getStorageSize();
    long npix4 = ds2size/2/sizeof(long);
    cerr << "sizes " << sizeof(int64) << " " << sizeof(PredType::NATIVE_INT64) << " " << sizeof(PredType::NATIVE_INT32) << " " << sizeof(long) << endl;
    cerr << "mask 2 has " << npix4-1 << " pixels" << endl;
    long *data4 = (long*) malloc( ds2size );
    dataset2.read(data4, PredType::NATIVE_INT64);
    int mask2_order = data4[0];
    int mask2_ordering = data4[1];
    if( mask2_ordering != RING ){
        cerr << "Problem, mask2 ordering is not RING but " << mask2_ordering << endl;
        return;
    }
    cerr << "mask 2 info " << mask2_order << " " << mask2_ordering << endl;
    vector<long> mask2_pixels;
    for( int i=2; i<2*npix4; i+=2 ){
            if( data4[i+1] > 0.5 ){
                mask2_pixels.push_back(data4[i]);
            }
    }
    free(data4);
    Healpix_Base2 mask2_base( mask2_order, RING );
    
    
    // file 3
    H5File file3( mapname3, H5F_ACC_RDONLY );
    firstType = file3.getObjTypeByIdx(0);
    if( firstType != H5G_GROUP ){
        cerr << "The map is not standard format;  the first element is not a group." << endl;
        return;
    }
    group_name = file3.getObjnameByIdx(0);
    if( group_name != "data" ){
        cerr << "Group name is " << group_name << ", not data" << endl;
        return;
    }
    group = file3.openGroup("data");
    dataset1_name = group.getObjnameByIdx(0);
    if( dataset1_name != "map" ){
        cerr << "Dataset 1 name is " << dataset1_name << ", not map" << endl;
        return;
    }
    dataset1 = group.openDataSet(dataset1_name);
    ds1size = dataset1.getStorageSize();
    int npix5 = ds1size/2/sizeof(double);
    double *data5 = (double*) malloc( ds1size );
    dataset1.read(data5, PredType::NATIVE_DOUBLE);
    int map3_order = data5[0];
    int map3_ordering = data5[1];
    cerr << "map 3 info " << map3_order << " " << map3_ordering << " " << npix5-1 << " pixels" << endl;
    
    dataset2_name = group.getObjnameByIdx(1);
    if( dataset2_name != "mask" ){
        cerr << "Dataset 2 name is " << dataset2_name << ", not mask" << endl;
        return;
    }
    dataset2 = group.openDataSet(dataset2_name);
    ds2size = dataset2.getStorageSize();
    long npix6 = ds2size/2/sizeof(long);
    cerr << "mask 3 has " << npix6-1 << " pixels" << endl;
    long *data6 = (long*) malloc( ds2size );
    dataset2.read(data6, PredType::NATIVE_INT64);
    int mask3_order = data6[0];
    int mask3_ordering = data6[1];
    if( mask3_ordering != RING ){
        cerr << "Problem, mask3 ordering is not RING but " << mask3_ordering << endl;
        return;
    }
    cerr << "mask 3 info " << mask3_order << " " << mask3_ordering << endl;
    vector<long> mask3_pixels;
    for( int i=2; i<2*npix6; i+=2 ){
            if( data6[i+1] > 0.5 ){
                mask3_pixels.push_back(data6[i]);
            }
    }
    free(data6);
    Healpix_Base2 mask3_base( mask3_order, RING );

    if( (map1_order != map2_order) || (map1_order != map3_order) ){
        cerr << "Map orders must be the same: are " << map1_order << ", " << map2_order << ", and " << map3_order << endl;
        return;
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
        for( unsigned int i=0; i<mask3_pixels.size(); ++i ){
            (*footprintMap)[ footprintMap->ang2pix(mask3_base.pix2ang(mask3_pixels[i])) ] = 1;
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

    if( OUTPUT_FITS ){
        system( "rm fp.fits" );
        fitshandle myfits;
        myfits.create("fp.fits");
        write_Healpix_map_to_fits(myfits, *footprintMap, PLANCK_FLOAT64);
        myfits.close();
    }

    // fill in the mask Partpix maps
    Partpix_Map2<int> *Mask1 = new Partpix_Map2<int>(mask1_order, *footprintMap);
    Mask1->fill(0);
    for( unsigned int i=0; i<mask1_pixels.size(); ++i ){
      if( ! (*footprintMap)[ footprintMap->ang2pix(mask1_base.pix2ang(mask1_pixels[i])) ] ){
          cerr << "Skipping low z mask pixel " << i << "!" << endl;
          continue;
      }
      (*Mask1)[mask1_pixels[i]] = 1;
    }

    mask1_pixels.clear();
    // make full-resolution versions of the masks so we don't keep having to check for mask indices
    Partpix_Map2<int> *MatchedMask1 = new Partpix_Map2<int>(order, *footprintMap);
    if( order > Mask1->Order() ){
        MatchedMask1->Import_upgrade( *Mask1, *footprintMap );
    }
    else if( order == Mask1->Order() ){
        MatchedMask1->Import_nograde( *Mask1, *footprintMap );
    }
    else{
        MatchedMask1->Import_degrade( *Mask1, *footprintMap );
    }

    cerr << "Finished making Partpix matched mask #1" << endl;
    cerr << "Imported from resolution " << mask1_order << " to " << order << endl;

    Partpix_Map2<int> *Mask2 = new Partpix_Map2<int>(mask2_order, *footprintMap);
    Mask2->fill(0);
    for( unsigned int i=0; i<mask2_pixels.size(); ++i ){
      if( ! (*footprintMap)[ footprintMap->ang2pix(mask2_base.pix2ang(mask2_pixels[i])) ] ){
          cerr << "Skipping mask2 pixel " << i << "!" << endl;
          continue;
      }
      (*Mask2)[mask2_pixels[i]] = 1;
    }
    mask2_pixels.clear();
    Partpix_Map2<int> *MatchedMask2 = new Partpix_Map2<int>(order, *footprintMap);
    if( order > Mask2->Order() ){
        MatchedMask2->Import_upgrade( *Mask2, *footprintMap );
    }
    else if( order == Mask2->Order() ){
        MatchedMask2->Import_nograde( *Mask2, *footprintMap );
    }
    else{
        MatchedMask2->Import_degrade( *Mask2, *footprintMap );
    }
    
    cerr << "Finished making Partpix matched mask #2" << endl;
    cerr << "Imported from resolution " << mask2_order << " to " << order << endl;

    Partpix_Map2<int> *Mask3 = new Partpix_Map2<int>(mask3_order, *footprintMap);
    Mask3->fill(0);
    for( unsigned int i=0; i<mask3_pixels.size(); ++i ){
      if( ! (*footprintMap)[ footprintMap->ang2pix(mask3_base.pix2ang(mask3_pixels[i])) ] ){
          cerr << "Skipping mask3 pixel " << i << "!" << endl;
          continue;
      }
      (*Mask3)[mask3_pixels[i]] = 1;
    }
    mask3_pixels.clear();
    Partpix_Map2<int> *MatchedMask3 = new Partpix_Map2<int>(order, *footprintMap);
    if( order > Mask3->Order() ){
        MatchedMask3->Import_upgrade( *Mask3, *footprintMap );
    }
    else if( order == Mask3->Order() ){
        MatchedMask3->Import_nograde( *Mask3, *footprintMap );
    }
    else{
        MatchedMask3->Import_degrade( *Mask3, *footprintMap );
    }

    cerr << "Finished making Partpix matched mask #3" << endl;
    cerr << "Imported from resolution " << mask3_order << " to " << order << endl;


    if( OUTPUT_FITS ){

        system( "rm Mask10.fits" );
        fitshandle myfits;
        myfits.create("Mask10.fits");
        Healpix_Map<int> lowzHMask = Mask1->to_Healpix( 0 );
        write_Healpix_map_to_fits(myfits, lowzHMask, PLANCK_INT32);
        myfits.close();
        
        if( order < 14 ){ 
            system( "rm Mask1.fits" );
            myfits = fitshandle();
            myfits.create("Mask1.fits");
            Healpix_Map<int> lowzHMask0 = MatchedMask1->to_Healpix( 0 );
            write_Healpix_map_to_fits(myfits, lowzHMask0, PLANCK_INT32);
            myfits.close();

            system( "rm Mask2.fits" );
            myfits = fitshandle();
            myfits.create("Mask2.fits");
            Healpix_Map<int> highzHMask = MatchedMask2->to_Healpix( 0 );
            write_Healpix_map_to_fits(myfits, highzHMask, PLANCK_INT32);
            myfits.close();
        }
        else{
            system( "rm Mask1.fits" );
            myfits = fitshandle();
            myfits.create("Mask1.fits");
            Partpix_Map2<int> lowzHMask1( 10, *footprintMap );
            lowzHMask1.Import_degrade( *MatchedMask1, *footprintMap );
            Healpix_Map<int> lowzHMask0 = lowzHMask1.to_Healpix( 0 );
            write_Healpix_map_to_fits(myfits, lowzHMask0, PLANCK_INT32);
            myfits.close();

            system( "rm Mask2.fits" );
            myfits = fitshandle();
            myfits.create("Mask2.fits");
            Partpix_Map2<int> highzHMask1( 10, *footprintMap );
            highzHMask1.Import_degrade( *MatchedMask2, *footprintMap );
            Healpix_Map<int> highzHMask = highzHMask1.to_Healpix( 0 );
            write_Healpix_map_to_fits(myfits, highzHMask, PLANCK_INT32);
            myfits.close();
        }
    }

    delete Mask1;
    delete Mask2;
    delete Mask3;


    // make the correlation maps (partpix)
    Partpix_Map2<float> *Map1 = new Partpix_Map2<float>(order, *footprintMap);
    Map1->fill(0.0);
    Partpix_Map2<float> *Map2 = new Partpix_Map2<float>(order, *footprintMap);
    Map2->fill(0.0);
    Partpix_Map2<float> *Map3 = new Partpix_Map2<float>(order, *footprintMap);
    Map3->fill(0.0);
    cerr << "Created Partpix maps for data, order " << order << endl;

#if USE_WEIGHTS
    Partpix_Map2<float> *WeightMap1 = new Partpix_Map2<float>(order, *footprintMap);
    WeightMap1->fill(1.0);  
    Partpix_Map2<float> *WeightMap2 = new Partpix_Map2<float>(order, *footprintMap);
    WeightMap2->fill(1.0);
    Partpix_Map2<float> *WeightMap3 = new Partpix_Map2<float>(order, *footprintMap);
    WeightMap3->fill(1.0);
    cerr << "Created Partpix maps for weights" << endl;
#endif

    for( int i=1; i<npix1; ++i ){
        // cerr << "trying low z pixel " << data1[2*i] << " value " << data1[2*i+1] << " i = " << i << endl;
        (*Map1)[data1[2*i]] = data1[2*i+1];
    }
    for( int i=1; i<npix3; ++i ){
        // cerr << "trying high z pixel " << data3[2*i] << " value " << data3[2*i+1] << " i = " << i << endl;
        (*Map2)[data3[2*i]] = data3[2*i+1];
    }
    for( int i=1; i<npix5; ++i ){
        // cerr << "trying high z pixel " << data3[2*i] << " value " << data3[2*i+1] << " i = " << i << endl;
        (*Map3)[data5[2*i]] = data5[2*i+1];
    }
    cerr << "Read in the data maps" << endl;

    if( OUTPUT_FITS ){
        Healpix_Map<float> lowzHMap = Map1->to_Healpix( 0. );
        system( "rm Map1.fits" );
        fitshandle myfits;
        myfits.create("Map1.fits");
        write_Healpix_map_to_fits(myfits, lowzHMap, PLANCK_FLOAT64);
        myfits.close();
        Healpix_Map<float> highzHMap = Map2->to_Healpix( 0. );
        system( "rm Map2.fits" );
        myfits = fitshandle();
        myfits.create("Map2.fits");
        write_Healpix_map_to_fits(myfits, highzHMap, PLANCK_FLOAT64);
        myfits.close();
    }

    // figure out the jackknife regions using the mask pixels
    // we want >=50 regions, bigger than the max radial bin.
    // insist on 90% completeness within each jackknife pixel.
    int jk_order = footprintMap->Order();
    vector<int> jkpixels;
    Partpix_Map<int> *jackknifeMap1 = new Partpix_Map<int>(jk_order, *footprintMap);
    jackknifeMap1->fill(0);
    double pixel_size = sqrt(41253./jackknifeMap1->Npix());
    while( pixel_size > r_23s[0] ){
      if( order < jackknifeMap1->Order() ){ // this is unlikely
          Partpix_Map2<int> *mask1 = new Partpix_Map2<int>(jk_order, *footprintMap);
          mask1->Import_upgrade(*MatchedMask1, *footprintMap);
          Partpix_Map2<int> *mask2 = new Partpix_Map2<int>(jk_order, *footprintMap);
          mask2->Import_upgrade(*MatchedMask2, *footprintMap);
          Partpix_Map2<int> *mask3 = new Partpix_Map2<int>(jk_order, *footprintMap);
          mask3->Import_upgrade(*MatchedMask3, *footprintMap);
          for( int i1=0; i1<jackknifeMap1->Npartpix(); ++i1 ){
              int i = jackknifeMap1->highResPix(i1);
              if( (*mask1)[i] == 1 && (*mask2)[i] == 1 && (*mask3)[i] == 1 ){
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
              jackknifeMap1 = new Partpix_Map<int>(jk_order, *footprintMap);
              jackknifeMap1->fill(0.);
              pixel_size = sqrt(41253./jackknifeMap1->Npix());
          }
          delete mask1;
          delete mask2;
      }
      else if( order == jackknifeMap1->Order() ){ // also unlikely
          for( int i1=0; i1<jackknifeMap1->Npartpix(); ++i1 ){
              int i = jackknifeMap1->highResPix(i1);
              if( (*MatchedMask1)[i] == 1 && (*MatchedMask2)[i] == 1 && (*MatchedMask3)[i] == 1 ){
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
              jackknifeMap1 = new Partpix_Map<int>(jk_order, *footprintMap);
              jackknifeMap1->fill(0);
              pixel_size = sqrt(41253./jackknifeMap1->Npix());
          }
      }
      else{ // data are higher resolution than the jackknife map
          for( int i1=0; i1<MatchedMask1->Npartpix(); ++i1 ){
              int64 i = MatchedMask1->highResPix(i1);
              if( (*MatchedMask1)[i] == 1 && (*MatchedMask2)[i] == 1 && (*MatchedMask3)[i] == 1){
                  (*jackknifeMap1)[jackknifeMap1->ang2pix(MatchedMask1->pix2ang(i))]++;
              }
          }
          double threshold = 0.;
          for( int j1=0; j1<jackknifeMap1->Npartpix(); ++j1 ){
              int j = jackknifeMap1->highResPix(j1);
              if( (*jackknifeMap1)[j] > threshold )
                  threshold = (*jackknifeMap1)[j];
          }
          threshold /= 2.0;
          //double threshold = 0.9*(MatchedMask2.Npix()/jackknifeMap1.Npix());
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
              jackknifeMap1 = new Partpix_Map<int>(jk_order, *footprintMap);
              jackknifeMap1->fill(0.);
              pixel_size = sqrt(41253./jackknifeMap1->Npix());
          }
      }
      delete jackknifeMap1;
    }
    if( pixel_size < r_12s[0] ){
      cerr << "Can not find a jackknife resolution with >50 regions and pixel size <" << 180/3.1415926*r_12s[r_12s.size()-1] << " degrees" << endl;
      return;
    }
    Healpix_Base jackknifeMap( jk_order, RING );
    cerr << "Constructed a jackknife map with order " << jk_order << ", "  << jkpixels.size() << " good pixels" << endl;


    int64 n_mask=0;
    for( int64 i1=0; i1<MatchedMask1->Npartpix(); ++i1 ){
      int64 i = MatchedMask1->highResPix(i1);
      if( (*MatchedMask1)[i] && (*MatchedMask2)[i] && (*MatchedMask3)[i] )
          ++n_mask;
    }
    float area_mask = 41252.962*n_mask/MatchedMask1->Npix();
    float area_jk = 41252.962*jkpixels.size()/jackknifeMap.Npix();
    float area_fraction = area_mask/area_jk;
    cerr << "Combined mask (" << area_mask << " sq deg) over jackknife (" << area_jk << " sq deg) area = " << area_fraction << " from " << n_mask << " mask pixels" << endl;
    

    // calculate a jkindex_array to save time in the loops
    vector<int> jkindex_array( jackknifeMap.Npix(), jkpixels.size() );
    for( unsigned int k=0; k<jkpixels.size(); ++k )
        jkindex_array[jkpixels[k]] = k;


    // make a copy of the original maps so that we can use degraded versions in the correlations
    Partpix_Map2<float> *Map1_orig = Map1;
    Partpix_Map2<float> *Map2_orig = Map2;
    Partpix_Map2<float> *Map3_orig = Map3;
    Partpix_Map2<int> *MatchedMask1_orig = MatchedMask1;
    Partpix_Map2<int> *MatchedMask2_orig = MatchedMask2;
    Partpix_Map2<int> *MatchedMask3_orig = MatchedMask3;

    cerr << "starting correlation stuff, with map resolutions " << Map1->Order() << endl;
    
    // for r_12 & r_23 pair
    // # pragma omp parallel for schedule(dynamic, 1)
    for( int rbin=0; rbin<(int)r_12s.size(); ++rbin ){

        float r_12 = r_12s[rbin];
        float r_23 = r_23s[rbin];
        
        vector<double> correlations_ijk(thetas.size(), 0.);
        vector<double> correlations_ij(thetas.size(), 0.);
        vector<double> correlations_jk(thetas.size(), 0.);
        vector<double> correlations_ik(thetas.size(), 0.);
        vector<double> correlations_q3(thetas.size(), 0.);

        vector< vector<double> > pss_ijk( thetas.size(), vector<double>() );
        vector< vector<double> > pss_ij( thetas.size(), vector<double>() );
        vector< vector<double> > pss_jk( thetas.size(), vector<double>() );
        vector< vector<double> > pss_ik( thetas.size(), vector<double>() );
        vector< vector<double> > pss_q3( thetas.size(), vector<double>() );
        vector< vector<double> > jkcorrs_ijk( thetas.size(), vector<double>() );
        vector< vector<double> > jkcorrs_ij( thetas.size(), vector<double>() );
        vector< vector<double> > jkcorrs_jk( thetas.size(), vector<double>() );
        vector< vector<double> > jkcorrs_ik( thetas.size(), vector<double>() );
        vector< vector<double> > jkcorrs_q3( thetas.size(), vector<double>() );
        
        vector<double> jk_means_ijk( thetas.size(), 0. );
        vector<double> jk_means_ij( thetas.size(), 0. );
        vector<double> jk_means_jk( thetas.size(), 0. );
        vector<double> jk_means_ik( thetas.size(), 0. );
        vector<double> jk_means_q3( thetas.size(), 0. );
        vector<double> jk_variance_ijk( thetas.size(), 0. );
        vector<double> jk_variance_ij( thetas.size(), 0. );
        vector<double> jk_variance_jk( thetas.size(), 0. );
        vector<double> jk_variance_ik( thetas.size(), 0. );
        vector<double> jk_variance_q3( thetas.size(), 0. );
        vector<int> jk_ns_ijk( thetas.size(), 0. );
        vector<int> jk_ns_ij( thetas.size(), 0. );
        vector<int> jk_ns_jk( thetas.size(), 0. );
        vector<int> jk_ns_ik( thetas.size(), 0. );
        vector<int> jk_ns_q3( thetas.size(), 0. );

        cerr << "Starting radii " << r_12 << " and " << r_23 << endl;

        for( unsigned int thbin=0; thbin<thetas.size(); ++thbin ){
            float theta = thetas[thbin];

            // make sure that the data resolution is necessary.
            // this is assuming that the bins are never decreasing in width with increasing radius!!
            // cerr << "Changing the map orders from " << Map1_orig->Order() << " to ";
            order = 1;
            while(1){
              float pixel_size = sqrt(41253./(12.0*pow(pow(2.0, order), 2.0)));
              if( 3.0*pixel_size < 57.3*r_23 && 3.0*pixel_size < 57.3*r_12 ) 
                  break;
              ++order;
            }
            // cerr << "order " << order << "... ";

            // DO NOT DEGRADE
            // order = Map1_orig->Order();
            
            if( order < Map1_orig->Order() && order > jackknifeMap.Order() ){
                Map1 = new Partpix_Map2<float>(order, *footprintMap);
                Map2 = new Partpix_Map2<float>(order, *footprintMap);
                Map3 = new Partpix_Map2<float>(order, *footprintMap);
                MatchedMask1 = new Partpix_Map2<int>(order, *footprintMap);
                MatchedMask2 = new Partpix_Map2<int>(order, *footprintMap);
                MatchedMask3 = new Partpix_Map2<int>(order, *footprintMap);
                Map1->Import_degrade(*Map1_orig, *footprintMap);
                Map2->Import_degrade(*Map2_orig, *footprintMap);
                Map3->Import_degrade(*Map3_orig, *footprintMap);
                MatchedMask1->Import_degrade(*MatchedMask1_orig, *footprintMap);
                MatchedMask2->Import_degrade(*MatchedMask2_orig, *footprintMap);
                MatchedMask3->Import_degrade(*MatchedMask3_orig, *footprintMap);
            }
            else if( order < Map1_orig->Order() && order < jackknifeMap.Order()+1 && Map1_orig->Order() > jackknifeMap.Order()+1 ){
                Map1 = new Partpix_Map2<float>(jackknifeMap.Order()+1, *footprintMap);
                Map2 = new Partpix_Map2<float>(jackknifeMap.Order()+1, *footprintMap);
                Map3 = new Partpix_Map2<float>(jackknifeMap.Order()+1, *footprintMap);
                MatchedMask1 = new Partpix_Map2<int>(jackknifeMap.Order()+1, *footprintMap);
                MatchedMask2 = new Partpix_Map2<int>(jackknifeMap.Order()+1, *footprintMap);
                MatchedMask3 = new Partpix_Map2<int>(jackknifeMap.Order()+1, *footprintMap);
                Map1->Import_degrade(*Map1_orig, *footprintMap);
                Map2->Import_degrade(*Map2_orig, *footprintMap);
                Map3->Import_degrade(*Map3_orig, *footprintMap);
                MatchedMask1->Import_degrade(*MatchedMask1_orig, *footprintMap);
                MatchedMask2->Import_degrade(*MatchedMask2_orig, *footprintMap);
                MatchedMask3->Import_degrade(*MatchedMask3_orig, *footprintMap);
            }
            else if( order == Map1_orig->Order() ){
                Map1 = Map1_orig;
                Map2 = Map2_orig;
                Map3 = Map3_orig;
                MatchedMask1 = MatchedMask1_orig;
                MatchedMask2 = MatchedMask2_orig;
                MatchedMask3 = MatchedMask3_orig;
            }
            else{
                cerr << "This combination of orders is not supported!  Map1.Order() = " << Map1->Order() << ", order = " << order << ", jackknifeMap.Order() = " << jackknifeMap.Order() << endl;
                throw;
            }
            order = Map1->Order();
            // cerr << "done!" << endl;

            vector<double> sumijk( jkpixels.size()+1, 0. );
            vector<double> sumij( jkpixels.size()+1, 0. );
            vector<double> sumjk( jkpixels.size()+1, 0. );
            vector<double> sumik( jkpixels.size()+1, 0. );
            vector< vector< vector<double> > > sumijk_injk( jkpixels.size()+1, vector< vector<double> >( jkpixels.size()+1, vector<double>( jkpixels.size()+1, 0. ) ) );
            vector< vector<double> > sumij_injk( jkpixels.size()+1, vector<double>( jkpixels.size()+1, 0. ) );
            vector< vector<double> > sumjk_injk( jkpixels.size()+1, vector<double>( jkpixels.size()+1, 0. ) );
            vector< vector<double> > sumik_injk( jkpixels.size()+1, vector<double>( jkpixels.size()+1, 0. ) );
            vector<double> nijk( jkpixels.size()+1, 0. );
            vector<double> nij( jkpixels.size()+1, 0. );
            vector<double> njk( jkpixels.size()+1, 0. );
            vector<double> nik( jkpixels.size()+1, 0. );
            vector< vector< vector<double> > > nijk_injk( jkpixels.size()+1, vector< vector<double> >( jkpixels.size()+1, vector<double>( jkpixels.size()+1, 0. ) ) );
            vector< vector<double> > nij_injk( jkpixels.size()+1, vector<double>( jkpixels.size()+1, 0. ) );
            vector< vector<double> > njk_injk( jkpixels.size()+1, vector<double>( jkpixels.size()+1, 0. ) );
            vector< vector<double> > nik_injk( jkpixels.size()+1, vector<double>( jkpixels.size()+1, 0. ) );

            // for pixel in high_res map
            double pixel_radius = Map1->max_pixrad();
            for( int64 i1=0; i1<Map1->Npartpix(); ++i1 ){
                int64 i = Map1->highResPix(i1);

                pointing pointing_i = Map1->pix2ang(i);
                if( ! (*MatchedMask1)[i] )
                    continue;
                int jkpix_index_i = jkindex_array[jackknifeMap.ang2pix( pointing_i )];

                // get pixel indices within r_low + annulus_width
                vector<int64> listpix_outer;
                pointing center = Map1->pix2ang(i);
                Map1->query_disc( center, r_12+pixel_radius, listpix_outer );
                
                // get pixel indices within r_low
                vector<int64> listpix_inner;
                Map1->query_disc( center, r_12-pixel_radius, listpix_inner );

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
                
                for( unsigned int j=0; j<listpix_annulus.size(); ++j ){

                    pointing pointing_j = Map2->pix2ang(listpix_annulus[j]);

                    // make sure this pixel is in the high-resolution area!
                    if( (*footprintMap)[footprintMap->ang2pix( pointing_j )] == 0 )
                        continue;

                    if( ! (*MatchedMask2)[listpix_annulus[j]] )
                        continue;

                    if( i == listpix_annulus[j] ){
                        continue;
                    }

                    int jkpix_index_j= jkindex_array[jackknifeMap.ang2pix( pointing_j )];

                    // 12
                    double mult3 = (*Map1)[i]*(*Map2)[listpix_annulus[j]];
                    double multw = 1.0;
                    
                    #if USE_WEIGHTS
                    mult3 *= (*WeightMap1)[i]*(*WeightMap2)[listpix_annulus[j]];
                    multw *= (*WeightMap1)[i]*(*WeightMap2)[listpix_annulus[j]];
                    #endif
                    
                    sumij_injk[jkpix_index_i][jkpix_index_j] += mult3;
                    nij_injk[jkpix_index_i][jkpix_index_j] += multw;

                    // find the 3rd point!
                    double a0 = theta - asin( pointing_i.theta * (sin(pointing_j.phi-pointing_i.phi)/r_12 ) );
                    pointing pointing_l;
                    double pi2mthetaj = 1.5707963-pointing_j.theta;
                    double sinpi2mthetaj = sin(pi2mthetaj);
                    double cospi2mthetaj = cos(pi2mthetaj);
                    double sinr23 = sin(r_23);
                    double cosr23 = cos(r_23);
                    double sina0 = sin(a0);
                    double cosa0 = cos(a0);
                    pointing_l.theta = 1.5707963 - asin( sinpi2mthetaj*cosr23 + cospi2mthetaj*sinr23*cosa0 );
                    pointing_l.phi = pointing_j.phi + atan2( sinr23*sina0, cosr23*cospi2mthetaj - cosa0*sinr23*sinpi2mthetaj );
                    if( pointing_l.phi < 0.0 ){
                        pointing_l.phi = pointing_l.phi + 6.2831852;
                    }
                    int pix_l = Map3->ang2pix( pointing_l );
                    
                    // double ra1 = pointing_i.phi*180.0/3.1415926;
                    // double dec1 = 90.0 - pointing_i.theta*180.0/3.1415926;
                    // double ra2 = pointing_j.phi*180.0/3.1415926;
                    // double dec2 = 90.0 - pointing_j.theta*180.0/3.1415926;
                    // double ra3 = pointing_l.phi*180.0/3.1415926;
                    // double dec3 = 90.0 - pointing_l.theta*180.0/3.1415926;
                    // cout << ra1 << " " << dec1 << " " << ra2 << " " << dec2 << " " << ra3 << " " << dec3 << " theta = " << theta <<  endl;
                    // cout << "distances " << distance_between( ra1, dec1, ra2, dec2 ) << " " << distance_between( ra1, dec1, ra3, dec3 ) << " " << distance_between( ra2, dec2, ra3, dec3 ) << endl;
                    
                    if( (*footprintMap)[footprintMap->ang2pix( pointing_l )] == 1 ){
                        if( (*MatchedMask3)[pix_l] ){
                            int jkpix_index_l = jkindex_array[jackknifeMap.ang2pix( pointing_l )];
                            
                            // 23
                            mult3 = (*Map2)[listpix_annulus[j]]*(*Map3)[pix_l];
                            multw = 1.0;
                            
                            #if USE_WEIGHTS
                            mult3 *= (*WeightMap2)[listpix_annulus[j]]*(*WeightMap3)[pix_l];
                            multw *= (*WeightMap2)[listpix_annulus[j]]*(*WeightMap3)[pix_l];
                            #endif
                            
                            sumjk_injk[jkpix_index_j][jkpix_index_l] += mult3;
                            njk_injk[jkpix_index_j][jkpix_index_l] += multw;
                            
                            // 13
                            mult3 = (*Map1)[i]*(*Map3)[pix_l];
                            multw = 1.0;
                            
                            #if USE_WEIGHTS
                            mult3 *= (*WeightMap1)[i]*(*WeightMap3)[pix_l];
                            multw *= (*WeightMap1)[i]*(*WeightMap3)[pix_l];
                            #endif
                            
                            sumik_injk[jkpix_index_i][jkpix_index_l] += mult3;
                            nik_injk[jkpix_index_i][jkpix_index_l] += multw;

                            // 123
                            mult3 = (*Map1)[i]*(*Map2)[listpix_annulus[j]]*(*Map3)[pix_l];
                            multw = 1.0;

                            #if USE_WEIGHTS
                            mult3 *= (*WeightMap1)[i]*(*WeightMap2)[listpix_annulus[j]]*(*WeightMap3)[pix_l];
                            multw *= (*WeightMap1)[i]*(*WeightMap2)[listpix_annulus[j]]*(*WeightMap3)[pix_l];
                            #endif

                            sumijk_injk[jkpix_index_i][jkpix_index_j][jkpix_index_l] += mult3;
                            nijk_injk[jkpix_index_i][jkpix_index_j][jkpix_index_l] += multw;
                        }
                    }
                    
                    // find the other 3rd point!
                    double theta2 = 6.2831852-theta;
                    a0 = theta2 - asin( pointing_i.theta * (sin(pointing_j.phi-pointing_i.phi)/r_12 ) );
                    pointing_l.theta = 3.1415926 - ( asin( sinpi2mthetaj*cosr23 + cospi2mthetaj*sinr23*cosa0 ) );
                    pointing_l.phi = pointing_j.phi + atan2( sinr23*sina0, cosr23*cospi2mthetaj - cosa0*sinr23*sinpi2mthetaj );
                    pix_l = Map3->ang2pix( pointing_l );
                    if( (*footprintMap)[footprintMap->ang2pix( pointing_l )] == 1 ){
                        if( (*MatchedMask3)[pix_l] ){
                            int jkpix_index_l = jkindex_array[jackknifeMap.ang2pix( pointing_l )];

                            // 23
                            mult3 = (*Map2)[listpix_annulus[j]]*(*Map3)[pix_l];
                            multw = 1.0;
                            
                            #if USE_WEIGHTS
                            mult3 *= (*WeightMap2)[listpix_annulus[j]]*(*WeightMap3)[pix_l];
                            multw *= (*WeightMap2)[listpix_annulus[j]]*(*WeightMap3)[pix_l];
                            #endif
                            
                            sumjk_injk[jkpix_index_j][jkpix_index_l] += mult3;
                            njk_injk[jkpix_index_j][jkpix_index_l] += multw;
                            
                            // 13
                            mult3 = (*Map1)[i]*(*Map3)[pix_l];
                            multw = 1.0;
                            
                            #if USE_WEIGHTS
                            mult3 *= (*WeightMap1)[i]*(*WeightMap3)[pix_l];
                            multw *= (*WeightMap1)[i]*(*WeightMap3)[pix_l];
                            #endif
                            
                            sumik_injk[jkpix_index_i][jkpix_index_l] += mult3;
                            nik_injk[jkpix_index_i][jkpix_index_l] += multw;

                            // 123
                            mult3 = (*Map1)[i]*(*Map2)[listpix_annulus[j]]*(*Map3)[pix_l];
                            multw = 1.0;

                            #if USE_WEIGHTS
                            mult3 *= (*WeightMap1)[i]*(*WeightMap2)[listpix_annulus[j]]*(*WeightMap3)[pix_l];
                            multw *= (*WeightMap1)[i]*(*WeightMap2)[listpix_annulus[j]]*(*WeightMap3)[pix_l];
                            #endif

                            sumijk_injk[jkpix_index_i][jkpix_index_j][jkpix_index_l] += mult3;
                            nijk_injk[jkpix_index_i][jkpix_index_j][jkpix_index_l] += multw;
                        }
                    }

                } // for j pixels in map 2
                
            }// for i pixels in map 1

            // now rearrange jackknife stuff to exclude the pixel
            for( unsigned int i=0; i<jkpixels.size()+1; ++i ){
                for( unsigned int j=0; j<jkpixels.size()+1; ++j ){
                    for( unsigned int l=0; l<jkpixels.size()+1; ++l ){
                        for( unsigned int k=0; k<jkpixels.size(); ++k ){
                            if( k != i && k != j && k != l ){
                                sumijk[k] += sumijk_injk[i][j][l];
                                nijk[k] += nijk_injk[i][j][l];
                            }
                            sumijk[jkpixels.size()] += sumijk_injk[i][j][l];
                            nijk[jkpixels.size()] += nijk_injk[i][j][l];
                        }
                        if( l != i && l != j ){
                            sumij[l] += sumij_injk[i][j];
                            nij[l] += nij_injk[i][j];
                            sumjk[l] += sumjk_injk[i][j];
                            njk[l] += njk_injk[i][j];
                            sumik[l] += sumik_injk[i][j];
                            nik[l] += nik_injk[i][j];
                        }
                        sumij[jkpixels.size()] += sumij_injk[i][j];
                        nij[jkpixels.size()] += nij_injk[i][j];
                        sumjk[jkpixels.size()] += sumjk_injk[i][j];
                        njk[jkpixels.size()] += njk_injk[i][j];
                        sumik[jkpixels.size()] += sumik_injk[i][j];
                        nik[jkpixels.size()] += nik_injk[i][j];
                    }
                }
            }
            
            
            // construct q3
            vector<double> q3( jkpixels.size()+1, 0. );
            for( unsigned int k=0; k<jkpixels.size()+1; ++k ){
                double w_h = sumij[k]*sumjk[k] + sumik[k]*sumjk[k] + sumik[k]*sumij[k];
                q3[k] = sumijk[k]/w_h;
            }
            
            
            if( nijk[jkpixels.size()] > 0 ){
                correlations_ijk[thbin] = sumijk[jkpixels.size()]/nijk[jkpixels.size()];
            }
            if( nij[jkpixels.size()] > 0 ){
                correlations_ij[thbin] = sumij[jkpixels.size()]/nij[jkpixels.size()];
            }
            if( njk[jkpixels.size()] > 0 ){
                correlations_jk[thbin] = sumjk[jkpixels.size()]/njk[jkpixels.size()];
            }
            if( nik[jkpixels.size()] > 0 ){
                correlations_ik[thbin] = sumik[jkpixels.size()]/nik[jkpixels.size()];
            }
            if( nijk[jkpixels.size()] > 0 && nij[jkpixels.size()] > 0 && njk[jkpixels.size()] > 0 && nik[jkpixels.size()] > 0 ){
                correlations_q3[thbin] = correlations_ijk[thbin] / (correlations_ij[thbin]*correlations_jk[thbin] + correlations_ik[thbin]*correlations_jk[thbin] + correlations_ik[thbin]*correlations_ij[thbin]);
            }
            
            // now do the jackknife
            for( unsigned int k=0; k<jkpixels.size(); ++k ){
                double corr_wo1 = 0.0;
                if( nijk[k] > 0.0 ){
                    corr_wo1 = sumijk[k]/nijk[k];
                }
                else{
                    cerr << "Warning, empty ijk jackknife pixel being used..." << endl;
                }
                double ps = jkpixels.size()*correlations_ijk[thbin] - (jkpixels.size()-1)*corr_wo1;
                pss_ijk[thbin].push_back(ps);
                jkcorrs_ijk[thbin].push_back(corr_wo1);
              
                corr_wo1 = 0.0;
                if( nij[k] > 0.0 ){
                    corr_wo1 = sumij[k]/nij[k];
                }
                else{
                    cerr << "Warning, empty ij jackknife pixel being used..." << endl;
                }
                ps = jkpixels.size()*correlations_ij[thbin] - (jkpixels.size()-1)*corr_wo1;
                pss_ij[thbin].push_back(ps);
                jkcorrs_ij[thbin].push_back(corr_wo1);
                
                corr_wo1 = 0.0;
                if( njk[k] > 0.0 ){
                    corr_wo1 = sumjk[k]/njk[k];
                }
                else{
                    cerr << "Warning, empty jk jackknife pixel being used..." << endl;
                }
                ps = jkpixels.size()*correlations_jk[thbin] - (jkpixels.size()-1)*corr_wo1;
                pss_jk[thbin].push_back(ps);
                jkcorrs_jk[thbin].push_back(corr_wo1);
                
                corr_wo1 = 0.0;
                if( nik[k] > 0.0 ){
                    corr_wo1 = sumik[k]/nik[k];
                }
                else{
                    cerr << "Warning, empty ik jackknife pixel being used..." << endl;
                }
                ps = jkpixels.size()*correlations_ik[thbin] - (jkpixels.size()-1)*corr_wo1;
                pss_ik[thbin].push_back(ps);
                jkcorrs_ik[thbin].push_back(corr_wo1);

                corr_wo1 = 0.0;
                if( nijk[k] > 0 && nij[k] > 0 && njk[k] > 0 && nik[k] > 0 ){
                    corr_wo1 = (sumijk[k]/nijk[k]) / ( (sumij[k]/nij[k])*(sumjk[k]/njk[k]) + (sumik[k]/nik[k])*(sumjk[k]/njk[k]) + (sumik[k]/nik[k])*(sumij[k]/nij[k]) );
                }
                else{
                    cerr << "Warning, empty q3 jackknife pixel being used..." << endl;
                }
                ps = jkpixels.size()*correlations_q3[thbin] - (jkpixels.size()-1)*corr_wo1;
                pss_q3[thbin].push_back(ps);
                jkcorrs_q3[thbin].push_back(corr_wo1);
            }

            // now average up the jackknife results for each radius
            for( unsigned int p=0; p<pss_ijk[thbin].size(); ++p ){
                jk_means_ijk[thbin] += jkcorrs_ijk[thbin][p];
                jk_means_ij[thbin] += jkcorrs_ij[thbin][p];
                jk_means_jk[thbin] += jkcorrs_jk[thbin][p];
                jk_means_ik[thbin] += jkcorrs_ik[thbin][p];
                jk_means_q3[thbin] += jkcorrs_q3[thbin][p];
            }
            if( pss_ijk[thbin].size()>2 ){
              jk_means_ijk[thbin] /= jkcorrs_ijk[thbin].size();
              jk_means_ij[thbin] /= jkcorrs_ij[thbin].size();
              jk_means_jk[thbin] /= jkcorrs_jk[thbin].size();
              jk_means_ik[thbin] /= jkcorrs_ik[thbin].size();
              jk_means_q3[thbin] /= jkcorrs_q3[thbin].size();

              for( unsigned int p=0; p<pss_ijk[thbin].size(); ++p ){
                  jk_variance_ijk[thbin] += (jkcorrs_ijk[thbin][p]-jk_means_ijk[thbin])*(jkcorrs_ijk[thbin][p]-jk_means_ijk[thbin]);
                  jk_variance_ij[thbin] += (jkcorrs_ij[thbin][p]-jk_means_ij[thbin])*(jkcorrs_ij[thbin][p]-jk_means_ij[thbin]);
                  jk_variance_jk[thbin] += (jkcorrs_jk[thbin][p]-jk_means_jk[thbin])*(jkcorrs_jk[thbin][p]-jk_means_jk[thbin]);
                  jk_variance_ik[thbin] += (jkcorrs_ik[thbin][p]-jk_means_ik[thbin])*(jkcorrs_ik[thbin][p]-jk_means_ik[thbin]);
                  jk_variance_q3[thbin] += (jkcorrs_q3[thbin][p]-jk_means_q3[thbin])*(jkcorrs_q3[thbin][p]-jk_means_q3[thbin]);
              }
              float fraction = (jkcorrs_ijk[thbin].size()-1.)/jkcorrs_ijk[thbin].size();
              jk_variance_ijk[thbin] *= (fraction/area_fraction);
              jk_variance_ij[thbin] *= (fraction/area_fraction);
              jk_variance_jk[thbin] *= (fraction/area_fraction);
              jk_variance_ik[thbin] *= (fraction/area_fraction);
              jk_variance_q3[thbin] *= (fraction/area_fraction);
            }
            else{
              jk_means_ijk[thbin] = 0.;
              jk_variance_ijk[thbin] = 0.;
              jk_means_ij[thbin] = 0.;
              jk_variance_ij[thbin] = 0.;
              jk_means_jk[thbin] = 0.;
              jk_variance_jk[thbin] = 0.;
              jk_means_ik[thbin] = 0.;
              jk_variance_ik[thbin] = 0.;
              jk_means_q3[thbin] = 0.;
              jk_variance_q3[thbin] = 0.;
            }

            cerr << 180.0/3.1415926*theta << " " << correlations_q3[thbin] << " " << jk_means_q3[thbin] << " " << sqrt(jk_variance_q3[thbin]) 
                 << " " << correlations_ijk[thbin] << " " << jk_means_ijk[thbin] << " " << sqrt(jk_variance_ijk[thbin]) 
                 << " " << correlations_ij[thbin] << " " << jk_means_ij[thbin] << " " << sqrt(jk_variance_ij[thbin]) 
                 << " " << correlations_jk[thbin] << " " << jk_means_jk[thbin] << " " << sqrt(jk_variance_ik[thbin]) 
                 << " " << correlations_ik[thbin] << " " << jk_means_ik[thbin] << " " << sqrt(jk_variance_ik[thbin]) << endl;

            if( MatchedMask1 != MatchedMask1_orig )
                delete MatchedMask1;
            if( MatchedMask2 != MatchedMask2_orig )
                delete MatchedMask2;
            if( MatchedMask3 != MatchedMask3_orig )
                delete MatchedMask3;
            if( Map1 != Map1_orig )
                delete Map1;
            if( Map2 != Map2_orig )
                delete Map2;
            if( Map3 != Map3_orig )
                delete Map3;

        } // for theta

        // write out results to file (3pt)
        stringstream corr_line;
        stringstream jk_line;
        stringstream err_line;
        corr_line << "[" << correlations_ijk[0];
        jk_line << "[" << jk_means_ijk[0];
        err_line << "[" << sqrt(jk_variance_ijk[0]);
        for( unsigned int i=1; i<correlations_ijk.size(); ++i ){
            corr_line << ", " << correlations_ijk[i];
            jk_line << ", " << jk_means_ijk[i];
            err_line << ", " << sqrt(jk_variance_ijk[i]);
        }
        corr_line << "]";
        jk_line << "]";
        err_line << "]";
        string corr_filename = "corr_3pt_r" + string_ize(rbin) + suffix;
        ofstream corr_file(corr_filename.c_str());
        corr_file << "[" << corr_line.str() << ", " << jk_line.str() << ", " << err_line.str() << "]" << endl;
        corr_file.close();

        // and q3
        corr_line.str("");
        jk_line.str("");
        err_line.str("");
        corr_line << "[" << correlations_q3[0];
        jk_line << "[" << jk_means_q3[0];
        err_line << "[" << sqrt(jk_variance_q3[0]);
        for( unsigned int i=1; i<correlations_q3.size(); ++i ){
            corr_line << ", " << correlations_q3[i];
            jk_line << ", " << jk_means_q3[i];
            err_line << ", " << sqrt(jk_variance_q3[i]);
        }
        corr_line << "]";
        jk_line << "]";
        err_line << "]";
        corr_filename = "corr_q3_r" + string_ize(rbin) + suffix;
        ofstream corrq3_file(corr_filename.c_str());
        corrq3_file << "[" << corr_line.str() << ", " << jk_line.str() << ", " << err_line.str() << "]" << endl;
        corrq3_file.close();
        

        // calculate the covariance matrix (3pt)
        vector< vector<double> > c_ijk( thetas.size(), vector<double>(thetas.size(), 0.) );
        cerr << "calculating covariance... " << endl;
        for( unsigned int i=0; i<thetas.size(); ++i ){
          for( unsigned int j=0; j<thetas.size(); ++j ){
            if( jkcorrs_ijk[i].size() != jkcorrs_ijk[j].size() ){
                cerr << "Problem, jackknife lists " << i << ", " << j << " are not the same size: " << jkcorrs_ijk[i].size() << ", " << jkcorrs_ijk[j].size() << endl;
                return;
            }
        for( unsigned int k=0; k<pss_ijk[i].size(); ++k ){
          c_ijk[i][j] += (jkcorrs_ijk[i][k]-jk_means_ijk[i])*(jkcorrs_ijk[j][k]-jk_means_ijk[j]);
        }
        float fraction = (jkcorrs_ijk[i].size()-1.)/(jkcorrs_ijk[i].size());
        c_ijk[i][j] *= fraction;
        c_ijk[i][j] /= area_fraction;

          }
        }
        
        string cov_filename = "cov_3pt_r" + string_ize(rbin) + suffix;
        ofstream cov_file(cov_filename.c_str());
        cov_file << "[";
        for( unsigned int k=0; k<thetas.size(); ++k ){
            if( k == 0 )
                cov_file << "[" << c_ijk[k][0];
            else
                cov_file << ", [" << c_ijk[k][0];
            for( unsigned int l=1; l<thetas.size(); ++l ){
                cov_file << ", " << c_ijk[k][l];
            }
            cov_file << "]";
        }
        cov_file << "]" << endl;
        cov_file.close();
        
        // and q3
        vector< vector<double> > c_q3( thetas.size(), vector<double>(thetas.size(), 0.) );
        cerr << "calculating covariance... " << endl;
        for( unsigned int i=0; i<thetas.size(); ++i ){
          for( unsigned int j=0; j<thetas.size(); ++j ){
            if( jkcorrs_q3[i].size() != jkcorrs_q3[j].size() ){
                cerr << "Problem, q3 jackknife lists " << i << ", " << j << " are not the same size: " << jkcorrs_q3[i].size() << ", " << jkcorrs_q3[j].size() << endl;
                return;
            }
            for( unsigned int k=0; k<pss_q3[i].size(); ++k ){
              c_q3[i][j] += (jkcorrs_q3[i][k]-jk_means_q3[i])*(jkcorrs_q3[j][k]-jk_means_q3[j]);
            }
            float fraction = (jkcorrs_q3[i].size()-1.)/(jkcorrs_q3[i].size());
            c_q3[i][j] *= fraction;
            c_q3[i][j] /= area_fraction;
          }
        }

        cov_filename = "cov_q3_r" + string_ize(rbin) + suffix;
        ofstream covq3_file(cov_filename.c_str());
        covq3_file << "[";
        for( unsigned int k=0; k<thetas.size(); ++k ){
            if( k == 0 )
                covq3_file << "[" << c_q3[k][0];
            else
                covq3_file << ", [" << c_q3[k][0];
            for( unsigned int l=1; l<thetas.size(); ++l ){
                covq3_file << ", " << c_q3[k][l];
            }
            covq3_file << "]";
        }
        covq3_file << "]" << endl;
        covq3_file.close();

    } // for rbin in r_12 and r_23 list
    
    delete MatchedMask1_orig;
    delete MatchedMask2_orig;
    delete MatchedMask3_orig;
    delete Map1_orig;
    delete Map2_orig;
    delete Map3_orig;
    delete footprintMap;

#if USE_WEIGHTS
    delete WeightMap1;
    delete WeightMap2;
    delete WeightMap3;
#endif
    
    return;
    
}


int main (int argc, char **argv){
    
    if( argc != 5 ){
        cerr << "Usage: correlate3pt map1 map2 map3 suffix" << endl;
        cerr << "       where the maps are hdf5 files produced by make_maps, that each have a map and a mask." << endl;
        return 1;
    }
    
    correlate3pt( argv[1], argv[2], argv[3], argv[4] );
    
    return 0;
    
}
