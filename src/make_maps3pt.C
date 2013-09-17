
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


void make_maps3pt( const char *catalog_filename_c, const char *mask_filename_c, float *angles_c, int n_angles, float *r12s_c, int n_r12s, float *r23s_c, int n_r23s, char *suffix_c ){

    cerr << setprecision(16);

    string catalog_filename(catalog_filename_c);
    string mask_filename(mask_filename_c);
    string suffix(suffix_c);
    
    int min_footprint_order = 3;
    
    // From the smallest angular bin width, figure out the necessary map resolution.
    float min_width = 1.e9;
    for( int i=0; i<n_r12s; ++i ){
        if( r12s_c[i] < min_width )
            min_width = r12s_c[i];
        if( r23s_c[i] < min_width )
            min_width = r23s_c[i];
    }
    int order = 1;
    while(1){
      float pixel_size = sqrt(41253./(12.0*pow(pow(2.0, order), 2.0)));
      if( 3.0*pixel_size < min_width )
          break;
      ++order;
    }
    cerr << "Using order " << order << " for minimum triangle side length " << min_width << endl;
    

    // read in the galaxy catalog
    vector<pointing> pointings;

    // is the input a healpix map rather than a catalog?
    size_t index1 = catalog_filename.rfind(".fits");
    if(index1 != string::npos){
        cerr << "The catalog is a Healpix Map!  Assuming it is a COUNTS map!!" << endl;

        Healpix_Map<double> map;
        read_Healpix_map_from_fits( catalog_filename_c, map );

        // now make a pointings vector as if this came from a catalog
        for( int i=0; i<map.Npix(); ++i ){
            int val = int(map[i]);
            for( int v=0; v<val; ++v ){
                pointings.push_back( map.pix2ang(i) );
            }
        }
    }
    
    else{

        FILE * file = fopen( catalog_filename.c_str(), "rb" );
        if( file == NULL ){
            cerr << "Problem opening " << catalog_filename << endl;
            return;
        }

        double invals[3];
        while( fread( invals, sizeof( double ), 3, file ) == 3 ){
            double phi = invals[0]*3.1415926/180.;
            double theta = (90.-invals[1])*3.1415926/180.;
            pointing mypointing(theta, phi);
            pointings.push_back( mypointing );
        }
        fclose( file );    
        cerr << "Read in " << pointings.size() << " objects" << endl;
    
    }
    
    // read in the mask
    char line_char[256];
    string scheme;
    int mask_order;
    vector<int> mask_pixels;
    Healpix_Base2 mask_base;
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
        FILE * maskfile = fopen( mask_filename.c_str(), "r" );
        if( maskfile == NULL ){
            cerr << "Problem opening " << mask_filename << endl;
            return;
        }
        bool header = true;
        while( fgets( line_char, 256, maskfile ) ){
            if( header ){
                char sch[16];
                int num = sscanf(line_char, "%s %d", sch, &mask_order );
                if( num != 2 ){
                    cerr << "Problem reading header from mask file, number of arguments read = " << num << ", not 2" << endl;
                    return;
                }
                scheme = string(sch);
                if( (!scheme.compare("RING")) && (!scheme.compare("NEST")) ){
                    cerr << "Scheme listed in mask header line is not RING or NEST, but " << scheme << endl;
                    return;
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
                    return;
                }
            }
            if( val > 0.5 )
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
    }
    cerr << "Finished with the mask" << endl;
    
    
    // from the masks, determine the best possible healpix footprint
    // requirements: must cover all of the mask area
    // must cover <80% of the area covered by the next lowest resolution
    // must be order <= 6 ?
    int footprint_order = min_footprint_order;
    double footprint_area = 50000.;
    Healpix_Map<int> *footprintMap = new Healpix_Map<int>(1, RING);
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
    Partpix_Map2<int> *dcMask = new Partpix_Map2<int>(mask_order, *footprintMap);
    dcMask->fill(0);
    for( unsigned int i=0; i<mask_pixels.size(); ++i ){
      if( ! (*footprintMap)[ footprintMap->ang2pix(mask_base.pix2ang(mask_pixels[i])) ] ){
          cerr << "Skipping mask pixel " << i << "!" << endl;
          continue;
      }
      (*dcMask)[mask_pixels[i]] = 1;
    }
    mask_pixels.clear();
    
    cerr << "Finished filling in the mask(s)" << endl;
  
    // system( "rm fp.fits" );
    // fitshandle myfits;
    // myfits.create("fp.fits");
    // write_Healpix_map_to_fits(myfits, footprintMap, PLANCK_FLOAT64);
    // myfits.close();
    
    // make the maps
        
    // always make the delta_counts map since we need that in order to make the delta_mag MASK
    Partpix_Map2<float> *dcMap = new Partpix_Map2<float>(order, *footprintMap);
    dcMap->fill(0.);
    for( unsigned int p=0; p<pointings.size(); ++p ){
        // check if objects are outside the footprint map
        // might happen if objects are outside the mask.
        int fppix = footprintMap->ang2pix(pointings[p]);
        if( (*footprintMap)[fppix] < 0.5 )
            continue;
        int64 pixnum = dcMap->ang2pix( pointings[p] );
        (*dcMap)[pixnum] += 1.0;
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
        cerr << "No counts in the unmasked part of the map??" << endl;
        return;
    }
    mean_counts /= n_pix;
    for( int64 i1=0; i1<dcMap->Npartpix(); ++i1 ){
        int64 i = dcMap->highResPix(i1);
        (*dcMap)[i] = (*dcMap)[i]/mean_counts - 1.0;
    }

    cerr << "Writing the map" << endl;

    // write out the maps (hdf5)
    {
    string outfilename = "dc_map3pt" + suffix + ".h5";
    H5::H5File *file = new H5::H5File( outfilename, H5F_ACC_TRUNC ); // clobber!
    H5::Group* group = new H5::Group( file->createGroup( "/data" ));
    H5::PredType datatype_map( H5::PredType::NATIVE_DOUBLE );
    H5::PredType datatype_mask( H5::PredType::NATIVE_INT64 );
    hsize_t dimsf[2];
    dimsf[1] = 2;
    int rank = 2;

    // the map dataset
    double *data = new double[2*dcMap->Npartpix()+2];
    data[0] = dcMap->Order();
    data[1] = RING;
    unsigned int index = 2;
    for( int64 i1=0; i1<dcMap->Npartpix(); ++i1 ){
        int64 i = dcMap->highResPix(i1);
        // only write it out if it's inside the mask
        if( (*dcMask)[dcMask->ang2pix(dcMap->pix2ang(i))] > 0.5 ){
            data[index] = i;
            ++index;
            data[index] = (*dcMap)[i];
            ++index;
        }
    }
    dimsf[0] = index/2; //dcMap->Npartpix()+1;
    H5::DataSpace *dataspace = new H5::DataSpace( rank, dimsf );
    string datasetname = "/data/map";
    H5::DataSet *dataset = new H5::DataSet( file->createDataSet( datasetname, datatype_map, *dataspace ) );
    dataset->write( data, H5::PredType::NATIVE_DOUBLE );
    delete[] data;
    delete dataspace;
    delete dataset;
    
    cerr << "Wrote " << dimsf[0]-1 << " pixels to the map" << endl;
    delete dcMap;
    
    // the mask dataset
    int64 *data_mask = new int64[2*dcMask->Npartpix()+2];
    data_mask[0] = dcMask->Order();
    data_mask[1] = RING;
    index = 2;
    for( int64 i1=0; i1<dcMask->Npartpix(); ++i1 ){
        int64 i = dcMask->highResPix(i1);
        data_mask[index] = i;
        ++index;
        data_mask[index] = (*dcMask)[i];
        ++index;
    }
    dimsf[0] = dcMask->Npartpix()+1;
    H5::DataSpace *dataspace_mask = new H5::DataSpace( rank, dimsf );
    datasetname = "/data/mask";
    H5::DataSet *dataset_mask = new H5::DataSet( file->createDataSet( datasetname, datatype_mask, *dataspace_mask ) );
    dataset_mask->write( data_mask, H5::PredType::NATIVE_INT64 );
    cerr << "Wrote " << dcMask->Npartpix() << " pixels to the mask" << endl;
    delete dataspace_mask;
    delete dataset_mask;
    delete[] data_mask;
    
    // the metadata
    H5::Group *metagroup = new H5::Group( file->createGroup( "/meta" ) );
    H5::PredType metadatatype( H5::PredType::NATIVE_CHAR );
    //H5::DataType metadatatype( H5T_STRING, 1 );
    //rank = 1;
    dimsf[0] = 1;
    H5::DataSpace *metadataspace = new H5::DataSpace( rank, dimsf );
    H5::DataSet *metadataset = new H5::DataSet( file->createDataSet( "/meta/meta", metadatatype, *metadataspace ) );

    string angles_line = "[";
    char buffer[32];
    sprintf( buffer, "%f", angles_c[0] );
    angles_line.append(buffer);
    for( int i=1; i<n_angles; ++i ){
        angles_line.append(", ");
        sprintf( buffer,"%f", angles_c[i] );
        angles_line.append(buffer);
    }
    angles_line.append("]");

    // Create new dataspace for attribute
    H5::DataSpace attr1_dataspace = H5::DataSpace( H5S_SCALAR );

    // Create new string datatype for attribute
    H5::StrType strdatatype( H5::PredType::C_S1, angles_line.size()+1 ); // of length 256 characters

    // Set up write buffer for attribute
    const H5std_string u_mean_buf( angles_line.c_str() );

    // Create attribute and write to it
    H5::Attribute u_mean_attr = metadataset->createAttribute("angles", strdatatype, attr1_dataspace);
    u_mean_attr.write(strdatatype, u_mean_buf);

    string r12s_line = "[";
    sprintf( buffer, "%f", r12s_c[0] );
    r12s_line.append(buffer);
    for( int i=1; i<n_r12s; ++i ){
        r12s_line.append(", ");
        sprintf( buffer,"%f", r12s_c[i] );
        r12s_line.append(buffer);
    }
    r12s_line.append("]");

    // Create new dataspace for attribute
    H5::DataSpace attr2_dataspace = H5::DataSpace( H5S_SCALAR );

    // Create new string datatype for attribute
    strdatatype = H5::StrType( H5::PredType::C_S1, r12s_line.size()+1 ); // of length 256 characters

    // Set up write buffer for attribute
    const H5std_string r12s_buf( r12s_line.c_str() );

    // Create attribute and write to it
    H5::Attribute r12s_attr = metadataset->createAttribute("r_12s", strdatatype, attr2_dataspace);
    r12s_attr.write(strdatatype, r12s_buf);

    string r23s_line = "[";
    sprintf( buffer, "%f", r23s_c[0] );
    r23s_line.append(buffer);
    for( int i=1; i<n_r23s; ++i ){
        r23s_line.append(", ");
        sprintf( buffer,"%f", r23s_c[i] );
        r23s_line.append(buffer);
    }
    r23s_line.append("]");

    // Create new dataspace for attribute
    H5::DataSpace attr20_dataspace = H5::DataSpace( H5S_SCALAR );

    // Create new string datatype for attribute
    strdatatype = H5::StrType( H5::PredType::C_S1, r23s_line.size()+1 ); // of length 256 characters

    // Set up write buffer for attribute
    const H5std_string r23s_buf( r23s_line.c_str() );

    // Create attribute and write to it
    H5::Attribute r23s_attr = metadataset->createAttribute("r_23s", strdatatype, attr20_dataspace);
    r23s_attr.write(strdatatype, r23s_buf);
    

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
    const H5std_string ftype_buf( "[\"counts\"]" );    
    // Create attribute and write to it
    H5::Attribute ftype_attr = metadataset->createAttribute("ftype", strdatatype, attr4_dataspace);
    ftype_attr.write(strdatatype, ftype_buf);
    
    delete metadataspace;
    delete metadataset;
    
    delete group;
    delete metagroup;
    delete file;
    }
    
    cerr << "Wrote delta counts map to file" << endl;
    delete dcMask;

    
    delete footprintMap;

    cerr << "Finished with make_maps!" << endl;

    return;
}




int main (int argc, char **argv){

    if( argc != 4 ){
        cerr << "Usage: make_maps metadata catalog mask" << endl;
        cerr << "Example input is given in the pxcorr/example_input directory." << endl;
        return 1;
    }
    
    string metadatafilename(argv[1]);
    string catalog_filename(argv[2]);
    string mask_filename(argv[3]);
    
    // read in the binning metadata
    ifstream metadatafile( metadatafilename.c_str() );
    vector<float> angles;
    vector<float> r_12s;
    vector<float> r_23s;
    string line;
    getline( metadatafile, line );
    while( !metadatafile.eof()){
        stringstream linestream(line);
        
        // look for ang_mean
        size_t index1 = line.find("angles");
        if(index1 != string::npos){
            ++index1;
            while( index1 < line.length() && line[index1] != '[' ){
                ++index1;
            }
            if( line[index1] != '[' ){
                cerr << "Problem parsing metadata line for angles" << endl << line << endl;
                return 1;
            }
            ++index1;
            while(1){
                size_t index2 = index1+1;
                while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
                    ++index2;
                }
                float val = atof(line.substr(index1, index2-index1).c_str());
                angles.push_back(val);
                if( line[index2] == ']' )
                    break;
                ++index2;
                index1 = index2;
            }
        }
        
        // look for r_12s
        index1 = line.find("r_12s");
        if(index1 != string::npos){
            ++index1;
            while( index1 < line.length() && line[index1] != '[' ){
                ++index1;
            }
            if( line[index1] != '[' ){
                cerr << "Problem parsing metadata line for r_12s" << endl << line << endl;
                return 1;
            }
            ++index1;
            while(1){
                size_t index2 = index1+1;
                while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
                    ++index2;
                }
                float val = atof(line.substr(index1, index2-index1).c_str());
                r_12s.push_back(val);
                if( line[index2] == ']' )
                    break;
                ++index2;
                index1 = index2;
            }
        }
        
        // look for r_23s
        index1 = line.find("r_23s");
        if(index1 != string::npos){
            ++index1;
            while( index1 < line.length() && line[index1] != '[' ){
                ++index1;
            }
            if( line[index1] != '[' ){
                cerr << "Problem parsing metadata line for r_23s" << endl << line << endl;
                return 1;
            }
            ++index1;
            while(1){
                size_t index2 = index1+1;
                while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
                    ++index2;
                }
                float val = atof(line.substr(index1, index2-index1).c_str());
                r_23s.push_back(val);
                if( line[index2] == ']' )
                    break;
                ++index2;
                index1 = index2;
            }
        }
        
        getline( metadatafile, line );
    }
    cerr << "angles: ";
    for( unsigned int i=0; i<angles.size(); ++i ){
        cerr << angles[i] << " ";
    }
    cerr << endl;
    cerr << "r_12s: ";
    for( unsigned int i=0; i<r_12s.size(); ++i ){
        cerr << r_12s[i] << " ";
    }
    cerr << endl;
    cerr << "r_23s: ";
    for( unsigned int i=0; i<r_23s.size(); ++i ){
        cerr << r_23s[i] << " ";
    }
    cerr << endl;
    if( r_12s.size() == 0 || r_23s.size() == 0 || r_12s.size() != r_23s.size() ){
        cerr << "Different length triangle side arrays??" << endl;
        return 1;
    }

    float angles_c[angles.size()];
    float r_12s_c[r_12s.size()];
    float r_23s_c[r_23s.size()];
    for( unsigned int i=0; i<angles.size(); ++i ){
        angles_c[i] = angles[i];
    }
    for( unsigned int i=0; i<r_12s.size(); ++i ){
        r_12s_c[i] = r_12s[i];
        r_23s_c[i] = r_23s[i];
    }

    char suffix[4] = "ahb";

    make_maps3pt( catalog_filename.c_str(), mask_filename.c_str(), angles_c, angles.size(), r_12s_c, r_12s.size(), r_23s_c, r_23s.size(), suffix );
    
    return 0;
    

}


