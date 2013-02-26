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

int main (int argc, char **argv){
    
    
    int min_footprint_order = 6;
    
    // read in the binning metadata
    string metadatafilename( argv[1] );
    ifstream metadatafile( metadatafilename.c_str() );
    vector<float> z_means;
    vector<float> z_widths;
    vector<float> ang_means;
    vector<float> ang_widths;
    string line;
    getline( metadatafile, line );
    while( !metadatafile.eof()){
        stringstream linestream(line);
        
        // look for z_mean
        size_t index1 = line.find("z_mean");
        if(index1 != string::npos){
            ++index1;
            while( index1 < line.length() && line[index1] != '[' ){
                ++index1;
            }
            if( line[index1] != '[' ){
                cerr << "Problem parsing metadata line for z_mean" << endl << line << endl;
                return 1;
            }
            ++index1;
            while(1){
                size_t index2 = index1+1;
                while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
                    ++index2;
                }
                float val = atof(line.substr(index1, index2-index1).c_str());
                z_means.push_back(val);
                if( line[index2] == ']' )
                    break;
                ++index2;
                index1 = index2;
            }
        }
        
        // look for z_width
        index1 = line.find("z_width");
        if(index1 != string::npos){
            ++index1;
            while( index1 < line.length() && line[index1] != '[' ){
                ++index1;
            }
            if( line[index1] != '[' ){
                cerr << "Problem parsing metadata line for z_mean" << endl << line << endl;
                return 1;
            }
            ++index1;
            while(1){
                size_t index2 = index1+1;
                while( index2 < line.length() && (line[index2] != ',' && line[index2] != ']') ){
                    ++index2;
                }
                float val = atof(line.substr(index1, index2-index1).c_str());
                z_widths.push_back(val);
                if( line[index2] == ']' )
                    break;
                ++index2;
                index1 = index2;
            }
        }
        
        // look for ang_mean
        index1 = line.find("ang_mean");
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
                cerr << "Problem parsing metadata line for z_mean" << endl << line << endl;
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
    cerr << "z_means: ";
    for( unsigned int i=0; i<z_means.size(); ++i )
        cerr << z_means[i] << " ";
    cerr << endl;
    cerr << "z_widths: ";
    for( unsigned int i=0; i<z_widths.size(); ++i )
        cerr << z_widths[i] << " ";
    cerr << endl;
    cerr << "ang_means: ";
    for( unsigned int i=0; i<ang_means.size(); ++i )
        cerr << ang_means[i] << " ";
    cerr << endl;
    cerr << "ang_widths: ";
    for( unsigned int i=0; i<ang_widths.size(); ++i )
        cerr << ang_widths[i] << " ";
    cerr << endl;
    if( z_means.size() == 0 || z_widths.size() == 0 || z_means.size() != z_widths.size() ){
        cerr << "Different length z_means and z_widths arrays??" << endl;
        return 1;
    }
    if( ang_means.size() == 0 || ang_widths.size() == 0 || ang_means.size() != ang_widths.size() ){
        cerr << "Different length ang_means and ang_widths arrays??" << endl;
        return 1;
    }
    
    // From the smallest angular bin width, figure out the necessary map resolution.
    float min_width = 1.e9;
    for( unsigned int i=0; i<ang_widths.size(); ++i ){
        if( ang_widths[i] < min_width )
            min_width = ang_widths[i];
    }
    int order = 1;
    while(1){
      float pixel_size = sqrt(41253./(12.0*pow(pow(2.0, order), 2.0)));
      if( 0.85*pixel_size < min_width )  // this 0.85 is a bit arbitrary, but it's ok to use pixels a bit bigger than the radial resolution.
          break;
      ++order;
    }
    cerr << "Using order " << order << " for minimum angular width " << min_width << endl;
    
    // read in the instructions about which maps to make
    string catalog_filename(argv[2]);
    string mask_filename(argv[3]);
    
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
    

    // read in the galaxy catalog
    FILE * file = fopen( catalog_filename.c_str(), "r" );
    if( file == NULL ){
        cerr << "Problem opening " << catalog_filename << endl;
        return 1;
    }
    char line_char[256];
    vector<pointing> pointings;
    vector<double> mags;
    while( fgets( line_char, 256, file ) ){
        double ra, dec, z, mag;
        int num = sscanf(line_char, "%lf %lf %lf %lf", &ra, &dec, &z, &mag );
        if( num != 4 ){
            cerr << "Problem reading from file, number of arguments read = " << num << ", not 4" << endl;
            return 1;
        }
        double phi = ra*3.1415926/180.;
        double theta = (90.-dec)*3.1415926/180.;
        pointing mypointing(theta, phi);
        pointings.push_back( mypointing );
        if( use_mags )
            mags.push_back( mag );
    }
    fclose( file );    
    cerr << "Read in " << pointings.size() << " objects" << endl;

    // read in the mask
    FILE * maskfile = fopen( mask_filename.c_str(), "r" );
    if( maskfile == NULL ){
        cerr << "Problem opening " << mask_filename << endl;
        return 1;
    }
    string scheme;
    int mask_order;
    vector<int> mask_pixels;
    bool header = true;
    while( fgets( line_char, 256, maskfile ) ){
        if( header ){
            char sch[16];
            int num = sscanf(line_char, "%s %d", sch, &mask_order );
            if( num != 2 ){
                cerr << "Problem reading header from mask file, number of arguments read = " << num << ", not 2" << endl;
                return 1;
            }
            scheme = string(sch);
            if( (!scheme.compare("RING")) && (!scheme.compare("NEST")) ){
                cerr << "Scheme listed in mask header line is not RING or NEST, but " << scheme << endl;
                return 1;
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
                return 1;
            }
        }
        if( val > 0.5 )
            mask_pixels.push_back(pix);
    }
    cerr << "Read in mask information, scheme " << scheme << ", order " << mask_order << ", " << mask_pixels.size() << " pixels" << endl;
    Healpix_Base2 mask_base( mask_order, RING );
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
    cerr << "Finished with the mask" << endl;
    
    
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
        for( unsigned int i=0; i<mask_pixels.size(); ++i ){
            footprintMap[ footprintMap.ang2pix(mask_base.pix2ang(mask_pixels[i])) ] = 1;
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
            footprint_area = footprint_area_new;
            break;
        }
    }
    cerr << "Using footprint order " << footprintMap.Order() << " with " << footprint_area << " sq degrees." << endl;

    // fill in the Partpix mask(s)
    Partpix_Map2<double> dcMask(mask_order, footprintMap);
    dcMask.fill(0);
    Partpix_Map2<double> dmMask;
    if( use_mags )
        dmMask = Partpix_Map2<double>(mask_order, footprintMap);
    for( unsigned int i=0; i<mask_pixels.size(); ++i ){
      if( ! footprintMap[ footprintMap.ang2pix(mask_base.pix2ang(mask_pixels[i])) ] ){
          cerr << "Skipping mask pixel " << i << "!" << endl;
          continue;
      }
      dcMask[mask_pixels[i]] = 1.0;
      if( use_mags )
        dmMask[mask_pixels[i]] = 1.0;
    }
    mask_pixels.clear();
    cerr << "Finished filling in the mask(s)" << endl;
  
    // system( "rm fp.fits" );
    // fitshandle myfits;
    // myfits.create("fp.fits");
    // write_Healpix_map_to_fits(myfits, footprintMap, PLANCK_FLOAT64);
    // myfits.close();
    
    // make the maps
    // one at a time, for kicks (maybe will parallelize later?)
    for( unsigned int z=0; z<z_means.size(); ++z ){
        
        // always make the delta_counts map since we need that in order to make the delta_mag MASK
        Partpix_Map2<double> dcMap(order, footprintMap);
        dcMap.fill(0.);
        Partpix_Map2<double> dmMap;
        double mean_mag = 0.;
        int n_mags = 0;
        if( use_mags ){
            dmMap = Partpix_Map2<double>(order, footprintMap);
            dmMap.fill(0.);
        }
        for( unsigned int p=0; p<pointings.size(); ++p ){
            int64 pixnum = dcMap.ang2pix( pointings[p] );
            dcMap[pixnum] += 1.0;
            if( use_mags ){
                int64 maskpix = dmMask.ang2pix(dmMap.pix2ang(pixnum));
                if( dmMask[maskpix] > 0.5 ){
                    dmMap[pixnum] += mags[p];
                    mean_mag += mags[p];
                    n_mags++;
                }
            }
        }
        
        // now make the delta_mag mask for pixels with counts=0
        // mask out the pixels that don't have any mag info.
        if( use_mags ){
            mean_mag /= n_mags;
            for( int i1=0; i1<dmMap.Npartpix(); ++i1 ){
                int i = dmMap.highResPix(i1);
                if( dcMap[i] > 0. ){
                    dmMap[i] = dmMap[i]/dcMap[i] - mean_mag;
                }
                else{
                    int64 maskpix = dmMask.ang2pix(dmMap.pix2ang(i));
                    dmMask[maskpix] = 0.0;
                }
            }
        }
        
        // now turn the counts map into a delta map
        double mean_counts = 0.0;
        double n_pix = 0.0;
        for( int64 i1=0; i1<dcMap.Npartpix(); ++i1 ){
            int64 i = dcMap.highResPix(i1);
            int64 maskpix = dcMask.ang2pix(dcMap.pix2ang(i));
            if( dcMask[maskpix] > 0.5 ){
                mean_counts += dcMap[i];
                n_pix += 1.0;
            }
        }
        if( n_pix == 0.0 ){
            cerr << "No counts in the unmasked part of the map??" << endl;
            return 1;
        }
        mean_counts /= n_pix;
        for( int64 i1=0; i1<dcMap.Npartpix(); ++i1 ){
            int64 i = dcMap.highResPix(i1);
            dcMap[i] = dcMap[i]/mean_counts - 1.0;
        }
        
        // counts:
        {
            // write out the maps (hdf5)
            string outfilename = "dc_maps.h5";
            H5File *file = new H5File( outfilename, H5F_ACC_TRUNC ); // clobber!
            Group* group = new Group( file->createGroup( "/data" ));
            PredType datatype( PredType::NATIVE_DOUBLE );
            hsize_t dimsf[2];
            dimsf[1] = 2;
            int rank = 2;

            // the map dataset
            double *data = (double*) malloc( sizeof(double) * 2*dcMap.Npartpix()+2 );
            data[0] = dcMap.Order();
            data[1] = RING;
            unsigned int index = 2;
            for( int64 i1=0; i1<dcMap.Npartpix(); ++i1 ){
                int64 i = dcMap.highResPix(i1);
                data[index] = i;
                ++index;
                data[index] = dcMap[i];
                ++index;
            }
            dimsf[0] = dcMap.Npartpix()+1;
            DataSpace *dataspace = new DataSpace( rank, dimsf );
            string datasetname = "/data/map";
            DataSet *dataset = new DataSet( file->createDataSet( datasetname, datatype, *dataspace ) );
            dataset->write( data, PredType::NATIVE_DOUBLE );
            delete dataspace;
            delete dataset;
            free( data );
        
            // the mask dataset
            double *data_mask = (double*) malloc( sizeof(double) * 2*dcMask.Npartpix()+2 );
            data_mask[0] = dcMask.Order();
            data_mask[1] = RING;
            index = 2;
            for( int64 i1=0; i1<dcMask.Npartpix(); ++i1 ){
                int64 i = dcMask.highResPix(i1);
                data_mask[index] = i;
                ++index;
                data_mask[index] = dcMask[i];
                ++index;
            }
            dimsf[0] = dcMask.Npartpix()+1;
            DataSpace *dataspace_mask = new DataSpace( rank, dimsf );
            datasetname = "/data/mask";
            DataSet *dataset_mask = new DataSet( file->createDataSet( datasetname, datatype, *dataspace_mask ) );
            dataset->write( data_mask, PredType::NATIVE_DOUBLE );
            delete dataspace_mask;
            delete dataset_mask;
            free( data_mask );
            
            // the metadata
            Group* metagroup = new Group( file->createGroup( "/meta" ) );
            DataType metadatatype( H5T_STRING, 1 );
            rank = 1;
            dimsf[0] = 1;
            DataSpace *metadataspace = new DataSpace( rank, dimsf );
            DataSet *metadataset = new DataSet( file->createDataSet( "/meta/meta", metadatatype, *metadataspace ) );

            stringstream ang_mean_line;
            ang_mean_line << "[" << ang_means[0];
            for( unsigned int i=1; i<ang_means.size(); ++i ){
                ang_mean_line << ", " << ang_means[i];
            }
            ang_mean_line << "]";

            // Create new dataspace for attribute
            DataSpace attr1_dataspace = DataSpace( H5S_SCALAR );
    
            // Create new string datatype for attribute
            StrType strdatatype( PredType::C_S1, ang_mean_line.str().size()+1 ); // of length 256 characters
    
            // Set up write buffer for attribute
            const H5std_string u_mean_buf( ang_mean_line.str().c_str() );
    
            // Create attribute and write to it
            Attribute u_mean_attr = metadataset->createAttribute("u_mean", strdatatype, attr1_dataspace);
            u_mean_attr.write(strdatatype, u_mean_buf);

            stringstream ang_widths_line;
            ang_widths_line << "[" << ang_widths[0];
            for( unsigned int i=1; i<ang_widths.size(); ++i ){
                ang_widths_line << ", " << ang_widths[i];
            }
            ang_widths_line << "]";

            // Create new dataspace for attribute
            DataSpace attr2_dataspace = DataSpace( H5S_SCALAR );
    
            // Create new string datatype for attribute
            strdatatype = StrType( PredType::C_S1, ang_widths_line.str().size()+1 ); // of length 256 characters
    
            // Set up write buffer for attribute
            const H5std_string u_widths_buf( ang_widths_line.str().c_str() );
    
            // Create attribute and write to it
            Attribute u_width_attr = metadataset->createAttribute("u_width", strdatatype, attr2_dataspace);
            u_width_attr.write(strdatatype, u_widths_buf);

            // Create new dataspace for attribute
            DataSpace attr3_dataspace = DataSpace( H5S_SCALAR );
            // Create new string datatype for attribute
            strdatatype = StrType( PredType::C_S1, 6 ); // of length 256 characters
            // Set up write buffer for attribute
            const H5std_string fourier_buf( "false" );
            // Create attribute and write to it
            Attribute fourier_attr = metadataset->createAttribute("fourier", strdatatype, attr3_dataspace);
            fourier_attr.write(strdatatype, fourier_buf);
            
            
            // Create new dataspace for attribute
            DataSpace attr4_dataspace = DataSpace( H5S_SCALAR );
            // Create new string datatype for attribute
            strdatatype = StrType( PredType::C_S1, 11 ); // of length 256 characters    
            // Set up write buffer for attribute
            const H5std_string ftype_buf( "[\"counts\"]" );    
            // Create attribute and write to it
            Attribute ftype_attr = metadataset->createAttribute("ftype", strdatatype, attr4_dataspace);
            ftype_attr.write(strdatatype, ftype_buf);
            
            
            delete metadataspace;
            delete metadataset;
            
            delete group;
            delete metagroup;
            delete file;
            
            cerr << "Wrote delta counts map (" << dcMap.Npartpix() << " + " << dcMask.Npartpix() << " entries) to file" << endl;
        }
        // mags:
        if( use_mags ){
            // write out the maps (hdf5)
            string outfilename = "dm_maps.h5";
            H5File *file = new H5File( outfilename, H5F_ACC_TRUNC ); // clobber!
            Group* group = new Group( file->createGroup( "/data" ));
            PredType datatype( PredType::NATIVE_DOUBLE );
            hsize_t dimsf[2];
            dimsf[1] = 2;
            int rank = 2;

            // the map dataset
            double *data = (double*) malloc( sizeof(double) * 2*dmMap.Npartpix()+2 );
            data[0] = dmMap.Order();
            data[1] = RING;
            unsigned int index = 2;
            for( int64 i1=0; i1<dmMap.Npartpix(); ++i1 ){
                int64 i = dmMap.highResPix(i1);
                data[index] = i;
                ++index;
                data[index] = dmMap[i];
                ++index;
            }
            dimsf[0] = dmMap.Npartpix()+1;
            DataSpace *dataspace = new DataSpace( rank, dimsf );
            string datasetname = "/data/map";
            DataSet *dataset = new DataSet( file->createDataSet( datasetname, datatype, *dataspace ) );
            dataset->write( data, PredType::NATIVE_DOUBLE );
            delete dataspace;
            delete dataset;
            free( data );
    
            // the mask dataset
            double *data_mask = (double*) malloc( sizeof(double) * 2*dmMask.Npartpix()+2 );
            data_mask[0] = dmMask.Order();
            data_mask[1] = RING;
            index = 2;
            for( int64 i1=0; i1<dmMask.Npartpix(); ++i1 ){
                int64 i = dmMask.highResPix(i1);
                data_mask[index] = i;
                ++index;
                data_mask[index] = dmMask[i];
                ++index;
            }
            dimsf[0] = dmMask.Npartpix()+1;
            DataSpace *dataspace_mask = new DataSpace( rank, dimsf );
            datasetname = "/data/mask";
            DataSet *dataset_mask = new DataSet( file->createDataSet( datasetname, datatype, *dataspace_mask ) );
            dataset->write( data_mask, PredType::NATIVE_DOUBLE );
            delete dataspace_mask;
            delete dataset_mask;
            free( data_mask );

            // the metadata
            Group* metagroup = new Group( file->createGroup( "/meta" ) );
            DataType metadatatype( H5T_STRING, 1 );
            rank = 1;
            dimsf[0] = 1;
            DataSpace *metadataspace = new DataSpace( rank, dimsf );
            DataSet *metadataset = new DataSet( file->createDataSet( "/meta/meta", metadatatype, *metadataspace ) );

            stringstream ang_mean_line;
            ang_mean_line << "[" << ang_means[0];
            for( unsigned int i=1; i<ang_means.size(); ++i ){
                ang_mean_line << ", " << ang_means[i];
            }
            ang_mean_line << "]";

            // Create new dataspace for attribute
            DataSpace attr1_dataspace = DataSpace( H5S_SCALAR );
            // Create new string datatype for attribute
            StrType strdatatype( PredType::C_S1, ang_mean_line.str().size()+1 ); // of length 256 characters
            // Set up write buffer for attribute
            const H5std_string u_mean_buf( ang_mean_line.str().c_str() );
            // Create attribute and write to it
            Attribute u_mean_attr = metadataset->createAttribute("u_mean", strdatatype, attr1_dataspace);
            u_mean_attr.write(strdatatype, u_mean_buf);

            stringstream ang_widths_line;
            ang_widths_line << "[" << ang_widths[0];
            for( unsigned int i=1; i<ang_widths.size(); ++i ){
                ang_widths_line << ", " << ang_widths[i];
            }
            ang_widths_line << "]";

            // Create new dataspace for attribute
            DataSpace attr2_dataspace = DataSpace( H5S_SCALAR );
            // Create new string datatype for attribute
            strdatatype = StrType( PredType::C_S1, ang_widths_line.str().size()+1 ); // of length 256 characters
            // Set up write buffer for attribute
            const H5std_string u_widths_buf( ang_widths_line.str().c_str() );
            // Create attribute and write to it
            Attribute u_width_attr = metadataset->createAttribute("u_width", strdatatype, attr2_dataspace);
            u_width_attr.write(strdatatype, u_widths_buf);

            
            // Create new dataspace for attribute
            DataSpace attr3_dataspace = DataSpace( H5S_SCALAR );
            // Create new string datatype for attribute
            strdatatype = StrType( PredType::C_S1, 6 ); // of length 256 characters
            // Set up write buffer for attribute
            const H5std_string fourier_buf( "false" );
            // Create attribute and write to it
            Attribute fourier_attr = metadataset->createAttribute("fourier", strdatatype, attr3_dataspace);
            fourier_attr.write(strdatatype, fourier_buf);
            
            
            // Create new dataspace for attribute
            DataSpace attr4_dataspace = DataSpace( H5S_SCALAR );
            // Create new string datatype for attribute
            strdatatype = StrType( PredType::C_S1, 8 ); // of length 256 characters    
            // Set up write buffer for attribute
            const H5std_string ftype_buf( "[\"mag\"]" );    
            // Create attribute and write to it
            Attribute ftype_attr = metadataset->createAttribute("ftype", strdatatype, attr4_dataspace);
            ftype_attr.write(strdatatype, ftype_buf);
            
            delete metadataspace;
            delete metadataset;
    
            delete group;
            delete metagroup;
            delete file;
            
            
            cerr << "Wrote delta mag map (" << dmMap.Npartpix() << " + " << dmMask.Npartpix() << " entries) to file" << endl;

        }
        
        
    } // for z bins

    return 0;
}
