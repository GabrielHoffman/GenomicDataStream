
#include <string>

#include <GenomicDataStream.h>
#include <vcfstream.h>
#include <bgenstream.h>
#include <vcfpp.h>

#include <armadillo>

using namespace std;
using namespace GenomicDataStreamLib;
 

// [[Rcpp::export]]
void fastLM( const arma::colvec& y, 
                const std::string &file,
                const std::string &field){

}

int main(int argc, char *argv[]){

    std::cout << "Start program"  << std::endl;
 
 	string file, field;

 	vcfpp::BcfReader *reader = new vcfpp::BcfReader( file );

 	int chunkSize = 4;
    Param param(file, field);

    // Initialise GenomicDataStream with file
    unique_ptr<GenomicDataStream> gdsStream;
    
    switch( getFileType(file) ){
        case VCF:
        case VCFGZ:
        case BCF:
            gdsStream = make_unique<vcfstream>( param );
            break;
        case BGEN:
            gdsStream = make_unique<bgenstream>( param );
            break;
        case OTHER:
            break;
    }


    return 0;
}