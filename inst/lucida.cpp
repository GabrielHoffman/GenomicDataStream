
#include <string>

#include <GenomicDataStream.h>

#include <armadillo>

using namespace std;
using namespace gds;
 

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
    
    unique_ptr<GenomicDataStream> gdsStream = createFileView( param );



    return 0;
}