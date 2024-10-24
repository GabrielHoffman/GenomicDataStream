
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

	// std::cout << "Start program"  << std::endl;
 
 	// parameters 
	string file = "test.vcf.gz";
	string field = "DS";    // read dosage field
	string region = "";     // no region filter
	string samples = "-";   // no samples filter
	int chunkSize = 4;      // each chunk will read 4 variants

	// initialize parameters
	Param param(file, region, samples, chunkSize);
	param.setField( field );

	// Initialise GenomicDataStream to read 
	// VCF/VCFGZ/BCF and BGEN with same interface
	unique_ptr<GenomicDataStream> gdsStream = createFileView( param );

	// declare DataChunk storing an Armadillo matrix for each chunk
	DataChunk<arma::mat, VariantInfo> chunk;

	// loop through chunks
	while( gdsStream->getNextChunk( chunk ) ){

	    // get data from chunk
	    // chunk.getData();

	    // get variant information
	    // chunk.getInfo();

	    // Do analysis with variants in this chunk
	}


	return 0;
}