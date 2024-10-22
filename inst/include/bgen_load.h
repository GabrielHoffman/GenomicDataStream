/** Adapted from load.cpp in rbgen v1.1.5
 * 
 * Changes: Remove dependency on Rcpp
 */

#include <sstream>
#include <map>
#include <set>
#include "genfile/bgen/View.hpp"
#include "genfile/bgen/IndexQuery.hpp"

using namespace std;

namespace {
	
	template< typename T >
	string atoi( const T &value ) {
	        std::ostringstream stream ;
	        stream << value ;
	        return stream.str() ;
	}

	struct DataSetter {
		DataSetter(
			vector<int>* ploidy,
			vector<int> const& ploidy_dimension,
			vector<double> * data,
			vector<int> const& data_dimension,
			vector<bool>* phased,
			size_t variant_i,
			map<size_t, size_t> const& requested_samples
		):
			m_ploidy( ploidy ),
			m_ploidy_dimension( ploidy_dimension ),
			m_data( data ),
			m_data_dimension( data_dimension ),
			m_phased( phased ),
			m_variant_i( variant_i ),
			m_requested_samples( requested_samples ),
			m_requested_sample_i( m_requested_samples.begin() ),
			m_storage_i( 0 ),
			m_order_type( genfile::eUnknownOrderType )
		{
			assert( m_data_dimension[0] == m_ploidy_dimension[0] ) ;
			assert( m_data_dimension[1] == m_ploidy_dimension[1] ) ;
			assert( m_data_dimension[1] == m_requested_samples.size() ) ;
			assert( m_variant_i < m_data_dimension[0] ) ;
			assert( m_data_dimension[2] >= 3 ) ;
			assert( m_phased->size() == m_data_dimension[0] ) ;
		}
	
		// Called once allowing us to set storage.
		void initialise( size_t number_of_samples, size_t number_of_alleles ) {
		}

		// If present with this signature, called once after initialise()
		// to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
		// This enables us to set up storage for the data ahead of time.
		void set_min_max_ploidy(
			genfile::bgen::uint32_t min_ploidy, genfile::bgen::uint32_t max_ploidy,
			genfile::bgen::uint32_t min_entries, genfile::bgen::uint32_t max_entries
		) {
			if( max_entries > m_data_dimension[2] ) {
				throw std::invalid_argument( "max_entries=" + atoi( max_entries )
					+ " (expected at most " + atoi( m_data_dimension[2] ) + ")" ) ;
			}
		}

		// Called once per sample to determine whether we want data for this sample
		bool set_sample( size_t i ) {
			if( m_requested_sample_i != m_requested_samples.end() && m_requested_sample_i->first == i ) {
				m_storage_i = m_requested_sample_i->second ;
				++m_requested_sample_i ;
// #if DEBUG
// 				std::cerr << "DataSetter::set_sample(): sample " << i << " has storage index " << m_storage_i << ".\n" ;
// #endif
				return true ;
			} else {
				// Don't want this sample
				return false ;
			}
		}

		// Called once per sample to set the number of probabilities that are present.
		void set_number_of_entries(
			size_t ploidy,
			size_t number_of_entries,
			genfile::OrderType order_type,
			genfile::ValueType value_type
		) {
			if( value_type != genfile::eProbability ) {
				throw std::invalid_argument(
					"value_type ("
					+ atoi( value_type ) + ", expected "
					+ atoi( genfile::eProbability ) + "=genfile::eProbability)"
				) ;
			}
			if( m_order_type == genfile::eUnknownOrderType ) {
				m_order_type = order_type ;
				// (*m_phased)( m_variant_i ) = ( m_order_type == genfile::ePerPhasedHaplotypePerAllele ) ;
				(*m_phased)[ m_variant_i ] = ( m_order_type == genfile::ePerPhasedHaplotypePerAllele ) ;
			} else {
				assert( order_type == m_order_type ) ;
			}
			int const index = m_variant_i + m_storage_i * m_ploidy_dimension[0] ;
			(*m_ploidy)[ index ] = ploidy ;
		}

		void set_value( genfile::bgen::uint32_t entry_i, double value ) {

			int const index = m_variant_i + m_storage_i * m_data_dimension[0] + entry_i * m_data_dimension[0] * m_data_dimension[1] ;
// #if DEBUG
// 			std::cerr << "Setting data for index " << m_variant_i << ", " << m_storage_i << ", " << entry_i << ": index " << index << "...\n" << std::flush ;
// #endif
			(*m_data)[ index ] = value ;
		}

		void set_value( genfile::bgen::uint32_t entry_i, genfile::MissingValue value ) {
			int const index = m_variant_i + m_storage_i * m_data_dimension[0] + entry_i * m_data_dimension[0] * m_data_dimension[1] ;
// #if DEBUG
// 			std::cerr << "Setting data for index " << m_variant_i << ", " << m_storage_i << ", " << entry_i << ": index " << index << "...\n" << std::flush ;
// #endif
			(*m_data)[ index ] = numeric_limits<double>::quiet_NaN() ;
		}

		void finalise() {
		// nothing to do
		}

	private:
		vector<int>* m_ploidy ;
		vector<int> const m_ploidy_dimension ;
		vector<double> * m_data ;
		vector<int> const m_data_dimension ;
		vector<bool> * m_phased ;

		size_t const m_variant_i ;

		map<size_t, size_t> const& m_requested_samples ;
		map<size_t, size_t>::const_iterator m_requested_sample_i ;
		size_t m_storage_i ;

		genfile::OrderType m_order_type ;
	} ;

	struct set_sample_names {
		typedef map< size_t, size_t > SampleIndexMap ;
		
		set_sample_names( vector<string>* result, SampleIndexMap* sample_indices ):
			m_result( result ),
			m_sample_indices( sample_indices ),
			m_index(0)
		{
			assert( result != 0 ) ;
			assert( sample_indices != 0 ) ;
			assert( sample_indices->size() == result->size() ) ;
		}

		void operator()( string const& value ) {
			m_sample_indices->insert( std::make_pair( m_index, m_index ) ) ;
			(*m_result)[m_index++] = value ;
		}
	private:
		vector<string>* m_result ;
		SampleIndexMap* m_sample_indices ;
		size_t m_index ;
	} ;
	
	struct set_requested_sample_names {
		typedef map<string, size_t> RequestedSamples ;
		typedef map<size_t, size_t> SampleIndexMap ;
		
		set_requested_sample_names(
			vector<string>* result,
			SampleIndexMap* sample_indices,
			RequestedSamples const& requested_samples
		):
			m_result( result ),
			m_sample_indices( sample_indices ),
			m_requested_samples( requested_samples ),
			m_index(0),
			m_value_index(0)
		{
			assert( result != 0 ) ;
			assert( sample_indices != 0 ) ;
			assert( sample_indices->size() == 0 ) ;
			assert( result->size() == requested_samples.size() ) ;
		}
		
		void operator()(string const& value ) {
			RequestedSamples::const_iterator where = m_requested_samples.find( value ) ;
			if( where != m_requested_samples.end() ) {
				(*m_result)[ where->second ] = value ;
				m_sample_indices->insert( make_pair( m_value_index, where->second ) ) ;
			}
			++m_value_index ;
		}
	private:
		vector<string>* m_result ;
		SampleIndexMap* m_sample_indices ;
		RequestedSamples const& m_requested_samples ;
		size_t m_index ;
		size_t m_value_index ;
	} ;
}

void get_all_samples(
	genfile::bgen::View const& view,
	size_t* number_of_samples,
	vector<string>* sampleNames,
	map<size_t, size_t>* requestedSamplesByIndexInDataIndex
) {
	*number_of_samples = view.number_of_samples() ;
	sampleNames->resize( *number_of_samples ) ;
	view.get_sample_ids( set_sample_names( sampleNames, requestedSamplesByIndexInDataIndex ) ) ;
}

void get_requested_samples(
	genfile::bgen::View const& view,
	vector<string> const& requestedSamples,
	size_t* number_of_samples,
	vector<string>* sampleNames,
	map<size_t, size_t>* requestedSamplesByIndexInDataIndex
) {
	// convert requested sample IDs to a map of requested indices.
	map<string, size_t> requestedSamplesByName ;
	for( size_t i = 0; i < requestedSamples.size(); ++i ) {
		requestedSamplesByName.insert( std::map< string, size_t >::value_type( string( requestedSamples[i] ), i )) ;
	}
	if( requestedSamplesByName.size() != requestedSamples.size() ) {
		throw std::invalid_argument(
			"load_unsafe(): requiredSamples: expected a list of unique samples with no repeats."
		) ;
	}

	*number_of_samples = requestedSamples.size() ;
	sampleNames->resize( requestedSamples.size() ) ;
	view.get_sample_ids( set_requested_sample_names( sampleNames, requestedSamplesByIndexInDataIndex, requestedSamplesByName ) ) ;

	// Check each requested sample has been asked for exactly once
	// We count distinct samples, among those requested, that we've found in the data
	// And we also count the min and max index of those samples.
	// If min = 0 and max = (#requested samples-1) and each sample was unique, we're ok.
	set<size_t> checkSamples ;
	size_t minIndex = std::numeric_limits< size_t >::max() ;
	size_t maxIndex = 0 ;
	
	for(
		map<size_t, size_t>::const_iterator p = requestedSamplesByIndexInDataIndex->begin();
		p != requestedSamplesByIndexInDataIndex->end();
		++p
	 ) {
		checkSamples.insert( p->second ) ;
		minIndex = std::min( minIndex, p->second ) ;
		maxIndex = std::max( maxIndex, p->second ) ;
	}
	if( checkSamples.size() != requestedSamples.size() || minIndex != 0 || maxIndex != (requestedSamples.size()-1) ) {
		// Huh.  To be most useful, let's print diagnostics
		// std::cerr << "!! Uh-oh: requested sample indices (data, request) are:\n" ;
		for(
			map<size_t, size_t>::const_iterator p = requestedSamplesByIndexInDataIndex->begin();
			p != requestedSamplesByIndexInDataIndex->end();
			++p
		 ) {
			// std::cerr << p->first << ", " << p->second << ".\n" ;
		}
		
		throw std::invalid_argument(
			"load_unsafe(): requiredSamples contains a sample not present in the data, or data contains a repeated sample ID."
		) ;
	}
}


