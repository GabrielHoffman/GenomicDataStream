/* Adapted from pgenlibr v0.5.0 by Christopher Chang
https://github.com/chrchang/plink-ng/tree/master/2.0/pgenlibr/src
*/

#ifndef __RPgenReader_H_
#define __RPgenReader_H_

// #include "include/pgenlib_ffi_support.h"
// #include "include/pgenlib_read.h"

#include "../pgenlibr/include/include/pvar_ffi_support.h"
#include "../pgenlibr/include/include/pgenlib_read.h"
#include "pgen/pvar.h"  // includes Rcpp

using namespace std;

class RPgenReader {
public:
  // similar to Python/pgenlib.pyx ; has a bit more functionality as of Feb
  // 2025
  RPgenReader();

#if __cplusplus >= 201103L
  RPgenReader(const RPgenReader&) = delete;
  RPgenReader& operator=(const RPgenReader&) = delete;
#endif

  void Load(const string &filename, RPvar *rp, 
            int raw_sample_ct,
            const vector<int> &sample_subset_1based);

  uint32_t GetRawSampleCt() const;

  uint32_t GetSubsetSize() const;

  uint32_t GetVariantCt() const;

  uint32_t GetAlleleCt(uint32_t variant_idx) const;

  uint32_t GetMaxAlleleCt() const;

  void ReadList( vector<double> &buf, const vector<int> &variant_subset, bool meanimpute);

  void Close();

  ~RPgenReader();

private:
  plink2::PgenFileInfo* _info_ptr;
  plink2::RefcountedWptr* _allele_idx_offsetsp;
  plink2::RefcountedWptr* _nonref_flagsp;
  plink2::PgenReader* _state_ptr;
  uintptr_t* _subset_include_vec;
  uintptr_t* _subset_include_interleaved_vec;
  uint32_t* _subset_cumulative_popcounts;
  plink2::PgrSampleSubsetIndex _subset_index;
  uint32_t _subset_size;

  plink2::PgenVariant _pgv;

  uintptr_t* _raregeno_buf;
  uint32_t* _difflist_sample_ids_buf;

  plink2::VecW* _transpose_batch_buf;
  // kPglNypTransposeBatch (= 256) variants at a time, and then transpose
  uintptr_t* _multivar_vmaj_geno_buf;
  uintptr_t* _multivar_vmaj_phasepresent_buf;
  uintptr_t* _multivar_vmaj_phaseinfo_buf;
  uintptr_t* _multivar_smaj_geno_batch_buf;
  uintptr_t* _multivar_smaj_phaseinfo_batch_buf;
  uintptr_t* _multivar_smaj_phasepresent_batch_buf;

  void SetSampleSubsetInternal(const vector<int> &sample_subset_1based);

  void ReadMaybeSparseHardcallsInternal(int variant_idx, int max_simple_difflist_len, uint32_t* difflist_common_geno_ptr, uint32_t* difflist_len_ptr);

  void ReadAllelesPhasedInternal(int variant_idx);
};


#endif
