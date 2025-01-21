#ifndef __RPGEN_READER_H__
#define __RPGEN_READER_H__

#include <include/pgenlib_ffi_support.h>
#include <include/pgenlib_read.h>
#include "pgen/pvar.h"  

#include <vector>
using namespace std;

class RPgenReader {
public:
  // imitates Python/pgenlib.pyx
  RPgenReader();

#if __cplusplus >= 201103L
  RPgenReader(const RPgenReader&) = delete;
  RPgenReader& operator=(const RPgenReader&) = delete;
#endif

  void Load(const string &filename, RPvar *rp, int raw_sample_ct,
            const vector<int> &sample_subset_1based);

  // uint32_t GetRawSampleCt() const;

  // uint32_t GetSubsetSize() const;

  // uint32_t GetVariantCt() const;

  // uint32_t GetAlleleCt(uint32_t variant_idx) const;

  // uint32_t GetMaxAlleleCt() const;

  // bool HardcallPhasePresent() const;

  // void ReadIntHardcalls(IntegerVector buf, int variant_idx, int allele_idx);

  // void ReadHardcalls(NumericVector buf, int variant_idx, int allele_idx);

  // void Read(NumericVector buf, int variant_idx, int allele_idx);

  // void ReadAlleles(IntegerMatrix acbuf,
  //                  Nullable<LogicalVector> phasepresent_buf, int variant_idx);

  // void ReadAllelesNumeric(NumericMatrix acbuf,
  //                         Nullable<LogicalVector> phasepresent_buf,
  //                         int variant_idx);

  // void ReadIntList(IntegerMatrix buf, IntegerVector variant_subset);

  void ReadList(vector<double> & buf, const vector<int> &variant_subset, bool meanimpute);

  // void FillVariantScores(NumericVector result, NumericVector weights, Nullable<IntegerVector> variant_subset);

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

  plink2::VecW* _transpose_batch_buf;
  // kPglNypTransposeBatch (= 256) variants at a time, and then transpose
  uintptr_t* _multivar_vmaj_geno_buf;
  uintptr_t* _multivar_vmaj_phasepresent_buf;
  uintptr_t* _multivar_vmaj_phaseinfo_buf;
  uintptr_t* _multivar_smaj_geno_batch_buf;
  uintptr_t* _multivar_smaj_phaseinfo_batch_buf;
  uintptr_t* _multivar_smaj_phasepresent_batch_buf;

  void SetSampleSubsetInternal(const vector<int> &sample_subset_1based);

  void ReadAllelesPhasedInternal(int variant_idx);
};

#endif  // __RPGEN_READER_H__