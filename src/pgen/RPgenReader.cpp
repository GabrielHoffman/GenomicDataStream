/* Adapted from pgenlibr v0.5.0 by Christopher Chang
https://github.com/chrchang/plink-ng/tree/master/2.0/pgenlibr/src
*/
#include "pgen/RPgenReader.h"  // includes Rcpp

using namespace std;

RPgenReader::RPgenReader() : _info_ptr(nullptr),
                             _allele_idx_offsetsp(nullptr),
                             _nonref_flagsp(nullptr),
                             _state_ptr(nullptr) {
}

void RPgenReader::Load(const string &filename, RPvar *rp, 
            int raw_sample_ct,
            const vector<int> &sample_subset_1based) {
  if (_info_ptr) {
    Close();
  }
  _info_ptr = static_cast<plink2::PgenFileInfo*>(malloc(sizeof(plink2::PgenFileInfo)));
  if (!_info_ptr) {
    stop("Out of memory");
  }
  plink2::PreinitPgfi(_info_ptr);
  uint32_t cur_sample_ct = UINT32_MAX;
  // if (raw_sample_ct.isNotNull()) {
  //   cur_sample_ct = as<int>(raw_sample_ct.get());
  // }
  cur_sample_ct = raw_sample_ct;
  uint32_t cur_variant_ct = UINT32_MAX;
  const char* fname = filename.c_str();
  plink2::PgenHeaderCtrl header_ctrl;
  uintptr_t pgfi_alloc_cacheline_ct;
  char errstr_buf[plink2::kPglErrstrBufBlen];
  if (plink2::PgfiInitPhase1(fname, nullptr, cur_variant_ct, cur_sample_ct, &header_ctrl, _info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) {
    stop(&(errstr_buf[7]));
  }
  const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
  if (rp != nullptr) {
    // if (pvar.isNotNull()) {

    // List pvarl = as<List>(pvar);
    // if (strcmp_r_c(pvarl[0], "pvar")) {
    //   stop("pvar is not a pvar object");
    // }
    // XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvarl[1]);
    if (rp->GetVariantCt() != raw_variant_ct) {
      stop("pvar and pgen have different variant counts");
    }
    _allele_idx_offsetsp = rp->GetAlleleIdxOffsetsp();
    if (_allele_idx_offsetsp) {
      _info_ptr->allele_idx_offsets = _allele_idx_offsetsp->p;
    }
    _info_ptr->max_allele_ct = rp->GetMaxAlleleCt();
  } else {
    if (header_ctrl & 0x30) {
      // no need to zero-initialize this
      _allele_idx_offsetsp = plink2::CreateRefcountedWptr(raw_variant_ct + 1);
      _info_ptr->allele_idx_offsets = _allele_idx_offsetsp->p;
      // _info_ptr->max_allele_ct updated by PgfiInitPhase2() in this case
    }
    _info_ptr->max_allele_ct = 2;
  }
  if ((header_ctrl & 0xc0) == 0xc0) {
    // todo: load this in pvar, to enable consistency check.  we use a
    // (manually implemented) shared_ptr in preparation for this.
    const uintptr_t raw_variant_ctl = plink2::DivUp(raw_variant_ct, plink2::kBitsPerWord);
    // no need to zero-initialize this
    _nonref_flagsp = plink2::CreateRefcountedWptr(raw_variant_ctl + 1);
    _info_ptr->nonref_flags = _nonref_flagsp->p;
  }
  const uint32_t file_sample_ct = _info_ptr->raw_sample_ct;
  unsigned char* pgfi_alloc = nullptr;
  if (plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
    stop("Out of memory");
  }
  uint32_t max_vrec_width;
  uintptr_t pgr_alloc_cacheline_ct;
  if (plink2::PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width, _info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) {
    if (pgfi_alloc && (!_info_ptr->vrtypes)) {
      plink2::aligned_free(pgfi_alloc);
    }
    stop(&(errstr_buf[7]));
  }
  if ((!_allele_idx_offsetsp) && (_info_ptr->gflags & 4)) {
    // Note that it's safe to be ignorant of multiallelic variants when
    // phase and dosage info aren't present; GetAlleleCt() then always returns
    // 2 when that isn't actually true, and all ALTs are treated as if they
    // were ALT1, but otherwise everything works properly.
    stop("Multiallelic variants and phase/dosage info simultaneously present; pvar required in this case");
  }
  _state_ptr = static_cast<plink2::PgenReader*>(malloc(sizeof(plink2::PgenReader)));
  if (!_state_ptr) {
    stop("Out of memory");
  }
  plink2::PreinitPgr(_state_ptr);
  plink2::PgrSetFreadBuf(nullptr, _state_ptr);
  const uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * plink2::kCacheline;
  const uintptr_t sample_subset_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec;
  const uintptr_t cumulative_popcounts_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerWord * plink2::kInt32PerVec) * plink2::kBytesPerVec;
  const uintptr_t genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kNypsPerVec) * plink2::kBytesPerVec;
  const uint32_t is_not_plink1_bed = (_info_ptr->vrtypes != nullptr);
  uintptr_t raregeno_byte_ct = 0;
  uintptr_t difflist_sample_ids_byte_ct = 0;
  if (is_not_plink1_bed) {
    const uint32_t max_stored_single_difflist_len = file_sample_ct / plink2::kPglMaxDifflistLenDivisor;
    raregeno_byte_ct = plink2::DivUp(2 * max_stored_single_difflist_len, plink2::kNypsPerVec) * plink2::kBytesPerVec;
    difflist_sample_ids_byte_ct = plink2::RoundUpPow2(3 * max_stored_single_difflist_len * sizeof(int32_t), plink2::kBytesPerVec);
  }
  const uintptr_t ac_byte_ct = plink2::RoundUpPow2(file_sample_ct * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
  const uintptr_t ac2_byte_ct = plink2::RoundUpPow2(file_sample_ct * 2 * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
  uintptr_t multiallelic_hc_byte_ct = 0;
  if (_info_ptr->max_allele_ct != 2) {
    multiallelic_hc_byte_ct = 2 * sample_subset_byte_ct + ac_byte_ct + ac2_byte_ct;
  }
  const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;
  unsigned char* pgr_alloc;
  if (plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglNypTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct + raregeno_byte_ct + difflist_sample_ids_byte_ct + multiallelic_hc_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes + 4 * (plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8), &pgr_alloc)) {
    throw logic_error("Out of memory");
  }
  plink2::PglErr reterr = plink2::PgrInit(fname, max_vrec_width, _info_ptr, _state_ptr, pgr_alloc);
  if (reterr != plink2::kPglRetSuccess) {
    if (!plink2::PgrGetFreadBuf(_state_ptr)) {
      plink2::aligned_free(pgr_alloc);
    }
    snprintf(errstr_buf, plink2::kPglErrstrBufBlen, "PgrInit() error %d", static_cast<int>(reterr));
    throw logic_error(errstr_buf);
  }
  unsigned char* pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct]);
  _subset_include_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _subset_include_interleaved_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);

#ifdef USE_AVX2
  _subset_include_interleaved_vec[-3] = 0;
  _subset_include_interleaved_vec[-2] = 0;
#endif
  _subset_include_interleaved_vec[-1] = 0;

  _subset_cumulative_popcounts = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
  _pgv.genovec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);
  if (is_not_plink1_bed) {
    _raregeno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[raregeno_byte_ct]);
    _difflist_sample_ids_buf = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[difflist_sample_ids_byte_ct]);
  } else {
    _raregeno_buf = nullptr;
    _difflist_sample_ids_buf = nullptr;
  }
  if (multiallelic_hc_byte_ct) {
    _pgv.patch_01_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.patch_01_vals = reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[ac_byte_ct]);
    _pgv.patch_10_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.patch_10_vals = reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[ac2_byte_ct]);
  } else {
    _pgv.patch_01_set = nullptr;
    _pgv.patch_01_vals = nullptr;
    _pgv.patch_10_set = nullptr;
    _pgv.patch_10_vals = nullptr;
  }
  _pgv.phasepresent = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _pgv.phaseinfo = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _pgv.dosage_present = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _pgv.dosage_main = reinterpret_cast<uint16_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct]);
  // if (sample_subset_1based.isNotNull()) {
  //   SetSampleSubsetInternal(sample_subset_1based.get());
  // } else {
  //   _subset_size = file_sample_ct;
  // }
  if (sample_subset_1based.size() > 0) {
    SetSampleSubsetInternal(sample_subset_1based);
  } else {
    _subset_size = file_sample_ct;
  }
  _transpose_batch_buf = reinterpret_cast<plink2::VecW*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglBitTransposeBufbytes]);
  _multivar_vmaj_geno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * genovec_byte_ct]);
  _multivar_vmaj_phasepresent_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
  _multivar_vmaj_phaseinfo_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
  _multivar_smaj_geno_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 4]);
  _multivar_smaj_phaseinfo_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
  _multivar_smaj_phasepresent_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  // pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
}

uint32_t RPgenReader::GetRawSampleCt() const {
  if (!_info_ptr) {
    throw logic_error("pgen is closed");
  }
  return _info_ptr->raw_sample_ct;
}

uint32_t RPgenReader::GetSubsetSize() const {
  return _subset_size;
}

uint32_t RPgenReader::GetVariantCt() const {
  if (!_info_ptr) {
    throw logic_error("pgen is closed");
  }
  return _info_ptr->raw_variant_ct;
}

uint32_t RPgenReader::GetAlleleCt(uint32_t variant_idx) const {
  if (!_info_ptr) {
    throw logic_error("pgen is closed");
  }
  if (variant_idx >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    throw logic_error(errstr_buf);
  }
  if (!_allele_idx_offsetsp) {
    return 2;
  }
  const uintptr_t* allele_idx_offsets = _allele_idx_offsetsp->p;
  return allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
}

uint32_t RPgenReader::GetMaxAlleleCt() const {
  if (!_info_ptr) {
    throw logic_error("pgen is closed");
  }
  return _info_ptr->max_allele_ct;
}


static const double kGenoRDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, NA_REAL);


void RPgenReader::ReadList(vector<double> &buf, const vector<int> & variant_subset, bool meanimpute) {
  if (!_info_ptr) {
    throw logic_error("pgen is closed");
  }

  // Check buffer size
  if( buf.capacity() == 0){
    throw logic_error("Buffer capacity is " + to_string(buf.capacity()));
  }

  if( buf.capacity() < variant_subset.size()*_subset_size){     
    throw logic_error("Buffer capacity is " + to_string(buf.capacity()));
  }

  // assume that buf has the correct dimensions
  const uintptr_t vsubset_size = variant_subset.size();
  const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
  double* buf_iter = &buf[0];
  for (uintptr_t col_idx = 0; col_idx != vsubset_size; ++col_idx) {
    uint32_t variant_idx = variant_subset[col_idx] - 1;
    if (static_cast<uint32_t>(variant_idx) >= raw_variant_ct) {
      char errstr_buf[256];
      snprintf(errstr_buf, 256, "variant_subset element out of range (%d; must be 1..%u)", variant_idx + 1, raw_variant_ct);
      throw logic_error(errstr_buf);
    }
    uint32_t dosage_ct;
    plink2::PglErr reterr = plink2::PgrGetD(_subset_include_vec, _subset_index, _subset_size, variant_idx, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
    if (reterr != plink2::kPglRetSuccess) {
      char errstr_buf[256];
      snprintf(errstr_buf, 256, "PgrGetD() error %d", static_cast<int>(reterr));
      throw logic_error(errstr_buf);
    }
    if (!meanimpute) {
      plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, buf_iter);
    } else {
      plink2::ZeroTrailingNyps(_subset_size, _pgv.genovec);
      if (plink2::Dosage16ToDoublesMeanimpute(_pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, buf_iter)) {
        char errstr_buf[256];
        snprintf(errstr_buf, 256, "variant %d has only missing dosages", variant_idx + 1);
        throw logic_error(errstr_buf);
      }
    }
    buf_iter = &(buf_iter[_subset_size]);
  }
}

void RPgenReader::Close() {
  // don't bother propagating file close errors for now
  if (_info_ptr) {
    CondReleaseRefcountedWptr(&_allele_idx_offsetsp);
    CondReleaseRefcountedWptr(&_nonref_flagsp);
    if (_info_ptr->vrtypes) {
      plink2::aligned_free(_info_ptr->vrtypes);
    }
    plink2::PglErr reterr = plink2::kPglRetSuccess;
    plink2::CleanupPgfi(_info_ptr, &reterr);
    free(_info_ptr);
    _info_ptr = nullptr;
  }
  if (_state_ptr) {
    plink2::PglErr reterr = plink2::kPglRetSuccess;
    plink2::CleanupPgr(_state_ptr, &reterr);
    if (PgrGetFreadBuf(_state_ptr)) {
      plink2::aligned_free(PgrGetFreadBuf(_state_ptr));
    }
    free(_state_ptr);
    _state_ptr = nullptr;
  }
  _subset_size = 0;
}

void RPgenReader::SetSampleSubsetInternal(const vector<int> & sample_subset_1based) {
  const uint32_t raw_sample_ct = _info_ptr->raw_sample_ct;
  const uint32_t raw_sample_ctv = plink2::DivUp(raw_sample_ct, plink2::kBitsPerVec);
  const uint32_t raw_sample_ctaw = raw_sample_ctv * plink2::kWordsPerVec;
  uintptr_t* sample_include = _subset_include_vec;
  plink2::ZeroWArr(raw_sample_ctaw, sample_include);
  const uint32_t subset_size = sample_subset_1based.size();
  if (subset_size == 0) {
    throw logic_error("Empty sample_subset is not currently permitted");
  }
  uint32_t sample_uidx = sample_subset_1based[0] - 1;
  uint32_t idx = 0;
  uint32_t next_uidx;
  while (1) {
    if (sample_uidx >= raw_sample_ct) {
      char errstr_buf[256];
      snprintf(errstr_buf, 256, "sample number out of range (%d; must be 1..%u)", static_cast<int>(sample_uidx + 1), raw_sample_ct);
      throw logic_error(errstr_buf);
    }
    plink2::SetBit(sample_uidx, sample_include);
    if (++idx == subset_size) {
      break;
    }
    next_uidx = sample_subset_1based[idx] - 1;

    // prohibit this since it implies that the caller expects genotypes to be
    // returned in a different order
    if (next_uidx <= sample_uidx) {
      throw logic_error("sample_subset is not in strictly increasing order");
    }
    sample_uidx = next_uidx;
  }
  plink2::FillInterleavedMaskVec(sample_include, raw_sample_ctv, _subset_include_interleaved_vec);
  const uint32_t raw_sample_ctl = plink2::DivUp(raw_sample_ct, plink2::kBitsPerWord);
  plink2::FillCumulativePopcounts(sample_include, raw_sample_ctl, _subset_cumulative_popcounts);
  plink2::PgrSetSampleSubsetIndex(_subset_cumulative_popcounts, _state_ptr, &_subset_index);
  _subset_size = subset_size;
}

void RPgenReader::ReadMaybeSparseHardcallsInternal(int variant_idx, int max_difflist_len, uint32_t* difflist_common_geno_ptr, uint32_t* difflist_len_ptr) {
  // Fills {_raregeno_buf, _difflist_sample_ids_buf, *difflist_common_geno_ptr,
  // *difflist_len_ptr} iff hardcalls are either (i) stored as a simple
  // difflist no longer than max_difflist_len, or (ii) stored as a list of
  // differences from an earlier variant, which is itself stored as a simple
  // difflist, and the sum of the two list lengths is <= max_difflist_len.
  //
  // Otherwise, fills _pgv.genovec and sets *difflist_common_geno_ptr to
  // UINT32_MAX.
  if (!_info_ptr) {
    throw logic_error("pgen is closed");
  }
  const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
  if (static_cast<uint32_t>(variant_idx) >= raw_variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    // stop(errstr_buf);
    throw logic_error(errstr_buf);
  }
  plink2::PglErr reterr = plink2::PgrGetDifflistOrGenovec(_subset_include_vec, _subset_index, _subset_size, max_difflist_len, variant_idx, _state_ptr, _pgv.genovec, difflist_common_geno_ptr, _raregeno_buf, _difflist_sample_ids_buf, difflist_len_ptr);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "PgrGetDifflistOrGenovec() error %d", static_cast<int>(reterr));
    // stop(errstr_buf);
    throw logic_error(errstr_buf);
  }
}

void RPgenReader::ReadAllelesPhasedInternal(int variant_idx) {
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    // stop(errstr_buf);
    throw logic_error(errstr_buf);
  }
  plink2::PglErr reterr = plink2::PgrGetMP(_subset_include_vec, _subset_index, _subset_size, variant_idx, _state_ptr, &_pgv);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "PgrGetMP() error %d", static_cast<int>(reterr));
    // stop(errstr_buf);
    throw logic_error(errstr_buf);
  }
}

RPgenReader::~RPgenReader() {
  Close();
}