#include "pgen/pvar.h"  // includes Rcpp

#include <stdexcept>

using namespace std;

RPvar::RPvar() {
  PreinitMinimalPvar(&_mp);
}

void RPvar::Load(const string &filename) {
  char errbuf[plink2::kPglErrstrBufBlen];
  plink2::PglErr reterr = LoadMinimalPvar(filename.c_str(), &_mp, errbuf);
  if (reterr != plink2::kPglRetSuccess) {
    if (reterr == plink2::kPglRetNomem) {
      // stop("Out of memory");
      throw logic_error("Out of memory");
    } else if (reterr == plink2::kPglRetReadFail) {
      // stop("File read failure");
      throw logic_error("File read failure");
    } else {
      // stop(&errbuf[7]);
      throw logic_error(&errbuf[7]);
    }
  }
}

uint32_t RPvar::GetVariantCt() const {
  return _mp.variant_ct;
}

const char* RPvar::GetVariantId(uint32_t variant_idx) const {
  if (variant_idx >= _mp.variant_ct) {
    char errbuf[256];
    if (_mp.variant_ct) {
      snprintf(errbuf, 256, "variant_num out of range (%d; must be 1..%d)", variant_idx + 1, _mp.variant_ct);
    } else {
      strcpy(errbuf, "pvar closed");
    }
    // stop(errbuf);
    throw logic_error(errbuf);
  }
  return _mp.variant_ids[variant_idx];
}

std::pair<std::multimap<const char*, int, classcomp>::iterator, std::multimap<const char*, int, classcomp>::iterator> RPvar::GetVariantsById(const char* id) {
  if (_nameToIdxs.empty()) {
    const uint32_t len = _mp.variant_ct;
    for (uint32_t variant_idx = 0; variant_idx != len; ++variant_idx) {
      _nameToIdxs.insert(std::pair<const char*, int>(_mp.variant_ids[variant_idx], variant_idx));
    }
  }
  return _nameToIdxs.equal_range(id);
}

uint32_t RPvar::GetAlleleCt(uint32_t variant_idx) const {
  if (variant_idx >= _mp.variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _mp.variant_ct);
    // stop(errstr_buf);
    throw logic_error(errstr_buf);
  }
  if (!_mp.allele_idx_offsetsp) {
    return 2;
  }
  const uintptr_t* allele_idx_offsets = _mp.allele_idx_offsetsp->p;
  return allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
}

const char* RPvar::GetAlleleCode(uint32_t variant_idx, uint32_t allele_idx) const {
  if (variant_idx >= _mp.variant_ct) {
    char errbuf[256];
    if (_mp.variant_ct) {
      snprintf(errbuf, 256, "variant_num out of range (%d; must be 1..%d)", variant_idx + 1, _mp.variant_ct);
    } else {
      strcpy(errbuf, "pvar closed");
    }
    // stop(errbuf);
    throw logic_error(errbuf);
  }
  uintptr_t allele_idx_offset_base = 2 * variant_idx;
  uint32_t allele_ct = 2;
  if (_mp.allele_idx_offsetsp) {
    const uintptr_t* allele_idx_offsets = _mp.allele_idx_offsetsp->p;
    allele_idx_offset_base = allele_idx_offsets[variant_idx];
    allele_ct = allele_idx_offsets[variant_idx + 1] - allele_idx_offset_base;
  }
  if (allele_idx >= allele_ct) {
    char errbuf[256];
    snprintf(errbuf, 256, "allele_num out of range (%d; must be 1..%d)", allele_idx + 1, allele_ct);
    // stop(errbuf);
    throw logic_error(errbuf);
  }
  return _mp.allele_storage[allele_idx_offset_base + allele_idx];
}

plink2::RefcountedWptr* RPvar::GetAlleleIdxOffsetsp() {
  if (_mp.allele_idx_offsetsp) {
    _mp.allele_idx_offsetsp->ref_ct += 1;
  }
  return _mp.allele_idx_offsetsp;
}

uint32_t RPvar::GetMaxAlleleCt() const {
  return _mp.max_allele_ct;
}

void RPvar::Close() {
  _nameToIdxs.clear();
  plink2::CleanupMinimalPvar(&_mp);
}

RPvar::~RPvar() {
  _nameToIdxs.clear();
  plink2::CleanupMinimalPvar(&_mp);
}
