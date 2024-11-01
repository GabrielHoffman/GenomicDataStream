
#include "pgenstream.h"
#include "pgen/RPgenReader.h"
#include "pgen/pvar.h"

using namespace Rcpp;


//' Loads variant IDs and allele codes from a .pvar or .bim file (which can be
//' compressed with gzip or Zstd).
//'
//' @param filename .pvar/.bim file path.
//' @return A pvar object, which can be queried for variant IDs and allele
//' codes.
//' @export
// [[Rcpp::export]]
SEXP NewPvar(String filename) {
  XPtr<class RPvar> pvar(new RPvar(), true);
  pvar->Load(filename);
  return List::create(_["class"] = "pvar", _["pvar"] = pvar);
}

//' Convert variant index to variant ID string.
//'
//' @param pvar Object returned by NewPvar().
//' @param variant_num Variant index (1-based).
//' @return The variant_numth variant ID string.
//' @export
// [[Rcpp::export]]
String GetVariantId(List pvar, int variant_num) {
  if (strcmp_r_c(pvar[0], "pvar")) {
    stop("pvar is not a pvar object");
  }
  XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar[1]);
  String ss(rp->GetVariantId(variant_num - 1));
  return ss;
}

//' Convert variant ID string to variant index(es).
//'
//' @param pvar Object returned by NewPvar().
//' @param id Variant ID to look up.
//' @return A list of all (1-based) variant indices with the given variant ID.
//' @export
// [[Rcpp::export]]
IntegerVector GetVariantsById(List pvar, String id) {
  if (strcmp_r_c(pvar[0], "pvar")) {
    stop("pvar is not a pvar object");
  }
  XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar[1]);
  std::pair<std::multimap<const char*, int, classcomp>::iterator, std::multimap<const char*, int, classcomp>::iterator> equal_range = rp->GetVariantsById(id.get_cstring());
  std::multimap<const char*, int, classcomp>::iterator i1 = equal_range.first;
  std::multimap<const char*, int, classcomp>::iterator i2 = equal_range.second;
  const uint32_t len = std::distance(i1, i2);
  IntegerVector iv = IntegerVector(len);
  for (uint32_t uii = 0; uii != len; ++uii) {
    iv[uii] = i1->second + 1;
    ++i1;
  }
  return iv;
}

//' Look up an allele code.
//'
//' @param pvar Object returned by NewPvar().
//' @param variant_num Variant index (1-based).
//' @param allele_num Allele index (1-based).
//' @return The allele_numth allele code for the variant_numth variant.
//' allele_num=1 corresponds to the REF allele, allele_num=2 corresponds to the
//' first ALT allele, allele_num=3 corresponds to the second ALT allele if it
//' exists and errors out otherwise, etc.
//' @export
// [[Rcpp::export]]
String GetAlleleCode(List pvar, int variant_num, int allele_num) {
  if (strcmp_r_c(pvar[0], "pvar")) {
    stop("pvar is not a pvar object");
  }
  XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar[1]);
  String ss(rp->GetAlleleCode(variant_num - 1, allele_num - 1));
  return ss;
}

//' Closes a pvar object, releasing memory.
//'
//' @param pvar Object returned by NewPvar().
//' @return No return value, called for side-effect.
//' @export
// [[Rcpp::export]]
void ClosePvar(List pvar) {
  if (strcmp_r_c(pvar[0], "pvar")) {
    stop("pvar is not a pvar object");
  }
  XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar[1]);
  rp->Close();
}




//' Opens a .pgen or PLINK 1 .bed file.
//'
//' @param filename .pgen/.bed file path.
//' @param pvar Object (see NewPvar()) corresponding to the .pgen's companion
//' .pvar; technically optional, but necessary for some functionality.  In
//' particular, at multiallelic variants, all ALT alleles may be collapsed
//' together when .pvar information is not available.
//' @param raw_sample_ct Number of samples in file; required if it's a PLINK 1
//' .bed file, otherwise optional.
//' @param sample_subset List of 1-based positions of samples to load;
//' optional, all samples are loaded if this is not specified.
//' @return A pgen object, which can be queried for genotype/dosage data.
//' @export
// [[Rcpp::export]]
SEXP NewPgen(String filename, Nullable<List> pvar = R_NilValue,
             Nullable<int> raw_sample_ct = R_NilValue,
             Nullable<IntegerVector> sample_subset = R_NilValue) {
  XPtr<class RPgenReader> pgen(new RPgenReader(), true);
  pgen->Load(filename, pvar, raw_sample_ct, sample_subset);
  return List::create(_["class"] = "pgen", _["pgen"] = pgen);
}

//' Returns the number of samples in the file.
//'
//' @param pgen Object returned by NewPgen().
//' @return Number of samples.
//' @export
// [[Rcpp::export]]
int GetRawSampleCt(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return rp->GetRawSampleCt();
}

//' Returns the number of variants in the file.
//'
//' @param pvar_or_pgen Object returned by NewPvar() or NewPgen().
//' @return Number of variants.
//' @export
// [[Rcpp::export]]
int GetVariantCt(List pvar_or_pgen) {
  const char* c_str = as<String>(pvar_or_pgen[0]).get_cstring();
  if (!strcmp(c_str, "pvar")) {
    XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar_or_pgen[1]);
    return rp->GetVariantCt();
  } else if (!strcmp(c_str, "pgen")) {
    XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pvar_or_pgen[1]);
    return rp->GetVariantCt();
  }
  stop("pvar_or_pgen is not a pvar or pgen object");
}

// Throughout this codebase, _num indicates R 1-based indexing.

//' Returns the effective number of alleles for a variant.  Note that if no
//' pvar was provided to the NewPgen() call, this function may return 2 even at
//' multiallelic variants, since the .pgen may not store allele-count
//' information.
//'
//' @param pvar_or_pgen Object returned by NewPvar() or NewPgen().
//' @param variant_num Variant index (1-based).
//' @return max(2, <number of alleles the variant_numth variant is known to
//' have>).  Note that if no
//' @export
// [[Rcpp::export]]
int GetAlleleCt(List pvar_or_pgen, int variant_num) {
  const char* c_str = as<String>(pvar_or_pgen[0]).get_cstring();
  const uint32_t variant_idx = variant_num - 1;
  if (!strcmp(c_str, "pvar")) {
    XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar_or_pgen[1]);
    return rp->GetAlleleCt(variant_idx);
  } else if (!strcmp(c_str, "pgen")) {
    XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pvar_or_pgen[1]);
    return rp->GetAlleleCt(variant_idx);
  }
  stop("pvar_or_pgen is not a pvar or pgen object");
}

//' Returns the maximum GetAlleleCt() value across all variants in the file.
//'
//' @param pvar_or_pgen Object returned by NewPvar() or NewPgen().
//' @return Maximum GetAlleleCt() value across all variants.
//' @export
// [[Rcpp::export]]
int GetMaxAlleleCt(List pvar_or_pgen) {
  const char* c_str = as<String>(pvar_or_pgen[0]).get_cstring();
  if (!strcmp(c_str, "pvar")) {
    XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar_or_pgen[1]);
    return rp->GetMaxAlleleCt();
  } else if (!strcmp(c_str, "pgen")) {
    XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pvar_or_pgen[1]);
    return rp->GetMaxAlleleCt();
  }
  stop("pvar_or_pgen is not a pvar or pgen object");
}

//' Returns whether explicitly phased hardcalls are present.
//'
//' @param pgen Object returned by NewPgen().
//' @return TRUE if the file contains at least one phased heterozygous
//' hardcall, FALSE otherwise.
//' @export
// [[Rcpp::export]]
bool HardcallPhasePresent(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return rp->HardcallPhasePresent();
}

//' Returns a numeric buffer that Read() or ReadHardcalls() can load to.
//'
//' @param pgen Object returned by NewPgen().
//' @return Numeric vector with appropriate length for Read() and
//' ReadHardcalls().
//' @export
// [[Rcpp::export]]
NumericVector Buf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return NumericVector(rp->GetSubsetSize());
}

//' Returns an empty two-row numeric matrix that ReadAlleles() can load to.
//'
//' @param pgen Object returned by NewPgen().
//' @return Numeric matrix with two rows, and appropriate number of columns for
//' ReadAlleles().
//' @export
// [[Rcpp::export]]
NumericVector AlleleCodeBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return NumericMatrix(2, rp->GetSubsetSize());
}

//' Returns an integer buffer that ReadHardcalls() can load to.
//'
//' @param pgen Object returned by NewPgen().
//' @return Integer vector with appropriate length for ReadHardcalls().
//' @export
// [[Rcpp::export]]
IntegerVector IntBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return IntegerVector(rp->GetSubsetSize());
}

//' Returns an empty two-row integer matrix that ReadAlleles() can load to.
//'
//' @param pgen Object returned by NewPgen().
//' @return Integer matrix with two rows, and appropriate number of columns for
//' ReadAlleles().
//' @export
// [[Rcpp::export]]
IntegerVector IntAlleleCodeBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return IntegerMatrix(2, rp->GetSubsetSize());
}

//' Returns a bool buffer that ReadAlleles() can load phasing information to.
//'
//' @param pgen Object returned by NewPgen().
//' @return Logical vector with appropriate length for ReadAlleles().
//' @export
// [[Rcpp::export]]
LogicalVector BoolBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return LogicalVector(rp->GetSubsetSize());
}

//' Loads the variant_numth variant, and then fills buf with \{0, 1, 2, NA\}
//' values indicating the number of copies of the first ALT (or user-specified)
//' allele each sample has.
//'
//' This function treats the data as diploid; you can divide by 2, and then
//' treat 0.5 as NA, if it's actually haploid.
//'
//' @param pgen Object returned by NewPgen().
//' @param buf Buffer returned by Buf() or IntBuf().
//' @param variant_num Variant index (1-based).
//' @param allele_num Allele index; 1 corresponds to REF, 2 to the first ALT
//' allele, 3 to the second ALT allele if it exists, etc.  Optional, defaults
//' 2.
//' @return No return value, called for buf-filling side-effect.
//' @export
// [[Rcpp::export]]
void ReadHardcalls(List pgen, SEXP buf, int variant_num, int allele_num = 2) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  if (Rf_isMatrix(buf)) {
    // otherwise the original buffer is not modified by Read[Int]Hardcalls
    stop("buf must be a non-matrix vector");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  const int variant_idx = variant_num - 1;
  const int allele_idx = allele_num - 1;
  if (TYPEOF(buf) == REALSXP) {
    rp->ReadHardcalls(buf, variant_idx, allele_idx);
  } else if (TYPEOF(buf) == INTSXP) {
    rp->ReadIntHardcalls(buf, variant_idx, allele_idx);
  } else {
    stop("Unsupported buf type");
  }
}

//' Loads the variant_numth variant, and then fills buf with numeric dosages
//' in [0, 2] indicating the dosages of the first ALT (or user-specified)
//' allele for each sample, with missing values represented by NA.
//'
//' This function treats the data as diploid; divide by 2 to obtain haploid
//' dosages.
//'
//' @param pgen Object returned by NewPgen().
//' @param buf Buffer returned by Buf().
//' @param variant_num Variant index (1-based).
//' @param allele_num Allele index; 1 corresponds to REF, 2 to the first ALT
//' allele, 3 to the second ALT allele if it exists, etc.  Optional, defaults
//' 2.
//' @return No return value, called for buf-filling side-effect.
//' @export
// [[Rcpp::export]]
void Read(List pgen, NumericVector buf, int variant_num, int allele_num = 2) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  if (Rf_isMatrix(buf)) {
    stop("buf must be a non-matrix vector");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  rp->Read(buf, variant_num - 1, allele_num - 1);
}

//' Loads the variant_numth variant, and then fills acbuf with integer allele
//' codes, where each column of the buffer corresponds to a sample.  An allele
//' code of 0 corresponds to the REF allele, 1 to the first ALT, 2 to the
//' second ALT, etc.  Missing hardcalls are represented by a pair of NA codes.
//'
//' This function treats the data as diploid.  If it's really haploid, you may
//' want to compare the two rows, and then treat samples where the allele codes
//' differ as missing values.
//'
//' @param pgen Object returned by NewPgen().
//' @param acbuf Buffer returned by AlleleCodeBuf() or IntAlleleCodeBuf().
//' @param variant_num Variant index (1-based).
//' @param phasepresent_buf Buffer returned by BoolBuf().  Optional; if
//' provided, elements are set to true when the sample has known phase.  Most
//' of these values will be TRUE even when the raw data is unphased, because
//' homozygous genotypes always have known phase.  (Missing genotypes are
//' considered to have unknown phase.)
//' @return No return value, called for acbuf-filling side-effect.
//' @export
// [[Rcpp::export]]
void ReadAlleles(List pgen, SEXP acbuf, int variant_num, Nullable<LogicalVector> phasepresent_buf = R_NilValue) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  const int variant_idx = variant_num - 1;
  // in this case, integer may be a more appropriate default than numeric?
  if (TYPEOF(acbuf) == INTSXP) {
    rp->ReadAlleles(acbuf, phasepresent_buf, variant_idx);
  } else if (TYPEOF(acbuf) == REALSXP) {
    rp->ReadAllelesNumeric(acbuf, phasepresent_buf, variant_idx);
  } else {
    stop("Unsupported acbuf type");
  }
}

//' Load hardcalls for multiple variants as an integer matrix.
//'
//' This function treats the data as diploid; you can divide by 2, and then
//' treat 0.5 as NA, if it's actually haploid.
//'
//' @param pgen Object returned by NewPgen().
//' @param variant_subset Integer vector containing 1-based indexes of variants
//' to load.
//' @return Integer matrix, where rows correspond to samples, columns
//' correspond to variant_subset, and values are in \{0, 1, 2, NA\} indicating
//' the number of hardcall ALT allele copies.  For multiallelic variants, all
//' ALT alleles are combined.
//' @export
// [[Rcpp::export]]
IntegerMatrix ReadIntList(List pgen, IntegerVector variant_subset) {
  // return value: rows = samples, columns = variants (from R's perspective)
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  IntegerMatrix result(rp->GetSubsetSize(), variant_subset.size());
  rp->ReadIntList(result, variant_subset);
  return result;
}

//' Load dosages for multiple variants as a numeric matrix.
//'
//' This function treats the data as diploid; divide by 2 to obtain haploid
//' dosages.
//'
//' @param pgen Object returned by NewPgen().
//' @param variant_subset Integer vector containing 1-based indexes of variants
//' to load.
//' @param meanimpute Optional; if true, missing values are mean-imputed
//' instead of being represented by NA.
//' @return Numeric matrix, where rows correspond to samples, and columns
//' correspond to variant_subset.  Values are in [0, 2] indicating ALT
//' allele dosages, or NA for missing dosages.  For multiallelic variants, all
//' ALT alelles are combined.
//' @export
// [[Rcpp::export]]
NumericMatrix ReadList(List pgen, IntegerVector variant_subset, bool meanimpute = false) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  NumericMatrix result(rp->GetSubsetSize(), variant_subset.size());
  rp->ReadList(result, variant_subset, meanimpute);
  return result;
}

//' Compute variant scores.
//'
//' This function treats the data as diploid; divide by 2 to obtain scores
//' based on a haploid dosage matrix.
//'
//' @param pgen Object returned by NewPgen().
//' @param weights Sample weights.
//' @param variant_subset Integer vector containing 1-based indexes of variants
//' to include in the dosage matrix.  Optional; by default, all variants are
//' included.
//' @return Numeric vector, containing product of sample-weight vector and the
//' specified subset of the dosage matrix.
//' @export
// [[Rcpp::export]]
NumericVector VariantScores(List pgen, NumericVector weights,
                            Nullable<IntegerVector> variant_subset = R_NilValue) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  int variant_ct;
  if (variant_subset.isNotNull()) {
    variant_ct = as<IntegerVector>(variant_subset).size();
  } else {
    variant_ct = rp->GetVariantCt();
  }
  NumericVector result(variant_ct);
  rp->FillVariantScores(result, weights, variant_subset);
  return result;
}

//' Closes a pgen object, releasing resources.
//'
//' @param pgen Object returned by NewPgen().
//' @return No return value, called for side-effect.
//' @export
// [[Rcpp::export]]
void ClosePgen(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  rp->Close();
}
