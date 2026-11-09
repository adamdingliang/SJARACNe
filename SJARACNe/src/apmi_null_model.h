// Built-in estimator-matched empirical PCHIP cutoff surface.
#ifndef APMI_NULL_MODEL_H__
#define APMI_NULL_MODEL_H__

namespace Apmi
{
// Opt-in cutoff representation, not the legacy-affine default. The bundled
// calibration uses unique rank permutations and npar=40. No extrapolation or
// cutoff clipping is performed. The computable domain reaches p=1e-7; this is
// not a claim of uniformly validated rare-tail probability accuracy.
// Errors use std::string, following the native command-line convention.
double pchipCutoff(int m, double p, int nparLimit = 40);
const char* pchipPrimarySha256();
const char* pchipDataSha256();
const char* pchipSummarySha256();
}

#endif
