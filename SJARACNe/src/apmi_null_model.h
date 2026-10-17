// Built-in estimator-matched empirical PCHIP cutoff surface.
#ifndef APMI_NULL_MODEL_H__
#define APMI_NULL_MODEL_H__

namespace Apmi
{
// Opt-in cutoff representation, not the legacy-affine default. The bundled
// calibration uses unique rank permutations and npar=40. No extrapolation or
// cutoff clipping is performed. Below p=0.0005 requires explicit consent to
// the less-supported extension, although the numerical table reaches p=1e-7.
// Errors use std::string, following the native command-line convention.
double pchipCutoff(int m, double p, int nparLimit = 40,
                   bool allowRareTail = false);
const char* pchipPrimarySha256();
const char* pchipDataSha256();
const char* pchipSummarySha256();
}

#endif
