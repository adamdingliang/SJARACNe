//------------------------------------------------------------------------------------
// Copyright (C) 2003, 2004  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017
//------------------------------------------------------------------------------------

#ifndef PARAM_H__
#define PARAM_H__

#include <string>
#include <vector>
#include <cstddef>

//------------------------------------------------------------------------------------

struct Parameter
{
   static const double default_threshold, default_pvalue, default_eps, default_sigma;
   static const int    default_sample;
   static const double default_percent, default_mean, default_cv, default_correction;
   static const int    default_nparLimit;
   static const int    default_seed;

   double threshold;  // mi threshold
   double pvalue;     // p-value for mi threshold
   double eps;        // DPI tolerance
   double sigma;      // gaussian kernel width
   int    sample;     // sample number
   double percent;    // high/low percentage for the conditional analysis
   double mean;       // filter mean
   double cv;         // filter coefficient of variance
   double correction; // coorection for noise
   int    nparLimit;  // maximum allowed value of npar
   int    seed;       // seed
   int    samplingPopulation; // observations eligible for resampling
   int    samplingSize;       // observations selected for MI calculation
   bool   thresholdSpecified; // whether -t was supplied, including -t 0
   bool   epsSpecified;       // whether fixed -e was explicitly supplied
   bool   dynamicDpi;         // experimental hub-size-dependent DPI tolerance
   bool   dynamicDpiParametersSpecified;
   int    dpiHref;            // reference count of reconstructed source rows
   double dpiSlope;           // epsilon increase per doubling above dpiHref
   double dpiMax;             // upper bound of experimental dynamic epsilon
   std::size_t dpiActualH;    // unique reconstructed rows, including empty rows

   std::string verbose, infile, outfile, adjfile, hub;
   std::string dpiWitnessFile;
   std::string subnetfile, annotfile, controlId, condition, home_dir;
   std::string subsampleSpec, samplingMethod;
   std::string nullModelFile, thresholdMethod;
   std::string pchipPrimarySha256, pchipDataSha256, pchipSummarySha256;
   int pchipActualM;

   std::vector<std::string> subnet, tf_list;

   Parameter()
      : threshold(default_threshold), pvalue(default_pvalue), eps(default_eps),
        sigma(default_sigma), sample(default_sample), percent(default_percent),
        mean(default_mean), cv(default_cv), correction(default_correction),
        nparLimit(default_nparLimit), seed(default_seed), samplingPopulation(0),
        samplingSize(0), thresholdSpecified(false), epsSpecified(false),
        dynamicDpi(false), dynamicDpiParametersSpecified(false), dpiHref(479),
        dpiSlope(0.05), dpiMax(0.30), dpiActualH(0),
        verbose("off"), infile(""), outfile(""), adjfile(""),
        hub(""), dpiWitnessFile(""), subnetfile(""), annotfile(""), controlId(""), condition(""),
        home_dir("./"), subsampleSpec(""), samplingMethod(""), nullModelFile(""),
        thresholdMethod("none"), pchipPrimarySha256(""), pchipDataSha256(""),
        pchipSummarySha256(""), pchipActualM(0),
        subnet(), tf_list() { }
};

//------------------------------------------------------------------------------------

bool equalIgnoreCase(std::string a, std::string b);
void checkParameter(Parameter &p);
void displayParameter(Parameter &p);
void createOutfileName(Parameter &p );
int resolveSubsampleSize(const std::string& spec, int populationSize);
double resolveDynamicDpiEpsilon(std::size_t sourceRows, int referenceHubs,
                                double slope, double maximum);

#endif
