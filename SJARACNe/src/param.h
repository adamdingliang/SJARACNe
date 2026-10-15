//------------------------------------------------------------------------------------
// Copyright (C) 2003, 2004  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017
//------------------------------------------------------------------------------------

#ifndef PARAM_H__
#define PARAM_H__

#include <string>
#include <vector>

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

   std::string verbose, infile, outfile, adjfile, hub;
   std::string subnetfile, annotfile, controlId, condition, home_dir;
   std::string subsampleSpec, samplingMethod;
   std::string nullModelFile, thresholdMethod;
   std::string nullModelFormat, nullModelKernelSchema, nullModelEstimator;
   std::string nullModelTailModel, nullModelGeneratorSha256;
   std::string nullModelCalibratorSchema, nullModelCalibratorSha256;
   std::string nullModelFitValuesSha256, nullModelValidationValuesSha256;
   int nullModelM, nullModelNparLimit;
   double nullModelSupportedPMin, nullModelSupportedPMax;
   double nullModelValidatedPMin, nullModelValidatedPMax;
   bool nullModelHasValidatedPMin, nullModelTailExtrapolated;

   std::vector<std::string> subnet, tf_list;

   Parameter()
      : threshold(default_threshold), pvalue(default_pvalue), eps(default_eps),
        sigma(default_sigma), sample(default_sample), percent(default_percent),
        mean(default_mean), cv(default_cv), correction(default_correction),
        nparLimit(default_nparLimit), seed(default_seed), samplingPopulation(0),
        samplingSize(0), thresholdSpecified(false), verbose("off"), infile(""), outfile(""), adjfile(""),
        hub(""), subnetfile(""), annotfile(""), controlId(""), condition(""),
        home_dir("./"), subsampleSpec(""), samplingMethod(""), nullModelFile(""),
        thresholdMethod("none"), nullModelFormat(""), nullModelKernelSchema(""),
        nullModelEstimator(""), nullModelTailModel(""), nullModelGeneratorSha256(""),
        nullModelCalibratorSchema(""), nullModelCalibratorSha256(""),
        nullModelFitValuesSha256(""), nullModelValidationValuesSha256(""),
        nullModelM(0), nullModelNparLimit(0), nullModelSupportedPMin(0.0),
        nullModelSupportedPMax(0.0), nullModelValidatedPMin(0.0),
        nullModelValidatedPMax(0.0),
        nullModelHasValidatedPMin(false), nullModelTailExtrapolated(false),
        subnet(), tf_list() { }
};

//------------------------------------------------------------------------------------

bool equalIgnoreCase(std::string a, std::string b);
void checkParameter(Parameter &p);
void displayParameter(Parameter &p);
void createOutfileName(Parameter &p );
int resolveSubsampleSize(const std::string& spec, int populationSize);

#endif
