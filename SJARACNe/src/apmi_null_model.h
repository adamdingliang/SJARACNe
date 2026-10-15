//------------------------------------------------------------------------------------
// Estimator-matched permutation-null tail model for adaptive-partitioning MI.
//------------------------------------------------------------------------------------

#ifndef APMI_NULL_MODEL_H__
#define APMI_NULL_MODEL_H__

#include <string>

struct ApmiNullModel
{
   std::string format;
   std::string kernelSchema;
   std::string estimator;
   std::string samplingNull;
   std::string rankPolicy;
   std::string tailModel;
   std::string rng;
   std::string generatorSha256;
   std::string fitValuesSha256;
   std::string validationValuesSha256;
   std::string scipyVersion;
   std::string calibrationStatus;
   std::string calibratorSchema;
   std::string calibratorSha256;
   std::string validationMethod;

   int m;
   int nparLimit;
   int fitDraws;
   int validationDraws;
   unsigned long fitSeed;
   unsigned long validationSeed;

   double tailThresholdQuantile;
   double tailThreshold;
   double tailProbability;
   double tailShape;
   double tailScale;
   bool hasTailEndpoint;
   double tailEndpoint;
   double supportedPMin;
   double supportedPMax;
   bool hasValidatedPMin;
   double validatedPMin;
   double validatedPMax;
   double defaultP;
   double defaultPCutoff;
   double stabilityProbability;
   double stabilityRelativeRange;
   double stabilityRelativeTolerance;
   double validationFamilyConfidence;
   double validationPointConfidence;

   static ApmiNullModel load(const std::string& filename);
   double cutoff(double pvalue) const;
};

#endif
