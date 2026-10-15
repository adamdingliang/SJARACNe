//------------------------------------------------------------------------------------
// Strict loader and evaluator for estimator-matched AP-MI null tail models.
//------------------------------------------------------------------------------------

#include "apmi_null_model.h"

#include <cerrno>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

namespace
{
const char *MODEL_FORMAT = "sjaracne-apmi-gpd-tail-v1";
const char *KERNEL_SCHEMA = "sjaracne-apmi-v1";
const char *ESTIMATOR = "sjaracne-adaptive-partitioning";
const char *SAMPLING_NULL = "independent-uniform-rank-permutation";
const char *RANK_POLICY = "unique-ordinal-ranks";
const char *TAIL_MODEL = "generalized-pareto-mle-floc0";
const char *RNG = "mt19937-rejection-fisher-yates-v1";
const char *CALIBRATION_STATUS = "accepted";
const char *CALIBRATOR_SCHEMA = "sjaracne-apmi-gpd-calibrator-v1";
const char *VALIDATION_METHOD = "independent-rank-permutation-stream";

std::string trim(const std::string& value)
{
   const std::string whitespace(" \t\n\r\f\v");
   const std::string::size_type first = value.find_first_not_of(whitespace);
   if (first == std::string::npos)
      return "";
   const std::string::size_type last = value.find_last_not_of(whitespace);
   return value.substr(first, last - first + 1);
}

std::string context(const std::string& filename, int lineNumber)
{
   std::ostringstream out;
   out << filename;
   if (lineNumber > 0)
      out << ":" << lineNumber;
   return out.str();
}

std::map<std::string, std::string> readFields(const std::string& filename)
{
   std::ifstream input(filename.c_str());
   if (!input.is_open())
      throw std::string("Unable to open AP-MI null model: ") + filename;

   const char *names[] = {
      "format", "kernel_schema", "estimator", "sampling_null", "rank_policy",
      "m", "npar_limit", "tail_model", "tail_threshold_quantile",
      "tail_threshold", "tail_probability", "tail_shape", "tail_scale",
      "tail_endpoint", "supported_p_min", "supported_p_max", "validated_p_min",
      "validated_p_max", "calibration_status", "calibrator_schema",
      "calibrator_sha256", "validation_method", "stability_probability",
      "stability_relative_range", "stability_relative_tolerance",
      "validation_family_confidence", "validation_point_confidence",
      "default_p", "default_p_cutoff", "fit_draws", "validation_draws",
      "fit_seed", "validation_seed", "rng", "generator_sha256",
      "fit_values_sha256", "validation_values_sha256", "scipy_version"
   };
   const std::set<std::string> allowed(names, names + sizeof(names) / sizeof(names[0]));
   std::map<std::string, std::string> fields;
   std::string line;
   int lineNumber = 0;

   while (std::getline(input, line))
   {
      lineNumber++;
      line = trim(line);
      if (line.empty() || line[0] == '#')
         continue;

      const std::string::size_type separator = line.find('=');
      if (separator == std::string::npos || line.find('=', separator + 1) != std::string::npos)
         throw std::string("Malformed AP-MI null-model record at ") +
               context(filename, lineNumber) + "; expected exactly one key=value pair.";

      const std::string key = trim(line.substr(0, separator));
      const std::string value = trim(line.substr(separator + 1));
      if (key.empty() || value.empty())
         throw std::string("Empty key or value in AP-MI null model at ") +
               context(filename, lineNumber) + ".";
      if (allowed.find(key) == allowed.end())
         throw std::string("Unknown AP-MI null-model field '") + key + "' at " +
               context(filename, lineNumber) + ".";
      if (!fields.insert(std::make_pair(key, value)).second)
         throw std::string("Duplicate AP-MI null-model field '") + key + "' at " +
               context(filename, lineNumber) + ".";
   }

   for (std::set<std::string>::const_iterator name = allowed.begin();
        name != allowed.end(); ++name)
      if (fields.find(*name) == fields.end())
         throw std::string("Missing required AP-MI null-model field '") + *name +
               "' in " + filename + ".";

   return fields;
}

double parseDouble(const std::map<std::string, std::string>& fields,
                   const std::string& key, const std::string& filename)
{
   const std::string& text = fields.find(key)->second;
   errno = 0;
   char *end = NULL;
   const double value = std::strtod(text.c_str(), &end);
   if (errno == ERANGE || end == text.c_str() || *end != '\0' || !std::isfinite(value))
      throw std::string("AP-MI null-model field '") + key +
            "' must be a finite number in " + filename + ".";
   return value;
}

long parseLong(const std::map<std::string, std::string>& fields,
               const std::string& key, long lower, long upper,
               const std::string& filename)
{
   const std::string& text = fields.find(key)->second;
   errno = 0;
   char *end = NULL;
   const long value = std::strtol(text.c_str(), &end, 10);
   if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
       value < lower || value > upper)
      throw std::string("AP-MI null-model field '") + key +
            "' is outside its valid integer range in " + filename + ".";
   return value;
}

unsigned long parseUint32(const std::map<std::string, std::string>& fields,
                          const std::string& key, const std::string& filename)
{
   const std::string& text = fields.find(key)->second;
   errno = 0;
   char *end = NULL;
   const unsigned long value = std::strtoul(text.c_str(), &end, 10);
   if (text.empty() || text[0] == '-' || errno == ERANGE || end == text.c_str() ||
       *end != '\0' || value > 0xffffffffUL)
      throw std::string("AP-MI null-model field '") + key +
            "' is outside the uint32 range in " + filename + ".";
   return value;
}

void requireEqual(const std::map<std::string, std::string>& fields,
                  const std::string& key, const std::string& expected,
                  const std::string& filename)
{
   const std::string& actual = fields.find(key)->second;
   if (actual != expected)
      throw std::string("Unsupported AP-MI null-model ") + key + " '" + actual +
            "' in " + filename + "; expected '" + expected + "'.";
}

bool isSha256(const std::string& value)
{
   if (value.length() != 64)
      return false;
   for (std::string::const_iterator c = value.begin(); c != value.end(); ++c)
      if (!(std::isdigit(static_cast<unsigned char>(*c)) ||
            (*c >= 'a' && *c <= 'f') || (*c >= 'A' && *c <= 'F')))
         return false;
   return true;
}

bool approximatelyEqual(double a, double b)
{
   const double scale = std::fmax(1.0, std::fmax(std::fabs(a), std::fabs(b)));
   return std::fabs(a - b) <= 1e-10 * scale;
}
}

ApmiNullModel ApmiNullModel::load(const std::string& filename)
{
   const std::map<std::string, std::string> fields = readFields(filename);
   requireEqual(fields, "format", MODEL_FORMAT, filename);
   requireEqual(fields, "kernel_schema", KERNEL_SCHEMA, filename);
   requireEqual(fields, "estimator", ESTIMATOR, filename);
   requireEqual(fields, "sampling_null", SAMPLING_NULL, filename);
   requireEqual(fields, "rank_policy", RANK_POLICY, filename);
   requireEqual(fields, "tail_model", TAIL_MODEL, filename);
   requireEqual(fields, "rng", RNG, filename);
   requireEqual(fields, "calibration_status", CALIBRATION_STATUS, filename);
   requireEqual(fields, "calibrator_schema", CALIBRATOR_SCHEMA, filename);
   requireEqual(fields, "validation_method", VALIDATION_METHOD, filename);

   ApmiNullModel model;
   model.format = fields.find("format")->second;
   model.kernelSchema = fields.find("kernel_schema")->second;
   model.estimator = fields.find("estimator")->second;
   model.samplingNull = fields.find("sampling_null")->second;
   model.rankPolicy = fields.find("rank_policy")->second;
   model.tailModel = fields.find("tail_model")->second;
   model.rng = fields.find("rng")->second;
   model.generatorSha256 = fields.find("generator_sha256")->second;
   model.fitValuesSha256 = fields.find("fit_values_sha256")->second;
   model.validationValuesSha256 = fields.find("validation_values_sha256")->second;
   model.scipyVersion = fields.find("scipy_version")->second;
   model.calibrationStatus = fields.find("calibration_status")->second;
   model.calibratorSchema = fields.find("calibrator_schema")->second;
   model.calibratorSha256 = fields.find("calibrator_sha256")->second;
   model.validationMethod = fields.find("validation_method")->second;
   model.m = static_cast<int>(parseLong(fields, "m", 2, INT_MAX, filename));
   model.nparLimit = static_cast<int>(parseLong(fields, "npar_limit", 1, INT_MAX, filename));
   model.fitDraws = static_cast<int>(parseLong(fields, "fit_draws", 1, INT_MAX, filename));
   model.validationDraws = static_cast<int>(parseLong(fields, "validation_draws", 1, INT_MAX, filename));
   model.fitSeed = parseUint32(fields, "fit_seed", filename);
   model.validationSeed = parseUint32(fields, "validation_seed", filename);
   model.tailThresholdQuantile = parseDouble(fields, "tail_threshold_quantile", filename);
   model.tailThreshold = parseDouble(fields, "tail_threshold", filename);
   model.tailProbability = parseDouble(fields, "tail_probability", filename);
   model.tailShape = parseDouble(fields, "tail_shape", filename);
   model.tailScale = parseDouble(fields, "tail_scale", filename);
   model.supportedPMin = parseDouble(fields, "supported_p_min", filename);
   model.supportedPMax = parseDouble(fields, "supported_p_max", filename);
   model.defaultP = parseDouble(fields, "default_p", filename);
   model.defaultPCutoff = parseDouble(fields, "default_p_cutoff", filename);
   model.stabilityProbability = parseDouble(fields, "stability_probability", filename);
   model.stabilityRelativeRange = parseDouble(fields, "stability_relative_range", filename);
   model.stabilityRelativeTolerance = parseDouble(fields, "stability_relative_tolerance", filename);
   model.validationFamilyConfidence = parseDouble(fields, "validation_family_confidence", filename);
   model.validationPointConfidence = parseDouble(fields, "validation_point_confidence", filename);

   const std::string endpoint = fields.find("tail_endpoint")->second;
   model.hasTailEndpoint = endpoint != "none";
   model.tailEndpoint = model.hasTailEndpoint
      ? parseDouble(fields, "tail_endpoint", filename) : 0.0;
   const std::string validated = fields.find("validated_p_min")->second;
   model.hasValidatedPMin = validated != "none";
   model.validatedPMin = model.hasValidatedPMin
      ? parseDouble(fields, "validated_p_min", filename) : 0.0;
   model.validatedPMax = parseDouble(fields, "validated_p_max", filename);

   if (!(model.tailThresholdQuantile > 0.0 && model.tailThresholdQuantile < 1.0))
      throw std::string("AP-MI null-model tail_threshold_quantile must be within (0,1) in ") + filename + ".";
   if (model.tailThreshold < 0.0)
      throw std::string("AP-MI null-model tail_threshold must be nonnegative in ") + filename + ".";
   const double theoreticalMaximum = std::log(static_cast<double>(model.m));
   if (model.tailThreshold > theoreticalMaximum)
      throw std::string("AP-MI null-model tail_threshold exceeds the theoretical AP-MI maximum log(m) in ") + filename + ".";
   if (!(model.tailProbability > 0.0 && model.tailProbability <= 1.0))
      throw std::string("AP-MI null-model tail_probability must be within (0,1] in ") + filename + ".";
   if (!(model.tailScale > 0.0))
      throw std::string("AP-MI null-model tail_scale must be positive in ") + filename + ".";
   if (!(model.stabilityProbability > 0.0 &&
         model.stabilityProbability <= model.tailProbability &&
         approximatelyEqual(model.stabilityProbability, model.defaultP)))
      throw std::string("AP-MI null-model stability_probability must equal default_p and lie within the fitted tail in ") + filename + ".";
   if (!(model.stabilityRelativeRange >= 0.0 &&
         model.stabilityRelativeRange <= model.stabilityRelativeTolerance &&
         model.stabilityRelativeTolerance > 0.0 &&
         model.stabilityRelativeTolerance < 1.0))
      throw std::string("AP-MI null-model stability diagnostics are not accepted in ") + filename + ".";
   if (!(model.validationFamilyConfidence > 0.0 &&
         model.validationFamilyConfidence < 1.0 &&
         model.validationPointConfidence >= model.validationFamilyConfidence &&
         model.validationPointConfidence < 1.0))
      throw std::string("AP-MI null-model validation confidence values are invalid in ") + filename + ".";
   if (!(model.supportedPMin > 0.0 && model.supportedPMin <= model.supportedPMax &&
         model.supportedPMax <= model.tailProbability))
      throw std::string("AP-MI null-model supported p range is invalid in ") + filename + ".";
   if (!approximatelyEqual(model.supportedPMin, model.defaultP) ||
       model.defaultP > model.supportedPMax)
      throw std::string("AP-MI null-model supported_p_min must equal default_p and lie within its supported p range in ") + filename + ".";
   if (!model.hasValidatedPMin)
      throw std::string("Accepted AP-MI null model must record a held-out validated_p_min in ") + filename + ".";
   if (model.validatedPMin < model.supportedPMin ||
       model.validatedPMin > model.validatedPMax ||
       model.validatedPMax > model.supportedPMax)
      throw std::string("AP-MI null-model validated_p_min lies outside its supported p range in ") + filename + ".";
   if (model.fitSeed == model.validationSeed)
      throw std::string("AP-MI null-model fit_seed and validation_seed must differ in ") + filename + ".";

   if (model.tailShape < 0.0)
   {
      const double expectedEndpoint = model.tailThreshold - model.tailScale / model.tailShape;
      if (!model.hasTailEndpoint || model.tailEndpoint <= model.tailThreshold ||
          !approximatelyEqual(model.tailEndpoint, expectedEndpoint))
         throw std::string("AP-MI null-model tail_endpoint is inconsistent with its negative-shape GPD in ") + filename + ".";
      if (model.tailEndpoint > theoreticalMaximum)
         throw std::string("AP-MI null-model tail_endpoint exceeds the theoretical AP-MI maximum log(m) in ") + filename + ".";
   }
   else if (model.hasTailEndpoint)
      throw std::string("AP-MI null-model tail_endpoint must be 'none' for a nonnegative-shape GPD in ") + filename + ".";

   if (!isSha256(model.calibratorSha256) || !isSha256(model.generatorSha256) ||
       !isSha256(model.fitValuesSha256) ||
       !isSha256(model.validationValuesSha256))
      throw std::string("AP-MI null-model SHA-256 provenance fields must contain 64 hexadecimal characters in ") + filename + ".";
   if (model.scipyVersion.empty())
      throw std::string("AP-MI null-model scipy_version must not be empty in ") + filename + ".";

   const double calculatedDefault = model.cutoff(model.defaultP);
   if (!approximatelyEqual(calculatedDefault, model.defaultPCutoff))
      throw std::string("AP-MI null-model default_p_cutoff is inconsistent with its GPD parameters in ") + filename + ".";
   if (model.defaultPCutoff > theoreticalMaximum)
      throw std::string("AP-MI null-model default_p_cutoff exceeds the theoretical AP-MI maximum log(m) in ") + filename + ".";

   return model;
}

double ApmiNullModel::cutoff(double pvalue) const
{
   if (!std::isfinite(pvalue) || pvalue < supportedPMin || pvalue > supportedPMax)
   {
      std::ostringstream out;
      out << std::setprecision(17)
          << "Requested MI p-value " << pvalue
          << " is outside the AP-MI null model's supported range ["
          << supportedPMin << ", " << supportedPMax << "].";
      throw out.str();
   }

   const double logRatio = std::log(pvalue / tailProbability);
   double excess;
   if (std::fabs(tailShape) < 1e-10)
      excess = -tailScale * logRatio;
   else
      excess = tailScale / tailShape * std::expm1(-tailShape * logRatio);

   double value = tailThreshold + excess;
   if (!std::isfinite(value) || value < tailThreshold)
      throw std::string("AP-MI null model produced a non-finite or sub-threshold cutoff.");
   if (hasTailEndpoint && value > tailEndpoint)
   {
      if (approximatelyEqual(value, tailEndpoint))
         value = tailEndpoint;
      else
         throw std::string("AP-MI null model produced a cutoff beyond its finite GPD endpoint.");
   }
   if (value > std::log(static_cast<double>(m)))
      throw std::string("AP-MI null model produced a cutoff beyond the theoretical AP-MI maximum log(m).");
   return value;
}
