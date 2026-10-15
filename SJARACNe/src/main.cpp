//------------------------------------------------------------------------------------
// Copyright (C) 2003  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017; this version assumes adaptive_partitioning
//------------------------------------------------------------------------------------

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include "apmi_null_model.h"
#include "matrix.h"
#include "parseargs.h"

//------------------------------------------------------------------------------------

const int NUM_OPTIONS = 22;

const char *option[NUM_OPTIONS] =
{
"-a <algorithm>     default: adaptive_partitioning",
"-c <+/-probeId %>  Conditional network reconstruction, default: NONE [****]\n"
"                   [format: \"+24 0.35\", \"-1973_s_at 0.4\"]",
"-e <tolerance>     DPI tolerance, default: 1",
"-f <mean> <cv>     Gene filter by the mean and coefficient of variance (cv) of\n"
"                   the expression values, default: mean=0, cv=0",
"-H <ARACNE_HOME>   Directory containing ARACNE configuration files,\n"
"                   default: current working directory",
"-h <probeId>       Hub gene (only MI w/ hub gene will be computed),\n"
"                   default: NONE",
"--h | --help       Display this help and exit",
"-i <file>          Input gene expression profile dataset (required)",
"-j <file>          Existing adjacency matrix (.adj) file",
"-k <kernel_width>  Gaussian kernel width (accurate method only),\n"
"                   default: determined by program",
"-l <file>          File containing a list of probes annotated as transcription\n"
"                   factors in the input dataset, default: NONE [***]",
"-M <file>          Estimator-matched AP-MI GPD-tail null model. The model must\n"
"                   exactly match the selected observation count and -N limit",
"-N <npar_limit>    Maximum allowed value of npar, default: 20",
"-S <Seed>	    Initial seed for random number generator, default: 1",
"-n <level>         Array measurement noise level, default: 0",
"-o <file>          Output file name (optional) [*]",
"-p <p-value>       P-value for MI threshold (e.g., 1e-7), default: 1 [**]",
"-r <sample_number> Legacy full-size bootstrap with replacement; retained only\n"
"                   for reproducibility, default: 0 (disabled)",
"-s <file>          File containing a list of probes for which a subnetwork will\n"
"                   be constructed, default: NONE",
"-t <threshold>     MI threshold, default: 0",
"-u <count|percent> Fixed-size subsampling without replacement, e.g. 80 or 80%;\n"
"                   default: disabled in the native command",
"-v <verbose>       on|off, default: off"
};

const int NUM_USAGE_NOTES = 4;

const char *usageNotes[NUM_USAGE_NOTES] =
{
"   [*] If no output file is specified by the user, an output will be\n"
"       automatically generated in the same directory as the input file by\n"
"       appending some of the parameter values, such as kernel width, MI\n"
"       threshold, tolerance and so on, at the end of the input file name, and\n"
"       changing the file extension to \".adj\".",
"  [**] If the '-t' option is supplied, it will enforce the program to use the\n"
"       specified MI threshold, therefore the '-p' option will be ignored.\n"
"       Otherwise, the program will automatically determine the MI threshold\n"
"       given the p-value. The default, p-value=1, will preserve all pairwise MI.",
" [***] This option is ideal for transcriptional network reconstruction. If\n"
"       provided, DPI will not remove any connection of a transcription factor (TF)\n"
"       by connections between two probes not annotated as TFs. This option is\n"
"       often used in conjunction with '-s', which specifies a list of probes that\n"
"       are either the same or a subset of the probes specified by '-l'.",
"[****] Conditional network reconstructs the network given a specified probe being\n"
"       most expressed or least expressed. In the format that follows '-c', probeId\n"
"       indicates the probe to be conditioned on; '+' or '-' specify whether the\n"
"       upper or lower tail of the probe's expression should be used as the\n"
"       condition, and '%' is a percentage between (0,1) specifying the proportion\n"
"       of samples used as the conditioning subset. Example usage: \"-c +24 0.35\",\n"
"       \"-c -1973_s_at 0.4\"."
};

//------------------------------------------------------------------------------------
// writeUsage() writes usage to stdout

void writeUsage(const char *progname)
{
   std::printf("Usage: %s [OPTION] ...\n\n", progname);

   for (int i = 0; i < NUM_OPTIONS; i++)
      std::printf("    %s\n", option[i]);

   std::printf("\n");

   for (int i = 0; i < NUM_USAGE_NOTES; i++)
      std::printf("%s\n", usageNotes[i]);
}

//------------------------------------------------------------------------------------

double parseFiniteNumericOption(const char *optionName, const char *argument)
{
   const std::string text = argument == NULL ? "" : argument;
   if (text.empty())
      throw std::string("Option '") + optionName + "' requires a finite numeric value.";

   errno = 0;
   char *end = NULL;
   const double value = std::strtod(text.c_str(), &end);
   if (errno == ERANGE || end == text.c_str() || *end != '\0' || !std::isfinite(value))
      throw std::string("Option '") + optionName + "' requires a finite numeric value, not '" +
            text + "'.";
   return value;
}

//------------------------------------------------------------------------------------

Parameter parseParameter(int argc, char *argv[])
{
   Parameter p; // initializes parameters to default values

   std::string temp;

   // parse arguments using macros defined in parseargs.h
   ARGBEGIN
   {
      case 'a': temp         = ARGF(); break;            // algorithm
      case 'c': temp         = ARGF();                   // condition
                p.condition  = temp.substr(0,1);
                p.controlId  = temp.substr(1);
                p.percent    = std::atof(ARGF());
            break;
      case 'e': p.eps        = std::atof(ARGF()); break; // DPI tolerance
      case 'f': p.mean       = std::atof(ARGF());        // mean
                p.cv         = std::atof(ARGF());        // coefficient of variance
                break;
      case 'H': p.home_dir   = ARGF(); break;            // ARACNE_HOME
      case 'h': p.hub        = ARGF(); break;            // hub gene
      case 'i': p.infile     = ARGF(); break;            // input file
      case 'j': p.adjfile    = ARGF(); break;            // adjacency matrix file
      case 'k': p.sigma      = std::atof(ARGF()); break; // gaussian kernel width
      case 'l': p.annotfile  = ARGF(); break;            // TF annotation file
      case 'M': p.nullModelFile = ARGF();                 // AP-MI null model
                if (p.nullModelFile == "")
                   throw std::string("Option '-M' requires a model file.");
                break;
      case 'N': p.nparLimit  = std::atoi(ARGF()); break; // max npar value
      case 'S': p.seed       = std::atoi(ARGF()); break; // seed
      case 'n': p.correction = std::atof(ARGF()); break; // correction for noise
      case 'o': p.outfile    = ARGF(); break;            // output file
      case 'p': p.pvalue     = std::atof(ARGF()); break; // p-value
      case 'r': p.sample     = std::atoi(ARGF()); break; // bootstrap sample number
      case 's': p.subnetfile = ARGF();                   // subset of probes
                if (p.subnetfile == "")
                   throw std::string("Option '-s' requires a subnetwork file.");
                break;
      case 't': p.threshold  = parseFiniteNumericOption("-t", ARGF()); // mi threshold
                p.thresholdSpecified = true;
                break;
      case 'u': p.subsampleSpec = ARGF();                // unique subsample
                if (p.subsampleSpec == "")
                   throw std::string("Option '-u' requires an exact count or "
                                     "percentage such as 80%.");
                break;
      case 'v': p.verbose    = ARGF(); break;            // verbose
      default : throw std::string("unknown parameter ") + ARGC();
   }
   ARGEND;

   checkParameter(p);

   std::cout << "Displaying parameters" << std::endl;
   displayParameter(p);

   return p;
}

//------------------------------------------------------------------------------------

void findLegacyThreshold(int n, Parameter& p)
{
   std::string filename = p.home_dir + "config_threshold.txt";

   std::ifstream infile(filename.c_str());
   if (!infile.is_open())
      throw "Unable to open " + filename;

   std::string line;

   std::getline(infile, line);
   while (line.length() > 0 && line[0] == '>')
      std::getline(infile, line);

   std::istringstream sin(line);

   if (!sin.good() || sin.peek() == EOF)
      throw "Configuration file format error: " + filename;

   double alpha, beta, gamma;
   sin >> alpha >> beta >> gamma;
   p.threshold = (alpha - std::log(p.pvalue)) / (-beta - gamma * n);

   infile.close();
}

//------------------------------------------------------------------------------------

void findEstimatorMatchedThreshold(int n, Parameter& p)
{
   const ApmiNullModel model = ApmiNullModel::load(p.nullModelFile);

   if (model.m != n)
   {
      std::ostringstream message;
      message << "AP-MI null model was calibrated for m=" << model.m
              << ", but this run uses exactly m=" << n
              << " observations. Exact-m models are not interpolated.";
      throw message.str();
   }
   if (model.nparLimit != p.nparLimit)
   {
      std::ostringstream message;
      message << "AP-MI null model was calibrated with npar_limit="
              << model.nparLimit << ", but this run uses -N " << p.nparLimit
              << ". The AP-MI estimator settings must match exactly.";
      throw message.str();
   }

   p.threshold = model.cutoff(p.pvalue);
   p.thresholdMethod = "estimator-matched AP-MI permutation-null GPD tail";
   p.nullModelFormat = model.format;
   p.nullModelKernelSchema = model.kernelSchema;
   p.nullModelEstimator = model.estimator;
   p.nullModelTailModel = model.tailModel;
   p.nullModelCalibratorSchema = model.calibratorSchema;
   p.nullModelCalibratorSha256 = model.calibratorSha256;
   p.nullModelGeneratorSha256 = model.generatorSha256;
   p.nullModelFitValuesSha256 = model.fitValuesSha256;
   p.nullModelValidationValuesSha256 = model.validationValuesSha256;
   p.nullModelM = model.m;
   p.nullModelNparLimit = model.nparLimit;
   p.nullModelSupportedPMin = model.supportedPMin;
   p.nullModelSupportedPMax = model.supportedPMax;
   p.nullModelHasValidatedPMin = model.hasValidatedPMin;
   p.nullModelValidatedPMin = model.validatedPMin;
   p.nullModelValidatedPMax = model.validatedPMax;
   p.nullModelTailExtrapolated =
      !model.hasValidatedPMin || p.pvalue < model.validatedPMin ||
      p.pvalue > model.validatedPMax;

   if (model.hasValidatedPMin && p.pvalue < model.validatedPMin)
      std::cout << "[THRESHOLD] WARNING: requested p=" << p.pvalue
                << " is below the model's held-out validated_p_min="
                << model.validatedPMin
                << "; this cutoff is a fitted-tail extrapolation."
                << std::endl;
   else if (model.hasValidatedPMin && p.pvalue > model.validatedPMax)
      std::cout << "[THRESHOLD] WARNING: requested p=" << p.pvalue
                << " is above the model's held-out validated_p_max="
                << model.validatedPMax
                << "; this cutoff lies outside the held-out validated range."
                << std::endl;
   else if (!model.hasValidatedPMin)
      std::cout << "[THRESHOLD] WARNING: this AP-MI null model records no "
                   "held-out validated p range; the requested cutoff is a "
                   "fitted-tail extrapolation."
                << std::endl;
}

//------------------------------------------------------------------------------------

void runStandard(int argc, char *argv[])
{
   Parameter p = parseParameter(argc, argv);
   std::srand(p.seed);
   Microarray_Set data;

   data.read(p.infile);

   if (p.mean > 0.0 || p.cv > 0.0)
      std::cout << data.filter(p.mean, p.cv)
                << " markers disabled due to lack of dynamic range."
                << std::endl << std::endl;

   std::vector<int> lower, upper, selected, *arrays = NULL;

   int controlId = -1;
   int nsample   = data.uarrays.size();

   if (p.controlId != "")
   {
      controlId = data.getProbeId(p.controlId);
      if (controlId == -1)
         throw "Cannot find marker: " + p.controlId;

      data.getHighLowPercent(p.percent, controlId, lower, upper);
      arrays = (p.condition == "+" ? &upper : &lower);
      nsample = arrays->size();
   }

   if (p.subsampleSpec != "")
   {
      p.samplingPopulation = nsample;
      p.samplingSize = resolveSubsampleSize(p.subsampleSpec, nsample);
      p.samplingMethod = "fixed-size without replacement";

      data.sampleWithoutReplacement(selected, p.samplingSize,
                                    static_cast<unsigned int>(p.seed), arrays);
      arrays = &selected;
      nsample = selected.size();

      std::cout << "[SAMPLING] Fixed-size sampling without replacement: selected "
                << nsample << " of " << p.samplingPopulation
                << " eligible observations (request: " << p.subsampleSpec
                << ")." << std::endl;

      if (equalIgnoreCase(p.verbose, "on"))
      {
         std::cout << "[SAMPLING] Selected original observation indices (0-based):";
         for (std::vector<int>::const_iterator id = selected.begin();
              id != selected.end(); ++id)
            std::cout << " " << *id;
         std::cout << std::endl;
      }
   }

   if (nsample < 2)
   {
      std::ostringstream s;
      s << "At least 2 observations are required for MI calculation; found "
        << nsample;

      if (arrays != NULL)
         s << " after conditional selection";

      s << ".";
      throw s.str();
   }

   std::cout << "Marker No: " << data.markerset.size()
             << " (" << data.Get_Num_Active_Markers() << " active)"
             << ", Array No: " << nsample << std::endl;

   if (p.thresholdSpecified)
   {
      p.thresholdMethod = "explicit -t";
      if (!p.nullModelFile.empty())
         std::cout << "[THRESHOLD] '-M' is ignored because an explicit '-t' "
                      "threshold was supplied."
                   << std::endl;
   }
   else if (p.pvalue != 1.0)
   {
      if (!p.nullModelFile.empty())
         findEstimatorMatchedThreshold(nsample, p);
      else
      {
         std::cout << "[THRESHOLD] WARNING: no estimator-matched model was "
                      "supplied with '-M'; using the legacy affine threshold "
                      "calibration for backward compatibility."
                   << std::endl;
         findLegacyThreshold(nsample, p);
         p.thresholdMethod = "legacy affine calibration";
      }

      std::cout << "MI threshold determined for p=" << p.pvalue << ": " << p.threshold
                << std::endl;
   }
   else if (!p.nullModelFile.empty())
      std::cout << "[THRESHOLD] '-M' was supplied but p=1 preserves all MI values; "
                   "the model is not loaded."
                << std::endl;

   std::vector<int> ids;

   if (p.hub != "")
   {
      std::string hub = p.hub;

      if (hub.length() == 0 || hub[0] != '_')
         hub = "_" + hub;

      int hubId = data.getProbeId(hub);
      if (hubId == -1)
         throw "Cannot find the hub probe: " + p.hub + ", nothing to be computed!";

      ids.push_back(hubId);
   }

   int numSubnets = p.subnet.size();
   int numMissingSubnets = 0;

   for (int i = 0; i < numSubnets; i++)
   {
      int gid = data.getProbeId(p.subnet[i]);
      if (gid == -1)
      {
         std::cout << "Cannot find probe: " << p.subnet[i] << " in \""
                   << p.subnetfile << "\" ... ignored." << std::endl;
         numMissingSubnets++;
      }
      else
         ids.push_back(gid);
   }

   if (p.subnetfile != "")
   {
      std::cout << "[SUBNETWORK] Requested: " << numSubnets
                << ", matched: " << ids.size()
                << ", missing: " << numMissingSubnets << std::endl;

      if (ids.size() == 0)
      {
         std::ostringstream s;
         s << "Subnetwork file \"" << p.subnetfile << "\" was supplied, but zero "
           << "probe IDs matched the first column of \"" << p.infile
           << "\" (requested: " << numSubnets
           << "). Refusing to construct an all-gene network.";
         throw s.str();
      }
   }

   if (p.correction != 0.0)
      data.computeMarkerVariance(arrays);

   // Adaptive-partitioning MI does not consume Marker::bandwidth.  Avoid the
   // unnecessary variance pass and per-marker sort performed by that calculation.

   Transfac transfac;

   int numTFs = p.tf_list.size();

   for (int i = 0; i < numTFs; i++)
   {
      int gid = data.getProbeId(p.tf_list[i]);
      if (gid == -1)
         std::cout << "Cannot find probe: " << p.tf_list[i] << " in \""
                   << p.annotfile << "\" ... ignored." << std::endl;
      else
     transfac[gid] = 1;
   }

   Matrix matrix;

   if (p.adjfile != "")
   {
      matrix.read(data, p);
      std::vector<int> missingAdjacencyRows;

      for (std::vector<int>::const_iterator id = ids.begin(); id != ids.end(); ++id)
         if (!matrix.hasAdjacencyRow(*id))
            missingAdjacencyRows.push_back(*id);

      if (!missingAdjacencyRows.empty())
      {
         std::ostringstream s;
         s << "Adjacency file \"" << p.adjfile
           << "\" does not contain a source row for requested hub";

         if (missingAdjacencyRows.size() > 1)
            s << "s";

         s << ": ";

         for (std::vector<int>::size_type i = 0;
              i < missingAdjacencyRows.size(); ++i)
         {
            if (i > 0)
               s << ", ";

            const Marker& marker = data.markerset[missingAdjacencyRows[i]];
            s << marker.accnum.substr(1);
         }

         s << ". Refusing to treat absent adjacency rows as empty networks.";
         throw s.str();
      }
   }
   else
   {
      std::vector<int> bs;

      if (p.sample > 0)
      {
         p.samplingPopulation = nsample;
         p.samplingSize = nsample;
         p.samplingMethod = "legacy bootstrap with replacement";
         std::cout << "[SAMPLING] WARNING: '-r' uses the deprecated legacy "
                      "full-size bootstrap with replacement. Use '-u' for "
                      "fixed-size sampling without replacement."
                   << std::endl;
         data.bootStrap(bs, arrays);
         arrays = &bs;
      }

      data.addNoise();

      data.createEdgeMatrix(nsample, matrix, p.threshold, controlId, p.correction,
                            p.nparLimit, ids, arrays);
   }

   if (p.eps != 1.0)
   {
      if (matrix.hasEnoughSourceRowsForDpi())
      {
         std::cout << "[NETWORK] Applying DPI ..." << std::endl;
         matrix.reduce(p.eps, ids, transfac);
      }
      else
         std::cout << "[NETWORK] Skipping DPI: fewer than two source rows are "
                      "available; no eligible triangle can be formed."
                   << std::endl;
   }

   if (p.outfile == "")
      createOutfileName(p);

   matrix.write(data, ids, p);
}

//------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   bool displayUsage = false;

   if (argc < 2)
      displayUsage = true;
   else
   {
      std::string firstOption = argv[1];

      if (firstOption == "--help" || firstOption == "--h")
         displayUsage = true;
   }

   if (displayUsage)
   {
      writeUsage(argv[0]);
      return 1;
   }

   try
   {
      runStandard(argc, argv);
   }
   catch (const std::string& s)
   {
      std::cerr << argv[0] << ": " << s << std::endl;
      return 1;
   }

   return 0;
}
