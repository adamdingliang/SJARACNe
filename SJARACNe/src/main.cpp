//------------------------------------------------------------------------------------
// Copyright (C) 2003  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017; this version assumes adaptive_partitioning
//------------------------------------------------------------------------------------

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "apmi.h"
#include "apmi_null_model.h"
#include "matrix.h"
#include "parseargs.h"

//------------------------------------------------------------------------------------

const int NUM_OPTIONS = 27;

const char *option[NUM_OPTIONS] =
{
"-a <algorithm>     default: adaptive_partitioning",
"-c <+/-probeId %>  Conditional network reconstruction, default: NONE [****]\n"
"                   [format: \"+24 0.35\", \"-1973_s_at 0.4\"]",
"-e <tolerance>     DPI tolerance, default: 1",
"--dynamic-dpi      Opt in to experimental hub-size-dependent DPI; requires -s,\n"
"                   excludes fixed -e, adjacency replay -j and conditional -c",
"--dpi-href <count> Dynamic DPI reference source rows, positive integer; default: 479",
"--dpi-slope <s>    Dynamic DPI increase per doubling of H above href;\n"
"                   finite positive value, default: 0.05",
"--dpi-max <value>  Dynamic DPI epsilon cap, finite in [0,1); default: 0.30",
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
"-M pchip           Opt in to the bundled experimental AP-MI PCHIP cutoff;\n"
"                   requires -N 40, m=50..10000 and p=1e-7..0.05; no tail flag",
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
"-v <verbose>       on|off, default: off",
"-W <file>          Benchmark diagnostic: per-source counts of pre-DPI edges\n"
"                   with at least 1, 2, 3, 5, or 10 qualifying witnesses"
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

int parseDynamicDpiReference(const std::string& argument)
{
   if (argument.empty())
      throw std::string("Option '--dpi-href' requires a positive integer.");
   for (std::string::const_iterator c = argument.begin(); c != argument.end(); ++c)
      if (*c < '0' || *c > '9')
         throw std::string("Option '--dpi-href' requires a positive integer, not '") + argument + "'.";
   errno = 0;
   char *end = NULL;
   const long value = std::strtol(argument.c_str(), &end, 10);
   if (errno == ERANGE || *end != '\0' || value < 1 || value > INT_MAX)
      throw std::string("Option '--dpi-href' must be an integer in [1,2147483647].");
   return static_cast<int>(value);
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
      case 'e': p.eps        = parseFiniteNumericOption("-e", ARGF()); // DPI tolerance
                p.epsSpecified = true;
                break;
      case '-':
      {
         // Long options are parsed only at an option boundary. ARGF still
         // consumes short-option values normally, even when a filename begins
         // with '--', preserving the historical single-character interface.
         const std::string longArgument = _arg;
         _arg = (char *)"";
         const std::size_t equal = longArgument.find('=');
         const std::string name = longArgument.substr(0, equal);
         const bool inlineValue = equal != std::string::npos;
         if (name == "dynamic-dpi")
         {
            if (inlineValue)
               throw std::string("Option '--dynamic-dpi' is a flag and takes no value.");
            p.dynamicDpi = true;
         }
         else if (name == "dpi-href" || name == "dpi-slope" || name == "dpi-max")
         {
            const std::string value = inlineValue ? longArgument.substr(equal + 1) : ARGF();
            p.dynamicDpiParametersSpecified = true;
            if (name == "dpi-href")
               p.dpiHref = parseDynamicDpiReference(value);
            else if (name == "dpi-slope")
               p.dpiSlope = parseFiniteNumericOption("--dpi-slope", value.c_str());
            else
               p.dpiMax = parseFiniteNumericOption("--dpi-max", value.c_str());
         }
         else
            throw std::string("unknown parameter --") + longArgument;
         break;
      }
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
                   throw std::string("Option '-M' requires the bundled model name 'pchip'.");
                break;
      case 'N': p.nparLimit  = std::atoi(ARGF()); break; // max npar value
      case 'S': p.seed       = std::atoi(ARGF()); break; // seed
      case 'n': p.correction = std::atof(ARGF()); break; // correction for noise
      case 'o': p.outfile    = ARGF(); break;            // output file
      case 'p': p.pvalue     = parseFiniteNumericOption("-p", ARGF()); break; // p-value
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
      case 'W': p.dpiWitnessFile = ARGF();               // DPI witness diagnostic
                if (p.dpiWitnessFile == "")
                   throw std::string("Option '-W' requires an output file.");
                break;
      default : throw std::string("unknown parameter ") + ARGC();
   }
   ARGEND;

   if (p.dynamicDpi && argc > 0)
      throw std::string("Unexpected positional argument '") + argv[0] +
            "' with '--dynamic-dpi'; this flag takes no value.";

   // The legacy parser stops at a positional token or '--'. Do not let a
   // new experimental option silently disappear in that unconsumed tail.
   // Inspect only remaining arguments: an already-consumed filename may
   // legitimately equal one of these option names.
   for (int argument = 0; argument < argc; ++argument)
   {
      const std::string remaining = argv[argument];
      const std::string name = remaining.substr(0, remaining.find('='));
      if (name == "--dynamic-dpi" || name == "--dpi-href" ||
          name == "--dpi-slope" || name == "--dpi-max")
         throw std::string("Dynamic DPI option '") + name +
               "' was not parsed; these options cannot follow positional arguments or '--'.";
   }

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
   if (std::string(adaptivePartitionKernelSchema()) != "sjaracne-apmi-v1")
      throw std::string("Bundled PCHIP data do not match the active AP-MI kernel schema.");
   p.threshold = Apmi::pchipCutoff(n, p.pvalue, p.nparLimit);
   p.thresholdMethod = "experimental PCHIP empirical cutoff";
   p.pchipPrimarySha256 = Apmi::pchipPrimarySha256();
   p.pchipDataSha256 = Apmi::pchipDataSha256();
   p.pchipSummarySha256 = Apmi::pchipSummarySha256();
   p.pchipActualM = n;
   std::cout << "[THRESHOLD] WARNING: PCHIP is experimental; uniform +/-10% "
                "null-rate accuracy and biological FDR are not established." << std::endl;
   if (n % 50 != 0)
      std::cout << "[THRESHOLD] WARNING: m=" << n
                << " is off the tested 50-step grid; this is unvalidated interpolation."
                << std::endl;
   if (p.pvalue < 0.0005)
      std::cout << "[THRESHOLD] WARNING: p is outside the principal assessment range "
                   "[0.0005,0.05]; this boundary is not a demonstrated failure threshold."
                << std::endl;
   if (p.pvalue <= 1e-5)
      std::cout << "[THRESHOLD] WARNING: extreme-tail Monte Carlo precision is limited; "
                   "10M validation draws expect only 100, 10 and 1 exceedances at "
                   "p=1e-5,1e-6,1e-7. Interpolation does not guarantee 10% accuracy."
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
      if (!p.pchipPrimarySha256.empty())
         std::cout << "[THRESHOLD] PCHIP cutoff binary64=" << std::setprecision(17)
                   << p.threshold << std::setprecision(6) << std::endl;
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

   bool dpiApplied = false;

   if (p.dynamicDpi)
   {
      p.dpiActualH = matrix.sourceRowCount();
      if (p.dpiActualH < 2)
         throw std::string("Experimental dynamic DPI requires at least two reconstructed "
                           "source rows; unresolved IDs do not count, but empty rows do.");
      p.eps = resolveDynamicDpiEpsilon(p.dpiActualH, p.dpiHref, p.dpiSlope, p.dpiMax);
      const std::streamsize savedPrecision = std::cout.precision();
      const std::ios::fmtflags savedFlags = std::cout.flags();
      std::cout << "[DPI] WARNING: dynamic epsilon is an experimental engineering policy; "
                   "biological accuracy and a universal anchor or slope are not established."
                << std::endl;
      std::cout << "[DPI] mode=experimental-dynamic H=" << p.dpiActualH
                << " H_ref=" << p.dpiHref
                << " slope=" << std::setprecision(17) << p.dpiSlope
                << " epsilon_max=" << p.dpiMax
                << " epsilon=" << p.eps
                << " epsilon_hex=" << std::hexfloat << p.eps << std::endl;
      std::cout.flags(savedFlags);
      std::cout.precision(savedPrecision);
   }

   if (p.eps != 1.0)
   {
      if (matrix.hasEnoughSourceRowsForDpi())
      {
         std::cout << "[NETWORK] Applying DPI ..." << std::endl;
         matrix.reduce(p.eps, ids, transfac);
         dpiApplied = true;
      }
      else
         std::cout << "[NETWORK] Skipping DPI: fewer than two source rows are "
                      "available; no eligible triangle can be formed."
                   << std::endl;
   }

   if (p.outfile == "")
      createOutfileName(p);

   if (!p.dpiWitnessFile.empty())
   {
      if (!dpiApplied)
         throw std::string("DPI witness diagnostic '-W' requires DPI to be "
                           "applied, but this network has fewer than two source "
                           "rows.");

      matrix.writeDpiWitnessDiagnostics(p, ids, transfac);
   }

   const DpiEdgeStatistics dpiStatistics = matrix.dpiEdgeStatistics(ids);
   matrix.write(data, ids, p);

   std::cout << "[DPI_STATS] pre_edges=" << dpiStatistics.preEdges
             << " pruned_edges=" << dpiStatistics.prunedEdges
             << " post_edges=" << dpiStatistics.postEdges
             << " dpi_applied=" << (dpiApplied ? 1 : 0) << std::endl;
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
