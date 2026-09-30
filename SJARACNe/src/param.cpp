//------------------------------------------------------------------------------------
// Copyright (C) 2003, 2004  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017
//------------------------------------------------------------------------------------

#include <cerrno>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include "param.h"

//------------------------------------------------------------------------------------

const double Parameter::default_threshold  =  0.00; // mi threshold
const double Parameter::default_pvalue     =  1.00; // p-value for mi threshold
const double Parameter::default_eps        =  1.00; // DPI tolerance
const double Parameter::default_sigma      = 99.00; // kernel width
const int    Parameter::default_sample     =  0;    // sample number
const double Parameter::default_percent    =  0.35; // high/low percentage
const double Parameter::default_mean       =  0.00; // filter mean
const double Parameter::default_cv         =  0.00; // filter standard deviation
const double Parameter::default_correction =  0.00; // array measurement noise level
const int    Parameter::default_nparLimit  = 20;    // max allowed value of npar
const int    Parameter::default_seed       = 1;     // Initial seed for random number generator
//------------------------------------------------------------------------------------

bool equalIgnoreCase(std::string a, std::string b)
{
   int len = a.length();
   if (b.length() != len)
      return false;

   int i = 0;
   while (i < len && std::toupper(a[i]) == std::toupper(b[i]))
      i++;

   return (i == len);
}

//------------------------------------------------------------------------------------

void checkParameter(Parameter &p)
{
   if (p.infile == "")
      throw std::string("No input file specified!");

   if (p.hub != "")
      p.hub = "_" + p.hub;

   if (p.subnetfile != "" && p.hub != "")
      throw std::string("Either supply one hub gene by '-h' or multiple genes in a "
                        "file by '-s', but not both!");

   if (p.sample < 0)
      throw std::string("Legacy bootstrap sample number '-r' must be nonnegative!");

   if (p.sample > 0 && p.subsampleSpec != "")
      throw std::string("Options '-r' (legacy bootstrap) and '-u' (unique "
                        "subsampling) cannot be used together!");

   if (p.adjfile != "" && (p.subsampleSpec != "" || p.sample > 0))
      throw std::string("Sampling options '-u' and '-r' cannot be used with an "
                        "existing adjacency matrix supplied by '-j'.");

   if (p.condition != "+" && p.condition != "-" && p.condition != "")
      throw std::string("Condition must be '+' or '-'!");

   if ((p.condition == "+" || p.condition == "-") && p.controlId == "")
      throw std::string("Control gene ID must be specified using '-c'!");

   if (p.sigma != Parameter::default_sigma && (p.sigma <= 0.0 || p.sigma >= 1.0))
      throw std::string("Kernel width '-k' must be within (0,1)!");

   if (p.threshold < 0.0)
      throw std::string("MI threshold '-t' must be nonnegative!");

   if (p.thresholdSpecified && p.pvalue != 1.0)
      std::cout << "P-value will not be used, since a threshold has been specified."
                << std::endl;

   // An explicit threshold makes the model irrelevant, including the valid
   // preserve-all request '-t 0'. Otherwise, both legacy replacement sampling and
   // adjacency replay violate the estimator/model contract.
   if (!p.nullModelFile.empty() && !p.thresholdSpecified && p.sample > 0)
      throw std::string("Estimator-matched AP-MI null models ('-M') cannot be used "
                        "with legacy replacement sampling '-r'; use fixed-size "
                        "sampling without replacement '-u' or an explicit '-t'.");

   if (!p.nullModelFile.empty() && !p.thresholdSpecified && !p.adjfile.empty())
      throw std::string("Estimator-matched AP-MI null models ('-M') cannot be used "
                        "while replaying an existing adjacency matrix with '-j'; "
                        "use the threshold recorded in that matrix or an explicit '-t'.");

   if (p.pvalue <= 0.0 || p.pvalue > 1.0)
      throw std::string("P-value '-p' must be in the range (0,1]!");

   if (p.eps < 0.0 || p.eps > 1.0)
      throw std::string("DPI tolerance '-e' must be within [0,1]!");

   if (p.percent <= 0.0 || p.percent >= 1.0)
      throw std::string("Percentage microarray must be within (0,1)!");

   if (p.mean < 0.0)
      throw std::string("Gene filter mean must be nonnegative!");

   if (p.cv < 0.0)
      throw std::string("Gene filter cv (coefficient of variance) must be "
                        "nonnegative!");

   if (p.correction < 0.0)
      throw std::string("Array measurement noise level '-n' must be nonnegative!");

   if (!equalIgnoreCase(p.verbose, "on") && !equalIgnoreCase(p.verbose, "off"))
      throw std::string("Verbose '-v' must be 'on' or 'off'!");

   if (p.nparLimit < 1)
      throw std::string("Maximum allowed value of npar must be positive!");

   if (p.home_dir != "./")
   {
      int len = p.home_dir.length();
      int b   = p.home_dir.find_last_of("/");

      if (b == std::string::npos || b < len - 1)
         p.home_dir += "/";
   }
}

//------------------------------------------------------------------------------------
// resolveSubsampleSize() converts either an exact count (for example, "80") or an
// explicit percentage (for example, "80%") to a fixed observation count.

int resolveSubsampleSize(const std::string& spec, int populationSize)
{
   if (populationSize < 1)
      throw std::string("Cannot subsample an empty observation population.");

   if (spec.empty())
      throw std::string("Unique subsampling '-u' requires a count or percentage.");

   const bool isPercent = spec[spec.length() - 1] == '%';
   const std::string value =
      (isPercent ? spec.substr(0, spec.length() - 1) : spec);

   if (value.empty())
      throw std::string("Unique subsampling '-u' must be an integer count or a "
                        "percentage such as 80%.");

   errno = 0;
   char *end = NULL;
   int sampleSize = 0;

   if (isPercent)
   {
      const double percentage = std::strtod(value.c_str(), &end);

      if (errno == ERANGE || end == value.c_str() || *end != '\0' ||
          !std::isfinite(percentage))
         throw std::string("Unique subsampling '-u' has an invalid percentage: ") +
                           spec;

      if (percentage <= 0.0 || percentage > 100.0)
         throw std::string("Unique subsampling percentage '-u' must be within "
                           "(0%,100%].");

      const double resolved =
         std::ceil(percentage * static_cast<double>(populationSize) / 100.0);

      if (resolved > INT_MAX)
         throw std::string("Unique subsample size is too large.");

      sampleSize = static_cast<int>(resolved);
   }
   else
   {
      const long parsed = std::strtol(value.c_str(), &end, 10);

      if (errno == ERANGE || end == value.c_str() || *end != '\0' ||
          parsed > INT_MAX || parsed < INT_MIN)
         throw std::string("Unique subsampling '-u' has an invalid observation "
                           "count: ") + spec;

      sampleSize = static_cast<int>(parsed);
   }

   if (sampleSize < 2)
      throw std::string("Unique subsampling must select at least 2 observations; "
                        "requested: ") + spec + ".";

   if (sampleSize > populationSize)
   {
      std::ostringstream message;
      message << "Unique subsampling requested " << sampleSize
              << " observations, but only " << populationSize
              << " are eligible.";
      throw message.str();
   }

   return sampleSize;
}

//------------------------------------------------------------------------------------
// appendProbeId() normalizes and stores one identifier from a hub-list record.

static void appendProbeId(std::string gid, bool first_record,
                          std::vector<std::string>& probe_list,
                          std::set<std::string>& seen,
                          int& duplicate_count)
{
   const std::string utf8_bom("\xEF\xBB\xBF", 3);

   if (first_record && gid.compare(0, utf8_bom.length(), utf8_bom) == 0)
      gid.erase(0, utf8_bom.length());

   const std::string whitespace(" \t\n\r\f\v");
   std::string::size_type first = gid.find_first_not_of(whitespace);

   if (first == std::string::npos)
      return;

   std::string::size_type last = gid.find_last_not_of(whitespace);
   gid = gid.substr(first, last - first + 1);
   gid = "_" + gid;

   if (seen.insert(gid).second)
      probe_list.push_back(gid);
   else
      duplicate_count++;
}

//------------------------------------------------------------------------------------
// readProbeList() reads and normalizes the list of nodes used by -s and -l.

static int readProbeList(const std::string& infilename,
                         std::vector<std::string>& probe_list,
                         int& duplicate_count)
{
   std::ifstream in(infilename.c_str(), std::ios::binary);
   if (!in.is_open())
      throw "Unable to open " + infilename;

   probe_list.clear();
   duplicate_count = 0;

   std::set<std::string> seen;
   std::string line;
   bool first_record = true;
   char c;

   while (in.get(c))
   {
      if (c == '\n' || c == '\r')
      {
         appendProbeId(line, first_record, probe_list, seen, duplicate_count);
         line.clear();
         first_record = false;

         // Consume LF as part of a CRLF delimiter. A lone CR is also a valid
         // delimiter, which keeps classic Mac and mixed-ending files readable.
         if (c == '\r' && in.peek() == '\n')
            in.get(c);
      }
      else
         line += c;
   }

   appendProbeId(line, first_record, probe_list, seen, duplicate_count);

   return probe_list.size();
}

//------------------------------------------------------------------------------------

void displayParameter(Parameter &p)
{
   std::cout << std::endl;

   std::cout << "[PARA] Input file:    " << p.infile  << std::endl;
   std::cout << "[PARA] Output file:   " << p.outfile << std::endl;

   if (p.thresholdSpecified)
      std::cout << "[PARA] MI threshold:  " << p.threshold << std::endl;
   else
      std::cout << "[PARA] MI P-value:    " << p.pvalue    << std::endl;

   if (!p.nullModelFile.empty())
      std::cout << "[PARA] AP-MI null model: " << p.nullModelFile << std::endl;

   std::cout << "[PARA] DPI tolerance: " << p.eps << std::endl;

   if (p.correction > 0.0)
      std::cout << "[PARA] Correction for MI estimation (array noise level: "
                << p.correction << ")" << std::endl;

   if (p.subnetfile != "")
   {
      int duplicate_count = 0;
      std::cout << "[PARA] Subset of probes to reconstruct: " << p.subnetfile
                << " (" << readProbeList(p.subnetfile, p.subnet, duplicate_count)
                << ")" << std::endl;

      if (duplicate_count > 0)
         std::cout << "[PARA] Duplicate subnetwork probes ignored: "
                   << duplicate_count << std::endl;
   }

   if (p.hub != "")
      std::cout << "[PARA] Hub probe to reconstruct: " << p.hub << std::endl;

   if (p.controlId != "")
   {
      std::cout << "[PARA] Control gene:  " << p.controlId << std::endl;
      std::cout << "[PARA] Condition:     " << p.condition << std::endl;
      std::cout << "[PARA] Percentage:    " << p.percent   << std::endl;
   }

   if (p.subsampleSpec != "")
      std::cout << "[PARA] Subsampling:   " << p.subsampleSpec
                << " without replacement" << std::endl;
   else if (p.sample > 0)
      std::cout << "[PARA] Sampling:      legacy bootstrap with replacement"
                << std::endl;

   if (p.annotfile != "")
   {
      int duplicate_count = 0;
      std::cout << "[PARA] TF annotation list: " << p.annotfile
                << " (" << readProbeList(p.annotfile, p.tf_list, duplicate_count)
                << ")" << std::endl;

      if (duplicate_count > 0)
         std::cout << "[PARA] Duplicate TF annotation probes ignored: "
                   << duplicate_count << std::endl;
   }

   if (p.mean != 0.0 || p.cv != 0.0)
   {
      std::cout << "[PARA] Filter mean:   " << p.mean << std::endl;
      std::cout << "[PARA] Filter CV:     " << p.cv   << std::endl;
   }

   std::cout << "[PARA] Npar limit:    " << p.nparLimit << std::endl;
}

//------------------------------------------------------------------------------------

static std::string getFileName(const std::string& matrixName)
{
   int b = matrixName.find_last_of("/");
   if (b == std::string::npos)
      b = matrixName.find_last_of("\\");

   std::string dirname(""), basename(matrixName);

   if (b != std::string::npos)
   {
      // Extract the directory and the filename if path is included
      basename = matrixName.substr(b + 1);
      dirname  = matrixName.substr(0, b) + "/";
   }

   int c = basename.find_last_of(".");
   if (c != std::string::npos)
       basename = basename.substr(0, c);

   return dirname + basename;
}

//------------------------------------------------------------------------------------

void createOutfileName(Parameter &p)
{
   std::string filename = getFileName(p.infile);

   if (p.hub != "")
      filename += "_h" + p.hub;

   if (p.controlId != "")
      filename += "_c" + p.controlId + (p.condition == "+" ? "H" : "L");

   char buffer[20];

   std::sprintf(buffer, "%0.3g", p.sigma);
   filename += std::string("_k") + buffer;

   if (p.threshold > 0.0)
   {
      std::sprintf(buffer, "%0.2g", p.threshold);
      filename += std::string("_t") + buffer;
   }

   if (p.eps < 1.0)
   {
      std::sprintf(buffer, "%0.2g", p.eps);
      filename += std::string("_e") + buffer;
   }

   if (p.sample > 0)
   {
      std::sprintf(buffer, "%03i", p.sample);
      filename += std::string("_r") + buffer;
   }

   if (p.subsampleSpec != "")
   {
      std::string samplingLabel = p.subsampleSpec;
      std::string::size_type percent = samplingLabel.find('%');
      if (percent != std::string::npos)
         samplingLabel.replace(percent, 1, "pct");
      filename += "_u" + samplingLabel;

      std::sprintf(buffer, "%i", p.seed);
      filename += std::string("_S") + buffer;
   }

   p.outfile = filename + ".adj";
}
