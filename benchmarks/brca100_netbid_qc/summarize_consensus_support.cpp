#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

struct EdgeMetric {
   std::string source;
   std::string target;
   double consensusMi;
   double miSum;
   int support;

   EdgeMetric(const std::string& sourceValue, const std::string& targetValue,
              double miValue)
      : source(sourceValue), target(targetValue), consensusMi(miValue),
        miSum(0.0), support(0) {}
};

std::vector<std::string> splitTabs(const std::string& line) {
   std::vector<std::string> fields;
   std::size_t start = 0;
   while (true) {
      const std::size_t end = line.find('\t', start);
      fields.push_back(line.substr(start, end == std::string::npos
                                           ? std::string::npos : end - start));
      if (end == std::string::npos)
         break;
      start = end + 1;
   }
   return fields;
}

bool endsWith(const std::string& value, const std::string& suffix) {
   return value.size() >= suffix.size() &&
          value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<std::string> adjacencyFiles(const std::string& directory) {
   DIR* handle = opendir(directory.c_str());
   if (handle == NULL)
      throw std::runtime_error("Cannot open adjacency directory: " + directory);
   std::vector<std::string> files;
   while (dirent* entry = readdir(handle)) {
      const std::string name(entry->d_name);
      if (endsWith(name, ".adj"))
         files.push_back(directory + "/" + name);
   }
   closedir(handle);
   std::sort(files.begin(), files.end());
   if (files.size() != 100)
      throw std::runtime_error("Expected exactly 100 adjacency files, got " +
                               std::to_string(files.size()));
   return files;
}

int internGene(const std::string& gene,
               std::unordered_map<std::string, int>& geneIds) {
   const std::unordered_map<std::string, int>::const_iterator found =
      geneIds.find(gene);
   if (found != geneIds.end())
      return found->second;
   const int id = static_cast<int>(geneIds.size());
   geneIds.insert(std::make_pair(gene, id));
   return id;
}

std::uint64_t edgeKey(int source, int target) {
   return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(source)) << 32) |
          static_cast<std::uint32_t>(target);
}

double parseFiniteDouble(const std::string& value, const std::string& context) {
   char* end = NULL;
   const double parsed = std::strtod(value.c_str(), &end);
   if (end == value.c_str() || *end != '\0' || !std::isfinite(parsed))
      throw std::runtime_error("Invalid numeric value in " + context + ": " + value);
   return parsed;
}

void stripCarriageReturn(std::string& line) {
   if (!line.empty() && line[line.size() - 1] == '\r')
      line.erase(line.size() - 1);
}

} // namespace

int main(int argc, char* argv[]) {
   try {
      if (argc != 4) {
         std::cerr << "Usage: summarize_consensus_support CONSENSUS_NCOL "
                      "ADJACENCY_DIR OUTPUT_TSV\n";
         return 2;
      }
      const std::string consensusPath(argv[1]);
      const std::string adjacencyDirectory(argv[2]);
      const std::string outputPath(argv[3]);

      std::ifstream consensus(consensusPath.c_str());
      if (!consensus)
         throw std::runtime_error("Cannot open consensus: " + consensusPath);
      std::string line;
      if (!std::getline(consensus, line))
         throw std::runtime_error("Missing consensus header");
      stripCarriageReturn(line);
      if (line != "source\ttarget\tsource.symbol\ttarget.symbol\tMI\tpearson\t"
                   "spearman\tslope\tp-value")
         throw std::runtime_error("Unexpected consensus header");

      std::unordered_map<std::string, int> geneIds;
      std::unordered_map<std::uint64_t, std::size_t> retained;
      std::vector<EdgeMetric> metrics;
      while (std::getline(consensus, line)) {
         stripCarriageReturn(line);
         const std::vector<std::string> fields = splitTabs(line);
         if (fields.size() != 9)
            throw std::runtime_error("Malformed consensus row");
         const int source = internGene(fields[0], geneIds);
         const int target = internGene(fields[1], geneIds);
         const std::uint64_t key = edgeKey(source, target);
         if (retained.find(key) != retained.end())
            throw std::runtime_error("Duplicate consensus edge");
         const double mi = parseFiniteDouble(fields[4], "consensus");
         retained.insert(std::make_pair(key, metrics.size()));
         metrics.push_back(EdgeMetric(fields[0], fields[1], mi));
      }
      if (metrics.empty())
         throw std::runtime_error("Consensus contains no edges");

      const std::vector<std::string> files = adjacencyFiles(adjacencyDirectory);
      std::string sourceText;
      std::string targetText;
      for (std::size_t fileIndex = 0; fileIndex < files.size(); ++fileIndex) {
         std::ifstream adjacency(files[fileIndex].c_str());
         if (!adjacency)
            throw std::runtime_error("Cannot open adjacency: " + files[fileIndex]);
         std::unordered_set<std::uint64_t> retainedEdgesSeen;
         while (std::getline(adjacency, line)) {
            stripCarriageReturn(line);
            if (line.empty() || line[0] == '>')
               continue;
            const std::size_t sourceEnd = line.find('\t');
            if (sourceEnd == std::string::npos)
               throw std::runtime_error("Malformed adjacency source row");
            sourceText.assign(line, 0, sourceEnd);
            const std::unordered_map<std::string, int>::const_iterator sourceFound =
               geneIds.find(sourceText);
            std::size_t position = sourceEnd + 1;
            while (position < line.size()) {
               const std::size_t targetEnd = line.find('\t', position);
               if (targetEnd == std::string::npos)
                  throw std::runtime_error("Missing adjacency MI field");
               targetText.assign(line, position, targetEnd - position);
               const std::size_t miStart = targetEnd + 1;
               const std::size_t miEnd = line.find('\t', miStart);
               if (sourceFound != geneIds.end()) {
                  const std::unordered_map<std::string, int>::const_iterator targetFound =
                     geneIds.find(targetText);
                  if (targetFound != geneIds.end()) {
                     const std::uint64_t key = edgeKey(sourceFound->second,
                                                       targetFound->second);
                     const std::unordered_map<std::uint64_t, std::size_t>::const_iterator
                        edgeFound = retained.find(key);
                      if (edgeFound != retained.end()) {
                         if (!retainedEdgesSeen.insert(key).second)
                            throw std::runtime_error(
                               "Duplicate retained edge in adjacency: " +
                               files[fileIndex]);
                         const std::string miText = line.substr(
                           miStart, miEnd == std::string::npos
                                       ? std::string::npos : miEnd - miStart);
                        EdgeMetric& metric = metrics[edgeFound->second];
                        metric.miSum += parseFiniteDouble(miText, "adjacency");
                        ++metric.support;
                     }
                  }
               }
               if (miEnd == std::string::npos)
                  break;
               position = miEnd + 1;
            }
         }
         if ((fileIndex + 1) % 10 == 0)
            std::cerr << "Scanned " << (fileIndex + 1) << "/100 adjacency files\n";
      }

      int mismatchCount = 0;
      for (std::size_t index = 0; index < metrics.size(); ++index) {
         const EdgeMetric& metric = metrics[index];
         if (metric.support <= 0 || metric.support > 100)
            throw std::runtime_error("Invalid support count for retained edge");
         const double meanMi = metric.miSum / metric.support;
         mismatchCount +=
            std::fabs(meanMi - metric.consensusMi) > 0.0000500001;
      }
      if (mismatchCount != 0)
         throw std::runtime_error("Consensus MI round-trip mismatch count: " +
                                  std::to_string(mismatchCount));

      const std::string partialOutputPath = outputPath + ".partial";
      std::ofstream output(partialOutputPath.c_str());
      if (!output)
         throw std::runtime_error("Cannot create output: " + partialOutputPath);
      output << "source\ttarget\tconsensus_MI\tsupport_count\tsupport_fraction\t"
                "mean_observed_MI\tconsensus_MI_roundtrip_match\n";
      output << std::setprecision(17);
      for (std::size_t index = 0; index < metrics.size(); ++index) {
         const EdgeMetric& metric = metrics[index];
         const double meanMi = metric.miSum / metric.support;
         const bool match = std::fabs(meanMi - metric.consensusMi) <= 0.0000500001;
         output << metric.source << '\t' << metric.target << '\t'
                << metric.consensusMi << '\t' << metric.support << '\t'
                << (static_cast<double>(metric.support) / 100.0) << '\t'
                 << meanMi << '\t' << (match ? 1 : 0) << '\n';
      }
      output.close();
      if (!output)
         throw std::runtime_error("Failed writing output: " + partialOutputPath);
      if (std::rename(partialOutputPath.c_str(), outputPath.c_str()) != 0)
         throw std::runtime_error("Cannot finalize output: " + outputPath);
      std::cerr << "Summarized " << metrics.size() << " retained edges\n";
      return 0;
   }
   catch (const std::exception& error) {
      std::cerr << "ERROR: " << error.what() << '\n';
      return 1;
   }
}
