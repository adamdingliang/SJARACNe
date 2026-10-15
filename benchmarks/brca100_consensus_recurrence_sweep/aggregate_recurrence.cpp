#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

const std::size_t kExpectedRuns = 100;
const std::uint16_t kMinimumSupport = 6;

struct EdgeStats {
   double miSum;
   std::uint16_t support;

   EdgeStats() : miSum(0.0), support(0) {}
};

struct RunStats {
   std::string filename;
   std::size_t edgeCount;

   RunStats(const std::string& filenameValue, std::size_t edgeCountValue)
      : filename(filenameValue), edgeCount(edgeCountValue) {}
};

bool endsWith(const std::string& value, const std::string& suffix) {
   return value.size() >= suffix.size() &&
          value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string joinPath(const std::string& directory, const std::string& filename) {
   if (!directory.empty() && directory[directory.size() - 1] == '/')
      return directory + filename;
   return directory + "/" + filename;
}

bool pathExists(const std::string& path) {
   struct stat status;
   return stat(path.c_str(), &status) == 0;
}

std::vector<std::string> adjacencyFiles(const std::string& directory) {
   DIR* handle = opendir(directory.c_str());
   if (handle == NULL)
      throw std::runtime_error("Cannot open adjacency directory: " + directory);

   std::vector<std::string> files;
   while (dirent* entry = readdir(handle)) {
      const std::string name(entry->d_name);
      if (endsWith(name, ".adj"))
         files.push_back(name);
   }
   closedir(handle);
   std::sort(files.begin(), files.end());
   if (files.size() != kExpectedRuns) {
      throw std::runtime_error(
         "Expected exactly 100 .adj files, got " + std::to_string(files.size()));
   }
   return files;
}

std::vector<std::string> splitTabs(const std::string& line) {
   std::vector<std::string> fields;
   std::size_t start = 0;
   while (true) {
      const std::size_t end = line.find('\t', start);
      fields.push_back(line.substr(
         start, end == std::string::npos ? std::string::npos : end - start));
      if (end == std::string::npos)
         break;
      start = end + 1;
   }
   return fields;
}

void stripCarriageReturn(std::string& line) {
   if (!line.empty() && line[line.size() - 1] == '\r')
      line.erase(line.size() - 1);
}

double parseFiniteDouble(const std::string& value, const std::string& context) {
   char* end = NULL;
   const double parsed = std::strtod(value.c_str(), &end);
   if (end == value.c_str() || *end != '\0' || !std::isfinite(parsed))
      throw std::runtime_error("Invalid numeric value in " + context + ": " + value);
   return parsed;
}

std::uint32_t internGene(
   const std::string& gene,
   std::unordered_map<std::string, std::uint32_t>& geneIds,
   std::vector<std::string>& geneNames) {
   const std::unordered_map<std::string, std::uint32_t>::const_iterator found =
      geneIds.find(gene);
   if (found != geneIds.end())
      return found->second;
   if (geneNames.size() >=
       static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
      throw std::runtime_error("Too many distinct gene identifiers");
   }
   const std::uint32_t id = static_cast<std::uint32_t>(geneNames.size());
   geneIds.insert(std::make_pair(gene, id));
   geneNames.push_back(gene);
   return id;
}

std::uint64_t edgeKey(std::uint32_t source, std::uint32_t target) {
   return (static_cast<std::uint64_t>(source) << 32) |
          static_cast<std::uint64_t>(target);
}

std::uint32_t sourceId(std::uint64_t key) {
   return static_cast<std::uint32_t>(key >> 32);
}

std::uint32_t targetId(std::uint64_t key) {
   return static_cast<std::uint32_t>(key & UINT64_C(0xffffffff));
}

std::string edgeDescription(
   std::uint64_t key, const std::vector<std::string>& geneNames) {
   return geneNames[sourceId(key)] + "----" + geneNames[targetId(key)];
}

void requireFreshOutput(const std::string& path) {
   if (pathExists(path))
      throw std::runtime_error("Refusing to overwrite existing output: " + path);
}

} // namespace

int main(int argc, char* argv[]) {
   try {
      if (argc != 5) {
         std::cerr
            << "Usage: aggregate_recurrence ADJACENCY_DIR EDGE_OUTPUT_TSV "
               "RUN_COUNTS_TSV SUMMARY_TSV\n";
         return 2;
      }

      const std::string adjacencyDirectory(argv[1]);
      const std::string edgeOutputPath(argv[2]);
      const std::string runOutputPath(argv[3]);
      const std::string summaryOutputPath(argv[4]);
      if (edgeOutputPath == runOutputPath || edgeOutputPath == summaryOutputPath ||
          runOutputPath == summaryOutputPath) {
         throw std::runtime_error("Output paths must be distinct");
      }
      requireFreshOutput(edgeOutputPath);
      requireFreshOutput(runOutputPath);
      requireFreshOutput(summaryOutputPath);

      const std::vector<std::string> files = adjacencyFiles(adjacencyDirectory);
      std::unordered_map<std::string, std::uint32_t> geneIds;
      std::vector<std::string> geneNames;
      std::unordered_map<std::uint64_t, EdgeStats> edges;
      std::vector<RunStats> runs;
      runs.reserve(files.size());

      std::string line;
      for (std::size_t fileIndex = 0; fileIndex < files.size(); ++fileIndex) {
         const std::string path = joinPath(adjacencyDirectory, files[fileIndex]);
         std::ifstream adjacency(path.c_str());
         adjacency.imbue(std::locale::classic());
         if (!adjacency)
            throw std::runtime_error("Cannot open adjacency file: " + path);

         std::unordered_set<std::uint64_t> seen;
         std::size_t lineNumber = 0;
         while (std::getline(adjacency, line)) {
            ++lineNumber;
            stripCarriageReturn(line);
            if (line.empty() || line[0] == '>')
               continue;

            const std::vector<std::string> fields = splitTabs(line);
            if (fields.size() < 3 || fields.size() % 2 == 0 || fields[0].empty()) {
               throw std::runtime_error(
                  "Malformed adjacency row at " + path + ":" +
                  std::to_string(lineNumber));
            }
            const std::uint32_t source =
               internGene(fields[0], geneIds, geneNames);
            for (std::size_t index = 1; index < fields.size(); index += 2) {
               if (fields[index].empty()) {
                  throw std::runtime_error(
                     "Empty target identifier at " + path + ":" +
                     std::to_string(lineNumber));
               }
               const std::uint32_t target =
                  internGene(fields[index], geneIds, geneNames);
               const std::uint64_t key = edgeKey(source, target);
               if (!seen.insert(key).second) {
                  throw std::runtime_error(
                     "Duplicate ordered edge " + edgeDescription(key, geneNames) +
                     " in adjacency file " + files[fileIndex]);
               }
               const double mi = parseFiniteDouble(
                  fields[index + 1], path + ":" + std::to_string(lineNumber));
               EdgeStats& stats = edges[key];
               stats.miSum += mi;
               if (stats.support >= kExpectedRuns)
                  throw std::runtime_error("Edge support exceeds bootstrap count");
               ++stats.support;
            }
         }
         if (!adjacency.eof())
            throw std::runtime_error("Failed while reading adjacency file: " + path);

         runs.push_back(RunStats(files[fileIndex], seen.size()));
         if ((fileIndex + 1) % 10 == 0) {
            std::cerr << "Scanned " << (fileIndex + 1) << "/100 adjacency files\n";
         }
      }

      std::vector<std::uint64_t> retained;
      retained.reserve(edges.size());
      for (std::unordered_map<std::uint64_t, EdgeStats>::const_iterator item =
              edges.begin();
           item != edges.end(); ++item) {
         if (item->second.support >= kMinimumSupport)
            retained.push_back(item->first);
      }
      std::sort(
         retained.begin(), retained.end(),
         [&geneNames](std::uint64_t left, std::uint64_t right) {
            const std::string& leftSource = geneNames[sourceId(left)];
            const std::string& rightSource = geneNames[sourceId(right)];
            if (leftSource != rightSource)
               return leftSource < rightSource;
            return geneNames[targetId(left)] < geneNames[targetId(right)];
         });

      std::ofstream edgeOutput(edgeOutputPath.c_str());
      std::ofstream runOutput(runOutputPath.c_str());
      std::ofstream summaryOutput(summaryOutputPath.c_str());
      edgeOutput.imbue(std::locale::classic());
      runOutput.imbue(std::locale::classic());
      summaryOutput.imbue(std::locale::classic());
      if (!edgeOutput)
         throw std::runtime_error("Cannot create edge output: " + edgeOutputPath);
      if (!runOutput)
         throw std::runtime_error("Cannot create run-count output: " + runOutputPath);
      if (!summaryOutput)
         throw std::runtime_error("Cannot create summary output: " + summaryOutputPath);

      edgeOutput << "source\ttarget\tmean_observed_MI\tconsensus_MI\t"
                    "support_count\tsupport_fraction\n";
      for (std::size_t index = 0; index < retained.size(); ++index) {
         const std::uint64_t key = retained[index];
         const EdgeStats& stats = edges.find(key)->second;
         const double meanMi = stats.miSum / static_cast<double>(stats.support);
         edgeOutput << geneNames[sourceId(key)] << '\t' << geneNames[targetId(key)]
                    << '\t' << std::defaultfloat << std::setprecision(17) << meanMi
                    << '\t' << std::fixed << std::setprecision(4) << meanMi << '\t'
                    << std::defaultfloat << static_cast<unsigned int>(stats.support)
                    << '\t' << std::setprecision(17)
                    << (static_cast<double>(stats.support) /
                        static_cast<double>(kExpectedRuns))
                    << '\n';
      }

      runOutput << "run_ordinal\tadjacency_file\tedge_count\n";
      for (std::size_t index = 0; index < runs.size(); ++index) {
         runOutput << (index + 1) << '\t' << runs[index].filename << '\t'
                   << runs[index].edgeCount << '\n';
      }

      summaryOutput << "metric\tvalue\n"
                    << "bootstrap_runs\t" << files.size() << '\n'
                    << "minimum_support\t" << kMinimumSupport << '\n'
                    << "union_edges\t" << edges.size() << '\n'
                    << "retained_edges\t" << retained.size() << '\n';

      edgeOutput.close();
      runOutput.close();
      summaryOutput.close();
      if (!edgeOutput)
         throw std::runtime_error("Failed writing edge output: " + edgeOutputPath);
      if (!runOutput)
         throw std::runtime_error("Failed writing run-count output: " + runOutputPath);
      if (!summaryOutput)
         throw std::runtime_error("Failed writing summary output: " + summaryOutputPath);

      std::cerr << "Aggregated " << edges.size() << " union edges; retained "
                << retained.size() << " with support >= " << kMinimumSupport
                << "\n";
      return 0;
   }
   catch (const std::exception& error) {
      std::cerr << "ERROR: " << error.what() << '\n';
      return 1;
   }
}
