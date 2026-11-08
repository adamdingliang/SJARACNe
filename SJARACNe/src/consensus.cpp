#include "consensus.h"
#include "matrix.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include <sys/stat.h>
#include <unordered_map>

namespace {

std::string canonicalPath(const std::string& path)
{
#ifdef _WIN32
   char *resolved = _fullpath(NULL, path.c_str(), 0);
#else
   char *resolved = realpath(path.c_str(), NULL);
#endif
   if (resolved == NULL)
      throw std::string("Cannot resolve consensus input '") + path + "': " + std::strerror(errno);
   const std::string result(resolved);
   std::free(resolved);
   return result;
}

std::vector<std::string> manifestPaths(const std::string& manifest)
{
   const std::string resolvedManifest = canonicalPath(manifest);
   const std::size_t separator = resolvedManifest.find_last_of("/\\");
   const std::string directory = separator == std::string::npos ? "" : resolvedManifest.substr(0, separator + 1);
   std::ifstream input(resolvedManifest.c_str());
   if (!input.is_open())
      throw std::string("Cannot open consensus manifest '") + manifest + "'.";
   std::vector<std::string> paths;
   std::set<std::string> unique;
   std::set<std::pair<std::uint64_t, std::uint64_t> > identities;
   std::string line;
   std::size_t lineNumber = 0;
   while (std::getline(input, line))
   {
      ++lineNumber;
      if (!line.empty() && line[line.size() - 1] == '\r')
         line.resize(line.size() - 1);
      if (line.empty() || line.find('\0') != std::string::npos)
      {
         std::ostringstream message;
         message << "Invalid empty/NUL path at consensus manifest line " << lineNumber << ".";
         throw message.str();
      }
      const bool absolute = line[0] == '/' || line[0] == '\\' ||
         (line.size() > 1 && line[1] == ':');
      const std::string path = canonicalPath(absolute ? line : directory + line);
      if (!unique.insert(path).second)
         throw std::string("Duplicate consensus input path: ") + path;
      struct stat identity;
      if (stat(path.c_str(), &identity) != 0)
         throw std::string("Cannot stat raw consensus input: ") + path;
      if (!identities.insert(std::make_pair(static_cast<std::uint64_t>(identity.st_dev),
                                            static_cast<std::uint64_t>(identity.st_ino))).second)
         throw std::string("Duplicate consensus input file identity (hardlink): ") + path;
      paths.push_back(path);
      if (paths.size() > static_cast<std::size_t>(INT_MAX))
         throw std::string("Too many consensus inputs; supported run count is at most INT_MAX.");
   }
   if (input.bad())
      throw std::string("Error reading consensus manifest '") + manifest + "'.";
   if (paths.empty())
      throw std::string("Consensus manifest contains no inputs.");
   return paths;
}

bool sameExistingFile(const std::string& first, const std::string& second)
{
   struct stat a, b;
   return !first.empty() && !second.empty() &&
      stat(first.c_str(), &a) == 0 && stat(second.c_str(), &b) == 0 &&
      a.st_dev == b.st_dev && a.st_ino == b.st_ino;
}

class RawReader
{
public:
   explicit RawReader(const std::string& name)
      : filename(name), input(name.c_str()), lineNumber(0), dataStarted(false), dpiHeaderSeen(false)
   {
      if (!input.is_open())
         throw std::string("Cannot open raw consensus input '") + filename + "': " + std::strerror(errno);
   }

   bool next(std::string& row)
   {
      while (std::getline(input, row))
      {
         ++lineNumber;
         if (!row.empty() && row[row.size() - 1] == '\r')
            row.resize(row.size() - 1);
         if (!dataStarted && !row.empty() && row[0] == '>')
         {
            std::istringstream header(row.substr(1));
            std::string key, subkey, value;
            header >> key >> subkey >> value;
            if (key == "DPI" && subkey == "tolerance")
            {
               if (dpiHeaderSeen)
                  fail("duplicate DPI tolerance header");
               dpiHeaderSeen = true;
               char *end = NULL;
               const double epsilon = std::strtod(value.c_str(), &end);
               if (end == value.c_str() || *end != '\0' || !std::isfinite(epsilon) || epsilon != 1.0)
                  fail("input is not an unpruned raw MI graph (DPI tolerance must be 1)");
            }
            continue;
         }
         if (!dpiHeaderSeen)
            fail("missing required DPI tolerance header before raw MI rows");
         dataStarted = true;
         if (row.empty() || row[0] == '>')
            fail("blank rows or headers after data are not canonical raw MI rows");
         return true;
      }
      if (input.bad())
         fail("I/O error while reading raw MI rows");
      return false;
   }

   void fail(const std::string& reason) const
   {
      std::ostringstream message;
      message << "Raw consensus input '" << filename << "', line " << lineNumber << ": " << reason << ".";
      throw message.str();
   }

private:
   std::string filename;
   std::ifstream input;
   std::size_t lineNumber;
   bool dataStarted, dpiHeaderSeen;
};

double parseMi(const std::string& token, const RawReader& reader)
{
   errno = 0;
   char *end = NULL;
   const double value = std::strtod(token.c_str(), &end);
   if (end == token.c_str() || *end != '\0' || errno == ERANGE ||
       !std::isfinite(value) || value < 0.0)
      reader.fail("MI must be a finite, nonnegative complete numeric token, got '" + token + "'");
   return value;
}

} // namespace

void validateNativeOutputPaths(const Parameter& p,
                               const std::vector<std::string>& rawInputs)
{
   std::vector<std::string> protectedInputs(rawInputs);
   protectedInputs.push_back(p.infile);
   protectedInputs.push_back(p.subnetfile);
   protectedInputs.push_back(p.annotfile);
   protectedInputs.push_back(p.consensusInputs);
   const std::string outputs[] = {p.outfile, p.dpiWitnessFile};
   for (std::size_t i = 0; i < 2; ++i)
      for (std::vector<std::string>::const_iterator input = protectedInputs.begin(); input != protectedInputs.end(); ++input)
         if (sameExistingFile(outputs[i], *input))
            throw std::string("Native output '") + outputs[i] + "' aliases a protected input '" + *input + "'.";
   if (sameExistingFile(p.outfile, p.dpiWitnessFile))
      throw std::string("Network and DPI diagnostic outputs must be different files.");
}

std::vector<std::string> buildConsensusFromManifest(
   Microarray_Set& data, const std::vector<int>& sources, Matrix& matrix, Parameter& p)
{
   const std::vector<std::string> paths = manifestPaths(p.consensusInputs);
   if (static_cast<std::size_t>(p.minRecurrence) > paths.size())
      throw std::string("Minimum recurrence exceeds the number of raw consensus inputs.");
   validateNativeOutputPaths(p, paths);
   std::vector<std::unique_ptr<RawReader> > readers;
   for (std::vector<std::string>::const_iterator path = paths.begin(); path != paths.end(); ++path)
      readers.push_back(std::unique_ptr<RawReader>(new RawReader(*path)));

   const std::size_t geneCount = data.markerset.size();
   std::set<int> uniqueSources;
   for (std::vector<int>::const_iterator source = sources.begin(); source != sources.end(); ++source)
      if (*source < 0 || static_cast<std::size_t>(*source) >= geneCount || !uniqueSources.insert(*source).second)
         throw std::string("Consensus requires distinct valid reconstructed source IDs.");
   std::unordered_map<std::string, int> geneIds;
   for (std::size_t i = 0; i < geneCount; ++i)
      geneIds.insert(std::make_pair(data.markerset[i].accnum.substr(1), static_cast<int>(i)));
   std::vector<double> sums(geneCount, 0.0);
   std::vector<std::uint32_t> counts(geneCount, 0), seen(geneCount, 0);
   matrix.nmv.assign(geneCount, AdjacencyRow());
   matrix.adjacencyRowsPresent.assign(geneCount, false);
   p.consensusRuns = paths.size();
   p.consensusUnionEdges = p.consensusCandidateEdges = 0;
   std::string row;
   for (std::vector<int>::const_iterator source = sources.begin(); source != sources.end(); ++source)
   {
      const std::string expected = data.markerset[*source].accnum.substr(1);
      std::fill(sums.begin(), sums.end(), 0.0);
      std::fill(counts.begin(), counts.end(), 0);
      std::fill(seen.begin(), seen.end(), 0);
      for (std::size_t run = 0; run < readers.size(); ++run)
      {
         RawReader& reader = *readers[run];
         if (!reader.next(row))
            reader.fail("missing requested source row '" + expected + "'");
         const std::size_t firstTab = row.find('\t');
         if (row.substr(0, firstTab) != expected)
            reader.fail("expected source row '" + expected + "' in supplied hub-list order");
         if (firstTab == std::string::npos)
            continue; // Explicit empty source row, never an absent row.
         std::size_t cursor = firstTab + 1;
         while (cursor < row.size())
         {
            const std::size_t split = row.find('\t', cursor);
            if (split == std::string::npos || split == cursor)
               reader.fail("malformed target/MI pair");
            const std::string target = row.substr(cursor, split - cursor);
            const std::size_t next = row.find('\t', split + 1);
            const std::string token = row.substr(split + 1, next == std::string::npos ? next : next - split - 1);
            const std::unordered_map<std::string, int>::const_iterator found = geneIds.find(target);
            if (found == geneIds.end())
               reader.fail("unknown target '" + target + "'");
            const int id = found->second;
            if (id == *source)
               reader.fail("self-edge for source '" + expected + "'");
            if (seen[id] == run + 1)
               reader.fail("duplicate target '" + target + "' within one resample row");
            seen[id] = static_cast<std::uint32_t>(run + 1);
            const double value = parseMi(token, reader);
            sums[id] += value;
            if (!std::isfinite(sums[id]))
               reader.fail("MI sum overflow");
            ++counts[id];
            if (next == std::string::npos)
            {
               cursor = row.size();
               break;
            }
            cursor = next + 1;
            if (cursor == row.size())
               reader.fail("trailing empty target/MI pair");
         }
         if (cursor == firstTab + 1)
            reader.fail("empty source rows must contain only their source ID");
      }
      matrix.markAdjacencyRow(*source);
      data.markerset[*source].isActive = true;
      AdjacencyRow& selected = matrix.nmv[*source];
      for (std::size_t target = 0; target < geneCount; ++target)
      {
         if (counts[target] == 0)
            continue;
         ++p.consensusUnionEdges;
         if (counts[target] >= static_cast<std::uint32_t>(p.minRecurrence))
         {
            selected.push_back(AdjacencyEdge(static_cast<int>(target), sums[target] / counts[target]));
            ++p.consensusCandidateEdges;
         }
      }
   }
   for (std::size_t run = 0; run < readers.size(); ++run)
      if (readers[run]->next(row))
         readers[run]->fail("extra source row after the supplied hub list");
   matrix.compactRows();
   return paths;
}
