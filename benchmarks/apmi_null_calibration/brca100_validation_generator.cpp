// Held-out BRCA expression permutation null for SJARACNe AP-MI.
//
// This executable deliberately lives under benchmarks.  It reads expression
// values, selects gene pairs and observation subsets, permutes one gene within
// each selected subset, constructs the same unique ordinal ranks used by the
// network executable, and calls the shared production AP-MI kernel.

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "apmi.h"

namespace
{

struct Options
{
   Options()
      : draws(0), nparLimit(40), seed(20260820), expression(), output(),
        mValues() { }

   std::uint64_t draws;
   int nparLimit;
   std::uint32_t seed;
   std::string expression;
   std::string output;
   std::vector<int> mValues;
};

struct ExpressionMatrix
{
   int observations;
   std::vector<std::vector<double> > rows;
   std::vector<std::size_t> eligibleRows;
   std::size_t nonconstantRows;
   std::size_t tiedNonconstantRows;
};

void usage(const char *program, std::ostream& out)
{
   out << "Usage: " << program
       << " --expression FILE --m INT[,INT...] --draws INT --output FILE\n"
       << "       [--npar INT] [--seed INT]\n\n"
       << "Generate a held-out, marginal-preserving expression permutation null.\n"
       << "Only nonconstant rows without exact ties are eligible, matching the\n"
       << "unique-ordinal-rank null model. Output is raw little-endian float64 in\n"
       << "draw-major, listed-m order; metadata is written to FILE.meta.\n";
}

std::string requireValue(int argc, char **argv, int& index, const char *option)
{
   if (++index >= argc)
      throw std::runtime_error(std::string(option) + " requires a value");
   return argv[index];
}

std::uint64_t parseUnsigned(const std::string& text, const char *option)
{
   if (text.empty() || text[0] == '-')
      throw std::runtime_error(std::string(option) + " requires a nonnegative integer");

   errno = 0;
   char *end = NULL;
   const unsigned long long value = std::strtoull(text.c_str(), &end, 10);
   if (errno == ERANGE || end == text.c_str() || *end != '\0')
      throw std::runtime_error(std::string("invalid value for ") + option + ": " + text);
   return static_cast<std::uint64_t>(value);
}

std::vector<int> parseMValues(const std::string& text)
{
   std::vector<int> result;
   std::size_t begin = 0;
   while (begin <= text.size())
   {
      const std::size_t comma = text.find(',', begin);
      const std::string token =
         text.substr(begin, comma == std::string::npos ? std::string::npos : comma - begin);
      const std::uint64_t value = parseUnsigned(token, "--m");
      if (value < 2 || value > static_cast<std::uint64_t>(INT_MAX))
         throw std::runtime_error("every --m value must be in [2, INT_MAX]");
      result.push_back(static_cast<int>(value));
      if (comma == std::string::npos)
         break;
      begin = comma + 1;
   }

   std::sort(result.begin(), result.end());
   if (std::adjacent_find(result.begin(), result.end()) != result.end())
      throw std::runtime_error("--m values must be unique");
   return result;
}

Options parseOptions(int argc, char **argv)
{
   Options options;
   for (int i = 1; i < argc; ++i)
   {
      const std::string argument = argv[i];
      if (argument == "--help" || argument == "-h")
      {
         usage(argv[0], std::cout);
         std::exit(0);
      }
      else if (argument == "--expression")
         options.expression = requireValue(argc, argv, i, "--expression");
      else if (argument == "--output" || argument == "-o")
         options.output = requireValue(argc, argv, i, argument.c_str());
      else if (argument == "--m")
         options.mValues = parseMValues(requireValue(argc, argv, i, "--m"));
      else if (argument == "--draws")
         options.draws = parseUnsigned(requireValue(argc, argv, i, "--draws"), "--draws");
      else if (argument == "--npar")
      {
         const std::uint64_t value =
            parseUnsigned(requireValue(argc, argv, i, "--npar"), "--npar");
         if (value < 1 || value > static_cast<std::uint64_t>(INT_MAX))
            throw std::runtime_error("--npar must be in [1, INT_MAX]");
         options.nparLimit = static_cast<int>(value);
      }
      else if (argument == "--seed")
      {
         const std::uint64_t value =
            parseUnsigned(requireValue(argc, argv, i, "--seed"), "--seed");
         if (value > std::numeric_limits<std::uint32_t>::max())
            throw std::runtime_error("--seed must fit in an unsigned 32-bit integer");
         options.seed = static_cast<std::uint32_t>(value);
      }
      else
         throw std::runtime_error("unknown option: " + argument);
   }

   if (options.expression.empty())
      throw std::runtime_error("--expression is required");
   if (options.output.empty())
      throw std::runtime_error("--output is required");
   if (options.mValues.empty())
      throw std::runtime_error("--m is required");
   if (options.draws == 0)
      throw std::runtime_error("--draws must be positive");
   return options;
}

void stripCarriageReturn(std::string& line)
{
   if (!line.empty() && line[line.size() - 1] == '\r')
      line.erase(line.size() - 1);
}

std::vector<std::string> splitTabs(const std::string& line)
{
   std::vector<std::string> fields;
   std::size_t begin = 0;
   while (true)
   {
      const std::size_t tab = line.find('\t', begin);
      fields.push_back(line.substr(begin, tab == std::string::npos
                                           ? std::string::npos : tab - begin));
      if (tab == std::string::npos)
         break;
      begin = tab + 1;
   }
   return fields;
}

double parseFiniteDouble(const std::string& text, std::size_t line,
                         std::size_t observation)
{
   errno = 0;
   char *end = NULL;
   const double value = std::strtod(text.c_str(), &end);
   while (end != NULL && *end != '\0' &&
          std::isspace(static_cast<unsigned char>(*end)))
      ++end;

   if (errno == ERANGE || end == text.c_str() || end == NULL || *end != '\0' ||
       !std::isfinite(value))
   {
      std::ostringstream message;
      message << "invalid or non-finite expression value at line " << line
              << ", observation " << observation;
      throw std::runtime_error(message.str());
   }
   return value;
}

ExpressionMatrix readExpression(const std::string& filename)
{
   std::ifstream input(filename.c_str());
   if (!input.is_open())
      throw std::runtime_error("unable to open expression file: " + filename);

   std::string line;
   if (!std::getline(input, line))
      throw std::runtime_error("expression file is empty: " + filename);
   stripCarriageReturn(line);
   const std::vector<std::string> header = splitTabs(line);
   if (header.size() < 4)
      throw std::runtime_error("expression matrix must contain two metadata columns and at least two observations");
   if (header.size() - 2 > static_cast<std::size_t>(INT_MAX))
      throw std::runtime_error("expression observation count exceeds INT_MAX");

   ExpressionMatrix matrix;
   matrix.observations = static_cast<int>(header.size() - 2);
   matrix.nonconstantRows = 0;
   matrix.tiedNonconstantRows = 0;

   std::size_t lineNumber = 1;
   while (std::getline(input, line))
   {
      ++lineNumber;
      stripCarriageReturn(line);
      if (line.empty())
      {
         std::ostringstream message;
         message << "blank expression row at line " << lineNumber;
         throw std::runtime_error(message.str());
      }

      const std::vector<std::string> fields = splitTabs(line);
      if (fields.size() != header.size())
      {
         std::ostringstream message;
         message << "incorrect expression dimensions at line " << lineNumber
                 << ": expected " << header.size() << " fields, found "
                 << fields.size();
         throw std::runtime_error(message.str());
      }

      std::vector<double> row(matrix.observations);
      for (int observation = 0; observation < matrix.observations; ++observation)
         row[observation] = parseFiniteDouble(fields[observation + 2], lineNumber,
                                              observation + 1);

      std::vector<double> sorted = row;
      std::sort(sorted.begin(), sorted.end());
      const bool nonconstant = sorted.front() != sorted.back();
      const bool hasTies =
         std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end();

      matrix.rows.push_back(row);
      if (nonconstant)
      {
         ++matrix.nonconstantRows;
         if (hasTies)
            ++matrix.tiedNonconstantRows;
         else
            matrix.eligibleRows.push_back(matrix.rows.size() - 1);
      }
   }

   if (!input.eof())
      throw std::runtime_error("failed while reading expression file: " + filename);
   if (matrix.eligibleRows.size() < 2)
      throw std::runtime_error("fewer than two nonconstant, tie-free expression rows are available");
   return matrix;
}

std::uint32_t randomBelow(std::mt19937& generator, std::uint32_t bound)
{
   const std::uint32_t threshold = static_cast<std::uint32_t>(-bound) % bound;
   while (true)
   {
      const std::uint32_t value = static_cast<std::uint32_t>(generator());
      if (value >= threshold)
         return value % bound;
   }
}

template <typename T>
void shuffleVector(std::vector<T>& values, std::mt19937& generator)
{
   for (std::size_t position = 0; position + 1 < values.size(); ++position)
   {
      const std::uint32_t remaining =
         static_cast<std::uint32_t>(values.size() - position);
      const std::size_t selected =
         position + static_cast<std::size_t>(randomBelow(generator, remaining));
      std::swap(values[position], values[selected]);
   }
}

struct RankedObservation
{
   RankedObservation(double inValue, int inPosition)
      : value(inValue), position(inPosition) { }
   double value;
   int position;
};

struct RankedObservationLess
{
   bool operator()(const RankedObservation& left,
                   const RankedObservation& right) const
   {
      // This is intentionally the same comparator as BuildRankCache.
      return left.value < right.value ||
             (left.value == right.value && left.position < right.position);
   }
};

void makeRanks(const std::vector<double>& values, std::vector<int>& ranks)
{
   std::vector<RankedObservation> observations;
   observations.reserve(values.size());
   for (std::size_t position = 0; position < values.size(); ++position)
      observations.push_back(RankedObservation(values[position],
                                               static_cast<int>(position)));
   std::sort(observations.begin(), observations.end(), RankedObservationLess());

   ranks.resize(values.size());
   for (std::size_t rank = 0; rank < observations.size(); ++rank)
      ranks[observations[rank].position] = static_cast<int>(rank) + 1;
}

void appendDoubleLittleEndian(std::vector<char>& output, double value)
{
   static_assert(sizeof(double) == sizeof(std::uint64_t),
                 "binary output requires 64-bit double");
   static_assert(std::numeric_limits<double>::is_iec559,
                 "binary output requires IEEE-754 double");
   std::uint64_t bits = 0;
   std::memcpy(&bits, &value, sizeof(bits));
   for (int byte = 0; byte < 8; ++byte)
      output.push_back(static_cast<char>((bits >> (8 * byte)) & 0xffu));
}

std::string joinMValues(const std::vector<int>& values)
{
   std::ostringstream output;
   for (std::size_t i = 0; i < values.size(); ++i)
   {
      if (i != 0)
         output << ',';
      output << values[i];
   }
   return output.str();
}

void generate(const Options& options, const ExpressionMatrix& expression)
{
   if (options.mValues.back() > expression.observations)
      throw std::runtime_error("an --m value exceeds the expression observation count");

   std::ofstream output(options.output.c_str(), std::ios::out | std::ios::binary);
   if (!output.is_open())
      throw std::runtime_error("unable to open output file: " + options.output);

   std::mt19937 generator(options.seed);
   std::vector<int> sampleOrder(expression.observations);
   std::vector<int> yOrder;
   std::vector<double> xValues;
   std::vector<double> yValues;
   std::vector<int> xRanks;
   std::vector<int> yRanks;
   std::vector<AdaptivePartitionWorkspace> workspaces(options.mValues.size());
   for (std::size_t i = 0; i < options.mValues.size(); ++i)
      workspaces[i].initialize(options.mValues[i], options.nparLimit);

   const std::size_t recordsPerBuffer = 8192;
   std::vector<char> buffer;
   buffer.reserve(recordsPerBuffer * sizeof(double));

   if (expression.eligibleRows.size() >
       static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
      throw std::runtime_error("eligible gene count exceeds the RNG's uint32 range");
   const std::uint32_t geneCount =
      static_cast<std::uint32_t>(expression.eligibleRows.size());
   for (std::uint64_t draw = 0; draw < options.draws; ++draw)
   {
      const std::uint32_t firstIndex = randomBelow(generator, geneCount);
      std::uint32_t secondIndex = randomBelow(generator, geneCount - 1);
      if (secondIndex >= firstIndex)
         ++secondIndex;
      const std::vector<double>& xRow =
         expression.rows[expression.eligibleRows[firstIndex]];
      const std::vector<double>& yRow =
         expression.rows[expression.eligibleRows[secondIndex]];

      std::iota(sampleOrder.begin(), sampleOrder.end(), 0);
      shuffleVector(sampleOrder, generator);

      for (std::size_t mIndex = 0; mIndex < options.mValues.size(); ++mIndex)
      {
         const int m = options.mValues[mIndex];
         yOrder.assign(sampleOrder.begin(), sampleOrder.begin() + m);
         shuffleVector(yOrder, generator);
         xValues.resize(m);
         yValues.resize(m);
         for (int position = 0; position < m; ++position)
         {
            xValues[position] = xRow[sampleOrder[position]];
            yValues[position] = yRow[yOrder[position]];
         }

         makeRanks(xValues, xRanks);
         makeRanks(yValues, yRanks);
         const double mi = computeAdaptivePartitionMI(
            &xRanks[0], &yRanks[0], m, options.nparLimit, workspaces[mIndex]);
         const double tolerance = 64 * std::numeric_limits<double>::epsilon();
         if (!std::isfinite(mi) || mi < -tolerance ||
             mi > std::log(static_cast<double>(m)) + tolerance)
            throw std::runtime_error("shared AP-MI kernel produced an invalid value");

         appendDoubleLittleEndian(buffer, mi);
         if (buffer.size() >= recordsPerBuffer * sizeof(double))
         {
            output.write(&buffer[0], static_cast<std::streamsize>(buffer.size()));
            buffer.clear();
         }
      }
   }

   if (!buffer.empty())
      output.write(&buffer[0], static_cast<std::streamsize>(buffer.size()));
   output.close();
   if (!output)
      throw std::runtime_error("failed while writing output file: " + options.output);

   const std::string metadataName = options.output + ".meta";
   std::ofstream metadata(metadataName.c_str());
   if (!metadata.is_open())
      throw std::runtime_error("unable to open metadata file: " + metadataName);
   metadata << "format=sjaracne-expression-permutation-null-binary-v1\n"
            << "kernel_schema=" << adaptivePartitionKernelSchema() << "\n"
            << "estimator=sjaracne-adaptive-partitioning\n"
            << "sampling_null=held-out-expression-within-subset-permutation\n"
            << "rank_policy=unique-ordinal-value-then-selected-position\n"
            << "expression=" << options.expression << "\n"
            << "observations=" << expression.observations << "\n"
            << "genes_total=" << expression.rows.size() << "\n"
            << "genes_nonconstant=" << expression.nonconstantRows << "\n"
            << "genes_tie_free_eligible=" << expression.eligibleRows.size() << "\n"
            << "genes_nonconstant_excluded_for_ties="
            << expression.tiedNonconstantRows << "\n"
            << "m_values=" << joinMValues(options.mValues) << "\n"
            << "draws=" << options.draws << "\n"
            << "npar_limit=" << options.nparLimit << "\n"
            << "seed=" << options.seed << "\n"
            << "rng=mt19937-rejection-fisher-yates-v1\n"
            << "gene_pair_policy=uniform-distinct-eligible-rows\n"
            << "observation_policy=uniform-without-replacement-prefix\n"
            << "permutation_policy=uniform-within-selected-subset\n"
            << "dtype=float64\nbyte_order=little\nrecord_bytes=8\n"
            << "record_layout=draw-major-m-values-order\n";
   metadata.close();
   if (!metadata)
      throw std::runtime_error("failed while writing metadata file: " + metadataName);
}

} // namespace

int main(int argc, char **argv)
{
   try
   {
      const Options options = parseOptions(argc, argv);
      const ExpressionMatrix expression = readExpression(options.expression);
      generate(options, expression);
      return 0;
   }
   catch (const std::exception& error)
   {
      std::cerr << "ERROR: " << error.what() << std::endl;
      usage(argv[0], std::cerr);
      return 2;
   }
}
