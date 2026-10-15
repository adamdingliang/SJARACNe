//------------------------------------------------------------------------------------
// Deterministic permutation-null generator for SJARACNe adaptive-partitioning MI.
//------------------------------------------------------------------------------------

#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include "apmi.h"

namespace
{

struct Options
{
   Options()
      : m(0), draws(0), drawsSpecified(false), enumerate(false), nparLimit(20),
        seed(1), format("tsv"), output("-") { }
   int m;
   std::uint64_t draws;
   bool drawsSpecified;
   bool enumerate;
   int nparLimit;
   std::uint32_t seed;
   std::string format;
   std::string output;
};

void usage(const char *program, std::ostream& out)
{
   out << "Usage: " << program
       << " --m INT --draws INT [--npar INT] [--seed INT]\n"
       << "       [--format tsv|binary] [--output FILE]\n"
       << "       " << program << " --m INT --enumerate [--npar INT]\n\n"
       << "Generate SJARACNe AP-MI values for identity ranks paired with uniformly\n"
       << "random rank permutations. TSV output defaults to stdout. Binary output is\n"
       << "raw little-endian float64 and requires FILE; metadata is written to FILE.meta.\n";
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

std::string requireValue(int argc, char **argv, int& index, const char *option)
{
   if (++index >= argc)
      throw std::runtime_error(std::string(option) + " requires a value");
   return argv[index];
}

Options parseOptions(int argc, char **argv)
{
   Options options;
   for (int i = 1; i < argc; i++)
   {
      const std::string argument = argv[i];
      if (argument == "--help" || argument == "-h")
      {
         usage(argv[0], std::cout);
         std::exit(0);
      }
      else if (argument == "--m")
      {
         const std::uint64_t value = parseUnsigned(requireValue(argc, argv, i, "--m"), "--m");
         if (value > static_cast<std::uint64_t>(INT_MAX))
            throw std::runtime_error("--m is too large");
         options.m = static_cast<int>(value);
      }
      else if (argument == "--draws")
      {
         options.draws = parseUnsigned(requireValue(argc, argv, i, "--draws"), "--draws");
         options.drawsSpecified = true;
      }
      else if (argument == "--enumerate")
         options.enumerate = true;
      else if (argument == "--npar")
      {
         const std::uint64_t value = parseUnsigned(requireValue(argc, argv, i, "--npar"), "--npar");
         if (value > static_cast<std::uint64_t>(INT_MAX))
            throw std::runtime_error("--npar is too large");
         options.nparLimit = static_cast<int>(value);
      }
      else if (argument == "--seed")
      {
         const std::uint64_t value = parseUnsigned(requireValue(argc, argv, i, "--seed"), "--seed");
         if (value > std::numeric_limits<std::uint32_t>::max())
            throw std::runtime_error("--seed must fit in an unsigned 32-bit integer");
         options.seed = static_cast<std::uint32_t>(value);
      }
      else if (argument == "--format")
         options.format = requireValue(argc, argv, i, "--format");
      else if (argument == "--output" || argument == "-o")
         options.output = requireValue(argc, argv, i, argument.c_str());
      else
         throw std::runtime_error("unknown option: " + argument);
   }

   if (options.m < 2)
      throw std::runtime_error("--m must be at least 2");
   if (options.enumerate)
   {
      if (options.drawsSpecified)
         throw std::runtime_error("--enumerate cannot be combined with --draws");
      if (options.m > 10)
         throw std::runtime_error("--enumerate supports --m at most 10");
      options.draws = 1;
      for (int value = 2; value <= options.m; value++)
         options.draws *= static_cast<std::uint64_t>(value);
   }
   else if (options.draws == 0)
      throw std::runtime_error("--draws must be positive unless --enumerate is used");
   if (options.nparLimit < 1)
      throw std::runtime_error("--npar must be positive");
   if (options.format != "tsv" && options.format != "binary")
      throw std::runtime_error("--format must be 'tsv' or 'binary'");
   if (options.format == "binary" && options.output == "-")
      throw std::runtime_error("--format binary requires --output FILE");
   return options;
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

void shuffleRanks(std::vector<int>& ranks, std::mt19937& generator)
{
   const std::size_t count = ranks.size();
   for (std::size_t position = 0; position + 1 < count; position++)
   {
      const std::uint32_t remaining = static_cast<std::uint32_t>(count - position);
      const std::size_t selected =
         position + static_cast<std::size_t>(randomBelow(generator, remaining));
      std::swap(ranks[position], ranks[selected]);
   }
}

class NullMiGenerator
{
public:
   NullMiGenerator(int m, int nparLimit, std::uint32_t seed, bool enumerate)
      : m_(m), nparLimit_(nparLimit), identity_(m), permutation_(m),
        generator_(seed), workspace_(), enumerate_(enumerate), draw_(0)
   {
      std::iota(identity_.begin(), identity_.end(), 1);
      permutation_ = identity_;
      workspace_.initialize(m_, nparLimit_);
   }

   double next()
   {
      if (enumerate_)
      {
         if (draw_ > 0 && !std::next_permutation(permutation_.begin(), permutation_.end()))
            throw std::runtime_error("complete permutation enumeration ended unexpectedly");
      }
      else
      {
         permutation_ = identity_;
         shuffleRanks(permutation_, generator_);
      }
      draw_++;
      return computeAdaptivePartitionMI(&identity_[0], &permutation_[0], m_,
                                        nparLimit_, workspace_);
   }

private:
   int m_;
   int nparLimit_;
   std::vector<int> identity_;
   std::vector<int> permutation_;
   std::mt19937 generator_;
   AdaptivePartitionWorkspace workspace_;
   bool enumerate_;
   std::uint64_t draw_;
};

void writeMetadata(const Options& options, std::ostream& out, const char *prefix,
                   const char *format)
{
   out << prefix << "format=" << format << "\n"
       << prefix << "kernel_schema=" << adaptivePartitionKernelSchema() << "\n"
       << prefix << "estimator=sjaracne-adaptive-partitioning\n"
       << prefix << "rng="
       << (options.enumerate ? "complete-lexicographic-permutation-enumeration-v1"
                             : "mt19937-rejection-fisher-yates-v1")
       << "\n"
       << prefix << "m=" << options.m << "\n"
       << prefix << "draws=" << options.draws << "\n"
       << prefix << "npar_limit=" << options.nparLimit << "\n"
       << prefix << "seed=" << options.seed << "\n";
}

void generateTsv(const Options& options, std::ostream& out)
{
   NullMiGenerator generator(options.m, options.nparLimit, options.seed,
                             options.enumerate);
   writeMetadata(options, out, "# ", "sjaracne-apmi-null-tsv-v1");
   out << "draw\tmi\n";
   out << std::setprecision(std::numeric_limits<double>::max_digits10);

   for (std::uint64_t draw = 0; draw < options.draws; draw++)
      out << draw << '\t' << generator.next() << '\n';
}

void appendDoubleLittleEndian(std::vector<char>& output, double value)
{
   static_assert(sizeof(double) == sizeof(std::uint64_t),
                 "binary null output requires 64-bit double");
   static_assert(std::numeric_limits<double>::is_iec559,
                 "binary null output requires IEEE-754 double");

   std::uint64_t bits = 0;
   std::memcpy(&bits, &value, sizeof(bits));
   for (int byte = 0; byte < 8; byte++)
      output.push_back(static_cast<char>((bits >> (8 * byte)) & 0xffu));
}

void generateBinary(const Options& options)
{
   std::ofstream output(options.output.c_str(), std::ios::out | std::ios::binary);
   if (!output.is_open())
      throw std::runtime_error("unable to open output file: " + options.output);

   NullMiGenerator generator(options.m, options.nparLimit, options.seed,
                             options.enumerate);
   const std::size_t recordsPerBuffer = 8192;
   std::vector<char> buffer;
   buffer.reserve(recordsPerBuffer * sizeof(double));

   for (std::uint64_t draw = 0; draw < options.draws; draw++)
   {
      appendDoubleLittleEndian(buffer, generator.next());
      if (buffer.size() == recordsPerBuffer * sizeof(double))
      {
         output.write(&buffer[0], static_cast<std::streamsize>(buffer.size()));
         buffer.clear();
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
   writeMetadata(options, metadata, "", "sjaracne-apmi-null-binary-v1");
   metadata << "dtype=float64\nbyte_order=little\nrecord_bytes=8\n"
            << "draw_index=implicit-zero-based\n";
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
      if (options.format == "binary")
         generateBinary(options);
      else if (options.output == "-")
         generateTsv(options, std::cout);
      else
      {
         std::ofstream output(options.output.c_str());
         if (!output.is_open())
            throw std::runtime_error("unable to open output file: " + options.output);
         generateTsv(options, output);
         if (!output.good())
            throw std::runtime_error("failed while writing output file: " + options.output);
      }
   }
   catch (const std::string& error)
   {
      std::cerr << argv[0] << ": " << error << '\n';
      return 1;
   }
   catch (const std::exception& error)
   {
      std::cerr << argv[0] << ": " << error.what() << '\n';
      return 1;
   }
   return 0;
}
