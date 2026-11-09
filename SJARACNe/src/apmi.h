//------------------------------------------------------------------------------------
// Shared adaptive-partitioning mutual-information kernel.
//------------------------------------------------------------------------------------

#ifndef APMI_H__
#define APMI_H__

#include <algorithm>
#include <cstddef>
#include <vector>

const char *adaptivePartitionKernelSchema();

// One workspace belongs to one AP-MI evaluation stream. Every buffer range that can
// be read is reset or overwritten before each MI pair or partition.
class AdaptivePartitionWorkspace
{
public:
   AdaptivePartitionWorkspace()
      : observationCount(0), partitionLimit(0), poc(), kon(), poradi(), marg(),
        apor(), quadrant(), NN(), amarg() { }

   void initialize(int N, int M)
   {
      observationCount = N;
      partitionLimit   = M;

      poc.resize(M);
      kon.resize(M);
      poradi.resize(N);
      marg.resize(4 * static_cast<std::size_t>(M));
      apor.resize(N);
      quadrant.resize(N);
      NN.resize(4);
      amarg.resize(16);
   }

   bool matches(int N, int M) const
   {
      return observationCount == N && partitionLimit == M;
   }

   void resetEdge()
   {
      std::fill(poc.begin(), poc.end(), 1);
      std::fill(kon.begin(), kon.end(), observationCount);
      std::fill(marg.begin(), marg.end(), 0);

      for (int i = 0; i < observationCount; i++)
         poradi[i] = i + 1;

      marg[0] = marg[partitionLimit] = 1;
      marg[2 * partitionLimit] = marg[3 * partitionLimit] = observationCount;
   }

   int observationCount;
   int partitionLimit;
   std::vector<int> poc;
   std::vector<int> kon;
   std::vector<int> poradi;
   std::vector<int> marg;
   std::vector<int> apor;
   std::vector<unsigned char> quadrant;
   std::vector<int> NN;
   std::vector<int> amarg;
};

// Computes SJARACNe's exact adaptive-partitioning MI from two 1-based rank vectors.
// The optional diagnostic preserves the network executable's legacy maximum-npar
// reporting without introducing shared mutable state into the reusable kernel.
double computeAdaptivePartitionMI(const int *xranks, const int *yranks, int N,
                                  int nparLimit,
                                  AdaptivePartitionWorkspace& workspace,
                                  int *maximumObservedNpar = NULL);

#endif
