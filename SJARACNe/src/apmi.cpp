//------------------------------------------------------------------------------------
// Shared adaptive-partitioning mutual-information kernel.
//------------------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <string>
#include "apmi.h"

const char *adaptivePartitionKernelSchema()
{
   return "sjaracne-apmi-v1";
}

double computeAdaptivePartitionMI(const int *xranks, const int *yranks, int N,
                                  int nparLimit,
                                  AdaptivePartitionWorkspace& workspace,
                                  int *maximumObservedNpar)
{
   const int M = nparLimit;

   int npar = 1;
   int run  = 0;

   if (maximumObservedNpar != NULL)
      *maximumObservedNpar = 1;

   double xcor = 0.0;

   if (!workspace.matches(N, M))
      throw std::string("Adaptive-partitioning workspace dimensions do not match MI input.");

   workspace.resetEdge();

   std::vector<int>& poc    = workspace.poc;
   std::vector<int>& kon    = workspace.kon;
   std::vector<int>& poradi = workspace.poradi;
   std::vector<int>& marg   = workspace.marg;
   std::vector<int>& apor   = workspace.apor;
   std::vector<unsigned char>& quadrant = workspace.quadrant;
   std::vector<int>& NN     = workspace.NN;
   std::vector<int>& amarg  = workspace.amarg;

   while (npar > 0)
   {
      run++;

      int np   = npar - 1;
      int apoc = poc[np];
      int akon = kon[np];
      int Nex  = akon - apoc + 1;

      for (int i = 0; i < Nex; i++)
         apor[i] = poradi[apoc + i - 1];

      int ave1 = std::floor((marg[np] + marg[np + 2 * M]) / 2);
      int ave2 = std::floor((marg[np + M] + marg[np + 3 * M]) / 2);

      std::fill(NN.begin(), NN.end(), 0);

      for (int i = 0; i < Nex; i++)
      {
         int k = apor[i] - 1;

         int j = (xranks[k] <= ave1 ? 0 : 2) + (yranks[k] <= ave2 ? 0 : 1);

         quadrant[i] = static_cast<unsigned char>(j);
         NN[j]++;
      }

      double c   = Nex / 4.0;
      double sum = 0.0;

      for (int i = 0; i < 4; i++)
      {
         double d = NN[i] - c;
         sum += d * d;
      }

      double tst = 4 * sum / Nex;

      if (tst > 7.8 || run == 1)
      {
         amarg[ 0] = amarg[ 1] = marg[np];
         amarg[ 2] = amarg[ 3] = ave1 + 1;
         amarg[ 4] = amarg[ 6] = marg[np + M];
         amarg[ 5] = amarg[ 7] = ave2 + 1;
         amarg[ 8] = amarg[ 9] = ave1;
         amarg[10] = amarg[11] = marg[np + 2 * M];
         amarg[12] = amarg[14] = ave2;
         amarg[13] = amarg[15] = marg[np + 3 * M];

         // Pack every non-leaf child in one stable pass. These offsets preserve
         // both quadrant and observation order.
         int writePosition[4] = { -1, -1, -1, -1 };
         int nextPosition = apoc - 1;

         for (int i = 0; i < 4; i++)
            if (NN[i] > 2)
            {
               writePosition[i] = nextPosition;
               nextPosition += NN[i];
            }

         for (int i = 0; i < Nex; i++)
         {
            int childQuadrant = quadrant[i];

            if (writePosition[childQuadrant] >= 0)
               poradi[writePosition[childQuadrant]++] = apor[i];
         }

         npar--;

         for (int i = 0; i < 4; i++)
            if (NN[i] > 2)
            {
               if (++npar > M)
                  throw std::string("Exceeded npar limit!");

               if (maximumObservedNpar != NULL && npar > *maximumObservedNpar)
                  *maximumObservedNpar = npar;

               akon = apoc + NN[i] - 1;

               int np = npar - 1;

               poc[np] = apoc;
               kon[np] = akon;

               for (int j = 0; j < 4; j++)
                  marg[np + j * M] = amarg[i + 4 * j];

               apoc = akon + 1;
            }
            else if (NN[i] > 0)
            {
               double Nx = amarg[i +  8] - amarg[i] + 1;
               double Ny = amarg[i + 12] - amarg[i + 4] + 1;

               xcor += NN[i] * std::log(NN[i] / (Nx * Ny));
            }
      }
      else
      {
         double Nx = marg[np + 2 * M] - marg[np] + 1;
         double Ny = marg[np + 3 * M] - marg[np + M] + 1;

         xcor += Nex * std::log(Nex / (Nx * Ny));

         npar--;
      }
   }

   return (xcor / N + std::log(N));
}
