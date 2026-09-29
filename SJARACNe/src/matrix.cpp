//------------------------------------------------------------------------------------
// Copyright (C) 2003  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017
//------------------------------------------------------------------------------------

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <ios>
#include <iterator>
#include <limits>
#include <new>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include "apmi.h"
#include "matrix.h"
#include "util.h"

int maxNpar = 0; // maximum observed value of npar

//------------------------------------------------------------------------------------

class ArrayValuePair
{
public:
   ArrayValuePair(int id, double v)
      : arrayId(id), value(v) { }

   int arrayId;  // getId(), set()
   double value; // getValue(), set()
};

typedef std::binary_function<ArrayValuePair, ArrayValuePair, bool>
   ArrayValuePairBinaryFunction;

class SortIncreasing_ArrayValuePair : ArrayValuePairBinaryFunction
{
public:
   bool operator()(const ArrayValuePair& a, const ArrayValuePair& b) const
   {
      return (a.value < b.value);
   }
};

class SortDecreasing_ArrayValuePair : ArrayValuePairBinaryFunction
{
public:
   bool operator()(const ArrayValuePair& a, const ArrayValuePair& b) const
   {
      return (a.value > b.value);
   }
};

//------------------------------------------------------------------------------------

class AdjacencyTargetLess
{
public:
   bool operator()(const AdjacencyEdge& left,
                   const AdjacencyEdge& right) const
   {
      return left.target < right.target;
   }

   bool operator()(const AdjacencyEdge& edge, int target) const
   {
      return edge.target < target;
   }
};

static void normalizeImportedRow(AdjacencyRow& row)
{
   if (row.size() < 2)
      return;

   bool targetSorted = true;
   for (std::size_t i = 1; i < row.size(); ++i)
      if (row[i - 1].target > row[i].target)
      {
         targetSorted = false;
         break;
      }

   // Stable ordering is required so the last retained occurrence of a
   // duplicate edge remains authoritative, matching the legacy map update.
   if (!targetSorted)
      std::stable_sort(row.begin(), row.end(), AdjacencyTargetLess());

   std::size_t write = 0;
   for (std::size_t first = 0; first < row.size();)
   {
      std::size_t last = first + 1;
      while (last < row.size() && row[last].target == row[first].target)
         ++last;

      row[write++] = row[last - 1];
      first = last;
   }

   row.resize(write);
}

static AdjacencyRow::iterator findAdjacencyEdge(AdjacencyRow& row, int target)
{
   AdjacencyRow::iterator edge =
      std::lower_bound(row.begin(), row.end(), target, AdjacencyTargetLess());

   return edge != row.end() && edge->target == target ? edge : row.end();
}

class DpiNeighborIndex
{
public:
   explicit DpiNeighborIndex(const Matrix& matrix)
      : rows(matrix.nmv), incomingForEmptyRows()
   {
      std::size_t nodeCount = std::max(matrix.nmv.size(),
                                      matrix.adjacencyRowsPresent.size());

      for (AdjacencyRows::const_iterator row = rows.begin(); row != rows.end(); ++row)
         if (!row->empty())
         {
            if (row->front().target < 0)
               throw std::string(
                  "Internal error: negative gene ID in adjacency matrix.");

            int largestTarget = row->back().target;
            nodeCount = std::max(nodeCount,
                                 static_cast<std::size_t>(largestTarget) + 1);
         }

      if (rows.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max()))
         throw std::string(
            "Internal error: adjacency matrix exceeds the supported gene-ID range.");

      incomingForEmptyRows.resize(nodeCount);

      for (std::size_t source = 0; source < rows.size(); ++source)
         for (AdjacencyRow::const_iterator edge = rows[source].begin();
              edge != rows[source].end(); ++edge)
         {
            int target = edge->target;
            bool targetHasRow = static_cast<std::size_t>(target) < rows.size() &&
                                !rows[target].empty();

            if (!targetHasRow)
               incomingForEmptyRows[target].push_back(
                  ArrayValuePair(static_cast<int>(source),
                                 edge->mutinfo));
         }
   }

   std::size_t nodeCount() const
   {
      return incomingForEmptyRows.size();
   }

   bool hasDirectRow(int geneId) const
   {
      return geneId >= 0 && static_cast<std::size_t>(geneId) < rows.size() &&
             !rows[geneId].empty();
   }

   const AdjacencyRow& directRow(int geneId) const
   {
      return rows[geneId];
   }

   const std::vector<ArrayValuePair>& incomingRow(int geneId) const
   {
      return incomingForEmptyRows[geneId];
   }

   std::size_t effectiveDegree(int geneId) const
   {
      if (hasDirectRow(geneId))
         return directRow(geneId).size();

      if (geneId >= 0 &&
          static_cast<std::size_t>(geneId) < incomingForEmptyRows.size())
         return incomingRow(geneId).size();

      return 0;
   }

private:
   const AdjacencyRows& rows;
   std::vector<std::vector<ArrayValuePair> > incomingForEmptyRows;
};

//------------------------------------------------------------------------------------

class RankedObservation
{
public:
   RankedObservation()
      : value(), position() { }

   RankedObservation(double inValue, int inPosition)
      : value(inValue), position(inPosition) { }

   double value;
   int position;
};

class Sort_RankedObservation
{
public:
   bool operator()(const RankedObservation& a, const RankedObservation& b) const
   {
      return (a.value < b.value ||
              (a.value == b.value && a.position < b.position));
   }
};

//------------------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, const Marker& m)
{
   out << "(" << m.label << "; " << m.accnum << "; " << m.idnum << "; "
       << (m.isActive ? "Active; " : "Not Active; ")
       << (m.isControl ? "Control" : "Not Control" ) << ")";
   return out;
}

//------------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const Probe& p)
{
   out << "(" << p.value << ", " << p.pvalue << ")";
   return out;
}

//------------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const Microarray_Set& ms)
{
   const char *delim = "\t";

   int numMarkers     = ms.markerset.size();
   int numMicroarrays = ms.uarrays.size();

   if (ms.header.size() > 0)
      std::copy(ms.header.begin(), ms.header.end(),
                std::ostream_iterator<std::string>(out, delim));
   else
   {
      // make up a col header
      out << "Id" << delim << "Desc" << delim;

      int col = 0;
      while (col < numMarkers - 1)
         out << "exp" << col++ << delim;

      out << "exp" << col;
   }

   out << "\n";

   for (int i = 0; i < numMarkers; i++)
   {
      out << ms.markerset[i].accnum << "\t" << ms.markerset[i].label;

      for (int j = 0; j < numMicroarrays; j++)
         out << delim << ms.uarrays[j][i].value;

      if (i < numMarkers - 1)
         out << "\n";
   }

   return out;
}

//------------------------------------------------------------------------------------

void Matrix::saveNode(int i, int j, double mi)
{
   if (!std::isfinite(mi))
      throw std::string("Refusing to store a non-finite MI value.");

   AdjacencyRow& row = nmv[i];

   if (row.empty() || row.back().target < j)
   {
      row.push_back(AdjacencyEdge(j, mi));
      return;
   }

   if (row.back().target == j)
   {
      row.back().mutinfo = mi;
      return;
   }

   AdjacencyRow::iterator edge =
      std::lower_bound(row.begin(), row.end(), j, AdjacencyTargetLess());

   if (edge == row.end() || edge->target != j)
      row.insert(edge, AdjacencyEdge(j, mi));
   else
      edge->mutinfo = mi;
}

//------------------------------------------------------------------------------------

void Matrix::read(Microarray_Set& data, const Parameter& p)
{
   std::ifstream in(p.adjfile.c_str());
   if (!in.is_open())
      throw "Unable to open " + p.adjfile;

   read(in, data, p);

   in.close();
}

//------------------------------------------------------------------------------------

static double parseMutualInformation(const std::string& text,
                                     const std::string& source,
                                     const std::string& target,
                                     const std::string& filename)
{
   const char *begin = text.c_str();
   char *end = NULL;
   errno = 0;

   double mi = std::strtod(begin, &end);
   const bool converted = (end != NULL && end != begin);

   while (end != NULL && *end != '\0' &&
          std::isspace(static_cast<unsigned char>(*end)))
      end++;

   std::ostringstream s;

   if (!converted || end == NULL || *end != '\0' || errno == ERANGE)
   {
      s << "Invalid MI value \"" << text << "\" for edge " << source
        << " -> " << target << " in adjacency file \"" << filename << "\".";
      throw s.str();
   }

   if (!std::isfinite(mi))
   {
      s << "Non-finite MI value \"" << text << "\" for edge " << source
        << " -> " << target << " in adjacency file \"" << filename << "\".";
      throw s.str();
   }

   return mi;
}

//------------------------------------------------------------------------------------

void Matrix::read(std::istream& in, Microarray_Set& data, const Parameter& p)
{
   // The expression matrix defines the complete gene-ID space.  Pre-sizing
   // prevents sparse -j inputs from leaving requested rows out of bounds.
   nmv.assign(data.markerset.size(), AdjacencyRow());
   adjacencyRowsPresent.assign(data.markerset.size(), false);

   std::string line;

   std::getline(in, line);
   while (line.length() > 0 && line[0] == '>')
      std::getline(in, line);

   while (in.good())
   {
      std::istringstream sin(line);
      std::string label;

      std::getline(sin, label, '\t');
      label = "_" + label;

      int geneId1 = data.getProbeId(label);
      if (geneId1 == -1)
         throw "Cannot find marker: " + label + " in the ADJ file!";

      data.markerset[geneId1].isActive = true;
      markAdjacencyRow(geneId1);

      std::getline(sin, label, '\t');
      label = "_" + label;

      while (sin.good())
      {
         std::string value;

         std::getline(sin, value, '\t');

         const std::string source = data.markerset[geneId1].accnum.substr(1);
         const std::string target = label.substr(1);
         double mi = parseMutualInformation(value, source, target, p.adjfile);

         if (mi >= p.threshold)
         {
            int geneId2 = data.getProbeId(label);
            if (geneId2 == -1)
               throw "Cannot find marker: " + label + " in the ADJ file!";

            // Imported adjacency rows can be unordered or repeated. Append in
            // file order, then stable-sort and deduplicate each row once after
            // parsing instead of shifting a vector on every unordered edge.
            nmv[geneId1].push_back(AdjacencyEdge(geneId2, mi));
         }

         std::getline(sin, label, '\t');
         label = "_" + label;
      }

      std::getline(in, line);
   }

   for (AdjacencyRows::iterator row = nmv.begin(); row != nmv.end(); ++row)
      normalizeImportedRow(*row);

   compactRows();
}

//------------------------------------------------------------------------------------

void Matrix::markAdjacencyRow(int geneId)
{
   if (geneId < 0)
      throw std::string("Internal error: cannot mark a negative adjacency row.");

   if (static_cast<std::size_t>(geneId) >= adjacencyRowsPresent.size())
      adjacencyRowsPresent.resize(geneId + 1, false);

   adjacencyRowsPresent[geneId] = true;
}

//------------------------------------------------------------------------------------

bool Matrix::hasAdjacencyRow(int geneId) const
{
   return geneId >= 0 &&
          static_cast<std::size_t>(geneId) < adjacencyRowsPresent.size() &&
          adjacencyRowsPresent[geneId];
}

//------------------------------------------------------------------------------------

bool Matrix::hasEnoughSourceRowsForDpi() const
{
   int sourceRows = 0;

   for (std::vector<bool>::const_iterator row = adjacencyRowsPresent.begin();
        row != adjacencyRowsPresent.end(); ++row)
      if (*row && ++sourceRows == 2)
         return true;

   return false;
}

//------------------------------------------------------------------------------------

void Matrix::writeGeneLine(std::ostream& out, const Microarray_Set& data, int geneId)
{
   AdjacencyRow& row = nmv[geneId];

   const Marker& marker = data.markerset[geneId];

   if (row.size() == 0 &&
       (!writeEmptyGenes || !marker.isActive && !marker.isControl))
      return;

   out << marker.accnum.substr(1);

   for (AdjacencyRow::iterator edge = row.begin(); edge != row.end(); ++edge)
   {
      int id = edge->target;

      if (writeTriangular && id <= geneId ||
          !writeReduced && edge->intermediate >= 0)
         continue;

      const Marker& marker = data.markerset[id];

      out << "\t" << marker.accnum.substr(1);

      if (edge->intermediate >= 0)
         out << "." << edge->intermediate;

      out << "\t" << edge->mutinfo;
   }

   out << std::endl;
}

//------------------------------------------------------------------------------------

void Matrix::writeGeneList(const Microarray_Set& data, const std::string& name,
                           int probeId)
{
   std::string filename = name + ".adj";

   std::ofstream out(filename.c_str());
   if (!out.is_open())
      throw "Unable to open " + filename;

   std::cout << "Writing gene list: "<< filename << std::endl;

   AdjacencyRow& row = nmv[probeId];

   const Marker& marker = data.markerset[probeId];

   if (row.size() == 0 &&
       (!writeEmptyGenes || !marker.isActive && !marker.isControl))
      return;

   for (AdjacencyRow::iterator edge = row.begin(); edge != row.end(); ++edge)
   {
      int id = edge->target;

      const Marker& marker = data.markerset[id];

      if (id != probeId)
         out << id << "\t" << marker.accnum << "t" << edge->mutinfo << std::endl;
   }

   out.close();
}

//------------------------------------------------------------------------------------

void Matrix::write(std::ostream& out, const Microarray_Set& data,
                   const std::vector<int>& ids)
{
   int numIds = ids.size();

   if (numIds == 0)
   {
      int numNodeMaps = nmv.size();

      for (int id = 0; id < numNodeMaps; id++)
         writeGeneLine(out, data, id);
   }
   else
      for (int i = 0; i < numIds; i++)
         writeGeneLine(out, data, ids[i]);
}

//------------------------------------------------------------------------------------

void Matrix::write(const Microarray_Set& data, const std::vector<int>& ids,
                   const Parameter& p, bool writeFull)
{
   if (!writeFull)
      return;

   std::ofstream out(p.outfile.c_str());
   if (!out.is_open())
      throw "Unable to open " + p.outfile;

   std::cout << "Writing matrix: " << p.outfile << std::endl;

   out << ">  Input file      " << p.infile     << std::endl;
   //out << ">  ADJ file        " << p.adjfile    << std::endl;
   //out << ">  Output file     " << p.outfile    << std::endl;
   //out << ">  Kernel width    " << p.sigma      << std::endl;
   out << ">  MI threshold    " << p.threshold  << std::endl;
   out << ">  MI P-value      " << p.pvalue     << std::endl;
   if (!p.nullModelFormat.empty())
   {
      out << ">  MI threshold method " << p.thresholdMethod << std::endl;
      out << ">  AP-MI null model file " << p.nullModelFile << std::endl;
      out << ">  AP-MI null model format " << p.nullModelFormat << std::endl;
      out << ">  AP-MI kernel schema " << p.nullModelKernelSchema << std::endl;
      out << ">  AP-MI estimator " << p.nullModelEstimator << std::endl;
      out << ">  AP-MI tail model " << p.nullModelTailModel << std::endl;
      out << ">  AP-MI calibrator schema " << p.nullModelCalibratorSchema
          << std::endl;
      out << ">  AP-MI calibrator SHA256 " << p.nullModelCalibratorSha256
          << std::endl;
      out << ">  AP-MI null model m " << p.nullModelM << std::endl;
      out << ">  AP-MI null model Npar " << p.nullModelNparLimit << std::endl;
      out << ">  AP-MI supported p range [" << p.nullModelSupportedPMin
          << "," << p.nullModelSupportedPMax << "]" << std::endl;
      if (p.nullModelHasValidatedPMin)
         out << ">  AP-MI validated p min " << p.nullModelValidatedPMin
             << std::endl;
      if (p.nullModelHasValidatedPMin)
         out << ">  AP-MI validated p max " << p.nullModelValidatedPMax
             << std::endl;
      else
         out << ">  AP-MI validated p min none" << std::endl;
      out << ">  AP-MI cutoff tail extrapolated "
          << (p.nullModelTailExtrapolated ? "yes" : "no") << std::endl;
      out << ">  AP-MI generator SHA256 " << p.nullModelGeneratorSha256
          << std::endl;
      out << ">  AP-MI fit values SHA256 " << p.nullModelFitValuesSha256
          << std::endl;
      out << ">  AP-MI validation values SHA256 "
          << p.nullModelValidationValuesSha256 << std::endl;
   }
   out << ">  DPI tolerance   " << p.eps        << std::endl;
   //out << ">  Correction      " << p.correction << std::endl;
   out << ">  Subnetwork file " << p.subnetfile << std::endl;

   if (p.samplingMethod != "")
   {
      out << ">  Sampling method " << p.samplingMethod << std::endl;
      if (p.subsampleSpec != "")
         out << ">  Sampling request " << p.subsampleSpec << std::endl;
      else
         out << ">  Sampling request legacy -r " << p.sample << std::endl;
      out << ">  Eligible observations " << p.samplingPopulation << std::endl;
      out << (p.subsampleSpec != "" ? ">  Sampled observations "
                                     : ">  Resample draws ")
          << p.samplingSize << std::endl;
   }

   //out << ">  Hub probe       " << p.hub        << std::endl;
   //out << ">  Control probe   " << p.controlId  << std::endl;
   //out << ">  Condition       " << p.condition  << std::endl;
   //out << ">  Percentage      " << p.percent    << std::endl;
   //out << ">  TF annotation   " << p.annotfile  << std::endl;
   //out << ">  Filter mean     " << p.mean       << std::endl;
   //out << ">  Filter CV       " << p.cv         << std::endl;

   write(out, data, ids);

   out.close();

   std::cout << "Maximum observed npar: " << maxNpar << std::endl;
}

//------------------------------------------------------------------------------------

void Matrix::createEntries(int numEntries)
{
   // create entries for all rows in the matrix

   for (int i = 0; i < numEntries; i++)
      nmv.push_back(AdjacencyRow());
}

//------------------------------------------------------------------------------------

void Matrix::addNode(int i, int j, double edgeValue, bool symmetric)
{
   saveNode(i, j, edgeValue);

   if (symmetric)
      saveNode(j, i, edgeValue);
}

//------------------------------------------------------------------------------------

void Matrix::compactRows()
{
   try
   {
      for (AdjacencyRows::iterator row = nmv.begin(); row != nmv.end(); ++row)
         if (row->capacity() != row->size())
            AdjacencyRow(*row).swap(*row);
   }
   catch (const std::bad_alloc&)
   {
      throw std::string("Unable to compact adjacency rows.");
   }
   catch (const std::length_error&)
   {
      throw std::string("Adjacency row is too large to compact.");
   }
}

//------------------------------------------------------------------------------------

double Matrix::getNodeMI(int geneId1, int geneId2)
{
   // this function returns a positive MI value if the node has been computed and
   // survived the thresholding; returns 0.0 if the node has been computed but did not
   // survive the thresholding; or returns -1.0 if the node has not yet been computed

   if (nmv.size() > geneId1 && nmv[geneId1].size() > 0)
   {
      AdjacencyRow& row = nmv[geneId1];

      AdjacencyRow::iterator edge = findAdjacencyEdge(row, geneId2);

      return (edge == row.end() ? 0.0 : edge->mutinfo);
   }
   else if (nmv.size() > geneId2 && nmv[geneId2].size() > 0)
   {
      AdjacencyRow& row = nmv[geneId2];

      AdjacencyRow::iterator edge = findAdjacencyEdge(row, geneId1);

      return (edge == row.end() ? 0.0 : edge->mutinfo);
   }
   else
      return -1.0;
}

//------------------------------------------------------------------------------------

static bool protectedByTFLogic(Transfac &transfac,
                               int geneId1, int geneId2, int geneId3)
{
   bool isA = (transfac.find(geneId1) != transfac.end());
   bool isB = (transfac.find(geneId2) != transfac.end());
   bool isC = (transfac.find(geneId3) != transfac.end());

   return ((isA || isB) && !isC);
}

//------------------------------------------------------------------------------------

static int strongerNeighborCount(const std::vector<ArrayValuePair>& miVector,
                                 int candidatePosition, double minMI)
{
   if (candidatePosition == 0 || miVector[0].value <= minMI)
      return 0;

   if (miVector[candidatePosition - 1].value > minMI)
      return candidatePosition;

   int lower = 0;
   int upper = candidatePosition;

   while (lower < upper)
   {
      int middle = lower + (upper - lower) / 2;

      if (miVector[middle].value > minMI)
         lower = middle + 1;
      else
         upper = middle;
   }

   return lower;
}

//------------------------------------------------------------------------------------

static void considerDpiIntermediate(int geneId3, double valueBC,
                                    int strongerCount, double minMI,
                                    int row_idx, int geneId1,
                                    const std::vector<int>& neighborPosition,
                                    Transfac& transfac, int& bestPosition)
{
   if (geneId3 < 0 ||
       static_cast<std::size_t>(geneId3) >= neighborPosition.size())
      return;

   int position = neighborPosition[geneId3];

   if (position < 0 || position >= strongerCount || position >= bestPosition ||
       valueBC <= minMI ||
       (transfac.size() > 0 &&
        protectedByTFLogic(transfac, row_idx, geneId1, geneId3)))
      return;

   bestPosition = position;
}

//------------------------------------------------------------------------------------

static void reduceOneNodeByIntersection(
   Matrix& matrix, int row_idx, double epsilon, Transfac& transfac,
   const DpiNeighborIndex& neighborIndex,
   std::vector<int>& neighborPosition)
{
   AdjacencyRow& row = matrix.nmv[row_idx];

   std::vector<ArrayValuePair> miVector;
   miVector.reserve(row.size());

   for (AdjacencyRow::iterator edge = row.begin(); edge != row.end(); ++edge)
      miVector.push_back(ArrayValuePair(edge->target, edge->mutinfo));

   SortDecreasing_ArrayValuePair sorter;
   std::sort(miVector.begin(), miVector.end(), sorter);

   int numPairs = miVector.size();

   for (int position = 0; position < numPairs; ++position)
      neighborPosition[miVector[position].arrayId] = position;

   for (int i = 0; i < numPairs; i++)
   {
      int    geneId1 = miVector[i].arrayId;
      double valueAB = miVector[i].value;
      double minMI   = valueAB / (1.0 - epsilon);

      int strongerCount = strongerNeighborCount(miVector, i, minMI);
      if (strongerCount == 0)
         continue;

      int bestPosition = strongerCount;
      std::size_t effectiveDegree = neighborIndex.effectiveDegree(geneId1);

      if (static_cast<std::size_t>(strongerCount) <= effectiveDegree)
      {
         // The stronger A-neighbor prefix is the smaller side of the
         // intersection. Preserve its MI-descending order and stop at the first
         // qualifying intermediate, exactly as the legacy loop did.
         for (int j = 0; j < strongerCount; ++j)
         {
            int geneId2 = miVector[j].arrayId;
            double valueBC = matrix.getNodeMI(geneId1, geneId2);

            if (valueBC > minMI && (transfac.size() == 0 ||
                !protectedByTFLogic(transfac, row_idx, geneId1, geneId2)))
            {
               bestPosition = j;
               break;
            }
         }
      }
      else if (neighborIndex.hasDirectRow(geneId1))
      {
         // geneId1 has a stored source row, so getNodeMI() consults that row
         // exclusively. Enumerate only its actual edges and intersect them with
         // the stronger A-neighbor prefix.
         const AdjacencyRow& directRow = neighborIndex.directRow(geneId1);

         for (AdjacencyRow::const_iterator edge = directRow.begin();
              edge != directRow.end(); ++edge)
         {
            considerDpiIntermediate(edge->target, edge->mutinfo,
                                    strongerCount, minMI, row_idx, geneId1,
                                    neighborPosition, transfac, bestPosition);

            if (bestPosition == 0)
               break;
         }
      }
      else
      {
         // A missing/empty geneId1 row makes getNodeMI() fall back to incoming
         // geneId2 -> geneId1 edges. The reverse index reproduces that behavior
         // without probing every stronger A-neighbor.
         const std::vector<ArrayValuePair>& row =
            neighborIndex.incomingRow(geneId1);

         for (std::vector<ArrayValuePair>::const_iterator edge = row.begin();
              edge != row.end(); ++edge)
         {
            considerDpiIntermediate(edge->arrayId, edge->value,
                                    strongerCount, minMI, row_idx, geneId1,
                                    neighborPosition, transfac, bestPosition);

            if (bestPosition == 0)
               break;
         }
      }

      if (bestPosition < strongerCount)
      {
         AdjacencyRow::iterator edge = findAdjacencyEdge(row, geneId1);

         if (edge != row.end())
            edge->intermediate = miVector[bestPosition].arrayId;
      }
   }

   for (int position = 0; position < numPairs; ++position)
      neighborPosition[miVector[position].arrayId] = -1;
}

//------------------------------------------------------------------------------------

void Matrix::reduceOneNode(int row_idx, double epsilon, Transfac& transfac)
{
   compactRows();
   DpiNeighborIndex neighborIndex(*this);
   std::vector<int> neighborPosition(neighborIndex.nodeCount(), -1);
   reduceOneNodeByIntersection(*this, row_idx, epsilon, transfac,
                               neighborIndex, neighborPosition);
}

//------------------------------------------------------------------------------------

void Matrix::reduce(double epsilon, const std::vector<int>& ids, Transfac& transfac)
{
   std::time_t t1, t2;
   std::time(&t1);

   compactRows();
   DpiNeighborIndex neighborIndex(*this);
   std::vector<int> neighborPosition(neighborIndex.nodeCount(), -1);

   if (ids.size() == 0)
   {
      int numMarkers = nmv.size();

      for (int i = 0; i < numMarkers; i++)
         reduceOneNodeByIntersection(*this, i, epsilon, transfac,
                                     neighborIndex, neighborPosition);
   }
   else
   {
      int numIds = ids.size();

      for (int i = 0; i < numIds; i++)
         reduceOneNodeByIntersection(*this, ids[i], epsilon, transfac,
                                     neighborIndex, neighborPosition);
   }

   std::time(&t2);
   std::cout << "DPI running time is: " << std::difftime(t2, t1) << "\n";
}

//------------------------------------------------------------------------------------

int Microarray_Set::Get_Num_Active_Markers() const
{
   int numMarkers = markerset.size();
   int numActive  = 0;

   for (int i = 0; i < numMarkers; i++)
      if (markerset[i].isActive)
         numActive++;

   return numActive;
}

//------------------------------------------------------------------------------------

bool Microarray_Set::isSameGene(int i, int j) const
{
   return (i == j || markerset[i].label == markerset[j].label &&
                     markerset[i].label != "---");
}

//------------------------------------------------------------------------------------

int Microarray_Set::getAccessionId(const std::string& accnum) const
{
   if (!accessionIdsCurrent)
      rebuildAccessionIndex();

   std::unordered_map<std::string, int>::const_iterator found =
      accessionIds.find(accnum);

   if (found != accessionIds.end())
      return found->second;

   return -1; // did not find it
}

//------------------------------------------------------------------------------------

void Microarray_Set::rebuildAccessionIndex() const
{
   accessionIds.clear();
   accessionIds.reserve(markerset.size());

   int numMarkers = markerset.size();

   for (int i = 0; i < numMarkers; i++)
      // emplace() does not replace an existing entry. This preserves the legacy
      // linear lookup behavior when an expression matrix contains duplicate
      // accession IDs: the first expression row remains authoritative.
      accessionIds.emplace(markerset[i].accnum, i);

   accessionIdsCurrent = true;
}

//------------------------------------------------------------------------------------

int Microarray_Set::getProbeId(const std::string& label) const
{
   int len = label.length();

   for (int i = 0; i < len; i++)
      if (!std::isdigit(label[i]))     // label is not a number
         return getAccessionId(label); // perhaps it is an accession ID

   return std::atoi(label.c_str());    // label is a number, convert it to int
}

//------------------------------------------------------------------------------------

void Microarray_Set::Set_ColHeader(int i, const std::string& hdr)
{
   while (i >= header.size())
      header.push_back("");

   header[i] = hdr;
}

//------------------------------------------------------------------------------------

void Microarray_Set::Set_Marker(int i, const Marker& m)
{
   while (i >= markerset.size())
      markerset.push_back(Marker());

   markerset[i] = m;
   accessionIdsCurrent = false;
}

//------------------------------------------------------------------------------------

void Microarray_Set::Set_Probe(int i, int j, Probe p)
{
   while (i >= uarrays.size())
      uarrays.push_back(Microarray());

   while (j >= uarrays[i].size())
      uarrays[i].push_back(Probe());

   uarrays[i][j] = p;
}

//------------------------------------------------------------------------------------
// Read a matrix data file into a Microarray_Set. The file consists of sets of lines
// that start with a description of a marker (including accession number and
// descriptive label). The remainder of the lines are probe values (the values only,
// or (value, pvalue) pairs) for each microarray for the marker that started the line.
// We assign each marker an ID number according to the line number it appeared on.
// Conceptually, the file can be thought of as a matrix of probe data in column major
// order.
//
// The program can automatically tell whether the input file consists of
// (value, pvalue) pairs or single expression values. It will read the data and set
// up the structure correctly according to the format of the file. (If the data
// contains only expression measurements, a p-value of 0.0 will be assigned to each
// expression value.)
//------------------------------------------------------------------------------------

int Microarray_Set::readMarkerWithPvalue(std::istream& in, const int arrno)
{
   // Read a marker line with (value, pvalue) format

   int markern = 0; // number of expression values of each probe

   in.exceptions(std::ios_base::badbit | std::ios_base::failbit);

   try
   {
      do
      {
         std::string sPVal;
         double val, pval;

         in >> val >> sPVal;

         switch (sPVal[0])
         {
            case 'A': pval = 0.7; break;
            case 'M': pval = 0.5; break;
            case 'P': pval = 0.1; break;
            default : pval = std::atof(sPVal.c_str());
         }

         if (!std::isfinite(val))
         {
            std::ostringstream s;
            s << "Non-finite expression value at gene row " << arrno + 1
              << ", observation " << markern + 1 << ".";
            throw s.str();
         }

         if (!std::isfinite(pval))
         {
            std::ostringstream s;
            s << "Non-finite expression p-value at gene row " << arrno + 1
              << ", observation " << markern + 1 << ".";
            throw s.str();
         }

         Set_Probe(markern++, arrno, Probe(val, pval));
      }
      while (in.good() && in.get() != '\015' && in.peek() != EOF &&
             in.peek() != '\015');
   }
   catch (std::ios_base::failure& f)
   {
      std::ostringstream s;
      s << "Could not read data at line no: " << arrno + 2;
      std::cout << s.str() << std::endl;
      throw s.str();
   }

   return markern;
}

//------------------------------------------------------------------------------------

int Microarray_Set::readMarkerNoPvalue(std::istream& in, const int arrno)
{
   // Read a marker line with no pvalues

   int markern = 0; // number of expression values of each probe

   in.exceptions(std::ios_base::badbit | std::ios_base::failbit);

   try
   {
      do
      {
         double val;

         in >> val;

         if (!std::isfinite(val))
         {
            std::ostringstream s;
            s << "Non-finite expression value at gene row " << arrno + 1
              << ", observation " << markern + 1 << ".";
            throw s.str();
         }

         Set_Probe(markern++, arrno, Probe(val, 0.0));
      }
      while (in.good() && in.get() != '\015' && in.peek() != EOF &&
             in.peek() != '\015');
   }
   catch (std::ios_base::failure& f)
   {
      std::ostringstream s;
      s << "Could not read data at line no: " << arrno;
      std::cout << s.str() << std::endl;
      throw s.str();
   }

   return markern;
}

//------------------------------------------------------------------------------------

int Microarray_Set::readHeader(std::istream& in)
{
   in.exceptions(std::ios_base::badbit | std::ios_base::failbit);

   try
   {
      do
      {
         std::string hdr;

         std::getline(in, hdr, '\t');

         Set_ColHeader(header.size(), hdr);
      }
      while (in.good() && in.peek() != '\015' && in.peek() != EOF);
   }
   catch (std::ios_base::failure& f)
   {
      throw std::string("Error while reading file headers"
                        "(win/*nux/mac end of line?).");
   }

   return header.size() - 2;
}

//------------------------------------------------------------------------------------

static void validateObservationCount(int expected, int actual, int lineNumber)
{
   if (actual != expected)
   {
      std::ostringstream s;
      s << "Incorrect expression dimensions at line " << lineNumber
        << ": expected " << expected << " observations, found " << actual << ".";
      throw s.str();
   }
}

//------------------------------------------------------------------------------------

void Microarray_Set::read(const std::string& filename)
{
   std::ifstream in(filename.c_str());
   if (!in.is_open())
      throw "Unable to open " + filename;

   read(in);

   in.close();
}

//------------------------------------------------------------------------------------

void Microarray_Set::read(std::istream& in)
{
   std::ios_base::iostate oldState = in.exceptions();
   in.exceptions(std::ios_base::badbit | std::ios_base::failbit);

   try
   {
      std::string line;
      std::getline(in, line);
      std::istringstream sin(line);
      int arrno = readHeader(sin);

      if (arrno < 2)
      {
         std::ostringstream s;
         s << "Expression matrix must contain at least 2 observation columns; found "
           << arrno << ".";
         throw s.str();
      }

      // we need to decide whether the input file contain only expression
      // or (value, pvalue) pairs

      int bypass_line_cnt = 0;
      std::getline(in, line);

      while (line.length() >= 11 && line.substr(0, 11) == "Description")
      {
         bypass_line_cnt++;
         std::getline(in, line);
      }

      std::cout << "\n[READ] " << bypass_line_cnt
                << " Description lines bypassed." << std::endl;

      const int firstDataLine = bypass_line_cnt + 2;

      std::istringstream pin(line);
      std::vector<std::string> firstprobe;
      std::string token;

      do
      {
         std::getline(pin, token, '\t');
         firstprobe.push_back(token);
      }
      while (pin.good() && pin.peek() != '\015' && pin.peek() != EOF);

      int valueNo = firstprobe.size() - 2;

      if (arrno == valueNo)
      {
         // single expression values
         std::cout << "[READ] P-value columns not found." << std::endl;

         int proben = 0; // probe number

         std::istringstream fin(line);
         std::string accnum, label;

         std::getline(fin, accnum, '\t');
         std::getline(fin, label,  '\t');

         accnum = "_" + accnum;

         Marker m(proben, accnum, label);
         Set_Marker(proben, m);
         int observations = readMarkerNoPvalue(fin, proben);
         validateObservationCount(arrno, observations, firstDataLine);
         proben++;

         while (in.good() && in.peek() != EOF)
         {
            std::getline(in, line);
            std::istringstream sin(line);

            std::getline(sin, accnum, '\t');
            std::getline(sin, label,  '\t');

            accnum = "_" + accnum;

            Marker m(proben, accnum, label);
            Set_Marker(proben, m);

            int observations = readMarkerNoPvalue(sin, proben);
            validateObservationCount(arrno, observations,
                                     firstDataLine + proben);

            proben++;
         }
      }
      else if (2 * arrno == valueNo)
      {
         // (value, pvalue) pairs
         std::cout << "[READ] (value, p-value) pairs found." << std::endl;

         int proben = 0; // probe number

         std::istringstream fin(line);
         std::string accnum, label;

         std::getline(fin, accnum, '\t');
         std::getline(fin, label,  '\t');

         accnum = "_" + accnum;

         Marker m(proben, accnum, label);
         Set_Marker(proben, m);
         int observations = readMarkerWithPvalue(fin, proben);
         validateObservationCount(arrno, observations, firstDataLine);
         proben++;

         while (in.good() && in.peek() != EOF)
         {
            std::getline(in, line);
            std::istringstream sin(line);

            std::getline(sin, accnum, '\t');
            std::getline(sin, label,  '\t');

            accnum = "_" + accnum;

            Marker m(proben, accnum, label);
            Set_Marker(proben, m);

            int observations = readMarkerWithPvalue(sin, proben);
            validateObservationCount(arrno, observations,
                                     firstDataLine + proben);

            proben++;
         }
      }
      else
      {
         std::ostringstream s;
         s << "Incorrect expression dimensions at line " << firstDataLine
           << ": header declares " << arrno << " observations, but the first row "
           << "contains " << valueNo << " data fields; expected " << arrno
           << " expression values or " << 2 * arrno
           << " expression/p-value fields.";
         throw s.str();
      }
   }
   catch (std::ios_base::failure& f)
   {
      throw std::string("Could not read data.(Last line empty?)");
   }

   in.exceptions(oldState); //reset excep state
}

//------------------------------------------------------------------------------------

int Microarray_Set::filter(const std::vector<int>& ids, double minMean,
                           double minSigma, int ctlid)
{
   int numIds         = ids.size();
   int numMicroarrays = uarrays.size();
   int numMarkers     = markerset.size();

   for (int i = 0; i < numMarkers; i++)
      markerset[i].isActive = false;

   int numDisabled = numMarkers;

   for (int i = 0; i < numMarkers; i++)
      if (i != ctlid)
      {
         double nx = 0.0, nxx = 0.0;

         for (int j = 0; j < numIds; j++)
         {
            double v = uarrays[ids[j]][i].value;
            nx  += v;
            nxx += v * v;
         }

         double mean     = nx / numMicroarrays;
         double variance = (numMicroarrays * nxx - nx * nx) /
                           (numMicroarrays * numMicroarrays);
         double stdev    = std::sqrt(variance);

         if (mean >= minMean && stdev >= mean * minSigma)
         {
            markerset[i].isActive = true;
            numDisabled--;
         }
      }

   return numDisabled;
}

//------------------------------------------------------------------------------------

int Microarray_Set::filter(double minMean, double minSigma, int ctlid)
{
   int numMicroarrays = uarrays.size();

   std::vector<int> v;

   for (int id = 0; id < numMicroarrays; id++)
      v.push_back(id);

   return filter(v, minMean, minSigma, ctlid);
}

//------------------------------------------------------------------------------------

void Microarray_Set::computeMarkerVariance(const std::vector<int> *arrays)
{
   int numMarkers = markerset.size();

   for (int i = 0; i < numMarkers; i++)
      markerset[i].var = variance(i, arrays);
}

//------------------------------------------------------------------------------------

void Microarray_Set::computeMarkerBandwidth(const std::vector<int> *arrays)
{
   computeMarkerVariance(arrays);

   int n = (arrays ? arrays->size() : uarrays.size());

   double *data = new double[n];

   int numMarkers = markerset.size();

   for (int i = 0; i < numMarkers; i++)
   {
      double prop = 1.06; // Gaussian
      int dim = 1;        // dimension of data

      double stdev = std::sqrt(markerset[i].var);

      for (int j = 0; j < n; j++)
         data[j] = uarrays[arrays ? arrays->at(j) : j][i].value;

      std::sort(data, data + n);

      double iqr = interQuartileRange(data, n);

      double iqrSig = 0.7413 * iqr; // estimate of sigma
      if (iqrSig == 0.0)
         iqrSig = stdev;

      double sig = std::min(stdev, iqrSig);

      markerset[i].bandwidth = prop * sig * std::pow(n, -1.0 / (dim + 4));
   }

   delete[] data;
}

//------------------------------------------------------------------------------------

double Microarray_Set::variance(int m, const std::vector<int> *arrays)
{
   // compute the variance of a gene expression vector

   int n = (arrays ? arrays->size() : uarrays.size());

   double s = 0.0, ss = 0.0;

   for (int i = 0; i < n; i++)
   {
      double v = uarrays[arrays ? arrays->at(i) : i][m].value;
      s  += v;
      ss += v * v;
   }

   return (ss - s * s / n) / (n - 1);
}

//------------------------------------------------------------------------------------

void Microarray_Set::getHighLowPercent(double x, int mId, std::vector<int>& lower,
                                       std::vector<int>& upper)
{
   std::vector<ArrayValuePair> sortArray;
   SortIncreasing_ArrayValuePair sorter;

   int numMicroarrays = uarrays.size();

   for (int id = 0; id < numMicroarrays; id++)
      sortArray.push_back(ArrayValuePair(id, uarrays[id][mId].value));

   std::sort(sortArray.begin(), sortArray.end(), sorter);

   int idPercNo = numMicroarrays * x;

   for (int id = 0; id < idPercNo; id++)
   {
      lower.push_back(sortArray[id].arrayId);
      upper.push_back(sortArray[numMicroarrays - idPercNo - 1 + id].arrayId);
   }
}

//------------------------------------------------------------------------------------

void Microarray_Set::bootStrap(std::vector<int>& boot, const std::vector<int> *arrays)
{
   boot.clear();

   int numIds = (arrays ? arrays->size() : uarrays.size());

   for (int id = 0; id < numIds; id++)
   {
      int r = std::rand() % numIds;
      boot.push_back(arrays ? arrays->at(r) : r);
   }
}

//------------------------------------------------------------------------------------

static std::uint32_t randomBelow(std::mt19937& generator, std::uint32_t bound)
{
   // Rejection sampling avoids the modulo bias of rand() % bound while retaining
   // a fully specified mt19937 bit stream across supported C++ runtimes.
   const std::uint32_t threshold = static_cast<std::uint32_t>(-bound) % bound;

   while (true)
   {
      const std::uint32_t value = static_cast<std::uint32_t>(generator());
      if (value >= threshold)
         return value % bound;
   }
}

//------------------------------------------------------------------------------------

void Microarray_Set::sampleWithoutReplacement(std::vector<int>& sample,
                                              int sampleSize,
                                              unsigned int seed,
                                              const std::vector<int> *arrays) const
{
   const int populationSize = (arrays ? arrays->size() : uarrays.size());

   if (sampleSize < 0 || sampleSize > populationSize)
      throw std::string("Without-replacement sample size is outside the eligible "
                        "observation population.");

   std::vector<int> positions(populationSize);
   std::iota(positions.begin(), positions.end(), 0);

   std::mt19937 generator(seed);

   // A partial Fisher-Yates shuffle chooses a uniform fixed-size subset in O(N)
   // memory and O(m) random draws.  Sorting the chosen positions restores their
   // order in the eligible population; it does not change the selected subset.
   for (int position = 0; position < sampleSize; position++)
   {
      const std::uint32_t remaining =
         static_cast<std::uint32_t>(populationSize - position);
      const int selected =
         position + static_cast<int>(randomBelow(generator, remaining));
      std::swap(positions[position], positions[selected]);
   }

   positions.resize(sampleSize);
   std::sort(positions.begin(), positions.end());

   sample.clear();
   sample.reserve(sampleSize);

   for (std::vector<int>::const_iterator position = positions.begin();
        position != positions.end(); ++position)
      sample.push_back(arrays ? arrays->at(*position) : *position);
}

//------------------------------------------------------------------------------------

void Microarray_Set::addNoise()
{
   int numMicroarrays = uarrays.size();
   int numMarkers     = markerset.size();

   for (int id = 0; id < numMicroarrays; id++)
      for (int mid = 0; mid < numMarkers; mid++)
      {
         double r = std::rand();
         double noise = (r / RAND_MAX) * 1e-10;

         uarrays[id][mid].value += noise;
      }
}

//------------------------------------------------------------------------------------

static void BuildRankCache(const Microarray_Set& data, int maNum,
                           const std::vector<int> *arrays,
                           const std::vector<bool>& markersNeeded,
                           std::vector<int>& rankRows,
                           std::vector<int>& rankCache)
{
   // Ranks belong to the exact selected/bootstrap observation sequence.  Position,
   // not original array ID, is the secondary sort key so repeated bootstrap
   // observations retain the same tie behavior as the legacy per-edge sorts.
   if (arrays != NULL && static_cast<int>(arrays->size()) != maNum)
      throw std::string("Rank-cache observation count does not match array selection.");

   const std::size_t numMarkers = data.markerset.size();

   if (markersNeeded.size() != numMarkers)
      throw std::string("Rank-cache marker selection has an invalid size.");

   const std::size_t numRankedMarkers =
      std::count(markersNeeded.begin(), markersNeeded.end(), true);

   if (maNum <= 0 || numRankedMarkers > rankCache.max_size() / maNum)
      throw std::string("Rank-cache dimensions are invalid or too large.");

   rankRows.assign(numMarkers, -1);
   rankCache.resize(numRankedMarkers * static_cast<std::size_t>(maNum));

   if (numRankedMarkers == 0)
      return;

   std::vector<RankedObservation> observations(maNum);
   Sort_RankedObservation sorter;
   int rankRow = 0;

   for (std::size_t marker = 0; marker < numMarkers; marker++)
   {
      if (!markersNeeded[marker])
         continue;

      rankRows[marker] = rankRow++;

      for (int position = 0; position < maNum; position++)
      {
         int arrayId = (arrays ? arrays->at(position) : position);
         observations[position] =
            RankedObservation(data.uarrays[arrayId][marker].value, position);
      }

      std::sort(observations.begin(), observations.end(), sorter);

      std::size_t offset = static_cast<std::size_t>(rankRows[marker]) * maNum;

      for (int rank = 0; rank < maNum; rank++)
         rankCache[offset + observations[rank].position] = rank + 1;
   }
}

//------------------------------------------------------------------------------------

double Microarray_Set::calculateMI(int maNum, int probeId1, int probeId2,
                                   double threshold, double noise2, int nparLimit,
                                   const std::vector<int>& rankCache,
                                   const std::vector<int>& rankRows,
                                   AdaptivePartitionWorkspace& workspace) const
{
   // compute mutual information between two gene expression vectors;
   // zero is returned if there is no connection

   if (isSameGene(probeId1, probeId2))
      return 0.0;

   if (rankRows[probeId1] < 0 || rankRows[probeId2] < 0)
      throw std::string("Internal error: MI requested for an uncached marker rank.");

   std::size_t offset1 = static_cast<std::size_t>(rankRows[probeId1]) * maNum;
   std::size_t offset2 = static_cast<std::size_t>(rankRows[probeId2]) * maNum;

   int pairMaxNpar = 1;
   double mi = computeAdaptivePartitionMI(&rankCache[offset1], &rankCache[offset2],
                                          maNum, nparLimit, workspace,
                                          &pairMaxNpar);

   if (pairMaxNpar > maxNpar)
      maxNpar = pairMaxNpar;

   if (!std::isfinite(mi))
   {
      std::ostringstream s;
      s << "Non-finite MI calculated for edge "
        << markerset[probeId1].accnum.substr(1) << " -> "
        << markerset[probeId2].accnum.substr(1) << " using " << maNum
        << " observations.";
      throw s.str();
   }

   if (mi < threshold)
      return 0.0;

   if (noise2 <= 0.0)
      return mi;

   double v1 = markerset[probeId1].var;
   double v2 = markerset[probeId2].var;

   double lambda = (v1 / (v1 - noise2)) * (v2 / (v2 - noise2));

   double value  = 1 + (std::exp(2 * mi) - 1) * (1 - 1 / lambda);

   double correctedMi = mi + 0.5 * std::log(value);

   if (!std::isfinite(correctedMi))
   {
      std::ostringstream s;
      s << "Non-finite noise-corrected MI calculated for edge "
        << markerset[probeId1].accnum.substr(1) << " -> "
        << markerset[probeId2].accnum.substr(1) << " using " << maNum
        << " observations.";
      throw s.str();
   }

   return correctedMi;
}

//------------------------------------------------------------------------------------

void Microarray_Set::computeOneRow(int maNum, Matrix& matrix, double threshold,
                                   int row_idx, int numMarkers, int controlId,
                                   bool half_matrix, bool symmetric, double noise2,
                                   int nparLimit,
                                   const std::vector<int>& rankCache,
                                   const std::vector<int>& rankRows,
                                   AdaptivePartitionWorkspace& workspace) const
{
   // this function computes one row of the adjacency matrix; it is called by
   // createEdgeMatrix(); note that since the adjacency matrix is symmetric,
   // only the upper right triangle is computed; so here this function computes only
   // the "upper right triangle" of the row

   for (int j = (half_matrix ? row_idx + 1 : 0); j < numMarkers; j++)
      if (j != controlId && markerset[j].isActive)
      {
         double edge = calculateMI(maNum, row_idx, j, threshold, noise2, nparLimit,
                                   rankCache, rankRows, workspace);
         if (edge != 0.0)
            matrix.addNode(row_idx, j, edge, symmetric);
      }
}

//------------------------------------------------------------------------------------

void Microarray_Set::createEdgeMatrix(int maNum, Matrix& matrix, double threshold,
                                      int controlId, double noise2, int nparLimit,
                                      const std::vector<int>& ids,
                                      const std::vector<int> *arrays) const
{
   // if controlId == -1, there is no constraint (use all arrays to compute the
   // mutual information; if ids.size == 0, all genes will be computed; otherwise,
   // only selected genes will be computed; arrays points to a vector of array ids
   // used for mutual information computation

   std::time_t t1, t2;
   std::time(&t1); // get a timestamp to assess how much time this will take

   int numMarkers = markerset.size();
   int count      = (ids.size() == 0 ? numMarkers : ids.size());
   int step       = std::ceil(0.1 * count);

   matrix.adjacencyRowsPresent.assign(numMarkers, false);

   std::vector<bool> markersNeeded(numMarkers, false);
   bool hasSource = false;

   if (ids.size() == 0)
   {
      for (int marker = 0; marker < numMarkers; marker++)
         if (marker != controlId && markerset[marker].isActive)
            hasSource = true;
   }
   else
      for (std::vector<int>::const_iterator id = ids.begin(); id != ids.end(); ++id)
         if (*id != controlId)
         {
            markersNeeded[*id] = true; // requested hubs may be inactive
            hasSource = true;
         }

   if (hasSource)
      for (int marker = 0; marker < numMarkers; marker++)
         if (marker != controlId && markerset[marker].isActive)
            markersNeeded[marker] = true;

   std::size_t numRankedMarkers =
      std::count(markersNeeded.begin(), markersNeeded.end(), true);
   double rankCacheMiB = numRankedMarkers * static_cast<double>(maNum) *
                         sizeof(int) / (1024.0 * 1024.0);

   std::vector<int> rankRows;
   std::vector<int> rankCache;
   AdaptivePartitionWorkspace workspace;

   try
   {
      BuildRankCache(*this, maNum, arrays, markersNeeded, rankRows, rankCache);
   }
   catch (const std::bad_alloc&)
   {
      std::ostringstream s;
      s << "Unable to allocate rank cache for " << numRankedMarkers << " genes x "
        << maNum << " observations (" << rankCacheMiB << " MiB).";
      throw s.str();
   }
   catch (const std::length_error&)
   {
      std::ostringstream s;
      s << "Rank cache is too large for " << numRankedMarkers << " genes x "
        << maNum << " observations (" << rankCacheMiB << " MiB).";
      throw s.str();
   }

   std::time(&t2);
   std::cout << "[MI] Cached ranks for " << numRankedMarkers << " genes x " << maNum
             << " observations in " << std::difftime(t2, t1) << " seconds ("
             << rankCacheMiB << " MiB)." << std::endl;

   if (hasSource)
      try
      {
         workspace.initialize(maNum, nparLimit);
      }
      catch (const std::bad_alloc&)
      {
         std::ostringstream s;
         s << "Unable to allocate adaptive-partitioning workspace for " << maNum
           << " observations and partition limit " << nparLimit << ".";
         throw s.str();
      }
      catch (const std::length_error&)
      {
         std::ostringstream s;
         s << "Adaptive-partitioning workspace is too large for " << maNum
           << " observations and partition limit " << nparLimit << ".";
         throw s.str();
      }

   if (ids.size() == 0) // all genes will be computed
   {
      matrix.createEntries(count);

      for (int i = 0; i < count; i++)
      {
         if (i != controlId && markerset[i].isActive)
          {
             matrix.markAdjacencyRow(i);
             computeOneRow(maNum, matrix, threshold, i, numMarkers, controlId,
                           true, true, noise2, nparLimit, rankCache, rankRows,
                           workspace);
          }

         if ((i + 1) % step == 0)
         {
            std::time(&t2);
            std::cout << 10 * (i + 1) / step << "%, time: "
                      << std::difftime(t2, t1) << std::endl;
         }
      }
   }
   else // only selected genes will be computed
      for (int i = 0; i < count; i++)
         if (ids[i] != controlId)
         {
            while (matrix.nmv.size() <= ids[i])
               matrix.nmv.push_back(AdjacencyRow());

            matrix.markAdjacencyRow(ids[i]);
            computeOneRow(maNum, matrix, threshold, ids[i], numMarkers, controlId,
                          false, false, noise2, nparLimit, rankCache, rankRows,
                          workspace);

            if ((i + 1) % step == 0)
            {
               std::time(&t2);
               std::cout << 10 * (i + 1) / step << "%, time: "
                         << std::difftime(t2, t1) << std::endl;
            }
         }

   matrix.compactRows();

   std::time(&t2);
   std::cout << "Gene: " << count << " Time: " << std::difftime(t2, t1) << std::endl;
}
