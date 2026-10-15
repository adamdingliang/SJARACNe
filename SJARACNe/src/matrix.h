//------------------------------------------------------------------------------------
// Copyright (C) 2003  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017
//------------------------------------------------------------------------------------

#ifndef MATRIX_H__
#define MATRIX_H__

#include <cstdint>
#include <iostream>
#include <map>
#include <unordered_map>
#include "param.h"

typedef std::map<int, int> Transfac;

class Microarray_Set;
class AdaptivePartitionWorkspace;

//------------------------------------------------------------------------------------
// The Matrix class represents the adjacency matrix data structure that consists a
// list of genes and their connections to other genes.  One can think of an adjacency
// matrix as a two dimensional array: each entry in the array is called a node;
// a node[i,j] will present if there is a connection between gene i and gene j, absent
// otherwise. The matrix is symmetric.
//------------------------------------------------------------------------------------

class AdjacencyEdge // represents an edge between two genes
{
public:
   AdjacencyEdge(int inTarget=-1, double mi=0.0)
      : mutinfo(mi), target(inTarget), intermediate(-1) { }

   double mutinfo;   // Get_MutualInfo(), Set_MutualInfo()
   int target;        // target gene ID
   int intermediate; // Get_Intermediate(), Set_Intermediate()
};

typedef std::vector<AdjacencyEdge> AdjacencyRow;
typedef std::vector<AdjacencyRow> AdjacencyRows;

//------------------------------------------------------------------------------------

struct DpiEdgeStatistics
{
   DpiEdgeStatistics()
      : preEdges(0), prunedEdges(0), postEdges(0) { }

   std::uint64_t preEdges;
   std::uint64_t prunedEdges;
   std::uint64_t postEdges;
};

//------------------------------------------------------------------------------------

class Matrix
{
public:
   Matrix()
      : nmv(), adjacencyRowsPresent(), writeTriangular(false),
        writeReduced(false), writeEmptyGenes(false) { }

   AdjacencyRows nmv;
   std::vector<bool> adjacencyRowsPresent;

   bool writeTriangular; // if true, only triangular half of matrix will be written
   bool writeReduced;    // if true, intermediate nodes will be written
   bool writeEmptyGenes; // if true, empty lines will be written

   void saveNode(int i, int j, double mi);

   void read(Microarray_Set& data, const Parameter& p);
   void read(std::istream& in, Microarray_Set& data, const Parameter& p);
   void markAdjacencyRow(int geneId);
   bool hasAdjacencyRow(int geneId) const;
   bool hasEnoughSourceRowsForDpi() const;

   void writeGeneLine(std::ostream& out, const Microarray_Set& data, int geneId);
   void writeGeneList(const Microarray_Set& data, const std::string& name,
                      int probeId);

   void write(std::ostream& out, const Microarray_Set& data,
              const std::vector<int>& ids);
   void write(const Microarray_Set& data, const std::vector<int>& ids,
              const Parameter& p, bool writeFull=true);

   void createEntries(int numEntries);
   void addNode(int i, int j, double edgeValue, bool symmetric);
   void compactRows();

   double getNodeMI(int geneId1, int geneId2);
   void reduceOneNode(int row_idx, double epsilon, Transfac& transfac);
   void reduce(double epsilon, const std::vector<int>& ids, Transfac& transfac);
   DpiEdgeStatistics dpiEdgeStatistics(const std::vector<int>& ids) const;
   void writeDpiWitnessDiagnostics(const Parameter& p,
                                   const std::vector<int>& ids,
                                   const Transfac& transfac) const;
};

//------------------------------------------------------------------------------------
// A Microarray dataset consists of multiple arrays, each of which measures the
// expression of thousands of Markers. Here we represent each of these measurements
// using a Probe.
//
// Microarrays are really silicon chips that contain probes to detect genes.  Eachi
// probe is a tuple containing values and pvalues, and there are n such tuples
// arranged in an array, indexed by the same numbers which index Markers in a set.
// Since Markers describe genes, the correspondence is natural.  Associated with each
// set of Microarrays is also a Marker_Set object.  If the cardinality of the marker
// set is m, then each Microarray will have m probes associated with it.
//------------------------------------------------------------------------------------

class Marker // describes a gene in a detailed way
{
public:
   Marker(int inId=0, const std::string& inAccnum="", const std::string& inLabel="",
          double inVar=0.0, double inBandwidth=0.0, bool inIsActive=true,
          bool inIsControl=false)
      : idnum(inId), accnum(inAccnum), label(inLabel), var(inVar),
        bandwidth(inBandwidth), isActive(inIsActive), isControl(inIsControl) { }

   int idnum;          // Get_ID(), Set_ID()
   std::string accnum; // Get_Accnum(), Set_Accnum(), hasAccession()
   std::string label;  // Get_Label(), Set_Label()
   double var;         // Get_Var(), Set_Var()
   double bandwidth;   // Get_Bandwidth(), Set_Bandwidth()
   bool isActive;      // Enable(), Disable(), Enabled()
   bool isControl;     // set_Control(), set_No_Control(), Is_Control()

   friend std::ostream& operator<<(std::ostream& out, const Marker& m);
};

typedef std::vector<Marker> Marker_Set; // array of markers, indexed by the ID
                                        // numbers contained by the markers
//------------------------------------------------------------------------------------

class Probe
{
public:
   Probe(double inValue=0.0, double inPValue=0.0)
      : value(inValue), pvalue(inPValue) { }

   double value;  // Get_Value(),  Set_Value()
   double pvalue; // Get_PValue(), Set_PValue()

   friend std::ostream& operator<<(std::ostream& out, const Probe& p);
};

typedef std::vector<Probe> Microarray;
typedef std::vector<Microarray> Microarray_Vector;

//------------------------------------------------------------------------------------

class Microarray_Set
{
public:
   Microarray_Set()
      : markerset(), uarrays(), header(), accessionIds(),
        accessionIdsCurrent(false) { }

   Marker_Set markerset;                 // Get_Num_Markers(), Get_Marker(),
                                         // Get_Marker_AffyId() gets marker's accnum,
                                         // getMarkerVariance() gets marker's variance
   Microarray_Vector uarrays;            // Get_Num_Microarrays(), Get_Microarray(),
                                         // Get_Probe(), Get_Value(), Get_PValue()
   std::vector<std::string> header;      // Get_Header(), Get_Array_Header()

   friend std::ostream& operator<<(std::ostream& out, const Microarray_Set& ms);

   int  Get_Num_Active_Markers() const;

   bool isSameGene(int i, int j) const;

   int  getAccessionId(const std::string& accnum) const; // formerly, get_Id()
   int  getProbeId(const std::string& label) const;

   void Set_ColHeader(int i, const std::string& hdr);
   void Set_Marker(int i, const Marker& m);
   void Set_Probe (int i, int j, Probe p);

   int readMarkerWithPvalue(std::istream& in, const int arrno);
   int readMarkerNoPvalue  (std::istream& in, const int arrno);
   int readHeader(std::istream& in);
   void read(const std::string& filename);
   void read(std::istream& in);

   int  filter(const std::vector<int>& ids, double minMean=50.0, double minSigma=20.0,
               int ctlid=-1);
   int  filter(double minMean=50.0, double minSigma=20.0, int ctlid=-1);

   void computeMarkerVariance (const std::vector<int> *arrays);
   void computeMarkerBandwidth(const std::vector<int> *arrays);

   double variance(int m, const std::vector<int> *arrays);

   void getHighLowPercent(double x, int mId, std::vector<int>& lower,
                          std::vector<int>& upper);

   void bootStrap(std::vector<int>& boot, const std::vector<int> *arrays);
   void sampleWithoutReplacement(std::vector<int>& sample, int sampleSize,
                                 unsigned int seed,
                                 const std::vector<int> *arrays) const;
   void addNoise();

private:
   void rebuildAccessionIndex() const;

   mutable std::unordered_map<std::string, int> accessionIds;
   mutable bool accessionIdsCurrent;

   double calculateMI(int maNum, int probeId1, int probeId2, double threshold,
                      double noise2, int nparLimit,
                      const std::vector<int>& rankCache,
                      const std::vector<int>& rankRows,
                      AdaptivePartitionWorkspace& workspace) const;
   void computeOneRow(int maNum, Matrix& matrix, double threshold, int row_idx,
                      int numMarkers, int controlId, bool half_matrix,
                      bool symmetric, double noise2, int nparLimit,
                      const std::vector<int>& rankCache,
                      const std::vector<int>& rankRows,
                      AdaptivePartitionWorkspace& workspace) const;

public:
   void createEdgeMatrix(int maNum, Matrix& matrix, double threshold, int controlId,
                         double noise2, int nparLimit, const std::vector<int>& ids,
                         const std::vector<int> *arrays) const;
};

#endif
