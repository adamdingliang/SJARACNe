#ifndef SJARACNE_CONSENSUS_H
#define SJARACNE_CONSENSUS_H

#include <string>
#include <vector>
#include "param.h"

class Microarray_Set;
class Matrix;

// Read canonical, unpruned source rows in manifest order. Additional working
// memory is O(genes + runs); only the K-selected directed graph is retained.
std::vector<std::string> buildConsensusFromManifest(
   Microarray_Set& data, const std::vector<int>& sources, Matrix& matrix,
   Parameter& parameters);

// Protect every scientific input, including aliases/hardlinks, from either
// native output. Called before aggregation and again once default output names
// have been resolved. It does not create or truncate an output file.
void validateNativeOutputPaths(const Parameter& parameters,
                               const std::vector<std::string>& rawInputs);

#endif
