#ifndef SJARACNE_CONSENSUS_H
#define SJARACNE_CONSENSUS_H

#include <string>
#include <vector>
#include "param.h"

class Microarray_Set;
class Matrix;

/**
 * Build the directed, K-supported mean-MI graph for the caller's final DPI pass.
 *
 * Each manifest entry must represent one distinct resample. Every raw graph
 * must contain exactly the supplied source rows, in the supplied order,
 * including source-only rows for empty regulons. A graph contributes at most
 * one observation per source-target pair. The mean divides by supporting runs,
 * not all runs; an emitted zero MI still counts as support.
 *
 * Numeric accumulators use O(G) storage and B input streams remain open. The
 * retained graph, expression data, row/stream buffers and path/ID strings also
 * contribute to memory use; this is not a total-memory O(G) guarantee.
 *
 * @param data Expression markers/IDs; matched source markers become active.
 * @param sources Distinct valid marker indices defining the raw row order.
 * @param matrix Receives the selected edges and explicit source-row presence.
 * @param parameters Supplies the manifest/K and receives run/edge counters.
 * @return Resolved raw input paths in manifest order, for later output checks.
 * @throws std::string If an input, row, MI value or protected output is invalid.
 *
 * No DPI or output writing occurs here. On failure, matrix/data/counters may
 * already be partially updated; the caller must abandon that reconstruction.
 */
std::vector<std::string> buildConsensusFromManifest(
   Microarray_Set& data, const std::vector<int>& sources, Matrix& matrix,
   Parameter& parameters);

/**
 * Reject existing output aliases to the supplied network inputs.
 *
 * Device/inode identity covers aliases and hardlinks to the expression file,
 * source/annotation lists, manifest and supplied raw graphs. This read-only
 * preflight does not create/truncate outputs or lock filesystem identities.
 * The caller repeats it after resolving default output names, before writing.
 *
 * @param parameters Configured inputs, network output and optional DPI witness.
 * @param rawInputs Raw graphs to protect; empty for raw-MI inference mode.
 * @throws std::string If an existing output aliases a protected input or the
 *     two outputs refer to the same existing file.
 */
void validateNativeOutputPaths(const Parameter& parameters,
                               const std::vector<std::string>& rawInputs);

#endif
