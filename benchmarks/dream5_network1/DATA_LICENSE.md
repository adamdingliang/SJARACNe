# DREAM5 Network 1 data provenance

Source: [DREAM5 Gene Network Inference challenge archive on Zenodo](https://zenodo.org/records/17854236),
DOI [10.5281/zenodo.17854236](https://doi.org/10.5281/zenodo.17854236),
`1_Challenge_Data_Supplement.zip` (published MD5
`fabfb9cbba5761d29d681e752a45df75`). Cite Marbach et al., *Wisdom of
crowds for robust gene network inference*, Nature Methods 9, 796–804 (2012).
Zenodo identifies the archive as [CC BY-ND 4.0](https://creativecommons.org/licenses/by-nd/4.0/).
No endorsement by the DREAM5 authors or data hosts is implied. Consult the
source terms before reusing or redistributing the data.

`data/Network1/input data/net1_expression_data.tsv`,
`data/Network1/input data/net1_transcription_factors.tsv`, and
`data/Network1/gold standard/DREAM5_NetworkInference_GoldStandard_Network1.tsv`
are byte-for-byte copies of the corresponding ZIP members, with SHA256 digests
pinned in `benchmark.py`. `data/Network1/dream5_net1.exp` is a technical,
value-preserving transposition of that expression table into SJARACNe's
gene-by-sample input format; its preparation code and SHA256 are recorded in
`benchmark.py` and `data/Network1/provenance.json`.

The `.adj` network predictions are **not** included here and are not part of
the original DREAM5 archive or its gold standard. The benchmark's inclusion
does not change the upstream dataset license or convert simulated network
edges into biological validation.
