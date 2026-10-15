# Representative NetBID2 recurrence-QC reports

These two single-file reports are the representative QC views for the proposed
BRCA100 operating points. Their inputs are the final \(K=6\) consensus networks
built from the frozen post-DPI recurrence sweep; generating the HTML did not
rerun SJARACNe inference.

| Network | \(p_b\) | AP-MI cutoff | \(K\) | Edges | Adjusted scale-free \(R^2\) | Report | SHA-256 |
|---|---:|---:|---:|---:|---:|---|---|
| TF | \(10^{-3}\) | 0.147322 | 6 | 416,408 | 0.879250 | [HTML](k006_tf_netbid2_qc.html) | `8fbaaf0f341e11a74ab1ba26e01c497b563fb6204ad3fa330235dd00065dd62b` |
| SIG | \(5\times10^{-4}\) | 0.164467 | 6 | 739,958 | 0.864118 | [HTML](k006_sig_netbid2_qc.html) | `d8bb181014ba2e1f9d33be6b58bfa8bd8bb4c3f3fc613ed0da97872ccf8761cf` |

## Provenance

- Frozen design fingerprint:
  `8948b3d2be1a9e1a69591131b0475b1c5265475b16665caca34be4f72d0a6e8e`
- Frozen design SHA-256:
  `1228de2bd5f2ae0e59bfaadefc43d8a9a080c68e6ae7313713f1bba762269216`
- R 4.4.3 (2025-02-28)
- NetBID2 2.2.0, remote commit
  `5defa454d600b94f5dd6d1f9f4428f99759a6821`
- igraph 2.3.3

The per-arm provenance directories retain the R Markdown source, exact TSV
outputs, environment record, network and summary manifests, and HTML-generation
manifest. The materialized 31 MB TF and 55 MB SIG consensus networks are
intentionally omitted; their hashes and byte sizes remain in the manifests.
The exact generation scripts are under `provenance/scripts/`.

## Viewing and limits

GitHub generally displays HTML as source. Download the raw file, verify it
against the package-level [SHA256SUMS](../SHA256SUMS), and open it locally. Plots are
embedded in each file. MathJax is loaded from `mathjax.rstudio.com`, so
equations may not typeset without network access.

These reports are topology QC. Adjusted \(R^2\) is descriptive and is not a
selection criterion, biological validation, edge-level p-value, or FDR
estimate.
