# Changelog

## v0.4.0

First SJARACNe2 release, following legacy SJARACNe v0.3.0.

- Use 80% resampling without replacement and PCHIP-calibrated MI cutoffs.
- Apply native C++ K20 consensus and supporting-run mean MI before DPI, with
  an adjustable dynamic-tolerance offset (default 0.025).
- Read raw MI files sequentially with buffered input; make hash re-verification
  optional to avoid an extra full read.
- Provide WDL/CWL workflows through one CLI, automatic native builds and locked
  environments. Require Python 3.11 or newer; CI tests 3.11 and 3.14.

The explicit legacy profile remains available. The command remains `sjaracne`.
