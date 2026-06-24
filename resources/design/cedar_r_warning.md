# CEDAR warning in Snakemake R jobs

## Summary

R rules run via Snakemake's `conda:` directive produce a cosmetic warning:

```
Warning messages:
1: No CEDAR controlled environment in .libPaths():
 - /scratch/.../scprocess_envs/<hash>_/lib/R/library
2: Aborted CEDAR session initialization, contact #help-cedar on slack
```

CEDAR is an in-house R package management system. The warning means CEDAR's startup hook ran but couldn't find its library path — which is the intended outcome of the `R_LIBS_SITE=""` fix in the scprocess CLI (lines 39-40 in `scprocess`). No CEDAR packages are loaded; the warning is cosmetic.

## Existing mitigation

The scprocess CLI sets `R_LIBS_SITE=""` and `R_LIBS_USER=""` in the environment passed to the Snakemake subprocess. This prevents conda R from falling back to system site-library paths.

## What we ruled out (2026-06-24)

| Source checked | Result |
|---|---|
| Compute node env vars (`R_PROFILE`, `R_PROFILE_SITE`, `R_LIBS_SITE`, etc.) | All empty — submitted SLURM job to verify |
| `Rprofile.site` in conda env | Does not exist |
| `~/.Rprofile` | Does not exist |
| Project `.Rprofile` (loads workflowr) | Tested directly on compute node — no CEDAR warning |
| Any file in the rlibs conda env | `grep -ri cedar` across entire env tree — nothing |
| Direct `Rscript` from conda env on compute node | No warning produced |

The warning **only** appears when R is launched via Snakemake's `conda:` directive, not when the same conda env's `Rscript` is invoked directly.

## Likely cause

Snakemake's internal conda activation mechanism sets environment variables or sources scripts that re-introduce a system path, which triggers CEDAR's hook. The `R_LIBS_SITE=""` set by the scprocess CLI may be overridden during Snakemake's conda env activation step.

## Potential fixes (not yet implemented)

- Set `R_PROFILE_SITE=""` alongside the existing `R_LIBS` vars in the scprocess CLI to suppress any site profile
- Investigate Snakemake's conda activation to find where `R_LIBS_SITE` is being re-set
- Add `R_LIBS_SITE=""` to `envmodules:` or shell preamble in individual rules (less desirable — per-rule rather than global)
