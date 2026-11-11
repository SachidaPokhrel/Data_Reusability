# Installation

Prerequisites: Bash, Python 3.10+ (standard library), network access to official package
metadata/downloads, and writable project storage. Linux glibc/WSL2 and supported macOS
architectures are detected; unsupported combinations fail with a specific message.

```bash
./bootstrap.sh --check
./bootstrap.sh
```

The bootstrap checks existing micromamba/Conda first. Otherwise a pinned micromamba 2.0.5
package is chosen from conda-forge's official API, using its platform metadata and SHA256.
No package lacking a trusted checksum is executed. Environments go under CONDA_ROOT (default
project .conda). No sudo or edits to ~/.bashrc occur. Environment specs and exact-build exports
are described in ../envs/README.md. Smoke checks cover Python, Snakemake, MAFFT, R and R packages.

Installation and the synthetic Snakemake DAG completed on the source Linux x86_64 HPC.
The published exact-build exports describe that platform; other platforms are unverified.

For offline HPC, install/cache on an approved transfer node using the same shared project and
CONDA_ROOT. Bootstrap prepares Snakemake environment caches without running analysis. On a
compute node reuse those environments. Use --offline only when verified manager/packages are
already available. If a site provides modules, select a supported Python/Conda module explicitly.
