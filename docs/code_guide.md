# Reading the code

Start with the [replay guide](reproduce_pilot.md) for commands and the [HTML report](../index.html) for the scientific result. The shell wrapper and Snakemake rules call the same Python stage functions. The original HPC pilot was run stage by stage.

| Task | Read first | What to look for |
|---|---|---|
| Download and verify | [acquisition.py](../src/pe_upstream/acquisition.py) | Pinned record, checksum validation, archive path checks |
| Map orthogroups to genes | [prepare_zenodo_v7.py](../scripts/qc/prepare_zenodo_v7.py) | Haplotype agreement through GFF3 parents and exclusion of reused gene assignments |
| Check inputs and choose loci | [quality.py](../src/pe_upstream/quality.py), [selection.py](../src/pe_upstream/selection.py) | All-19 core definition; six isolates selected by genotype distance; sampling before upstream measurements |
| Extract and align | [core.py](../src/pe_upstream/core.py), [sequences.py](../src/pe_upstream/sequences.py) | Zero-based, half-open coordinates, strand orientation, CDS phase, and extraction flags |
| Classify variation | [variation.py](../src/pe_upstream/variation.py), [repeats.py](../src/pe_upstream/repeats.py) | Matched callable pairs; separate coverage and identity; unresolved status; repeat overlap rather than insertion history |
| Summarize results | [analyze.R](../scripts/statistics/analyze.R), [figures.R](../scripts/visualization/figures.R) | Locus-level denominators, unresolved bounds, and conditional chromosome-block bootstrap |
| Review the candidate | [review_chr5_38.py](../scripts/qc/review_chr5_38.py) | Project the original windows onto the larger alignment; compare occupied columns before interpreting divergence |
| Track execution | [cli.py](../src/pe_upstream/cli.py), [common.py](../src/pe_upstream/common.py) | Stage contracts, file hashes, and audit checks |

The code comments explain decisions that can change the interpretation. For exact column meanings and thresholds, see the [data dictionary](data_dictionary.md) and [statistical methods](statistical_methods.md). The core and pilot definitions are separate: a core gene is present in all 19 source assemblies, while the divergence measurements compare six selected isolates.
