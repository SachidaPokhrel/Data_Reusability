# Data sources and source review

- Paper: https://doi.org/10.1371/journal.pbio.3003596
- Versioned record: https://zenodo.org/records/17977529
- REST API: https://zenodo.org/api/records/17977529
- BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA772192
- Authors' code: https://github.com/TeamMGE/Skiadas2025_pangenome_haplotypes

The biological pilot on HPC downloaded and checksum-verified the versioned record 17977529.
See the [verified files](../publication/verified_downloads.tsv). Archive inspection confirmed
19 assembly and gene annotation files and 18 repeat annotation files. Pe6 repeat data are absent.
The haplotype membership table was reconciled through transcript-to-gene links in the deposited
GFF3; the resulting pilot configuration is versioned separately.

Acquisition starts from the API and saves the entire response to data/metadata/zenodo_record.json.
Download links are read from API fields; required entries must have HTTPS URLs, sizes and
recognized checksum algorithms. Raw archives remain unchanged. The PHYLIP genotype matrix is
additionally requested for diversity-based pilot selection; graph data are optional.

Genome FASTA paths, gene/repeat annotation formats and gene identifier conventions must be
inspected before populating config/isolate_files.tsv. A version change requires a new record,
configuration review and revalidation, rather than silently following a concept DOI.

The paper's methods describe an initial annotation followed by a comparative common annotation.
Select the deposited annotation set that actually matches orthogroup identifiers; filename similarity
alone is insufficient. See the [primary methods](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3003596).
