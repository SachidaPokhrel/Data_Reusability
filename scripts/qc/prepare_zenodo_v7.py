"""Build the pilot manifest and strict gene-level memberships for Zenodo record 17977529.

Run after extracting the verified assemblies and GFF3 files. The raw membership table
is preserved. This script writes derived TSVs and a separate pilot configuration.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / "src"))
from pe_upstream.common import load_config, path, rel, table, write_json
from pe_upstream.core import features

config = load_config(root / "config/pilot.yaml")
data = path(config, "data")
isolates = [f"Pe{i}" for i in range(1, 20)]
samples, parents = [], {}

for s in isolates:
    files = {
        "fasta": data / "external/genome_assembly" / f"{s}_assembly.fasta",
        "gff": data / "external/gene_annotation" / f"{s}_annotation.gff3",
        "repeats": data / "external/repeat_annotation" / f"{s}_repeat_annotation.gff",
    }
    for p in files.values():
        if not p.is_file() and not (s == "Pe6" and p == files["repeats"]):
            raise FileNotFoundError(p)
    genes, links = set(), defaultdict(set)
    for f in features(files["gff"]):
        if f.kind == "gene" and "ID" in f.attrs:
            genes.add(f.attrs["ID"])
        if f.kind in {"transcript", "mRNA"} and {"ID", "Parent"} <= f.attrs.keys():
            links[f.attrs["ID"]].update(f.attrs["Parent"].split(","))
    parents[s] = {tx: next(iter(ids)) for tx, ids in links.items()
                  if len(ids) == 1 and ids <= genes}
    samples.append({"isolate": s, "haplotype": "assembly_unphased",
                    "assembly_role": "haploid_primary",
                    "repeat_status": "available" if files["repeats"].is_file() else "missing_from_release",
                    **{k: rel(p) if p.is_file() else "" for k, p in files.items()}})

groups, assignments = {}, defaultdict(set)
with (data / "raw/Pe_orthogroups.tsv").open(newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader)
    expected = {f"{s}.{h}" for s in isolates for h in ("h1", "h2")}
    if len(header) != 39 or set(header[1:]) != expected:
        raise ValueError("Unexpected haplotype table header")
    columns = {name: i for i, name in enumerate(header)}
    for row in reader:
        if len(row) != len(header) or not row[0] or row[0] in groups:
            raise ValueError(f"Invalid orthogroup on line {reader.line_num}")
        group, agreed = row[0], {}
        for s in isolates:
            pair = []
            for h in ("h1", "h2"):
                call = row[columns[f"{s}.{h}"]].strip()
                if not call:
                    continue
                suffix = f"-{h}"
                gene = parents[s].get(call.removesuffix(suffix)) if call.endswith(suffix) else None
                if gene is None:
                    raise ValueError(f"Unresolved transcript: {group}, {s}.{h}, {call}")
                pair.append(gene)
                assignments[s, gene].add(group)
            if len(pair) == 2 and pair[0] == pair[1]:
                agreed[s] = pair[0]
        groups[group] = agreed

# Check reuse against all rows, including rows that fail the all-isolate filter.
memberships, audit = [], []
for group, agreed in groups.items():
    reused = [s for s, gene in agreed.items() if len(assignments[s, gene]) > 1]
    status = ("incomplete_or_discordant" if len(agreed) != 19 else
              "cross_group_assignment" if reused else "retained")
    audit.append({"orthogroup": group, "agreeing_isolates": len(agreed),
                  "reused_isolates": ";".join(reused), "status": status})
    if status == "retained":
        memberships.extend({"orthogroup": group, "isolate": s, "gene_id": gene}
                           for s, gene in agreed.items())

if not memberships:
    raise ValueError("No unambiguous all-19 groups retained")

manifest = root / "config/isolate_files.zenodo_v7.tsv"
table(manifest, samples)
table(data / "processed/Pe_gene_memberships.tsv", memberships)
table(data / "processed/Pe_gene_mapping_audit.tsv", audit)
config.pop("_config_path", None)
config["samples"] = rel(manifest)
config["orthology"] = {
    "file": "../processed/Pe_gene_memberships.tsv",
    "format": "long", "group_column": "orthogroup",
    "isolate_column": "isolate", "gene_column": "gene_id", "id_map": None,
}
write_json(root / "config/pilot.zenodo_v7.yaml", config)

print(f"Isolates mapped: {len(samples)}")
print(f"Retained all-19 groups: {len(memberships) // 19}")
print("Prepared config/pilot.zenodo_v7.yaml")
