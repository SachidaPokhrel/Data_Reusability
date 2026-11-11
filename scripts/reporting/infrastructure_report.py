#!/usr/bin/env python3
"""Regenerate the PARTIAL report from observed infrastructure/stage outcomes only."""
import html
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from pe_upstream.common import ROOT,read_json,write_json,atomic_text,now,git_hash,table,digest,rel,markdown_page


def main():
    folder=ROOT/'reports/provenance/stages'
    records=[read_json(p) for p in sorted(folder.glob('*.json'))]
    lines=['# PARTIAL: Conserved Genes, Variable Neighborhoods','',
        '## Abstract','',
        'The repository and a synthetic analysis harness were implemented and tested. Installation '
        'and biological data acquisition failed because DNS resolution was unavailable. No biological '
        'isolates, genes or candidate loci have been analyzed; those quantities are not estimable.',
        '', '## Background','',
        'Core membership and coding conservation do not guarantee upstream conservation. '
        'The proposed re-use tests the scope of gene-centered conservation labels.',
        '', '## Original data source','',
        'Skiadas et al. (2026), PLOS Biology 24:e3003596. '
        'https://doi.org/10.1371/journal.pbio.3003596; pinned dataset '
        'https://doi.org/10.5281/zenodo.17977529; PRJNA772192.',
        '', '## New data-reuse question','',
        'Among coding-conserved, all-19 single-copy core genes, how often does upstream DNA fail '
        'prespecified conservation criteria?',
        '', '## Methods','',
        'The implemented stages validate membership and annotation mapping, select diversity-based '
        'pilot isolates and seeded random loci, extract strand-aware CDS/upstream sequence, align '
        'and summarize callable distances, annotate repeats, classify uncertainty, and generate reports. '
        'R/MAFFT/Snakemake integration code remains unexecuted here. Synthetic tests use an explicitly '
        'TEST-only alignment stub and standard-library descriptive fallback when tools are absent.',
        '', '## Data and QC summary','',
        'No raw biological files were retrieved. Actual source schemas, isolate identifiers, haplotype '
        'representation and archive contents remain unverified. No biological QC summary exists.',
        '', '## Results','', '| Stage | Status | Observed outcome |','|---|---|---|']
    for row in records:
        lines.append(f'| {row["stage"]} | {row["status"]} | {row["details"]} |')
    tests=ROOT/'reports/provenance/test_results.json'
    if tests.exists():
        t=read_json(tests);lines+=['',f'Automated tests: {t["tests_run"]} run; '
            f'{t["passed"]} passed; {t["skipped"]} skipped; {t["failed"]} failed. '
            'These are software validation outcomes, not biological results.']
    lines+=['','## Statistical analysis','',
        'Biological prevalence, effect sizes, confidence intervals and candidate counts are unavailable. '
        'The R analysis code specifies conditional prevalence, unresolved bounds, chromosome-block '
        'bootstrap and prespecified threshold sensitivity; none has run on biological data.',
        '', '## Figures','',
        'No scientific figures were generated from biological data. No mock figures are substituted.',
        '', '## Interpretation','',
        'Synthetic validation supports selected extraction, accounting, transfer and audit behaviors. '
        'It does not validate source compatibility, actual tool integrations or hypotheses H1–H3.',
        '', '## Limitations','',
        'DNS, absent software and unavailable source data prevent scientific execution. Further '
        'limitations include annotation-derived anchors, assembly paths versus diploid variation, '
        'fixed-window homology, repeat-class uncertainty and related, race-selected isolates.',
        '', '## Conclusions','',
        'Infrastructure and offline tests are available; the biological research question remains unanswered.',
        '', '## Future directions','',
        'Restore permitted network access; bootstrap and smoke-test environments; retrieve REST '
        'metadata/files; inspect schemas and populate reviewed isolate/ID mappings; run and validate '
        'the six-isolate pilot; expand only after the pilot gate passes.',
        '', '## Software and execution','',
        'See reports/provenance/software_versions.tsv, bootstrap_status.json, repository_validation.json '
        'and test_results.json. Unavailable software is labeled UNAVAILABLE, never an installed version.',
        '', '## References','',
        '- Skiadas et al. (2026). https://doi.org/10.1371/journal.pbio.3003596',
        '- Zenodo dataset. https://doi.org/10.5281/zenodo.17977529','']
    text='\n'.join(lines)
    report=ROOT/'reports/pilot'
    md=atomic_text(report/'pilot_report.md',text)
    ht=atomic_text(report/'pilot_report.html',markdown_page(text,'PARTIAL project report'))
    atomic_text(report/'week12_presentation_summary.md','# Week 12 — PARTIAL project status\n\n'
        'Question: Does conservation of single-copy core coding DNA extend upstream?\n\n'
        'The workflow and synthetic tests are implemented. Real-data analysis is blocked by DNS and '
        'missing tools. No biological findings or candidate loci can yet be presented.\n\n'
        'Next: verified API retrieval, reviewed sample mapping, six-isolate pilot and candidate validation.\n')
    manifest={'status':'partial','timestamp_utc':now(),'git_commit_at_generation':git_hash(),
        'biological_analysis_completed':False,'isolates_analyzed':None,'genes_analyzed':None,
        'high_confidence_candidates':None,'stages':records,
        'outputs':{rel(p):digest(p) for p in [md,ht]}}
    write_json(ROOT/'reports/provenance/run_manifest.json',manifest)
    table(ROOT/'reports/provenance/input_checksums.tsv',
        [{'path':rel(p),'sha256':digest(p),'role':'configuration_or_status_not_biological_data'}
         for p in sorted((ROOT/'config').glob('*.yaml'))+sorted(folder.glob('*.json'))])
    print('Regenerated PARTIAL report from recorded outcomes; no biological results invented.')

if __name__=='__main__':main()
