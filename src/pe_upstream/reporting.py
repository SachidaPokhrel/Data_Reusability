"""Reproducible reports from actual tables; no simulated biological conclusions."""
import html
import json
import shutil
import subprocess
from collections import Counter,defaultdict
from pathlib import Path
from .common import ROOT,path,rows,table,read_json,write_json,atomic_text,software,now,digest,rel,git_hash,markdown_page


def statistics(c):
    """Run R locus-level summaries, or the explicitly synthetic fallback if R is unavailable."""
    out=path(c,'results')/'statistics';out.mkdir(parents=True,exist_ok=True)
    master=path(c,'results')/'candidates/master_locus_table.tsv'
    if c['mode']=='test' and not shutil.which('Rscript'):
        counts=[];summary=[]
        data=rows(master)
        for w in c['upstream_windows_bp']:
            dd=[r for r in data if int(r['window'])==w];ct=Counter(r['candidate_status'] for r in dd)
            counts.extend({'window':w,'candidate_status':k,'n':ct[k]} for k in 'ABCDE')
            denom=ct['A']+ct['B']
            summary.append({'window':w,'sampled_loci':len(dd),'unresolved':ct['E'],
                'high_quality_coding_conserved':denom,'candidate_screen_count':ct['B'],
                'conditional_candidate_fraction':ct['B']/denom if denom else None,
                'bootstrap_status':'TEST_ONLY_R_UNAVAILABLE'})
        return [table(out/'prevalence.tsv',summary),table(out/'classification_counts.tsv',counts),
                write_json(out/'statistics_status.json',{'status':'synthetic_descriptive_only','R':'unavailable'})]
    subprocess.run(['Rscript',str(ROOT/'scripts/statistics/analyze.R'),str(master),
        str(out/'sensitivity_locus.tsv'),str(out)],check=True)
    return [out/name for name in ['prevalence.tsv','classification_counts.tsv','sensitivity_summary.tsv',
        'repeat_association_coefficients.tsv','model_status.tsv','R_sessionInfo.txt','statistics_status.json']]


def figures(c):
    """Render the main R figures from measured tables and prespecified thresholds."""
    directory=path(c,'figures')
    if c['mode']=='test' and not shutil.which('Rscript'):
        return [write_json(directory/'figure_status.json',{'status':'skipped_R_unavailable',
            'reason':'No mock scientific figures generated'})]
    subprocess.run(['Rscript',str(ROOT/'scripts/visualization/figures.R'),
        str(path(c,'results')/'candidates/master_locus_table.tsv'),
        str(path(c,'results')/'statistics'),str(directory),str(c['_config_path'])],check=True)
    return list(directory.glob('figure*'))


def extra_figures(c):
    """Plot isolate-level differences and nominated screen hits without implying validation."""
    if c['mode']=='test':
        try:import matplotlib
        except ImportError:
            return [write_json(path(c,'figures')/'candidate_figure_status.json',
                {'status':'skipped_matplotlib_unavailable','figures':[]})]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from .core import fasta,load_genes
    from .repeats import repeat_records
    from .sequences import file_id
    directory=path(c,'figures');out=path(c,'results')
    master=rows(out/'candidates/master_locus_table.tsv');selection=read_json(out/'orthogroups/selection.json')
    outputs=[];reviews=[]
    sample_manifest={r['isolate']:r for r in rows(out/'qc/isolate_manifest.tsv')}
    coordinate_rows=rows(out/'upstream/upstream_coordinates.tsv')
    model_cache={};repeat_cache={}
    for w in c['upstream_windows_bp']:
        reliable=[r for r in master if int(r['window'])==w and r['candidate_status']!='E']
        top=sorted(reliable,key=lambda r:float(r['upstream_divergence']),reverse=True)[:20]
        if top:
            heat=[]
            for r in top:
                obj=read_json(out/'variation/intermediate'/f'{file_id(r["orthogroup_id"])}.json')
                pp=obj[f'upstream_{w}']['pairs'];means=[]
                for s in selection['isolates']:
                    vals=[p['combined_divergence'] for p in pp if p['valid'] and s in (p['isolate_a'],p['isolate_b'])]
                    means.append(sum(vals)/len(vals) if vals else float('nan'))
                heat.append(means)
            fig,ax=plt.subplots(figsize=(8,max(4,len(top)*.25)))
            image=ax.imshow(np.ma.masked_invalid(np.array(heat)),aspect='auto',cmap='magma')
            ax.set_xticks(range(len(selection['isolates'])),selection['isolates'],rotation=90)
            ax.set_yticks(range(len(top)),[r['orthogroup_id'] for r in top],fontsize=7)
            ax.set_title(f'Mean upstream divergence to other isolates; {w} bp')
            fig.colorbar(image,ax=ax,label='Combined alignment divergence');fig.tight_layout()
            for ext in ['png','pdf']:
                p=directory/f'figure4_heatmap_{w}.{ext}';fig.savefig(p,dpi=300);outputs.append(p)
            plt.close(fig)
        candidates=sorted([r for r in reliable if r['candidate_status']=='B'],
            key=lambda r:float(r['upstream_divergence']),reverse=True)[:c['max_candidate_loci']]
        # Aligned bases form a common coordinate framework; gaps are shown explicitly.
        for rank,row in enumerate(candidates,1):
            og=row['orthogroup_id'];aligned=fasta(out/'upstream/alignments'/f'{file_id(og)}.upstream_{w}.fa')
            seqs=[aligned[s] for s in selection['isolates']];colors={'A':0,'C':1,'G':2,'T':3,'-':4,'N':5}
            matrix=np.array([[colors.get(ch,5) for ch in seq] for seq in seqs])
            from matplotlib.colors import ListedColormap
            fig,(locus,ax)=plt.subplots(2,1,figsize=(12,7),gridspec_kw={'height_ratios':[1,1.2]})
            ax.imshow(matrix,aspect='auto',interpolation='nearest',
                cmap=ListedColormap(['#4daf4a','#377eb8','#ff7f00','#e41a1c','#eeeeee','#777777']),vmin=0,vmax=5)
            local=[r for r in coordinate_rows if r['orthogroup']==og and int(r['window'])==w]
            by_isolate={r['isolate']:r for r in local}
            gene_limit=w*.35
            for y,isolate in enumerate(selection['isolates']):
                info=by_isolate[isolate];a,b=int(info['start']),int(info['end'])
                locus.plot([-int(info['actual_length']),0],[y,y],color='gray',lw=2)
                locus.arrow(0,y,gene_limit,0,width=.04,head_width=.22,head_length=w*.04,
                    color='steelblue',length_includes_head=True)
                if isolate not in repeat_cache:
                    repeat_cache[isolate]=repeat_records(sample_manifest[isolate]['repeats'],c['repeat_te_classes'])
                for rep in repeat_cache[isolate]:
                    if rep['chrom']!=info['chrom']:continue
                    start=max(a,rep['start']);end=min(b,rep['end'])
                    if start>=end:continue
                    anchor=int(info['gene_start']) if info['strand']=='+' else int(info['gene_end'])
                    left,right=(start-anchor,end-anchor) if info['strand']=='+' else (anchor-end,anchor-start)
                    locus.plot([left,right],[y,y],color='darkorange' if rep['is_TE'] else 'dimgray',lw=6)
            locus.set_yticks(range(len(seqs)),selection['isolates'])
            locus.set_xlim(-w,gene_limit*1.1)
            locus.set(xlabel='Gene-boundary coordinates (bp); gene arrows truncated, not to scale',
                title=f'{og}: transcription-oriented locus; orange classified TE, gray other repeat')
            reviews.append({'orthogroup':og,'window':w,'rank':rank,'evidence':'moderate',
                'synteny':'conserved_immediate_core_neighbors','strand':'transcription_oriented',
                'assembly_continuity':'passed_window_QC','gene_model':'passed_CDS_QC',
                'TE_insertion':'not_established','graph_support':'not_run_by_main_workflow','manual_review':'pending'})
            ax.set_yticks(range(len(seqs)),selection['isolates'])
            ax.set(xlabel='Upstream multiple-alignment column (transcription-oriented)',
                title=f'{og}, {w} bp: assembly alignment; gaps light gray, ambiguity dark gray')
            fig.text(.5,.01,'Candidate screen; TE insertion history and regulatory function are not established.',ha='center',fontsize=8)
            fig.tight_layout(rect=[0,.04,1,1])
            for ext in ['png','pdf']:
                p=directory/f'figure5_candidate_{rank}_{w}.{ext}';fig.savefig(p,dpi=300);outputs.append(p)
            plt.close(fig)
    outputs.append(table(out/'candidates/detailed_candidate_review.tsv',reviews,
        ['orthogroup','window','rank','evidence','synteny','strand','assembly_continuity',
         'gene_model','TE_insertion','graph_support','manual_review']))
    outputs.append(write_json(directory/'candidate_figure_status.json',{'status':'complete',
        'files':[rel(p) for p in outputs],'no_candidate_plot_meaning':'No qualifying screen candidates or no valid alignments'}))
    return outputs


def report(c):
    """Write the pilot report and record the gate decision without starting full mode."""
    out=path(c,'results');reports=path(c,'reports');figs=path(c,'figures')
    selection=read_json(out/'orthogroups/selection.json');master=rows(out/'candidates/master_locus_table.tsv')
    stats=rows(out/'statistics/prevalence.tsv');evidence=rows(out/'synteny/locus_evidence.tsv')
    label='TEST — SYNTHETIC DATA' if c['mode']=='test' else c['mode'].upper()
    lines=[f'# {label}: Conserved Genes, Variable Neighborhoods','',
        '## Abstract','',f'This {label.lower()} run screened {len(selection["orthogroups"])} loci in '
        f'{len(selection["isolates"])} assembly representations for coding versus upstream sequence conservation. '
        'Candidate calls are computational screens; regulatory function and TE insertion history are not established.',
        '', '## Background','',
        'A shared gene need not have a conserved coding sequence or upstream neighborhood. '
        'This reanalysis tests whether gene-centered summaries conceal adjacent sequence variation.',
        '', '## Original data source','',
        ('Synthetic fixtures generated by this repository; no biological findings.' if c['mode']=='test' else
         'Skiadas et al. (2026), PLOS Biology, doi:10.1371/journal.pbio.3003596. '
         'Pinned Zenodo record 17977529; NCBI BioProject PRJNA772192.'),
        '', '## New data-reuse question','',
        'How often does high coding conservation in a single-copy core orthogroup fail to extend into upstream DNA?',
        '', '## Methods','',
        f'Core definition: {c["core_population_size"]} manifest isolates. Selection: {selection["isolate_selection_basis"]}. '
        f'Gene sampling seed {c["random_seed"]}; sampled before upstream QC. CDS divergence threshold '
        f'{c["max_cds_divergence"]}; variable-upstream threshold {c["min_upstream_divergence"]}. '
        f'Callable fraction threshold {c["min_pair_callable_fraction"]}; minimum matched-pair fraction '
        f'{c["min_valid_pair_fraction"]}. Annotation-derived windows: {c["upstream_windows_bp"]} bp. '
        'CDS and upstream are compared on matching valid isolate pairs.',
        '', '## Data and QC summary','',
        f'All-manifest single-copy core groups: {selection["all_manifest_core_genes"]}. '
        f'CDS/synteny-eligible groups: {selection["eligible_cds_synteny_genes"]}. '
        f'Sampled groups: {len(selection["orthogroups"])}. Upstream failures remain unresolved.',
        '', '## Results','', '| Window | Sampled loci | Unresolved | High-quality coding-conserved | Candidate screens | Conditional fraction |',
        '|---|---|---|---|---|---|']
    for r in stats:
        lines.append('| '+' | '.join(str(r.get(k,'')) for k in ['window','sampled_loci','unresolved',
            'high_quality_coding_conserved','candidate_screen_count','conditional_candidate_fraction'])+' |')
    lines+=['','## Statistical analysis','',
        'Statistical outputs include denominator-specific prevalence and threshold sensitivity. '
        'Chromosome-block bootstrap intervals, when supported, are conditional on sampled isolates. '
        'Exploratory regressions require adequate locus/chromosome counts and full-rank predictors. '
        'Pairwise comparisons are not treated as independent replicates. See prevalence.tsv, '
        'model_status.tsv and coefficient tables for effect sizes and BH-adjusted p-values when estimated.']
    for r in stats:
        if r.get('ci_lower') and r.get('ci_upper'):
            lines.append(f'\nWindow {r["window"]}: conditional chromosome-block 95% interval '
                         f'{r["ci_lower"]}–{r["ci_upper"]}.')
    lines+=['','## Figures','']
    import os
    for p in sorted(figs.glob('*.png')):
        lines.append(f'![{p.stem}]({os.path.relpath(p,reports)})\n')
    lines+=['','## Interpretation','',
        'Class B is a candidate conserved-CDS/variable-upstream pattern, not proof of altered regulation. '
        'Assembly-derived indels with conserved flanks support structural differences, but their ancestral '
        'direction and association with particular TE events require additional validation.',
        '', '## Limitations','',
        'Gene boundaries are not verified transcription start sites. Equal-length windows may truncate '
        'homologous distal sequence after insertions. Noncoding and coding regions experience different '
        'constraints. Haploid assembly representations omit within-isolate variation. Filters enrich for '
        'well-annotated, collinear genes. Race-selected isolates are not random population samples. '
        'Graph validation and expression follow-up are optional, not presumed complete.',
        '', '## Conclusions','',
        'The tables quantify the specified sequence patterns conditional on the input assemblies and QC. '
        'Strong biological conclusions require manual candidate review and orthogonal validation.',
        '', '## Future directions','',
        'Review 3–5 loci using syntenic flanks, graph paths or local genome alignments; validate annotation '
        'boundaries and repeat classes; use comparable expression data only if its design supports inference.',
        '', '## Software and provenance','']
    sw=software(c)
    for r in rows(sw):lines.append(f'- {r["software"]}: {r["version"]}')
    lines+=['','## References','',
        '- Skiadas et al. (2026). https://doi.org/10.1371/journal.pbio.3003596',
        '- Dataset: https://doi.org/10.5281/zenodo.17977529','']
    text='\n'.join(lines);name=f'{c["mode"]}_report'
    md=atomic_text(reports/f'{name}.md',text)
    ht=atomic_text(reports/f'{name}.html',markdown_page(text,label+' report'))
    summary=atomic_text(reports/'week12_presentation_summary.md',
        f'# {label} presentation summary\n\nQuestion: Does coding conservation extend upstream?\n\n'
        f'Analyzed {len(selection["isolates"])} assembly representations and {len(selection["orthogroups"])} sampled loci.\n\n'
        'Use the generated prevalence table with explicit denominators. Candidate screens require manual validation.\n\n'
        +('Synthetic validation only; no biological findings.\n' if c['mode']=='test' else 'No expression or causal regulatory claims are supported by this analysis alone.\n'))
    gate={'status':'pass' if c['mode']=='pilot' and len(selection['isolates'])==6 and
        len(selection['orthogroups'])>=50 and any(int(r['unresolved'])<int(r['sampled_loci'])/2 for r in stats) else 'not_passed',
        'mode':c['mode'],'config_sha256':digest(c['_config_path']),
        'manifest_sha256':digest(out/'qc/isolate_manifest.tsv'),
        'input_files':{rel(r[k]):digest(r[k]) for r in rows(out/'qc/isolate_manifest.tsv') for k in ['fasta','gff','repeats'] if r[k]},
        'reason':'requires real six-isolate pilot, >=50 sampled loci, <50% unresolved in at least one window',
        'timestamp_utc':now()}
    gp=write_json(reports/'analysis_gate.json',gate)
    return [md,ht,summary,gp,sw]
