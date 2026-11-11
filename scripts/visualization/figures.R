#!/usr/bin/env Rscript
# Reproduce PNG/PDF plots from locus tables and the recorded thresholds.
# Figure 1 omits class E for readability; Figure 2 retains those loci.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4)stop('Usage: figures.R MASTER_TSV STATISTICS_DIR FIGURES_DIR WINDOW_CONFIG')
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(ggplot2);library(jsonlite)})
d <- read_tsv(args[1],show_col_types=FALSE,na=c('','NA'));stats <- args[2];out <- args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
config <- fromJSON(args[4])
theme_set(theme_bw(base_size=11))
save_plot <- function(p,name,width=8,height=5) {
 ggsave(file.path(out,paste0(name,'.png')),p,width=width,height=height,dpi=300)
 ggsave(file.path(out,paste0(name,'.pdf')),p,width=width,height=height)
}
labels <- c(A='Conserved CDS / conserved upstream',B='Conserved CDS / variable upstream',
 C='Variable CDS / conserved upstream',D='Variable CDS / variable upstream',E='Unresolved / intermediate')
p <- ggplot(d %>% filter(candidate_status!='E'),aes(CDS_divergence,upstream_divergence,color=upstream_repeat_fraction))+
 geom_point(alpha=.6,size=1.5)+facet_wrap(~window,labeller=label_both)+scale_color_viridis_c(na.value='gray')+
 geom_vline(xintercept=config$max_cds_divergence,linetype=2)+geom_hline(yintercept=config$min_upstream_divergence,linetype=2)+
 labs(x='CDS divergence (substitutions + gap bases)',y='Upstream divergence',color='Repeat fraction',
 caption='Unresolved loci are retained in Figure 2; windows are not independent.')
save_plot(p,'figure1_CDS_vs_upstream')
classes <- read_tsv(file.path(stats,'classification_counts.tsv'),show_col_types=FALSE)
p <- ggplot(classes,aes(factor(window),n,fill=candidate_status))+geom_col(position='fill')+
 scale_fill_brewer(palette='Set2',labels=labels,drop=FALSE)+labs(x='Upstream window (bp)',y='Proportion of sampled loci',fill='Class',
 caption='E includes extraction failures, uncertain homology, low coverage and intermediate upstream divergence.')
save_plot(p,'figure2_conservation_classes',10,5)
p <- ggplot(d %>% filter(candidate_status!='E'),aes(upstream_repeat_fraction,upstream_divergence))+
 geom_point(alpha=.5,color='darkcyan',size=1.5)+facet_wrap(~window,labeller=label_both)+
 labs(x='Mean upstream repeat fraction',y='Upstream divergence',caption='Association does not establish TE insertion or causation.')
save_plot(p,'figure3_repeat_association')
sens <- read_tsv(file.path(stats,'sensitivity_summary.tsv'),show_col_types=FALSE)
p <- ggplot(sens,aes(factor(cds_threshold),factor(upstream_threshold),fill=candidate_fraction))+
 geom_tile()+facet_wrap(~window,labeller=label_both)+scale_fill_viridis_c(na.value='gray80')+
 labs(x='CDS conservation threshold',y='Upstream variability threshold',fill='Candidate fraction',
 caption='QC is held fixed. Fractions conditional on callable, homology-supported loci.')
save_plot(p,'figure6_threshold_sensitivity')
write_json(list(status='complete',figures=c(1,2,3,6),optional_figures=c(4,5)),
 file.path(out,'figure_status.json'),pretty=TRUE,auto_unbox=TRUE)
