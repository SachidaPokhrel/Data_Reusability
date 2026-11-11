#!/usr/bin/env Rscript
# One observation per locus/window. Isolate pairs are never independent replicates.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop('Usage: analyze.R MASTER_TSV SENSITIVITY_TSV OUTPUT_DIRECTORY')
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(ggplot2)
  library(sandwich);library(lmtest);library(jsonlite)})
master <- read_tsv(args[1], show_col_types=FALSE, na=c('', 'NA'))
sens <- read_tsv(args[2], show_col_types=FALSE)
out <- args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(2026)
summary <- master %>% group_by(window) %>% summarise(
  sampled_loci=n(), unresolved=sum(candidate_status=='E'),
  coding_conserved=sum(coding_conserved_independent==1),
  coding_conserved_unresolved=sum(coding_conserved_independent==1 & candidate_status=='E'),
  high_quality_coding_conserved=sum(candidate_status %in% c('A','B')),
  candidate_screen_count=sum(candidate_status=='B'),.groups='drop') %>% mutate(
  conditional_candidate_fraction=ifelse(high_quality_coding_conserved>0,candidate_screen_count/high_quality_coding_conserved,NA_real_),
  all_coding_conserved_lower_bound=ifelse(coding_conserved>0,candidate_screen_count/coding_conserved,NA_real_),
  all_coding_conserved_upper_bound=ifelse(coding_conserved>0,(candidate_screen_count+coding_conserved_unresolved)/coding_conserved,NA_real_))
# Block bootstrap conditional on the sampled isolates and QC-passing locus set.
intervals <- lapply(sort(unique(master$window)),function(w) {
  d <- master %>% filter(window==w,candidate_status %in% c('A','B'))
  blocks <- split(d,d$chrom)
  if (length(blocks)<5) return(data.frame(window=w,ci_lower=NA_real_,ci_upper=NA_real_,bootstrap_status='fewer_than_five_chromosomes'))
  boot <- replicate(2000,{
    b <- bind_rows(blocks[sample(seq_along(blocks),length(blocks),replace=TRUE)])
    mean(b$candidate_status=='B')})
  data.frame(window=w,ci_lower=unname(quantile(boot,.025)),ci_upper=unname(quantile(boot,.975)),
             bootstrap_status='chromosome_block_percentile_conditional_on_isolates')
}) %>% bind_rows()
summary <- left_join(summary,intervals,by='window')
write_tsv(summary,file.path(out,'prevalence.tsv'))
classes <- master %>% count(window,candidate_status) %>% complete(window,candidate_status=c('A','B','C','D','E'),fill=list(n=0))
write_tsv(classes,file.path(out,'classification_counts.tsv'))
sens_summary <- sens %>% group_by(window,cds_threshold,upstream_threshold) %>% summarise(
  resolved_loci=sum(resolved),coding_conserved=sum(coding_conserved),candidates=sum(candidate),.groups='drop') %>%
  mutate(candidate_fraction=ifelse(coding_conserved>0,candidates/coding_conserved,NA_real_))
write_tsv(sens_summary,file.path(out,'sensitivity_summary.tsv'))
coefficients <- list();model_status <- list()
for (w in sort(unique(master$window))) {
  d <- master %>% filter(window==w,candidate_status!='E')
  formula <- upstream_divergence ~ CDS_divergence + upstream_repeat_fraction +
    log1p(intergenic_length) + local_repeat_density + log1p(distance_to_contig_end) +
    missing_data_fraction + factor(chrom)
  d <- d[complete.cases(model.frame(formula,d,na.action=na.pass)),]
  status <- 'not_enough_loci_or_chromosomes'
  if (nrow(d)>=50 && length(unique(d$chrom))>=8) {
    fit <- lm(formula,data=d)
    if (fit$rank==length(coef(fit)) && nrow(d)>=5*fit$rank) {
      test <- coeftest(fit,vcov.=vcovCL(fit,cluster=d$chrom,type='HC1'),df=length(unique(d$chrom))-1)
      crit <- qt(.975,df=length(unique(d$chrom))-1)
      coefficients[[length(coefficients)+1]] <- data.frame(window=w,term=rownames(test),
        effect_size=test[,1],standard_error=test[,2],p_value=test[,4],
        ci_lower=test[,1]-crit*test[,2],ci_upper=test[,1]+crit*test[,2],
        n_loci=nrow(d),n_chromosomes=length(unique(d$chrom)))
      status <- 'exploratory_OLS_chromosome_clustered'
    } else status <- 'rank_deficiency_or_insufficient_observations_per_parameter'
  }
  model_status[[length(model_status)+1]] <- data.frame(window=w,status=status)
}
if (length(coefficients)>0) {
  coefs <- bind_rows(coefficients) %>% mutate(p_value_BH=p.adjust(p_value,method='BH'))
  write_tsv(coefs,file.path(out,'repeat_association_coefficients.tsv'))
} else write_tsv(data.frame(window=integer(),term=character(),effect_size=numeric(),p_value_BH=numeric()),file.path(out,'repeat_association_coefficients.tsv'))
write_tsv(bind_rows(model_status),file.path(out,'model_status.tsv'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
write_json(list(status='complete',inference='exploratory; no expression or causal regulatory inference',
  bootstrap_seed=2026,bootstrap_iterations=2000),file.path(out,'statistics_status.json'),pretty=TRUE,auto_unbox=TRUE)
