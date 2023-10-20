AOPEP_gwas<-subset(MiBioGen_clean,CHR==9 & POS>94726699-100000 & POS<95150224+100000)
AOPEP_dat<-format_data(AOPEP_gwas,
                          snp_col = "rsID",
                          beta_col = "beta",
                          se_col = "SE",
                          pval_col = "P",
                          samplesize_col = "N",
                          effect_allele_col = "eff.allele",
                          other_allele_col = "re.allele")
AOPEP_ALM_gwas <- extract_outcome_data(snps = AOPEP_dat$SNP,
                                       outcomes = 'ebi-a-GCST90000025')
dat <- harmonise_data(AOPEP_dat,AOPEP_ALM_gwas)
res3 <- mr(dat)
res3

mr_heterogeneity(dat)

run_mr_presso(harm_rt,NbDistribution = 1000)

mr_pleiotropy_test(dat)

singlesnp_res<- mr_singlesnp(dat)
View(singlesnp_res)
singlesnpOR=generate_odds_ratios(singlesnp_res)

sen_res<- mr_leaveoneout(dat)
View(sen_res)

p1 <- mr_scatter_plot(res3, dat)
p1[[1]]

p2 <- mr_forest_plot(singlesnp_res)
p2[[1]]

p3 <- mr_leaveoneout_plot(sen_res)
p3[[1]]

res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(singlesnp_res)
p4[[1]]

expo_rt<- read_exposure_data(
  filename = "Lachnospiraceae.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "eff.allele",
  other_allele_col = "re.allele",
  pval_col = "P",
  samplesize_col = "N")


expo_rt<- expo_rt[expo_rt$pval.exposure < 1e-06,]
expo_rt <- clump_data(expo_rt,clump_kb = 100,clump_r2 = 0.01)

outc_rt <- read_outcome_data(
  snps = expo_rt$SNP,
  filename = "AOPEP.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "Freq",
  pval_col = "p")

dat <- harmonise_data(expo_rt,outc_rt)
res3 <- mr(dat)
res3

het <- mr_heterogeneity(dat)

ple <- mr_pleiotropy_test(dat)

singlesnp_res<- mr_singlesnp(dat)
singlesnpOR=generate_odds_ratios(singlesnp_res)

sen_res<- mr_leaveoneout(dat)

p1 <- mr_scatter_plot(res3, dat)
p1[[1]]

p2 <- mr_forest_plot(singlesnp_res)
p2[[1]]

p3 <- mr_leaveoneout_plot(sen_res)
p3[[1]]

res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(singlesnp_res)
p4[[1]]

AOPEP_gwas<-subset(MiBioGen_clean,CHR==9 & POS>94726699-100000 & POS<95150224+100000)
AOPEP_dat<-format_data(AOPEP_gwas,
                       snp_col = "rsID",
                       beta_col = "beta",
                       se_col = "SE",
                       pval_col = "P",
                       samplesize_col = "N",
                       effect_allele_col = "eff.allele",
                       other_allele_col = "re.allele")
AOPEP_ALM_gwas <- extract_outcome_data(snps = AOPEP_dat$SNP,
                                       outcomes = 'ebi-a-GCST90000025')
dat <- harmonise_data(AOPEP_dat,AOPEP_ALM_gwas)
res3 <- mr(dat)
res3

het <- mr_heterogeneity(dat)

run_mr_presso(harm_rt,NbDistribution = 1000)

ple <- mr_pleiotropy_test(dat)


singlesnp_res<- mr_singlesnp(dat)
View(singlesnp_res)
singlesnpOR=generate_odds_ratios(singlesnp_res)

sen_res<- mr_leaveoneout(dat)
View(sen_res)

p1 <- mr_scatter_plot(res3, dat)
p1[[1]]

p2 <- mr_forest_plot(singlesnp_res)
p2[[1]]

p3 <- mr_leaveoneout_plot(sen_res)
p3[[1]]

res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(singlesnp_res)
p4[[1]]

expo_rt<- read_exposure_data(
  filename = "Lachnospiraceae.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "eff.allele",
  other_allele_col = "re.allele",
  pval_col = "P",
  samplesize_col = "N")

expo_rt<- expo_rt[expo_rt$pval.exposure < 1e-06,]
expo_rt <- clump_data(expo_rt,clump_kb = 100,clump_r2 = 0.01)  

MiBioGen_Clean_ALM_gwas <- extract_outcome_data(snps = expo_rt$SNP,
                                       outcomes = 'ebi-a-GCST90000025')


dat <- harmonise_data(expo_rt,MiBioGen_Clean_ALM_gwas)
res3 <- mr(dat)
res3

het <- mr_heterogeneity(dat)


ple <- mr_pleiotropy_test(dat)

singlesnp_res<- mr_singlesnp(dat)
singlesnpOR=generate_odds_ratios(singlesnp_res)

sen_res<- mr_leaveoneout(dat)

p1 <- mr_scatter_plot(res3, dat)
p1[[1]]

p2 <- mr_forest_plot(singlesnp_res)
p2[[1]]

p3 <- mr_leaveoneout_plot(sen_res)
p3[[1]]

res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(singlesnp_res)
p4[[1]]
