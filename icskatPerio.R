# Apply interval-censored SKAT test to loose teeth outcome from questionnaire 
# the UK Biobank. 
# Only chr 5 for now - 2277 genes.

library(tidyr)
library(SeqArray)
library(data.table)
library(dplyr)
library(magrittr)
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/ICSKAT_fit_null.R")
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/ICSKAT.R")
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/make_IC_dmat.R")

# input
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
chr <- as.numeric(args[2])

# output name
outRoot <- "icskatPerio_"
# output directory
outputDir <- "/rsrch3/home/biostatistics/rsun3/ukb/analysis/output/periodontitis"

# max gene size?
max_gene_size <- 2000

# load the gene information
setwd("/rsrch3/home/biostatistics/rsun3/github/LungCancerAssoc/Data")
load(file="ensembl_refgene_hg19_20180109.rda")
all_genes <- ensembl_refgene_hg19_20180109 %>%
	filter(Notes == 0 & Chr == chr) %>%
	mutate(Diff = txEnd - txStart) %>% 
	arrange(Diff)

# Pick a section of genes to analyze based on aID argument
genes_per_chunk <- 10
nchunks <- ceiling(table(all_genes$Chr) / genes_per_chunk)
chunk_pts <- cumsum(nchunks)
if (aID <= chunk_pts[1]) {
  #chr <- 1
  start_row <- genes_per_chunk*(aID - 1) + 1
  end_row <- genes_per_chunk*aID
} else {
  #chr <- min(which(chunk_pts >= aID))
  start_row <- (aID - chunk_pts[chr-1] -1) * genes_per_chunk + 1
  end_row <- (aID - chunk_pts[chr-1]) * genes_per_chunk
}
buffer <- 5000

# Cut to that chromosome
#gene_info <- all_genes %>% filter(Chr == chr) %>%
gene_info <- all_genes %>%
  slice(start_row:min(nrow(.), end_row))

# get loose teeth outcome
setwd("/rsrch3/home/biostatistics/rsun3/ukb/phenotypes")
outcomeDat <- fread("periodontitisUKBIC.txt") %>%
  select(eid, Lage, Rage) %>%
  set_colnames(c("ID", "leftDays", "rightDays")) %>%
  filter(!is.na(ID) & !is.na(leftDays) & !is.na(rightDays)) %>%
  mutate(leftTime = leftDays / 365.25) %>%
  mutate(rightTime = ifelse(rightDays == 99999, 999, rightDays / 365.25))

# open the age and sex table
setwd("/rsrch3/home/biostatistics/rsun3/ukb/phenotypes")
ageSexDat <- fread("ukbAgeSex.txt", header=T, data.table=F) %>%
	set_colnames(c("ID", "Age", "Age2", "Age3", "Age4", "Sex")) %>% 
	select(ID, Age, Sex) %>% 
	filter(!is.na(ID) & !is.na(Age) & !is.na(Sex)) %>%
	mutate(AgeSq = Age^2) %>%
	mutate(AgeSqSex = AgeSq * Sex) %>%
	mutate(AgeSex = Age * Sex)

# open the eigenvectors
setwd("/rsrch3/home/biostatistics/rsun3/ukb/phenotypes")
evecsDat <- fread("evecs_n502524_20200702.txt") %>%
	set_colnames(c("ID", paste0("evec", 1:40))) %>%
	select(1:21) %>%
	filter(!is.na(ID) & !is.na(evec1))

# merge covariates
allCov <- merge(ageSexDat, evecsDat, by="ID") %>%
	merge(., outcomeDat, by="ID") %>%
	mutate(ID = paste0(ID, "_", ID)) 

# open GDS file
setwd("/rsrch3/home/biostatistics/rsun3/ukb/gwas_data/imputed_gds")
gdsfile <- seqOpen(paste0("ukb_imp_chr", chr, "_v3_qc.gds"), allow.duplicate = TRUE)

# get the sample IDs
sampleID <- seqGetData(gdsfile, "sample.id")

# get all the SNPs for this run
snpDF <- data.frame(variant_id = seqGetData(gdsfile, "variant.id"), chromosome = seqGetData(gdsfile, "chromosome"),
  position = seqGetData(gdsfile, "position"), allele = seqGetData(gdsfile, "allele"),
	RS = seqGetData(gdsfile, "annotation/id")) %>% 
	filter(chromosome == chr)

# results data frame
resultsDF <- data.frame(gene = gene_info$HGNC_name, start=gene_info$txStart,
                        end=gene_info$txEnd, q=NA, cleanq=NA, rareq=NA, skatp=NA,
                        burdenp=NA, complex=NA, Wskatp=NA, Wburdenp=NA, Wcomplex=NA,
                        Rskatp=NA, Rburdenp=NA, Rcomplex=NA, RWskatp=NA, RWburdenp=NA,
                        RWcomplex=NA)

# loop through nFilters (make sure this is whole number)
for (gene_it in 1:nrow(gene_info)) {

	# pick out the indices
  start_idx <- min(which(snpDF$position >= gene_info$txStart[gene_it] - buffer))
  end_idx <- max(which(snpDF$position <= gene_info$txEnd[gene_it] + buffer))

  # stop if only one SNP
  q <- end_idx - start_idx + 1
  resultsDF$q[gene_it] <- q
  if (q <= 1) (next)

	# reset the filter
	seqResetFilter(gdsfile)
	# set the filter
	snpsToFilter <- snpDF$variant_id[start_idx:end_idx]
	seqSetFilter(gdsfile, variant.id=snpsToFilter)
 
	# get genotypes 
	seqData <- SeqArray::seqGetData(gdsfile, "annotation/format/DS")
	geno <- data.frame(seqData$data) %>%
		mutate(ID = sampleID) %>%
		merge(., allCov, by="ID")

	# drop NA rows - we filtered SNPs with too much missingness so shouldn't drop too many
	cleanG <- geno %>% drop_na() %>%
		select(paste0("X", 1:q)) %>%
		as.matrix(.)

	# some will have MAF 0
  MAFs <- apply(cleanG, 2, mean) / 2
  MAFpos <- which(MAFs > 0)
  cleanG <- cleanG[, MAFpos]
  MAFs <- MAFs[MAFpos]
  if (length(MAFs) < 2 | length(MAFs) > max_gene_size)   {next}
  # weights by beta
  weights <- dbeta(MAFs, 1, 25)
  GW <- cleanG %*% diag(weights)

	# make design matrices
  dmats <- make_IC_dmat(X=as.matrix(geno %>% select("Sex", paste0("evec", 1:10))),
      lt=geno$leftTime, rt=geno$rightTime)

  # fit null model
  obs_ind <- as.numeric(geno$rightTime < 999)
  tpos_ind <- as.numeric(geno$leftTime > 0)
  null_fit <- ICSKAT_fit_null(init_beta=c(rep(0, 14)), lt=geno$leftTime, rt=geno$rightTime,
                              left_dmat=dmats$left_dmat,
                              right_dmat=dmats$right_dmat,
                              obs_ind=obs_ind, tpos_ind=tpos_ind)
  # get pvalue
  skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                        right_dmat=dmats$right_dmat, G=cleanG, lt=geno$leftTime, rt=geno$rightTime,
                        null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)
  # get weighted p-value
  #weighted_skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
  #                                 right_dmat=dmats$right_dmat, G=GW, lt=geno$leftTime, rt=geno$rightTime,
  #                                 null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)
	
	# regular skat results
  resultsDF$skatp[gene_it] <- skat_output$p_SKAT
  resultsDF$burdenp[gene_it] <- skat_output$p_burden
  resultsDF$complex[gene_it] <- skat_output$complex
  # weighted results
  #resultsDF$Wskatp[gene_it] <- weighted_skat_output$p_SKAT
  #resultsDF$Wburdenp[gene_it] <- weighted_skat_output$p_burden
  #resultsDF$Wcomplex[gene_it] <- weighted_skat_output$complex

	# checkpoint
	cat("Done with ", gene_it, " genes \n")
}

# close gds
seqClose(gdsfile)
# write results
setwd(outputDir)
write.table(resultsDF, paste0(outRoot, aID, ".txt"), append=F, quote=F, row.names=F, col.names=T) 





