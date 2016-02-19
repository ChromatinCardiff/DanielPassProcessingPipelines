library("sleuth")
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
# Source https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html
###############
## Load data ##
###############
base_dir <- "/home/sbi6dap/Projects/NON-ATHAL/BMR"
sample_id <- dir(file.path(base_dir,"kal"))
sample_id
# Make directory links
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "kal", id))
kal_dirs

## Load meta data
s2c <- read.table(file.path(base_dir, "run_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c
# Merge paths (column called paths)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

## Make pairwise inputs
s2c.miR221 <- subset(s2c, condition=="miR-221" | condition=="Neg-ctrl")

#############
## Ok, Go! ##
#############
# Now the “sleuth object” can be constructed. This requires three commands that:
# (1) load the kallisto processed data into the object 
# (2) estimate parameters for the sleuth response error measurement model and
# (3) perform differential analysis (testing). 
# On a laptop the three steps should take about 2 minutes altogether.

## GENE ANNOTATED RUN ##
# Annotate transcripts with gene names from biomart before running
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "btaurus_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
so <- sleuth_fit(so)
# models(so)  ## Look at the models available
so <- sleuth_wt(so, which_beta = 'conditionmiR-221')



## STANDARD RUN ##
so <- sleuth_prep(s2c, ~ condition)
so <- sleuth_fit(so)
so <- sleuth_wt(so, 'conditionmiR-221')


#############################
## Output to shiny or docs ##
#############################
# Make shiny
sleuth_live(so)

# Make tables
results_table <- sleuth_results(so, 'conditionmiR-221')



