install('~/dev/sleuth/')
library('sleuth')

# install.packages('cowplot')
# if it isn't installed
library('cowplot')

base_dir <- '..'

# get the sample to covariate mapping
s2c <- read.table(file.path(base_dir, 'metadata', 'sample_info.tsv'),
  header = TRUE, stringsAsFactors = FALSE)

# get the sample id names
sample_id <- dir(file.path(base_dir, "results", "paired"))
sample_id

# create a vector that points to all the kallisto results
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results",
  "paired", id, "kallisto"))
kal_dirs

# add the kallisto results to the experimental setup
s2c <- mutate(s2c, path = kal_dirs)
s2c

# get the gene names using biomaRt
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# if you want to re-create the analysis off-line, make sure to save the
# transcript to gene mapping below
#
# saveRDS(t2g, '../metadata/t2g.rds')
# t2g <- readRDS('../metadata/t2g.rds')


# create the sleuth object
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
# fit the full model
so <- sleuth_fit(so)
# perform the test
so <- sleuth_wt(so, which_beta = 'conditionscramble')

# see which models and test have been performed
models(so)

# open up the shiny object for interactive analysis
sleuth_live(so)
