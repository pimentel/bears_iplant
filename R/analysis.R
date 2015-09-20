library('sleuth')

# install.pakages('cowplot')
# if it isn't installed
library('cowplot')

base_dir <- '..'

s2c <- read.table(file.path(base_dir, 'metadata', 'sample_info.tsv'),
  header = TRUE, stringsAsFactors = FALSE)

sample_id <- dir(file.path(base_dir,"results"))
sample_id

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id, "kallisto"))
kal_dirs

so <- sleuth_prep(kal_dirs, s2c, ~ condition)
so <- sleuth_fit(so)
so <- sleuth_test(so, which_beta = 'conditionscramble')

models(so)

sleuth_live(so)
