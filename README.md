# kallisto and sleuth walkthrough

This is a walk-through for [kallisto](http://pachterlab.github.io/kallisto) and
[sleuth](http://pachterlab.github.io/sleuth). 

### Preliminaries --- THE NEXT FEW CODES NEED TO BE DONE ON YOUR INFRASTRUCTURE ACCOUNT

## data organization

In principle, you can organize your data however you please. However, I strongly
recommend following a format such as the one we present here.

Login to your infrastructure account.
First, download the walk-through along some of the auxiliary materials:

```
mkdir -p ~/analysis
cd ~/analysis
wget -O bears_iplant.zip https://github.com/pimentel/bears_iplant/archive/master.zip
unzip bears_iplant.zip
mv bears_iplant-master bears_iplant
cd bears_iplant
```

This will put the walkthrough at `~/analysis/bears_iplant` and place you in that
directory.

Next, you will need to download the data, this takes about 10-20 minutes:

```
wget http://de.iplantcollaborative.org/dl/d/777051B5-1339-457B-A776-5C57B99D7EC0/cuffdiff2_data_res.tar.gz
tar -xvf cuffdiff2_data_res.tar.gz
```

# kallisto overview

The kallisto pipeline is quite simple. There are basically two steps:

1. build an index (once per organism or annotation)
2. quantify against the index (once per each experimental sample)

You may have noticed that there is no alignment step. This is because part of
the kallisto algorithm performs a very fast "alignment" which we call the
pseudoalignment. This simplifies things from the user's point of view since
there are no extra intermediate files.

## index

Before you can quantify with kallisto, you must create an
index from an annotation file. For RNA-Seq, an annotation file is the set of
cDNA transcripts in FASTA format (the "transcriptome"). We provide some common
transcriptomes [here](http://bio.math.berkeley.edu/kallisto/transcriptomes/).
Ensembl has a collection of many annotations, including the one that we will use today.
If you are not working with a model organism, one option is to do a de novo
assembly.

An `index` operation creates the target [de Bruijn](http://thegenomefactory.blogspot.com/2013/08/how-to-pronounce-de-bruijn.html?m=1) graph structure along with a
few other data structures that we use in kallisto for fast pseudoalignment. It
simply takes a FASTA file and outputs an index in a binary format that we
designed for kallisto. There is no default extension for the index, though I
like to use `.kidx`.

Let's start off by downloading the annotation.  This download takes several minutes:

```bash
mkdir ~/analysis/bears_iplant/annotation
curl -o ~/analysis/bears_iplant/annotation/human_trans.fa.gz  http://bio.math.berkeley.edu/kallisto/transcriptomes/Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz
```

Next, let's build an index from the provided Ensembl human transcriptome:

```bash
kallisto index -i ~/analysis/bears_iplant/annotation/human_trans.kidx ~/analysis/bears_iplant/annotation/human_trans.fa.gz
```

This creates an index at `~/analysis/bears_iplant/annotation/human_trans.kidx` from the annotation
`~/analysis/bears_iplant/annotation/human_trans.fa.gz` with the default k-mer size (k = 31). Note, that
if you have very short reads (e.g. 35bp), you should change `k` to something
smaller (e.g. `-k 21`).

## quantification

Once you index an annotation, you can quantify any number of samples against it.
The quantification step includes pseudoalignment as well as running the EM
algorithm to do estimation of transcript level abundances.

The parameters are pretty minimal. You must supply an index, an output location,
and a set of reads. In the case of single-end data, you must specify that they
are single-end, as well as a fragment length distribution. In the case of
paired-end data, we can infer the fragment length distribution from the data.

There is also one other important parameter: the number of bootstrap iterations.
By default, kallisto runs zero bootstrap iterations. If you do not plan to run
sleuth for differential expression analysis, this is okay. But if you plan to
run sleuth, you must provide a nonzero number of bootstraps. In general, this
number should be at least 30.

A basic quantification example for running sleuth afterwards looks like this
(though, the command below won't work because this data doesn't actually exist):

```bash
kallisto quant -i ~/analysis/bears_iplant/annotation/some_index.kidx -b 30 -t 2 -o ~/analysis/bears_iplant/results/sample_id \
  ~/analysis/bears_iplant/data/sample_id/sample_id_1.fastq.gz ~/analysis/bears_iplant/data/sample_id/sample_id_2.fastq.gz
```

where:

- `-i ~/bears_iplant/analysis/annotation/some_index.kidx` indicates to use this `kallisto index`
- `-b 30` indicates to generate 30 bootstrap samples
- `-t 2` indicates to use 2 threads
- `sample_id` is some sample identifier that is unique to each sample
- `-o ~/bears_iplant/analysis/results/sample_id` indicates to place output in directory `~/analysis/results/sample_id` 
-       (e.g. `~/bears_iplant/analysis/results/sample_id/abundance.h5`)
- the final two arguments are the 'left' and 'right' reads, respectively

### running kallisto

Run the following commands to start kallisto.  Note how fast it is compared to other counting programs.

```
mkdir ~/bears_iplant/analysis/results/paired/SRR493366
kallisto quant -i ~/bears_iplant/analysis/annotation/human_trans.kidx -b 30 --bias -t 2 \
  -o ~/bears_iplant/analysis/results/paired/SRR493366/kallisto \
  ~/bears_iplant/analysis/data/SRR493366/SRR493366_1.fastq.gz ~/analysis/data/SRR493366/SRR493366_2.fastq.gz
```

### output

After quantification, you will get a number of files in the output directory.

- `run_info.json` -  some high-level information about the run, including the command and versions of kallisto used to generate the output
- `abundance.tsv` - a plain text file with transcript level abundance estimates. This file can be read into R or any other statistical language easily (e.g. `read.table('abundance.tsv')`)
- `abundance.h5` - a HDF5 file containing all of the quantification information including bootstraps and other auxiliary information from the run. This file is read by sleuth

# sleuth

In the wild, a sleuth is a pack of bears. In the context of RNA-Seq, a sleuth is
a pack of kallisto. The job of sleuth is to perform aggregate analysis of many
related samples at once. Currently, the main role of sleuth is to do
differential expression analysis at the transcript level. Sleuth differs from most other differential
expression tools by modeling the technical variance due to the transcript abundance
estimation along with the biological variability between samples. Most methods
only model the biological variability.

To explain how to use __sleuth__ we provide an example based on the data in the "Cuffdiff2 paper":

* [Differential analysis of gene regulation at transcript resolution with RNA-seq](http://www.nature.com/nbt/journal/v31/n1/full/nbt.2450.html)	by Cole Trapnell, David G Henderickson, Martin Savageau, Loyal Goff, John L Rinn and Lior Pachter, Nature Biotechnology __31__, 46--53 (2013).

The human fibroblast RNA-Seq data for the paper is available on GEO at accession [GSE37704](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37704). The samples to be analyzed are the six samples LFB_scramble_hiseq_repA, LFB_scramble_hiseq_repB, LFB_scramble_hiseq_repC, LFB_HOXA1KD_hiseq_repA, LFB_HOXA1KD_hiseq_repA, and LFB_HOXA1KD_hiseq_repC. These are three biological replicates in each of two conditions (scramble and HoxA1 knockdown) that will be compared with __sleuth__.

## Now we need to download our data onto our desktop to continue

Make a directory on your desktop called "bears_iplant"
Copy your output using WinSCP or Cyberduck into this directory.

Start up RStudio and navigate to `R` subdirectory in the directory we've been
working in.

```r
setwd('~/Desktop/bears_iplant/R')
```

First, let's install `sleuth` and `biomaRt`, a tool that we will use later for
getting the gene names:

```r
source('http://bioconductor.org/biocLite.R')
biocLite("devtools")
biocLite("pachterlab/sleuth")
biocLite("biomaRt")
```

Open a new file if you would like to type the commands and add them as we go
along, or you can simply open `analysis.R` and follow along. You can execute a
line in Rstudio using `ctrl + enter`.

Next, load sleuth:

```r
library('sleuth')
```

Though not required, I also suggest loading a package called `cowplot` which makes the `ggplot`
default much more aesthetically pleasing:

```r
# install.packages('cowplot')
# if it isn't installed
library('cowplot')
```

Let's also set the base working directory:

```r
base_dir <- '..'
```

From here on, all the commands will be in R unless otherwise specified.

## preparing your data

The main requirements of sleuth are:

- sample to covariates mapping
- output from kallisto

### sample to covariate mapping

The sample to covariate mapping is simply a table that describes the experiment. possibly the most
challenging part of sleuth is organizing your data into a way that it can be
read and easily. The only real requirement is that there is at least one column
labeled 'sample'. The remaining columns are important as they describe your
experiment, but the column names can be pretty much any valid string.

Our data is pretty simple in that there is only one covariate here: the
experimental condition.

This is what the file looks like (from the terminal --- NOT INSIDE R):

```bash
cat metadata/sample_info.tsv
sample  condition
SRR493366       scramble
SRR493367       scramble
SRR493368       scramble
SRR493369       HOXA1KD
SRR493370       HOXA1KD
SRR493371       HOXA1KD
```

Let's load this file in R:

```r
s2c <- read.table(file.path(base_dir, 'metadata', 'sample_info.tsv'),
  header = TRUE, stringsAsFactors = FALSE)
```

### locating kallisto output

Next, we have to tell `sleuth` where all the kallisto output is. If you've
organized your data like we did in the snake file, this is quite easy. Sleuth is
simply expecting a character array pointing to all the directories.

Next get the list of sample IDs with

```r
sample_id <- dir(file.path(base_dir, "results", "paired"))
```

The result can be displayed by typing
```r
sample_id
## [1] "SRR493366" "SRR493367" "SRR493368" "SRR493369" "SRR493370" "SRR493371"
```

In the box above, lines beginning with ## show the output of the command (in
what follows we include the output that should appear with each command).

A list of paths to the kallisto results indexed by the sample IDs is collated with

```r
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results",
  "paired", id, "kallisto"))
```

Finally, we have to add these directories to the sample to covariates data.frameto:

```r
s2c <- mutate(s2c, path = kal_dirs)
```

Be sure to look at the data frame to ensure consistency between the sample ids as well as the path that was provided:

```r
s2c
```

### getting gene names

Since the gene names are not automatically in the annotation, we need to get
them from elsewhere. `biomaRt` provide a good way of getting the mappings for
Ensembl annotations.

```r
# get the gene names using biomaRt
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
```

### fitting the model

Now the “sleuth object” can be constructed. This requires three commands that
(1) load the kallisto processed data into the object (2) estimate parameters for
the sleuth response error measurement model and (3) perform differential analyis
(testing). On a laptop the three steps should take about 2 minutes altogether.

First type

```r
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
## reading in kallisto results
## ......
## normalizing est_counts
## 42193 targets passed the filter
## normalizing tpm
## normalizing bootstrap samples
```

then

```r
so <- sleuth_fit(so)
## summarizing bootstraps
## fitting measurement error models
## shrinkage estimation
## computing variance of betas
```

and finally

```r
so <- sleuth_wt(so, which_beta = 'conditionscramble')
```

In general, one can see the possible tests that could be performed using the
`which_beta` parameter in `sleuth_wt` and examining the coefficients:

```r
models(so)
## [  full  ]
## formula:  ~condition
## coefficients:
##  (Intercept)
##      conditionscramble
```

### interactive analysis

Sleuth provides many different ways visualize your data. Most visualizations are
prefixed by `plot_`. While this is true, we think the best way to analyze your
data is using sleuth live. Sleuth live gives you an interactive visualization
along with all the differential expression analysis together. You can execute
sleuth live with:

```r
sleuth_live(so)
```

Let's chat about what sort of things to look out for.
