# kallisto and sleuth walkthrough

# Preliminaries

There are many different ways to run these commands in a pipeline. I personally prefer to use
[`snakemake`](https://bitbucket.org/johanneskoester/snakemake/wiki/Home).
A more complete analysis of using this data using `snakemake` can be [found
here](https://github.com/pachterlab/bears_analyses) (along with other
examples).

Provided is a very simplified `Snakefile` that will help us iterate
over the several different `fastq` files.

I will not do more complicated analyses like get the sample ID numbers and pull
down the particular data sets.

## data organization

In principle you can organize your data however you please. However I strongly
suggest you organize it in the way of it we expect in the `Snakefile`. Please
download your data in the following organization:

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

Next, you will need to download the data.

```
TODO: download link once it's done
```

All that the following scripts expect is that you are working in this directory and the raw data conforms to this format (\*.fastq.gz files).

# kallisto overview

The kallisto pipeline is quite simple. There are basically two steps:

1. build an index (once per organism)
2. quantify against the index (once per each experimental sample)

You may have noticed that there is no alignment step. This is because part of
the kallisto algorithm performs a very fast "alignment" which we call the
pseudoalignment. This simplifies things from the user's point of view since
there are no extra intermediate files.

So the basic pipeline look something like this:

** TODO:  insert figure **

## index

Before you can use kallisto for any type of quantification, you must create an
index from an annotation file. For RNA-Seq, an annotation file is the set of
cDNA transcripts in a FASTA format (the "transcriptome").  For model organisms
these can commonly be found online. We provide some common transcriptomes
[here](http://bio.math.berkeley.edu/kallisto/transcriptomes/). Ensembl has a
collection of many of them including one that we will use today. If you do not
have a model organism one option might be to do a de novo assembly. I won't
cover that today but others might.

An index operation creates the target de Bruijn graph structure along with a few
other data structures that we use in kallisto for fast pseudoalignment. It
simply takes a fasta file and outputs a index in a binary format that we specify
specifically for kallisto. There is no default extension for the index, though I
like to use `.kidx`.

Let's start off by downloading the annotation (note: this URL might change in the future):

```{sh}
mkdir annotation
curl -o annotation/human_trans.fa.gz http://bio.math.berkeley.edu/kallisto/transcriptomes/Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz
```

Next, let's build an index from the provided ensemble human transcriptome:

```{sh}
kallisto index -i annotation/human_trans.kidx annotation/human_trans.fa.gz
```

## quantification

Once you index an annotation, you can quantify any number of samples against it.
The quantification step includes pseudoalignment as well as running the EM
algorithm to do estimation of the transcript level abundances.

The parameters are pretty minimal. You must supply an index, an output location,
and a set of reads. In the case of single and data, you must specify that they
are single and as well as a fragment length distribution. In the case of paired
and data we can simply infer the fragment length distribution.

There is also one other important parameter: the number of bootstrap iterations.
By default, kallisto runs zero bootstrap iterations. If you do not plan to run
sleuth, this is okay. But if you plan to run sleuth for differential expression
analysis, you must provide a nonzero number of bootstraps. In general, this
number should be at least 30.

A basic quantification example for running sleuth afterwards looks like this
(though, the command below won't work because this data doesn't actually exist):

```{sh}
kallisto quant -i annotation/some_index.kidx -b 30 -t 2 -o results/sample_id \
  data/sample_id/sample_id_1.fastq.gz data/sample_id/sample_id_2.fastq.gz
```

- `-i annotation/some_index.kidx` indicates to use this `kallisto index`
- `-b 30` indicates to generate 30 bootstrap samples
- `-t 2` indicates use 2 threads
- `sample_id` is some sample identifier that is unique to each sample
- `-o results/sample_id` indicates to place output starting at `results/sample_id` (e.g. `results/sample_id/abundance.h5`)
- the final two arguments are the 'left' and 'right' reads, respectively

### executing snakemake

TODO: explain why snakemake is important

Open up `Snakefile` and have a let's talk about what it's doing.

Since we have a `Snakefile`, we can simply run the `snakemake` command:

```{sh}
snakemake -p -j 2
```


- `-p` prints out the actual command that will execute
- `-j 2` specifies that there are two available processors

### output

After quantification, you will get a number of files in the output directory.

- `run_info.json` -  some high-level information about the run including the command and versions of kallisto used to generate the output
- `abundance.tsv` - a plain text file with transcript level abundance estimates.This file can be read into our or any other statistical language easily (e.g. `read.tabe('abundance.tsv')`)
- `abundance.h5` - a HDF5 file containing all of the quantification information including bootstraps and other auxiliary information from the run. This file is imported into sleuth

# sleuth

In the wild, a sleuth is a pack of bears. In the context of RNA-Seq, a sleuth is
a pack of kallisto. The job of sleuth is to do aggregate analysis of many
related samples at once. Currently, the main role of sleuth is to do differential
transcription expression analysis. Sleuth differs from most other differential
expression tools by modeling the technical error due to the transcript abundance
estimation along with the biological variability between samples. Most methods
only model the biological variability.

## preliminaries

start up RStudio and  load up sleuth:

```{r}
library('sleuth')
```

Though not required, I also suggest loading a package called `cowplot` which makes the `ggplot`
default much more aesthetically pleasing:

```{r}
library('cowplot')
```

From here on, all the commands will be in our unless otherwise specified.

## preparing your data

The main requirements of sleuth are:

- sample to covariates mapping
- output from kallisto


### sample to covariate mapping

This is simply a table that describes the experiment. possibly the most
challenging part of sleuth is organizing your data into a way that it can be
read and easily. The only real requirement is that there is at least one column
labeled 'sample'. The remaining columns are important as they describe your
experiment but the column names can realistically be pretty much anything.

Our data is pretty simple in that there is only one covariate here: the
experimental condition.

This is what the file looks like (from the terminal):

```{sh}
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

```{r}
s2c <- read.table('metadata/hiseq_info.tsv', header = TRUE,
  stringsAsFactors = FALSE)
```

### locating kallisto output

Next, we have to tell `sleuth` where all the kallisto output is. If you've
organized your data like we did in the snake file, this is quite easy. Sleuth is
simply expecting a character array pointing to all the directories.

Next get the list of sample IDs with

```{r}
base_dir <- '.'
sample_id <- dir(file.path(base_dir,"results"))
```

The result can be displayed by typing
```{r}
sample_id
## [1] "SRR493366" "SRR493367" "SRR493368" "SRR493369" "SRR493370" "SRR493371"
```

In the box above, lines beginning with ## show the output of the command (in
what follows we include the output that should appear with each command).

A list of paths to the kallisto results indexed by the sample IDs is collated with

```{r}
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id, "kallisto"))
kal_dirs
```
TODO: put kal_dirs output

### fitting the model

Now the “sleuth object” can be constructed. This requires three commands that
(1) load the kallisto processed data into the object (2) estimate parameters for
the sleuth response error measurement model and (3) perform differential analyis
(testing). On a laptop the three steps should take about 2 minutes altogether.

First type

```{r}
so <- sleuth_prep(kal_dirs, s2c, ~ condition)
## reading in kallisto results
## ......
## normalizing est_counts
## 42193 targets passed the filter
## normalizing tpm
## normalizing bootstrap samples
```

then

```{r}
so <- sleuth_fit(so)
## summarizing bootstraps
## fitting measurement error models
## shrinkage estimation
## computing variance of betas
```

and finally

```{r}
so <- sleuth_test(so, which_beta = 'conditionscramble')
```

In general, one can see the possible tests that could be performed using the
`which_beta` parameter in `sleuth_test` and examining the coefficients:

```{r}
models(so)
## [  full  ]
## formula:  ~condition
## coefficients:
##  (Intercept)
##      conditionscramble
## tests:
##  conditionscramble
```

### interactive analysis

Sleuth provides many different ways visualize your data. Most visualizations are
prefixed by `plot_`. While this is true, we think the best way to analyze your
data is using sleuth live. Sleuth live gives you an interactive visualization
along with all the differential expression analysis together. You can execute
sleuth live with:

```{r}
sleuth_live(so)
```

Let's chat about what sort of things to look out for.
