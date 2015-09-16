# kallisto and sleuth walkthrough

# Preliminaries

There are many different ways to run these commands in a pipeline. I personally prefer to use
[`snakemake`](https://bitbucket.org/johanneskoester/snakemake/wiki/Home).
A more complete analysis of using this data using `snakemake` can be [found
here](https://github.com/pachterlab/bears_analyses) (along with other
examples). I am going to use a very simplified `Snakefile` that will help us iterate
over the several different fastq files.

I will not do more complicated analyses like get the sample ID numbers and pull
down the particular data sets.

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

Let's build an index from the provided ensemble human transcriptome:

```{sh}
kallisto index -i
```


## Quantification
