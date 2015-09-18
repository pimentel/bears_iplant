N_THREADS = 2

ANNO_PRE = "annotation/human_trans"
ANNO_FA = ANNO_PRE + ".fa.gz"
KAL_IDX = ANNO_PRE + ".kidx"

SAMPLES = ['SRR493366', 'SRR493367', 'SRR493368', 'SRR493369', 'SRR493370', 'SRR493371']

rule all:
    input:
        expand('results/paired/{id}/kallisto/abundance.h5', id = SAMPLES)

rule kallisto_paired:
    input:
        'data/{id}/{id}_1.fastq.gz',
        'data/{id}/{id}_2.fastq.gz',
        KAL_IDX
    output:
        'results/paired/{id}/kallisto',
        'results/paired/{id}/kallisto/abundance.h5'
    threads: N_THREADS
    shell:
        'kallisto quant '
        '-i {KAL_IDX} '
        '-b 30 '
        '--bias '
        '-t {threads} '
        '-o {output[0]} '
        '{input[0]} {input[1]}'
