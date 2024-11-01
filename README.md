[![Version][version-shield]][version-url]
[![Python versions][python-shield]][python-url]
[![Platforms][platform-shield]][python-url]
[![DOI](https://zenodo.org/badge/265117246.svg)](https://zenodo.org/badge/latestdoi/265117246)

# Easy-Prime: an optimized prime editor gRNA design tool based on gradient boosting trees

Easy-Prime provides optimized pegRNA and ngRNA combinations for efficient PE design. Visit: http://easy-prime.cc/

# Summary

PE design involves carefully choosing a standard sgRNA, a RT template that contains the desired edits, a PBS that primes the RT reaction, and a ngRNA that nicks the non-edit strand. Usually thousands of combinations are available for one single disired edit. Therefore, it is overwhelming to select the most likely high-efficient candidate from the huge number of combinations.

Easy-Prime applies a machine learning model (i.e., XGboost) that learns important PE design features from multiple published PE data sources to help researchers selecting the best candidate.

# Docker usage

The simplest way to use easy-prime is to use our dockerized version.

```
docker pull liyc1989/easy_prime

docker run -p 80:80 easy_prime

```

ref: https://github.com/YichaoOU/easy_prime_docker

# Installation

The most easiest way to install Easy-Prime is via conda (version >=4.9). 

```

conda create -n genome_editing -c cheng_lab easy_prime

source activate genome_editing

easy_prime -h

easy_prime_vis -h

```

See https://easy-prime.readthedocs.io/en/latest/content/Installation.html for step-by-step installation screenshots.

# Usage

```

## Make sure you have installed Easy-Prime before running the commands below

git clone https://github.com/YichaoOU/easy_prime

cd easy_prime/test

easy_prime -h

easy_prime --version

## Please update the genome_fasta in config.yaml, otherwise an error may occur!

easy_prime -c config.yaml -f test.vcf

## Will output results to a folder

## to output visualization, go to this output folder and run:

## Please make sure you have provided the correct path to genome fasta using -s

easy_prime_vis -f topX_pegRNAs.csv -s ~/Data/Human/hg19/fasta/hg19.fa

```

Easy-Prime also provides a dash application. 

Please have dash installed before running the dash application.

```

## Make sure you have installed Easy-Prime before running the commands below

git clone https://github.com/YichaoOU/easy_prime

cd easy_prime/dash_app

## Please update the genome_fasta in config.yaml, otherwise an error may occur!

## Please also update the genome_fasta variable in line 56 in utils2.py

python application.py

```

![screenshot](screenshot.png)

# Easy-Prime on AWS

Please use this URL: http://easy-prime.cc/

New to the AWS version, we added linker sequence from `pegLIT` prediction. It may take several minutes to get the linker sequence. See the example below.

![linker](easy_prime_AWS_updates.PNG)

# Tutorial

See https://easy-prime.readthedocs.io/en/latest/content/AWS.html for step-by-step tutorial with screenshots.

## Input

1. vcf input example

VCF headers will be ignored. Only the first 5 columns from the vcf file will be used; they are: chr, pos, name/id, ref, alt.

```
## comment line, will be ignored
chr9	110184636	FIG5G_HEK293T_HEK3_6XHIS	G	GCACCATCATCACCATCAT
chr1	185056772	FIG5E_U2OS_RNF2_1CG	G	C
chr1	173878832	rs5878	T	C
chr11	22647331	FIG3C_FANCF_7AC_PE3B	T	G
chr19	10244324	EDFIG5B_DNMT1_dPAM	G	T

```

2. fasta input example

To specify reference and alternative allele, you need two fasta sequences; `_ref` is a keyword that will be recognized as the reference allele and `_alt` is a keyword for target mutations.

```
>rs2251964_ref
GTTACCAAAGCAAATGACATCTTGTGAAAGGGGAGGTCTGAAAAAAAAAAACAAGTGGGTGGGTTTTTTCAAAGTAGGCCACCGGGCCTGAGATGACCAGAATTCAAATTAGGATGACAGTGTAGTAGGGGAAGCAACCAGAATCGGACCT
>rs2251964_alt
GTTACCAAAGCAAATGACATCTTGTGAAAGGGGAGGTCTGAAAAAAAAAAACAAGTGGGTGGGTTTTTTCAAAGTAGGCCACCGGGCCTGAGATAACCAGAATTCAAATTAGGATGACAGTGTAGTAGGGGAAGCAACCAGAATCGGACCT
```

The PrimeDesign format input is only supported in the Easy-Prime web server. The vcf file separated by space is only supported in the Easy-Prime web server.

## Parameters

Genome: only support hg19 for now. 

## Results

The web output contain two parts:


1. pegRNA table

In this result table, each predicted sgRNA/ngRNA/RTT/PBS configuration will be provided in 4 rows, they will have the same variant ID and predicted efficiency.


2. Sequence visualization

By default, the top prediction will be shown automatically. 


# Input

A vcf file containing at least 5 columns. See `test/test.vcf` for examples.


## Searching parameters for PE design

Default values are shown in the following yaml files.

```yaml

genome_fasta: /path/to/genome.fa
scaffold: GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC
debug: 0
n_jobs: 4
min_PBS_length: 8
max_PBS_length: 17
min_RTT_length: 10
max_RTT_length: 25
min_distance_RTT5: 3
max_ngRNA_distance: 100
max_target_to_sgRNA: 10
sgRNA_length: 20
offset: -3
PAM: NGG

```

# Output

The output folder contains:

- topX_pegRNAs.csv
- rawX_pegRNAs.csv.gz
- X_p_pegRNAs.csv.gz
- summary.csv

The top candidates are provided in `topX_pegRNAs.csv`. This is a rawX format file. 

## rawX format

X means the input to machine learning models. Here, rawX basically means the file before machine learning featurization. Specifically, rawX contains 11 + 1 columns. The first 5 columns are from the input vcf file: sample_ID, chr, pos, ref, alt, where sample_ID ends with `_candidate_xxx`, this indicates the N-th combination. The next 6 columns are genomic coordinates: type, seq, chr, start, end, strand, where the `type` could be sgRNA, PBS, RTT, or ngRNA. Since for one PE design, it has to have these 4 components, which means that for one unique `sample_ID`, it has 4 rows specifying the sequences for each of them. The 12-th column, which is optional, is the predicted efficiency; in other words, the Y for machine learning.

Both `topX_pegRNAs.csv` and `rawX_pegRNAs.csv.gz` use this format.

## X format

X format is the numeric representation of rawX. `X_p` format appends the predicted efficiency to the last column of X.

## Main results

The main results, which is the top condidates, is provided in `topX_pegRNAs.csv`.

# PE design visualization

Users can visualize the predicted combinations using:

```bash

easy_prime_vis -f topX_pegRNAs.csv -s /path/to/genome_fasta.fa

```

This will output pdf files to a result dir. 


## FAQ

#### Can't read input file

```
Error: The requested fasta database file (/Users/yli11/Data/hg19.fa) could not be opened. Exiting!
Reading fasta file: test.vcf
{}
no valid sequences in: test.vcf
Exit...
Can't read test.vcf as vcf or fasta. Please check input. Exit...

```

This error, as it outputs, is caused by the `genome_fasta` parameter in the `config.yaml` file. Make sure you have the correct path to the input genome fasta.


[version-shield]: https://img.shields.io/conda/v/cheng_lab/easy_prime.svg
[version-url]: https://anaconda.org/cheng_lab/easy_prime
[python-shield]: https://img.shields.io/pypi/pyversions/easy_prime.svg
[python-url]: https://pypi.python.org/pypi/easy_prime
[platform-shield]: https://anaconda.org/cheng_lab/easy_prime/badges/platforms.svg

