# Unorganized notes on how we got the amplicon sequence and editing efficiency.

## Summary

Observation:
 - a lot amplicon sequences don't contain gRNAs
 - reads can't mapped to amplicon sequences
 
1. To identify the correct amplicon sequences for each fastq file, we can first map fastq to the human genome. Then we will know what the amplicon sequence look like.

2. get the correct pegRNA for each fastq file, there are only limited number of pegRNAs provided in the paper, so we can just directly map pegRNA to the amplicon sequence.

3. to get PBS RTS sequnece, this one has a one to one relationship, so it is not hard for this part.

4. to get nick gRNA, again, one to one

5. to get target mutation (we need target mutation to convert it to vcf file.), we can just map PBS RTS sequence to the amplicon sequence. 
        - to get the editing frequency, we just need the relative position
        
        - to get the vcf file, we need the absolution position, using amplicons coordinates


## Data Summary

The following data is removed:

- PE1
- ABEmax
- HDR
- etc

RNA-seq data is paired end, will also be filtered out.

Hard to separater PE2 and PE3 using just the labels. So we are going to analyze them together.


### confused labelling examples

EDFIG   =  FIGED


HBB_EDFIG4G_2TTOA this figure not found in the data         this is PE1

PRNP_EDFIG8 this figure not found in the data      this is RNA-seq

HEXA_EDFIG8 this figure not found in the data        this is RNA-seq

let just used FigED10, in the SRA data summary, this has been divided to abvcdef

HBB_EDFIG10_4ATOT this figure not found in the data

PRNP_EDFIG10_6GTOT this figure not found in the data

﻿sgRNA and pegRNA sequences can be found in Supplementary Table 3 under the Figure 5 (should be EDfig6) heading



## get targeted mutation

1. extract 3' extend sequence, to fasta, to fastq bwa mem mapping, to vcf, automatically infer targeted mutation, some could be 293T SNPs.

## method to calculate editing frequency


no need to get frequency right now, we can first parse amplicon.fasta and compare 3' extsion sequence and get mutation pos


Given RTS and PBS sequence together, in David Liu's paper, they call it 3' extension sequencing.

check where is your PBS_RTS sequence map in the amplicon sequence (Nucleotide_frequency_table.txt), find the differenct positions

should be only 1 difference. 

Check if it is the same specified in the file name

if more than 1, check which one is the same

if not any similar, print error

if yes, get the editing frequency from (Nucleotide_percentage_table.txt)

if count less than 10, print error

py2

>>> from skbio.alignment import local_pairwise_align_ssw
>>> alignment = local_pairwise_align_ssw(
...                 "ACTAAGGCTCTCTACCCCTCTCAGAGA",
...                 "ACTAAGGCTCCTAACCCCCTTTTCTCAGA"
...             )
>>> print alignment
>query
ACTAAGGCTCTC-TACCC----CTCTCAGA
>target
ACTAAGGCTC-CTAACCCCCTTTTCTCAGA


>>> df = pd.read_csv("Quantification_window_nucleotide_percentage_table.txt",sep="\t",index_col=0)
>>> df
          C         G         T       G.1         A       T.1       G.2  \
A  0.000092  0.000246  0.049252  0.000092  0.998138  0.000092  0.000185
C  0.998846  0.000062  0.000231  0.000108  0.000215  0.000077  0.000046
G  0.000323  0.999523  0.001062  0.999507  0.000569  0.000431  0.999431
T  0.000739  0.000169  0.949455  0.000292  0.001077  0.999400  0.000339
N  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
-  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000

        G.3       C.1       A.1    ...          C.6       C.7       T.5  \
A  0.000185  0.000154  0.998661    ...     0.001431  0.000893  0.000123
C  0.000139  0.999338  0.000354    ...     0.998261  0.998492  0.000308
G  0.999523  0.000277  0.000847    ...     0.000031  0.000031  0.000215
T  0.000154  0.000231  0.000139    ...     0.000277  0.000585  0.999354
N  0.000000  0.000000  0.000000    ...     0.000000  0.000000  0.000000
-  0.000000  0.000000  0.000000    ...     0.000000  0.000000  0.000000

        C.8       C.9       A.8      G.11       A.9      G.12      G.13
A  0.000292  0.000539  0.997861  0.000185  0.997645  0.000631  0.000754
C  0.999277  0.999123  0.000262  0.000077  0.000954  0.000031  0.000062
G  0.000000  0.000062  0.000754  0.999631  0.000939  0.999030  0.999000
T  0.000431  0.000277  0.001108  0.000077  0.000446  0.000277  0.000185
N  0.000000  0.000000  0.000015  0.000031  0.000015  0.000031  0.000000
-  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000

[6 rows x 40 columns]
>>> seq = [x.split(".")[0] for x in df.columns]
>>> seq
['C', 'G', 'T', 'G', 'A', 'T', 'G', 'G', 'C', 'A', 'G', 'A', 'G', 'G', 'A', 'A', 'A', 'G', 'G', 'A', 'A', 'G', 'C', 'C', 'C', 'T', 'G', 'C', 'T', 'T', 'C', 'C', 'T', 'C', 'C', 'A', 'G', 'A', 'G', 'G']
>>> seq = "".join(seq)
>>> seq
'CGTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGG'
>>>           TCTGCCATCTCGTGCTCA
KeyboardInterrupt
>>> TGAGCACGAGATGGCAGA




### correct crispresso output

/research/rgs01/project_space/chenggrp/blood_regulome/chenggrp/gRNAs/David_liu_pegRNA/analysis/samples_with_no_alginment/given_amplicon_rerun_yli11_2020-04-11

## file

pegRNA.tsv is the original file copied from supp file

fixed problems for these two gRNAs
HEK4_2b_18 GGCACTGCGGCTGGAGGTGG TCTCCGCTTTAACCCCAACCTCCAGCCGC
HEK4_2b_19 GGCACTGCGGCTGGAGGTGG GCTCCGCTTTAACCCCAACCTCCAGCCGC




## debug crisprEsso error

Note that the errors were not caused by errors in the pipeline, but simply caused by the confusing info in the paper. So we need to find the correct pairing between data and amplicon sequence, pegRNA, PBS RTS, and nick-gRNAs.

1. Most of the data should have alignments. Let's find out those wrong data.

	error data have been identified and reanalyzed

2. now, given amplicon sequence, find target mutation

assume PBS RTS sequence is correct, match it to amplicon sequence, find differences.

all_amplicon.fa
all_RTS.fa





4-11 - 4-12 work in folder find_target_mutation



## construct input for pegRNA

replace CATGGTGCACCTGACTCCTG to CATGGTGCATCTGACTCCTG

(TAACGGCAGACTTCTCCAC)
(TAACGGCAGACTTCTCCTC)


replace some primers


the following data has limited number of reads

/home/yli11/dirs/David_liu_pegRNA/sra_download_yli11_2020-03-29/SRR10284412.fastq.gz	EDFIG6D_HEK3_PE3_A8_REP2	ATGTGGGCTGCCTAGAAAGG	CCCAGCCAAACTTGTCAACC	GCCCAGACTGAGCACGTGA

alginment problem is not caused by wrong primer
/home/yli11/dirs/David_liu_pegRNA/sra_download_yli11_2020-03-29/SRR10284415.fastq.gz	EDFIG6D_HEK3_PE3_A5_REP2	ATGTGGGCTGCCTAGAAAGG	CCCAGCCAAACTTGTCAACC	GCCCAGACTGAGCACGTGA

alginment problem is not caused by wrong primer
/home/yli11/dirs/David_liu_pegRNA/sra_download_yli11_2020-03-29/SRR10284471.fastq.gz	EDFIG6A_HEK3_PE3_C4_REP2	ATGTGGGCTGCCTAGAAAGG	CCCAGCCAAACTTGTCAACC	GCCCAGACTGAGCACGTGA

[yli11@nodecn204 analysis]$ run_lsf.py -f update_input.list --user_lsf human_SNV_analysis.lsf -j updated_human_SNV
2020-04-08 14:56:28,731 - INFO - main - The job id is: updated_human_SNV
2020-04-08 14:56:28,961 - INFO - submit_pipeline_jobs - BaseE has been submitted; JobID: 100726252
2020-04-08 14:56:29,101 - INFO - submit_pipeline_jobs - email has been submitted; JobID: 100726254



2. now, given amplicon sequence, find target mutation

assume PBS RTS sequence is correct, match it to amplicon sequence, find differences.

all_amplicon.fa
all_RTS.fa


note that our RTS assignment may not be correct.


## note

I could have included some off-target fastq files here in my study, now I have to remove them, the way to find them is through primer

CRISPResso -r1 $rep2 -a $original -g $guide -o rep2_to_match_figure



## steps


bash run.sh

python get_mutation.py 

manually check results



##

manually fixed FIG4H_HEK3_1CTTINS_2GC problem rerun

also changed david_liu_pegRNA_final_table_everything.tsv corresponding expected amplicon sequnece and target mutation correction VCF pos ref alt columns

removed EMX1 67, 56, 567 combinations.

raw_table_XY_for_training.csv contains all the data for traning



