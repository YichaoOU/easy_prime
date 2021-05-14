====================================================================
Easy-Prime: *Machine learning based pegRNA design*
====================================================================


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   Installation

Ask questions here
^^^^^^^^^^^^^^^^^

https://github.com/YichaoOU/easy_prime

Summary
^^^^^^

PE design involves carefully choosing a standard sgRNA, a RT template that contains the desired edits, a PBS that primes the RT reaction, and a ngRNA that nicks the non-edit strand. Usually thousands of combinations are available for one single disired edit. Therefore, it is overwhelming to select the most likely high-efficient candidate from the huge number of combinations.

Easy-Prime applies a machine learning model (i.e., XGboost) that learned important PE design features from public PE amplicon sequencing data to help researchers selecting the best candidate.


Installation
^^^^^^^^^^^^^^


::

	conda create -n genome_editing -c cheng_lab easy_prime

	source activate genome_editing

	easy_prime -h


For detailed installation with screenshots, see: :doc:`Installation <Installation>`

Input
^^^^




1. vcf input example

VCF headers will be ignored. Only the first 5 columns from the vcf file will be used; they are: chr, pos, name/id, ref, alt.

::

	## comment line, will be ignored
	chr9	110184636	FIG5G_HEK293T_HEK3_6XHIS	G	GCACCATCATCACCATCAT
	chr1	185056772	FIG5E_U2OS_RNF2_1CG	G	C
	chr1	173878832	rs5878	T	C
	chr11	22647331	FIG3C_FANCF_7AC_PE3B	T	G
	chr19	10244324	EDFIG5B_DNMT1_dPAM	G	T



2. fasta input example

To specify reference and alternative allele, you need two fasta sequences; `_ref` is a keyword that will be recognized as the reference allele and `_alt` is a keyword for target mutations.


::


	>test_ref
	AAAAAAAAAAAAAAAAAAAAAAAAAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
	>test_alt
	AAAAAAAAAAAAAAAAAAAAAAAAAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA





Config file
^^^^^^^

Default values are shown in the following yaml files.

::

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




Output
^^^^^^^^

The output folder contains:

- topX_pegRNAs.csv
- rawX_pegRNAs.csv.gz
- X_p_pegRNAs.csv.gz
- summary.csv

The top candidates are provided in `topX_pegRNAs.csv`. This is a rawX format file. 

rawX format
------

X means the input to machine learning models. Here, rawX basically means the file before machine learning featurization. Specifically, rawX contains 11 + 1 columns. The first 5 columns are from the input vcf file: sample_ID, chr, pos, ref, alt, where sample_ID ends with `_candidate_xxx`, this indicates the N-th combination. The next 6 columns are genomic coordinates: type, seq, chr, start, end, strand, where the `type` could be sgRNA, PBS, RTT, or ngRNA. Since for one PE design, it has to have these 4 components, which means that for one unique `sample_ID`, it has 4 rows specifying the sequences for each of them. The 12-th column, which is optional, is the predicted efficiency; in other words, the Y for machine learning.

Both `topX_pegRNAs.csv` and `rawX_pegRNAs.csv.gz` use this format.

X format
------

X format is the numeric representation of rawX. `X_p` format appends the predicted efficiency to the last column of X.

Main results
-----

The main results, which is the top condidates, is provided in `topX_pegRNAs.csv`.

PE design visualization
----------

Users can visualize the predicted combinations using:


::

	easy_prime_vis -f topX_pegRNAs.csv -s /path/to/genome_fasta.fa



This will output pdf files to a result dir. 



Usage
^^^^^^


::


	git clone https://github.com/YichaoOU/easy_prime

	cd easy_prime/test

	easy_prime -h

	easy_prime --version

	## Please update the genome_fasta in config.yaml 

	easy_prime -c config.yaml -f test.vcf

	## Will output results to a folder


DASH application
----------------

Easy-Prime also provides a dash application. 

Please have dash installed before running the dash application.

::

	git clone https://github.com/YichaoOU/easy_prime

	cd easy_prime/dash_app

	python main.py









