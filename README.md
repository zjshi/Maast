# Catalogue

* [Maast](https://github.com/zjshi/Maast#maast)
* [What Maast does](https://github.com/zjshi/Maast#what-maast-does)
* [How to cite](https://github.com/zjshi/Maast#how-to-cite)
* [Installation](https://github.com/zjshi/Maast#installation)
* [Conda Installation](https://github.com/zjshi/Maast#conda-installation)
* [How to use](https://github.com/zjshi/Maast#how-to-use)
* * [Type SNPs from a set of whole genome assemblies and sequencing reads from beginning to end in one single command line](https://github.com/zjshi/Maast#type-snps-from-a-set-of-whole-genome-assemblies-and-sequencing-reads-from-beginning-to-end-in-one-single-command-line)
* * [Genotype SNPs step by step](https://github.com/zjshi/Maast#genotype-snps-step-by-step)
* * * [Step 1a: Call SNP with a collection of whole genome assemblies](https://github.com/zjshi/Maast#step-1a-call-snp-with-a-collection-of-whole-genome-assemblies)
* * * [Step 1b: Call SNPs from a set of whole genomes without redundancy reduction](https://github.com/zjshi/Maast#step-1b-call-snps-from-a-set-of-whole-genomes-without-redundancy-reduction)
* * * [Step 1c: Call SNPs with customized minimum prevalence and minor allele frequency (MAF) thresholds](https://github.com/zjshi/Maast#step-1c-call-snps-with-customized-minimum-prevalence-and-minor-allele-frequency-maf-thresholds)
* * * [Step 2: Build SNP covering k-mer database](https://github.com/zjshi/Maast#step-2-build-snp-covering-k-mer-database)
* * * [Step 3: Genotype whole genome assemblies, sequencing reads or both](https://github.com/zjshi/Maast#step-3-genotype-whole-genome-assemblies-sequencing-reads-or-both)
* * * [Construct a SNP tree with Maast genotypes (optional)](https://github.com/zjshi/Maast#construct-a-snp-tree-with-maast-genotypes-optional)
* * * [More helper text and arguments](https://github.com/zjshi/Maast#more-helper-text-and-arguments)
* [Example tutorial](https://github.com/zjshi/Maast#example)
* * [Download and decompress test dataset](https://github.com/zjshi/Maast#download-and-decompress-test-dataset)
* * [Genotype SNPs from begin to end in one single command line with the test dataset](https://github.com/zjshi/Maast#genotype-snps-from-begin-to-end-in-one-single-command-line-with-the-test-dataset)
* * [Genotype SNPs step by step with the test dataset](https://github.com/zjshi/Maast#genotype-snps-step-by-step-with-the-test-dataset)
* * * [Step 1: Call SNPs with whole genome assemblies](https://github.com/zjshi/Maast#step-1-call-snps-with-whole-genome-assemblies)
* * * [Step 2: Build SNP covering k-mer database](https://github.com/zjshi/Maast#step-2-build-snp-covering-k-mer-database-1)
* * * [Step 3: Genotype whole genome assemblies, sequencing reads or both](https://github.com/zjshi/Maast#step-3-genotype-whole-genome-assemblies-sequencing-reads-or-both-1)
* * * [Construct a SNP tree with Maast genotypes (optional)](https://github.com/zjshi/Maast#construct-a-snp-tree-with-maast-genotypes-optional-1)

# Maast

Maast for microbial agile accurate SNP typing  

## What Maast does

Recent spikes in available whole-genome sequences have greatly expanded intra-species diversity especially for prevalent species. As the number of genomes per species grows, it becomes computationally challenging to perform whole-genome alignment and call single nucleotide polymorphisms (SNPs). Furthermore, the genomes from some species are highly similar and hence redundant for SNP discovery. These trends are irreversible and worse over time. To address the challenge, we present Maast, a tool for discovering core-genome SNPs and genotyping these SNPs in conspecific genomes, contigs, or unassembled reads. Maast runs orders of magnitude faster than existing tools and uses less RAM because it is free of read alignment and assembly. Maast is also comparably accurate and recovers more core-genome SNPs compared to other the-state-of-art tools.

## How to cite

The publication of Maast is in preparation. Please cite this GitHub repo as alternative for now. 

## Installation

<b>Python requirement</b>
* Python3 (>=3.6.9)

<b>Required Python libraries</b>
* [NumPy] (https://numpy.org/install/) (>=1.19.5)
* [SciPy] (https://scipy.org/install/) (>=1.5.4)
* [Biopython] (https://biopython.org/wiki/Download) (>=1.79)
* [NetworkX] (https://pypi.org/project/networkx/) (>=2.5.1)

Note: the following installation command line might be helpful
`pip install numpy biopython`

<b>Required external programs</b>
* [Mash](https://github.com/marbl/Mash) (>=v2.2)
* [MUMmer4](https://github.com/mummer4/mummer) (>=v4.0.0)

<b>Optional installation</b>
* [FastTreeMP](http://www.microbesonline.org/fasttree/FastTreeMP) (>= v2.1.11) (Optional; only required when tree subcommand is run)  
* [pigz](https://zlib.net/pigz/) (Optional; A parallel implementation of gzip for modern multi-processor, multi-core machines)
* [lbzip2](http://lbzip2.org/) (Optional; A free, multi-threaded compression utility with support for bzip2 compressed file format)
* [lz4](http://www.lz4.org) (Optional; Extremely Fast Compression algorithm)

Note: the optional dependencies are not required for essential features of Maast, but they are recommended to be installed for better performance or additional features.  

First, retrieve a copy of Maast to your local computing environment   

`git clone https://github.com/zjshi/Maast.git`

Change your current working directory into where you put Maast   
`cd /path/to/Maast/`

Type in the command line to compile the source code of Maast   
`make`

Type in the command line to make GT-Pro ready to execute   
`chmod 755 maast`

The main program (`maast`) should be found in the same directory as `/path/to/Maast/`. This location can be added to the system path so that the main program can be accessed from anywhere. Reference through full path is also allowed.

Type in the command line to display help text   

`./maast -h`  

<b>Notes for C++ compiler</b>   

Maast requires a C++ compiler that is compatible with C++ 11 standards to work properly. All the tests have been done and passed with clang-900.0.38, but it should be compatible for GNU C Compiler (newer than 5.4.0). We have not tested Maast with older compilers, but we expect it to run similiarly as long as it compiles successfully.


## Conda Installation

<b>Create a new conda environment</b>   
`conda create -n maast`

<b>Activate the environment just created</b>   
`conda activate maast`

<b>Conda automatic installation with all dependencies</b>   
`conda install -c conda-forge -c bioconda maast`

<b>Quick installation verification</b>   
`maast -h`  

## How to use 

### Type SNPs from a set of whole genome assemblies and sequencing reads from beginning to end in one single command line
 
`maast end_to_end --in-dir /path/to/directory/containing/genomes/reads/or/both --out-dir /path/Maast/output/`

Note:  

Input directory must have a number of whole genome assemblies in FASTA format. 

Maast can automatically identify file types with supported file suffix: whole genome assemblies (.fa, .fsa, .fna and .fasta) and sequencing reads (.fq and .fastq). Files compressed with popular algorithms, including .gz, .lz4 and .bz2, are also supported.

The running of end_to_end subcomand is equavalent to the running of genomes, db and genotype subcommand with default settings in a row.

### Genotype SNPs step by step

#### Step 1a: Call SNP with a collection of whole genome assemblies
`maast genomes --fna-dir /path/to/genomes/ --rep-fna /path/to/rep_genome.fna --out-dir /path/Maast/output/`  

Note:  
By default, Maast first collapsed redundancy in the input genomes and then call common SNPs from a subset of tag genomes. It also automatically identifies a centroid-genome and use it for the representative genome.

Upon a successful run, this step will produce several important files that are required for downstream steps.
* reference.fna (Reference genome that provides genomic coordinate for SNPs)
* core_snps.vcf (SNP catalog)
* tag_paths.list (Selected tag genomes)

#### Step 1b: Call SNPs from a set of whole genomes without redundancy reduction

`maast genomes --fna-dir /path/to/genomes/ --rep-fna /path/to/rep_genome.fna --out-dir /path/Maast/output/ --skip-centroid --keep-redundancy`  

#### Step 1c: Call SNPs with customized minimum prevalence and minor allele frequency (MAF) thresholds

`maast genomes --fna-dir /path/to/genomes/ --rep-fna /path/to/rep_genome.fna --out-dir /path/Maast/output/ --min-prev 0.95 --snp-freq 0.001`  

#### Step 2: Build SNP covering k-mer database

`maast db --ref-genome /path/to/reference.fna --vcf /path/to/core_snps.vcf --msa /path/to/tag_msa.fna --tag-fna-list /path/to/tag_paths.list --fna-dir /path/to/genomes/ --out-dir /path/Maast/output/`

Note:  

Upon a successful run, this step will produce a SNP covering k-mer database that is required for genotyping sequencing reads.
* kmer_db.bin (SNP covering k-mer database)

#### Step 3: Genotype whole genome assemblies, sequencing reads or both

`maast genotype --in-dir /path/to/directory/containing/genomes/reads/or/both --ref-genome /path/to/reference.fna --db /path/to/kmer_db.bin --vcf /path/to/core_snps.vcf --out-dir /path/Maast/output/`

#### Construct a SNP tree with Maast genotypes (optional)

`maast tree --input-list /path/to/Maast/genotypes.input.tsv --out-dir /path/Maast/output/`

#### More helper text and arguments

`maast end_to_end|genomes|db|genotype|tree -h`

## Example

### Download and decompress test dataset 

`wget --content-disposition https://fileshare.czbiohub.org/s/TwGJAsAZ6dQsM49/download`

`tar xzvf 101346.tar.gz`

Note: after running the two command line above, one directory named 101346 can be found in the current directory. In the directory 101346, there are 300 whole genome assemblie in FASTA format (.fna) and 8 gzipped files of WGS sequencing reads in FASTQ format (.fastq.gz).

### Genotype SNPs from begin to end in one single command line with the test dataset

`maast end_to_end --in-dir ./101346 --out-dir ./101346_out`

Note: after running the above command line, one directory name 101346_out can be found in the currently directory, which contains all resulting files and directories. 

The files include
* reference.fna (selected reference genome)
* tag_paths.list (list of selected tag genomes)
* tag_msa.fna (multiple sequence alignment of tag genomes)
* coords.tsv (coordinates of consensus genome)
* core_snps.vcf (called SNPs in VCF format)
* nr_kmer_set.tsv (raw SNP-covering k-mers)
* check_fna_paths.list (a list of genomes used for validating SNP-covering k-mers)
* kmer_prof.tsv (hit profile of SNP-covering k-mers)
* selected_kmers.tsv (validated SNP-covering k-mers)
* kmer_db.bin (optimized database of SNP-covering k-mers)

The directories include
* gt_results (SNP genotyping results)
* temp (tempory directory for hosting )

### Genotype SNPs step by step with the test dataset

#### Step 1: Call SNPs with whole genome assemblies

`maast genomes --fna-dir ./101346 --out-dir ./101346_out`

Note: upon a successful run of the first step, the output files include
* reference.fna (selected reference genome)
* tag_paths.list (list of selected tag genomes)
* tag_msa.fna (multiple sequence alignment of tag genomes)
* coords.tsv (coordinates of consensus genome)
* core_snps.vcf (called SNPs in VCF format)

#### Step 2: Build SNP covering k-mer database

`maast db --ref-genome ./101346_out/reference.fna --vcf ./101346_out/core_snps.vcf --msa ./101346_out/tag_msa.fna --tag-fna-list ./101346_out/tag_paths.list --fna-dir ./101346/ --out-dir ./101346_out/`

Note: all the required input files can be found from the output files of the first step. 

Upon a successful run of the second step, the output files include
* nr_kmer_set.tsv (raw SNP-covering k-mers)
* check_fna_paths.list (a list of genomes used for validating SNP-covering k-mers)
* kmer_prof.tsv (hit profile of SNP-covering k-mers)
* selected_kmers.tsv (validated SNP-covering k-mers)
* kmer_db.bin (optimized database of SNP-covering k-mers)

Among them, kmer_db.bin is the database file that will be used in the next step along with a few other required files from the first step.

#### Step 3: Genotype whole genome assemblies, sequencing reads or both

`maast genotype --in-dir ./101346/ --ref-genome ./101346_out/reference.fna --db ./101346_out/kmer_db.bin --vcf ./101346_out/core_snps.vcf --out-dir ./101346_out/`

Note: Files to genotype should be supplied in a directory with --in-dir. Supported file types including FASTA and FASTQ formats. Input files can be all FASTAs, FASTQs or a mixture of both.

all other required input files could be found from the output files of two previous steps. 

The main output files are the SNP genotypes that can be found in the a directory named "gt_results" in the designated output directory, ./101346_out/ in this case.

It has seven fields as the following:

1. Contig: string type with arbitary length which specifies the contig of a representative genome where a SNP is from
2. Local Pos: up to seven digits which specifies the local position of a SNP on a contig
3. Global Pos: up to seven digits which specifies the global position of a SNP in a species, served as sort of ID
4. Allele 1: single character, A, C, G or T, which specifies allele 1 of a SNP
5. Allele 2: similiar as Ref allele but specifies allele 2 of a SNP
6. Allele 1 Cnt: an integer specifying the count of detected allele 1 in a metagenome
7. Allele 2 Cnt: an integer specifying the count of detected allele 2 in a metagenome

An example of such looks like the following:

| Contig                                     | Local Pos     | Global Pos     | Allele 1       | Allele 2       | Allele 1 Cnt   | Allele 2 Cnt   |
| :---                                       |    :----:     |     :----:     |    :----:      |    :----:      |    :----:      |    :----:      |
| NODE_10_length_179788_cov_11.0000_ID_43085 | 15829         | 349759         | C              | T              | 65             | 0              |
| NODE_10_length_179788_cov_11.0000_ID_43085 | 15863         | 20713          | C              | T              | 62             | 1              |
| NODE_10_length_179788_cov_11.0000_ID_43085 | 15889         | 131457         | C              | A              | 62             | 0              |
| NODE_10_length_179788_cov_11.0000_ID_43085 | 15907         | 4457           | G              | A              | 59             | 0              |
| NODE_10_length_179788_cov_11.0000_ID_43085 | 15910         | 4553           | C              | A              | 59             | 0              |
| NODE_10_length_179788_cov_11.0000_ID_43085 | 15937         | 151893         | C              | T              | 56             | 0              |
| NODE_10_length_179788_cov_11.0000_ID_43085 | 15940         | 101338         | C              | T              | 55             | 0              |
| ...                                        | ...           | ...            | ...            |  ...           |  ...           |  ...           |

#### Construct a SNP tree with Maast genotypes (optional)

This is an optional step that helps take advantage of genotyped SNPs for a quick application - SNP tree building

`paste <(find ./101346_out/gt_results/ -name '*tsv' | sort) <(find ./101346_out/gt_results/ -name '*tsv' | sort | cut -d'/' -f4 | cut -d'.' -f1) > 101346_genotypes.input.tsv`

Note: the step above generates a list of input pairs. Each pair per row contains a path to a genotype result file generated from Maast genotype command and a unique name of the file. The path and name are separated by a tab, like the following
/file/path/1	name1
/file/path/2	name2
/file/path/3    name3
...


The first three rows of 101346_genotypes.input.tsv in this example look like
./101346_out/gt_results/GUT_GENOME000400.fna.tsv	GUT_GENOME000400
./101346_out/gt_results/GUT_GENOME000466.fna.tsv	GUT_GENOME000466
./101346_out/gt_results/GUT_GENOME000688.fna.tsv	GUT_GENOME000688


`maast tree --input-list ./101346_genotypes.input.tsv --out-dir ./101346_out/`

Note: upon the successful completion of this command, the following three output can be found:
* concat_allele.aln.fasta (concatenated allele sequences with genotyped SNPs)
* concat_allele.aln.mat (Pairwise genomic distances between concatenated allele sequences)
* concat_allele.aln.tre (Phylogenetic tree built with concatenated allele sequences)
