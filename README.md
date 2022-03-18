# Maast

Maast for microbial agile accurate SNP typing  

## What Maast does

Recent spikes in available whole-genome sequences have greatly expanded intra-species diversity especially for prevalent species. However, the increasing availability also introduced high-level of redundancy which imposed computing burden for core-genome SNP calling. The issue will exacerbate as the trend is irreversible. Maast presented here is a tool free of read alignment and assembly, which detects subspecies structure in conspecific genomes and picks tag and reference genomes for rapid and sensitive core-genome SNP calling in both sequencing reads and whole genomes. Maast runs orders of magnitude faster with less RAM use and recovers more core-genome SNPs comparing to other the-state-of-art tools.

## How to cite

The publication of Maast is in preparation. Please cite this GitHub repo as alternative for now. 

## Installation

<b>Dependencies</b>

* Python3
* numpy
* Biopython

* [Mash](https://github.com/marbl/Mash) (>= v2.2)
* [MUMmer4](https://github.com/mummer4/mummer) (>= v4.0.0)

* [FastTreeMP](http://www.microbesonline.org/fasttree/FastTreeMP) (>= v2.1.11) (Optional; only required when tree subcommand is run)  
* [pigz](https://zlib.net/pigz/) (Optional; A parallel implementation of gzip for modern multi-processor, multi-core machines)
* [lbzip2](http://lbzip2.org/) (Optional; A free, multi-threaded compression utility with support for bzip2 compressed file format)
* [lz4](http://www.lz4.org) (Optional; Extremely Fast Compression algorithm)


First, retrieve a copy of Maast to your local computing environment

`git clone https://github.com/zjshi/Maast.git`

Change your current working directory into where you put Maast
`cd /path/to/Maast/`

Type in the command line to compile the source code of Maast
`make`

Type in the command line to make GT-Pro ready to execute
`chmod 755 maast`

The main program (`maast`) should be found in the same directory as `/path/to/Maast/`. The maast can be added to the system path so that the main program can be accessed from anywhere. Reference through full path is also allowed.

Type in the command line to display help text

`./maast -h`  

<b>Notes for C++ compiler</b>

Maast requires a C++ compiler that is compatible with C++ 11 standards to work properly. All the tests have been done and passed with clang-900.0.38, but it should be compatible for GNU C Compiler (newer than 5.4.0). We have not tested Maast with older compilers, but we expect it to run similiarly as long as it compiles successfully.

## Examples

### Type SNPs from a set of whole genome assemblies and sequencing reads from beginning to end in one single command line
 
`maast end_to_end --in-dir /path/to/directory/containing/genomes/reads/or/both --out-dir /path/Maast/output/`

Note:  

Input directory must have a number of whole genome assemblies in FASTA format. 

Maast can automatically identify file types with supported file suffix: whole genome assemblies (.fa, .fsa, .fna and .fasta) and sequencing reads (.fq and .fastq). Files compressed with popular algorithms, including .gz, .lz4 and .bz2, are also supported.

The running of end_to_end subcomand is equavalent to the running of genomes, db and genotype subcommand with default settings in a row.

### Call SNPs from a set of whole genomes

#### Default SNP calling
`maast genomes --fna-dir /path/to/genomes/ --rep-fna /path/to/rep_genome.fna --out-dir /path/Maast/output/`  

Note:  
By default, Maast first collapsed redundancy in the input genomes and then call common SNPs from a subset of tag genomes. It also automatically identifies a centroid-genome and use it for the representative genome.

Upon the successful running, this step will produce several important files that are required for downstream steps.
* reference.fna (Reference genome that provides genomic coordinate for SNPs)
* core_snps.vcf (SNP catalog)
* tag_paths.list (Selected tag genomes)

#### Call SNPs from a set of whole genomes without redundancy reduction

`maast genomes --fna-dir /path/to/genomes/ --rep-fna /path/to/rep_genome.fna --out-dir /path/Maast/output/ --skip-centroid --keep-redundancy`  

#### Call SNPs with customized Min. prevalence and MAF thresholds

`maast genomes --fna-dir /path/to/genomes/ --rep-fna /path/to/rep_genome.fna --out-dir /path/Maast/output/ --min-prev 0.95 --snp-freq 0.001`  

### Build SNP covering k-mer database

`maast db --ref-genome /path/to/reference.fna --vcf /path/to/core_snps.vcf --msa /path/to/tag_msa.fna --tag-fna-list /path/to/tag_paths.list --fna-dir /path/to/genomes/ --out-dir /path/Maast/output/`

Note:  

Upon the successful running, this step will produce a SNP covering k-mer database that is required for genotyping sequencing reads.
* kmer_db.bin (SNP covering k-mer database)

### Genotype whole genome assemblies, sequencing reads or both

`maast genotype --in-dir /path/to/directory/containing/genomes/reads/or/both --ref-genome /path/to/reference.fna --db /path/to/kmer_db.bin --vcf /path/to/core_snps.vcf --out-dir /path/Maast/output/`

### Construct a SNP tree with Maast genotypes

`maast tree --input-list /path/to/Maast/genotypes.input.tsv --out-dir /path/Maast/output/`
