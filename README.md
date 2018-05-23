# QIIME2 18S rRNA tag-sequencing pipeline - Caron Lab

We're slowly updating to QIIME2 for analyzing our tag-sequencing data. For our QIIME1 pipeline see another [github repo](https://github.com/shu251/https://github.com/shu251/V4_tagsequencing_18Sdiversity_q1) repo. Upstream from this pipeline is our [DNA/RNA extraction protocol](https://www.protocols.io/view/rna-and-optional-dna-extraction-from-environmental-hk3b4yn) and library prep for amplifying and sequencing the V4 region of the 18S rRNA gene, [here](https://www.protocols.io/view/18s-v4-tag-sequencing-pcr-amplification-and-librar-hdmb246).

This protocol is specific to analyzing microbial eukaryotic diversity by way of 18S rRNA gene tag sequencing. Here, we use [V4 Stoeck et al. 2010 primers](https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-294X.2009.04480.x). The product is an approximately 400 bp region. [more](https://onlinelibrary.wiley.com/doi/full/10.1111/jeu.12217)

* Contributors to this pipeline: Sarah Hu & Zhenfeng Liu.
* Alternate protocol for [subsampled open-reference OTU clustering](http://qiime.org/tutorials/open_reference_illumina_processing.html) using QIIME1 can be found at step 10 at another [github repo](https://github.com/shu251/https://github.com/shu251/V4_tagsequencing_18Sdiversity_q1) repo. This protocol follows [Rideout et al. 2014 PeerJ](https://peerj.com/articles/545/).

## Requirements:
* Contents of this repo
* [QIIME2](https://docs.qiime2.org/2018.4/) version 2018.4
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* R
* *See section below* Reference database to be used for downstream sequence clustering & taxonomy assignment. For 18S (microbial eukaryotic) work, I prefer [PR2](https://github.com/vaulot/pr2_database/wiki)
* To follow step by step instructions below, follow along with test files provided from here: zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1236641.svg)](https://doi.org/10.5281/zenodo.1236641)

## Prep reference database

For several steps (including: closed OTU clusering, chimera detection, and taxonomy assignment) you will need a reference database that is imported as an artifact into qiime2.

Use below commands to set up this reference database with accompanying taxonomy information. Here, I'm using [PR2](https://github.com/vaulot/pr2_database/wiki) v4.10 (May 2018).

You will need to download both the fasta file and a taxonomy file. Then you'll need to import as artifacts (*create a .qza file*), [see instructions here](https://docs.qiime2.org/2018.4/tutorials/feature-classifier/).

Replace the below $PWD/pr2.fasta and pr2_tax.txt with appropriate path and reference fasta and taxonomy text file. 

You can import both the fasta DB file and taxonomy file as QIIME2 artifacts:
```
# activate the qiime2 environment
source activate qiime2-2018.4

# First import the database
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path $PWD/db/pr2.fasta \
  --output-path $PWD/db/pr2.qza

# Then the taxonomy file
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --source-format HeaderlessTSVTaxonomyFormat \
  --input-path $PWD/db/pr2_tax.txt \
  --output-path $PWD/db/pr2_tax.qza

```
Then you have the option to select the region within the fasta db specific to the primers you used and train the classifer:
```
# Select V4 region from the PR2 database
# Use appropriate forward and reverse primers
qiime feature-classifier extract-reads \
  --i-sequences $PWD/db/pr2.qza \
  --p-f-primer CCAGCASCYGCGGTAATTCC \
  --p-r-primer ACTTTCGTTCTTGATYRA \
  --p-trunc-len 150 \
  --o-reads $PWD/db/v4_extracts.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $PWD/db/v4_extracts.qza \
  --i-reference-taxonomy $PWD/db/pr2_tax.qza \
  --o-classifier $PWD/db/pr2_classifier.qza

# tip: make sure you version the databases and taxonomy files you're using. These are often updated so you want o keep them current, but also be able to match the appropriate fasta and taxonomy file.
  
```

## Prep directories and sample IDs to import into QIIME2 environment.

This protocol assumes that you have demultiplexed paired-end data. If you have data that has not been demultiplexed, you can add a demultiplex step at the begining of the protocol, which can also be done in [QIIME2](https://docs.qiime2.org/2018.4/). QIIME2 can also handle single-end data, but the steps are not described in this protocol, as paired-end data are recommended for V4 tag sequencing analysis.

If you would like to follow along to steps below, download 'test_fastq_large.zip' files from here: zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1236641.svg)](https://doi.org/10.5281/zenodo.1236641).
These sequences are a full set of raw sequences (R1 and R2 per sample) from the SPOT station. This provides the opportunity to run the below tutorial with a full dataset.

```
unzip test_fastq_large.zip # download from zenodo
mkdir raw_seqs_dir
mv Test*.fastq* raw_seqs_dir/
```

First, create a data manifest file in a text editor, replicate the format below. Use "manifest.txt" if using the test dataset. Every line must be comma separated, with sample-id, followed by the path of fastq file, followed by either "forward" or "reverse". 

```
# Example:
sample1,$PWD/raw_seqs_dir/Test01_full_L001_R1_001.fastq.gz,forward
sample1,$PWD/raw_seqs_dir/Test01_full_L001_R2_001.fastq.gz,reverse
sample2,$PWD/raw_seqs_dir/Test02_full_L001_R1_001.fastq.gz,forward
sample2,$PWD/raw_seqs_dir/Test02_full_L001_R2_001.fastq.gz,reverse
```
* Replace $PWD with your path
* The fastq files can be gziped. 
* List all of your fastq files. 
* Save the file and name it "manifest.txt".

## Run the run_trim.py script.
This script does two things:
(1) Runs trimmomatic as a preliminary trimming step. However, this does not remove PCR primers (e.g. Stoeck et al V4 region primers). However, the trimmomatic command can be used to remove excess barcodes or indices and run a preliminary trim of low quality bps.
(2) Generates three output files:
* trimmed_manifest.txt    A manifest file with trimmed fastq files
* trimming.log            A log file generated by Trimmomatic detailing the number of original and trimmed sequences
* trimmed_seq_dir         A directory containing all the trimmed fastq files
```
python $PWD/run_trim.py manifest.txt
```

## Start QIIME environment & import data

### Activate qiime2 environment:
```
source activate qiime2-2018.4
# This may depend on how you've installed QIIME2
```

### Load files into qiime formatted '.qza'
```
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path trimmed_manifest.txt --output-path demux.qza --source-format PairedEndFastqManifestPhred33
# trimmed_manifest.txt file was generated from the 'run_trim.py' command
```

### Remove primers
Below script is specific for using Stoeck et al. V4 primers. Make sure to change the '--p-front-f' and '--p-front-r' input sequences if you're using other primers.
```
qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza --p-cores 8 --p-front-f CCAGCASCYGCGGTAATTCC --p-front-r ACTTTCGTTCTTGATYRA --o-trimmed-sequences demux_trimmed.qza


# Output:
## QIIME is caching your current deployment for improved performance. This may take a few moments and should only happen once per deployment.
## Saved SampleData[PairedEndSequencesWithQuality] to: demux_trimmed.qza
```

### Option to summarize and view:
QIIME2 allows you to open an interactive summary. For each '.qza' file you can generate a related '.qzv' file. This file type can be opened via a view port (web server). [See this github for more information](https://github.com/leylabmpi/qiime_2_ll_quick_viewer).

```
# Generate the '.qzv' visual file from your '.qza' file.
qiime demux summarize --i-data demux_trimmed.qza --o-visualization demux_trimmed.qzv
# Open via local web server.
qiime_2_ll_quick_viewer --filename demux_trimmed.qzv
```

The output '.qza' (QIIME artifact) will be used directly for dada2 ASV calling (see below) or will require additional QC steps (directly below) for other OTU clustering approaches.

## Join paired sequences and quality filter

Join paired end reads with vsearch in qiime. After merging, final sequence has to be > 300 bps in length or it will be discarded ('--p-minmergelen'). With 18S sequencing, we've found it is more ideal (and you get more phylogenetic information) to rely on longer sequences. Thus, we use 250 bp x 250 bp PE MiSeq sequencing and then the overlap will only be a portion of the read. 
```
qiime vsearch join-pairs --i-demultiplexed-seqs demux_trimmed.qza --o-joined-sequences demux-joined.qza --p-minmergelen 300
# Output: Saved SampleData[JoinedSequencesWithQuality] to: demux-joined.qza

# Option to summarize merged read outputs and visualize:
qiime demux summarize --i-data demux-joined.qza --o-visualization demux-joined.qzv
qiime_2_ll_quick_viewer --filename demux-joined.qzv
```

### Quality filter & dereplicate sequences:
Here, the minimum Q-score is set to 20. Change '--p-min-quality' if desired. Sequences are dereplicated here to reduce run time.
```
# Filter by Q score
qiime quality-filter q-score-joined \
	--i-demux demux-joined.qza \
	--o-filtered-sequences demux-joined-filtered.qza \
	--o-filter-stats demux-joined-filter-stats.qza \
	--p-min-quality 20

# Option to summarize and view after the QC step:
qiime demux summarize --i-data demux-joined-filtered.qza --o-visualization demux-joined-filtered.qzv
qiime_2_ll_quick_viewer --filename demux-joined-filtered.qzv

# Dereplicate sequences:
qiime vsearch dereplicate-sequences --i-sequences demux-joined-filtered.qza --o-dereplicated-table derep_table.qza --o-dereplicated-sequences derep_seqs.qza
#Output:
# Saved FeatureTable[Frequency] to: derep_table.qza
# Saved FeatureData[Sequence] to: derep_seqs.qza
```
Now you have final QC'ed reads (chimera checking is after OTU/ASV clustering below), that are merged and filtered.

## Three flavors of grouping sequences by a defined parameter (e.g. Operational Taxonomic Units or Amplicon Sequence Variants) 
[See QIIME2 documentation](https://docs.qiime2.org/2018.4/tutorials/otu-clustering/)

### OTU clustering by closed reference or *de novo*

Below are preferred approaches for closed or *de novo* OTU clustering specific for 18S tag-sequencing (how much fun is protistan diversity?). [For open-reference OTU clustering, I still prefer to use QIIME1](https://github.com/shu251/V4_tagsequencing_18Sdiversity_q1). But check back as we continue migrating over to QIIME2.

## **Closed reference OTU clustering**

Generate a new directory to store output files from closed reference OTU clustering and then run. Here, I'm using 0.97 percent identity. Replace '--i-reference-sequences' with location of choice (and trained!) database.

```
mkdir closedref
qiime vsearch cluster-features-closed-reference \
	--i-table derep_table.qza \
	--i-sequences derep_seqs.qza \
	--i-reference-sequences $PWD/db/pr2.qza \
	--o-clustered-table closedref/table_closed_97.qza \
	--o-unmatched-sequences closedref/unmatched_seqs.qza \
	--o-clustered-sequences closedref/rep-seqs_closed_97.qza \
	--p-perc-identity 0.97 \
	--p-threads 8
# Output:
## Saved FeatureTable[Frequency] to: closedref/table_closed_97.qza
## Saved FeatureData[Sequence] to: closedref/rep-seqs_closed_97.qza
## Saved FeatureData[Sequence] to: closedref/unmatched_seqs.qza
```

## ***De novo OTU clustering***

```
mkdir denovo
qiime vsearch cluster-features-de-novo \
	--i-table derep_table.qza \
	--i-sequences derep_seqs.qza \
	--o-clustered-table denovo/table_dn_97.qza \
	--o-clustered-sequences denovo/rep-seqs_dn_97.qza \
	--p-perc-identity 0.97 \
	--p-threads 8
# Output:
## Saved FeatureTable[Frequency] to: denovo/table_dn_97.qza
## Saved FeatureData[Sequence] to: denovo/rep-seqs_dn_97.qza
```

## Chimera filtering after OTU clustering:

This happens in two steps
1. Table and rep-seqs artifact files are searched for chimeric sequences by either *de novo* or reference-based filtering
2. Then the table and rep-seqs files need to be filtered (removed) out.

Input files:

| OTU clustering  |     Table           | Representative sequences|
|-----------------|:-------------------:|------------------------:|
|Closed reference | table_closed_97.qza | rep-seqs_closed_97.qza  |
| *De novo*       | table_dn_97.qza     | rep-seqs_dn_97.qza      |

**Reference-based chimera checking**
Here, use the reference database you may have used for OTU clustering and the one you'll use for taxonomy assignment. [QIIME 2 documentation](https://docs.qiime2.org/2018.4/plugins/available/vsearch/uchime-ref/)
```
# Replace paths and .qza files where appropriate
qiime vsearch uchime-ref \
	--i-table $PWD/table.qza \
	--i-sequences $PWD/rep-seqs.qza \
	--i-reference-sequences $PWD/pr2.qza \
	--p-threads 8 \
	--output-dir chimera_ref_dir

# Remove chimeric sequences from table
qiime feature-table filter-features \
	--i-table $PWD/table.qza \
	--m-metadata-file chimera_ref_dir/nonchimeras.qza \
	--o-filtered-table $PWD/table_ref_nc.qza

# Remove chimeric sequences from representative sequences
qiime feature-table filter-seqs \
	--i-data $PWD/rep-seqs.qza \
	--m-metadata-file chimera_ref_dir/nonchimeras.qza \
	--o-filtered-data $PWD/rep-seqs_ref_nc.qza
```	

***De novo*** **chimera checking**
This will take longer than the reference-based chimera sequences. [QIIME2 documentation](https://docs.qiime2.org/2018.4/tutorials/chimera/?highlight=chimera)
```
# Run de novo chimera checking
qiime vsearch uchime-denovo \
	--i-table $PWD/table.qza \
	--i-sequences $PWD/rep-seqs.qza \
	--output-dir uchime_denovo_dir

# Remove chimeric sequences from table
qiime feature-table filter-features \
	--i-table $PWD/table.qza \
	--m-metadata-file uchime_denovo_dir/nonchimeras.qza \
	--o-filtered-table $PWD/table_dn_nc.qza

# Remove chimeric sequences from representative sequences
qiime feature-table filter-seqs \
	--i-data $PWD/rep-seqs.qza \
	--m-metadata-file uchime_denovo_dir/nonchimeras.qza \
	--o-filtered-data $PWD/rep-seqs_dn_nc.qza
```
Option to summarize chimera removal steps
```
qiime feature-table summarize \
	--i-table $PWD/table_dn_nc.qza \
	--o-visualization table_dn_nc.qzv
```

## Amplicon Sequence Variants (ASVs)

Since the dada2 ASV employs its own qc steps, error rate estimation and chimera detection/removal, use the .qza artifact files generated after the cutadapt step above. [Documentation](https://docs.qiime2.org/2018.4/plugins/available/dada2/denoise-paired/?highlight=dada2)

Below command includes error estimation rate of 2 (this is also the default) and pooled *de novo* chimera detection, meaning it conducts *de novo* chimeria sequence filtering among all samples, rather than for each individual sample. You may need to modify '--p-max-ee' and '--p-n-reads-learn' depending on your sequencing depth and quality.

```
mkdir ASVs
qiime dada2 denoise-paired \
	--i-demultiplexed-seqs demux_trimmed.qza \
	--o-table ASVs/table \
	--o-representative-sequences ASVs/rep-seqs \
	--p-n-threads 8 \
	--p-trunc-len-f 240 \
	--p-trunc-len-r 240 \
	--p-max-ee 2 \
	--p-n-reads-learn 1000000 \
	--p-chimera-method pooled \
	--o-denoising-stats ASVs/stats-dada2.qza

# Example output:
## Saved FeatureTable[Frequency] to: ASVs/table.qza
## Saved FeatureData[Sequence] to: ASVs/rep-seqs.qza
```

## Assign taxonomy

Use rep-seqs*nc.qza (OTU clustering) or rep-seqs.qza (ASV clustering) files to assign taxonomy to the representative sequence.
See above section on preparing reference databases.

*Two options:*
[Vsearch](https://peerj.com/articles/2584/)
```
qiime feature-classifier classify-consensus-vsearch \
	--i-query rep-seqs_closed_ref_nc.qza \
	--i-reference-reads $PWD/db/pr2.qza \
	--i-reference-taxonomy $PWD/db/pr2_tax.qza \
	--p-perc-identity 0.9 \
	--p-threads 8 \
	--o-classification tax_vsearch.qza
```

Alternatively, you can pre-train a classifier to assign taxonomy:
[Documentation here](https://docs.qiime2.org/2018.4/plugins/available/feature-classifier/classify-sklearn/?highlight=sklearn). For this one you can pre-fit your reference database using the primers for your amplicon [see here](https://docs.qiime2.org/2018.4/tutorials/feature-classifier/).
```
qiime feature-classifier classify-sklearn \
	--i-classifier $PWD/db/pr2_classifier.qza \
	--i-reads rep-seqs.qza \
	--o-classification tax_sklearn.qza
```

## Generate output tables
### OTU or ASV tables:
From clustering steps above you should have some version of these files:
* table*.qza
* rep-seqs*.qza

Use these files to generate tabulated table with rep seqs. **Change OTU to ASV if you're using the dada2 pipeline.
```
qiime tools export table-non-chimeric.qza \
	--output-dir export_dir # Change to ASV if using dada2 pipeline

biom convert \
	-i export_dir/feature-table.biom \
	-o export_dir/feature-table.tsv \	
	--to-tsv
```

### Export taxonomy information
Export to table:
```
qiime tools export tax.qza \
	--output-dir tax_dir
# tabulated file tax_dir/taxonomy.tsv
```

## Compile OTU and taxonomy tables for downstream analysis in R

Import .tsv OTU table and reformat so that column names are sample names. Compiled count information with taxonomy assignments. See script: compile_counts_tax.r

```
# Start R environment
# Import data from biom convert output
count<-read.delim("*.tsv", header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")

# Re-classify and re-format
colnames(count)<-as.character(unlist(count[2,]))
count1<-count[-c(1,2),]
colnames(count1)[1]<-"Feature.ID"
head(count1[1:2,]) # check dataframe structure
x<-dim(count1)[2];x # number of columns
# Convert class in count1 to numerics
count2<-count1
x<-dim(count2)[2];x # number of columns
count2[2:x] <- lapply(count2[2:x], function(x) as.numeric(as.character(x)))

# Get base stats and generate a table called "dataset_info"
seq_total<-apply(count2[2:x],2,sum) #number of sequences per sample
OTU_count<-colSums(count2[2:x]>0) # OTUs per sample
OTU_single<-colSums(count2[2:x]==1) # OTUs with only 1 seq
OTU_double<-colSums(count2[2:x]==2) # OTUs that have only 2 seqs
OTU_true<-colSums(count2[2:x]>2) # Number of OTUs with >2 seqs
dataset_info<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true)

# Totals for whole dataset:
sum(seq_total) # total number of sequences over whole dataset
dim(count2)[1] # total unique OTUs over whole dataset

# Option to remove global singletons
rowsum<-apply(count2[2:x],1,sum) # Number of sequences per OTU (for whole dataset)
count.no1 = count2[ rowsum>1, ] # Remove rowsums =< 1 (only one sequence for 1 OTU)
dim(count2)[1] - dim(count.no1)[1] # Total number of global singletons removed from data
## Switch values above to remove global doubletons as well!

# Import taxonomy file and join with count data:
taxname<-read.delim("$PWD/tax_dir/taxonomy.tsv", header=TRUE, row.names=NULL)
library(plyr) # Load plyr library into R environment
counts_wtax<-join(count.no1, taxname, by="Feature.ID", type="left", match="first")

# Write output table with OTU or ASV cluster information and taxonomy IDs:
write.table(counts_wtax, file="OutputTable_wtax.txt", quote=FALSE, sep="\t", row.names=FALSE)
```


### Get sequence run stats

In the above QC steps, whenever a .qza file was worked on, you had the option to run this:
qiime demux summarize --i-data [XXX.qza] --o-visualization XXX.qzv
This generates a 'visualization' version of the same file (.qzv). Then you can run the quick viewer:
qiime_2_ll_quick_viewer --filename [XXX.qzv]. This will open a page where you can download a .csv file to your local computer

* demux.qzv - get input sequence information
* demux-trimmed.qzv - trimmed sequences (primer trimmed)
* demux-joined.qzv - merged paired end reads
* demux-joined-filtered.qzv - quality score filtered

### Last updated Sarah Hu 05-08-2018

