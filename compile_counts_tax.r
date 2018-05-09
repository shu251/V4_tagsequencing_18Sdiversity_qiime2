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
