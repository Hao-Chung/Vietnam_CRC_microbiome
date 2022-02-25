library("dada2")
packageVersion("dada2")
library(stringr)
library(Biostrings)
library(phyloseq)
library(ggplot2)

path <- "./trimmed_fastq/"
list.files(path)

## Forward and reverse fastq file names have format ${x}_1.fastq.gz and ${x}_2.fastq.gz
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names=TRUE))
## Extract sample names, assuming file names have consistent format
sample.names <- sapply(strsplit(basename(fnFs), "_"),`[`,1)

## Have a look of the data
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[5:8])

## Remove the primer on forward and reverse independently, using a fixed cut-off in length for each sequence.
filt.fwd <- str_c(sample.names, "_noprimer_1.fastq.gz", sep="")
filt.rev <- str_c(sample.names, "_noprimer_2.fastq.gz", sep="")
## Truncate the length of forward and reverse reads to retain only high quality reads, as well as removing primers
for (x in 1:length(sample.names)){
      filterAndTrim(fwd=fnFs[x], filt=filt.fwd[x], rev=fnRs[x], filt.rev = filt.rev[x], 
                    truncLen = c(240,180), maxN=0, trimLeft = c(19,20), rm.phix = TRUE, multithread = TRUE)
      print(cat('Finish trimming primers for sample', sample.names[x], sep=" "))
}
## After this step, the number of paired reads in each sample is greatly reduced

path2 <- "./noprimers_fastq/"
list.files(path2)
## Denote the new path for the fastq files that have removed the primers
fwdfin <- sort(list.files(path2, pattern="_1.fastq.gz", full.names=TRUE))
revfin <- sort(list.files(path2, pattern="_2.fastq.gz", full.names=TRUE))

plotQualityProfile(fwdfin[1:2])
plotQualityProfile(revfin[20:22])

## Learn the error rates using the trimmed sequences
errF <- learnErrors(fwdfin, multithread = TRUE)
errR <- learnErrors(revfin, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

## Sample inference, default is without pooling
dadaFs <- dada(fwdfin, err = errF, multithread = TRUE)
dadaRs <- dada(revfin, err = errR, multithread = TRUE)

## Merging paired reads
## Allow for maximum of 2 mismatches in the overlap region
mergers <- mergePairs(dadaFs, fwdfin, dadaRs, revfin, verbose=TRUE, maxMismatch = 2, returnRejects = FALSE)

seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab))) # most of the unique sequences still fall into the bacterial 16S expected size (251-256)

## Extra steps to ensure that no 16S rDNA was left in the pre-processing
shortseq <- getSequences(seqtab)[which(nchar(getSequences(seqtab)) %in% 245:249)]
shortseq.taxa <- assignTaxonomy(shortseq, "./silva_nr99_v138_train_set.fa.gz", multithread = TRUE, outputBootstraps = TRUE)
rownames(shortseq.taxa$boot) <- NULL
shortseq.taxa$boot
rownames(shortseq.taxa$tax) <- NULL
shortseq.taxa$tax
## There are two sequences of length 249 bp classified correctly to Alloprevotella and Solobacterium
# so consider to include 249 bp as the minimum cutoff of 16S rRNA length. 
## Now for merged sequences longer than 256 bp.
longseq <- getSequences(seqtab)[which(nchar(getSequences(seqtab)) %in% 257:260)]
longseq.taxa <- assignTaxonomy(longseq, "./silva_nr99_v138_train_set.fa.gz", multithread = TRUE, outputBootstraps = TRUE)
rownames(longseq.taxa$tax) <- NULL
longseq.taxa$tax
rownames(longseq.taxa$boot) <- NULL
longseq.taxa$boot
# No bacterial 16S sequences were found in sequences greater than 256 bp. 

## Checking to see if there is any 16S sequence fall within the 245-249 bp length
merged.245_249 <- DNAStringSet(getSequences(seqtab)[which(nchar(getSequences(seqtab)) %in% 245:249)], use.names = TRUE)
merged.245_249 <- merged.245_249[order(width(merged.245_249))]
names(merged.245_249) <- 1:65
writeXStringSet(merged.245_249,filepath = "./merged_sequences_245_249.fasta",format = "fasta")

########################################################
### Remove most contaminated sequences ###############
###################################################

## Based on previous interpretation, the minimum cutoff for 16S length is 249 bp. 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 249:256]
sum(seqtab2)/sum(seqtab) #~89% sequences were retained.

# remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

sum(seqtab.nochim)/sum(seqtab2) # the total number of reads is still mainly retained

## checking the steps of processing 
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
orig_depth <- read.table(file="./Original_fastqc_sequencing_depth.txt", sep="\t", header = FALSE)
rownames(orig_depth) <- orig_depth$V1
trim_depth <- read.table(file="./Trimmed_adapter_fastqc_sequencing_depth.txt", sep='\t', header=FALSE)
rownames(trim_depth) <- trim_depth$V1

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged", "trimlen", "nonchim")
rownames(track) <- sample.names
head(track)

## Combine track with original sample depth
rownames(track)
orig_depth <- orig_depth[match(rownames(track), rownames(orig_depth)),]
trim_depth <- trim_depth[match(rownames(track), rownames(trim_depth)),]
identical(rownames(track), rownames(orig_depth))
identical(rownames(track), rownames(trim_depth))
track.new <- cbind(orig_depth$V2, trim_depth$V2, track)
colnames(track.new)[c(1,2)] <- c("original", "trim_adapter")
head(track.new)
track.new <- as.data.frame(track.new)

## Assign taxonomy to the created nonchimera seqtab
taxa <- assignTaxonomy(seqtab.nochim, "./silva_nr99_v138_train_set.fa.gz", multithread=TRUE, outputBootstraps = TRUE)
taxa.species <- addSpecies(taxa$tax,"silva_species_assignment_v138.fa.gz")

taxa.print <- taxa.species
rownames(taxa.print) <- NULL
write.table(taxa.print, file="./27EN_210str_dada2_taxonomy_wcontam.csv",quote = FALSE,sep = ",",col.names = TRUE,row.names = TRUE)

taxa$boot
taxa.confidence <- taxa$boot
rownames(taxa.confidence) <- NULL
write.table(taxa.confidence, file="./27EN_210str_dada2_taxonomy_wcontam_bootstrap.csv", quote=FALSE, sep=",", col.names=TRUE,
            row.names=TRUE)

## Doing some exploration 
which(is.na(taxa.print[,1]))
which(!(taxa.print[,1] %in% c("Bacteria", "Archaea")))
taxa.print[which(taxa.print[,6] == "Fusobacterium"),]

## Output the all sequences that pass the chimera removing steps
nochim.all.seqs <- DNAStringSet(getSequences(seqtab.nochim), use.names = TRUE)
names(nochim.all.seqs) <- str_c('seq',1:4078, sep = "")
writeXStringSet(nochim.all.seqs,filepath = "./27EN_4078seqs_wcontam.fasta",format = "fasta")

### Evaluating the DADA2's accuracy on the mock community (zymo control)
rownames(seqtab.nochim) <- sample.names
unqs.mock <- seqtab.nochim["Zymo",]
unqs.mock <- sort(unqs.mock[unqs.mock > 0], decreasing = TRUE) ## Drop ASVs absent in the mock

##  The output sequences of 4078 ASVs were aligned and then used for tree building. Rough phylogeny showed that 
# there are 170 sequences having abnormal branch lengths
contam1 <- read.table(file="./pasta_alignment/27EN_contam_base_ontree.ID.txt", header=FALSE, stringsAsFactors = FALSE)
contam1 <- unique(contam1$V1)
rownames(taxa.print) <- str_c('seq',1:4078, sep = "")
taxa.print[sort(match(contam1, rownames(taxa.print))),]

rownames(taxid) <- str_c('seq',1:4078, sep = "")
taxid[sort(match(contam1, rownames(taxid))),]

## Recheck the rough phylogeny after removing these 170 ASVs
asv.tre <- read.newick(file="./pasta_alignment/27EN_4078seqs_wcontam_pasta.output.newick")
asv.trim.tre <- drop.tip(asv.tre, tip =contam1)
plot(asv.trim.tre)
write.tree(asv.trim.tre, file="./pasta_alignment/27EN_3908seqs_trim.newick")
retain.asv <- asv.trim.tre$tip.label[order(as.numeric(str_extract(asv.trim.tre$tip.label,pattern = "[0-9]+")), decreasing = FALSE)]
write.table(retain.asv, file="./27EN_3908sesqs_postfilter_ID.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)
## The trimmed tree showed that there is no abnormal branch length left, indicating that the human contaminants have been removed

## Checking the sequences that are not listed in contam1
filt1 <- which(!(rownames(taxa.print) %in% contam1))
which(is.na(taxa.print[filt1, 1])) ## For the remaining, no NA/unclassified at the kingdom level 
taxa.print[match(names(which(is.na(taxa.print[filt1, 2]))), rownames(taxa.print)),]

######################################
# Hand it off to phyloseq object #### 
######################################

## Transform the dada2 output into phyloseq object
samples.out <- rownames(seqtab.nochim)
write.table(samples.out, file="./27EN_sampledata_rough.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)

## Input the sample data table
meta.rough <- read.table(file="./27EN_sampledata_rough.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
rownames(meta.rough) <- meta.rough$sampleID
meta.rough <- meta.rough[,-c(1)]

# Building the phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows =FALSE),
         sample_data(meta.rough), tax_table(taxa.species))
taxa.dna <- DNAStringSet(taxa_names(ps))
names(taxa.dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("seq",seq(ntaxa(ps)))

## removing contaminant taxa from the phyloseq object
ps <- prune_taxa(asv.trim.tre$tip.label,ps)

tax_table(ps)
tax_table(ps)[which(otu_table(ps)['Zymo',] > 10 ),]
otu_table(ps)[,'seq1']
rowSums(otu_table(ps))

# Doing a few exploration
tax_table(ps)[which(otu_table(ps)['T52',] > 0 ),]

### Blank 
tax_table(ps)[which(otu_table(ps)['NTC-Primary',] > 0 ),]
otu_table(ps)['NTC-Primary', which(otu_table(ps)['NTC-Primary',] > 0 )]
tax_table(ps)[which(otu_table(ps)['NTC-index',] > 0 ),]
otu_table(ps)['NTC-index', which(otu_table(ps)['NTC-index',] > 0 )]

## DNA extraction kit
tax_table(ps)[which(otu_table(ps)['FastDNA-S',] > 0 ),]
otu_table(ps)['FastDNA-S', which(otu_table(ps)['FastDNA-S',] > 0 )]
kitome.fast1 <- taxa_names(ps)[which(otu_table(ps)['FastDNA-S', ] > 0 )]

tax_table(ps)[which(otu_table(ps)['FastDNA-L',] > 0 ),]
otu_table(ps)['FastDNA-L', which(otu_table(ps)['FastDNA-L',] > 0 )]
kitome.fast2 <- taxa_names(ps)[which(otu_table(ps)['FastDNA-L', ] > 0 )]

tax_table(ps)[which(otu_table(ps)['Relia-TE',] > 0 ),]
otu_table(ps)['Relia-TE', which(otu_table(ps)['Relia-TE',] > 0 )]
kitome.relia1 <- taxa_names(ps)[which(otu_table(ps)['Relia-TE', ] > 0 )]

tax_table(ps)[which(otu_table(ps)['Relia-SW',] > 0 ),]
otu_table(ps)['Relia-SW', which(otu_table(ps)['Relia-SW',] > 0 )]
kitome.relia2 <- taxa_names(ps)[which(otu_table(ps)['Relia-SW', ] > 0 )]

# Combine and reorder
kitome.id <- unique(c(kitome.fast1, kitome.fast2, kitome.relia1, kitome.relia2))
kitome.id <- kitome.id[order(as.numeric(str_extract(kitome.id, "[0-9]+")), decreasing = FALSE)]

kitome.presence <- c()
kitome.abund <- list()
for (x in kitome.id) {
      kitome.presence[x] <- length(which(otu_table(ps)[,x]>0))
      kitome.abund[[x]] <- summary(otu_table(ps)[which(otu_table(ps)[,x] >0), x])
}
## Kitome presence in only one sample is automatically seen as contamination
kitome.abund[names(which(kitome.presence == 1))]
kitome.abund[names(which(kitome.presence > 1 & kitome.presence < 5))]
tax_table(ps)[names(which(kitome.presence > 1 & kitome.presence < 5)),]
# seq476 (Prevotella), seq732(Rothia) are not contamination, the others are view as contamination
# seq783 (Ramlibacter), seq1034(Aquabacterium), seq1176 (Comamonadaceae), seq1258(Microbacterium), seq1261, 
# seq1576 (Herbaspirillum), seq1772 (Massilia), seq1847 (Pseudomonas), seq1959 and others are contamination

kitome.abund[names(which(kitome.presence >=5))]
tax_table(ps)[names(which(kitome.presence >=5)),]
# Contamination identified: seq37 (Salmonella), seq70 (Lactobacillus fermentum), seq132 (Bacillus subtilis),
# seq161 (Staphylococcus), seq178 (Listeria), seq174 (Pseudomonas), seq65(Enterococcus)??
# seq93 (Ralstonia), seq497 (Burkholderia group), seq750 (Ralstonia), seq550 (Anaeroglobus), seq791 (Pseudomonas),
# seq935 (Paucibacter), seq1039(Acinetobacter), seq1154(Ananerobacillus), seq2140 (Enhydrobacter)
otu_table(ps)[,'seq65']

kitome <- c(names(which(kitome.presence ==1)), names(which(kitome.presence > 1 & kitome.presence < 5)),
  'seq37', 'seq70', 'seq132', 'seq161', 'seq178', 'seq174', 'seq93', 'seq497', 'seq750', 'seq550', 'seq791',
  'seq935', 'seq1039', 'seq1154', 'seq2140')
## The identified kitome, including species in the mock community
kitome.fin <- kitome[which(!(kitome %in% c("seq476",'seq732')))]
kitome.abund[kitome.fin]

############################################
## Looking at mock community first 
mock.only <- prune_samples('Zymo',ps)
mock.only <- prune_taxa(taxa_sums(mock.only)>0, mock.only)
mock.only <- prune_taxa(taxa_sums(mock.only)>100, mock.only )
otu_table(mock.only)
plot_bar(mock.only, fill='Genus') + geom_bar(aes(colour=Genus, fill=Genus), stat='identity', position='stack') + 
      scale_fill_manual(values = brewer.pal(n=8, name='Set3')) + 
      scale_color_manual(values=brewer.pal(n=8, name='Set3')) + theme_bw() + 
      xlab("Mock community") + ylab("Abundance")

#################################################################
############## Checking the retainment of reads after each filtering step
#############################################################
library(reshape2)
require(RColorBrewer)
identical(rownames(meta.rough), rownames(track.new))# TRUE
track.new$sample <- meta.rough$sample_type
track.new.m <- melt(track.new[,c(1,2,5,6,7,8)],id.vars = c("sample"))
head(track.new.m)
colnames(track.new.m) <- c("sample_type", "Processing", "value")
track.new.m$sample_type <- factor(track.new.m$sample_type, levels=c("saliva","biopsy", "polyp", 
                                                                    "healthy_tissue", "tumour", "mock", "blank"))
ggplot(track.new.m, aes(x=sample_type, y=value, fill=Processing)) + geom_boxplot(width=0.7) + 
      scale_fill_manual(labels=c("original", "trim_adapter", "post-merging", "cut-length", "remove-chimera"),
                        values=c("azure",brewer.pal(n=8, name="Accent")[c(1:4)])) + theme_bw() + 
      xlab("Sample type") + ylab("Sequencing depth") + 
      theme(axis.title = element_text(size=12), axis.text.x=element_text(size=10, angle=50, hjust=1))

## Exclude mock community, negative controls, Archaea taxa and kitome taxa from the data
samples.only <- sample_names(ps)[which(!(sample_data(ps)$Group %in% c("mock", "kit", "blank")))]
crc.all <- prune_samples(samples.only, ps)
crc.all <- prune_taxa(taxa_names(ps)[which(!(taxa_names(ps) %in% kitome.fin))], crc.all)
crc.all <- prune_taxa(taxa_names(ps)[which(tax_table(ps)[,1] != 'Archaea')], crc.all)
crc.all <- prune_taxa(taxa_sums(crc.all)>0, crc.all)

######################
## Do some filtering to remove low abundant taxa
###########

length(taxa_names(crc.all)[which(taxa_sums(crc.all) >=10 )])

## removing taxa that are lower than 10 counts, keeping singleton 
keep <- genefilter_sample(crc.all, filterfun_sample(function(x) x > 10), A = 1)
crc.keep <- prune_taxa(keep, crc.all)
crc.keep

