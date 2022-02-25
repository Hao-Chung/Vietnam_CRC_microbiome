# Vietnam_CRC_microbiome
This page contains source data (RData, csv, newick files) and Rmarkdown codes required to replicate the analysis on 16S rRNA sequencing of colorectal microbiome (saliva and gut mucosa) from Vietnamese patients

## Source data:

1/ 27EN_compiled_metadata.csv: anonymous patient data.
    
2/ 27EN_saliva_fil_865seqs.snp.fa.contree.midroot.newick: IQ-tree constructed phylogeny of 865 ASVs detected in the salivary microbiome.
    
3/ tissue_fil_taxa_names_1073seqs_pasta.aln.contree.reroot.newick: IQ-tree constructed phylogeny of 1,073 ASVs detected in the gut mucosal microbiome.
    
4/ 27EN_phyloseq_CRC_keep.RData: phyloseq data file containing the otu_table and taxonomy of all ASVs filtered in 203 microbiome samples. 

## Analysis codes: 

1/ dada2_processing.R: preprocessing of raw fastq files into phyloseq object. Raw fastq files available under the Bioproject PRJNA791834. 

2/ 27EN_CRC_analysis_updated.Rmd: analysis of microbiome data, including codes for generating figures.
