if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

BiocManager::install(c("biomaRt", "Biostrings"))
install.packages(c('rentrez', 'dplyr', 'readr', 'seqinr', 'stringr'))

library(biomaRt)
library(Biostrings)
library(rentrez)
library(dplyr)
library(readr)
library(seqinr)
library(stringr)

setwd("~/R")

mart_finder <- function(file_name_1) {
  list_1 = listMarts(host='https://www.ensembl.org')
  write.csv(list_1, file_name_1, row.names=FALSE)
}

database_finder <- function(mart_name, file_name_2) {
  list_2 = useMart(biomart=mart_name, host='https://www.ensembl.org')
  list_2_new = listDatasets(list_2)
  write.csv(list_2_new, file_name_2, row.names=FALSE)
}

filters_attributes <- function(type, species, file_1, file_2) {
  species_dataset = useEnsembl(biomart=type, dataset=species)
  list_1 = listFilters(species_dataset)
  list_2 =listAttributes(species_dataset)
  write.csv(list_1, file_1, row.names=FALSE)
  write.csv(list_2, file_2, row.names=FALSE)
}

dataset_retrieve <- function(type, species, chrom, file_name) {
  species_dataset = useEnsembl(biomart=type, dataset=species)
  species_query <- getBM(attributes=c('refseq_mrna', 'refseq_peptide', 'ensembl_gene_id',     
                                      'external_gene_name', 'description', 'start_position', 'end_position', 'strand',  
                                      'chromosome_name', 'name_1006'), filters =
                           'chromosome_name', values =chrom, mart = species_dataset)
  write.csv(species_query, file_name, row.names=FALSE)
  species_csv = read.csv(file_name, na.strings = c("", "NA"))
  species_csv = species_csv %>% na.omit()
  write.csv(species_csv, file_name, row.names=FALSE)
}

gene_list <- function(species, chrom, species_2_id, species_2_gene_name, file_name) {
  species_dataset = useEnsembl(biomart="ensembl", dataset=species)
  gene_list_query <- getBM(attributes=c('ensembl_gene_id','external_gene_name',
                                        species_2_id, species_2_gene_name), filters =
                             'chromosome_name', values =chrom, mart = species_dataset)
  write.csv(gene_list_query, file_name, row.names=FALSE)
  genes_csv = read.csv(file_name, na.strings = c("", "NA"))
  genes_csv = genes_csv %>% na.omit()
  write.csv(genes_csv, file_name, row.names=FALSE)
}

gene_list_dataset_1_filter <- function(species, gene_list, species_filter, filter_gene) {
  dataset = read.csv(species)
  genes = read.csv(gene_list)
  list_1 = genes %>% select(external_gene_name)
  list_1_column = unique(list_1)
  list_1_vector = unlist(list_1_column)
  query_1 = dataset[dataset$external_gene_name %in% list_1_vector, ]
  write.csv(query_1, species_filter, row.names=FALSE)
  
  list_2 =  query_1 %>% select(external_gene_name)
  list_2_column = unique(list_2)
  list_2_vector = unlist(list_2_column)
  query_2 = genes[genes$external_gene_name %in% list_2_vector, ]
  write.csv(query_2, filter_gene, row.names=FALSE)
}

gene_list_dataset_2_filter <- function(species, gene_list, column_name, file_name_1, file_name_2) {
  dataset = read.csv(species)
  genes = read.csv(gene_list)
  list_1 = genes %>% select(column_name)
  list_1_column = unique(list_1)
  list_1_vector = unlist(list_1_column)
  query_1 = dataset[dataset$ensembl_gene_id %in% list_1_vector, ]
  write.csv(query_1, file_name_1, row.names=FALSE)
  list_2 =  query_1 %>% select(ensembl_gene_id)
  list_2_column = unique(list_2)
  list_2_vector = unlist(list_2_column)
  query_2 = genes[genes[, column_name] %in% list_2_vector, ]
  write.csv(query_2, file_name_2, row.names=FALSE)
}

dataset_1_final_filter<- function(species, gene_list, file_name) {
  dataset = read.csv(species)
  genes = read.csv(gene_list)
  list_1 = genes %>% select(external_gene_name)
  list_1_column = unique(list_1)
  list_1_vector = unlist(list_1_column)
  query_1 = dataset[dataset$external_gene_name %in% list_1_vector, ]
  write.csv(query_1, file_name, row.names=FALSE)
}  

gene_ontology_filter <- function(file, go_term, go_name_filter) {
  filtered_species = read.csv(file)
  query = filtered_species[filtered_species$name_1006 == go_term,]
  write.csv(query, go_name_filter, row.names=FALSE)
}

ref_seq_list <- function(file_name, column_name, gene_name, name) {
  file = read.csv(file_name)
  selected = file[, c(column_name, 'external_gene_name')]
  new_list = unique(selected)
  query = new_list[new_list$external_gene_name %in% c(gene_name),]
  write.csv(query, name, row.names=FALSE)
}

ref_seq_sequence <- function(db_type, id, file_name) {
  net_handle <- entrez_fetch(db=db_type, id=id, rettype="fasta", retmode='text')
  write(net_handle, file = file_name)
}

pairwise_alignment <- function(file_1, file_2, matrix, open_gap, extend_gap, file_name) {
  species_1 <- read.fasta(file_1)
  species_1_character <- unlist(species_1)
  species_1_upper <- lapply(species_1_character, toupper)
  species_1_unlist <- unlist(species_1_upper)
  species_1_string <- toString(species_1_unlist)
  species_1_comma = str_replace_all(species_1_string,",","")
  species_1_space = str_replace_all(species_1_comma," ","")
  
  
  species_2 <- read.fasta(file_2)
  species_2_character <- unlist(species_2)
  species_2_upper <- lapply(species_2_character, toupper)
  species_2_unlist <- unlist(species_2_upper)
  species_2_string <- toString(species_2_unlist)
  species_2_comma = str_replace_all(species_2_string,",","")
  species_2_space = str_replace_all(species_2_comma," ","")
  
  
  alignment <- pairwiseAlignment(species_1_space, species_2_space, type="global",
                                 substitutionMatrix = matrix,
                                 gapOpening = open_gap,
                                 gapExtension = extend_gap,
                                 scoreOnly = FALSE)
  
  
  
  writePairwiseAlignments(alignment, file=file_name, Matrix = matrix, block.width=60)
}

mart_finder('mart_list_R.csv') 

database_finder('ENSEMBL_MART_ENSEMBL', 'database_list_R.csv')

filters_attributes('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', 'h_filter_R.csv', 'h_attrib_R.csv')

filters_attributes('ENSEMBL_MART_ENSEMBL', 'mmusculus_gene_ensembl', 'm_filter_R.csv', 'm_attrib_R.csv')

dataset_retrieve('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', '5', 'species_1_R.csv')

dataset_retrieve('ENSEMBL_MART_ENSEMBL','mmusculus_gene_ensembl', '18', 'species_2_R.csv')

gene_list('hsapiens_gene_ensembl', '5', 'mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_associated_gene_name', 'genes_R.csv')

gene_list_dataset_1_filter('species_1_R.csv',  'genes_R.csv','species_1_filter_R.csv',
                           'filtered_gene_R.csv') 

gene_list_dataset_2_filter('species_2_R.csv', 'filtered_gene_R.csv', 'mmusculus_homolog_ensembl_gene',  
                           'species_2_filter_final_R.csv', 'filtered_gene_final_R.csv')

dataset_1_final_filter('species_1_filter_R.csv', 'filtered_gene_final_R.csv', 'species_1_filter_final_R.csv')

gene_ontology_filter('species_1_filter_final_R.csv', 'plasma membrane', 'species_1_go_R.csv')

gene_ontology_filter('species_2_filter_final_R.csv', 'plasma membrane', 'species_2_go_R.csv')

ref_seq_list('species_1_go_R.csv', 'refseq_peptide', 'APC', 'species_1_ref_R.csv')

ref_seq_list('species_2_go_R.csv', 'refseq_peptide', 'Apc', 'species_2_ref_R.csv')

ref_seq_sequence('protein', 'NP_001394379', 'H_APC_ref_seq_R.fasta')

ref_seq_sequence('protein', 'NP_001347909', 'M_Apc_ref_seq_R.fasta')

pairwise_alignment('H_APC_ref_seq_R.fasta', 'M_Apc_ref_seq_R.fasta', 'BLOSUM62', -10, -0.5, 'alignment_R.txt')
