# R Pipeline Tutorial to Query and Manipulate Data from BioMart and RefSeq in Order to Perform a Pairwise Alignment Between Two Species 

### Overview

The following tutorial gives a step by step guide on how to successfully use the pipeline in order to get the desired results. The R script without the explanations has also been provided and can be found above.

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.17")
```
The following chunk installs the package BiocManger. This package originates from BioConductor and is required to use other packages that originate from there.

```R
BiocManager::install(c("biomaRt", "Biostrings"))
install.packages(c('rentrez', 'dplyr', 'readr', 'seqinr', 'stringr'))
```
The following chunk installs all the necessary packages needed to run the pipeline. The first line of chunk installs all of the packages that originated from BioConductor and the second line of chunk installs all of the packages that originated from CRAN. The packages from BioConductor are installed through the use of BiocManager::install() function and the packages from CRAN are installed through the use of install.packages() function. If more than one package is being installed, then the c() function is used within the previously mentioned functions.

```R
library(biomaRt)
library(Biostrings)
library(rentrez)
library(dplyr)
library(readr)
library(seqinr)
library(stringr)
```
The following chunk loads all the packages that will be used for the pipeline. The function library() allows for successfully downloaded packages to be called and used.

```R
setwd("~/R")
```
The following chunk sets the directory as to where files can be saved to or accessed from.

```R
mart_finder <- function(file_name_1) {
  list_1 = listMarts(host='https://www.ensembl.org')
  write.csv(list_1, file_name_1, row.names=FALSE)
}
```
The following chunk creates a function that will return a csv that contains a list of all the possible BioMart Ensembl marts. The function has one argument. The argument file_name_1 refers to the name of the csv file for all possible marts. The list of marts was saved as list_1 and converted to a csv file. The list of all possible marts can be used to determine what mart will be entered for the mart_name argument of the database_finder function.

```R
database_finder <- function(mart_name, file_name_2) {
   list_2 = useMart(biomart=mart_name, host='https://www.ensembl.org')
   list_2_new = listDatasets(list_2)
   write.csv(list_2_new, file_name_2, row.names=FALSE)
}
```
The following chunk creates a function that will return a csv that contains a list of all the possible BioMart Ensembl databases. The function has two arguments. The first argument mart_name refers to the name of the selected mart. This mart can be found from the list that was created from the mart_finder function. The second argument file_name_2 refers to the name of the csv file for possible datasets. The list of datasets is saved as list_2 and converted to a csv file. The list of all possible databases can be used to determine what species will be entered for the species argument of the dataset_retrieve function.

```R
filters_attributes <- function(type, species, file_1, file_2) {
  species_dataset = useEnsembl(biomart=type, dataset=species)
  list_1 = listFilters(species_dataset)
  list_2 =listAttributes(species_dataset)
  write.csv(list_1, file_1, row.names=FALSE)
  write.csv(list_2, file_2, row.names=FALSE)
}
```
The following chunk creates a function that will return csvs containing all the possible attributes and filters for a specified dataset. This function has four arguments. These arguments are type, species, file_1, and file_2. The argument type refers to the mart that will be used and the argument species refers to the selected database. The arguments of file_1 and file_2 refer to what the csvs will be called. The species_dataset variable is set equal to the useEnsembl subpackage, which requires a biomart and dataset name to access a specified database. The biomart variable is set equal to the type of mart that will be used and the dataset variable is set equal to the species argument. The attributes and filters for species_dataset variable can be found by using the functions listFilters() and listAttributes() and adding species_dataset inside the parenthesis. The filters and attributes are then saved as list_1 and list_2. These lists are then written to csvs through the function write.csv(). 

```R
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
```
The following chunk creates a function that will return a csv file containing queried data from BioMart. To get this file, four arguments must be entered. These arguments are type, species, chrom, and file_name. The argument of type refers to the mart that will be used. The species argument refers to the database from which data will be queried. Furthermore, this database contains genetic information related to a particular species. A list of possible species can be found using the database_finder function mentioned above. The chrom argument refers to the chromosomal location in which the genes are located. The file_name argument refers to what the output csv file will be named. The species_dataset refers to the variable that is set equal to the useEnsembl subpackage, which requires a biomart and dataset name to access a specified database. The biomart variable is set equal to the type of mart that will be used and the dataset variable is set equal to the species argument. The species_query refers to the variable that is set equal to the queries that will be retrieved based on the specified attributes and filters. The queries are called through the function getBM(), which contains all of the filters and attributes used for the query. The queries are filtered by chromosomal location, which is set equal to the argument of chrom. All the attributes and filters can be found using the filters_attributes function mentioned above. The species_query is then saved to a csv through the argument file_name. This csv is then read back in and set equal to the variable species_csv. In addition to this, the reading in of the csv also replaces any blanks with NA. The species_csv variable is called again and the na.omit() function is used to remove any NAs present in the csv. The removal of NA allows for easier manipulation of the data. The species_csv is then converted to a csv file and given a name based on the file_name argument.

```R
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
```
The following chunk creates a function that lists the genes of a particular species and its corresponding homologs to another specified species. This information is then saved to a csv file. This function has five arguments that must be entered. The first argument is the name of the database from which the information will be queried. The second argument is the chromosomal location in which the genes are located. The third argument is the gene ids for the second specified species, which is used as an attribute when querying the data. The fourth is the gene names for the second specified species, which is also used as an attribute when querying the data. The fifth argument is file_name, which is used to name the newly created csv. To query them, it accesses a database as previously mentioned in the dataset_retrieve function. While this function is like the dataset_retrieve function, the attributes for this function are different. In particular, the attributes list the gene id and gene name for the specified species and homologs. All the queries are filtered by chromosomal location, which is set equal to the argument of chrom. All of this information is set equal to the variable of gene_list_query. The gene_list_query is then saved to a csv and given a name through the argument file_name. This csv is then read back in and set equal to the variable species_csv. In addition to this, the reading in of the csv also replaces any blanks with NA. The species_csv variable is called again and the na.omit() function is used to remove any NAs present in the csv. The removal of NA allows for easier manipulation of the data. The species_csv is then converted to a csv file and given a name based on the file_name argument.

```R
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
```
The following chunk creates a function that filters out queries that don't appear on the gene list for the first species dataset. This function has four arguments, which are species, gene_list, species_filter, and species_gene. The first argument refers to the filename of the species dataset that is read in and the second argument refers to the filename of the gene list that is read in. The argument species_filter refers to the file that is created after absent genes are filtered out of the dataset. The argument filter_gene refers to the file that is created after genes from the gene_list are filtered out if they are absent from the dataset. The gene list and the dataset for the first species are read in and set equal to the variables of dataset and genes. To get all the unique genes for the first species, a variable called list_1_column is created. This variable is a list that is made up of genes from the gene variable. The variable list_1_column is created, which removes any duplicate gene names from list_1 and helps to make the overall list be as small as possible. The variable list_1_vector is created to convert the list_1_column to a vector. The variable query_1 is created and returns queries from the dataset variable that have an external gene name that is present in list_1_vector. These queries are then saved to a csv file and given a name with the species_filter argument. The gene list then must be updated to filter out genes that were not present in the dataset. A variable called list_2 is created and makes a list of all the gene names from the filtered dataset. Like list_1, a variable called list_2_column is created to remove duplicate names that might be present in the list_2. The variable list_2_vector is created to convert the list_2_column to a vector. The variable query_2 is created and returns queries from the genes variable that have an external gene name that is present in list_2_vector. These queries are then saved to a csv file and given a name through the use of the filter_gene argument.

```R
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
```
The following chunk creates a function that filters out queries that don't appear on the gene list for the second species dataset. This function has five arguments, which are species, gene_list, column_name, file_name_1, and file_name_2. Like the arguments in the previously mentioned function, the first argument refers to the filename of the species dataset that is read in. The second argument refers to the filename of the filtered gene list that is read in. This filtered gene list was created from the gene_filter variable in the previous function. The argument column_name refers to the column name that corresponds to the homologs. To get all the unique genes for the second species, a variable called list_1_column is created. This variable is a list that is made up of genes from the gene variable. The variable list_1_column is created and removes any duplicate gene names from list_1 and helps to make the overall list be as small as possible. The variable list_1_vector is created to convert the list_1_column to a vector. The variable query_1 is created and returns queries from the dataset that have a gene name that is present in list_1_vector. These queries are then saved to a csv file and given a name with the file_name_1 argument. For the filtered gene list to reflect the genes that are present for species_2, it must be queried against the queried dataset. To do this, all the unique gene ids for the second species must be retrieved. This is done through the creation of a variable called list_2, which creates a list of unique gene ids from the query_1 variable. The variable list_2_column is created, which removes any duplicate gene ids from list_2 and helps to make the list as small as possible. The variable list_2_vector is then created to convert the list_2_column to a vector. Since the name of the homolog can vary depending on the species that is chosen for analysis, the argument of column_name allows for different parameters to be entered for the function. The variable query_2 is created and returns queries from the genes variable that have a gene id that is present in list_2_vector based on the parameter entered for column_name. These queries are then saved to a csv file and given a name with the file_name_2 argument.

```R
dataset_1_final_filter<- function(species, gene_list, file_name) {
  dataset = read.csv(species)
  genes = read.csv(gene_list)
  list_1 = genes %>% select(external_gene_name)
  list_1_column = unique(list_1)
  list_1_vector = unlist(list_1_column)
  query_1 = dataset[dataset$external_gene_name %in% list_1_vector, ]
  write.csv(query_1, file_name, row.names=FALSE)
}  
```
The following chunk creates a function that updates the first dataset for the first species. As shown in the gene_list_dataset_2_filter function, the gene list was updated to remove any genes that were not present in the dataset for the second species. This process also removed the corresponding homologs for the first species. Due to the removal of the corresponding homologs, the dataset for the first species is updated for a final time to reflect the filtering that occurred for the dataset of the second species. This function has three arguments, which are species, gene_list, and file_name. The first argument refers to the filename of the filtered dataset that is read in for the first species. The second argument refers to the filename of the filtered gene list that is read in. This filtered gene list was created in the previous function and shows all the genes present for both datasets. The third argument refers to what the filtered dataset will be called. After reading in the csvs for the genes and dataset variables, a variable called list_1 is created. The variable creates a list of unique gene names from the dataset variable. The variable list_1_column is created and removes any duplicate gene names from list_1 and helps to make the overall list be as small as possible. The variable list_1_vector is created to convert the list_1_column to a vector. The variable query_1 is created and returns queries from the dataset that have a gene name that is present in list_1_vector. These queries are then saved to a csv file and given a name with the file_name argument.

```R
gene_ontology_filter <- function(file, go_term, go_name_filter) {
  filtered_species = read.csv(file)
  query = filtered_species[filtered_species$name_1006 == go_term,]
  write.csv(query, go_name_filter, row.names=FALSE)
}
```
The following chunk creates a function that returns queries that have a particular gene ontology term. The function has three arguments, which are file, go_term, and go_name filter. The first argument refers to the filtered species file that is read in and the second parameter refers to the selected gene ontology term. The third argument refers to the name that the output csv file is given. The filtered species file is read in and saved as filtered_species. The query variable is then created, which returns queries if the specified go_term is present in the GO_term_name column. These queries are then saved to a csv file and given a name with the go_name_filter argument.

```R
ref_seq_list <- function(file_name, column_name, gene_name, name) {
  file = read.csv(file_name)
  selected = file[, c(column_name, 'external_gene_name')]
  new_list = unique(selected)
  query = new_list[new_list$external_gene_name %in% c(gene_name),]
  write.csv(query, name, row.names=FALSE)
}
```
The following chunk creates a function that returns refseqs for a particular gene. The function has four arguments, which are file_name, column_name, gene, and name. The first argument refers to the selected file that contains queries filtered by gene ontology name for a specified species. The second argument column_name refers to the column in which the desired refseqs are. The refseqs can be selected from either the RefSeq_mRNA_ID or RefSeq_peptide_ID column. These columns contain the sequences for a specified gene in mRNA and protein form. The third argument gene refers to the gene of interest, and the fourth argument refers to what the output csv will be called. The selected variable is created by selecting one of the two previously mentioned columns and the column Gene_name. This creates a list that only has information related to the two chosen columns. The variable new_list is created, which drops any duplicate refseqs that may be present. The variable query is created which returns queries from new_list based on whether or not they match the gene_name argument. The query variable is then saved to a csv file and given a name with the name argument.

```R
ref_seq_sequence <- function(db_type, id, file_name) {
  net_handle <- entrez_fetch(db=db_type, id=id, rettype="fasta", retmode='text')
  write(net_handle, file = file_name)
}
```
The following chunk creates a function that outputs a refseq based on its id and database type. The function has three arguments, which are db_type, id, and file. The first parameter refers to the database type, which can either be nucleotide or protein. The second parameter id refers to the refseq id that is entered. The third parameter file_name refers to what the output file will be named. A variable called net_handle is created and retrieves information for a specified refseq depending on the type of database and refseq id that is entered. The net handle variable is then written to a specified file based on the name given with the file_name argument. It should be noted, however, that while rettype is set to fasta, it could be changed to 'gb' which is for GenBank.

```R
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
```
The following chunk creates a function that returns a pairwise alignment for two species along with its score. The function has six arguments, which are file_1, file_2, matrix, open_gap, extend_gap, snd file_name. The arguments of file_1 and file_2 refer to the fasta files that read in for species_1 and species_2. The argument matrix refers to the type of substitution matrix that will be used for the alignment. The fourth and fifth arguments refer to the penalties for the open and extended gaps. It is important to note that all values for these penalties must be entered as negative numbers. The sixth argument refers to what the written alignment will be named when it is made into a file. The variable species_1 refers to a fasta file that is read in for species_1. The variable species_1_character is created, which converts species_1 to a vector. The variable species_1_upper is created, which capitalizes all of the values for species_1_character. The variable species_1_unlist is created, which converts the species_1_upper to a vector again. The variable species_1_string is created, which converts species_1_unlist to a string. The variable species_1_comma is created, which removes the commas in the character string. The variable species_1_space is created, which removes spaces from the character string. Similar steps are taken to deal with the data for species_2 The variable species_2 refers to a fasta file that is read in for species_2. The variable species_2_character is created, which converts species_2 to a vector. The variable species_2_upper is created, which capitalizes all of the values for species_1_character. The variable species_2_unlist is created, which converts the species_2_upper to a vector again. The variable species_2_string is created, which converts species_2_unlist to a string. The variable species_2_comma is created, which removes the commas in the character string. The variable species_2_space is created, which removes spaces from the character string. The variable alignment is created, which creates the pairwise alignment. The pairwise alignment requires the sequences for species_1 and species_2, the type of alignment, the substitution matrix, the gap opening score, the extended gap score, and whether or not the score will be the only thing returned. Species_1 and species_2 are indicated by species_1_space and species_2_space. While the type of alignment is set to global, the other possibilities are local, overlap, global-local, and local-global. While the pre-defined substitution matrix is set to BLOSUM62, other matrices can be used. These matrices include BLOSUM45, BLOSUM50, BLOSUM80, BLOSUM100, PAM30, PAM40, PAM70, PAM120, and PAM250. In addition to this, it is also possible to make up a substitution matrix, that is not covered in this tutorial. Overall, the type of matrix is dependent on whether or not the sequences are in amino acid or nucleotide form. The function writePairwiseAlignments is then used to write the pairwise alignment to a file. This function requires the pairwise alignment, the name of the file, the matrix being used, and block width, which dermines how the alignment is arranged.

### Function Arguments

```R
mart_finder('mart_list_R.csv') 
```

```R
database_finder('ENSEMBL_MART_ENSEMBL', 'database_list_R.csv')
```

```R
filters_attributes('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', 'h_filter_R.csv', 'h_attrib_R.csv')
```

```R
filters_attributes('ENSEMBL_MART_ENSEMBL', 'mmusculus_gene_ensembl', 'm_filter_R.csv', 'm_attrib_R.csv')
```

```R
dataset_retrieve('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', '5', 'species_1_R.csv')
```

```R
dataset_retrieve('ENSEMBL_MART_ENSEMBL','mmusculus_gene_ensembl', '18', 'species_2_R.csv')
```

```R
gene_list('hsapiens_gene_ensembl', '5', 'mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_associated_gene_name', 'genes_R.csv')
```

```R
gene_list_dataset_1_filter('species_1_R.csv',  'genes_R.csv','species_1_filter_R.csv',
                           'filtered_gene_R.csv') 
```

```R
gene_list_dataset_2_filter('species_2_R.csv', 'filtered_gene_R.csv', 'mmusculus_homolog_ensembl_gene',  
                           'species_2_filter_final_R.csv', 'filtered_gene_final_R.csv')
```

```R
dataset_1_final_filter('species_1_filter_R.csv', 'filtered_gene_final_R.csv', 'species_1_filter_final_R.csv')
```

```R
gene_ontology_filter('species_1_filter_final_R.csv', 'plasma membrane', 'species_1_go_R.csv')
```

```R
gene_ontology_filter('species_2_filter_final_R.csv', 'plasma membrane', 'species_2_go_R.csv')
```

```R
ref_seq_list('species_1_go_R.csv', 'refseq_peptide', 'APC', 'species_1_ref_R.csv')
```

```R
ref_seq_list('species_2_go_R.csv', 'refseq_peptide', 'Apc', 'species_2_ref_R.csv')
```

```R
ref_seq_sequence('protein', 'NP_001394379', 'H_APC_ref_seq_R.fasta')
```

```R
ref_seq_sequence('protein', 'NP_001347909', 'M_Apc_ref_seq_R.fasta')
```

```R
pairwise_alignment('H_APC_ref_seq_R.fasta', 'M_Apc_ref_seq_R.fasta', 'BLOSUM62', -10, -0.5, 'alignment_R.txt')
```
The following chunks call all the previously mentioned functions. The order in which the functions are called and the order in which the arguments are entered into the functions is important. If the functions are called in the wrong order, this could cause the chunk to fail and possibly give inaccurate results as well as extra or unwanted files. The entered information for each function corresponds to a specified argument within the function. Entering wrong or extra information into the arguments will interfere with the pipeline process, as it will cause the functions to return an error message. Some of the functions are called more than once, since they deal with both species. Other functions deal with one species, so their functionality is more limited in comparison to the ones that can deal with both species.
