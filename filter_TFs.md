Script for filtering transcription factors in MaxQuant proteingroups
table
================

Load necessary libraries.

``` r
library(data.table)
library(xlsx)
library(biomaRt)
```

Select rat BioMart database for later use.

``` r
mart_rat <- useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl")
```

Generate list of transcription factors based on “DNA-binding
transcription factor activity (<GO:0003700>)” and “RNA polymerase II
cis-regulatory region sequence-specific DNA binding (<GO:0000978>)”
categories with mammalian filter from www.geneontology.org (data
downloaded 16.03.2020).

``` r
#Read in the Gene Ontology files.
TFs <- fread("20200316_GENEONTOLOGY_TRANSCRIPTION_FACTOR_MAMMALIAN.txt", header = FALSE)
RNAPII <- fread("20200316_GENEONTOLOGY_RNAPOLII_SEQUENCE-SPECIFIC_DNA_BINDING_MAMMALIAN.txt", header = FALSE)

#Convert gene symbols uppercase and keep only unique records.
combined_TFs <- unique(c(toupper(TFs$V2), toupper(RNAPII$V2)))
```

Main function for filtering rows containing transcription factors.

``` r
Filter_TFs <- function(input, sheetname) {
  
#Read the indicated sheet from the input xlsx file and convert it to data.table format.
dt <- as.data.table(read.xlsx(input,sheetName = sheetname))

#Split Protein IDs groups, making new row for each individual Protein ID.
dt_protein_IDs <- dt[,list(UniprotIDs_split = unlist(strsplit(as.character(Protein.IDs), ";"))), by=eval(colnames(dt))]
#Convert Uniprot IDs to gene symbols using BioMart.
UniprotIDs <- dt_protein_IDs$UniprotIDs_split
ID_to_name <- getBM(attributes = c('external_gene_name', "uniprot_gn_id"), filters = "uniprot_gn_id", values = UniprotIDs, mart=mart_rat)
dt_protein_IDs <- merge(dt_protein_IDs, ID_to_name, by.x ="UniprotIDs_split", by.y = "uniprot_gn_id")
#Find the values of Protein IDs column in rows where the converted gene name is contained in the transcription factor list.
split_IDs_TF <- dt_protein_IDs[toupper(external_gene_name) %in% combined_TFs]$Protein.IDs

#Split Gene names groups, making new row for each individual Gene name.
dt_gene_symbol <- dt[,list(Genes_split = unlist(strsplit(as.character(Gene.names), ";"))), by=eval(colnames(dt))]
#Find the values of Gene names column in rows where the split gene name is contained in the transcription factor list.
split_names_TF <- dt_gene_symbol[toupper(Genes_split) %in% combined_TFs]$Gene.names

#Based on the lists in split_IDs_TF and split_names_TF, keep only the rows from the original table that contain at least one transcription factor.
TF_list <- dt[Protein.IDs %in% split_IDs_TF | Gene.names %in% split_names_TF]

return(TF_list)
}
```

Apply the filtering function on the MaxQuant proteingroups tables.

``` r
Filtered_replicate1 <- Filter_TFs("Supplementary File 3.xlsx", "Replicate1")
Filtered_replicate2 <- Filter_TFs("Supplementary File 3.xlsx", "Replicate2")

#Write the filtered tables to files.
fwrite(Filtered_replicate1, "TFs_replicate1.csv")
fwrite(Filtered_replicate2, "TFs_replicate2.csv")
```

Information about the R session.

``` r
sessionInfo()
```

    ## R version 4.0.1 (2020-06-06)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 17763)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=Estonian_Estonia.1257  LC_CTYPE=Estonian_Estonia.1257   
    ## [3] LC_MONETARY=Estonian_Estonia.1257 LC_NUMERIC=C                     
    ## [5] LC_TIME=Estonian_Estonia.1257    
    ## system code page: 1252
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] biomaRt_2.44.4    xlsx_0.6.4.2      data.table_1.13.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] progress_1.2.2       tidyselect_1.1.0     xfun_0.16           
    ##  [4] purrr_0.3.4          rJava_0.9-13         vctrs_0.3.4         
    ##  [7] generics_0.1.0       htmltools_0.5.0      stats4_4.0.1        
    ## [10] BiocFileCache_1.12.1 yaml_2.2.1           blob_1.2.1          
    ## [13] XML_3.99-0.5         rlang_0.4.7          pillar_1.4.6        
    ## [16] glue_1.4.2           DBI_1.1.0            rappdirs_0.3.1      
    ## [19] BiocGenerics_0.34.0  bit64_4.0.5          dbplyr_1.4.4        
    ## [22] lifecycle_0.2.0      stringr_1.4.0        memoise_1.1.0       
    ## [25] evaluate_0.14        Biobase_2.48.0       knitr_1.30          
    ## [28] IRanges_2.22.2       curl_4.3             parallel_4.0.1      
    ## [31] AnnotationDbi_1.50.3 xlsxjars_0.6.1       Rcpp_1.0.5          
    ## [34] openssl_1.4.3        S4Vectors_0.26.1     bit_4.0.4           
    ## [37] hms_0.5.3            askpass_1.1          digest_0.6.25       
    ## [40] stringi_1.4.6        dplyr_1.0.2          tools_4.0.1         
    ## [43] magrittr_1.5         RSQLite_2.2.0        tibble_3.0.3        
    ## [46] crayon_1.3.4         pkgconfig_2.0.3      ellipsis_0.3.1      
    ## [49] xml2_1.3.2           prettyunits_1.1.1    assertthat_0.2.1    
    ## [52] rmarkdown_2.5        httr_1.4.2           R6_2.5.0            
    ## [55] compiler_4.0.1
