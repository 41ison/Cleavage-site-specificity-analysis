# Protease cleavage site specificity analysis - PICS

This is a small code for cleavage site specificity analysis

If you have any problems running the code you can check the following:
- pay attention to the conflicts between the packages used here and those you have appended.
- of course, check the names of your files **peptide.tsv** and **protein.fas**. Sometimes people change names.
- ⚠️when importing your **psm.tsv** files, make sure you are using the `clean_names()` function from `janitor` package.

### Following you find the functions to help with the analysis

```
# to calculate the frequency of each amino acid in each column in percentage
aa_freq <- function(x) {
    table(x) / length(x) * 100
}

# to impute missing amino acids with zero and reorder if necessary
complete_and_reorder_amino_acids <- function(element) {
# Complete missing amino acids with zero
  for (amino_acid in twenty_amino_acids) {
    if (!amino_acid %in% names(element)) {
      element[[amino_acid]] <- 0
    }
  }

# Reorder amino acids
  element <- element[match(twenty_amino_acids, names(element))]
  
  return(element)
}
```

The simplest way to spilt the extended_peptide sequence keeping the 4 amino acids before and after the first and second dot into two separate columns. The first will be fingerprint_Nterm and the second fingerprint_Cterm. Once you imported the psm.tsv file you can create a fingerprint for the N-term and C-term with 4 amino acids on each side of the cleavage site. If you are using the latest version of **FragPipe**, you don't need to map the peptides to proteins in fasta file anymore. The following code works well:

```
psm_file <- psm_tsv %>%
    dplyr::mutate(
        fingerprint_Nterm = case_when(
            stringr::str_detect(extended_peptide, "^\\.") ~ "NA",
            TRUE ~ substr(extended_peptide, 2, 16)
            ),
        fingerprint_Cterm = substr(extended_peptide, nchar(extended_peptide) - 15, nchar(extended_peptide) - 2),
        fingerprint_Nterm = str_extract(fingerprint_Nterm, ".{4}\\..{4}"),
        fingerprint_Nterm = str_remove_all(fingerprint_Nterm, "\\."),
        fingerprint_Cterm = str_extract(fingerprint_Cterm, ".{4}\\..{4}"),
        fingerprint_Cterm = str_remove_all(fingerprint_Cterm, "\\.")
        ) %>%
        relocate(extended_peptide, .before = fingerprint_Nterm)
```

Extract the matrix of amino acids and calculate the frequency of each residue in each position (from P4 to P4').

```
# create a list of peptide sequences including the N-term and C-term fingerprints
fingerprint_protease <- c(psm_file$fingerprint_Nterm, psm_file$fingerprint_Cterm) %>% 
        na.omit() %>%
        strsplit("")
 
# create a matrix with the list of peptide sequences
mat_aa <- matrix(unlist(fingerprint_protease), ncol = 8, byrow = TRUE)

# remove the rows containing "B" or any other unwanted amino acids in the matrix
# mat_aa <- mat_aa[!apply(mat_aa, 1, function(x) any(x == "B")), ]

# List of 20 amino acids
twenty_amino_acids <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

# calculate the frequency of each amino acid in each column and plot a heatmap
mat_aa_freq <- apply(mat_aa, 2, aa_freq)

new_list <- lapply(mat_aa_freq, complete_and_reorder_amino_acids)

final_matrix <- matrix(unlist(new_list,), ncol = 8, byrow = FALSE)

colnames(final_matrix) <- c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
row.names(final_matrix) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
```

Plotting the result and save your heatmap

```
PICS_plot <- final_matrix %>%
    ComplexHeatmap::pheatmap(
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        fontsize = 20,
        show_rownames = TRUE,
        show_colnames = TRUE,
        color = colorRampPalette(c("grey95", "steelblue"))(20),
        main = "Cleavage Site specificity",
        gaps_col = 4,
        angle_col = "0"
    )

ggsave("PICS_plot.png",
    path = "plots",
    PICS_plot, width = 8,
    height = 8, units = "in", dpi = 350)
```

You will find the Rscript for build a matrix of frequency for each amino acid per position.

![cleavage](https://github.com/41ison/Cleavage-site-specificity-analysis/assets/108031197/8b08e17d-29b1-4051-83c4-39c9e97cb7ce)
![seqlogo](https://github.com/41ison/Cleavage-site-specificity-analysis/assets/108031197/0832882a-a29f-41bc-a140-0ce9ee83d11c)
