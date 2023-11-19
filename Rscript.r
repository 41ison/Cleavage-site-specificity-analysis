## Protease cleavage site specificity analysis - PICS
# If you have any problem running the code you can check the following:
# 1 - pay attention to the conflicts between the packages used here and those ones you have appended
# 2 - run the sessionInfo() to see the versions of each package used
# 3 - of course, check the names of your files peptide.tsv and protein.fas. Sometimes people change names

# output of the sessionInfo() in my workflow
# other attached packages:
# [1] ggseqlogo_0.1   janitor_2.2.0   pheatmap_1.0.12 lubridate_1.9.3
# [5] forcats_1.0.0   stringr_1.5.0   dplyr_1.1.3     purrr_1.0.2    
# [9] readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.4  
#[13] tidyverse_2.0.0 seqinr_4.2-30 

library(seqinr)
library(tidyverse)
library(pheatmap)
library(janitor)
library(ggseqlogo)

# load the peptide file from FragPipe output as peptides_df
peptides_df <- read_tsv("peptide.tsv") %>%
    clean_names() %>%
    rename(name = protein)

# read the fasta file and store sequences
fasta_sequences <- read.fasta("protein.fas",
    seqtype = "AA",
    as.string = TRUE
)

# extract names from fasta file
fasta_names <- sapply(fasta_sequences, attr, "name")

# function to match names and extract sequences
matched_sequences <- lapply(peptides_df$name, function(name) {
    idx <- match(name, fasta_names)
    if (!is.na(idx)) {
        return(as.character(fasta_sequences[[idx]]))
    } else {
        return("Name not found in fasta file")
    }
})

# add matched sequences to the data frame
peptides_df$protein_sequence <- matched_sequences

# match peptide sequence, extract previous and next 4 amino acids, and add to new columns
peptides_df <- peptides_df %>%
    mutate(previous_4_aa = sapply(1:nrow(.), function(row_idx) {
        peptide <- peptides_df$peptide[row_idx]
        protein <- peptides_df$protein_sequence[row_idx]
        match_idx <- regexpr(peptide, protein)
        if (match_idx > 0) {
            start_idx <- max(match_idx - 4, 1)
            previous_aa <- substr(protein, start_idx, match_idx - 1)
            return(previous_aa)
        } else {
            return("Peptide not found in protein")
        }
    }),
    next_4_aa = sapply(1:nrow(.), function(row_idx) {
        peptide <- peptides_df$peptide[row_idx]
        protein <- peptides_df$protein_sequence[row_idx]
        match_idx <- regexpr(peptide, protein)
        if (match_idx > 0) {
            end_idx <- match_idx + attr(match_idx, "match.length")
            next_aa <- substr(protein, end_idx, end_idx + 3)
            return(next_aa)
        } else {
            return("Peptide not found in protein")
        }
    }))

# Add new columns with the previous and the next 4 amino acids + the peptide fragment
# hint 1 - if you did your search in MaxQuant you already have this column in your peptide file
# hint 2 - if you are using the latest version of FragPipe and selecting your peptide sequences from psm.tsv
# file instead of peptide.tsv file, you already have this column as extended peptide sequence
peptides_df <- peptides_df %>%
    dplyr::mutate(
        fingerprint_Nterm = paste(previous_4_aa, peptide, sep = ""),
            fingerprint_Cterm = paste(peptide, next_4_aa, sep = "")
            ) %>%
    filter(
        !grepl("Peptide not found in protein", fingerprint_Nterm),
        previous_4_aa != "M"
    ) %>%
    filter(nchar(fingerprint_Nterm) >= 8) %>% # remove peptides with less than 8 amino acids in the N-term
    dplyr::mutate(
        fingerprint_Nterm = substr(fingerprint_Nterm, 1, 8),
        fingerprint_Cterm = substr(fingerprint_Cterm, nchar(fingerprint_Cterm) - 7, nchar(fingerprint_Cterm)),
        fingerprint_Cterm = case_when(
            nchar(next_4_aa) < 4 ~ NA_character_,
            TRUE ~ fingerprint_Cterm
        )
        )
# if the column next4_aa is of lenght 0, it means that the peptide is at the end of the protein and the next_4_aa are not available.
# So the peptide is removed from the fingerprint_Cterm column

view(peptides_df)

# at this point you should be able to plot the sequence logo of the fingerprint in the N-term
ggseqlogo(
    peptides_df$fingerprint_Nterm,
    method = "bits",
    seq_type = "AA") +
    theme_bw() +
    geom_vline(xintercept = 4.5, color = "black", linetype = "dashed") +
    labs(
        x = "Position",
        y = "Frequency (bits)",
        title = "Sequence logo of the N-terminal fingerprint"
    ) +
    theme(
        text = element_text(size = 25)
    )
 
# and the sequence logo of the fingerprint in the C-term
ggseqlogo(
    peptides_df$fingerprint_Cterm %>% na.omit(),
    method = "bits",
    seq_type = "AA",
    na_col = "grey50") +
    theme_bw() +
    geom_vline(xintercept = 4.5, color = "black", linetype = "dashed") +
    labs(
        x = "Position",
        y = "Frequency (bits)",
        title = "Sequence logo of the C-terminal fingerprint"
    ) +
    theme(
        text = element_text(size = 25)
    )

# create a list of peptide sequences including the N-term and C-term fingerprints
fingerprint_protease <- c(peptides_df$fingerprint_Nterm, peptides_df$fingerprint_Cterm) %>% 
        na.omit() %>%
        strsplit("")
 
# create a matrix with the list of peptide sequences
mat_aa <- matrix(unlist(fingerprint_protease), ncol = 8, byrow = TRUE)

# remove the rows containing "B" or any other unwanted amino acids in the matrix
# mat_aa <- mat_aa[!apply(mat_aa, 1, function(x) any(x == "B")), ]

# function to calculate the frequency of each amino acid in each column in percentage
aa_freq <- function(x) {
    table(x) / length(x) * 100
}

# calculate the frequency of each amino acid in each column and plot a heatmap
mat_aa_freq <- apply(mat_aa, 2, aa_freq)

# List of all 20 amino acids
# twenty_amino_acids <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

# Function to complete missing amino acids with zero and reorder if necessary
# you can skip this step if you have a complete matrix
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

# If you used the function above, you want to apply it to each column of the matrix
# new_list <- lapply(mat_aa_freq, complete_and_reorder_amino_acids)

# If you didn't use the function above, you can jump to this step
# final_matrix <- matrix(unlist(new_list,), ncol = 8, byrow = FALSE)

colnames(mat_aa_freq) <- c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
# row.names(final_matrix) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

PICS_map <- mat_aa_freq %>%
    pheatmap(
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        fontsize = 20,
        show_rownames = TRUE,
        show_colnames = TRUE,
        color = colorRampPalette(c("grey95", "#A73030FF"))(20),
        main = "Cleavage Site specificity",
        gaps_col = 4,
        angle_col = 0
    )

ggsave("PICS_map.png",
        PICS_map,
        width = 10, height = 12, 
        units = "in", dpi = 300)
