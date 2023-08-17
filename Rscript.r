# Cleavage site specificity analysis

library(seqinr)
library(tidyverse)
library(pheatmap)
library(janitor)
library(ggseqlogo)

# assuming your data frame is named 'key_df' and the key column is named 'key_column'
key_df <- read_tsv("/Users/chaves/Documents/Kidney_TMT16plex/Fixed_TMT_3h/peptide.tsv") %>%
    clean_names() %>%
    rename(name = protein)

# read the fasta file and store sequences
fasta_sequences <- read.fasta("/Users/chaves/Documents/Kidney_TMT16plex/Fixed_TMT_3h/protein.fas",
    seqtype = "AA",
    as.string = TRUE
)

# extract names from fasta file
fasta_names <- sapply(fasta_sequences, attr, "name")

# function to match names and extract sequences
matched_sequences <- lapply(key_df$name, function(name) {
    idx <- match(name, fasta_names)
    if (!is.na(idx)) {
        return(as.character(fasta_sequences[[idx]]))
    } else {
        return("Name not found in fasta file")
    }
})

# add matched sequences to the data frame
key_df$protein_sequence <- matched_sequences

# match peptide sequence, extract previous 4 amino acids, and add to new column
key_df <- key_df %>%
    mutate(previous_4_aa = sapply(1:nrow(.), function(row_idx) {
        peptide <- key_df$peptide[row_idx]
        protein <- key_df$protein_sequence[row_idx]
        match_idx <- regexpr(peptide, protein)
        if (match_idx > 0) {
            start_idx <- max(match_idx - 4, 1)
            previous_aa <- substr(protein, start_idx, match_idx - 1)
            return(previous_aa)
        } else {
            return("Peptide not found in protein")
        }
    }))

# mutate a new column with the previous 4 amino acids and the peptide
key_df <- key_df %>%
    mutate(fingerprint_protease = paste(previous_4_aa, peptide, sep = "")) %>%
    filter(
        !grepl("Peptide not found in protein", fingerprint_protease),
        previous_4_aa != "M"
    )

# remove the rows in which the previous 4 amino acids are less than 4
key_df <- key_df %>%
    filter(nchar(previous_4_aa) == 4)

# keep the first 8 amino acids of the fingerprint
key_df$fingerprint_protease <- substr(key_df$fingerprint_protease, 1, 8)

view(key_df)

# at this point you should be able to plot the sequence logo of the fingerprint
ggseqlogo(
    key_df$fingerprint_protease,
    method = "bits",
    seq_type = "AA") +
    theme_bw()

# create a list of peptide sequences
fingerprint_list <- strsplit(key_df$fingerprint_protease, "")
head(fingerprint_list)
 
 # create a matrix with the list of peptide sequences
mat_aa <- matrix(unlist(fingerprint_list), ncol = 8, byrow = TRUE)

# remove the rows containing "B" or any other unwanted amino acids in the matrix
mat_aa <- mat_aa[!apply(mat_aa, 1, function(x) any(x == "B")), ]

# calculate the frequency of each amino acid in each column and plot a heatmap
mat_aa_freq <- apply(mat_aa, 2, function(x) table(x) / length(x)) * 100
colnames(mat_aa_freq) <- c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")

# if for some weird reason you have an uneven number of amino acids in the matrix, you can use this code to add a row of zeros
# mat_aa_freq <- rbind(mat_aa_freq, rep(0, ncol(mat_aa_freq)))

mat_aa_freq %>%
    pheatmap(
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        fontsize = 20,
        show_rownames = TRUE,
        show_colnames = TRUE,
        color = colorRampPalette(c("grey95", "#A73030FF"))(20),
        main = "Cleavage Site specificity",
        gaps_col = 4
    )
