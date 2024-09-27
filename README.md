# Protease cleavage site specificity analysis - PICS

This is a small code for cleavage site specificity analysis

If you have any problems running the code you can check the following:
- pay attention to the conflicts between the packages used here and those you have appended.
- of course, check the names of your files **peptide.tsv** and **protein.fas**. Sometimes people change names.
- ⚠️when importing your **psm.tsv** files, make sure you are using the `clean_names()` function from `janitor` package.

### Following you find the packages and functions to help with the analysis

R packages you need to have installed and loaded:

```
library(tidyverse)
library(here)
library(janitor)
library(ComplexHeatmap)
```

Functions you need to follow the code in this repository:

```
# to calculate the frequency of each amino acid in each column in percentage
aa_freq <- function(x) {
    table(x) / length(x) * 100
}

# to impute missing amino acids with zero and reorder if necessary
complete_and_reorder_amino_acids <- function(element) {
# List of 20 amino acids
twenty_amino_acids <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

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

# FragPipe search results usage

The simplest way to spilt the extended_peptide sequence keeping the 4 amino acids before and after the first and second dot into two separate columns. The first will be _fingerprint_Nterm_ and the second _fingerprint_Cterm_. Once you imported the **psm.tsv** file you can create a fingerprint for the N-term and C-term with 4 amino acids on each side of the cleavage site. If you are using the latest version of **FragPipe**, you don't need to map the peptides to proteins in fasta file anymore, because you have the column extended peptide.

❗Remember that the **psm.tsv** file contains all the PSM scored. Only a fraction of these PSM will be transfered to the **peptide.tsv** file. So, you may want to apply a filter in _hyperscore_ column to select the best PSMs.

```
psm_file <- read_tsv("psm.tsv") %>%
    janitor::clean_names() %>%
    dplyr::mutate(
        fingerprint_Nterm = case_when(
            str_detect(extended_peptide, "^\\.") ~ "NA",
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
fingerprint_protease <- c(psm_file$fingerprint_Nterm,
            psm_file$fingerprint_Cterm) %>% 
        na.omit() %>%
        strsplit("")
 
# create a matrix with the list of peptide sequences
mat_aa <- matrix(unlist(fingerprint_protease),
                        ncol = 8, byrow = TRUE)

# remove the rows containing "B" or any other unwanted amino acids in the matrix
# mat_aa <- mat_aa[!apply(mat_aa, 1,
                    function(x) any(x == "B")), ]

# calculate the frequency of each amino acid in each column and plot a heatmap
mat_aa_freq <- apply(mat_aa, 2, aa_freq)

new_list <- lapply(mat_aa_freq, complete_and_reorder_amino_acids)

final_matrix <- matrix(unlist(new_list,), ncol = 8, byrow = FALSE)

colnames(final_matrix) <- c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
row.names(final_matrix) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
```

Plotting the result and save your heatmap

```
png("pics_plot_heatmap.png",
    width = 6, height = 6,
    units = "in", res = 300)
final_matrix %>%
    ComplexHeatmap::pheatmap(
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        fontsize = 10,
        show_rownames = TRUE,
        show_colnames = TRUE,
        color = colorRampPalette(c("grey95", "steelblue"))(20),
        main = "Cleavage Site specificity",
        gaps_col = 4,
        angle_col = "0"
    )
dev.off()
```
<p align="center">
<img src="https://github.com/41ison/Cleavage-site-specificity-analysis/blob/main/PICS_complexheatmap.png" width="500">
</p>

### Alternatively, you can plot using ggplot2

```
pics_plot <- final_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "residue") %>%
    pivot_longer(cols = -residue, names_to = "position",
                    values_to = "frequency") %>%
    dplyr::mutate(
        position = factor(position, c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")),
        residue = factor(residue, c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))
        ) %>%
    ggplot(aes(x = position, y = residue, fill = frequency)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "grey90", high = "dodgerblue4") +
    theme_void() +
    labs(
        title = "Cleavage Site Specificity",
        x = "Position",
        y = "Amino acid residue",
        fill = "Frequency (%)"
    ) +
    theme(text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(),
        legend.title = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"),
        legend.title.position = "top")

# save the plot
ggsave("PICS_plot_ggplot.png", 
    plot = pics_plot, 
    width = 6, height = 6, bg = "white",
    units = "in", dpi = 300)
```
<p align="center">
<img src="https://github.com/41ison/Cleavage-site-specificity-analysis/blob/main/PICS_ggplot2.png" width="500">
</p>

##Making your data more interesting and informative with seqLogos for N-term and C-term

Additional packages

```
library(ggseqlogo)    # to create the seqLogos
library(patchwork)    # to merge the figures in a nice panel
```

### Make the seqlogo for the N-term and C-term

```
Nterm_seqLogo_plot <- psm_file$fingerprint_Nterm %>%
    na.omit() %>%
    ggseqlogo::ggseqlogo(
  method = "bits",
  seq_type = "AA"
  ) +
  geom_hline(yintercept = 0, 
        color = "black", linetype = "dashed") +
    geom_vline(xintercept = 4.5, 
        color = "black", linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    text = element_text(size = 15, color = "black"),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "SeqLogo of the N-terminal fingerprint",
       x = "Amino acid position",
       y = "Bits")

ggsave("Nterm_seqLogo_plot.png", 
    plot = Nterm_seqLogo_plot, 
    width = 4, height = 5, bg = "white",
    units = "in", dpi = 300)

Cterm_seqLogo_plot <- psm_file$fingerprint_Cterm %>%
    na.omit() %>%
    ggseqlogo::ggseqlogo(
  method = "bits",
  seq_type = "AA"
  ) +
  geom_hline(yintercept = 0, 
        color = "black", linetype = "dashed") +
    geom_vline(xintercept = 4.5, 
        color = "black", linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    text = element_text(size = 15, color = "black"),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "SeqLogo of the C-terminal fingerprint",
       x = "Amino acid position",
       y = "Bits")
```

## Since you are here, you can use the `patchwork` package to merge the figures like this:

```
panel_plot <- (
    (
(Nterm_seqLogo_plot | Cterm_seqLogo_plot) + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
    ) / pics_plot) +
    plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(size = 40, face = "bold"))) +
    plot_layout(heights = c(1,4))

ggsave("panel_plot.png", 
        panel_plot, 
        width = 20, height = 25, 
        dpi = 300)
```

<p align="center">
<img src="https://github.com/41ison/Cleavage-site-specificity-analysis/blob/main/panel_plot.png" width="500">
</p>
