## Protease cleavage site specificity analysis - PICS
This is a small code for cleavage site specificity analysis

If you have any problems running the code you can check the following:
1 - pay attention to the conflicts between the packages used here and those you have appended.
2 - run the `sessionInfo()` to see the versions of each package used.
3 - of course, check the names of your files `peptide.tsv` and `protein.fas`. Sometimes people change names.

The simplest way to spilt the extended_peptide sequence keeping the 4 amino acids before and after the first and second dot into two separate columns. The first will be fingerprint_Nterm and the second fingerprint_Cterm.

```
figerprint_df <- peptide_file %>%
    dplyr::mutate(
        extended_sequence = stringr::str_replace_all(extended_peptide, "\\.", ""),
        fingerprint_Nterm = substr(extended_sequence, 5, 12),
        fingerprint_Cterm = substr(extended_sequence, nchar(extended_sequence) - 11, nchar(extended_sequence) - 4)) %>%
        dplyr::select(peptide, fingerprint_Nterm, fingerprint_Cterm) %>%
        distinct()
```

You will find the Rscript for build a matrix of frequency for each amino acid per position.
![cleavage](https://github.com/41ison/Cleavage-site-specificity-analysis/assets/108031197/8b08e17d-29b1-4051-83c4-39c9e97cb7ce)
![seqlogo](https://github.com/41ison/Cleavage-site-specificity-analysis/assets/108031197/0832882a-a29f-41bc-a140-0ce9ee83d11c)
