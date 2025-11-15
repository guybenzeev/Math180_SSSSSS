#Please run this function to install required packages for the project
    out_clean <- out |>
        filter(
        !grepl("-", ref_codon),  # no gaps anywhere in the codon
        !is.na(ref_aa)           # must translate to an amino acid
        )

  return(out_clean)