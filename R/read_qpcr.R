read_qpr <- function(datafolder,
                     metadatafile = "./metadata.csv",
                     norm){
  # Turn off the readr warnings
  options(readr.num_columns = 0)
  # Define Biorad CSV headers
  bhead <- c("File", "Well", "Target", "Content", "Replicate", "Sample", "Cq")
  # Check if norm is set and unquote it
  if(missing(norm)){stop("Parameter norm is missing. You need to specify a housekeeping gene")}
  norm <- enquo(norm)
  # Read the Bio-Rad CSVs. The first 19 rows in a Bio-Rad CSV export are the header, ignore them
  message("Loading the Bio-Rad CSV files.")
  qpcr_data <- lapply(list.files(path = datafolder, full.names = TRUE, pattern = "\\.csv$"), read_csv, skip = 19)
  qpcr_data <- bind_rows(qpcr_data, .id = "File")
  # Make sure it's the right format and that there is actually the specified calibrator
  if(!all(names(qpcr_data) == bhead)){
    stop("The data files need to have the format: ", bhead)
    }
  if(!(quo_name(norm) %in% qpcr_data$Target)){
    stop("The specified calibrator does not exist.")
  }
  # Read the Metadata and make sure it's the right format.
  # It has to contain at least the columns Sample and Genotype
  message("Loading the Metadata.")
  metadata <- read_csv(metadatafile)
  if(!all(c("Sample", "Genotype") %in% names(metadata))){
    stop("The metadata file needs to contain at least a Sample and Genotype column.")
  }
  # Collapse technical replicates and calculate statistics.
  message("Collapsing technical replicates and calculating statistics.")
  qpcr_data_collapsed <- qpcr_data %>%
    group_by(Target, Sample, Content) %>%
    summarize(Cq_mean = mean(Cq),
              Technical_replicates = n(),
              Spread = max(Cq) - min(Cq),
              SD = sd(Cq)) %>%
    ungroup()
  # Do the 2^-dCT calculations
  message("Perform 2^dCT calculations.")
  qpcr_data_calculated <- qpcr_data_collapsed %>%
    select(Sample, Target, Cq_mean) %>%
    # Remove rows without a sample name (these are most likely NTCs)
    filter(!is.na(Sample)) %>%
    spread(key = "Target", value = "Cq_mean") %>%
    mutate_if(is.double, funs(dCt      = . - !!norm,
                              `2^-dCt` = 2 ^(!!norm - .))) %>%
    left_join(metadata, by = "Sample")
  # Return the two tibbles
  message("Finished.")
  return(list(qpcr_data_collapsed, qpcr_data_calculated))
}
