#' @title Read and process qPCR data
#' @description Reads an exported Bio-Rad qPCR result file, calculates quality
#'   controls summaries and performs 2^-ddCt calculations.
#' @param datafolder The directory, which contains one or multiple Biorad
#'   CSVs. Important: all CSV files in that location will be read and cause
#'   problems if they are not in the right format!
#' @param metadatafile A CSV file that contains sample metadata. At least a
#'   "Sample" and"Genotype" field have to be in that file
#' @param norm Reference gene to use
#' @return A list of two data frames: The first one contains the data
#'   collapsed by technical replicates and contains quality control indicators.
#'   The second one contains the 2^-ddCt calculations.
#' @examples
#' \dontrun{
#' read_qpcr("./raw_data/", metadatafile = "./metadata.csv", norm = Gapdh)
#' }
#' @import magrittr
#' @export
read_qpcr <- function(datafolder,
                     metadatafile = "./metadata.csv",
                     norm){
  # Turn off the readr warnings
  options(readr.num_columns = 0)
  # Define Biorad CSV headers
  bhead <- c("File", "Well", "Target", "Content", "Replicate", "Sample", "Cq")
  # Check if norm is set and unquote it
  if(missing(norm)){stop("Parameter norm is missing. You need to specify a housekeeping gene")}
  norm <- dplyr::enquo(norm)
  # Read the Bio-Rad CSVs. The first 19 rows in a Bio-Rad CSV export are the header, ignore them
  message("Loading the Bio-Rad CSV files.")
  qpcr_data <- lapply(list.files(path = datafolder, full.names = TRUE, pattern = "\\.csv$"),
                      readr::read_csv, skip = 19)
  qpcr_data <- dplyr::bind_rows(qpcr_data, .id = "File")
  # Make sure it's the right format and that there is actually the specified calibrator
  if(!all(names(qpcr_data) == bhead)){
    stop("The data files need to have the format: ", bhead)
    }
  if(!(dplyr::quo_name(norm) %in% qpcr_data$Target)){
    stop("The specified calibrator does not exist.")
  }
  # Read the Metadata and make sure it's the right format.
  # It has to contain at least the columns Sample and Genotype
  message("Loading the Metadata.")
  metadata <- readr::read_csv(metadatafile)
  if(!all(c("Sample", "Genotype") %in% names(metadata))){
    stop("The metadata file needs to contain at least a Sample and Genotype column.")
  }
  # Collapse technical replicates and calculate statistics.
  message("Collapsing technical replicates and calculating statistics.")
  qpcr_data_collapsed <- qpcr_data %>%
    dplyr::group_by(Target, Sample, Content) %>%
    dplyr::summarize(Cq_mean = mean(Cq),
                     Technical_replicates = n(),
                     Spread = max(Cq) - min(Cq),
                     SD = stats::sd(Cq)) %>%
    dplyr::ungroup()
  # Do the 2^-dCT calculations
  message("Perform 2^dCT calculations.")
  qpcr_data_calculated <- qpcr_data_collapsed %>%
    dplyr::select(Sample, Target, Cq_mean) %>%
    # Remove rows without a sample name (these are most likely NTCs)
    dplyr::filter(!is.na(Sample)) %>%
    tidyr::spread(key = "Target", value = "Cq_mean") %>%
    dplyr::mutate_if(is.double, dplyr::funs(dCt      = . - !!norm,
                                            `2^-dCt` = 2 ^(!!norm - .))) %>%
    dplyr::left_join(metadata, by = "Sample")
  # Return the two tibbles
  message("Finished.")
  return(list(qpcr_data_collapsed, qpcr_data_calculated))
}
