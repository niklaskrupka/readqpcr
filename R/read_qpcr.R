#' @title Read Bio-Rad qPCR data
#' @description Reads an exported Bio-Rad qPCR result file and calculates quality
#'   controls summaries
#' @param datafiles A character vector of Bio-Rad files.
#' @return A tibble with collapsed by technical replicates and quality control indicators.
#' @examples
#' \dontrun{
#' read_qpcr("./raw_data/")
#' }
#' @import magrittr
#' @export
read_qpcr <- function(datafiles){
  # Turn off the readr warnings
  options(readr.num_columns = 0)
  # Define Biorad CSV headers
  bhead <- c("File", "Well", "Target", "Content", "Replicate", "Sample", "Cq")
  # Read the Bio-Rad CSVs. The first 19 rows in a Bio-Rad CSV export are the header, ignore them
  message("Loading the Bio-Rad CSV files.")
  qpcr_data <- lapply(datafiles, readr::read_csv, skip = 19)
  qpcr_data <- dplyr::bind_rows(qpcr_data, .id = "File")
  # Make sure it's the right format
  if(!all(names(qpcr_data) == bhead)){
    stop("The data files need to have the format: ", bhead)
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
  attr(qpcr_data_collapsed, "type") <- "read_qpcr"
  return(qpcr_data_collapsed)
}
