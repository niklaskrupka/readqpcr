#' @title Plot normalized qPCR data in a nice format
#' @description Takes normalized qPCR data and turns it into a ggplot2 object for visualization
#' @param the_data A data frame containing normalized qPCR data. If read_qpcr is used this would be the second
#'   element of the returned list. Make sure the data frame contains at least the columns "Sample" and "Genotype".
#' @param norm Reference gene to use
#' @return A ggplot2 object
#' @examples
#' \dontrun{
#' plot_qpcr(my_pcrdata[[2]], norm = Gapdh)
#' }
#' @import magrittr
#' @export
plot_qpcr <- function(the_data,
                      norm){
  # Unquote the normalization factor, which is provided so it isn't plotted.
  norm <- dplyr::enquo(norm)
  # Select Columns that contain calculated values, remove normalization gene and plot.
  the_data %>%
    dplyr::select(Sample, Genotype,
                  dplyr::contains("2^-dCt"),
                  -dplyr::contains(dplyr::quo_name(norm))) %>%
    tidyr::gather(-Sample, -Genotype, key = "Gene", value = "Relative expression") %>%
    tidyr::separate(col = Gene, into = c("Gene", "Method"), sep = "_") %>%
    ggplot2::ggplot(ggplot2::aes(x = Genotype, y = `Relative expression`, color = Genotype)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x  = ggplot2::element_blank()) +
    ggplot2::scale_colour_manual(values = c("#333333", "#d55e00")) +
    ggplot2::stat_summary(geom = "bar", fun.y = "mean", fill = NA, width = 0.8) +
    ggplot2::stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.2) +
    ggplot2::geom_jitter(width = 0.2, size = 2) +
    ggplot2::facet_wrap(~Gene, scales = "free_y")
}
