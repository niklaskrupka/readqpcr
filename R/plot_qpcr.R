plot_qpcr <- function(the_data,
                      norm){
  # Unquote the normalization factor, which is provided so it isn't plotted.
  norm <- enquo(norm)
  # Select Columns that contain calculated values, remove normalization gene and plot.
  the_data %>%
    select(Sample, Genotype, contains("2^-dCt"), -contains(quo_name(norm))) %>%
    gather(-Sample, -Genotype, key = "Gene", value = "Relative expression") %>%
    separate(col = Gene, into = c("Gene", "Method"), sep = "_") %>%
    ggplot(aes(x = Genotype, y = `Relative expression`, color = Genotype)) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x  = element_blank()) +
    scale_colour_manual(values = c("#333333", "#d55e00")) +
    stat_summary(geom = "bar", fun.y = "mean", fill = NA, width = 0.8) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.2) +
    geom_jitter(width = 0.2, size = 2) +
    facet_wrap(~Gene, scales = "free_y")
}
