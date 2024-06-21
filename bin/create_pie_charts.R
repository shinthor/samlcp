#!/usr/bin/env Rscript

library(ggplot2)
require(stringr)
args <- commandArgs(trailingOnly=TRUE)
input_file_path <- args[1]
project_name <- args[2]
output_dir <- args[5]

results_data <- read.table(input_file_path, header = TRUE, sep = "\t")
results_data[5] <- log((results_data[5] * 100 + 1))

ggplot(results_data,
       aes(x = .data[[colnames(results_data)[5]]] / 2,
           y = .data[[colnames(results_data)[6]]],
           fill = .data[[colnames(results_data)[3]]],
           width = .data[[colnames(results_data)[5]]])) +
  labs(fill = colnames(results_data)[3],
       y = colnames(results_data)[1],
       x = colnames(results_data)[2]) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(vars(.data[[colnames(results_data)[2]]]), vars(.data[[colnames(results_data)[1]]])) +
  coord_polar("y") +
  theme(axis.text.x = element_blank())

column3_name <- colnames(results_data)[3]
output_dir <- file.path(output_dir)
output_file_path <- file.path(output_dir, paste0(project_name, "_", column3_name, "_pie.png"))
ggsave(output_file_path, width = 12, height = 12, dpi = 300, create.dir = TRUE)
ggsave(file.path(".", paste0(project_name, "_", column3_name, "_pie.png")), width = 20, height = 20, dpi = 300, create.dir = TRUE)