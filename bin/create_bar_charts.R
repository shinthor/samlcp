library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
print(args)
input_file_path <- args[1]
project_name <- args[2]
output_dir <- args[5]

results_data <- read.table(input_file_path, header = TRUE, sep = "\t")

# Create unique cell types list
unique_cell_types <- unique(results_data[[args[4]]])

for (cell_type in unique_cell_types) {
  filtered_data <- subset(results_data, results_data[[args[4]]] == cell_type)
  
  # Order the data by age and the column 3 category
  filtered_data <- filtered_data[order(filtered_data[[args[3]]], filtered_data[[colnames(filtered_data)[1]]]), ]
  
  p <- ggplot(filtered_data, aes(x = factor(.data[[colnames(filtered_data)[3]]]), y = .data[[colnames(filtered_data)[6]]], fill = factor(.data[[colnames(filtered_data)[1]]]))) +
    geom_bar(position = "dodge", stat = "identity") +
    labs(fill = colnames(filtered_data)[1], y = "Proportion", x = colnames(filtered_data)[3]) +
    # theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle(paste0("Cell Type: ", cell_type))
  column3_name <- colnames(results_data)[3]
  output_file_path <- file.path(output_dir, paste0(project_name, cell_type, "_",  column3_name, "_bar_chart.pdf"))
  ggsave(output_file_path, plot = p, width = 12, height = 6, dpi = 300, device = "pdf")
}