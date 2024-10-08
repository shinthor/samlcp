#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    # Check if the required arguments are provided
    if len(sys.argv) < 6:
        print("Usage: python create_bar_charts.py <input_file_path> <project_name> <column1_name> <column2_name> <output_dir>")
        return

    input_file_path = sys.argv[1]
    project_name = sys.argv[2]
    column1_name = sys.argv[3]
    column2_name = sys.argv[4]
    output_dir = sys.argv[5]

    y_axis_column_pos = 6
    x_axis_column_pos = 0
    x_axis_grouping_column_pos = 2

    # Read the input file into a DataFrame
    results_data = pd.read_csv(input_file_path, sep='\t')
    # Check if the input file has 7 columns. If 7, it has an extra identity column in position 3(0-indexed). If less, it does not 
    if len(results_data.columns) < 7:
        y_axis_column_pos = 5
    # Create a unique cell types list
    unique_cell_types = results_data[column2_name].unique()

    for cell_type in unique_cell_types:
        # Filter the data by cell type
        filtered_data = results_data[results_data[column2_name] == cell_type]
        sorted_hue_categories = np.sort(results_data[results_data.columns[x_axis_grouping_column_pos]].unique())
        # Order the data by age and column 3 category
        filtered_data = filtered_data.sort_values(by=[results_data.columns[x_axis_grouping_column_pos]])

        # Create a bar chart using Seaborn
        plt.figure(figsize=(12, 6))
        sns.barplot(x=results_data.columns[x_axis_column_pos],
                    y=results_data.columns[y_axis_column_pos],
                    hue=results_data.columns[x_axis_grouping_column_pos],
                    hue_order=sorted_hue_categories,
                    data=filtered_data)

        # Set plot title and labels
        plt.title(f"Cell Type: {cell_type}")
        plt.xlabel(results_data.columns[x_axis_column_pos])
        plt.ylabel("Proportion")
        plt.legend(title=results_data.columns[x_axis_grouping_column_pos], bbox_to_anchor=(1.005, 1.005), loc="upper left", borderaxespad=0)

        # Rotate x-axis labels
        plt.xticks(rotation=90)

        # Save the plot to a file
        output_file_path = f"{output_dir}/{project_name}_{cell_type}_{results_data.columns[2]}_bar_chart.svg"
        plt.savefig(output_file_path, dpi=300, bbox_inches='tight')

        # Close the plot to free up memory
        plt.close()

if __name__ == "__main__":
    main()
