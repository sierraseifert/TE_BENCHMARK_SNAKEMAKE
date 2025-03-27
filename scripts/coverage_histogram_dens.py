# this .py script is called from terminal using the following command (ensure
# you're in appropriate directory):
# python3 ./calculate_coverage.py path_to_reference_csv path_to_EG_input_csv path_to_RM2_input_csv path_to_EDTA_input_csv seq_len "Model" "Nested"

#Example: if we want to calculate percent coverage over nested elements in arabidopsis & plot the histogram
# python Downloads/calculate_coverage.py Desktop/araTha_minor_nest.csv 
# Desktop/araTha/araTha_wg_EG.csv Desktop/araTha/araTha_wg_RM2.csv Desktop/araTha/araTha_wg_EDTA.csv 
# 100000000 "Arabidopsis" "Internal Nested Elements" arabidopsis_internal_elements.pdf

#updated on February 27, 2025 from Sierra's original script!

import sys
import csv
import matplotlib.pyplot as plt
from textwrap import wrap
import seaborn as sns

def build_list(csv, seq_len):
    # returns a list that's the same length as the background genomic
    # sequence, st every element represents a position in the genomic
    # sequence. if the element is 1, then a TE covers this position.
    # if the element is 0, then a TE does not cover this position.
     
    list = [0] * seq_len

    with open(csv, 'r') as csv_fh:        
        for line in csv_fh:
            elem = line.rstrip().split(",")
            start = int(elem[2])
            end = int(elem[3])

            for i in range(start - 1, end):
                list[i] += 1

    return list

def calc_cov_bins(ref, query_csv, seq_len, test_pipeline):
    # Build the query list that indicates positions covered by TEs in the query file
    query_list = build_list(query_csv, seq_len)

    with open(ref, 'r') as csv_fh:  
        csv_reader = csv.reader(csv_fh)

        coverage_count = 0
        perfect_coverage_count = 0

        elem_coverage_list = []  # List to hold coverage percentages for the pipeline

        # Iterate over each element in the reference file
        for line in csv_reader:
            ref_start = int(line[2])  # Start position of the element
            ref_end = int(line[3])    # End position of the element

            # Calculate the perfect coverage for the element (total length)
            elem_perfect_coverage_count = ref_end - ref_start + 1
            perfect_coverage_count += elem_perfect_coverage_count

            # Count the positions covered by the query (TEs detected by the test pipeline)
            elem_coverage_count = 0
            for i in range(ref_start - 1, ref_end):
                if query_list[i] == 1:
                    elem_coverage_count += 1
        
            coverage_count += elem_coverage_count

            # Calculate the coverage percentage for this element
            coverage_p = elem_coverage_count / elem_perfect_coverage_count
            elem_coverage_list.append(coverage_p)  # Store the coverage percentage

        # Calculate the average coverage percentage for the test pipeline
        average_coverage_p = coverage_count / perfect_coverage_count

        # Store the list of coverage percentages for the pipeline in counts_dict
        counts_dict[test_pipeline] = elem_coverage_list  # Store actual coverage percentages

        # Print summary information
        print()
        print(f"AVERAGE COVERAGE BY {test_pipeline}: {average_coverage_p}")
        print(f"Meaning that on average, TEs in input reference csv are covered by TEs in input test csv "
              f"{average_coverage_p * 100:.2f}%.\n")
        print()
        print("COVERAGE BREAKDOWN:")
        # Bin the data into categories based on coverage percentage (but we're not using this for the KDE plot)
        c_0_20_list = []
        c_20_40_list = []
        c_40_60_list = []
        c_60_80_list = []
        c_80_100_list = []
        
        for c in elem_coverage_list:
            if c <= 0.2:
                c_0_20_list.append(c)
            elif c <= 0.4:
                c_20_40_list.append(c)
            elif c <= 0.6:
                c_40_60_list.append(c)
            elif c <= 0.8:
                c_60_80_list.append(c)
            elif c <= 1.0:
                c_80_100_list.append(c)

        # Print the number of elements in each coverage range (for information)
        print(f"{len(c_0_20_list)} TEs in input reference csv are 0-20% covered by {test_pipeline}")
        print(f"{len(c_20_40_list)} TEs in input reference csv are 20-40% covered by {test_pipeline}")
        print(f"{len(c_40_60_list)} TEs in input reference csv are 40-60% covered by {test_pipeline}")
        print(f"{len(c_60_80_list)} TEs in input reference csv are 60-80% covered by {test_pipeline}")
        print(f"{len(c_80_100_list)} TEs in input reference csv are 80-100% covered by {test_pipeline}")

    # No return needed, as counts_dict is a global variable and updated directly

    # P all results for all test pipelines
import seaborn as sns
import matplotlib.pyplot as plt

def plot_kde(counts_dict, model, nest_status, save_file):
    plt.figure(figsize=(8, 6))  # Set the figure size for better visibility

    # For each pipeline, use the actual coverage percentages for KDE plot
    sns.kdeplot(counts_dict['EDTA'], shade=True, color='red', label='EDTA', linewidth=2)
    sns.kdeplot(counts_dict['EarlGrey'], shade=True, color='green', label='EarlGrey', linewidth=2)
    sns.kdeplot(counts_dict['RepeatModeler2'], shade=True, color='blue', label='RepeatModeler2', linewidth=2)

    # Customize plot labels and title
    plt.xlabel('Coverage Percentage')  # Label the x-axis
    plt.ylabel('Density')              # Label the y-axis
    plot_sub = model
    plt.title('\n'.join(wrap(plot_sub, width=50)), size=12, style='italic')
    plt.suptitle(nest_status + ' Coverage Distribution by Pipeline for ')
    plt.legend(title='Pipeline', loc='upper right')  # Add a legend

    # Save the plot to a file
    plt.savefig(save_file)

if __name__ == "__main__":
    ref = sys.argv[1]
    query_csvs = sys.argv[2:5]
    seq_len = int(sys.argv[5])
    model = sys.argv[6]
    nest_status = sys.argv[7]
    save_file = sys.argv[8]
    
    test_pipelines = ["EarlGrey", "RepeatModeler2", "EDTA"]

    counts_dict = {
	    "EarlGrey": [0, 0, 0, 0, 0],
	    "RepeatModeler2": [0, 0, 0, 0, 0],
	    "EDTA": [0, 0, 0, 0, 0]
    }
	
for query_csv, test_pipeline in zip(query_csvs, test_pipelines):
	calc_cov_bins(ref, query_csv, seq_len, test_pipeline)

#plot_hist(counts_dict, model, nest_status, save_file)
plot_kde(counts_dict, model, nest_status, save_file)
