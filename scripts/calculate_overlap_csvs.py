import sys
import csv
import os

def build_list(csv, seq_len):
    list = [0] * seq_len

    with open(csv, 'r') as csv_fh:        
        for line in csv_fh:
            elem = line.rstrip().split(",")
            start = int(elem[2])
            end = int(elem[3])

            for i in range(start - 1, end):
                list[i] += 1

    return list


def calculate_overlap(input_csv, list, minor_nest_csv, major_nest_csv):
    minor_nest_rows = []
    major_nest_rows = []

    with open(input_csv, 'r') as csv_fh:  
        csv_reader = csv.reader(csv_fh)

        for line in csv_reader:
            # default no overlap
            overlap = "no"
            p_overlap = 0
            num_overlaps = 0
            
            start = int(line[2])
            end = int(line[3])

            for i in range(start - 1, end):
                if list[i] >= 2:
                    if overlap == "no":
                        overlap = "yes"
                    num_overlaps += 1
            
            if overlap == "yes":
                p_overlap = num_overlaps / (end - start + 1)

            # append overlap and p_overlap to row
            line.append(overlap)
            line.append(p_overlap)

            # Check which file the row should go into
            if p_overlap == 1.0:
                minor_nest_rows.append(line)
            elif p_overlap > 0 and p_overlap < 1.0:
                major_nest_rows.append(line)
        
    # Write rows with overlap == 1.0 to their own CSV
    with open(minor_nest_csv, 'w', newline='') as csv_fh:
        csv_writer = csv.writer(csv_fh)
        csv_writer.writerows(minor_nest_rows)

    # Write rows with positive but not 1.0 overlap to their own CSV
    with open(major_nest_csv, 'w', newline='') as csv_fh:
        csv_writer = csv.writer(csv_fh)
        csv_writer.writerows(major_nest_rows)                


def main():
    # Get command-line arguments: input CSV, sequence length, and output CSV filenames
    input_csv = sys.argv[1]
    seq_len = int(sys.argv[2])
    minor_nest = sys.argv[3]  # Output CSV for overlap == 1.0
    major_nest = sys.argv[4]  # Output CSV for positive but not 1.0 overlap

    # Build the list of overlaps from the input CSV
    overlap_list = build_list(input_csv, seq_len)

    # Call the function to calculate overlap and write to the specified CSV files
    calculate_overlap(input_csv, overlap_list, minor_nest, major_nest)


if __name__ == "__main__":
    main()
