# this .py script is called from terminal using the following command (ensure
# you're in appropriate directory):
# python3 ./calculate_coverage.py path_to_garlic_csv path_to_rm2_csv path_to_eg_csv path_to_edta_csv seq_len

import sys
import csv

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

def get_coverage(ref_csv, test_csv, seq_len):
    query_list = build_list(test_csv, seq_len)

    with open(ref_csv, 'r') as csv_fh:  

        csv_reader = csv.reader(csv_fh)

        coverage_count = 0
        perfect_coverage_count = 0

        elem_coverage_list = []

        # for every element in ref
        for line in csv_reader:
            
            ref_start = int(line[2])
            ref_end = int(line[3])

            # a perfect pipeline would cover the whole elem
            elem_perfect_coverage_count = ref_end - ref_start + 1
            perfect_coverage_count += elem_perfect_coverage_count

            # how many positions between and including ref_start and
            # ref_end are identified as a TE by the test pipeline?
            elem_coverage_count = 0
            for i in range(ref_start - 1, ref_end):
                if query_list[i] == 1:
                    elem_coverage_count += 1
        
            coverage_count += elem_coverage_count

            coverage_p = elem_coverage_count / elem_perfect_coverage_count
            elem_coverage_list.append(coverage_p)
    

        average_coverage_p = coverage_count / perfect_coverage_count

        # note that 0.2 coverage will go in c_0_20_list, 0.4 coverage will go in c_20_40_list, etc
        # ie if a coverage value is exactly equal to a cutoff value, it will be binned in the lower list
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
        
        num_0_20 = len(c_0_20_list)
        num_20_40 = len(c_20_40_list)
        num_40_60 = len(c_40_60_list)
        num_60_80 = len(c_60_80_list)
        num_80_100 = len(c_80_100_list)

        #############################################################        
        print()
        print("AVERAGE COVERAGE: " + str(average_coverage_p))
        print("Meaning that on average, TEs in input reference csv are covered by TEs in input test csv " 
              + str(average_coverage_p * 100) + "%.\n")
        print()
        print("COVERAGE BREAKDOWN:")
        print(str(num_0_20) + " TEs in input reference csv are 0-20" + "%" + " covered")
        print(str(num_20_40) + " TEs in input reference csv are 20-40" + "%" + " covered")
        print(str(num_40_60) + " TEs in input reference csv are 40-60" + "%" + " covered")
        print(str(num_60_80) + " TEs in input reference csv are 60-80" + "%" + " covered")
        print(str(num_80_100) + " TEs in input reference csv are 80-100" + "%" + " covered")

        return [num_0_20, num_20_40, num_40_60, num_60_80, num_80_100]

if __name__ == "__main__":
    garlic_csv = sys.argv[1]
    rm2_csv = sys.argv[2]
    eg_csv = sys.argv[3]
    edta_csv = sys.argv[4]
    seq_len = int(sys.argv[5])

    rm2_coverage = get_coverage(garlic_csv, rm2_csv, seq_len)
    eg_coverage = get_coverage(garlic_csv, eg_csv, seq_len)
    edta_coverage = get_coverage(garlic_csv, edta_csv, seq_len)

    # write coverage to a text file
    # st [# TEs 0-20% covered, # TEs 20-40% covered, # TEs 40-60% covered, # TEs 60-80% covered, # TEs 80-100% covered]
    output_file = "output/coverage.txt" 
    with open(output_file, "w") as f:
        f.write("RepeatModeler2:")
        f.write(", ".join(map(str, rm2_coverage)) + "\n")  

        f.write("Earl Grey:")
        f.write(", ".join(map(str, eg_coverage)) + "\n")

        f.write("EDTA:")
        f.write(", ".join(map(str, edta_coverage)) + "\n")

    print(f"Stats written to {output_file}")
