import sys
import math

def build_list(csv, seq_len):
    # P = 1. N = 0.

    list = [0] * seq_len

    with open(csv, 'r') as csv_fh:
        # check if file is empty
        lines = csv_fh.readlines()
        if not lines:
            return list
        
        # reset file pointer to the beginning
        csv_fh.seek(0)
        
        for line in csv_fh:
            elem = line.rstrip().split(",")
            start = int(elem[2])
            end = int(elem[3])

            for i in range(start - 1, end):
                list[i] = 1

    return list


def get_stats(ref_csv, test_csv, seq_len):
    ref_list = build_list(ref_csv, seq_len)
    test_list = build_list(test_csv, seq_len)

    TP = 0
    FP = 0
    FN = 0
    TN = 0

    for i in range(0, seq_len):
        # TP
        if ref_list[i] == 1 and test_list[i] == 1:
            TP += 1
        
        # FP
        elif ref_list[i] == 0 and test_list[i] == 1:
            FP += 1

        # FN
        elif ref_list[i] == 1 and test_list[i] == 0:
            FN += 1
        
        # TN
        else:
            TN += 1
    
    mcc = (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0 else 0
    sens = TP / (TP + FN)
    spec = TN / (FP + TN)
    accu = (TP + TN) / (TP + TN + FP + FN)
    # check if no positive predictions
    if (TP + FP) == 0:
        prec = 0
    else: 
        prec = TP / (TP + FP)
    fdr = 1 - prec
    f1 = (2 * TP) / ((2 * TP) + FP + FN)

    # for downstream use, generate stats s.t. ['mcc', 'sens', 'spec', 'accu', 'prec', 'FDR', 'F1']
    return([mcc, sens, spec, accu, prec, fdr, f1])

if __name__ == "__main__":
    garlic_csv = sys.argv[1]
    rm2_csv = sys.argv[2]
    eg_csv = sys.argv[3]
    edta_csv = sys.argv[4]
    seq_len = int(sys.argv[5])

    rm2_stats = get_stats(garlic_csv, rm2_csv, seq_len)
    eg_stats = get_stats(garlic_csv, eg_csv, seq_len)
    edta_stats = get_stats(garlic_csv, edta_csv, seq_len)

    # write stats to a text file
    output_file = "output/stats.txt" 
    with open(output_file, "w") as f:
        f.write("stats RM2:\n")
        f.write(", ".join(map(str, rm2_stats)) + "\n\n")  

        f.write("stats EG:\n")
        f.write(", ".join(map(str, eg_stats)) + "\n\n")

        f.write("stats EDTA:\n")
        f.write(", ".join(map(str, edta_stats)) + "\n\n")

    print(f"Stats written to {output_file}")