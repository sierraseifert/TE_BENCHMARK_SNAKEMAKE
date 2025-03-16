# Load the config file
configfile: "config.yaml"

# Rule to ensure all CSV outputs and the statistics file are generated
rule all:
    input:
        config["output_dir"] + "/" + config["seq_name"] + "_Garlic.csv",
        config["output_dir"] + "/" + config["seq_name"] + "_RM2.csv",
        config["output_dir"] + "/" + config["seq_name"] + "_EG.csv",
        config["output_dir"] + "/" + config["seq_name"] + "_EDTA.csv",
        config["output_dir"] + "/" "stats.txt",
        config["output_dir"] + "/" + "plot_stats.png",
        config["output_dir"] + "/" + "plot_distribution.png"

# Rule for Garlic
rule garlic_to_csv:
    input:
        gff=config["Garlic_gff"] 
    output:
        csv="{output_dir}/{seq_name}_Garlic.csv" 
    script:
        "scripts/Garlic_gff_to_csv.R"  

# Rule for RM2
rule rm2_to_csv:
    input:
        gff=config["RM2_gff"]  
    output:
        csv="{output_dir}/{seq_name}_RM2.csv"  
    script:
        "scripts/RM2_gff_to_csv.R" 

# Rule for EG
rule eg_to_csv:
    input:
        gff=config["EG_gff"] 
    output:
        csv="{output_dir}/{seq_name}_EG.csv" 
    script:
        "scripts/EG_gff_to_csv.R" 

# Rule for EDTA
rule edta_to_csv:
    input:
        gff=config["EDTA_gff"] 
    output:
        csv="{output_dir}/{seq_name}_EDTA.csv"  
    script:
        "scripts/EDTA_gff_to_csv.R"  
        
# Rule to calculate statistics
rule calculate_statistics:
    input:
        garlic=config["output_dir"] + "/" + config["seq_name"] + "_Garlic.csv",
        rm2=config["output_dir"] + "/" + config["seq_name"] + "_RM2.csv",
        eg=config["output_dir"] + "/" + config["seq_name"] + "_EG.csv",
        edta=config["output_dir"] + "/" + config["seq_name"] + "_EDTA.csv",
    output:
        stats=config["output_dir"] + "/" + "stats.txt" 
    params:
        seq_len=config["seq_len"]
    shell:
        """
        python3 ./scripts/generate_stats.py {input.garlic} {input.rm2} {input.eg} {input.edta} {params.seq_len}
        """

rule generate_plot_stats:
    input:
        stats=config["output_dir"] + "/" + "stats.txt"
    output:
        plot=config["output_dir"] + "/" + "plot_stats.png"  
    shell:
        """
        jupyter nbconvert --to notebook --execute --inplace --allow-errors scripts/view_stats.ipynb --ExecutePreprocessor.timeout=600
        """

rule calculate_distribution:
    input:
        garlic=config["output_dir"] + "/" + config["seq_name"] + "_Garlic.csv",
        rm2=config["output_dir"] + "/" + config["seq_name"] + "_RM2.csv",
        eg=config["output_dir"] + "/" + config["seq_name"] + "_EG.csv",
        edta=config["output_dir"] + "/" + config["seq_name"] + "_EDTA.csv",
    output:
        stats=config["output_dir"] + "/" + "coverage.txt" 
    params:
        seq_len=config["seq_len"]
    shell:
        """
        python3 ./scripts/calculate_coverage.py {input.garlic} {input.rm2} {input.eg} {input.edta} {params.seq_len}
        """

rule generate_plot_distribution:
    input:
        stats=config["output_dir"] + "/" + "coverage.txt"
    output:
        plot=config["output_dir"] + "/" + "plot_distribution.png"  
    script:
        "scripts/view_distribution.R"  
