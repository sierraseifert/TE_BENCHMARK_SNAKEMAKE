# Load the config file
configfile: "config.yaml"

# Rule to ensure all CSV outputs and the statistics file are generated
rule all:
    input:
        config["output_dir"] + "/" "stats.txt",
        config["output_dir"] + "/" + "plot_stats.pdf",
        config["output_dir"] + "/" + "plot_distribution.pdf",
        expand(config["output_dir"] + "/" + config["seq_name"] + "_{prog}.csv",prog=config["programs_short"]),
        config["output_dir"] + "/" + config["seq_name"] + "_minnest.csv",
        config["output_dir"] + "/" + config["seq_name"] + "_hostnest.csv",
        config["output_dir"] + "/" + "plot_distribution_minnest.pdf",
        config["output_dir"] + "/" + "plot_distribution_hostnest.pdf"

rule garlic_to_csv:
    input:
        gff=config["Garlic_gff"]
    output:
        csv=config["output_dir"] + "/" + config["seq_name"] + "_Garlic.csv"
    script:
        "scripts/Garlic_gff_to_csv.R"

rule rm_to_csv:
    input:
        gff=config["RM2_gff"]
    output:
        csv=config["output_dir"] + "/" + config["seq_name"] + "_RM2.csv"
    script:
        "scripts/RM2_gff_to_csv.R"

rule eg_to_csv:
    input:
        gff=config["EG_gff"]
    output:
        csv=config["output_dir"] + "/" + config["seq_name"] + "_EG.csv"
    script:
        "scripts/EG_gff_to_csv.R"

rule edta_to_csv:
    input:
        gff=config["EDTA_gff"]
    output:
        csv=config["output_dir"] + "/" + config["seq_name"] + "_EDTA.csv"
    script:
        "scripts/EDTA_gff_to_csv.R"

rule nest_analysis:
    input:
        "{output_dir}/{seq_name}_Garlic.csv"
    output:
        csv1="{output_dir}/{seq_name}_minnest.csv",
        csv2="{output_dir}/{seq_name}_hostnest.csv"
    params:
        seq_len=config["seq_len"]
    shell:
        """
        python3 ./scripts/calculate_overlap_csvs.py {input} {params.seq_len} {output.csv1} {output.csv2}
        """

rule calculate_statistics:
    input:
        garlic=config["output_dir"] + "/" + config["seq_name"] + "_Garlic.csv",
        rm2=config["output_dir"] + "/" + config["seq_name"] + "_RM2.csv",
        eg=config["output_dir"] + "/" + config["seq_name"] + "_EG.csv",
        edta=config["output_dir"] + "/" + config["seq_name"] + "_EDTA.csv"
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
	    plot=config["output_dir"] + "/" + "plot_stats.pdf"
    params:
        model=config["Model"],
        coverage="Whole_Genome"
    shell:
	    """
	    python3 ./scripts/view_stats.py {input.stats} {output.plot} {params.model} {params.coverage}
        """

rule generate_plot_distribution_overall:
    input:
        garlic=config["output_dir"] + "/" + config["seq_name"] + "_Garlic.csv",
        rm2=config["output_dir"] + "/" + config["seq_name"] + "_RM2.csv",
        eg=config["output_dir"] + "/" + config["seq_name"] + "_EG.csv",
        edta=config["output_dir"] + "/" + config["seq_name"] + "_EDTA.csv"
    output:
        plot=config["output_dir"] + "/" + "plot_distribution.pdf"
    params:
        model=config["Model"],
        seq_len=config["seq_len"],
        nest_status="All Elements"
    shell:
        """
        python3 ./scripts/coverage_histogram_dens.py {input.garlic} {input.eg} {input.rm2} {input.edta} {params.seq_len} {params.model:q} {params.nest_status:q} {output.plot}
        """

rule generate_plot_distribution_minnest:
    input:
        garlic=config["output_dir"] + "/" + config["seq_name"] + "_minnest.csv",
        rm2=config["output_dir"] + "/" + config["seq_name"] + "_RM2.csv",
        eg=config["output_dir"] + "/" + config["seq_name"] + "_EG.csv",
        edta=config["output_dir"] + "/" + config["seq_name"] + "_EDTA.csv"
    output:
        plot=config["output_dir"] + "/" + "plot_distribution_minnest.pdf"
    params:
        model=config["Model"],
        seq_len=config["seq_len"],
        nest_status="Nested TEs"
    shell:
        """
        python3 ./scripts/coverage_histogram_dens.py {input.garlic} {input.eg} {input.rm2} {input.edta} {params.seq_len} {params.model:q} {params.nest_status:q} {output.plot}
        """

rule generate_plot_distribution_hostnest:
    input:
        garlic=config["output_dir"] + "/" + config["seq_name"] + "_hostnest.csv",
        rm2=config["output_dir"] + "/" + config["seq_name"] + "_RM2.csv",
        eg=config["output_dir"] + "/" + config["seq_name"] + "_EG.csv",
        edta=config["output_dir"] + "/" + config["seq_name"] + "_EDTA.csv"
    output:
        plot=config["output_dir"] + "/" + "plot_distribution_hostnest.pdf"
    params:
        model=config["Model"],
        seq_len=config["seq_len"],
        nest_status="Host TEs"
    shell:
        """
        python3 ./scripts/coverage_histogram_dens.py {input.garlic} {input.eg} {input.rm2} {input.edta} {params.seq_len} {params.model:q} {params.nest_status:q} {output.plot}
        """
