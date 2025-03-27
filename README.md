[![Snakemake](https://img.shields.io/badge/snakemake-≥3.13.3-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

# Welcome to the Snakemake workflow associated with the manuscript xxx.
The following is a quick guide on how to download and test the snakemake workflow. You can then adapt the config.yml file to suit your data!

## To Use on local computer or cluster via terminal with Conda installed
* We are unable to provide specific download instructions beyond those below. There are many resources available online if you run into computing environments beyond those described here.

### 1. Clone the repository to your local computer (simple download or git clone) or your cluser (git clone or wget)
* If you have git installed on your local computer or cluster (recommended)
> ```
> git clone https://github.com/sierraseifert/TE_BENCHMARK_SNAKEMAKE.git
> ```
* If you do not have git installed, you can perform a basic download on your local computer or cluster
  * Click the down arrow next to Code at the top of the repository and downlaod the zip file.
* If you do not have git insyalled, you can also use wget on your cluster if wget is installed
> ```
> wget https://github.com/sierraseifert/TE_BENCHMARK_SNAKEMAKE/archive/refs/heads/main.zip
> gunzip TE_BENCHMARK_SNAKEMAKE
> ```
### 2. Generate snakemake conda environment
> ```bash
> conda env create -f snakemake.yml
> ```
### 3. Remove the outputs from `/output/` directory
> ```bash
> rm output/*
> ```
### 4. Activate your snakemake environment
> ```bash
> conda activate snakemake
> ```
* We provide the expected outputs on our GitHub repository so you have a reference to compare your outputs to.
### 5. Check that the Snakefile rules are ready to go
> ```bash
> snakemake -n -r
> ```
* It should run with no errors, and you should see the following at the end of the terminal output:
> ```
> Job counts:
> count	jobs
> 1	all
> 1	calculate_statistics
> 1	edta_to_csv
> 1	eg_to_csv
> 1	garlic_to_csv
> 1	generate_plot_distribution_hostnest
> 1	generate_plot_distribution_minnest
> 1	generate_plot_distribution_overall
> 1	generate_plot_stats
> 1	nest_analysis
> 1	rm_to_csv
> 11
> ```
### 6. Then, run snakemake
> ```bash
> snakemake
> ```
* The job takes ~5 minutes to run without asking for additional compute resources.
* It should finsih all of the jobs, and you should see this at the end of the terminal output:
> ```
> Finished job x.
> 11 of 11 steps (100%) done.
> ```
* You should also see all 11 output files in the `/output` directory, which you can compare to those housed on this repository.

### Updates:
To come: an optional step to perform analyses on individual TE types (ie. LINE/SINE/LTR/etc) will be added in the near future!

### References:
This workflow was built using the Snakemake resources provided by

[Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.](https://f1000research.com/articles/10-33/v2)
