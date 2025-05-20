# Visualize differential gene module expression.

Requires samples to be pre-processed using Kallisto pseudo-aligner.

## Preprocessing (before running ImmunomeAnalyzeR)

To create the Kallisto index:
```
kallisto index -i [outputIndexName] -t [# of threads] [reference_cDNA_fasta_file]
```

To run Kallisto quantification:
```
kallisto quant -i [index] -o [output] pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq
```

## Running ImmunomeAnalyzeR

### Description
Compares two groups (ex. treatment vs. control, healthy vs. disease, etc.)
by comparing expression of groups of related genes.

### Details
The config.yml file must be configured prior to running the script.

* 'study_design_file': Path to the tab-delimited text file that contains two
   columns: a sample identifier (ID) and a group. See
   'study_design_template.txt'.

* 'treatment_name': The name of the group of interest in the study design
   file. Usually this is 'patient'. The script will generate a plot for each
   treatment sample.

* 'data_dir': Path to the directory with Kallisto results. The directory must
  be structured as the raw Kallisto output: one sub-directory for each
  patient ID with an 'abundance.h5' file inside:
  data_dir
  |---SampleID1
  | |---abundance.h5
  |---SampleID2
  | |---abundance.h5
  |---SampleID3
  | |---abundance.h5
  ...

* 'gene_sets': path to gene modules CSV (disease/molecular pathway related
  modules)

* 'cell_sets': Path to the gene modules CSV (cell-related modules)

* 'visualized_control_samples': Boolean. If false only treatment samples are
  visualized.

* 'generate_bar_plot': Boolean. Generate gene-sets module visualizations

* 'generate_cell_plots': Boolean. Generate cell-sets module visualizations.

* 'generate_spider_plot': Boolean. Generate the above plots as a spider plot

### Output
Generates gene module expression plots for every treatment (and control)
sample for every gene module in the 'generated' directory. Also, a
composition plot is generated that visualizes the log2 expression composition
for each sample.

### To run
Using R version 4.5 or later, navigate to the ImmunomeAnalyzeR directory and run the following command:
```
Rscript immunomeAnalyzeR.R <CONFIG_SETTING>
```
Where <CONFIG_SETTING> is the name of the settings group to use (ex. 'default' or 'production').

### Questions?
Email me at benmckgreen@gmail.com
