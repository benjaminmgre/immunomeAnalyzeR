#' Visualize differential gene module expression 
#' 
#' Created by Benjamin Green.
#' "benmckgreen@gmail.com"
#' 
#' @description 
#' Compares two groups (ex. treatment vs. control, healthy vs. disease, etc.) 
#' by comparing expression of groups of related genes.
#' 
#' @details
#' The config.yml file must be configured prior to running the script.
#' 
#' * 'study_design_file': Path to the tab-delimited text file that contains two
#'    columns: a sample identifier (ID) and a group. See 
#'    'study_design_template.txt'.
#'    
#' * 'treatment_name': The name of the group of interest in the study design
#'    file. Usually this is 'patient'. The script will generate a plot for each 
#'    treatment sample.
#'    
#' * 'data_dir': Path to the directory with Kallisto results. The directory must
#'   be structured as the raw Kallisto output: one sub-directory for each 
#'   patient ID with an 'abundance.h5' file inside:
#'   data_dir
#'   |---SampleID1
#'   | |---abundance.h5
#'   |---SampleID2
#'   | |---abundance.h5
#'   |---SampleID3
#'   | |---abundance.h5
#'   ...
#' 
#' * 'gene_sets': path to gene modules CSV (disease/molecular pathway related 
#'   modules)
#' 
#' * 'cell_sets': Path to the gene modules CSV (cell-related modules)
#' 
#' * 'visualized_control_samples': Boolean. If false only treatment samples are 
#'   visualized.
#'   
#' * 'generate_bar_plot': Boolean. Generate gene-sets module visualizations
#'   
#' * 'generate_cell_plots': Boolean. Generate cell-sets module visualizations.
#' 
#' * 'generate_spider_plot': Boolean. Generate the above plots as a spider plot
#'
#' @returns
#' Generates gene module expression plots for every treatment (and control) 
#' sample for every gene module in the 'generated' directory. Also, a 
#' composition plot is generated that visualizes the log2 expression composition
#' for each sample.
#' 


print("Welcome to ImmunomeAnalyzeR")
print("Created by Benjamin Green")
print("Last updated: May 12th, 2025")

# Ensure correct environment is loaded
if (Sys.getenv("RENV_PROJECT") == "") {
  source("renv/activate.R")
}

if (!renv::status()$synchronized) {
  message("Restoring packages from renv.lock ...")
  renv::restore(prompt = FALSE)
}

# Import packages ----
packages <- c("tidyverse", "config", "ensembldb", "EnsDb.Hsapiens.v86", 
              "tximport", "edgeR", "cowplot", "matrixStats", "scales", "fmsb")

for (p in packages) {
  suppressPackageStartupMessages( suppressWarnings(
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  ))
}

# Config ----
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
mode <- "production" # DEBUG
print(paste0("Running with '", mode, "' configuration settings."))
Sys.setenv(R_CONFIG_ACTIVE = mode)
config <- config::get()

# Create output directory if needed
if (!dir.exists("./generated/")) {
  dir.create("./generated/")
}

# Import and read Kallisto results ----
targets <- read.table(config$study_design_file, sep = '\t', header = TRUE)
sample.ids <- targets$Identifier

# patients with condition, controls are 'normal'
patient.ids <- subset(targets, Group == config$treatment_name)$Identifier
control.ids <- subset(targets, Group != config$treatment_name)$Identifier
paths <- file.path(config$data_dir, targets$Identifier, "abundance.tsv")

# Get annotations to convert transcripts to genes
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
# Change column names to make compatible
Tx <- dplyr::rename(Tx, target_id=tx_id)
# Remove unnecessary columns. Transcript ID has to be the first column
Tx <- dplyr::select(Tx, "target_id", "gene_name")

# Import the Kallisto results
# Ignore Tx Version and Ignore After Bar to correctly summarize to gene level
# Use browseVignettes("tximport") for help
print("Loading abundances from Kallisto")
Txi_gene <- suppressMessages(tximport(paths, 
                     type="kallisto", 
                     tx2gene=Tx, 
                     txOut=FALSE, 
                     countsFromAbundance = "lengthScaledTPM", 
                     ignoreTxVersion = TRUE,
                     ignoreAfterBar = TRUE,
))

# CPM normalization and filtering, generating composition plots ----
dge.unmodified <- DGEList(Txi_gene$abundance, group = targets$Group)

keep <- filterByExpr(dge.unmodified)
dge.filtered <- dge.unmodified[keep, , keep.lib.sizes = FALSE]

dge.normalized <- calcNormFactors(dge.filtered, method = "TMM")

log2.cpm <- cpm(dge.normalized, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
log2.cpm.df <- as_tibble(log2.cpm, rownames="gene_id")

# Prepare unmodified abundances for plotting
plot.counts <- as_tibble(dge.filtered$counts)
colnames(plot.counts) <- sample.ids
plot.counts$gene_id <- rownames(dge.filtered$counts)
plot.counts.pivot <- pivot_longer(plot.counts,
                                  cols=all_of(sample.ids),
                                  names_to = "sample",
                                  values_to = "expression")

# Prepare filtered, transformed, normalized for plotting
colnames(log2.cpm.df) <- c("gene_id", sample.ids)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, 
                                  cols=all_of(sample.ids),
                                  names_to = "sample",
                                  values_to = "expression")

composition.p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=sample, y=expression, fill=sample) +
  geom_violin(trim=FALSE, show.legend=FALSE) +
  stat_summary(fun="median",
               geom="point",
               shape=95,
               size=10,
               color="black",
               show.legend=FALSE) +
  labs(y="log2 expression", x="Sample",
       title="Filtered, TMM Normalized, CPM Normalized, and Log2 Transformed Expression"
  ) +
  theme_bw()

comp.plot.title <- paste("generated/composition_plots/composition_", Sys.time(), sep='')
comp.plot.title <- paste(comp.plot.title, ".png", sep='')
comp.plot.title <- gsub(" ", "_", comp.plot.title)
if (!dir.exists("generated/composition_plots/")) {
  dir.create("generated/composition_plots/")
}

png(comp.plot.title, width=700, height=500)
# plot_grid(composition.p1, composition.p2, labels=c('A', 'B'), label_size = 12)
suppressMessages(plot(composition.p1))
invisible(dev.off())

# Calculate scores ----
gene_sets <- read_csv(config$gene_sets, show_col_types = FALSE)

pathway.columns <- c("pathway_name", sample.ids)
pathway.scores <- data.frame(matrix(nrow=0, ncol=length(pathway.columns)))
colnames(pathway.scores) <- pathway.columns


for (set.name in colnames(gene_sets)) {
  gene.set <- gene_sets[,set.name]
  gene.set <- gene.set[!is.na(gene.set)]
  
  expression.df <- log2.cpm.df %>%
    dplyr::filter(gene_id %in% gene.set)
  
  genes <- expression.df$gene_id
  
  control.expression.df <- expression.df %>%
    dplyr::select(all_of(control.ids))
  control.expression <- data.matrix(control.expression.df)
  control.mean <- rowMeans(control.expression)
  control.sd <- rowSds(control.expression)
  
  expression <- subset(data.matrix(expression.df), select=-gene_id)
  
  z.scores <- sweep(expression, 1, control.mean, "-")
  z.scores <- sweep(z.scores, 1, control.sd, "/")
  
  mean.z <- colMeans(z.scores)
  mean.z <- as.numeric(mean.z)
  new.row <- c(set.name, mean.z)
  pathway.scores[nrow(pathway.scores) + 1,] <- new.row
}

# Cell type sets
cell_sets <- read_csv(config$cell_sets, show_col_types = FALSE)

cell.pathway.columns <- c("pathway_name", sample.ids)
cell.pathway.scores <- data.frame(matrix(nrow=0, ncol=length(cell.pathway.columns)))
colnames(cell.pathway.scores) <- cell.pathway.columns

for (set.name in colnames(cell_sets)) {
  gene.set <- cell_sets[,set.name]
  gene.set <- gene.set[!is.na(gene.set)]
  
  expression.df <- log2.cpm.df %>%
    dplyr::filter(gene_id %in% gene.set)
  
  genes <- expression.df$gene_id
  
  control.expression.df <- expression.df %>%
    dplyr::select(all_of(control.ids))
  control.expression <- data.matrix(control.expression.df)
  control.mean <- rowMeans(control.expression)
  control.sd <- rowSds(control.expression)
  
  expression <- subset(data.matrix(expression.df), select=-gene_id)
  
  z.scores <- sweep(expression, 1, control.mean, "-")
  z.scores <- sweep(z.scores, 1, control.sd, "/")
  
  mean.z <- colMeans(z.scores)
  mean.z <- as.numeric(mean.z)
  new.row <- c(set.name, mean.z)
  cell.pathway.scores[nrow(cell.pathway.scores) + 1,] <- new.row
}

# GRAPHING PARAMETERS ----
# Define the immunogram reference limits
y.green.min <- -2
y.green.max <- 2
y.yellow.max <- 3
y.yellow.min <- -3
y.red.min <- -6
y.red.max <- 6

# GRAPH DIRECTORY FUNCTION ----
immune.png <- function(id, tag) {
  data.dir <- file.path('generated', id)
  if (!dir.exists(data.dir)) {
    dir.create(data.dir)
  }
  file.path(data.dir, paste(id, '_', tag, '.png', sep=''))
}

# BARPLOT FUNCTION ----
immune.barplot <- function(data, id=NULL) {
  n.pathways <- length(data$pathway_name)
  # Format data
  data.long <- pivot_longer(data, 
                            cols=-c("pathway_name"), 
                            names_to= "sample",
                            values_to= "expression")
  data.long$expression <- sapply(data.long$expression, as.numeric)
  # # Order data descending
  # data.long <- data.long[order(data.long$expression, decreasing = TRUE),]
  # data.long$pathway_name <- factor(data.long$pathway_name, levels = data.long$pathway_name[order(data.long$expression)])
  
  p <- ggplot(data.long, fill="black") +
    geom_bar(aes(x=pathway_name, y=expression), stat="identity", show.legend = FALSE) + 
    scale_y_continuous(breaks=c(-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10), limits = c(-6, 6), oob=squish) +
    annotate("rect", xmin=0, xmax=n.pathways + 1, ymin=y.yellow.min, ymax=y.green.min, alpha=0.3, fill="yellow") +
    annotate("rect", xmin=0, xmax=n.pathways + 1, ymin=y.green.max, ymax=y.yellow.max, alpha=0.3, fill="yellow") +
    annotate("rect", xmin=0, xmax=n.pathways + 1, ymin=y.red.min, ymax=y.yellow.min, alpha=0.2, fill="red") +
    annotate("rect", xmin=0, xmax=n.pathways + 1, ymin=y.yellow.max, ymax=y.red.max, alpha=0.2, fill="red") + 
    labs(y="Expression z-score", x="Pathway", title=paste(id, "Immune Pathway Analysis", sep=" ")) +
    coord_flip()
  
  invisible(print(p))
  
  img.path <- immune.png(id, "Barplot")
  suppressMessages(invisible(ggsave(img.path, plot = p)))
}

# SPIDERPLOT FUNCTION ----
immune.spiderplot <- function(data, id=NULL) {
  # Process data
  data.transpose <- data.table::transpose(data)
  
  # Set the first row as the column names
  colnames(data.transpose) <- data.transpose[1,]
  data.transpose <- data.transpose[-1,]
  
  # https://r-graph-gallery.com/142-basic-radar-chart.html
  # Add max and min lines to the dataframe
  n.columns <- dim(data.transpose)[2]
  data.transpose <- rbind(rep(6, n.columns), rep(-6, n.columns), data.transpose)
  data.transpose <- data.transpose %>% mutate_all(as.numeric)
  
  img.path <- immune.png(id, "Spiderplot")
  png(img.path)
  
  # Grey area ON
  # data.transpose <- rbind(data.transpose, rep(2, n.columns))
  # p <- radarchart(data.transpose,
  #            cglty=1,
  #            cglcol = 'black',
  #            axislabcol = 'blue',
  #            pcol=c('red', 'grey'),
  #            pfcol=c(NA, "#99999980"),
  #            plty = c(1, 1),
  #            pty=c(16, 32),
  #            seg=6,
  #            axistype = 1,
  #            caxislabels = c(-6, -4, -2, 0, 2, 4, 6),
  #            title=paste(id, "Immune Pathway Analysis", sep=" "))
  # Grey area OFF
  invisible(radarchart(data.transpose,
              cglty=1,
              cglcol = 'black',
              axislabcol = 'blue',
              pcol='red',
              plty = 1,
              pty=16,
              seg=6,
              axistype = 1,
              caxislabels = c(-6, -4, -2, 0, 2, 4, 6),
              title=paste(id, "Immune Pathway Analysis", sep=" ")))

  dev.off()
}

# Plotting loop ----
print("Generating Plots")
for (sample in sample.ids) {
  if (sample %in% patient.ids | config$visualize_control_samples) {
    sample.data <- pathway.scores[, c("pathway_name", sample)]
    if (config$generate_bar_plot) {
      immune.barplot(sample.data, id = sample)
    }
    if (config$generate_spider_plot) {
      immune.spiderplot(sample.data, id = sample)
    }
    
    # Cell Type Plots
    if (config$generate_cell_plots) {
      cell.data <- cell.pathway.scores[, c("pathway_name", sample)]
      if (config$generate_bar_plot) {
        immune.barplot(cell.data, 
                                 id = paste(sample, "Cell Type", sep=' '))
      }
      if (config$generate_spider_plot) {
        immune.spiderplot(cell.data, 
                                    id = paste(sample, "Cell Type", sep=' '))
      }
    }
  }
}


print("Analysis Complete!")
