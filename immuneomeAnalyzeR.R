# Immunome AnalyzeR
# Created by Benjamin Green

print("Welcome to ImmunomeAnalyzeR")
print("Created by Benjamin Green")
print("Last updated: August 28th, 2023")

# Import packages ----
library(tidyverse)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(matrixStats)
library(cowplot)
library(config)
library(ggplot2)
library(scales)
library(plyr)
library(data.table)
library(fmsb)

# Config ----
Sys.setenv(R_CONFIG_ACTIVE = "production")
config <- config::get()

# Check if "generated" dir exists
if (!dir.exists("generated/")) {
  dir.create("generated/")
}

# Import and read Kallisto results ----
targets <- read_tsv(config$study_design_file)
sample.ids <- targets$Identifier
patient.ids <- subset(targets, Group == 'patient')$Identifier
control.ids <- subset(targets, Group == 'control')$Identifier
paths <- file.path(config$data_dir, targets$Identifier, "abundance.tsv")

# Get annotations to convert transcripts to genes
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
# Change column names to make compatible
Tx <- dplyr::rename(Tx, target_id=tx_id)
# Remove unessisary columns. Transcript ID has to be the first column
Tx <- dplyr::select(Tx, "target_id", "gene_name")

# Import the Kallisto results
# Ignore Tx Version and Ignore After Bar to correctly summarize to gene level
# Use browseVignettes("tximport") for help
Txi_gene <- tximport(paths, type="kallisto", 
                     tx2gene=Tx, 
                     txOut=FALSE, 
                     countsFromAbundance = "lengthScaledTPM", 
                     ignoreTxVersion = TRUE,
                     ignoreAfterBar = TRUE,
)

# CPM normalization and filtering, generating composition plots ----

# Non-filtered, non-normalized
myDGEList <- DGEList(Txi_gene$counts)
cpm <- cpm(myDGEList)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames="gene_id")
colnames(log2.cpm.df) <- c("gene_id", sample.ids)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, 
                                  cols=all_of(sample.ids),
                                  names_to = "sample",
                                  values_to = 'expression')
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
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
  ) +
  theme_bw()

# Filtered, non-normalized
keepers <- rowSums(cpm>1)>=2
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames='gene_id')
colnames(log2.cpm.filtered.df) <- c('gene_id', sample.ids)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df,
                                           cols=sample.ids,
                                           names_to='sample',
                                           values_to='expression')

composition.p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=sample, y=expression, fill=sample) +
  geom_violin(trim=FALSE, show.legend=FALSE) +
  stat_summary(fun="median",
               geom="point",
               shape=95,
               size=10,
               color="black",
               show.legend=FALSE) +
  labs(y="log2 expression", x="Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="Filtered, non-normalized",
  ) +
  theme_bw()

# Filtered, normalized
myDGEList.norm <- calcNormFactors(myDGEList.filtered, method="TMM")
log2.cpm.norm <- cpm(myDGEList.norm, log=TRUE)
log2.cpm.norm.df <- as_tibble(log2.cpm.norm, rownames="gene_id")
colnames(log2.cpm.norm.df) <- c("gene_id", sample.ids)
log2.cpm.norm.df.pivot <- pivot_longer(log2.cpm.norm.df,
                                       cols=sample.ids,
                                       names_to = "sample",
                                       values_to = "expression")

composition.p3 <- ggplot(log2.cpm.norm.df.pivot) +
  aes(x=sample, y=expression, fill=sample) +
  geom_violin(trim=FALSE, show.legend=FALSE) +
  stat_summary(fun="median",
               geom="point",
               shape=95,
               size=10,
               color="black",
               show.legend=FALSE) +
  labs(y="log2 expression", x="Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="Normalized",
  ) +
  theme_bw()

comp.plot.title <- paste("generated/composition_plots/composition_", Sys.time(), sep='')
comp.plot.title <- paste(comp.plot.title, ".png", sep='')
comp.plot.title <- gsub(" ", "_", comp.plot.title)
if (!dir.exists("generated/composition_plots/")) {
  dir.create("generated/composition_plots/")
}

png(comp.plot.title, width=700, height=500)
plot_grid(composition.p1, composition.p2, composition.p3, labels=c('A', 'B', 'C'), label_size=12) 
dev.off()

# Calculate scores ----

# Gene sets
gene_sets <- read_csv(config$gene_sets)

pathway.columns <- c("pathway_name", sample.ids)
pathway.scores <- data.frame(matrix(nrow=0, ncol=length(pathway.columns)))
colnames(pathway.scores) <- pathway.columns


for (set.name in colnames(gene_sets)) {
  gene.set <- gene_sets[,set.name]
  gene.set <- gene.set[!is.na(gene.set)]
  
  expression.df <- log2.cpm.norm.df %>%
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
cell_sets <- read_csv(config$cell_sets)

cell.pathway.columns <- c("pathway_name", sample.ids)
cell.pathway.scores <- data.frame(matrix(nrow=0, ncol=length(cell.pathway.columns)))
colnames(cell.pathway.scores) <- cell.pathway.columns


for (set.name in colnames(cell_sets)) {
  gene.set <- cell_sets[,set.name]
  gene.set <- gene.set[!is.na(gene.set)]
  
  expression.df <- log2.cpm.norm.df %>%
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
  
  print(p)
  
  img.path <- immune.png(id, "Barplot")
  ggsave(img.path)
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
  p <- radarchart(data.transpose,
                  cglty=1,
                  cglcol = 'black',
                  axislabcol = 'blue',
                  pcol='red',
                  plty = 1,
                  pty=16,
                  seg=6,
                  axistype = 1,
                  caxislabels = c(-6, -4, -2, 0, 2, 4, 6),
                  title=paste(id, "Immune Pathway Analysis", sep=" "))
  
  
  print(p)
  dev.off()
}

# Plotting loop ----
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
        immune.barplot(cell.data, id = paste(sample, "Cell Type", sep=' '))
      }
      if (config$generate_spider_plot) {
        immune.spiderplot(cell.data, id = paste(sample, "Cell Type", sep=' '))
      }
    }
  }
}


print("Analysis Complete!")
