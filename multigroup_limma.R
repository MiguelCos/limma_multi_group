### Script for running multi-group limma analysis ---- 
## Miguel Cosenza 02.07.2020 

## Input parameters ----

### 1. Please provide a meaningful name for the dataset/experiment ----

exper_code <- "Experiment XY"

### 2. How many top significant hits do you want to plot (boxplots group vs intersity)? 

n_top_hits <- 12

## Required packages ----

packages <- c("dplyr", "here", "tidyr", "ggplot2", "rmarkdown", "knitr")

biopackgs <- c("limma")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
          install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (length(setdiff(biopackgs, rownames(installed.packages()))) > 0){
          if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
          
          BiocManager::install(setdiff(biopackgs, rownames(installed.packages())))
          
}

library(dplyr)
library(limma)
library(rmarkdown)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)

## Load data ----

expr_dat <- read.delim("Data/input_limma.txt")

annot_dat <- read.delim("Data/annotation.txt")

## Define the design matrix ----

groups <- as.factor(annot_dat$Group)

design <- model.matrix(~groups)

row.names(design) <- annot_dat$Sample_ID

## Prep expression data matrix ----

tomat <- dplyr::select(expr_dat,
                       -ID, row.names(design)) %>% as.data.frame()

row.names(tomat) <- expr_dat$ID

mat <- as.matrix(tomat)

## Execute the linear model / limma ----

fit <- lmFit(mat, design = design)

fit <- eBayes(fit)

output_limma <- topTable(fit, adjust.method = "BH", number = Inf)

output_limma$Protein <- row.names(output_limma)


## Generate output ----

if(!dir.exists("Output")){dir.create("Output")}

write.table(output_limma,
            file = "Output/tab_output_multigroup_limma.txt",
            row.names = FALSE, col.names = TRUE)


sig_hits <- dplyr::filter(output_limma, 
                               adj.P.Val <= 0.05) %>% row.names(.)

n_significant <- length(sig_hits)

top_n_hits <- top_n(output_limma, n = n_top_hits, wt = `F`)

## Prep some boxplots for the top proteins with lowest P-values ----

slim_expr <- pivot_longer(expr_dat, cols = colnames(mat),
                          names_to = "Sample_ID",
                          values_to = "Abundance") 

slim_expr_g <- left_join(slim_expr, annot_dat,
                         by = "Sample_ID")

slim_expr_sig <- filter(slim_expr_g,
                        )

hits_expr <- filter(slim_expr_g,
                    ID %in% row.names(top_n_hits))


boxplots <- ggplot(hits_expr,
                   aes(x = Group, y = Abundance)) +
          geom_boxplot()+
          facet_wrap(ID ~ .)

rmarkdown::render(input = here::here("renderReport.R"),
                  output_file = paste0("Output/limma_anova_report_",exper_code))
