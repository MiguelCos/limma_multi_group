### Script for running multi-group limma analysis ---- 
## Miguel Cosenza 02.07.2020 

### 1. Please provide a meaningful name for the dataset/experiment ----

exper_code <- "Experiment XY"

### 2. How many top significant hits do you want to plot (boxplots group vs intersity)? 

n_top_hits <- 12

### 3. Which groups do you want to compare?  

## Required packages ----

packages <- c("dplyr", "here", "tidyr", "ggplot2", "rmarkdown", "knitr", "reshape")

biopackgs <- c("limma")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
          install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (length(setdiff(biopackgs, rownames(installed.packages()))) > 0){
          if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
          
          BiocManager::install(setdiff(biopackgs, rownames(installed.packages())))
          
}

library(dplyr)
library(stringr)
library(limma)
library(rmarkdown)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(reshape)
library(plotly)

## Load data ----

expr_dat <- read.delim("Data/input_limma.txt")

annot_dat <- read.delim("Data/annotation.txt")

## Define the design matrix ----

groups <- as.factor(annot_dat$Group)

design <- model.matrix(~0+groups)

row.names(design) <- annot_dat$Sample_ID

colnames(design) <- colnames(design) %>% str_remove(., "groups") %>% str_trim()

# 3. DEFINE WHICH GROUPS YOU WISH TO COMPARE ----

contrast.matrix <- makeContrasts(C-B, B-A, C-A, levels=design) # MODIFY THIS LINE

## Prep expression data matrix ----
n_contrasts <- dim(contrast.matrix)[2]

tomat <- dplyr::select(expr_dat,
                       -ID, row.names(design)) %>% as.data.frame()

row.names(tomat) <- expr_dat$ID

mat <- as.matrix(tomat)

## Execute the linear model / limma ----

fit <- lmFit(mat, design = design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

output_limma2 <- topTable(fit2, adjust.method = "BH", number = Inf)

output_limma2$Protein <- row.names(output_limma2)


## Generate output ----

## output tabular list ----

list_tabular <- list()

for (i in 1:n_contrasts){
          
          outlim <- topTable(fit2, coef = i, adjust.method = "BH", number = Inf)
          outlim$Protein <- row.names(outlim)
          outlim$Contrast <- colnames(contrast.matrix)[i]
          list_tabular[[i]] <- outlim
}

names(list_tabular) <- colnames(contrast.matrix)

if(!dir.exists("Output")){dir.create("Output")}

write.table(output_limma2,
            file = "Output/tab_output_general_Ftest_limma.txt",
            row.names = FALSE, col.names = TRUE)

sig_hits <- dplyr::filter(output_limma2, 
                               adj.P.Val <= 0.05) %>% row.names(.)

n_significant <- length(sig_hits)

top_n_hits <- top_n(output_limma2, n = n_top_hits, wt = `F`)


## Extract basic information about each contrasts ----

getinfo <- function(x){
          x %>% dplyr::filter(adj.P.Val <= 0.05) %>% dim(.) %>% .[1]
}

## Prep some boxplots for the top proteins with lowest P-values ----

slim_expr <- pivot_longer(expr_dat, cols = colnames(mat),
                          names_to = "Sample_ID",
                          values_to = "Abundance") 

slim_expr_g <- left_join(slim_expr, annot_dat,
                         by = "Sample_ID")

hits_expr <- filter(slim_expr_g,
                    ID %in% row.names(top_n_hits))


boxplots <- ggplot(hits_expr,
                   aes(x = Group, y = Abundance)) +
          geom_boxplot()+
          facet_wrap(ID ~ .)

## Prep volcano plots for contrasts ----

merged_limmacontr <- reshape::merge_all(list_tabular)

write.table(merged_limmacontr,
            file = "Output/tab_output_per_contrasts_limma.txt",
            row.names = FALSE, col.names = TRUE)


tovolc <- merged_limmacontr %>% mutate(Differentially_expressed = case_when(adj.P.Val <= 0.05 ~ TRUE,
                                                         TRUE ~ FALSE))

volcanoes <- ggplot(data = tovolc, 
                 mapping = aes(x = logFC,
                               y = -log10(adj.P.Val),
                               color = Differentially_expressed)) + 
                    geom_point(alpha = 0.5) + 
                    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
                    facet_wrap(~ Contrast, ncol = 3)+
                    theme(legend.position = "none")


rmarkdown::render(input = here::here("renderReport.R"),
                  output_file = paste0("Output/limma_anova_report_",exper_code))
