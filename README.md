# `limma`-based global F-test for protein expression data for more than 2 groups

The script would take one expression data file and an annotation file and execute a general F-test (equivalent to ANOVA), to evaluate which proteins are affected in general by the grouping variable.

## How to use:

###  Input requirements:

- 1 File named `input_limma.txt`. It should be tab separated and should contain normalized protein expression data. The first column should be named `ID` and should contain protein IDs (i.e. Uniprot or Symbol). Every other column should correspond to an experimental sample / individual, and should contain the corresponding protein expression values. These column should have a distinctive name.

- 1 File named `annotation.txt` that would map the names of your experimental samples / individuals to each experimental group. The first column should be called `Sample_ID` and should contain exactly the same names/sample codes given in the columns of the expression data file. The second column should be named `Group` and should contain the group that corresponds to each individual / sample.

### Step-by-step:

1. Download this repository in your local computer and initialize it as an R Project in RStudio.
2. Place the two input files in the `Data` folder.
3. Open the `multigroup_limma.R` script. 
4. Answer the questions on line `8` and line `12` of the script.
5. Define the paired contrasts that you are interested to evaluate after general F-test, by modifiying `contrast.matrix` input at line `60`. (Check below how to set up your contrast matrix)
6. Click `Source` in the top right corner of the script.

The script will extract the experimental design information from the annotation file and match it to the columns in the expression data file. Then it will perform a general F-test using `limma` to detect which proteins are in general affected by the grouping variable. Using the first group in your annotation file as a baseline.


### Modifying your contrast matrix (line `60`)

The contrast matrix would define which groups should be compared against each other.

In line `60` we define the contrast matrix for our comparisons such as this example:

```
contrast.matrix <- makeContrasts(C-B, B-A, C-A, levels=design)
```

In this example, we are testing for 3 groups (A, B and C), and are doing all pairwise comparisons between each other.

Following on that, we are defining our contrast matrix by setting `makeContrasts(C-B, B-A, C-A, levels = design)`, which indicates that we want to test C vs B, B vs A and C vs A. 

If you want to test only C vs A, then you need to modify the function input accordingly; such as `makeContrasts(C-A, levels = design)`

Notice that every contrast is defined by `X-Y`, with `X` beign the numerator of the comparison. Every interesting contrast should be separated by a comma `,`. Nothing else on the line `60` should be modified. 


### Output  

The script will generate an `Output` folder containing: 

- An HTML summary report showing:
  - The summary of the grouping variable (number of samples per group)
  - The number of proteins significantly affected by the grouping variable.
  - A set of boxplots showing the mean values of protein abundance for the top protein hits (those with larger values values for the _F statistic_) for each group.
  - Volcano plots for each pairwise comparison tested.
- A tabular file with the output statistics after `limma` analysis, including the p-values after testing for each protein.
- A tabular file with the output statistics after `limma` analysis for pair-wise comparisons. This contains a `Contrast` variable for each contrast, and includes log2-FC values.

