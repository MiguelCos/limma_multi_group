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
5. Click `Source` in the top right corner of the script.

The script will extract the experimental design information from the annotation file and match it to the columns in the expression data file. Then it will perform a general F-test using `limma` to detect which proteins are in general affected by the grouping variable. Using the first group in your annotation file as a baseline.

### Output  

The script will generate an `Output` folder containing: 

- An HTML summary report showing:
  - The summary of the grouping variable (number of samples per group)
  - The number of proteins significantly affected by the grouping variable.
  - A set of boxplots showing the mean values of protein abundance for the top protein hits (those with larger values values for the _F statistic_) for each group.
- A tabular file with the output statistics after `limma` analysis, including the p-values after testing for each protein.

