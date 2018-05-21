# BindSpace
TFBS models trained with HT-SELEX data using StarSpace embedding method

## installation

### install dependencies
The model file is > 200MB, make sure you have git lfs installed in order to successfully download the it. Follow this instruction to install git lfs:  
https://github.com/git-lfs/git-lfs/wiki/Installation

To install this R packages, please make sure you have installed the most updated versions of the following packages:
* data.table
* Matrix
* PRROC
* Biostrings
* ComplexHeatmap

These packages can be installed in R:  
install.packages(c("data.table","Matrix","PRROC"))  
source("https://bioconductor.org/biocLite.R")  
biocLite(c("ComplexHeatmap", "Biostrings"))  

### download BindSpace to local directory
git clone https://github.com/hy395/BindSpace.git

### install in R
install.packages("BindSpace/", repos=NULL, type="source")  
NOTE: This package can't be installed in R via install_github function because of the large model file.

## test
example code can be found under BindSpace/tests/test_package.R
