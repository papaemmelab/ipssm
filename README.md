# ipssm

R Package for the Molecular International Prognostic Scoring System (IPSS-M) for Myelodysplastic Syndromes.


## :rocket: Installation instructions

```R
# install devtools if you don't have it already for easy installation
# install.packages("devtools")
library(devtools)
install_github("papaemmelab/ipssm")
```


## :boom: Usage

The worflow below consists of 4 simple steps, namely 1) Read your input data file and perform some validation on the data, 2) Process the variables in a suitable format for the model, 2) Calculate the IPSS-M risk score and risk category (under the worst, mean, and best scenario to account for missing data), 4) Annotate the results.

```R
# load the ipssm library
library("ipssm")
# path to your input data file
path.file <- system.file("extdata", "IPSSMexample.csv", package = "ipssm") 
#path.file <- system.file("extdata", "IPSSMexample.xlsx", package = "ipssm") # equivalent

# 1) Read and Validate File
dd <- IPSSMread(path.file)

# 2) Process User Input Variables into Model Variables
dd.process <- IPSSMprocess(dd)

# 3) Calculate IPSS-M
dd.res <- IPSSMmain(dd.process)

# 4) Annotate Results
dd.annot <- IPSSMannotate(dd.res)
```



## :question: Question

Any questions feel free to contact [ElsaB](https://elsab.github.io/).
