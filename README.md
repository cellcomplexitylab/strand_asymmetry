# Introduction

The present repository contains the code and the instructions to reproduce the results of the article [Strand
asymmetry influences mismatch resolution during single-strand
annealing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02665-3).


# Requirements

To reproduce the analyses, you will need a Linux distribution, [GNU make](https://www.gnu.org/software/make/)
and [Docker](https://www.docker.com/). The code has only been tested on Ubuntu 20.04.


# Steps to reproduce the analyses
```
git clone https://github.com/cellcomplexitylab/strand_asymmetry
cd strand_asymmetry
make
```

The process will take approximately 2 days and will require 300 GB of available disk space. The data
is downloaded from [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJNA592699?show=reads), which is a
mirror of the repository on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141211). The
data regarding the mapping of the reporters is created in the directory `mapping`, the data regarding
the calls of the repair events is created in the directory `mismatches` and the figures are created
in the repository called `figures`.

# License
The content of the repository is covered by the [MIT license](https://opensource.org/licenses/MIT).
