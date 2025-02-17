# Lsd1 ko/ki whole bone marrow scRNA-Seq

Our study (Staehle et al., currently in revision at Nature Communcations - future link) examines the enzymatic and non-enzymatic functions of KDM1A (Lsd1) in hematopoiesis. To this end we performed scRNA-Seq of whole bone marrow in Lsd1 ko, Lsd1 ko BMT (secondary bone marrow transplant) and Lsd1 ei (enzymatically inactive knock-in) mice as well as their respective littermates. The final analysis includes a total of 28.869 cells across n = 13 murine donors.

This repository serves the dual purpose of making our scRNA-Seq analysis workflow transparent, and should also allow to reproduce the main figures of the manuscript (if that's not the case - get in touch). The raw counts and/or .fastq files required to run these analyses, as well as processed objects that include metadata can be found under [GSE286398](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286396).

## Running this yourself:

### R packages
Install and activate renv, then run:
   ```r
   renv::restore()
   ```
### conda environments
Create these conda environments
   ```r
   conda env create -f environment_scanpy_sce.yaml
   conda env create -f environment_scvi-tools.yaml
   ```

## Main analysis tools

### scanpy
https://github.com/scverse/scanpy
### scvi-tools
https://github.com/scverse/scvi-tools
### seurat
https://github.com/satijalab/seurat

For a full reference list, please consult the manuscript.