# LymphoSeq2
Adaptive Immune Receptor Repertoire Sequencing (AIRR-seq) provides a unique opportunity to interrogate the adaptive immune repertoire under various clinical conditions. The utility offered by this technology has quickly garnered interest from a community of clinicians and researchers investigating the immunological landscapes of a large spectrum of health and disease states. LymphoSeq2 is a toolkit that allows users to import, manipulate and visualize AIRR-Seq data from various AIRR-Seq assays such as Adaptive ImmunoSEQ and BGI-IRSeq, with support for 10X VDJ sequencing coming soon. The platform also supports the importing of AIRR-seq data processed using the MiXCR pipeline. LymphoSeq2 builds on the BioConductor package LymphoSeq and introduces tidyverse syntax.

## Installation:
```{r} 
devtools::install_github("shashidhar22/LymphoSeq2", build_vignette=TRUE)
```
