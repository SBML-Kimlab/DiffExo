# DiffExo  

## The DiffExo pipeline was developed to perform the following tasks: 
1) Normalization of DNA-binding protein binding intensity
2) Identification of differentially binding sites between two or more sample groups using DEseq2
3) Statistics analysis of normalized binding intensity 
4) Generation of plots including box plots, scatter plots, heat maps, and volcano plots.

Two dataset files (CSV), including information on the data directories and statistical groups (control or experiment), are formatted before starting the DiffExo pipeline. For preprocessing of ChIP-mini data for the DiffExo pipeline, the sequencing raw files (FASTQ) are converted into alignment files (BAM) and visualization files (GFF) using the cloud-based ChIP-exo analysis pipeline, ChEAP. Followed by visualization, a Deep-learning optimized ChIP-exo peak calling suite (DEOCSU) is utilized to predict binding peaks of target DNA-binding protein using GFF files. Since all the read count data from each ChIP-mini sequencing library are not normally distributed as they are in an RNA-seq library, a negative binomial distribution-based algorithm was adopted to estimate differentially binding sites. 

**1) localPNC.py**: calculate the normalized binding intensity of DNA-binding protein using an RPPM unit and generate a pickle file to convey the normalized data to the next process.  
**2) StatAnalysis.py**: proceed with statistical test (parametric or non-parametric test) using the RPPM.pickle file.  
**3) DEseq2Tools**: convert overlapping binding sites into reference file for DEseq2, and differential binding sites are calculated in R script (DEseq2_cal). In addition, the raw result file from DEseq2 merges with RPPM and information on binding sites to generate the final Deseq2 result file (CSV) using DEseq2Tools.  
**4) PlotTools.py**: generate four types of plots: box plots, heat maps, scatter plots for visualizing binding intensity, and volcano plots for visualizing DEseq2 results.  

![Figure_S3](https://github.com/SBML-Kimlab/DiffExo/assets/67301306/f4740828-9aad-4adc-b63e-f8ffff917542)

## List of PIP dependencies
* **numpy** (1.21.6)  
* **pandas** (1.2.3)  
* **pysam** (0.16.0.1)  
* **scipy** (1.7.3)  
* **matplotlib** (3.5.2)  
* **bioinfokit** (2.1.0)  

## List of R dependencies
* **GenomicFeatures** (1.34.8)  
* **SummarizedExperiment** (1.12.0)  
* **DESeq2** (1.22.2)  

## Bowtie mapping option
```
bowtie -p [core] -S [bt_index] [FASTQ_R1.fastq]
```

## Citation information
ChIP-mini: a low-input ChIP-exo protocol for elucidating DNA-binding protein dynamics in intracellular pathogens. Park JY#, Jang M#, Choi E, Lee SM, Bang I, Woo J, Kim S, Lee EJ*, Kim D*. Nucleic Acids Res. 2025 Jan 27. 


## References
Love, M.I., Huber, W., and Anders, S.J.G.b. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. 15, 1-21.  
Eder, T., and Grebien, F.J.G.B. (2022). Comprehensive assessment of differential ChIP-seq tools guides optimal algorithm selection. 23, 119.  
Huang, Y., Wang, J.Y., Wei, X.M., Hu, B.J.A.M., and Materials (2014). Bioinfo-Kit: a sharing software tool for bioinformatics. 472, 466-469.  

If you have any questions, please feel free to contact the author via email. [jjypark56@gmail.com]
