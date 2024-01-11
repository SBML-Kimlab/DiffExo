# DiffExo

## The DiffExo pipeline was developed to perform the following tasks: 
1) Normalization of TF binding intensity
2) Identification of differentially binding sites between two or more sample groups using DEseq2
3) Statistics analysis of normalized binding intensity 
4) Generation of plots including box plots, scatter plots, heat maps, and volcano plots.

Two dataset files (CSV), including information on the data directories and statistical groups (control or experiment), are formatted before starting the DiffExo pipeline. For preprocessing of ChIP-mini data for the DiffExo pipeline, the sequencing raw files (FASTQ) are converted into alignment files (BAM) and visualization files (GFF) using the cloud-based ChIP-exo analysis pipeline, ChEAP. Followed by visualization, a Deep-learning optimized ChIP-exo peak calling suite (DEOCSU) is utilized to predict binding peaks of target TF using GFF files. Since all the read count data from each ChIP-mini sequencing library are not normally distributed as they are in an RNA-seq library, a negative binomial distribution-based algorithm was adopted to estimate differentially binding sites. 

**1) localPNC.py**: calculates the normalized binding intensity of TF using an RPPM unit and generates a pickle file to convey the normalized data to the next process.  
**2) StatAnalysis**.py: proceeds with statistical test (parametric or non-parametric test) using the RPPM.pickle file.  
**3) DEseq2Tools**: converts overlapping binding sites into reference file for DEseq2, and differential binding sites are calculated in R script (DEseq2_cal). In addition, the raw result file from DEseq2 merges with RPPM and information on binding sites to generate the final Deseq2 result file (CSV) using DEseq2Tools.  
**4) PlotTools.py**: generates four types of plots: box plots, heat maps, scatter plots for visualizing binding intensity, and volcano plots for visualizing DEseq2 results.  

![Figure_S3](https://github.com/SBML-Kimlab/DiffExo/assets/67301306/f4740828-9aad-4adc-b63e-f8ffff917542)
