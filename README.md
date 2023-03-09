# rnaSeqPreprocessing

[Under development] Comprehensive Nextflow-based pipeline for RNA-seq analysis of cancer samples. 

Tools:
* QC with FastQC
* Read alignment with STAR
* Gene expression quantifiaction with HTseq
* Cancer type prediction using kNN and The Cancer Genome Atlas as a reference
* Alignment post processing for variant calling with GATK4
* Variant annotation with VEP
* Conversion of annotated variants from VCF to MAF format
* Isoform-level gene expression quantification with Kallisto
* RNA-based fusion gene prediction with Pizzly
