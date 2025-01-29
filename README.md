# ecig Spatial Transcriptome
**Spatial transcriptome of prenatal e-cig induced transcriptional changes in rat brain.**

This repository contains scripts and data to reproduce the analysis and figures used in the paper.
Various types of data such as spatial transcriptome (stRNA-seq), single nucles snRNA-seq and snATAC-seq.
Analysis of each data types are saved in separate folders.

## Date
January, 2025.

## Contributors

## Data Analysis Overview
There are 5 main folders that constitute analysis steps of different data types.

## **stRNAseq** 
Spatial transcriptomics data
1. Preprocess the data
2. Run velocity analysis
3. Visualize results

## **stRNAseq_StereoSeq**

## **snRNAseq**

## **RNA Velocity**

### Overview
This section describes the in-depth RNA velocity analysis on spatial transcriptomic data obtained from 10x Genomics Visium assays. As explained, we focus on four distinct rat brain tissue datasets: MA (Male Anterior), FA (Female Anterior), MP (Male Posterior), and FP (Female Posterior). Each dataset includes samples from both control groups and specimens exposed to E-cigarette smoke (during prenatal development), allowing for comparative analyses of gene expression dynamics under different treatment conditions. 

### Data Cleaning and Integration
The analysis through our RNA velocity pipeline begins with meticulous data cleaning to ensure high-quality input for downstream analyses. Using the scvi library, each dataset undergoes filtering to remove low-quality cells based on minimum gene counts and mitochondrial gene expression percentages. Specifically, cells with fewer than 200 genes or high mitochondrial content are excluded to maintain data integrity. After initial filtering, the datasets are integrated using the SCVI (Single-cell Variational Inference) model, which effectively harmonizes data from different samples by accounting for batch effects and other confounding variables. This integration process results in a unified AnnData object that serves as the foundation for subsequent RNA velocity computations.

### Spatial Coordinate Alignment
Accurate spatial mapping is crucial for spatial transcriptomics analysis. The pipeline includes steps to align and, when necessary, rotate spatial coordinates to ensure consistency across all samples. This alignment facilitates the precise spatial mapping of gene expression patterns, enabling the visualization of transcriptional dynamics within the anatomical context of the brain tissues. The rotation of coordinates, applied to specific samples, ensures that spatial data from different experiments are comparable and correctly oriented for integrated analysis.

### RNA Velocity Computation
RNA velocity analysis is performed using the scvelo library, which leverages spliced and unspliced RNA counts to infer dynamic transcriptional changes. The integrated AnnData object undergoes normalization and calculation of moments to prepare the data for velocity estimation. RNA velocity vectors are then computed, providing insights into the directional flow of cellular states and the temporal progression of gene expression. This dynamic information is crucial for understanding the underlying biological processes affected by E-cigarette exposure.

### Visualization
The analysis includes various visualization techniques to interpret RNA velocity results effectively. Velocity embeddings are generated to display the directional flow of cellular states on a UMAP (Uniform Manifold Approximation and Projection) plot, colored by spatial regions. Additionally, latent time scatter plots illustrate the temporal progression of gene expression, while velocity pseudotime plots highlight the inferred trajectory of cells over time. These visualizations collectively offer a comprehensive view of the transcriptional dynamics and spatial organization within the brain tissues.

### Differential Gene Kinetics Analysis
To identify genes with significant velocity changes between control and E-cigarette exposed groups, a differential gene kinetics analysis is conducted. This involves statistical testing using the Mann-Whitney U test and the diffxpy library for robust differential expression analysis. The analysis focuses on specific brain regions, identifying genes that exhibit significant changes in RNA velocity under treatment conditions. The results are visualized through plots such as velocity pseudotime and confidence metrics, providing insights into the molecular mechanisms influenced by E-cigarette smoke exposure.

### Results and Interpretation
The integration of spatial context with dynamic gene expression data reveals nuanced patterns of gene regulation and cellular responses in the rat brain tissues. Differential analyses highlight specific genes and pathways affected by E-cigarette exposure, with variations observed across different brain regions and between sexes. These findings contribute to a deeper understanding of the biological impact of environmental factors on brain function and may inform future studies on gene regulation and neurological health.

### Conclusion
This comprehensive RNA velocity analysis on spatial transcriptomic datasets offers valuable insights into the dynamic transcriptional changes induced by E-cigarette exposure in rat brain tissues. By integrating high-quality data cleaning, robust integration methods, and advanced visualization techniques, the study elucidates the complex interplay between environmental factors and gene expression dynamics. The findings pave the way for future research on the molecular mechanisms underlying environmental impacts on neurological health.

### Dependencies and Setup
The analysis pipeline leverages several Python libraries, including scvi, scvelo, scanpy, diffxpy, numpy, pandas, and matplotlib. Ensure all dependencies are installed and properly configured before running the analysis scripts. Detailed instructions for setting up the environment and executing the analysis are provided in the subsequent sections of this README.

### Usage
To replicate the analysis, follow the steps outlined in the provided scripts. Begin with data cleaning and integration, proceed to RNA velocity computation, and conclude with differential gene kinetics analysis. Visualizations are generated at each stage to facilitate the interpretation of results. Refer to the code comments and function descriptions for detailed guidance on each processing step.

## Data Integration & Deconvolution

## **enrichmentHuman**
1. Download human PsyhcENCODE data
2. Preprocess the data and get rat orthologous data
3. Run analysis
4. Visualize results

