Description of how the R code in this repository was used

Files output from the sequence processing were used along with a meta data file for the analysis in this code

Relative Abundance
FF2_diet_rel_abundances.R, FF2_inoc_rel_abundances.R and FF2_RS_but_stacked_bar.R
were all used for determining the relative abundances of all taxa at various taxonomic levels either across inoccula or samples and displayed in the stacked bar chart and table within the manuscript.

Differential Abundance
FF2_ANCOM.R, FF2_DESEQ2.R and FF2_LEfSe.R
were used to perform differential abundance analysis by each of these methods. For LEfSe the output files were used at the Huttenhower lab Galaxy web server and run with default parameters with "Inoculum" input as "Individual". 

SCFA changes
FF2_SCFA_analysis.R
was used to determine significant changes in SCFA between conditions and generate the bubble plots in the paper.

Diversity Analysis
FF2_diversity.R
was used to measure and compare alpha (Shannon, Inverse Simpson, Faith) and beta (Bray-Curtis, Aitchison, Weighted UniFrac) diversity across treatments and between treatments and controls. 

Correlation Analysis
(1) FF2_but_SCFA_corr.R and (2) FF2_Microbiome_Ferm_corr.R
were (1) used to determine correlations between butyrate and other fermentation products and pH 
and (2) used to determine correlations between butyrate producing organisms and butyrate levels, lactate levels, pH and RS degrading organisms

