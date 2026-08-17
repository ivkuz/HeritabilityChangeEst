# HeritabilityChangeEst

R code used in the study "Genetic influences on educational outcomes during and after the Soviet era: Revisiting evidence from Estonia".

This repository contains the analyses used to estimate and compare genetic associations with educational attainment (EA), occupational status (OS), height, and BMI between Soviet-era and post-Soviet groups in Estonia. 

Sample preparation (ancestry filtering, relatedness) was accomplished as in https://github.com/ivkuz/GeneticMigrationStructureEstonia.

Find a short description of what each of the scripts from the "scripts" directory does below:


  Data preparation and sample definition
  
    filterUnrelated.R  -  Prepares lists of individuals for different analyses and applies relatedness-based filtering. It also creates subsets by birth year, era, participation wave, age threshold, sex, and settlement type.
    makeKinshipMatrixOfRelated.R  -  Processes genomic relationship matrices (GRMs), including thresholding relationship coefficients and constructing GRM files for analyses involving related individuals and G×E comparisons.
    transformEA_INT.R  -  Prepares educational attainment and related phenotype variables, including transformations of categorical measures into inverse-normalised variables.
    EA_YoB_distribution.R  -  Examines the distribution of educational attainment across year-of-birth cohorts and produces diagnostic plots used to assess cohort-related patterns and thresholds.
    manova.R  -  Performs multivariate and univariate comparisons of educational attainment and occupational status between Soviet and post-Soviet groups, including analyses by participation wave.

    
  Heritability analyses
  
    figures_h2_2_3_4.R  -  Generates the main heritability figures, including comparisons between Soviet and post-Soviet groups, analyses by participation wave, family-inclusive model estimates, and analyses using inverse-normalised EA.
    figures_h2_S1_S2.R  -  Generates Supplementary figures, including alternative heritability analyses using GCTA, binary EA, a different relatedness cutoff, and subgroup analyses by participation wave, sex, and settlement type.
    figures_h2_decades_2_3_S3.R  -  Generates the birth-cohort heritability analyses. Heritability is examined across overlapping decade-long birth-year cohorts.
    plot_h2_origSample.R  -  Generates the heritability analysis for the original sample, reported as a supplementary analysis.
    variance_gen_env.R  -  Examines genetic and environmental variance components across birth-year cohorts for educational attainment, occupational status, height, and BMI.

    
  Polygenic-score analyses
  
    pgsR2.R  -  Main analysis script for the proportion of trait variance explained by polygenic scores (R2). The script includes:
      data preparation;
      the main R2 analysis;
      weighting;
      matching of phenotype distributions.
      pgsR2_origSample.R  -  Performs the polygenic-score variance-explained analysis in the original sample.
    figure_R2_5.R  -  Generates figures, showing trait variance explained by polygenic scores in Soviet and post-Soviet groups, including analyses of inverse-normalised traits.
    figures_R2_S4_S7.R  -  Generates Supplementary Figures, including analyses by sex, settlement type, and other phenotype/subsample definitions.
    figures_R2_decades_5_S5_S6.R  -  Examines variance explained by polygenic scores across birth cohorts and generates the corresponding main and supplementary figures.
    plot_R2_origSample.R  -  Generates the supplementary analysis of polygenic-score variance explained in the original sample.
    plot_R2_matching.R  -  Generates figures for analyses in which subsamples from the Soviet and post-Soviet groups are selected to match the distribution of educational attainment in another subcohort.
    figure_R2_weighting.R  -  Generates the figure and evaluates polygenic-score variance explained after applying different weighting schemes, including within-group and overall weighting.

      
  G×E and interaction analyses
  
    regressionInteractions.R  -  Performs gene-by-environment analyses using linear regression. The analyses examine interactions between polygenic scores and the Soviet/post-Soviet era and generate corresponding figures.


  Composite figures
  
    figures_composite_2_3_5.R  -  Combines previously generated figure components into the final composite versions.

    
  Bayesian analyses
  
    bayesian_estimates.R  -  Calculates Bayesian posterior estimates for differences in heritability, differences in polygenic-score variance explained, and PGS × era interaction effects. The script uses normal priors and produces posterior means, posterior standard deviations, 95% credible intervals, and probabilities of positive effects.
    
Main analyses

The repository therefore covers several complementary approaches to testing whether genetic influence on social outcomes changed following the end of the Soviet era:

Heritability comparisons
Comparison of genetic and environmental contributions to variation in EA and other traits between Soviet and post-Soviet groups.
Birth-cohort analyses
Examination of changes in heritability and polygenic-score variance explained across overlapping birth-year cohorts.
Polygenic-score analyses
Estimation of the proportion of phenotypic variance explained by polygenic scores and comparison of these estimates between eras.
Distribution matching
Analyses designed to reduce the possibility that differences in the phenotype distributions between groups explain differences in PGS-based estimates.
Weighting
Analyses using weighted subsamples to account for differences in group composition.
Original-sample analyses
Additional analyses using the original sample definitions as robustness checks.
G×E analyses
Direct tests of interactions between genetic propensity, represented by polygenic scores, and the social environment/era.
Additional robustness analyses
Analyses using alternative phenotypic transformations, relatedness thresholds, participation waves, sex, settlement type, and other subgroup definitions.
Data and computational requirements

The repository contains the R analysis code but does not contain the underlying individual-level data or intermediate genetic-analysis results.
Because these data are not included in the repository, the scripts are intended to document and reproduce the analyses within the corresponding research environment rather than provide a completely self-contained pipeline.
Some scripts also contain alternative or commented-out analyses. These are retained to document analysis decisions and intermediate approaches.

License

See the repository LICENSE file for licensing information.
