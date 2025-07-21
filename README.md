# HeritabilityChangeEst
R code used in the study "No evidence of a change in genetic influence on social outcomes following the end of the Soviet era in Estonia".

Sample preparation (ancestry filtering, relatedness) was accomplished as in https://github.com/ivkuz/GeneticMigrationStructureEstonia.


Find a short description of what each of the scripts from the "R_code" directory does below:


plot_h2_main.R -  Fig. 1 (cutoff 15), SF4 (cutoff 10) and SF 5 (GCTA and binary EA):
                  Heritability in the post-Soviet and Soviet groups and the groups additionally divided by the wave of participation
                  
plot_h2_origSample - SF1: Heritability in the post-Soviet and Soviet groups in the original sample

plot_h2_decades - Fig. 4 and SF 8: Variance explained by PGS in the birth cohorts by half-decade


pgsR2.R - All the calculations with R2:
          1. Prepare the data for the analysis 
          2. Main R2 analysis 
          3. Weighting 
          4. Matching distributions 
          5. Original sample 


plot_R2_main.R - Plot Fig. 3 and SF 7: Trait variance explained by PGS in the post-Soviet and Soviet groups

plot_R2_origSample.R -SF3: Variance explained by PGS in the post-Soviet and Soviet groups in the original sample

plot_R2_decades.R - Fig. 2  and SF 6: Heritability across different birth year cohorts 


plot R2_matching.R - Fig. 5 and SF 9:
                    Variance of EA explained by PGSEA in subsamples from the groups (defined by wave participation and era) 
                    selected based on the distribution of EA in another subcohort

plot_R2_weighting.R - Fig. 6: Trait variance explained by PGS in weighted subsamples 
