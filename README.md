# HeritabilityChangeEst
R code used in the study "No evidence of a change in genetic influence on social outcomes following the end of the Soviet era in Estonia".

Sample preparation (ancestry filtering, relatedness) was accomplished as in https://github.com/ivkuz/GeneticMigrationStructureEstonia.


Find a short description of what each of the scripts from the "scripts" directory does below:


        figures_h2_2_3_4.R generates fragments of Fig. 2 and 3, and Fig. 4; figures_h2_S1_S2 generates Supplementary Fig. 1 and 2: Heritability in the post-Soviet and Soviet groups and the groups additionally divided by the wave of participation
                          
        figures_h2_decades_2_3_S3.R generates fragments of Fig. 2 and 3, and Supplementary Fig. 3: Heritability across decade-long birth year cohorts 
        
        plot_h2_origSample generates Supplementary Figure 13: Heritability in the post-Soviet and Soviet groups in the original sample
        

        pgsR2.R - All the calculations with R2:
                  1. Prepare the data for the analysis 
                  2. Main R2 analysis 
                  3. Weighting 
                  4. Matching distributions 
                  5. Original sample 


        figure_R2_5.R generates Fig. 5; figures_R2_S4_S7.R generates Supplementary Fig. 4 and 7: Trait variance explained by PGS in the post-Soviet and Soviet groups
        
        figures_R2_decades_5_S5_S6.R generates a fragment of Fig. 5,  Supplementary Fig. 5 and 6: Variance explained by PGS in the birth cohorts by half-decade 
        
        plot_R2_origSample.R generates Supplementary Fig. 14 and 15: Variance explained by PGS in the post-Soviet and Soviet groups in the original sample


        plot R2_matching.R generates Fig. 6 and Supplementary Fig. 8: Variance of EA explained by PGSEA in subsamples from the groups (defined by wave participation and era) 
                            selected based on the distribution of EA in another subcohort

        plot_R2_weighting.R generates Fig. 7: Trait variance explained by PGS in weighted subsamples


        regressionInteractions.R analyses GxE through linear regression, generates Fig. 8 and Supplementary Fig. 9-11
