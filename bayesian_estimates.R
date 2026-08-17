library(data.table)

bayes_normal <- function(estimate,
                         se,
                         prior_mean = 0,
                         prior_sd = 0.10) {
  
  # variances
  prior_var <- prior_sd^2
  likelihood_var <- se^2
  
  # posterior variance
  post_var <- 1 / (1/prior_var + 1/likelihood_var)
  post_sd <- sqrt(post_var)
  
  # posterior mean
  post_mean <- post_var * (
    prior_mean/prior_var +
      estimate/likelihood_var
  )
  
  # credible interval
  ci_lower <- post_mean - 1.96 * post_sd
  ci_upper <- post_mean + 1.96 * post_sd
  
  # probability of positive effect
  prob_positive <- pnorm(0,
                         mean = post_mean,
                         sd = post_sd,
                         lower.tail = FALSE)
  
  data.frame(
    estimate = estimate,
    SE = se,
    posterior_mean = post_mean,
    posterior_SD = post_sd,
    CI95_lower = ci_lower,
    CI95_upper = ci_upper,
    P_positive = prob_positive
  )
}

###########
# h2 ######
###########

h2_results <- fread("~/EA_heritability/figures/paper/revision/h2_estimates_main_pedigree.tsv")

h2_res <- h2_results[Trait=="EA" & Cutoff=="15", ]

h2_w12_delta <- h2_res[Wave=="All" & Era=="PS", h2] - h2_res[Wave=="All" & Era=="S", h2]
h2_w12_se <- sqrt(h2_res[Wave=="All" & Era=="PS", se]^2 + h2_res[Wave=="All" & Era=="S", se]^2)
bayes_h2 <- bayes_normal(estimate = h2_w12_delta, se = h2_w12_se,
                         prior_mean = 0, prior_sd = 0.05)

h2_w2_delta <- h2_res[Wave=="2" & Era=="PS", h2] - h2_res[Wave=="2" & Era=="S", h2]
h2_w2_se <- sqrt(h2_res[Wave=="2" & Era=="PS", se]^2 + h2_res[Wave=="2" & Era=="S", se]^2)
bayes_normal(estimate = h2_w2_delta, se = h2_w2_se,
             prior_mean = 0.19, prior_sd = 0.146)


###########
# R2 ######
###########

r2_results <- fread("~/EA_heritability/results/revision/r2_adj_EduYears_1000.tsv")

r2_w12_dif <- r2_results[, ps_r2_inc] - r2_results[, s_r2_inc]
r2_w12_delta <- mean(r2_w12_dif)
r2_w12_se <- sd(r2_w12_dif)
bayes_r2 <- bayes_normal(estimate = r2_w12_delta, se = h2_w12_se,
                         prior_mean = 0, prior_sd = 0.05)

r2_w2_dif <- r2_results[, ps_r2_inc_p2] - r2_results[, s_r2_inc_p2]
r2_w2_delta <- mean(r2_w2_dif)
r2_w2_se <- sd(r2_w2_dif)
bayes_normal(estimate = r2_w2_delta, se = h2_w2_se,
             prior_mean = 0.015, prior_sd = 0.03)


###########
# PGSxEra #
###########

regr_results <- fread("~/EA_heritability/figures/paper/revision/regres_interact_EA_withYoB_15_all.tsv")

regr_beta <- regr_results[V1=="I(PRS * era)", Estimate]
regr_se <- regr_results[V1=="I(PRS * era)", `Std. Error`]
bayes_int_withYoB <- bayes_normal(estimate = regr_beta, se = regr_se,
                                  prior_mean = 0, prior_sd = 0.036) # based on SE for PGS beta

regr_results <- fread("~/EA_heritability/figures/paper/revision/regres_interact_EA_noYoB_15_all.tsv")

regr_beta <- regr_results[V1=="I(PRS * era)", Estimate]
regr_se <- regr_results[V1=="I(PRS * era)", `Std. Error`]
bayes_int_noYoB <- bayes_normal(estimate = regr_beta, se = regr_se,
             prior_mean = 0, prior_sd = 0.032) # based on SE for PGS beta


bayesian_estimates <- rbind(bayes_h2,
                            bayes_r2,
                            bayes_int_withYoB,
                            bayes_int_noYoB)
write.table(bayesian_estimates, "~/EA_heritability/figures/paper/revision/bayesian_estimates.tsv",
            row.names = F, quote = F, sep = "\t")
