library(reshape2)
library(reshape)
library(tidyverse)
library(MASS)
library(mice)


# Function for LOCF
coalesce <- function(...) {
  apply(cbind(...), 1, function(x) {
    x[which(!is.na(x))[1]]
  })
}

# Function for subset
subset_trt <- function(data, value) {
  data1 <- data %>% dplyr::filter(Trt == value)
  return(data1)
}

# Function for Variance of ratio using Delta Method
delta_var <-
  function(baseline, postbaseline, type = 0, base, post, basev, postv, covar) {
    if (type == 0) { 
      dvar <- (((var(postbaseline) / (mean(baseline) ** 2 )) +
                  ((var(baseline) * (mean(postbaseline) ** 2)) / (mean(baseline) ** 4)) -
                  (2 * mean(postbaseline) * cov(baseline, postbaseline) /
                     (mean(baseline)**3)
                  )) * (100 ** 2)) / ss
    } else if (type == 1) {
      dvar <- (((mean(postv) / (mean(base) ** 2)) +
                  ((mean(basev) * (mean(post) ** 2)) / (mean(base) ** 4)) -
                  (2 * mean(post) * mean(covar) / (mean(base) ** 3)
                  )) * (100 ** 2)) / ss
    }
    return(dvar)
  }
# Function for Delta estimate
delta_mean_mi <- function() {
  ((mean(delta_estimatet$post_t) / mean(delta_estimatet$base_t) ) -
     (mean(delta_estimatec$post_c) / mean(delta_estimatec$base_c))
  ) * 100
}
# Function for Delta estimate non MI
delta_mean <- function(data, postb) {
  ((mean(eval(substitute(postb), subset_trt(data, 1)))/
      mean(subset_trt(data, 1)$baseline)) -
     (mean(eval(substitute(postb), subset_trt(data, 0)))/
        mean(subset_trt(data, 0)$baseline))
  ) * 100
}
# Function for Delta SD MI
delta_se_mi <- function() {
  sqrt(
    delta_var(type=1, base=delta_estimatet$base_t, 
              post=delta_estimatet$post_t, basev=delta_estimatet$var_base_t,
              postv=delta_estimatet$var_post_t, covar=delta_estimatet$cov_t
    ) +
      delta_var(type=1, base=delta_estimatec$base_c, 
                post=delta_estimatec$post_c, basev=delta_estimatec$var_base_c,
                postv=delta_estimatec$var_post_c, covar=delta_estimatec$cov_c
      )
  )
}
# Function for Delta SD non MI
delta_se <- function(data, postb) {
  sqrt(
    delta_var(subset_trt(data, 1)$baseline, 
              eval(substitute(postb), subset_trt(data, 1))) +
      delta_var(subset_trt(data, 0)$baseline,
                eval(substitute(postb), subset_trt(data, 0))
      )
  )
}

# Test Drug mean 
mu1 <- c(0.7223, 0.72454, 0.724, 0.72147)
# Test var-cov matrix
sigma1 <- matrix(
  c(
    0.006035, 0.005872, 0.006363, 0.006749,
    0.005872, 0.005979, 0.00636, 0.006746,
    0.006363, 0.00636, 0.007052, 0.007406,
    0.006749, 0.006746, 0.007406, 0.008171
  ),
  ncol = 4, byrow = TRUE
)

# Comparator Drug mean
mu2 <- c(0.73455, 0.72935, 0.73068, 0.72844)
# Comparator var-cov matrix
sigma2 <- matrix(
  c(
    0.009362, 0.009508, 0.009203, 0.009055,
    0.009508, 0.009971, 0.00949, 0.009376,
    0.009203, 0.00949, 0.009355, 0.009207,
    0.009055, 0.009376, 0.009207, 0.009285
  ),
  ncol = 4, byrow = TRUE
)
# expected diff
expected1 <-
  (((0.72147 - 0.7223) / 0.7223) - ((0.72844 - 0.73455) / 0.73455)) * 100
#function to generate and analyze data
simulation_analysis <- function(sim, rho1 = 0.8) {
  for (i in 1:sim) {
    #Generate data for test drug
    set.seed(123789 + i)
    bvn1 <- data.frame(mvrnorm(ss,
                               mu = mu1,
                               Sigma = sigma1))
    colnames(bvn1) <- c("X1", "X4", "X3", "X2")
    bvn1 <- bvn1 %>% dplyr::select(X1, X2)
    bvn1$pch <- ((bvn1$X2 / bvn1$X1) - 1) * 100
    median_t <- median(bvn1$pch)
    bvn1$rand10 <-
      ifelse(bvn1$pch < median_t, rbinom(ss, 1, 0.8), 1)
    bvn1$rand15 <-
      ifelse(bvn1$pch < median_t, rbinom(ss, 1, 0.7), 1)
    bvn1$rand20 <-
      ifelse(bvn1$pch < median_t, rbinom(ss, 1, 0.6), 1)
    bvn1$rand25 <-
      ifelse(bvn1$pch < median_t, rbinom(ss, 1, 0.5), 1)
    bvn1$pchmis10 <-
      ifelse(bvn1$pch < median_t & bvn1$rand10 == 0, NA, bvn1$pch)
    bvn1$pchmis15 <-
      ifelse(bvn1$pch < median_t & bvn1$rand15 == 0, NA, bvn1$pch)
    bvn1$pchmis20 <-
      ifelse(bvn1$pch < median_t & bvn1$rand20 == 0, NA, bvn1$pch)
    bvn1$pchmis25 <-
      ifelse(bvn1$pch < median_t & bvn1$rand25 == 0, NA, bvn1$pch)
    bvn1$X2_10 <-
      ifelse(bvn1$pch < median_t & bvn1$rand10 == 0, NA, bvn1$X2)
    bvn1$X2_15 <-
      ifelse(bvn1$pch < median_t & bvn1$rand15 == 0, NA, bvn1$X2)
    bvn1$X2_20 <-
      ifelse(bvn1$pch < median_t & bvn1$rand20 == 0, NA, bvn1$X2)
    bvn1$X2_25 <-
      ifelse(bvn1$pch < median_t & bvn1$rand25 == 0, NA, bvn1$X2)
    #Generate data for comparator drug
    set.seed(123789 + i)
    bvn2 <- data.frame(mvrnorm(ss,
                               mu = mu2,
                               Sigma = sigma2))
    colnames(bvn2) <- c("X1", "X4", "X3", "X2")
    bvn2 <- bvn2 %>% dplyr::select(X1, X2)
    bvn2$pch <- ((bvn2$X2 / bvn2$X1) - 1) * 100
    median_c <- median(bvn2$pch)
    bvn2$rand10 <-
      ifelse(bvn2$pch < median_c, rbinom(ss, 1, 0.8), 1)
    bvn2$rand15 <-
      ifelse(bvn2$pch < median_c, rbinom(ss, 1, 0.7), 1)
    bvn2$rand20 <-
      ifelse(bvn2$pch < median_c, rbinom(ss, 1, 0.6), 1)
    bvn2$rand25 <-
      ifelse(bvn2$pch < median_c, rbinom(ss, 1, 0.5), 1)
    bvn2$pchmis10 <-
      ifelse(bvn2$pch < median_c & bvn2$rand10 == 0, NA, bvn2$pch)
    bvn2$pchmis15 <-
      ifelse(bvn2$pch < median_c & bvn2$rand15 == 0, NA, bvn2$pch)
    bvn2$pchmis20 <-
      ifelse(bvn2$pch < median_c & bvn2$rand20 == 0, NA, bvn2$pch)
    bvn2$pchmis25 <-
      ifelse(bvn2$pch < median_c & bvn2$rand25 == 0, NA, bvn2$pch)
    bvn2$X2_10 <-
      ifelse(bvn2$pch < median_c & bvn2$rand10 == 0, NA, bvn2$X2)
    bvn2$X2_15 <-
      ifelse(bvn2$pch < median_c & bvn2$rand15 == 0, NA, bvn2$X2)
    bvn2$X2_20 <-
      ifelse(bvn2$pch < median_c & bvn2$rand20 == 0, NA, bvn2$X2)
    bvn2$X2_25 <-
      ifelse(bvn2$pch < median_c & bvn2$rand25 == 0, NA, bvn2$X2)
    #combine test drug data and comparator drug data in one dataset
    predata0 <-
      rbind(cbind(bvn1, Treat = 1), cbind(bvn2, Treat = 0))
    #generate missing data
    for (j in c(20, 30, 40, 50)) {
      if (j == 20) {
        predata <- predata0 %>% dplyr::select(X1, pchmis10, X2_10, Treat)
        colnames(predata) <-
          c("baseline", "response", "postbase", "Trt")
      }
      if (j == 30) {
        predata <- predata0 %>% dplyr::select(X1, pchmis15, X2_15, Treat)
        colnames(predata) <-
          c("baseline", "response", "postbase", "Trt")
      }
      if (j == 40) {
        predata <- predata0 %>% dplyr::select(X1, pchmis20, X2_20, Treat)
        colnames(predata) <-
          c("baseline", "response", "postbase", "Trt")
      }
      if (j == 50) {
        predata <- predata0 %>% dplyr::select(X1, pchmis25, X2_25, Treat)
        colnames(predata) <-
          c("baseline", "response", "postbase", "Trt")
      }
      # multiple imputation
      imp <-
        mice(
          data = predata,
          m = 10,
          method = "pmm",
          maxit = 50,
          print = FALSE
        )
      # First, turn the datasets into long format
      imp_long <-
        mice::complete(imp, action = "long", include = TRUE)
      # Convert treatment variable to factor
      imp_long$Trt <- with(imp_long,
                           as.factor(imp_long$Trt))
      # Convert back to mids type - mice can work with this type
      imp_long_mids <- as.mids(imp_long)
      # Regression
      fitimp <- with(imp_long_mids,
                     lm(response ~ Trt))
      delta_estimatet <- imp_long %>%
        filter(.imp > 0 & Trt == 1) %>%
        group_by(factor(.imp)) %>%
        summarize(
          post_t = mean(postbase),
          base_t = mean(baseline),
          var_post_t = var(postbase),
          var_base_t = var(baseline),
          cov_t = cov(postbase, baseline)
        )
      delta_estimatec <- imp_long %>%
        filter(.imp > 0 & Trt == 0) %>%
        group_by(factor(.imp)) %>%
        summarize(
          post_c = mean(postbase),
          base_c = mean(baseline),
          var_post_c = var(postbase),
          var_base_c = var(baseline),
          cov_c = cov(postbase, baseline)
        )
      #Accumulating data from different simulations
      if (i < 2) {
        if (j == 20) {
          flagdat1_20 <- cbind(
            summary(pool(fitimp))[2, 2], summary(pool(fitimp))[2, 3],
            delta_mean_mi(), delta_se_mi()
          )
          colnames(flagdat1_20) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        } else if (j == 30) {
          flagdat1_30 <- cbind(
            summary(pool(fitimp))[2, 2], summary(pool(fitimp))[2, 3],
            delta_mean_mi(), delta_se_mi()
          )
          colnames(flagdat1_30) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        } else if (j == 40) {
          flagdat1_40 <- cbind(
            summary(pool(fitimp))[2, 2], summary(pool(fitimp))[2, 3],
            delta_mean_mi(), delta_se_mi()
          )
          colnames(flagdat1_40) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        } else if (j == 50) {
          flagdat1_50 <- cbind(
            summary(pool(fitimp))[2, 2], summary(pool(fitimp))[2, 3],
            delta_mean_mi(), delta_se_mi()
          )
          colnames(flagdat1_50) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        }
      } else {
        if (j == 20) {
          flagdat1_20 <- 
            rbind(flagdat1_20,
                  cbind(summary(pool(fitimp))[2, 2],summary(pool(fitimp))[2, 3],
                        delta_mean_mi(), delta_se_mi()
                  ))
        } else if (j == 30) {
          flagdat1_30 <- 
            rbind(flagdat1_30,
                  cbind(summary(pool(fitimp))[2, 2],summary(pool(fitimp))[2, 3],
                        delta_mean_mi(), delta_se_mi()
                  ))
        } else if (j == 40) {
          flagdat1_40 <- 
            rbind(flagdat1_40,
                  cbind(summary(pool(fitimp))[2, 2],summary(pool(fitimp))[2, 3],
                        delta_mean_mi(), delta_se_mi()
                  ))
        } else if (j == 50) {
          flagdat1_50 <- 
            rbind(flagdat1_50,
                  cbind(summary(pool(fitimp))[2, 2],summary(pool(fitimp))[2, 3],
                        delta_mean_mi(), delta_se_mi()
                  ))
        }
      }
      # completers set
      predata1 <- predata %>% filter(!is.na(response))
      fit2 <- aov(response ~ Trt, predata1)
      #Accumuating results from different simulations
      if (i < 2) {
        if (j == 20) {
          flagdat2_20 <- cbind(
            coef(summary.lm(fit2))[2, 1], coef(summary.lm(fit2))[2, 2],
            delta_mean(predata1, postbase), delta_se(predata1, postbase)
          )
          colnames(flagdat2_20) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        } else if (j == 30) {
          flagdat2_30 <- cbind(
            coef(summary.lm(fit2))[2, 1], coef(summary.lm(fit2))[2, 2],
            delta_mean(predata1, postbase), delta_se(predata1, postbase)
          )
          colnames(flagdat2_30) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        } else if (j == 40) {
          flagdat2_40 <- cbind(
            coef(summary.lm(fit2))[2, 1], coef(summary.lm(fit2))[2, 2],
            delta_mean(predata1, postbase), delta_se(predata1, postbase)
          )
          colnames(flagdat2_40) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        } else if (j == 50) {
          flagdat2_50 <- cbind(
            coef(summary.lm(fit2))[2, 1], coef(summary.lm(fit2))[2, 2],
            delta_mean(predata1, postbase), delta_se(predata1, postbase)
          )
          colnames(flagdat2_50) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        }
      } else {
        if (j == 20) {
          flagdat2_20 <- rbind(flagdat2_20, cbind(
            coef(summary.lm(fit2))[2, 1], coef(summary.lm(fit2))[2, 2],
            delta_mean(predata1, postbase), delta_se(predata1, postbase)
          ))
        } else if (j == 30) {
          flagdat2_30 <- rbind(flagdat2_30, cbind(
            coef(summary.lm(fit2))[2, 1], coef(summary.lm(fit2))[2, 2],
            delta_mean(predata1, postbase), delta_se(predata1, postbase)
          ))
        } else if (j == 40) {
          flagdat2_40 <- rbind(flagdat2_40, cbind(
            coef(summary.lm(fit2))[2, 1], coef(summary.lm(fit2))[2, 2],
            delta_mean(predata1, postbase), delta_se(predata1, postbase)
          ))
        } else if (j == 50) {
          flagdat2_50 <- rbind(flagdat2_50, cbind(
            coef(summary.lm(fit2))[2, 1], coef(summary.lm(fit2))[2, 2],
            delta_mean(predata1, postbase), delta_se(predata1, postbase)
          ))
        }
      }
      # LOCF
      predata2 <- predata
      predata2$response1 <- coalesce(predata$response, 0)
      predata2$postbase1 <- coalesce(predata$postbase, predata$baseline)
      fit3 <- aov(response1 ~ Trt, predata2)
      if (i < 2) {
        if (j == 20) {
          flagdat3_20 <- cbind(
            coef(summary.lm(fit3))[2, 1], coef(summary.lm(fit3))[2, 2],
            delta_mean(predata2, postbase1), delta_se(predata2, postbase1)
          )
          colnames(flagdat3_20) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        } else if (j == 30) {
          flagdat3_30 <- cbind(
            coef(summary.lm(fit3))[2, 1], coef(summary.lm(fit3))[2, 2],
            delta_mean(predata2, postbase1), delta_se(predata2, postbase1)
          )
          colnames(flagdat3_30) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        } else if (j == 40) {
          flagdat3_40 <- cbind(
            coef(summary.lm(fit3))[2, 1], coef(summary.lm(fit3))[2, 2],
            delta_mean(predata2, postbase1), delta_se(predata2, postbase1)
          )
          colnames(flagdat3_40) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        } else if (j == 50) {
          flagdat3_50 <- cbind(
            coef(summary.lm(fit3))[2, 1], coef(summary.lm(fit3))[2, 2],
            delta_mean(predata2, postbase1), delta_se(predata2, postbase1)
          )
          colnames(flagdat3_50) <-
            c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
        }
      } else {
        if (j == 20) {
          flagdat3_20 <- rbind(flagdat3_20, cbind(
            coef(summary.lm(fit3))[2, 1], coef(summary.lm(fit3))[2, 2],
            delta_mean(predata2, postbase1), delta_se(predata2, postbase1)
          ))
        } else if (j == 30) {
          flagdat3_30 <- rbind(flagdat3_30, cbind(
            coef(summary.lm(fit3))[2, 1], coef(summary.lm(fit3))[2, 2],
            delta_mean(predata2, postbase1), delta_se(predata2, postbase1)
          ))
        } else if (j == 40) {
          flagdat3_40 <- rbind(flagdat3_40, cbind(
            coef(summary.lm(fit3))[2, 1], coef(summary.lm(fit3))[2, 2],
            delta_mean(predata2, postbase1), delta_se(predata2, postbase1)
          ))
        } else if (j == 50) {
          flagdat3_50 <- rbind(flagdat3_50, cbind(
            coef(summary.lm(fit3))[2, 1], coef(summary.lm(fit3))[2, 2],
            delta_mean(predata2, postbase1), delta_se(predata2, postbase1)
          ))
        }
      }
    }
    # whole data
    predata_orig <- predata0 %>% dplyr::select(X1, X2, pch, Treat)
    colnames(predata_orig) <-
      c("baseline", "postbase", "response", "Trt")
    fit1 <- aov(response ~ Trt, predata_orig)
    if (i < 2) {
      flagdat <- cbind(
        coef(summary.lm(fit1))[2, 1], coef(summary.lm(fit1))[2, 2],
        delta_mean(predata_orig, postbase), delta_se(predata_orig, postbase)
      )
      colnames(flagdat) <-
        c("Estimate", "Std_Error", "D_estimate", "D_Std_Error")
    } else {
      flagdat <- rbind(flagdat, cbind(
        coef(summary.lm(fit1))[2, 1], coef(summary.lm(fit1))[2, 2],
        delta_mean(predata_orig, postbase), delta_se(predata_orig, postbase)
      ))
    }
  }
  # record data in a single list for output
  return(
    list(
      # Actual complete data estimates
      Orig_mean = mean(data.frame(flagdat)$Estimate),
      Orig_stderr = sqrt(mean(data.frame(flagdat)$Std_Error ** 2)),
      Orig_dmean = mean(data.frame(flagdat)$D_estimate),
      Orig_dstderr = sqrt(mean(data.frame(flagdat)$D_Std_Error ** 2)),
      # Completers set estimates ignoring the missing
      COMP10_mean = mean(data.frame(flagdat2_20)$Estimate),
      COMP10_stderr = sqrt(mean(data.frame(flagdat2_20)$Std_Error ** 2)),
      COMP10_dmean = mean(data.frame(flagdat2_20)$D_estimate),
      COMP10_dstderr = sqrt(mean(data.frame(flagdat2_20)$D_Std_Error ** 2)),
      COMP15_mean = mean(data.frame(flagdat2_30)$Estimate),
      COMP15_stderr = sqrt(mean(data.frame(flagdat2_30)$Std_Error ** 2)),
      COMP15_dmean = mean(data.frame(flagdat2_30)$D_estimate),
      COMP15_dstderr = sqrt(mean(data.frame(flagdat2_30)$D_Std_Error ** 2)),
      COMP20_mean = mean(data.frame(flagdat2_40)$Estimate),
      COMP20_stderr = sqrt(mean(data.frame(flagdat2_40)$Std_Error ** 2)),
      COMP20_dmean = mean(data.frame(flagdat2_40)$D_estimate),
      COMP20_dstderr = sqrt(mean(data.frame(flagdat2_40)$D_Std_Error ** 2)),
      COMP25_mean = mean(data.frame(flagdat2_50)$Estimate),
      COMP25_stderr = sqrt(mean(data.frame(flagdat2_50)$Std_Error ** 2)),
      COMP25_dmean = mean(data.frame(flagdat2_50)$D_estimate),
      COMP25_dstderr = sqrt(mean(data.frame(flagdat2_50)$D_Std_Error ** 2)),
      # Estimates imputing the missing with LOCF method
      LOCF10_mean = mean(data.frame(flagdat3_20)$Estimate),
      LOCF10_stderr = sqrt(mean(data.frame(flagdat3_20)$Std_Error ** 2)),
      LOCF10_dmean = mean(data.frame(flagdat3_20)$D_estimate),
      LOCF10_dstderr = sqrt(mean(data.frame(flagdat3_20)$D_Std_Error ** 2)),
      LOCF15_mean = mean(data.frame(flagdat3_30)$Estimate),
      LOCF15_stderr = sqrt(mean(data.frame(flagdat3_30)$Std_Error ** 2)),
      LOCF15_dmean = mean(data.frame(flagdat3_30)$D_estimate),
      LOCF15_dstderr = sqrt(mean(data.frame(flagdat3_30)$D_Std_Error ** 2)),
      LOCF20_mean = mean(data.frame(flagdat3_40)$Estimate),
      LOCF20_stderr = sqrt(mean(data.frame(flagdat3_40)$Std_Error ** 2)),
      LOCF20_dmean = mean(data.frame(flagdat3_40)$D_estimate),
      LOCF20_dstderr = sqrt(mean(data.frame(flagdat3_40)$D_Std_Error ** 2)),
      LOCF25_mean = mean(data.frame(flagdat3_50)$Estimate),
      LOCF25_stderr = sqrt(mean(data.frame(flagdat3_50)$Std_Error ** 2)),
      LOCF25_dmean = mean(data.frame(flagdat3_50)$D_estimate),
      LOCF25_dstderr = sqrt(mean(data.frame(flagdat3_50)$D_Std_Error ** 2)),
      # Estimates imputing the missing with MI method
      MI10_mean = mean(data.frame(flagdat1_20)$Estimate),
      MI10_stderr = sqrt(mean(data.frame(flagdat1_20)$Std_Error ** 2)),
      MI10_dmean = mean(data.frame(flagdat1_20)$D_estimate),
      MI10_dstderr = sqrt(mean(data.frame(flagdat1_20)$D_Std_Error ** 2)),
      MI15_mean = mean(data.frame(flagdat1_30)$Estimate),
      MI15_stderr = sqrt(mean(data.frame(flagdat1_30)$Std_Error ** 2)),
      MI15_dmean = mean(data.frame(flagdat1_30)$D_estimate),
      MI15_dstderr = sqrt(mean(data.frame(flagdat1_30)$D_Std_Error ** 2)),
      MI20_mean = mean(data.frame(flagdat1_40)$Estimate),
      MI20_stderr = sqrt(mean(data.frame(flagdat1_40)$Std_Error ** 2)),
      MI20_dmean = mean(data.frame(flagdat1_40)$D_estimate),
      MI20_dstderr = sqrt(mean(data.frame(flagdat1_40)$D_Std_Error ** 2)),
      MI25_mean = mean(data.frame(flagdat1_50)$Estimate),
      MI25_stderr = sqrt(mean(data.frame(flagdat1_50)$Std_Error ** 2)),
      MI25_dmean = mean(data.frame(flagdat1_50)$D_estimate),
      MI25_dstderr = sqrt(mean(data.frame(flagdat1_50)$D_Std_Error ** 2))
    )
  )
}
ss <- 400
out_matrix1 <- data.frame(simulation_analysis(sim = 1000))
output1 <- data.frame(t(out_matrix1))

output21 <- data.frame(cbind(output1, rownames(output1)))
colnames(output21) <- c("value", "stat")
row.names(output21) <- NULL

output22 <-
  data.frame(do.call("rbind", strsplit(as.character(output21$stat), "_",
                                       fixed = TRUE)))
output23 <- cbind(output22, value = output21$value)
cast_data1 <- cast(output23, X1 ~ X2)
# cast_data1$lsmean <- -1 * cast_data1$lsmean
cast_data1$lci <- cast_data1$mean - 1.96*cast_data1$stderr
cast_data1$uci <- cast_data1$mean + 1.96*cast_data1$stderr
cast_data1$dlci <- cast_data1$dmean - 1.96*cast_data1$dstderr
cast_data1$duci <- cast_data1$dmean + 1.96*cast_data1$dstderr
cast_data1$meanse <- paste0(round(cast_data1$mean,3), " (", round(cast_data1$stderr,3),")")
cast_data1$ci <- paste0("(", round(cast_data1$lci,3), ",",round(cast_data1$uci,3),")")
cast_data1$dmeanse <- paste0(round(cast_data1$dmean,3), " (", round(cast_data1$dstderr,3),")")
cast_data1$dci <- paste0("(", round(cast_data1$dlci,3), ",",round(cast_data1$duci,3),")")
cast_data1$lbias <- abs(cast_data1$mean - expected1)
cast_data1$dbias <- abs(cast_data1$dmean - expected1)
cast_data1
