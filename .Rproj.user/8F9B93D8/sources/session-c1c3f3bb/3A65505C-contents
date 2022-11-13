# install.packages("tmvtnorm")
# install.packages("ggplot2")
# install.packages("moments")
# install.packages("tidyverse")
library(tidyverse)
library(moments)
library(tmvtnorm)
library(ggplot2)
library(MASS)
library(dplyr)

# Sample size required per arm for 2 treatment comparision
power <- 0.85
alpha <- 0.05
# effect_size <-  0.787/4.838  # (mu1 - mu2)/sigma
# effect_size <-  1.57/3.17 #1.57/3.1688  # (mu1 - mu2)/sigma
effect_size <-  0.37/1.10565 #1.57/3.1688  # (mu1 - mu2)/sigma
Npergrp <- (2*((qnorm(1 - (alpha/2)) + qnorm(power))**2))/effect_size**2
ceiling(Npergrp)


# size <- c(22, 35, 45, 55, 65, 75, 85, 100, 150, 233, 250, 275, 300, 350, 400, 450, 500, 680, 1000)
size <- 161 #74 #243 #195 #120 #680
# rho1 <- c(-0.7, -0.5, -0.3, 0.3, 0.5, 0.7)
rho1 <- 0.7
s1 <- 1
s2 <- 1
s1*s2- s1*s2*rho1*s1*s2*rho1
base <- 0.7
t_pbase <- 1.029 #2.5 #0.3174
c_pbase <- 0.77 #1.4 #-0.2201
expected <- (((t_pbase - base)/base) - ((c_pbase - base)/base))*100
# for (rho in rho1){
v1 <- NULL
samplemean <- data.frame(v1)
for (ss in size){
  for(i in 1:5000){
    # set.seed(123789+i)
    # test <- data.frame(cbind(rtmvnorm(n=800, mean=c(compeffect, treateffect), 
    #                                   sigma=matrix(c(s1, s1*s2*rho, s1*s2*rho, s2),2 )
    #                                   # , lower=c(-0.5, -0.3), upper=c(1.5, 1.7), algorithm="rejection"
    # ), rho))
    bvn1 <- data.frame(mvrnorm(ss, mu = c(base, t_pbase), Sigma = matrix(c(s1, s1*s2*rho1, s1*s2*rho1, s2),2 ))) # from MASS package
    colnames(bvn1) <- c("X1","X2")
    bvn1$pch <- ((bvn1$X2/bvn1$X1)-1)*100
    Median_t<-median(bvn1$pch);
    bvn1$rand10 <- ifelse(bvn1$pch < Median_t, rbinom(ss,1,0.8), 1)
    bvn1$rand15 <- ifelse(bvn1$pch < Median_t, rbinom(ss,1,0.7), 1)
    bvn1$rand20 <- ifelse(bvn1$pch < Median_t, rbinom(ss,1,0.6), 1)
    bvn1$rand25 <- ifelse(bvn1$pch < Median_t, rbinom(ss,1,0.5), 1)
    bvn1$pchmis10 <- ifelse (bvn1$pch < Median_t & bvn1$rand10 == 0,NA, bvn1$pch);
    bvn1$pchmis15 <- ifelse (bvn1$pch < Median_t & bvn1$rand15 == 0,NA, bvn1$pch);
    bvn1$pchmis20 <- ifelse (bvn1$pch < Median_t & bvn1$rand20 == 0,NA, bvn1$pch);
    bvn1$pchmis25 <- ifelse (bvn1$pch < Median_t & bvn1$rand25 == 0,NA, bvn1$pch);
    bvn1$X2_10 <- ifelse (bvn1$pch < Median_t & bvn1$rand10 == 0,NA, bvn1$X2);
    bvn1$X2_15 <- ifelse (bvn1$pch < Median_t & bvn1$rand15 == 0,NA, bvn1$X2);
    bvn1$X2_20 <- ifelse (bvn1$pch < Median_t & bvn1$rand20 == 0,NA, bvn1$X2);
    bvn1$X2_25 <- ifelse (bvn1$pch < Median_t & bvn1$rand25 == 0,NA, bvn1$X2);
    
    # set.seed(123789+i)
    bvn2 <- data.frame(mvrnorm(ss, mu = c(base, c_pbase), Sigma = matrix(c(s1, s1*s2*rho1, s1*s2*rho1, s2),2 ))) # from MASS package
    colnames(bvn2) <- c("X1","X2")
    bvn2$pch <- ((bvn2$X2/bvn2$X1)-1)*100
    Median_c<-median(bvn2$pch);
    bvn2$rand10 <- ifelse(bvn2$pch < Median_c, rbinom(ss,1,0.8), 1)
    bvn2$rand15 <- ifelse(bvn2$pch < Median_c, rbinom(ss,1,0.7), 1)
    bvn2$rand20 <- ifelse(bvn2$pch < Median_c, rbinom(ss,1,0.6), 1)
    bvn2$rand25 <- ifelse(bvn2$pch < Median_c, rbinom(ss,1,0.5), 1)
    bvn2$pchmis10 <- ifelse (bvn2$pch < Median_c & bvn2$rand10 == 0,NA, bvn2$pch);
    bvn2$pchmis15 <- ifelse (bvn2$pch < Median_c & bvn2$rand15 == 0,NA, bvn2$pch);
    bvn2$pchmis20 <- ifelse (bvn2$pch < Median_c & bvn2$rand20 == 0,NA, bvn2$pch);
    bvn2$pchmis25 <- ifelse (bvn2$pch < Median_c & bvn2$rand25 == 0,NA, bvn2$pch);
    bvn2$X2_10 <- ifelse (bvn2$pch < Median_c & bvn2$rand10 == 0,NA, bvn2$X2);
    bvn2$X2_15 <- ifelse (bvn2$pch < Median_c & bvn2$rand15 == 0,NA, bvn2$X2);
    bvn2$X2_20 <- ifelse (bvn2$pch < Median_c & bvn2$rand20 == 0,NA, bvn2$X2);
    bvn2$X2_25 <- ifelse (bvn2$pch < Median_c & bvn2$rand25 == 0,NA, bvn2$X2);
    # sum(is.na(bvn1$pchmis25))
    
    mean <- cbind(
      # summary from comparator drug
      meant=mean(bvn1$pch), sqrt(var(bvn1$pch)), ss,mean(bvn1$X1),mean(bvn1$X2), 
      deltat=((mean(bvn1$X2)/mean(bvn1$X1))-1)*100, 
      meant_10 = mean(bvn1$pchmis10, na.rm=TRUE),
      meant_15 = mean(bvn1$pchmis15, na.rm=TRUE),
      meant_20 = mean(bvn1$pchmis20, na.rm=TRUE),
      meant_25 = mean(bvn1$pchmis25, na.rm=TRUE),
      # summary from comparator drug
      meanc=mean(bvn2$pch), sqrt(var(bvn2$pch)), mean(bvn2$X1),mean(bvn2$X2), 
      deltac=((mean(bvn2$X2)/mean(bvn2$X1))-1)*100, 
      meanc_10 = mean(bvn2$pchmis10, na.rm=TRUE),
      meanc_15 = mean(bvn2$pchmis15, na.rm=TRUE),
      meanc_20 = mean(bvn2$pchmis20, na.rm=TRUE),
      meanc_25 = mean(bvn2$pchmis25, na.rm=TRUE),
      #treatment difference
      effect=mean(bvn1$pch) - mean(bvn2$pch), 
      effect_10=mean(bvn1$pchmis10, na.rm=TRUE) - mean(bvn2$pchmis10, na.rm=TRUE),
    effect_15=mean(bvn1$pchmis15, na.rm=TRUE) - mean(bvn2$pchmis15, na.rm=TRUE),
    effect_20=mean(bvn1$pchmis20, na.rm=TRUE) - mean(bvn2$pchmis20, na.rm=TRUE),
    effect_25=mean(bvn1$pchmis25, na.rm=TRUE) - mean(bvn2$pchmis25, na.rm=TRUE),
      # Delta method
    d_effect=(((mean(bvn1$X2)/mean(bvn1$X1))-1) - ((mean(bvn2$X2)/mean(bvn2$X1))-1))*100,
    d_effect_10=(((mean(bvn1$X2_10, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_10, na.rm=TRUE)/mean(bvn2$X1))-1))*100,
    d_effect_15=(((mean(bvn1$X2_15, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_15, na.rm=TRUE)/mean(bvn2$X1))-1))*100,
    d_effect_20=(((mean(bvn1$X2_20, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_20, na.rm=TRUE)/mean(bvn2$X1))-1))*100,
    d_effect_25=(((mean(bvn1$X2_25, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_25, na.rm=TRUE)/mean(bvn2$X1))-1))*100
    )
    samplemean <- rbind(samplemean, mean)
  }
  samplemean <- samplemean %>% slice(501:5000)
}

samplemean %>% 
  group_by(factor(ss)) %>% 
  summarize(mean = mean(effect))

samplemean %>%
  group_by(factor(ss)) %>% 
  summarize(mean = mean(d_effect))

samplemean %>% 
  group_by(factor(ss)) %>% 
  summarize(mean = mean(effect_10))

samplemean %>%
  group_by(factor(ss)) %>% 
  summarize(mean = mean(d_effect_10))


samplemean %>%
  group_by(factor(ss)) %>% 
  summarize(mean = mean(meant))
samplemean %>%
  group_by(factor(ss)) %>% 
  summarize(mean = mean(meanc))

samplemean %>%
  group_by(factor(ss)) %>% 
  summarize(mean = mean(V2))

samplemean %>%
  group_by(factor(ss)) %>% 
  summarize(mean = mean(V4))

samplemean %>%
  group_by(factor(ss)) %>% 
  summarize(mean = mean(V5))

samplemean %>%
  group_by(factor(ss)) %>% 
  summarize(mean = mean(V9))

samplemean %>%
  group_by(factor(ss)) %>% 
  summarize(mean = mean(V10))
            